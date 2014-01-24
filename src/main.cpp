#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <utility>

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iterator/zip_iterator.hpp>

#include "defs.h"
#include "walkinggraph.h"
#include "particle.h"
#include "utils.h"

using std::cout;
using std::endl;

int main()
{
  simsys::WalkingGraph g;

  const std::string Commenter("//");

  // Read in vertices for walking graph.
  {
    std::ifstream fin("../data/node.txt");
    std::string line;
    while (std::getline(fin, line)) {
      boost::algorithm::trim_left(line);
      if (line.compare(0, Commenter.size(), Commenter) == 0) continue;
      std::stringstream ss(line);
      int id, type;
      double x, y;
      ss >> id >> x >> y >> type;
      g.add_vertex(id, x, y, (simsys::vertex_color_enum) type);
    }
    fin.close();
  }

  // Read in edge for walking graph and construct anchor points.
  // Anchor points are constructed based on the assumption that
  // walking graph edge is parallel to X or Y axis.
  {
    std::ifstream fin("../data/edge.txt");
    std::string line;
    while (std::getline(fin, line)) {
      boost::algorithm::trim_left(line);
      if (line.compare(0, Commenter.size(), Commenter) == 0) continue;
      std::stringstream ss(line);
      int v1, v2;
      ss >> v1 >> v2;
      g.add_edge(v1, v2);
    }
    fin.close();
  }

  // Read in RFID readers.
  {
    std::ifstream fin("../data/rfid.txt");
    std::string line;
    while (std::getline(fin, line)) {
      boost::algorithm::trim_left(line);
      if (line.compare(0, Commenter.size(), Commenter) == 0) continue;
      std::stringstream ss(line);
      int id, v1, v2;
      double x, y;
      ss >> id >> x >> y >> v1 >> v2;
      g.add_reader(x, y, v1, v2);
    }
    fin.close();
  }

  // Read in rooms configuration.
  double xmax = 0.0, ymax = 0.0;
  std::vector<simsys::IsoRect_2> rooms;
  {
    std::ifstream fin("../data/room.txt");
    std::string line;
    while (std::getline(fin, line)) {
      boost::algorithm::trim_left(line);
      if (line.compare(0, Commenter.size(), Commenter) == 0) continue;
      std::stringstream ss(line);
      double x0, y0, x1, y1;
      int id;
      ss >> id >> x0 >> y0 >> x1 >> y1;
      rooms.push_back(simsys::IsoRect_2(x0, y0, x1, y1));
      if (x1 > xmax) xmax = x1;
      if (y1 > ymax) ymax = y1;
    }
    fin.close();
  }

  // Read in halls configuration
  std::vector<std::pair<simsys::IsoRect_2, int> > halls;
  {
    std::ifstream fin("../data/hall.txt");
    std::string line;
    while (std::getline(fin, line)) {
      boost::algorithm::trim_left(line);
      if (line.compare(0, Commenter.size(), Commenter) == 0) continue;
      std::stringstream ss(line);
      double x0, y0, x1, y1;
      int dir;
      ss >> x0 >> y0 >> x1 >> y1 >> dir;
      halls.push_back(std::make_pair(simsys::IsoRect_2(x0, y0, x1, y1), dir));
    }
    fin.close();
  }

  g.build_index(UNIT_LENGTH);

  // for (auto it = boost::vertices(g()); it.first != it.second; ++it.first) {
  //   cout << g.indices()[*(it.first)] << ": ";
  //   for (auto out = boost::out_edges(*(it.first), g()); out.first != out.second; ++out.first) {
  //     cout << "(" << g.indices()[source(*(out.first), g())]
  //          << "," << g.indices()[target(*(out.first), g())] << ") ";
  //   }
  //   cout << endl;
  // }

  // Initialize and run each particle along the graph for *DURATION*.
  std::vector<simsys::Particle> objects;
  for (int i = 0; i < NUM_OBJECT; ++i) {
    objects.push_back(simsys::Particle(g, i));
    objects[i].advance(g, DURATION);
    // objects[i].print(g);
  }

  // Generate RFID readings for each object.
  std::vector<std::vector<int> > readings = detect(g, objects);

  // for (size_t i = 0; i < readings.size(); ++i) {
  //   for (size_t j = 0; j < readings[i].size(); ++j)
  //     std::cout << j << " (" << objects[i].pos(g, j) << ") " << readings[i][j] << std::endl;
  //   std::cout << std::endl;
  // }

  const std::vector<double> ratios = {
    0.01, 0.03, 0.1,
    0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9
  };

  std::vector<double> hitrates(ratios.size());

  boost::random::uniform_real_distribution<> unifd(50.0, DURATION);
  for (int i = 0; i < NUM_TIMESTAMP; ++i) {

    int timestamp = unifd(gen);

    AnchorMap anchors;
    simsys::Tree objecttree;
    {
      std::vector<simsys::Point_2> points;
      std::vector<int> indices;
      for (int j = 0; j < NUM_OBJECT; ++j) {
        predict(g, objects[j].id(), readings[j], timestamp, anchors);
        points.push_back(objects[j].pos(g));
        indices.push_back(objects[j].id());
      }
      objecttree.insert(
          boost::make_zip_iterator(boost::make_tuple(points.begin(), indices.begin())),
          boost::make_zip_iterator(boost::make_tuple(points.end(), indices.end())));
    }

    for (size_t j = 0; j < ratios.size(); ++j) {
      for (int test = 0; test < NUM_TEST_PER_TIMESTAMP; ++test) {
        // First, adjust the query window.
        simsys::IsoRect_2 win = random_window(ratios[j], xmax, ymax);
        std::vector<std::pair<simsys::Fuzzy_iso_box, double> >
            win_rooms = intersect_room(win, rooms),
            win_halls = intersect_hall(win, halls);

        // Do the query on real data as well as fake data.
        std::vector<simsys::Point_and_int> real_results;
        std::vector<simsys::Point_and_int> enclosed_anchors;

        for (size_t k = 0; k < win_rooms.size(); ++k) {
          objecttree.search(std::back_inserter(real_results), win_rooms[k].first);
          std::vector<simsys::Point_and_int> tmp = g.anchors(win_rooms[k].first);
          enclosed_anchors.insert(enclosed_anchors.end(), tmp.begin(), tmp.end());
        }

        for (size_t k = 0; k < win_halls.size(); ++k) {
          objecttree.search(std::back_inserter(real_results), win_halls[k].first);
          std::vector<simsys::Point_and_int> tmp = g.anchors(win_halls[k].first);
          enclosed_anchors.insert(enclosed_anchors.end(), tmp.begin(), tmp.end());
        }

        std::map<int, double> fake_results;
        for (size_t k = 0; k < enclosed_anchors.size(); ++k) {
          int ind = boost::get<1>(enclosed_anchors[k]);
          for (auto it = anchors[ind].cbegin(); it != anchors[ind].cend(); ++it)
            fake_results[it->first] += it->second;
        }

        std::set<int> real;
        for (size_t k = 0; k < real_results.size(); ++k)
          real.insert(boost::get<1>(real_results[k]));

        int hit = 0;
        for (auto it = fake_results.cbegin(); it != fake_results.end(); ++it)
          if (real.end() != real.find(it->first)) ++hit;

        hitrates[j] += 1.0 * hit / real.size();
      }
    }

    // std::map<int, int> test;

    // for (auto it = anchors.begin(); it != anchors.end(); ++it) {
    //   for (auto i = it->second.begin(); i != it->second.end(); ++i) {
    //     test[i->first] += i->second;
    //   }
    // }

    // for (auto it = test.begin(); it != test.end(); ++it)
    //   cout << it->first << ' ' << it->second << endl;

  }

  for (size_t i = 0; i < hitrates.size(); ++i)
    hitrates[i] /= NUM_TEST_PER_TIMESTAMP * NUM_TIMESTAMP;

  std::ofstream of("hitrate.txt");
  for (size_t i = 0; i < hitrates.size(); ++i)
    of << ratios[i] << ' ' << hitrates[i] << endl;
  of.close();

  return 0;
}
