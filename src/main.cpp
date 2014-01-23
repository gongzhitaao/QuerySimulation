#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <sstream>
#include <vector>
#include <utility>

#include <CGAL/Range_segment_tree_traits.h>
#include <CGAL/Range_tree_k.h>

#include <boost/graph/random.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/algorithm/string.hpp>

#include "defs.h"
#include "walkinggraph.h"
#include "particle.h"

using std::cout;
using std::endl;

typedef std::map<std::map<int, double> > AnchorMap;

typedef CGAL::Range_tree_map_traits_2<simsys::K, int> Traits;
typedef CGAL::Range_tree_2<Traits> RangeTree;
typedef Traits::Key Key;
typedef Traits::Interval Interval;

// Generate readings for each objects
std::vector<std::vector<int> > detect(simsys::WalkingGraph &g,
                                      const std::vector<simsys::Particle> &particles)
{
  boost::random::uniform_real_distribution<> unifd(0, 1);
  std::vector<std::vector<int> > readings;
  for (size_t i = 0; i < particles.size(); ++i) {
    std::vector<int> tmp;
    for (int j = 0; j < DURATION; ++j) {
      if (unifd(gen) > HIT_RATE) tmp.push_back(-1);
      else tmp.push_back(g.detected(particles[i].pos(g, j), RADIUS));
    }
    readings.push_back(tmp);
  }
  return readings;
}

bool predict(simsys::WalkingGraph &g, int id, const std::vector<int> &reading,
             double t, AnchorMap &anchors, int limit = 2)
{
  // The number of valid readings, i.e. reading >= 0, in [start, end]
  // is *limit*.
  int end = t;
  int start = end;
  {
    int last = -1;
    int count = 0;
    for (/* empty */; count < limit && start >= 0; --start) {
      if (reading[start] >= 0 && reading[start] != last) {
        ++count;
        last = reading[start];
      }
    }

    // Not enough observation
    if (count < limit) return false;

    // The last decrement is uncalled for.
    ++start;
  }

  // Initialize subparticles
  std::list<simsys::Particle> subparticles;
  for (int i = 0; i < NUM_PARTICLE; ++i)
    subparticles.push_back(simsys::Particle(g, id, RADIUS, reading[start]));

  // This is the filter process where we eliminate those that missed
  // the reader.  This is NOT particle filter, just an extention of
  // symbolic model IMO.
  for (int i = start + 1; i <= end; ++i) {
    for (auto it = subparticles.begin(); it != subparticles.end(); /* empty */) {
      simsys::Point_2 p = it->advance(g);
      if (reading[i] >= 0 && g.detected(p, RADIUS, reading[i]) < 0)
        it = subparticles.erase(it);
      else ++it;
    }
  }

  // Predicting.  During the *remain*, the object's position is
  // unknown, which is exactly what we'd like to predict.
  double remain = t - end;
  int total = subparticles.size();
  for (auto it = subparticles.begin(); it != subparticles.end(); ++it) {
    it->advance(g, remain);
  }

  return true;
}

simsys::IsoRect_2
random_window(double ratio, double xmax, double ymax)
{
  double r = 1 - std::sqrt(ratio);
  boost::random::uniform_real_distribution<> unifx(0, xmax * r), unify(0, ymax * r);
  double xmin = unifx(gen), ymin = unify(gen);
  return simsys::IsoRect_2(xmin, ymin, xmin + xmax * (1 - r), ymin + ymax * (1 - r));
}

std::vector<std::pair<Interval, double> >
intersect_room(const simsys::IsoRect_2 &win,
               const std::vector<simsys::IsoRect_2> &rooms)
{
  std::vector<std::pair<Interval, double> > results;
  for (size_t i = 0; i < rooms.size(); ++i) {
    auto res = CGAL::intersection(win, rooms[i]);
    // if the query intersects with a room, then the intersected part
    // extends to the whole room.
    if (res) {
      const simsys::IsoRect_2 tmp = *boost::get<simsys::IsoRect_2>(&*res);
      const simsys::IsoRect_2 room = rooms[i];
      results.push_back(std::make_pair(Interval(room.min(), room.max()),
                                       tmp.area() / room.area()));
    }
  }
  return results;
}

std::vector<std::pair<Interval, double> >
intersect_hall(const simsys::IsoRect_2 &win,
               const std::vector<std::pair<simsys::IsoRect_2, int> > &halls)
{
  std::vector<std::pair<Interval, double> > results;
  for (size_t i = 0; i < halls.size(); ++i) {
    auto res = CGAL::intersection(win, halls[i].first);
    if (res) {
      const simsys::IsoRect_2 tmp = *boost::get<simsys::IsoRect_2>(&*res);
      const simsys::IsoRect_2 hall = halls[i].first;
      if (0 == halls[i].second)
        results.push_back(std::make_pair(
            Interval(simsys::Point_2(tmp.xmin(), hall.ymin()),
                     simsys::Point_2(tmp.xmax(), hall.ymax())),
            (tmp.ymax() - tmp.ymin()) / (hall.ymax() - hall.ymax())));
      else
        results.push_back(std::make_pair(
            Interval(simsys::Point_2(hall.xmin(), tmp.ymin()),
                     simsys::Point_2(hall.xmax(), tmp.ymax())),
            (tmp.xmax() - tmp.xmin()) / (hall.xmax() - hall.xmin())));
    }
  }
  return results;
}

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
    objects[i].print(g);
  }

  // Generate RFID readings for each object.
  std::vector<std::vector<int> > readings = detect(g, objects);

  for (size_t i = 0; i < readings.size(); ++i) {
    for (size_t j = 0; j < readings[i].size(); ++j)
      std::cout << j << " (" << objects[i].pos(g, j) << ") " << readings[i][j] << std::endl;
    std::cout << std::endl;
  }

  /*

  const std::vector<double> ratios = {
    0.01, 0.03, 0.1,
    0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9
  };

  std::vector<int> hitrates(ratios.size());

  boost::random::uniform_real_distribution<> unifd(50.0, DURATION);
  for (int i = 0; i < NUM_TIMESTAMP; ++i) {
    int timestamp = unifd(gen);
    Anchor anchors;
    RangeTree objecttree;
    {
      std::vector<Key> inputs;
      for (int j = 0; j < NUM_OBJECT; ++j) {
        predict(g, objects[j].id(), readings[j], timestamp, anchors);
        inputs.push_back(Key(objects[j].pos(g), objects[j].id()));
      }
      objecttree.make_tree(inputs.begin(), inputs.end());
    }

    for (size_t j = 0; j < ratios.size(); ++j) {
      for (int test = 0; test < NUM_TESTS; ++test) {
        // First, adjust the query window.
        simsys::IsoRect_2 win = random_window(ratios[j], xmax, ymax);
        std::vector<std::pair<Interval, int> >
            win_rooms = intersect_room(win, rooms),
            win_halls = intersect_hall(win, halls);

        // Do the query on real data as well as fake data.
        std::vector<Key> real_results;
        std::vector<simsys::AnchorKey> tmp_results;

        for (size_t k = 0; k < win_rooms.size(); ++k) {
          objecttree.window_query(win_rooms[k].first, std::back_inserter(real_results));
          g.anchortree().window_query(win_rooms[k].first, std::back_inserter(tmp_results));
        }

        for (size_t k = 0; k < win_halls.size(); ++k) {
          objecttree.window_query(win_halls[k].first, std::back_inserter(real_results));
          g.anchortree().window_query(win_halls[k].first, std::back_inserter(tmp_results));
        }

        std::map<int, int> fake_results;
        for (auto it = tmp_results.cbegin(); it != tmp_results.end(); ++it) {
          std::pair<int, int> ind = it->second;
          for (auto it = anchors[ind].cbegin(); it != anchors[ind].cend(); ++it)
            fake_results[it->first] += it->second;
        }

        std::set<int> real;
        for (size_t k = 0; k < real_results.size(); ++k)
          real.insert(real_results[k].second);

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
    hitrates[i] /= NUM_TESTS * NUM_TIMESTAMP;

  std::ofstream of("hitrate.txt");
  for (size_t i = 0; i < hitrates.size(); ++i)
    of << ratios[i] << ' ' << hitrates[i] << endl;
  of.close();

  */

  return 0;
}
