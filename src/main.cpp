#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <utility>

#include <boost/algorithm/string.hpp>

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

  std::vector<double> hitrates =
      range_query_hitrate_vs_windowsize(g, objects, readings, rooms, halls, xmax, ymax);

  std::ofstream of("hitrate.txt");
  for (size_t i = 0; i < hitrates.size(); ++i)
    of << WINDOW_RATIOS[i] << ' ' << hitrates[i] << endl;
  of.close();

  return 0;
}
