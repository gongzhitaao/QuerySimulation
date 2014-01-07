#include <fstream>
#include <vector>

#include <boost/graph/random.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#include "defs.h"
#include "walkinggraph.h"
#include "particle.h"

using simsys::gen;
extern boost::random::mt19937 gen;

const double DURATION   = 100.0; // length of simulation
const double HIT_RATE   = 0.95;  // probability of an object in detection range being detected
const int NUM_OBJECT    = 10;    // number of moving objects under observation
const int NUM_PARTICLE  = 32;    // number of sub-particles each object generate
const int NUM_SNAPSHOTS = 10;    // number of timestamps to query
const double RADIUS     = 2.0;   // detection range of RFID readers
const double RATE       = 1.0;   // reader's reading rate, i.e. reading per second

std::vector<std::vector<int> > detect(simsys::WalkingGraph &g,
                                      const std::vector<simsys::Particle> &particles,
                                      double unit, int count)
{
  boost::random::uniform_real_distribution<> unifd(0, 1);

  std::vector<std::vector<int> > readings;
  for (size_t i = 0; i < particles.size(); ++i) {
    std::vector<int> tmp;
    for (int j = 0; j < count; ++j) {
      if (unifd(gen) > HIT_RATE) tmp.push_back(-1);
      else tmp.push_back(g.detected(particles[i].pos(g, j * unit), RADIUS));
    }
    readings.push_back(tmp);
  }

  return readings;
}

int main()
{
  // number of readings
  const int NUM_READING = (int) (DURATION * RATE);

  simsys::WalkingGraph g;

  // Read in vertices for walking graph.
  {
    std::ifstream fin("../data/node.txt");
    int id, type;
    double x, y;
    while (fin >> id >> x >> y >> type)
      g.add_vertex(id, x, y, (simsys::vertex_color_enum) type);
    fin.close();
  }

  // Read in edge for walking graph.
  {
    std::ifstream fin("../data/edge.txt");
    int v1, v2;
    while (fin >> v1 >> v2)
      g.add_edge(v1, v2);
    fin.close();
  }

  // Read in RFID readers.
  {
    std::ifstream fin("../data/rfid.txt");
    int id, v1, v2;
    double x, y;
    while (fin >> id >> x >> y >> v1 >> v2)
      g.add_reader(x, y, v1, v2);
    fin.close();
  }

  g.build_index();

  std::vector<simsys::Particle> objects;
  // Run each particle along the graph for *DURATION*.
  for (int i = 0; i < NUM_OBJECT; ++i) {
    simsys::Particle p(g);
    p.advance(g, DURATION);
    objects.push_back(p);
    // p.print(g);
  }

  // Generate random timestamps for testing.
  std::vector<double> timestamps;
  {
    boost::random::uniform_real_distribution<> unifd(50.0, DURATION);
    for (int i = 0; i < NUM_SNAPSHOTS; ++i)
      timestamps.push_back(unifd(gen));
  }

  std::vector<std::vector<int> > readings = detect(g, objects, 1.0 / RATE, NUM_READING);

  // for (size_t i = 0; i < readings.size(); ++i) {
  //   for (size_t j = 0; j < readings[i].size(); ++j)
  //     std::cout << readings[i][j] << ' ';
  //   std::cout << std::endl;
  // }

  return 0;
}
