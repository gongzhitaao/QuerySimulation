#include <fstream>
#include <vector>

#include <boost/graph/random.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#include "walkinggraph.h"
#include "particle.h"

int main()
{
  const double RADIUS = 1.0;    // detection range of RFID readers
  const double DURATION = 30.0; // length of simulation
  const int NUM_OBJECT = 10;    // number of moving objects under observation
  const int NUM_PARTICLE = 32;  // number of sub-particles each object generate
  const int NUM_SNAPSHOTS = 10; // number of timestamps to query

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
    int id;
    double x, y;
    while (fin >> id >> x >> y)
      g.add_reader(x, y, RADIUS);
    fin.close();
  }

  // Run each particle along the graph for *DURATION*.
  std::vector<simsys::Particle> objects(NUM_OBJECT, simsys::Particle(g));
  for (int i = 0; i < NUM_OBJECT; ++i)
    objects[i].advance(g, DURATION);

  // Generate random timestamps for testing.
  boost::random::mt19937 gen(time(0));
  boost::random::uniform_real_distribution<> unifd(0, DURATION);
  std::vector<double> timestamps;
  for (int i = 0; i < NUM_SNAPSHOTS; ++i)
    timestamps.push_back(unifd(gen));

  return 0;
}
