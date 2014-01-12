#include <fstream>
#include <list>
#include <vector>
#include <utility>

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
const int NUM_OBJECT    = 1;     // number of moving objects under observation
const int NUM_PARTICLE  = 128;   // number of sub-particles each object generate
const int NUM_SNAPSHOTS = 10;    // number of timestamps to query
const double RADIUS     = 2.0;   // detection range of RFID readers
const double RATE       = 1.0;   // reader's reading rate, i.e. reading per second
const double UNIT_LENGTH = 1.0;   // distance between anchor points along each axis.

typedef std::map<std::pair<int, int>, std::map<int, double> > Anchor;

// Generate readings for each objects
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

bool predict(simsys::WalkingGraph &g, int id, const std::vector<int> &reading,
             double t, Anchor &anchors, int limit = 2)
{
  // The number of valid readings, i.e. reading >= 0, in [start, end]
  // is *limit*.
  int end = t * RATE;
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
    for (auto it = subparticles.begin(); it != subparticles.end(); /* empt */) {
      simsys::Point_2 p = it->advance(g);
      if (reading[i] >= 0 && g.detected(p, RADIUS, reading[i]) < 0)
        it = subparticles.erase(it);
      else ++it;
    }
  }

  // Predicting.  During the *remain*, the object's position is
  // unknown, which is exactly what we'd like to predict.
  double remain = t - end / RATE;
  int total = subparticles.size();
  for (auto it = subparticles.begin(); it != subparticles.end(); ++it) {
    simsys::Point_2 p = it->advance(g, remain);
    anchors[simsys::key(p, UNIT_LENGTH)][id] += 1.0 / total;
  }

  return true;
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

  // Read in edge for walking graph and construct anchor points.
  // Anchor points are constructed based on the assumption that
  // walking graph edge is parallel to X or Y axis.
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

  // Read in rooms configuration.
  {
    std::ifstream fin("../data/room.txt");
    double x0, y0, x1, y1;
    int id, vertex;
    while (fin >> id >> x0 >> y0 >> x1 >> y1 >> vertex)
      g.add_room(x0, y0, x1, y1, vertex);
    fin.close();
  }

  // Read in halls configuration
  {
    std::ifstream fin("../data/hall.txt");
    double x0, y0, x1, y1;
    int dir;
    while (fin >> x0 >> y0 >> x1 >> y1 >> dir)
      g.add_hall(x0, y0, x1, y1, dir);
    fin.close();
  }

  g.build_index(UNIT_LENGTH);

  // Initialize and run each particle along the graph for *DURATION*.
  simsys::Particle::set_unit(1.0 / RATE);
  std::vector<simsys::Particle> objects;
  for (int i = 0; i < NUM_OBJECT; ++i) {
    simsys::Particle p(g);
    p.advance(g, DURATION);
    objects.push_back(p);
    // p.print(g);
  }

  // Generate RFID readings for each object.
  std::vector<std::vector<int> > readings = detect(g, objects, 1.0 / RATE, NUM_READING);
  // for (size_t i = 0; i < readings.size(); ++i) {
  //   for (size_t j = 0; j < readings[i].size(); ++j)
  //     std::cout << j / RATE << " (" << objects[i].pos(g, j / RATE) << ") " << readings[i][j] << std::endl;
  //   std::cout << std::endl;
  // }

  // Align subparticles to anchor points.
  Anchor anchors;
  {
    boost::random::uniform_real_distribution<> unifd(50.0, DURATION);
    double time = unifd(gen);
    for (int i = 0; i < NUM_OBJECT; ++i)
      predict(g, i, readings[i], time, anchors);

    simsys::IsoRect_2 win = g.random_window(0.2);
    // std::cout << win << std::endl;
  }

  return 0;
}
