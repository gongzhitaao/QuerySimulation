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
const double UNIT_LENGH = 1.0;   // distance between anchor points along each axis.

typedef CGAL::Range_tree_map_traits_2<simsys::K, std::pair<int, int> > Traits;
typedef CGAL::Range_tree_2<Traits> RangeTree_2;
typedef Traits::Key Node;
typedef Traits::Interval Window;

typedef std::map<std::pair<int, int>, std::map<int, double> > Anchor;

inline std::pair<int, int> key(const simsys::Point_2 p)
{
  return std::make_pair((int)(p.x() / UNIT_LENGH), (int)(p.y() / UNIT_LENGH));
}

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
    anchors[key(p)][id] += 1.0 / total;
  }

  return true;
}

int main()
{
  // number of readings
  const int NUM_READING = (int) (DURATION * RATE);
  const simsys::Line_2 Horizontal(simsys::Point_2(0, 0), simsys::Point_2(UNIT_LENGH, 0));
  const simsys::Line_2 Vertical(simsys::Point_2(0, 0), simsys::Point_2(0, UNIT_LENGH));

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
  RangeTree_2 anchortree;
  {
    std::ifstream fin("../data/edge.txt");
    int v1, v2;
    std::vector<Node> inputs;
    while (fin >> v1 >> v2) {
      simsys::Point_2 start, end;
      boost::tie(start, end) = g.add_edge(v1, v2);
      simsys::Vector_2 v = end - start;
      if (CGAL::parallel(simsys::Line_2(simsys::Point_2(0, 0), v), Vertical)) {
        int y0 = std::ceil(start.y() / UNIT_LENGH);
        int count = (end.y() - y0 * UNIT_LENGH) / UNIT_LENGH;
        for (int dy = 0; dy < count; ++dy) {
          simsys::Point_2 p = simsys::Point_2(start.x(), y0 + dy * UNIT_LENGH);
          inputs.push_back(Node(p, key(p)));
        }
      } else {
        int x0 = std::ceil(start.x() / UNIT_LENGH);
        int count = (end.x() - x0 * UNIT_LENGH) / UNIT_LENGH;
        for (int dx = 0; dx < count; ++dx) {
          simsys::Point_2 p = simsys::Point_2(x0 + dx * UNIT_LENGH, start.y());
          inputs.push_back(Node(p, key(p)));
        }
      }
      anchortree.make_tree(inputs.begin(), inputs.end());
    }
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
  }

  return 0;
}
