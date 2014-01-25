#include "utils.h"

#include <list>
#include <map>
#include <set>

#include <boost/algorithm/string.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/iterator/zip_iterator.hpp>

#include "defs.h"

// debugging
#include <fstream>

void build_graph(simsys::WalkingGraph &g)
{
  const std::string Commenter("//");
  const std::string DataRoot("/home/gongzhitaao/Documents/simsystem/data/jiao/");

  // Read in vertices for walking graph.
  {
    std::ifstream fin(DataRoot + "node.txt");
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
    std::ifstream fin(DataRoot + "edge.txt");
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
    std::ifstream fin(DataRoot + "rfid2.txt");
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
  {
    std::ifstream fin(DataRoot + "room.txt");
    std::string line;
    while (std::getline(fin, line)) {
      boost::algorithm::trim_left(line);
      if (line.compare(0, Commenter.size(), Commenter) == 0) continue;
      std::stringstream ss(line);
      double x0, y0, x1, y1;
      int id;
      ss >> id >> x0 >> y0 >> x1 >> y1;
      g.add_room(id, x0, y0, x1, y1);
    }
    fin.close();
  }

  // Read in halls configuration
  {
    std::ifstream fin(DataRoot + "hall.txt");
    std::string line;
    while (std::getline(fin, line)) {
      boost::algorithm::trim_left(line);
      if (line.compare(0, Commenter.size(), Commenter) == 0) continue;
      std::stringstream ss(line);
      double x0, y0, x1, y1;
      int dir;
      ss >> x0 >> y0 >> x1 >> y1 >> dir;
      g.add_hall(x0, y0, x1, y1, dir);
    }
    fin.close();
  }

  g.build_index(UNIT_LENGTH);
}

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
             double t, AnchorMap &anchors, int limit)
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

    int size = subparticles.size();

    if (0 == size) return false;

    if (size < NUM_PARTICLE) {
      boost::random::uniform_int_distribution<> unifi(0, subparticles.size() - 1);
      int left = NUM_PARTICLE - subparticles.size();
      for (int i = 0; i < left; ++i)
        subparticles.push_back(*boost::next(subparticles.begin(), unifi(gen)));
    }
  }

  // Predicting.  During the *remain*, the object's position is
  // unknown, which is exactly what we'd like to predict.
  double remain = t - end;
  int total = subparticles.size();
  for (auto it = subparticles.begin(); it != subparticles.end(); ++it) {
    simsys::Point_2 p = it->advance(g, remain);
    anchors[g.align(p)][it->id()] += 1.0 / total;
  }

  return true;
}

std::vector<double>
range_query_hitrate_vs_windowsize(
    simsys::WalkingGraph &g,
    const std::vector<simsys::Particle> &objects,
    const std::vector<std::vector<int> > &readings)
{
  std::vector<double> hitrates(WINDOW_RATIOS.size());

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

    for (size_t j = 0; j < WINDOW_RATIOS.size(); ++j) {
      for (int test = 0; test < NUM_TEST_PER_TIMESTAMP; ++test) {
        // First, adjust the query window.

        std::vector<std::pair<simsys::Fuzzy_iso_box, double> >
            wins = g.random_window(WINDOW_RATIOS[j]);

        // Do the query on real data as well as fake data.
        std::vector<simsys::Point_and_int> real_results;
        std::vector<simsys::Point_and_int> enclosed_anchors;

        for (size_t k = 0; k < wins.size(); ++k) {
          objecttree.search(std::back_inserter(real_results), wins[k].first);
          std::vector<simsys::Point_and_int> tmp = g.anchors(wins[k].first);
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

        hitrates[j] += real.size() > 0 ? 1.0 * hit / real.size() : 0.0;
      }
    }
  }

  for (size_t i = 0; i < hitrates.size(); ++i)
    hitrates[i] /= NUM_TEST_PER_TIMESTAMP * NUM_TIMESTAMP;

  return hitrates;
}
