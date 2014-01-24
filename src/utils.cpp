#include "utils.h"

#include <list>

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#include "defs.h"

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

simsys::IsoRect_2
random_window(double ratio, double xmax, double ymax)
{
  double r = 1 - std::sqrt(ratio);
  boost::random::uniform_real_distribution<> unifx(0, xmax * r), unify(0, ymax * r);
  double xmin = unifx(gen), ymin = unify(gen);
  return simsys::IsoRect_2(xmin, ymin, xmin + xmax * (1 - r), ymin + ymax * (1 - r));
}

std::vector<std::pair<simsys::Fuzzy_iso_box, double> >
intersect_room(const simsys::IsoRect_2 &win,
               const std::vector<simsys::IsoRect_2> &rooms)
{
  std::vector<std::pair<simsys::Fuzzy_iso_box, double> > results;
  for (size_t i = 0; i < rooms.size(); ++i) {
    auto res = CGAL::intersection(win, rooms[i]);
    // if the query intersects with a room, then the intersected part
    // extends to the whole room.
    if (res) {
      const simsys::IsoRect_2 tmp = *boost::get<simsys::IsoRect_2>(&*res);
      const simsys::IsoRect_2 room = rooms[i];
      results.push_back(std::make_pair(simsys::Fuzzy_iso_box(room.min(), room.max()),
                                       tmp.area() / room.area()));
    }
  }
  return results;
}

std::vector<std::pair<simsys::Fuzzy_iso_box, double> >
intersect_hall(const simsys::IsoRect_2 &win,
               const std::vector<std::pair<simsys::IsoRect_2, int> > &halls)
{
  std::vector<std::pair<simsys::Fuzzy_iso_box, double> > results;
  for (size_t i = 0; i < halls.size(); ++i) {
    auto res = CGAL::intersection(win, halls[i].first);
    if (res) {
      const simsys::IsoRect_2 tmp = *boost::get<simsys::IsoRect_2>(&*res);
      const simsys::IsoRect_2 hall = halls[i].first;
      if (0 == halls[i].second)
        results.push_back(std::make_pair(
            simsys::Fuzzy_iso_box(
                simsys::Point_2(tmp.xmin(), hall.ymin()),
                simsys::Point_2(tmp.xmax(), hall.ymax())),
            (tmp.ymax() - tmp.ymin()) / (hall.ymax() - hall.ymax())));
      else
        results.push_back(std::make_pair(
            simsys::Fuzzy_iso_box(
                simsys::Point_2(hall.xmin(), tmp.ymin()),
                simsys::Point_2(hall.xmax(), tmp.ymax())),
            (tmp.xmax() - tmp.xmin()) / (hall.xmax() - hall.xmin())));
    }
  }
  return results;
}
