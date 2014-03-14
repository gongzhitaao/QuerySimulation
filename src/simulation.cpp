#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/mpl/inserter.hpp>

#include "simulation.h"
#include "global.h"

namespace simulation {

using std::cout;
using std::endl;

// Query window intersects with rooms and hall, generating windows
// with probability accordingly.
static std::vector<std::pair<IsoRect_2, double> > wins_;

// Readings for each objects
static std::vector<std::vector<int> > readings_;

static boost::unordered_map<
  int, boost::unordered_map<int, double> > probs_;

// Simulation duration
static double duration_;

// Objects point set
static Point_set_2 objectset_;

void
Simulation_impl_::initialize()
{
  g_.enter_room(enter_room_prob_);
  g_.knock_door(knock_door_prob_);

  for (int i = 0; i < num_object_; ++i)
    objects_.push_back(Particle(g_, i));
}

// Generate readings for each objects
void
Simulation_impl_::detecting()
{
  readings_.clear();
  boost::random::uniform_real_distribution<> unifd(0, 1);
  for (size_t i = 0; i < objects_.size(); ++i) {
    std::vector<int> tmp;
    for (int j = 0; j < duration_; ++j) {
      if (unifd(gen) > success_rate_) tmp.push_back(-1);
      else tmp.push_back(g_.detected_by(objects_[i].pos(g_, j),
                                        radius_));
    }
    readings_.push_back(tmp);
  }
}

landmark_t
Simulation_impl_::random_inside_reader(int i) const
{
  landmark_t pos = g_.reader_pos(i);
  Particle p(g_, -1, pos);
  return p.advance(g_, radius_);
}

void
Simulation_impl_::random_window(double ratio) const
{
  wins_ = g_.random_window(ratio);
}

void
Simulation_impl_::snapshot(double t)
{
  std::vector<std::pair<int, landmark_t> > infos;

  // construct objects point set for range query
  {
    std::vector<Point_2> points;
    for (auto it = objects_.begin(); it != objects_.end(); ++it) {
      landmark_t pos = it->pos(g_, t);
      Point_2 p = g_.coord(pos);
      infos.push_back(std::make_pair(it->id(), pos));
      points.push_back(p);
    }

    objectset_.insert(
        boost::make_zip_iterator(boost::make_tuple(points.begin(),
                                                   infos.begin())),
        boost::make_zip_iterator(boost::make_tuple(points.end(),
                                                   infos.end())));
  }

  // insert object into graph for nearest neighbor query
  {
    std::vector<landmark_t> tmp(infos.size());
    std::transform(infos.begin(), infos.end(), tmp.begin(),
                   [](const std::pair<int, landmark_t> &p) {
                     return p.second;
                   });
    g_.insert_objects(tmp);
  }

  for (int i = 0; i < num_object_; ++i)
    predicting(i, t);
}

void
Simulation_impl_::reset()
{
  probs_.clear();
  objectset_.clear();
  g_.clear_objects();
}

boost::unordered_set<int>
Simulation_impl_::range_query()
{
  std::vector<vertex_handle> tmp;
  for (auto w = wins_.begin(); w != wins_.end(); ++w) {
    objectset_.range_search(w->first.vertex(0), w->first.vertex(1),
                            w->first.vertex(2), w->first.vertex(3),
                            std::back_inserter(tmp));
  }

  std::vector<int> results;
  std::transform(tmp.begin(), tmp.end(), std::back_inserter(results),
                 [](const vertex_handle vh) {
                   return vh->info().first;
                 });
  return boost::unordered_set<int>(results.begin(), results.end());
}

boost::unordered_map<int, double>
Simulation_impl_::range_query_pred()
{
  boost::unordered_map<int, double> results;
  for (auto w = wins_.begin(); w != wins_.end(); ++w) {
    std::vector<int> anchors = g_.anchors_in_win(w->first);
    for (auto ac = anchors.begin(); ac != anchors.end(); ++ac) {
      boost::unordered_map<int, double> &prob = probs_[*ac];
      for (auto it = prob.begin(); it != prob.end(); ++it)
        results[it->first] += it->second * w->second;
    }
  }
  return results;
}

boost::unordered_set<int>
Simulation_impl_::nearest_neighbors(int id, int k)
{
  return g_.nearest_neighbors(id, k);
}

boost::unordered_map<int, double>
Simulation_impl_::nearest_neighbors_pred(int id, int k)
{
  boost::unordered_map<int, double> results =
      g_.nearest_neighbors(id, k, probs_);
  return results;
}

void
Simulation_impl_::run(double duration)
{
  duration_ = duration;
  for (auto it = objects_.begin(); it != objects_.end(); ++it)
    it->advance(g_, duration);
}

int
Simulation_impl_::random_object()
{
  boost::random::uniform_int_distribution<> unifi(0, num_object_);
  return unifi(gen);
}

bool
Simulation_impl_::predicting(int obj, double t, int limit)
{
  const std::vector<int> &reading = readings_[obj];

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
  std::list<Particle> subparticles;
  for (int i = 0; i < num_particle_; ++i)
    subparticles.push_back(
        Particle(g_, obj, random_inside_reader(reading[start])));

  // This is the filter process where we eliminate those that missed
  // the reader.  This is NOT particle filter, just an extention of
  // symbolic model IMO.
  for (int i = start + 1; i <= end; ++i) {
    for (auto it = subparticles.begin(); it != subparticles.end();
         /* empty */) {
      landmark_t p = it->advance(g_);
      if (reading[i] >= 0 && g_.detected_by(p, radius_) != reading[i])
        it = subparticles.erase(it);
      else ++it;
    }

    int sz = subparticles.size();

    if (0 == sz) return false;

    if (sz < num_particle_) {
      boost::random::uniform_int_distribution<> unifi(0, sz - 1);
      int left = num_particle_ - sz;
      for (int i = 0; i < left; ++i)
        subparticles.push_back(
            *boost::next(subparticles.begin(), unifi(gen)));
    }
  }

  // Predicting.  During the *remain*, the object's position is
  // unknown, which is exactly what we'd like to predict.
  double remain = t - end;
  double prob = 1.0 / num_particle_;
  for (auto it = subparticles.begin();
       it != subparticles.end(); ++it) {
    landmark_t p = it->advance(g_, remain);
    probs_[g_.align(p)][it->id()] += prob;
  }

  return true;
}

}
