#ifndef SRC_SIMSYSTEM_H_
#define SRC_SIMSYSTEM_H_

#pragma once

#include <set>
#include <string>
#include <utility>
#include <vector>

#include "defs.h"
#include "walkinggraph.h"
#include "particle.h"

namespace simsys {


bool operator< (const std::pair<int, bool> &a, std::pair<int, bool> &b);

class SimSystem
{
 public:
  SimSystem();

  void run(const WalkingGraph &g, double duration, int num_object);
  // void snapshot(double t);

  // std::set<int> window_query_real(Window win, double t);
  // std::set<int> window_query_prediction(Window win, double t);

  // std::set<int> nearest_neighbor_real(int id, int k, double t);
  // std::set<int> nearest_neighbor_prediction(int id, int k, double t);

  // std::set<std::pair<int, bool> > window_query_real(Window win, double start, double duration);
  // std::set<std::pair<int, bool> > window_query_prediction(Window win, double start, double duration);

  // std::set<std::pair<int, bool> > nearest_neighbor_real(int id, int k, double start, double duration);
  // std::set<std::pair<int, bool> > nearest_neighbor_prediction(int id, int k, double start, double duration);

 private:
  // Objects under observation.  Note that in particle filter process,
  // each particle will generate sub-particles for prediction.
  std::vector<Particle> particles_;
};

}

#endif  // SRC_SIMSYSTEM_H_
