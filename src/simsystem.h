#ifndef SRC_SIMSYSTEM_H_
#define SRC_SIMSYSTEM_H_

#pragma once

#include <fstream>

#include <vector>
#include <string>
#include <utility>

#include "defs.h"

namespace simsys {

class WalkingGraph;

class SimSystem
{
 public:
  SimSystem() {}
  ~SimSystem() {}

  void run(const WalkingGraph *wg, double duration, int num_object);

 private:
  // Objects under observation.  Note that in particle filter process,
  // each particle will generate sub-particles for prediction.
  std::vector<Particle> particles_;
};

}

#endif  // SRC_SIMSYSTEM_H_
