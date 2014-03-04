#ifndef SRC_PARTICLE_H_
#define SRC_PARTICLE_H_

#pragma once

#include <iostream>
#include <vector>
#include <utility>  // std::pair

#include "walkinggraph.h"

namespace simulation {

class Particle
{
 public:
  Particle(const WalkingGraph &g, int id, landmark_t pos = {0, 0, -1});
  Particle(const Particle &other);

  landmark_t
  advance(const WalkingGraph &g, double duration = -1.0);

  landmark_t
  pos(const WalkingGraph &g, double t = -1) const;

  static void
  set_unit(double unit)
  { unit_ = unit; }

  static double
  get_unit() { return unit_; }

  void
  print(std::ostream &os) const;

  const int
  id() const { return id_; }

 private:

  static double unit_;

  const int id_;

  landmark_t pos_;
  double velocity_;

  std::vector<std::pair<double, int> > history_;
};

}

#endif  // SRC_PARTICLE_H_
