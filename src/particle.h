#ifndef SRC_PARTICLE_H_
#define SRC_PARTICLE_H_

#include <iostream>
#include <vector>
#include <utility>  // std::pair

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#include "walkinggraph.h"

namespace simsys {

typedef std::vector<std::pair<double, Point_2> > Trace;

class Particle
{
 public:
  Particle(const WalkingGraph &g, int id = -1, double radius = -1.0, int reader = -1);
  Particle(const Particle &other);

  Point_2 advance(const WalkingGraph &g, double duration = -1.0);
  Point_2 pos(const WalkingGraph &g, double t = -1) const;
  Trace pos(const WalkingGraph &g, double start, double duration, int count) const;

  static void set_unit(double unit) { unit_ = unit; }
  static double get_unit() { return unit_; }

  void print(const WalkingGraph &g) const;

  const int id() const { return id_; }

 private:
  Vertex random_next(const Vertex v, const WalkingGraph &g, const Vertex u = NullVertex) const;

  static double unit_;

  const int id_;

  Vertex source_;
  Vertex target_;
  double velocity_;
  double p_;

  std::vector<std::pair<double, Vertex> > history_;
};

}

#endif  // SRC_PARTICLE_H_
