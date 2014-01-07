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
  Particle(const WalkingGraph &g, int id = -1, int reader = -1);

  Point_2 advance(const WalkingGraph &g, double t = -1);
  Point_2 pos(const WalkingGraph &g, double t = -1) const;
  Trace pos(const WalkingGraph &g, double start, double duration, int count) const;

  void print(const WalkingGraph &g) const;

  const int id() const { return id_; }

 private:
  Vertex random_next(Vertex v, const WalkingGraph &g) const;

  static const double TimeUnit;

  Vertex source_;
  Vertex target_;
  double p_;
  double velocity_;

  const int id_;

  std::vector<std::pair<double, Vertex> > history_;
};

}

#endif  // SRC_PARTICLE_H_
