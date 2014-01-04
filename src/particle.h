#ifndef SRC_PARTICLE_H_
#define SRC_PARTICLE_H_

#include <vector>
#include <utility>  // std::pair

#include "defs.h"
#include "walkinggraph.h"

namespace simsys {

class Particle
{
 public:
  Particle(const WalkingGraph &g);

  Point_2 advance(const WalkingGraph &g, double t = -1);
  Point_2 pos(const WalkingGraph &g, double t = -1) const;
  Trace pos(const WalkingGraph &g, double start, double duration, int count) const;

 private:
  Vertex random_next(Vertex v, const WalkingGraph &g) const;

  static const double TimeUnit;

  Vertex source_;
  Vertex target_;
  double p_;
  double velocity_;

  std::vector<std::pair<double, Vertex> > history_;
};

}

#endif  // SRC_PARTICLE_H_
