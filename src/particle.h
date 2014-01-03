#ifndef SRC_PARTICLE_H_
#define SRC_PARTICLE_H_

#include <vector>
#include <utility>  // std::pair

#include "def.h"

class WalkingGraph;

class Particle
{
 public:
  Particle(const WalkingGraph *wg);

  Point_2 advance(const WalkingGraph *wg, double t = -1);
  Point_2 pos(const WalkingGraph *wg, double t = -1) const;
  trace_t pos(const Walkinggraph *wg, double start, double duration, double step = 0.1) const;

 private:
  Vertex random_next(Vertex v, const WalkingGraph *wg) const;

  static const double TimeUnit;

  Vertex source_;
  Vertex target_;
  double p_;
  double velocity_;

  std::vector<std::pair<double, Vertex> > history_;
};

#endif  // SRC_PARTICLE_H_
