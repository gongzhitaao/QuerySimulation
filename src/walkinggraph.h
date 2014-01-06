#ifndef SRC_WALKINGGRAPH_H_
#define SRC_WALKINGGRAPH_H_

#pragma once

#include <map>
#include <vector>

#include "defs.h"

namespace simsys {

class WalkingGraph
{
  typedef CGAL::Range_tree_map_traits_2<K, int> Traits;
  typedef CGAL::Range_tree_2<Traits> range_tree_t;
  typedef Traits::Key Key;
  typedef Traits::Interval Interval;

  struct Segment {
    Vertex source, target;
    double ratio, prob;
  };

  struct Reader {
    Circle_2 r;
    std::vector<Segment> segments;
  };

 public:
  WalkingGraph();

  void add_vertex(int id, double x, double y, vertex_color_enum c = HALL);
  void add_edge(int src, int des);
  void add_reader(double x, double y, double r);
  void finalize();

  const UndirectedGraph &operator() () const { return g_; }
  const CoordMap &coords() const { return coords_; }
  const WeightMap &weights() const { return weights_; }

 private:
  UndirectedGraph g_;

  ColorMap colors_;
  CoordMap coords_;
  WeightMap weights_;

  std::map<int, Vertex> vertices_;
  std::vector<Reader> readers_;

  // readers
  range_tree_t r_;
};

}

#endif  // SRC_WALKINGGRAPH_H_
