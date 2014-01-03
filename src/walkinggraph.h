#ifndef SRC_WALKINGGRAPH_H_
#define SRC_WALKINGGRAPH_H_

#pragma once

#include <cmath>
#include <map>
#include <vector>

#include "defs.h"

#include <CGAL/Range_segment_tree_traits.h>
#include <CGAL/Range_tree_k.h>

namespace simsys {

class WalkingGraph
{
  typedef CGAL::Range_tree_map_traits_2<K, int> Traits;
  typedef CGAL::Range_tree_2<Traits> range_tree_t;
  typedef Traits::Key Key;
  typedef Traits::Interval Interval;

  struct Reader {
    Circle_2 r;
    std::vector<index_t> segments;
  };

 public:
  WalkingGraph()
      : g_(UndirectedGraph())
      , colors_(boost::get(boost::vertex_color, g_))
      , coords_(boost::get(vertex_coord_t(), g_))
      , weights_(boost::get(boost::edge_weight, g_))
  {}

  ~WalkingGraph() {}

  void add_vertex(int id, double x, double y, vertex_color_enum c = HALL) {
    Vertex u = boost::add_vertex(g_);
    colors_[u] = c;
    coords_[u] = Point_2(x, y);
    vertices_[id] = u;
  }

  void add_edge(int src, int des) {
    Vertex u = vertices_[src];
    Vertex v = vertices_[des];
    Edge e = boost::add_edge(u, v, g_).first;
    weights_[e] = std::sqrt(CGAL::squared_distance(coords_[u], coords_[v]));
  }

  void add_reader(double x, double y, double r) {
    Reader rd;
    rd.r = Circle_2(Point_2(x, y), r);
    rd.segments = std::vector<index_t>();
    readers_.push_back(rd);
  }

  void finalize() {
    std::vector<Key> vec;
    for (std::map<int, Vertex>::const_iterator it = vertices_.begin(); it != vertices_.end(); ++it)
      vec.push_back(Key(coords_[it->second], it->first));
    v_.make_tree(vec.begin(), vec.end());

    vec.clear();
    for (size_t i = 0; i < readers_.size(); ++i) {
      const Reader &r = readers_[i];
      vec.push_back(Key(r.r.center(), i));

      double x = r.r.center().x(),
             y = r.r.center().y(),
           rad = std::sqrt(r.r.squared_radius());
      Interval win(Point_2(x - rad, x + rad), Point_2(y - rad, y + rad));
      std::vector<Key> res;
      v_.window_query(win, std::back_inserter(res));
    }
    r_.make_tree(vec.begin(), vec.end());
  }

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

  // vertices
  range_tree_t v_;

  // readers
  range_tree_t r_;
};

}

#endif  // SRC_WALKINGGRAPH_H_
