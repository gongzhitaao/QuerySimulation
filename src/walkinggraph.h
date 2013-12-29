#ifndef SRC_WALKINGGRAPH_H_
#define SRC_WALKINGGRAPH_H_

#pragma once

#include <cmath>
#include <map>
#include <vector>

#include "defs.h"

namespace simsys {

class WalkingGraph
{
 public:
  WalkingGraph()
      : g_(UndirectedGraph())
      , colors_(boost::get(boost::vertex_color, g_))
      , coords_(boost::get(vertex_coord(), g_))
      , weights_(boost::get(boost::edge_weight, g_))
      , readers_(boost::get(reader_list(), g_))
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

  void add_reader(int id, double radius, double x, double y, int v1,  int v2) {
    Edge e = boost::edge(vertices_[v1], vertices_[v2], g_);
    readers_[id] = std::make_pair(Circle_2(Point_2(x, y), radius * radius), e);
    readerlists_[e].push_back(id);
  }

  const UndirectedGraph &operator() () const { return g_; }
  const CoordMap &coords() const { return coords_; }
  const WeightMap &weights() const { return weights_; }

 private:
  UndirectedGraph g_;

  ColorMap colors_;
  CoordMap coords_;
  WeightMap weights_;
  ReaderListMap readerlists_;

  std::map<int, Vertex> vertices_;
  std::map<int, std::pair<Circle_2, Edge> > readers_;
};

}

#endif  // SRC_WALKINGGRAPH_H_
