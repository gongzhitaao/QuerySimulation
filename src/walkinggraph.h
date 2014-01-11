#ifndef SRC_WALKINGGRAPH_H_
#define SRC_WALKINGGRAPH_H_

#pragma once

#include <map>
#include <vector>
#include <utility>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Range_segment_tree_traits.h>
#include <CGAL/Range_tree_k.h>

namespace simsys {

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Circle_2 Circle_2;
typedef K::Vector_2 Vector_2;
typedef K::Line_2 Line_2;

enum vertex_color_enum { HALL, DOOR, ROOM, VERTEX_COLOR_ENUM };

struct vertex_coord_t { typedef boost::vertex_property_tag kind; };
typedef boost::property<boost::vertex_index1_t, int> VertexIndexProperty;
typedef boost::property<vertex_coord_t, Point_2, VertexIndexProperty> CoordProperty;
typedef boost::property<boost::vertex_color_t, vertex_color_enum, CoordProperty> VertexProperty;

typedef boost::property<boost::edge_weight_t, double> EdgeProperty;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                              VertexProperty, EdgeProperty> UndirectedGraph;
typedef boost::graph_traits<UndirectedGraph>::vertex_descriptor Vertex;
typedef boost::graph_traits<UndirectedGraph>::edge_descriptor Edge;
const Vertex NullVertex = boost::graph_traits<UndirectedGraph>::null_vertex();

typedef boost::property_map<UndirectedGraph, boost::vertex_index1_t>::type VertexIndexMap;
typedef boost::property_map<UndirectedGraph, boost::vertex_color_t>::type ColorMap;
typedef boost::property_map<UndirectedGraph, vertex_coord_t>::type CoordMap;
typedef boost::property_map<UndirectedGraph, boost::edge_weight_t>::type WeightMap;

struct Reader {
  Point_2 center;
  Vertex source, target;
  double ratio;
};

class WalkingGraph
{
  typedef CGAL::Range_tree_map_traits_2<K, int> Traits;
  typedef CGAL::Range_tree_2<Traits> RangeTree_2;
  typedef Traits::Key Node;
  typedef Traits::Interval Window;

 public:
  WalkingGraph();

  void add_vertex(int id, double x, double y, vertex_color_enum c = HALL);
  std::pair<Point_2, Point_2> add_edge(int src, int des);
  void add_reader(double x, double y, int v1, int v2);
  void build_index();

  int detected(const Point_2 &p, double r, int id = -1);

  UndirectedGraph &operator() () { return g_; }
  const UndirectedGraph &operator() () const { return g_; }

  const VertexIndexMap &indices() const { return vertindices_; }
  const CoordMap &coords() const { return coords_; }
  const WeightMap &weights() const { return weights_; }
  const Reader &reader(int i) const { return readers_[i]; }

 private:
  UndirectedGraph g_;

  VertexIndexMap vertindices_;
  ColorMap colors_;
  CoordMap coords_;
  WeightMap weights_;

  std::map<int, Vertex> vertices_;
  std::vector<Reader> readers_;

  RangeTree_2 readertree_;
};

}

#endif  // SRC_WALKINGGRAPH_H_
