#ifndef SRC_WALKINGGRAPH_H_
#define SRC_WALKINGGRAPH_H_

#pragma once

#include <map>
#include <vector>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <CGAL/Cartesian.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/property_map.h>
#include <CGAL/Fuzzy_iso_box.h>

namespace simsys {

typedef CGAL::Cartesian<double> K;
typedef K::Point_2 Point_2;
typedef CGAL::Iso_rectangle_2<K> IsoRect_2;

enum vertex_color_enum { HALL, DOOR, ROOM, VERTEX_COLOR_ENUM };

struct vertex_coord_t { typedef boost::vertex_property_tag kind; };
typedef boost::property<boost::vertex_index1_t, int> VertexLabelProperty;
typedef boost::property<vertex_coord_t, Point_2, VertexLabelProperty> CoordProperty;
typedef boost::property<boost::vertex_color_t, vertex_color_enum, CoordProperty> VertexProperty;

typedef boost::property<boost::edge_weight_t, double> EdgeProperty;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                              VertexProperty, EdgeProperty> UndirectedGraph;
typedef boost::graph_traits<UndirectedGraph>::vertex_descriptor Vertex;
typedef boost::graph_traits<UndirectedGraph>::edge_descriptor Edge;

typedef boost::property_map<UndirectedGraph, boost::vertex_index1_t>::type VertexLabelMap;
typedef boost::property_map<UndirectedGraph, boost::vertex_color_t>::type ColorMap;
typedef boost::property_map<UndirectedGraph, vertex_coord_t>::type CoordMap;
typedef boost::property_map<UndirectedGraph, boost::edge_weight_t>::type WeightMap;

const Vertex NullVertex = boost::graph_traits<UndirectedGraph>::null_vertex();

typedef boost::tuple<Point_2, int> Point_and_int;
typedef CGAL::Search_traits_2<K> Traits_base;
typedef CGAL::Search_traits_adapter<
  Point_and_int, CGAL::Nth_of_tuple_property_map<0, Point_and_int>, Traits_base> Traits;

typedef CGAL::Orthogonal_k_neighbor_search<Traits> K_neighbor_search;
typedef K_neighbor_search::Tree Tree;
typedef CGAL::Fuzzy_iso_box<Traits> Fuzzy_iso_box;

struct Reader {
  Point_2 center;
  Vertex source, target;
  double ratio;
};

class WalkingGraph
{
 public:
  WalkingGraph();

  void add_vertex(int id, double x, double y, vertex_color_enum c = HALL);
  void add_edge(int src, int des);
  void add_reader(double x, double y, int v1, int v2);

  void build_index(double unit);
  int detected(const Point_2 &p, double r, int id = -1);

  int align(const Point_2 &p);
  std::vector<Point_and_int> anchors(const Fuzzy_iso_box &win);
  std::vector<Vertex> path(Vertex source, Vertex target) const;

  UndirectedGraph &operator() () { return g_; }
  const UndirectedGraph &operator() () const { return g_; }

  int label(Vertex v) const { return labels_[v]; }
  const Point_2 &coord(Vertex v) const { return coords_[v]; }
  double weight(Vertex u, Vertex v) const { return weights_[boost::edge(u, v, g_).first]; }
  const Reader &reader(int i) const { return readers_[i]; }
  vertex_color_enum color(Vertex v) const { return colors_[v]; }

 private:
  UndirectedGraph g_;

  VertexLabelMap labels_;
  ColorMap colors_;
  CoordMap coords_;
  WeightMap weights_;

  std::map<int, Vertex> vertices_;
  std::vector<Reader> readers_;

  Tree readertree_;
  Tree anchortree_;
};

}

#endif  // SRC_WALKINGGRAPH_H_
