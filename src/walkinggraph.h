#ifndef SRC_WALKINGGRAPH_H_
#define SRC_WALKINGGRAPH_H_

#pragma once

#include <map>
#include <vector>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/random.hpp>

#include <boost/tuple/tuple.hpp>
#include <boost/unordered_map.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Point_set_2.h>

namespace simulation {

typedef boost::tuple<int, int, double> landmark_t;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Iso_rectangle_2 IsoRect_2;
typedef K::Circle_2 Circle_2;

typedef CGAL::Triangulation_vertex_base_with_info_2<
  std::pair<int, landmark_t>, K> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
typedef CGAL::Point_set_2<K, Tds> Point_set_2;
typedef Point_set_2::Vertex_handle  Vertex_handle;

enum vertex_color_enum { HALL, DOOR, ROOM, VERTEX_COLOR_ENUM };

struct vertex_coord_t { typedef boost::vertex_property_tag kind; };
typedef boost::property<
  boost::vertex_name_t, int,
  boost::property<
    vertex_coord_t, Point_2,
    boost::property<
      boost::vertex_color_t, vertex_color_enum> > > VertexProperty;

typedef boost::property<
  boost::edge_index_t, int,
  boost::property<boost::edge_weight_t, double> > EdgeProperty;

typedef boost::undirected_graph<
  VertexProperty, EdgeProperty> UndirectedGraph;
typedef boost::graph_traits<UndirectedGraph>::vertex_descriptor Vertex;
typedef boost::graph_traits<UndirectedGraph>::edge_descriptor Edge;

typedef boost::property_map<UndirectedGraph,
                            boost::vertex_name_t>::type name_map_t;
typedef boost::property_map<UndirectedGraph,
                            vertex_coord_t>::type coord_map_t;
typedef boost::property_map<UndirectedGraph,
                            boost::vertex_color_t>::type color_map_t;
typedef boost::property_map<UndirectedGraph,
                            boost::edge_index_t>::type index_map_t;
typedef boost::property_map<UndirectedGraph,
                            boost::edge_weight_t>::type weight_map_t;

typedef boost::unordered_map<
  int, std::vector<std::pair<int, double> > > anchor_map_t;

template <typename VertexNameMap>
struct remove_anchor
{
  remove_anchor() { }
  remove_anchor(const VertexNameMap &names)
      : names_(names) { }

  template <typename Vertex>
  bool operator () (const Vertex &v) const
  { return names_[v] >= 0; }

  const VertexNameMap &names_;
};

struct all
{
};

Point_2
linear_interpolate(const Point_2 &p0, const Point_2 &p1, double a);

class WalkingGraph
{
 public:

  WalkingGraph();

  const Point_2 &
  coord(int v) const
  { return coords_[vertices_.at(v)]; }

  double
  weight(int u, int v) const
  { return weights_[boost::edge(vertices_.at(u),
                                vertices_.at(v), g_).first]; }

  vertex_color_enum
  color(int v) const
  { return colors_[vertices_.at(v)]; }

  template <typename Generator>
  int
  random_vertex(Generator gen) const
  { return names_[boost::random_vertex(g_, gen)]; }

  int
  random_next(int cur, int pre = -1) const;

  landmark_t
  random_pos() const;

  landmark_t
  reader_pos(int i) const { return readermap_.at(i); }

  std::vector<std::pair<IsoRect_2, double> >
  random_window(double ratio) const;

  int
  detected_by(const landmark_t &pos, double radius);

  UndirectedGraph
  operator () ()
  { return g_; }

  void
  print(std::ostream &os) const;

 protected:
  void
  initialize();

  void
  insert_anchors(double unit = 20.0);

  UndirectedGraph g_;
  boost::filtered_graph<
    UndirectedGraph, boost::keep_all, remove_anchor<name_map_t> > fg_;

  name_map_t names_;
  coord_map_t coords_;
  color_map_t colors_;
  index_map_t indices_;
  weight_map_t weights_;

  anchor_map_t anchors_;

  std::map<int, Vertex> vertices_;

  Point_set_2 readerset_;
  boost::unordered_map<int, landmark_t> readermap_;

  std::vector<IsoRect_2> rooms_;
  std::vector<IsoRect_2> halls_;
  std::vector<int> dirs_;

  double xmax_, ymax_;
};

}  // namespace simsys

#endif  // SRC_WALKINGGRAPH_H_
