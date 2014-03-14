#ifndef SRC_WALKINGGRAPH_H_
#define SRC_WALKINGGRAPH_H_

#pragma once

#include <vector>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/random.hpp>

#include <boost/tuple/tuple.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Point_set_2.h>

namespace simulation {

// <source, target, perscentage>
typedef boost::tuple<int, int, double> landmark_t;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Iso_rectangle_2 IsoRect_2;
typedef K::Circle_2 Circle_2;

typedef CGAL::Triangulation_vertex_base_with_info_2<int, K> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
typedef CGAL::Point_set_2<K, Tds> Point_set_2;
typedef Point_set_2::Vertex_handle vertex_handle;

enum vertex_color_enum { HALL, DOOR, ROOM, VERTEX_COLOR_ENUM };

struct VP {
  Point_2 coord;
  vertex_color_enum color;
};

struct EP {
  double weight;
  std::vector<std::pair<int, double> > anchors;
};

typedef boost::property<boost::vertex_name_t, int> VertexProperty;
typedef boost::property<boost::edge_name_t, int> EdgeProperty;

typedef boost::undirected_graph<
  VertexProperty, EdgeProperty> UndirectedGraph;
typedef boost::graph_traits<UndirectedGraph>::vertex_descriptor Vertex;
typedef boost::graph_traits<UndirectedGraph>::edge_descriptor Edge;

typedef boost::property_map<UndirectedGraph,
                            boost::vertex_name_t>::type vertex_name_t;
typedef boost::property_map<UndirectedGraph,
                            boost::edge_name_t>::type edge_name_t;

Point_2
linear_interpolate(const Point_2 &p0, const Point_2 &p1, double a);

class WalkingGraph
{
  friend class InsertAnchor;
  friend class RangeQuery;

  enum { ANCHORID = 100, OBJECTID = 1000};

 public:

  WalkingGraph();

  void
  enter_room(double p)
  { enter_room_ = p; }

  void
  knock_door(double p)
  { knock_door_ = p; }

  Point_2
  coord(int v) const
  { return vp_.at(v).coord; }

  Point_2
  coord(const landmark_t &pos) const
  {
    return linear_interpolate(coord(pos.get<0>()),
                              coord(pos.get<1>()),
                              pos.get<2>());
  }

  double
  weight(int u, int v) const
  {
    return ep_.at(
      boost::get(enames_,
                 boost::edge(vertices_.at(u),
                             vertices_.at(v), wg_).first)).weight;
  }

  vertex_color_enum
  color(int v) const
  { return vp_.at(v).color; }

  template <typename Generator>
  int
  random_vertex(Generator gen) const
  { return boost::get(vnames_, boost::random_vertex(wg_, gen)); }

  int
  random_next(int to, int from = -1) const;

  landmark_t
  random_pos() const;

  landmark_t
  reader_pos(int i) const { return readermap_.at(i); }

  std::vector<std::pair<IsoRect_2, double> >
  random_window(double ratio) const;

  int
  detected_by(const landmark_t &pos, double radius);

  int
  align(const landmark_t &p);

  std::vector<int>
  anchors_in_win(const IsoRect_2 &w);

  UndirectedGraph
  operator () ()
  { return wg_; }

  void
  print(std::ostream &os) const;

 protected:
  void
  initialize();

  void
  insert_anchors(double unit = 20.0);

  // walkinggraph, anchorgraph
  UndirectedGraph wg_, ag_;

  vertex_name_t vnames_;
  edge_name_t enames_;

  // vertex property
  boost::unordered_map<int, VP> vp_;
  // edge property
  boost::unordered_map<int, EP> ep_;
  // anchor poiht property
  boost::unordered_map<int, landmark_t> ap_;

  boost::unordered_map<int, Vertex> vertices_;

  Point_set_2 anchorset_;

  Point_set_2 readerset_;
  boost::unordered_map<int, landmark_t> readermap_;

  std::vector<IsoRect_2> rooms_;
  std::vector<IsoRect_2> halls_;
  std::vector<int> dirs_;

  double xmax_, ymax_;

  double enter_room_, knock_door_;
};

}  // namespace simulation

#endif  // SRC_WALKINGGRAPH_H_
