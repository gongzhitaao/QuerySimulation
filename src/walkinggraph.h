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
typedef Point_set_2::Vertex_handle  vertex_handle;

enum vertex_color_enum { HALL, DOOR, ROOM, VERTEX_COLOR_ENUM };

struct vertex_coord_t { typedef boost::vertex_property_tag kind; };
typedef boost::property<
  boost::vertex_name_t, int,
  boost::property<
    vertex_coord_t, Point_2,
    boost::property<
      boost::vertex_color_t, vertex_color_enum> > > VertexProperty;

typedef boost::property<
  boost::edge_name_t, int,
  boost::property<boost::edge_weight_t, double> > EdgeProperty;

typedef boost::undirected_graph<
  VertexProperty, EdgeProperty> UndirectedGraph;
typedef boost::graph_traits<UndirectedGraph>::vertex_descriptor Vertex;
typedef boost::graph_traits<UndirectedGraph>::edge_descriptor Edge;

typedef boost::property_map<UndirectedGraph,
                            boost::vertex_name_t>::type vertex_name_t;
typedef boost::property_map<UndirectedGraph,
                            vertex_coord_t>::type coord_map_t;
typedef boost::property_map<UndirectedGraph,
                            boost::vertex_color_t>::type color_map_t;
typedef boost::property_map<UndirectedGraph,
                            boost::edge_name_t>::type edge_name_t;
typedef boost::property_map<UndirectedGraph,
                            boost::edge_weight_t>::type weight_map_t;

typedef boost::unordered_map<
  int, std::vector<std::pair<int, double> > > anchor_map_t;

template<typename NameMap>
struct positive_index
{
  positive_index() { }
  positive_index(NameMap names)
      : names_(names) { }

  template <typename T>
  bool operator () (const T &v) const
  { return names_[v] >= 0; }

  NameMap names_;
};

Point_2
linear_interpolate(const Point_2 &p0, const Point_2 &p1, double a);

class WalkingGraph
{
 public:

  WalkingGraph();

  const Point_2 &
  coord(int v) const
  { return boost::get(coords_, vertices_.at(v)); }

  double
  weight(int u, int v) const
  { return boost::get(weights_,
                      boost::edge(vertices_.at(u),
                                  vertices_.at(v), g_).first); }

  vertex_color_enum
  color(int v) const
  { return boost::get(colors_, vertices_.at(v)); }

  template <typename Generator>
  int
  random_vertex(Generator gen) const
  { return boost::get(vnames_, boost::random_vertex(fg_, gen)); }

  int
  random_next(int cur, int pre = -1) const;

  landmark_t
  random_pos() const;

  landmark_t
  reader_pos(int i) const { return readermap_.at(i); }

  std::vector<std::pair<IsoRect_2, double> >
  random_window(double ratio) const;

  void
  insert_objects(const std::vector<landmark_t> &objects);

  void
  clear_objects();

  int
  detected_by(const landmark_t &pos, double radius);

  std::vector<int>
  nearest_neighbors(int object, int k);

  int
  align(const landmark_t &p);

  UndirectedGraph
  operator () ()
  { return g_; }

  void
  print(std::ostream &os) const;

 protected:
  enum { OBJECT_START=1000 };

  void
  initialize();

  void
  insert_anchors(double unit = 20.0);

  UndirectedGraph g_;

  vertex_name_t vnames_;
  coord_map_t coords_;
  color_map_t colors_;
  edge_name_t enames_;
  weight_map_t weights_;

  boost::filtered_graph<UndirectedGraph,
                        positive_index<edge_name_t>,
                        positive_index<vertex_name_t> > fg_;

  anchor_map_t anchors_;
  anchor_map_t objects_;

  boost::unordered_map<int, Edge> edges_;
  boost::unordered_map<int, Vertex> vertices_;

  Point_set_2 readerset_;
  boost::unordered_map<int, landmark_t> readermap_;

  std::vector<IsoRect_2> rooms_;
  std::vector<IsoRect_2> halls_;
  std::vector<int> dirs_;

  double xmax_, ymax_;
};

}  // namespace simsys

#endif  // SRC_WALKINGGRAPH_H_
