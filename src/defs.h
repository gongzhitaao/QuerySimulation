#ifndef SRC_DEFS_H_
#define SRC_DEFS_H_

#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

namespace simsys {

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Circle_2 Circle_2;

enum vertex_color_enum { ROOM, DOOR, HALL, VERTEX_COLOR_ENUM };

struct vertex_coord { typedef boost::vertex_property_tag kind; };
typedef boost::property<vertex_coord, Point_2> CoordProperty;
typedef boost::property<boost::vertex_color_t, vertex_color_enum, CoordProperty> VertexProperty;

// edge property
struct reader_list { typedef boost::vertex_property_tag kind; };
typedef boost::property<reader_list, std::vector<int> > ReaderProperty;
typedef boost::property<boost::edge_weight_t, double, ReaderProperty> EdgeProperty;

// graph type
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                              VertexProperty, EdgeProperty> UndirectedGraph;

// vertex type
typedef boost::graph_traits<UndirectedGraph>::vertex_descriptor Vertex;

// edge type
typedef boost::graph_traits<UndirectedGraph>::edge_descriptor Edge;

// get property by descriptor
typedef boost::property_map<UndirectedGraph, boost::vertex_color_t>::type ColorMap;
typedef boost::property_map<UndirectedGraph, vertex_coord>::type CoordMap;
typedef boost::property_map<UndirectedGraph, boost::edge_weight_t>::type WeightMap;
typedef boost::property_map<UndirectedGraph, reader_list>::type ReaderListMap;

// pair: (timestamp, <x, y>)
typedef std::vector<std::pair<double, Point_2> > trace_t;

}

#endif  // SRC_DEFS_H_
