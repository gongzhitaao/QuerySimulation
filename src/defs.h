#ifndef SRC_DEFS_H_
#define SRC_DEFS_H_

#pragma once

#include <vector>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Range_segment_tree_traits.h>
#include <CGAL/Range_tree_k.h>

namespace simsys {

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Circle_2 Circle_2;

typedef CGAL::Range_tree_map_traits_2<K, int> Traits;
typedef CGAL::Range_tree_2<Traits> RangeTree_2;
typedef Traits::Key Node;
typedef Traits::Interval Window;

enum vertex_color_enum { HALL, DOOR, ROOM, VERTEX_COLOR_ENUM };

struct vertex_coord_t { typedef boost::vertex_property_tag kind; };
typedef boost::property<vertex_coord_t, Point_2> CoordProperty;
typedef boost::property<boost::vertex_color_t, vertex_color_enum, CoordProperty> VertexProperty;

// Edge property.
typedef boost::property<boost::edge_weight_t, double> EdgeProperty;

// Graph type
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                              VertexProperty, EdgeProperty> UndirectedGraph;

// Vertex type
typedef boost::graph_traits<UndirectedGraph>::vertex_descriptor Vertex;

// Edge type
typedef boost::graph_traits<UndirectedGraph>::edge_descriptor Edge;

// Get property by descriptor
typedef boost::property_map<UndirectedGraph, boost::vertex_color_t>::type ColorMap;
typedef boost::property_map<UndirectedGraph, vertex_coord_t>::type CoordMap;
typedef boost::property_map<UndirectedGraph, boost::edge_weight_t>::type WeightMap;

// Pair: (timestamp, <x, y>)
typedef std::vector<std::pair<double, Point_2> > Trace;

}

#endif  // SRC_DEFS_H_
