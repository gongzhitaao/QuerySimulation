#include <cmath>

#include "walkinggraph.h"

namespace simsys {

WalkingGraph::WalkingGraph()
    : g_(UndirectedGraph())
    , colors_(boost::get(boost::vertex_color, g_))
    , coords_(boost::get(vertex_coord_t(), g_))
    , weights_(boost::get(boost::edge_weight, g_))
{
}


void WalkingGraph::add_vertex(int id, double x, double y, vertex_color_enum c) {
  Vertex u = boost::add_vertex(g_);
  colors_[u] = c;
  coords_[u] = Point_2(x, y);
  vertices_[id] = u;
}

void WalkingGraph::add_edge(int src, int des) {
  Vertex u = vertices_[src];
  Vertex v = vertices_[des];
  Edge e = boost::add_edge(u, v, g_).first;
  weights_[e] = std::sqrt(CGAL::squared_distance(coords_[u], coords_[v]));
}

void WalkingGraph::add_reader(double x, double y, double r) {
  Reader rd;
  rd.r = Circle_2(Point_2(x, y), r);
  rd.segments = std::vector<Segment>();
  readers_.push_back(rd);
}

void WalkingGraph::finalize() {
  std::vector<Node> vec;
  range_tree_t verts;
  for (auto it = vertices_.cbegin(); it != vertices_.cend(); ++it)
    vec.push_back(Key(coords_[it->second], it->first));
  verts.make_tree(vec.begin(), vec.end());

  double r = std::sqrt(readers_[0].squared_radius());
}

}
