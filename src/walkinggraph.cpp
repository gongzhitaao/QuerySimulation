#include <cmath>

#include "walkinggraph.h"

namespace simsys {

int WalkingGraph::anchor_id = 0;

WalkingGraph::WalkingGraph()
    : g_(UndirectedGraph())
    , vertindices_(boost::get(boost::vertex_index1, g_))
    , colors_(boost::get(boost::vertex_color, g_))
    , coords_(boost::get(vertex_coord_t(), g_))
    , weights_(boost::get(boost::edge_weight, g_))
{
}


void WalkingGraph::add_vertex(int id, double x, double y, vertex_color_enum c)
{
  Vertex u = boost::add_vertex(g_);
  colors_[u] = c;
  coords_[u] = Point_2(x, y);
  vertices_[id] = u;
  vertindices_[u] = id;
}

void WalkingGraph::add_edge(int src, int des)
{
  Vertex u = vertices_[src];
  Vertex v = vertices_[des];
  Edge e = boost::add_edge(u, v, g_).first;
  weights_[e] = std::sqrt(CGAL::squared_distance(coords_[u], coords_[v]));
}

void WalkingGraph::add_reader(double x, double y, int v1, int v2)
{
  Reader rd;
  rd.center = Point_2(x, y);
  rd.source = vertices_[v1];
  rd.target = vertices_[v2];
  rd.ratio = std::sqrt(CGAL::squared_distance(rd.center, coords_[rd.source]) /
                       CGAL::squared_distance(coords_[rd.target], coords_[rd.source]));
  readers_.push_back(rd);
}

void WalkingGraph::build_index()
{
  std::vector<Node> inputs;
  for (size_t i = 0; i < readers_.size(); ++i)
    inputs.push_back(Node(readers_[i].center, i + 1));
  readertree_.make_tree(inputs.begin(), inputs.end());
}

int WalkingGraph::detected(const Point_2 &p, double r)
{
  double rx = r + 0.1;
  Window win(Point_2(p.x() - rx, p.y() - rx), Point_2(p.x() + rx, p.y() + rx));
  std::vector<Node> result;
  readertree_.window_query(win, std::back_inserter(result));

  // There is at most one element, if any, in the result set since
  // readers' detection ranges don't overlap.
  if (result.size() > 0)
    return result[0].second;

  return -1;
}

}
