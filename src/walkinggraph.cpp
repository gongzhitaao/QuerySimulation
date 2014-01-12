#include <cmath>
#include <limits>

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#include "walkinggraph.h"

namespace simsys {

extern boost::random::mt19937 gen;

WalkingGraph::WalkingGraph()
    : g_(UndirectedGraph())
    , vertindices_(boost::get(boost::vertex_index1, g_))
    , colors_(boost::get(boost::vertex_color, g_))
    , coords_(boost::get(vertex_coord_t(), g_))
    , weights_(boost::get(boost::edge_weight, g_))
    , xmax_(std::numeric_limits<double>::min())
    , ymax_(std::numeric_limits<double>::min())
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

void WalkingGraph::add_room(double x0, double y0, double x1, double y1, int node)
{
  rooms_.push_back(std::make_pair(IsoRect_2(x0, y0, x1, y1), vertices_[node]));
  if (x1 > xmax_) xmax_ = x1;
  if (y1 > ymax_) ymax_ = y1;
}

void WalkingGraph::add_hall(double x0, double y0, double x1, double y1, int dir)
{
  halls_.push_back(std::make_pair(IsoRect_2(x0, y0, x1, y1), dir));
}

void WalkingGraph::build_index(double unit)
{
  // Build the readers' index.
  {
    typedef ReaderTraits::Key Key;
    std::vector<Key> inputs;
    for (size_t i = 0; i < readers_.size(); ++i)
      inputs.push_back(Key(readers_[i].center, i));
    readertree_.make_tree(inputs.begin(), inputs.end());
  }

  // Build the anchor points' index.
  {
    typedef K::Line_2 Line_2;
    const Line_2 Horizontal(Point_2(0, 0), Point_2(unit, 0));
    const Line_2 Vertical(Point_2(0, 0), Point_2(0, unit));

    typedef AnchorTraits::Key Key;
    std::vector<Key> inputs;
    boost::graph_traits<UndirectedGraph>::edge_iterator it, end;
    for (boost::tie(it, end) = boost::edges(g_); it != end; ++it) {
      Point_2 start = coords_[boost::source(*it, g_)],
                end = coords_[boost::target(*it, g_)];
      Vector_2 v = end - start;
      if (CGAL::parallel(Line_2(Point_2(0, 0), v), Vertical)) {
        int y0 = std::ceil(start.y() / unit);
        int count = (end.y() - y0 * unit) / unit;
        for (int dy = 0; dy < count; ++dy) {
          Point_2 p = Point_2(start.x(), y0 + dy * unit);
          inputs.push_back(Key(p, key(p, unit)));
        }
      } else {
        int x0 = std::ceil(start.x() / unit);
        int count = (end.x() - x0 * unit) / unit;
        for (int dx = 0; dx < count; ++dx) {
          Point_2 p = Point_2(x0 + dx * unit, start.y());
          inputs.push_back(Key(p, key(p, unit)));
        }
      }
    }
    anchortree_.make_tree(inputs.begin(), inputs.end());
  }
}

int WalkingGraph::detected(const Point_2 &p, double r, int id)
{
  if (id > 0)
    return CGAL::squared_distance(p, readers_[id - 1].center) <= r * r ? id - 1 : -1;

  typedef ReaderTraits::Key Key;
  typedef ReaderTraits::Interval Window;

  double rx = r + 0.1;
  Window win(Point_2(p.x() - rx, p.y() - rx), Point_2(p.x() + rx, p.y() + rx));
  std::vector<Key> result;
  readertree_.window_query(win, std::back_inserter(result));

  // There is at most one element, if any, in the result set since
  // readers' detection ranges don't overlap.
  for (size_t i = 0; i < result.size(); ++i) {
    int id = result[i].second;
    if (CGAL::squared_distance(p, readers_[id].center) <= r * r)
      return id + 1;
  }

  return -1;
}

IsoRect_2 WalkingGraph::random_window(double ratio) const
{
  double r = 1 - std::sqrt(ratio);
  boost::random::uniform_real_distribution<> unifx(0, xmax_ * r), unify(0, ymax_ * r);
  double minx = unifx(gen), miny = unify(gen);
  return IsoRect_2(minx, miny, minx + xmax_ * (1 - r), miny + ymax_ * (1 - r));
}

}
