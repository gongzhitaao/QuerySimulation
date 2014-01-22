#include <algorithm>
#include <cmath>
#include <fstream>

#include <boost/iterator/counting_iterator.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/astar_search.hpp>
#include <boost/iterator/zip_iterator.hpp>

#include <CGAL/Fuzzy_iso_box.h>

#include "walkinggraph.h"
#include "defs.h"

namespace simsys {

WalkingGraph::WalkingGraph()
    : g_(UndirectedGraph())
    , labels_(boost::get(boost::vertex_index1, g_))
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
  labels_[u] = id;
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
  Point_2 pa = coords_[rd.source];
  rd.ratio = std::sqrt(CGAL::squared_distance(pa, rd.center)) /
             weights_[boost::edge(rd.source, rd.target, g_).first];
  readers_.push_back(rd);
}

void WalkingGraph::build_index(double unit)
{
  // Build the readers' index.
  {
    std::vector<Point_2> points;
    std::vector<int> indices;
    for (size_t i = 0; i < readers_.size(); ++i) {
      points.push_back(readers_[i].center);
      indices.push_back(i);
    }
    readertree_.clear();
    readertree_.insert(
        boost::make_zip_iterator(boost::make_tuple(points.begin(), indices.begin())),
        boost::make_zip_iterator(boost::make_tuple(points.end(), indices.end())));
  }

  // Build the anchor points' index.
  {
    typedef K::Segment_2 Segment_2;
    const Segment_2 Vertical(Point_2(0, 0), Point_2(0, 1.0));

    std::vector<Point_2> points;
    std::vector<int> indices;

    boost::graph_traits<UndirectedGraph>::edge_iterator it, end;
    for (boost::tie(it, end) = boost::edges(g_); it != end; ++it) {

      Point_2 start = coords_[boost::source(*it, g_)],
                end = coords_[boost::target(*it, g_)];

      if (CGAL::parallel(Vertical, Segment_2(start, end))) {
        if (start.y() > end.y()) {
          Point_2 p = end;
          end = start;
          start = p;
        }
        int y0 = std::ceil(1.0 * start.y() / unit);
        int count = (end.y() - y0 * unit) / unit;
        for (int dy = 0; dy <= count; ++dy) {
          indices.push_back(points.size());
          points.push_back(Point_2(start.x(), y0 + dy * unit));
        }
      } else {
        if (start.x() > end.x()) {
          Point_2 p = end;
          end = start;
          start = p;
        }
        int x0 = std::ceil(1.0 * start.x() / unit);
        int count = (end.x() - x0 * unit) / unit;
        for (int dx = 0; dx <= count; ++dx) {
          indices.push_back(points.size());
          points.push_back(Point_2(x0 + dx * unit, start.y()));
        }
      }
    }
    anchortree_.clear();
    anchortree_.insert(
        boost::make_zip_iterator(boost::make_tuple(points.begin(), indices.begin())),
        boost::make_zip_iterator(boost::make_tuple(points.end(), indices.end())));
  }
}

int WalkingGraph::detected(const Point_2 &p, double r, int id)
{
  K_neighbor_search search(readertree_, p, 1);

  if (search.begin() == search.end()) return -1;

  int ind = boost::get<1>(search.begin()->first);

  return (CGAL::squared_distance(p, readers_[ind].center) <= r * r) ? ind + 1 : -1;
}

struct found_goal {};

template <typename Vertex>
class astar_goal_visitor : public boost::default_astar_visitor
{
 public:
  astar_goal_visitor(Vertex goal) : goal_(goal) {}

  template <typename Graph>
  void examine_vertex(Vertex &u, const Graph &g) {
    if (u == goal_) throw found_goal();
  }
 private:
  Vertex goal_;
};

template<typename Graph, typename CoordMap>
class heuristic : public boost::astar_heuristic<Graph, int>
{
 public:
  typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;

  heuristic(Vertex &goal, const CoordMap &coord)
      :  goal_(goal), coord_(coord) {}

  int operator()(Vertex u) {
    return std::sqrt(CGAL::squared_distance(coord_[goal_], coord_[u]));
  }

 private:
  Vertex goal_;
  const CoordMap &coord_;
};

std::vector<Vertex> WalkingGraph::path(Vertex source, Vertex target) const
{
  std::vector<Vertex> p(boost::num_vertices(g_));
  astar_goal_visitor<Vertex> visitor(source);
  try {
    boost::astar_search(
        g_, target, heuristic<UndirectedGraph, CoordMap>(target, coords_),
        boost::predecessor_map(boost::make_iterator_property_map(
            p.begin(), boost::get(boost::vertex_index, g_))).
        visitor(visitor));
  } catch (found_goal fg) {
    return p;
  }

  return std::vector<Vertex>();
}

std::vector<int> WalkingGraph::anchors(const std::pair<Point_2, Point_2> &win)
{
  // typedef CGAL::Fuzzy_iso_box<Traits> Fuzzy_iso_box;
  // std::vector<std::pair<int, Point_2> > tmp;
  // anchortree_.search(std::back_inserter(tmp), Fuzzy_iso_box(win.first, win.second));
  std::vector<int> res;
  // std::transform(tmp.begin(), tmp.end(), std::back_inserter(res),
  //                [] (const std::pair<int, Point_2> &i) { return i.first; });
  return res;
}

}
