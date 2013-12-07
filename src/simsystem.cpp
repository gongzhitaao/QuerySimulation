#include <cmath>
#include <map>

#include <boost/property_map/property_map.hpp>
#include <boost/graph/astar_search.hpp>
#include <boost/graph/random.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#include "simsystem.h"

struct found_goal {};

template <typename Vertex>
class astar_goal_visitor : public boost::default_astar_visitor
{
public:
  astar_goal_visitor(Vertex goal) : goal_(goal) {}

  template <typename Graph>
  void examine_vertex(Vertex u, Graph &g) {
    if (u == goal_) throw found_goal();
  }
private:
  Vertex goal_;
};

template<typename Graph, typename Vec2>
class heuristic : public boost::astar_heuristic<Graph, double>
{
public:
  typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;

  heuristic(Vertex goal, Vec2 &coord)
    :  goal_(goal), coord_(coord) {}

  double operator()(Vertex u) {
    double dx = coord_[goal_].x - coord_[u].x;
    double dy = coord_[goal_].y - coord_[u].y;
    return sqrt(dx * dx + dy * dy);
  }

private:
  Vertex goal_;
  Vec2 coord_;
};

void SimSystem::run(double duration, int num_object, int reader_id)
{
  boost::random::mt19937 gen(time(0));

  // Human average walking speed ranging roughly from 4.5 to 5.4 km/h
  // See: https://en.wikipedia.org/wiki/Walking
  boost::random::normal_distribution<> norm(1.5, 0.3);

  // predecessor in shortest distance tree
  std::vector<Vertex> p(boost::num_vertices(g_));
  // shotest distance
  std::vector<double> d(boost::num_vertices(g_));

  // records of everyone's postions for every second.  First/row index
  // is time, second/column index is people's id'.  This is only
  // needed for particle filter.
  for (int id = 0; id < num_object; ++id) {

    Vertex start = boost::random_vertex(g_, gen);
    double elapsed = 0;
    double velocity = norm(gen);

    trace tr;
    tr.push_back(Vec3(vertex_coord_[start]));

    while (elapsed < duration) {
      Vertex goal = boost::random_vertex(g_, gen);

      astar_goal_visitor<Vertex> visitor(start);
      try {
	boost::astar_search(g_, goal, heuristic<UndirectedGraph, CoordMap>(goal, vertex_coord_),
			    boost::predecessor_map(boost::make_iterator_property_map(p.begin(), boost::get(boost::vertex_index, g_))).
			    distance_map(boost::make_iterator_property_map(d.begin(), boost::get(boost::vertex_index, g_))).
			    visitor(visitor));
      } catch (found_goal fg) {

	for (Vertex u = start, v; u != p[u] && elapsed < duration; u = p[u]) {
	  v = p[u];
	  Edge e = (boost::edge(u, v, g_)).first;
	  double dist = edge_weight_[e];
	  elapsed += dist / velocity;
	  tr.push_back(Vec3(vertex_coord_[v], elapsed));
	}
      }

      start = goal;
    }
    record_.push_back(tr);
  }
}
