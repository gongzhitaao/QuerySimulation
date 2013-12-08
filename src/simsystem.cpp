#include <cmath>
#include <map>

#include <boost/property_map/property_map.hpp>
#include <boost/graph/astar_search.hpp>
#include <boost/graph/random.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#include "simsystem.h"

namespace simsys {

using namespace std;
using namespace boost;

struct found_goal {};

template <typename Vertex>
class astar_goal_visitor : public default_astar_visitor
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
class heuristic : public astar_heuristic<Graph, double>
{
 public:
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;

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
  random::mt19937 gen(time(0));

  // Human average walking speed ranging roughly from 4.5 to 5.4 km/h
  // See: https://en.wikipedia.org/wiki/Walking
  random::normal_distribution<> norm(1.5, 0.3);

  vector<Vertex> p(num_vertices(g_));
  vector<double> d(num_vertices(g_));

  for (int id = 0; id < num_object; ++id) {

    Vertex start = random_vertex(g_, gen);
    double elapsed = 0;
    double velocity = norm(gen);  // very unlikely to be zero.

    trace tr;

    random::uniform_int_distribution<> unifi(0, in_degree(start, g_) - 1);
    graph_traits<UndirectedGraph>::in_edge_iterator it = (in_edges(start, g_)).first;
    it += unifi(gen);
    Vertex pre = source(*it, g_);

    random::uniform_real_distribution<> unifd(0, 1);
    double pre_elapsed = edge_weight_[*it] / velocity;
    elapsed = pre_elapsed * unifd(gen);

    tr.push_back(Vec3(vertex_coord_[pre], elapsed - pre_elapsed));
    tr.push_back(Vec3(vertex_coord_[start], elapsed));

    while (elapsed < duration) {
      Vertex goal = random_vertex(g_, gen);

      astar_goal_visitor<Vertex> visitor(start);
      try {
        astar_search(g_, goal, heuristic<UndirectedGraph, CoordMap>(goal, vertex_coord_),
                     predecessor_map(make_iterator_property_map(p.begin(), get(vertex_index, g_))).
                     distance_map(make_iterator_property_map(d.begin(), get(vertex_index, g_))).
                     visitor(visitor));
      } catch (found_goal fg) {

        for (Vertex u = start, v; u != p[u] && elapsed < duration; u = p[u]) {
          v = p[u];
          Edge e = (edge(u, v, g_)).first;
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

}
