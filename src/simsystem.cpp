#include <cmath>
#include <map>

#include <boost/property_map/property_map.hpp>
#include <boost/graph/astar_search.hpp>
#include <boost/graph/random.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#include "simsystem.h"
#include "walkinggraph.h"

namespace simsys {

void SimSystem::run(const WalkingGraph *wg, double duration, int num_object)
{
  records_.clear();

  const WalkingGraph &g = *wg;
  random::mt19937 gen(time(0));

  // Human average walking speed ranging roughly from 4.5 to 5.4 km/h
  // See: https://en.wikipedia.org/wiki/Walking
  random::normal_distribution<> norm(1.5, 0.3);

  vector<Vertex> p(num_vertices(g()));
  vector<double> d(num_vertices(g()));

  for (int id = 0; id < num_object; ++id) {

    Vertex start = random_vertex(g(), gen);
    double elapsed = 0;
    double velocity = norm(gen);  // very unlikely to be zero.

    trace_t trace;

    random::uniform_int_distribution<> unifi(0, in_degree(start, g()) - 1);
    graph_traits<UndirectedGraph>::in_edge_iterator it = (in_edges(start, g())).first;
    it += unifi(gen);
    Vertex pre = source(*it, g());

    random::uniform_real_distribution<> unifd(0, 1);
    double pre_elapsed = g.weights()[*it] / velocity;
    elapsed = pre_elapsed * unifd(gen);

    trace.push_back(std::make_pair(elapsed - pre_elapsed, g.coords()[pre]));
    trace.push_back(std::make_pair(elapsed, g.coords()[start]));

    while (elapsed < duration) {
      Vertex goal = random_vertex(g(), gen);

      astar_goal_visitor<Vertex> visitor(start);
      try {
        astar_search(g(), goal, heuristic<UndirectedGraph, CoordMap>(goal, g.coords()),
                     predecessor_map(make_iterator_property_map(p.begin(), get(vertex_index, g()))).
                     distance_map(make_iterator_property_map(d.begin(), get(vertex_index, g()))).
                     visitor(visitor));
      } catch (found_goal fg) {
        for (Vertex u = start, v; u != p[u] && elapsed < duration; u = p[u]) {
          v = p[u];
          Edge e = edge(u, v, g()).first;
          double dist = g.weights()[e];
          elapsed += dist / velocity;
          trace.push_back(std::make_pair(elapsed, g.coords()[v]));
        }
      }
      start = goal;
    }
    records_.push_back(trace);
  }
}

}
