#include <boost/property_map/property_map.hpp>
#include <boost/graph/astar_search.hpp>
#include <boost/graph/random.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#include "particle.h"
#include "walkinggraph.h"

namespace simsys {

inline Point_2 linear_interpolate(const Point_2 &p1, const Point_2 &p2, double r)
{
  return p1 + (p2 - p1) * r;
}

const double Particle::TimeUnit = 1.0;

Particle::Particle(const WalkingGraph &g)
{
  boost::random::mt19937 gen(time(0));

  // Human average walking speed ranging roughly from 4.5 to 5.4 km/h
  // See: https://en.wikipedia.org/wiki/Walking
  boost::random::normal_distribution<> norm(1.5, 0.3);
  velocity_ = norm(gen);

  target_ = boost::random_vertex(g(), gen);
  source_ = random_next(target_, g);

  boost::random::uniform_real_distribution<> unifd(0, 1);
  double p_ = unifd(gen);

  double pre_elapsed = g.weights()[boost::edge(source_, target_, g()).first] * p_ / velocity_;

  history_.push_back(std::make_pair(-pre_elapsed, source_));
}

Vertex Particle::random_next(Vertex v, const WalkingGraph &g) const
{
  boost::random::mt19937 gen(time(0));
  boost::random::uniform_int_distribution<> unifi(0, boost::in_degree(v, g()) - 1);
  boost::graph_traits<UndirectedGraph>::in_edge_iterator it =
      boost::next((boost::in_edges(target_, g())).first, unifi(gen));
  return boost::source(*it, g());
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
class heuristic : public boost::astar_heuristic<Graph, double>
{
 public:
  typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;

  heuristic(Vertex &goal, const CoordMap &coord)
      :  goal_(goal), coord_(coord) {}

  double operator()(Vertex u) {
    return std::sqrt(CGAL::squared_distance(coord_[goal_], coord_[u]));
  }

 private:
  Vertex goal_;
  const CoordMap &coord_;
};

Point_2 Particle::advance(const WalkingGraph &g, double t)
{
  if (t <= 0) {
    // Particle advances for a unit of time and may not go backwards
    // unless there is a deadend.

    double elapsed = history_.back().first;
    double left = (1 - p_) * g.weights()[boost::edge(source_, target_, g()).first];
    double dist = velocity_ * TimeUnit;

    dist -= left;

    while (dist > 0) {
      elapsed += left / velocity_;
      history_.push_back(std::make_pair(elapsed, source_));

      source_ = target_;
      target_ = random_next(source_, g);
      left = g.weights()[boost::edge(source_, target_, g()).first];
      dist -= left;
    }

    p_ = 1 + dist / g.weights()[boost::edge(source_, target_, g()).first];

  } else {
    // Particle advances for a period of time set by the user, during
    // which, the particle will randomly choose its destination and
    // advance for it.  Upon reaching its destination, the particle
    // will repeat the above process till it runs out of time.  This
    // is slightly different than the *t < 0* case.

    boost::random::mt19937 gen(time(0));
    double pre_elapsed = history_.back().first;
    double elapsed = (1 - p_) * g.weights()[boost::edge(source_, target_, g()).first] / velocity_;

    std::vector<Vertex> p(num_vertices(g()));
    std::vector<double> d(num_vertices(g()));

    while (elapsed < t) {
      Vertex goal = boost::random_vertex(g(), gen);
      astar_goal_visitor<Vertex> visitor(source_);
      try {
        boost::astar_search(
            g(), goal, heuristic<UndirectedGraph, CoordMap>(goal, g.coords()),
            boost::predecessor_map(boost::make_iterator_property_map(
                p.begin(), boost::get(boost::vertex_index, g()))).
            distance_map(boost::make_iterator_property_map(
                d.begin(), boost::get(boost::vertex_index, g()))).
            visitor(visitor));
      } catch (found_goal fg) {
        for (Vertex u = source_, v; u != p[u] && elapsed < t; u = p[u]) {
          v = p[u];
          Edge e = edge(u, v, g()).first;
          double dist = g.weights()[e];
          elapsed += dist / velocity_;
          history_.push_back(std::make_pair(elapsed + pre_elapsed, u));

          source_ = u;
          target_ = v;
        }
      }
    }

    p_ = (elapsed - t) * velocity_ / g.weights()[boost::edge(source_, target_, g()).first];

  }

  return linear_interpolate(g.coords()[source_], g.coords()[target_], p_);

}

Point_2 Particle::pos(const WalkingGraph &g, double t) const
{
  Point_2 p1, p2;
  double ratio;

  if (t >= 0) {
    int low = 0, high = history_.size() - 1, mid;

    if (t >= history_[high].first) return g.coords()[history_[high].second];

    while (low < high) {
      mid = (low + high) / 2;
      if (history_[mid].first > t) high = mid;
      else low = mid + 1;
    }

    p1 = g.coords()[history_[low - 1].second];
    p2 = g.coords()[history_[low].second];
    ratio = (t - history_[low - 1].first) / (history_[low].first - history_[low - 1].first);

  } else {
    p1 = g.coords()[source_];
    p2 = g.coords()[target_];
    ratio = p_;
  }

  return linear_interpolate(p1, p2, ratio);
}

Trace Particle::pos(const WalkingGraph &g, double start, double duration, int count) const
{
  double step = duration / count;
  Trace res;
  for (int i = 0; i < count; ++i)
    res.push_back(std::make_pair(start + i * step, pos(g, start + i * step)));
  return res;
}

}
