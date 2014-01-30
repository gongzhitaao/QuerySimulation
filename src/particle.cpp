#include "particle.h"

#include <iostream>

#include <boost/graph/random.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#include "param.h"

extern boost::random::mt19937 gen;

namespace simsys {

double Particle::unit_ = 1.0;

inline Point_2
linear_interpolate(const Point_2 &p1, const Point_2 &p2, double r)
{
  return p1 + (p2 - p1) * r;
}

Particle::Particle(const WalkingGraph &g, int id, double radius, int reader)
    : id_(id)
{
  boost::random::normal_distribution<> norm(80, 10);
  velocity_ = norm(gen);

  if (reader > 0) {
    const Reader &rd = g.reader(reader - 1);
    boost::random::uniform_real_distribution<> unifd(0, 1);
    if (unifd(gen) > 0.5) {
      source_ = rd.source;
      target_ = rd.target;
      p_ = rd.ratio;
    } else {
      source_ = rd.target;
      target_ = rd.source;
      p_ = 1 - rd.ratio;
    }

    history_.push_back(std::make_pair(-p_ * g.weight(source_, target_) / velocity_, source_));

    advance(g, unifd(gen) * radius);
    history_.clear();
    history_.push_back(std::make_pair(-p_ * g.weight(source_, target_) / velocity_, source_));

  } else {
    source_ = boost::random_vertex(g(), gen);
    target_ = random_next(source_, g);
    boost::random::uniform_real_distribution<> unifd(0.1, 0.9);
    p_ = unifd(gen);
    history_.push_back(std::make_pair(-p_ * g.weight(source_, target_) / velocity_, source_));
  }
}

Particle::Particle(const Particle &other)
    : id_(other.id_)
    , source_(other.source_)
    , target_(other.target_)
    , p_(other.p_)
    , history_(other.history_)
{
  boost::random::normal_distribution<> norm(80, 10);
  velocity_ = norm(gen);
}

// Randomly choose next vertex to advance to.  If u which is where the
// particle came from is present, then the randomly chosen vertex may
// not be u unless the out degree of v is only one in which we have no
// choice.  In this way, the particle preserves its direction.
Vertex
Particle::random_next(Vertex cur, const WalkingGraph &g, Vertex pre) const
{
  std::set<Vertex> hall, door, room;
  auto pairit = boost::out_edges(cur, g());
  for (auto it = pairit.first; it != pairit.second; ++it) {
    Vertex v = boost::target(*it, g());
    switch (g.color(v)) {
      case HALL: hall.insert(v); break;
      case DOOR: door.insert(v); break;
      case ROOM: room.insert(v); break;
      default: break;
    }
  }

  if (NullVertex == pre) {
    if (hall.size() > 0) {
      boost::random::uniform_int_distribution<> unifi(0, hall.size() - 1);
      return *(boost::next(hall.begin(), unifi(gen)));
    }
    return *(door.begin());
  }

  boost::random::uniform_real_distribution<> unifd(0, 1);

  if (room.size() > 0 && g.color(pre) != ROOM &&
      unifd(gen) < ENTER_ROOM)
    return *(room.begin());

  if (door.size() > 0
      && (g.color(cur) == ROOM ||
          unifd(gen) < KNOCK_DOOR))
      return *(door.begin());

  if (hall.size() > 1) hall.erase(pre);

  boost::random::uniform_int_distribution<> unifi(0, hall.size() - 1);
  return *(boost::next(hall.begin(), unifi(gen)));
}

Point_2
Particle::advance(const WalkingGraph &g, double duration)
{
  double w = g.weight(source_, target_);
  double elapsed = history_.back().first + w * p_ / velocity_,
            left = w * (1 - p_),
            dist = duration <= 0 ? unit_ * velocity_ - left : duration * velocity_ - left;

  while (dist >= 0) {
    elapsed += left / velocity_;
    history_.push_back(std::make_pair(elapsed, target_));

    Vertex pre = source_;
    source_ = target_;
    target_ = random_next(source_, g, pre);
    left = g.weight(source_, target_);
    dist -= left;
  }

  p_ = 1 + dist / g.weight(source_, target_);

  return linear_interpolate(g.coord(source_), g.coord(target_), p_);
}

Point_2
Particle::pos(const WalkingGraph &g, double t) const
{

  Point_2 p1, p2;
  double ratio;

  if (t >= 0) {
    int low = 0, high = history_.size() - 1, mid;

    if (t >= history_[high].first) {

      double left = (t - history_[high].first) * velocity_,
                w = g.weight(source_, target_);

      if (left > p_ * w) return g.coord(history_[high].second);
      else {
        p1 = g.coord(source_);
        p2 = g.coord(target_);
        ratio = left / w;
      }
    } else {

      while (low < high) {
        mid = (low + high) / 2;
        if (history_[mid].first > t) high = mid;
        else low = mid + 1;
      }

      Vertex u = history_[low - 1].second,
             v = history_[low].second;
      p1 = g.coord(u);
      p2 = g.coord(v);
      ratio = (t - history_[low - 1].first) * velocity_ / g.weight(u, v);
    }
  } else {
    p1 = g.coord(source_);
    p2 = g.coord(target_);
    ratio = p_;
  }

  return linear_interpolate(p1, p2, ratio);
}

std::vector<std::pair<double, Point_2> >
Particle::pos(const WalkingGraph &g, double start, double duration, int count) const
{
  double step = duration / count;
  std::vector<std::pair<double, Point_2> > res;
  for (int i = 0; i < count; ++i)
    res.push_back(std::make_pair(start + i * step, pos(g, start + i * step)));
  return res;
}

void
Particle::print(const WalkingGraph &g) const
{
  for (size_t i = 0; i < history_.size(); ++i)
    std::cout << history_[i].first << ' ' << g.label(history_[i].second) << std::endl;
}

}
