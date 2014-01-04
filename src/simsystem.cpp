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

bool operator< (const std::pair<int, bool> &a, std::pair<int, bool> &b)
{
  return a.first < b.first;
}

SimSystem::SimSystem()
{
}

void SimSystem::run(const WalkingGraph &g, double duration, int num_object)
{
  particles_.clear();
  for (int i = 0; i < num_object; ++i) {
    Particle p(g);
    p.advance(g, duration);
    particles_.push_back(p);
  }
}

}
