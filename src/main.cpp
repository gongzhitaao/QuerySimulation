#include <vector>
#include <iostream>
#include <fstream>

#include <boost/random/normal_distribution.hpp>

#include "simulation.h"
#include "global.h"

using std::cout;
using std::endl;

int main()
{
  const double DURATION = 200.0;
  const int NUM_TIMESTAMP = 100;
  const double RATIO = 0.02;

  using namespace simulation::param;
  simulation::Simulation sim(_num_object=200);

  sim.run(DURATION);

  boost::random::uniform_real_distribution<> unifd(100.0, 200.0);
  for (int ts = 0; ts < NUM_TIMESTAMP; ++ts) {

    sim.snapshot(unifd(gen));

    sim.random_window(RATIO);

    std::vector<int> real = sim.range_query();
    boost::unordered_map<int, double> fake = sim.range_query_pred();

    sim.reset();
  }

  return 0;
}
