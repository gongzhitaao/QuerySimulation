#include <vector>
#include <iostream>
#include <fstream>

// #include "param.h"
#include "simulation.h"

using std::cout;
using std::endl;

int main()
{
  using namespace simulation::param;

  using namespace simulation::param;
  simulation::Simulation sim(_num_object=200);

  sim.run(200);
  sim.snapshot(100);

  sim.random_window(0.02);
  std::vector<int> results = sim.range_query();
  for (size_t i = 0; i < results.size(); ++i)
    cout << results[i] << ' ';
  cout << endl;

  sim.reset();

  return 0;
}
