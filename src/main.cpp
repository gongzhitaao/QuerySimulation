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

  simulation::Simulation sim(_num_object=200);

  sim.run(200);
  sim.snapshot(100);

  int object = sim.random_object();

  std::vector<int> nn = sim.nearest_neighbors(object, 5);
  for (size_t i = 0; i < nn.size(); ++i)
    cout << nn[i] << ' ';
  cout << endl;

  return 0;
}
