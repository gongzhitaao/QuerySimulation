#include <gtest/gtest.h>

#include <fstream>
#include <iostream>

#include "param.h"
#include "simulation.h"

using namespace std;

class SimulationTest : public ::testing::Test {
 protected:
  static simulation::WalkingGraph g;
};

simulation::WalkingGraph SimulationTest::g;

TEST_F(SimulationTest, walkinggraph)
{
  ofstream fout("edge.txt");
  g.print(fout);
  fout.close();
}

TEST_F(SimulationTest, particle)
{
  simulation::Particle p(g, -1);
  p.advance(g, 100);
  p.print(cout);
}

TEST_F(SimulationTest, random_pos)
{
  for (int i = 0; i < 100; ++i) {
    simulation::landmark_t pos = g.random_pos();
    cout << pos.get<0>() << ' ' << pos.get<1>() << endl;
    cout << g.weight(pos.get<0>(), pos.get<1>()) << endl;
  }
}

TEST_F(SimulationTest, knn)
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

  sim.reset();
}

TEST_F(SimulationTest, range_query)
{
  using namespace simulation::param;
  simulation::Simulation sim(_num_object=200);

  sim.run(200);
  sim.snapshot(100);

  sim.random_window(0.01);
  std::vector<int> results = sim.range_query();
  for (size_t i = 0; i < results.size(); ++i)
    cout << results[i] << ' ';
  cout << endl;

  sim.reset();
}

int main(int argc, char** argv) {
  ::testing::GTEST_FLAG(filter) = "*range_query";
  // This allows the user to override the flag on the command line.
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
