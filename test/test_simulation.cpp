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
  ofstream fout("vertex.txt");
  g.print(fout);
  fout.close();
}

TEST_F(SimulationTest, particle)
{
  simulation::Particle p(g);
  p.advance(g, 100);
  p.print(cout);
}

TEST_F(SimulationTest, nearest_neighbors)
{
  using namespace simulation;
  simulation::Simulation sim(_num_object=10);
}

int main(int argc, char** argv) {
  ::testing::GTEST_FLAG(filter) = "*particle";
  // This allows the user to override the flag on the command line.
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
