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

int main(int argc, char** argv) {
  ::testing::GTEST_FLAG(filter) = "*random_pos";
  // This allows the user to override the flag on the command line.
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
