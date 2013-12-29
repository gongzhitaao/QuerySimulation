#include <gtest/gtest.h>

#include <cstdio>
#include <vector>
#include <iostream>

#include "walkinggraph.h"
#include "simsystem.h"

using namespace std;

class SimSystemTest : public ::testing::Test {
 protected:
  static void SetUpTestCase() {
    const int num_node = 10;
    for (int i = 0; i < num_node; ++i)
      g.add_vertex(i, 0, i);

    for (int i = 0; i < num_node - 1; ++i)
      g.add_edge(i, i + 1);
  }

  static simsys::WalkingGraph g;
  static simsys::SimSystem sim;
};

simsys::WalkingGraph SimSystemTest::g;
simsys::SimSystem SimSystemTest::sim;

TEST_F(SimSystemTest, run)
{
  sim.run(&g, 1.0, 1);
  sim.save_trace("test.txt");
}

TEST_F(SimSystemTest, snapshot0)
{
  sim.run(&g, 20.0, 3);
  double t = 0.5;
  vector<simsys::Point_2> result = sim.get_snapshot(t);
  char msg[64];
  cout << "ID TIME X Y" << endl;
  for (unsigned i = 0; i < result.size(); ++i) {
    sprintf(msg, "%d %.1f %.3f %.3f", i, t, result[i].x(), result[i].y());
    cout << msg << endl;
  }
}

TEST_F(SimSystemTest, snapshot1)
{
  sim.run(&g, 20.0, 3);
  vector<simsys::trace_t> result = sim.get_snapshot(2, 4, 0.5);
  char msg[64];
  cout << "ID TIME X Y" << endl;
  for (unsigned i = 0; i < result.size(); ++i) {
    for (unsigned j = 0; j < result[i].size(); ++j) {
      sprintf(msg, "%d %.1f %.3f %.3f",
              i, result[i][j].first, result[i][j].second.x(), result[i][j].second.y());
      cout << msg << endl;
    }
    cout << "====" << endl;
  }
}

int main(int argc, char** argv) {
  ::testing::GTEST_FLAG(filter) = "SimSystemTest.*";
  // This allows the user to override the flag on the command line.
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
