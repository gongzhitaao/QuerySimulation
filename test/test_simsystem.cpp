#include <gtest/gtest.h>

#include <vector>
#include <iostream>

#include "simsystem.h"

using namespace std;

class SimSystemTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    const int num_node = 10;
    for (int i = 0; i < num_node; ++i)
      sim.add_node(i, 0, i, simsys::SimSystem::HALL);

    for (int i = 0; i < num_node - 1; ++i)
      sim.add_edge(i, i + 1);;
  }

  simsys::SimSystem sim;
};

TEST_F(SimSystemTest, run)
{
  sim.run(1, 1);
  sim.save_trace("test.txt");
}

TEST_F(SimSystemTest, snapshot0)
{
  sim.run(20, 3);
  double t = 0.5;
  vector<simsys::Vec3> result = sim.get_snapshot(t);
  for (unsigned i = 0; i < result.size(); ++i)
    cout << result[i].x << ' ' << result[i].y << ' ' << result[i].t << endl;
}

TEST_F(SimSystemTest, snapshot1)
{
  sim.run(20, 3);
  vector<simsys::trace> result = sim.get_snapshot(2, 4, 0.5);
  for (unsigned i = 0; i < result.size(); ++i) {
    for (unsigned j = 0; j < result[i].size(); ++j)
      cout << result[i][j].x << ' ' << result[i][j].y << ' ' << result[i][j].t << endl;
    cout << "====" << endl;
  }
}
