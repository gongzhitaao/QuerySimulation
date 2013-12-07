#include <gtest/gtest.h>

#include "simsystem.h"

class SimSystemTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    const int num_node = 10;
    for (int i = 0; i < num_node; ++i)
      sim_ring.add_node(i, 0, i, SimSystem::HALL);

    for (int i = 0; i < num_node - 1; ++i)
      sim_ring.add_edge(i, i + 1);;
  }

  SimSystem sim_ring;
};

TEST_F(SimSystemTest, add_node)
{
  sim_ring.run(20, 1);
  sim_ring.save_trace("test.txt");
}
