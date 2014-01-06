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

    g.finalize();
  }

  static simsys::WalkingGraph g;
};

simsys::WalkingGraph SimSystemTest::g;

int main(int argc, char** argv) {
  ::testing::GTEST_FLAG(filter) = "";
  // This allows the user to override the flag on the command line.
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
