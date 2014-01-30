#include <gtest/gtest.h>

#include <cstdio>
#include <vector>
#include <fstream>
#include <iostream>

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#include "param.h"
#include "walkinggraph.h"
#include "particle.h"
#include "utils.h"

using namespace std;

class SimSystemTest : public ::testing::Test {
 protected:
  static void SetUpTestCase() {
    build_graph(g);

    std::vector<simsys::Particle> objects;
    for (int i = 0; i < NUM_OBJECT; ++i) {
      objects.push_back(simsys::Particle(g, i));
      objects[i].advance(g, DURATION);
      // objects[i].print(g);
    }

    // Generate RFID readings for each object.
    readings = detect(g, objects);
  }

  static simsys::WalkingGraph g;
  static std::vector<simsys::Particle> objects;
  static std::vector<std::vector<int> > readings;
};

simsys::WalkingGraph SimSystemTest::g;
std::vector<simsys::Particle> SimSystemTest::objects;
std::vector<std::vector<int> > SimSystemTest::readings;

TEST_F(SimSystemTest, particle)
{
  boost::random::uniform_real_distribution<> unifd(50.0, DURATION);
  char fname[32];
  for (int i = 0; i < 10; ++i) {
    AnchorMap anchors;
    double t = unifd(gen);
    predict(g, objects[0].id(), readings[0], t, anchors);

    {
      sprintf(fname, "%d-0.txt", i);
      ofstream fout(fname);
      fout << objects[0].pos(g, t) << endl;
      fout.close();
    }

    {
      sprintf(fname, "%d-1.txt", i);
      ofstream fout(fname);
      for (auto it = anchors.begin(); it != anchors.end(); ++it)
        fout << it->first << endl;
      fout.close();
    }
  }
}

int main(int argc, char** argv) {
  ::testing::GTEST_FLAG(filter) = "*particle";
  // This allows the user to override the flag on the command line.
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
