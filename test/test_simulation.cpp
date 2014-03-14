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

TEST_F(SimulationTest, particle_advance)
{
  simulation::Particle p(g, -1);
  p.advance(g, 100);
  p.print(cout);
}

TEST_F(SimulationTest, particle_pos)
{
  simulation::Particle p(g, -1);
  p.advance(g, 100);
  for (double t = 1.0; t < 100.0; t += 1.0)
    cout << g.coord(p.pos(g, t)) << endl;
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

  boost::unordered_set<int> nn = sim.nearest_neighbors(object, 5);
  for (auto i = nn.begin(); i != nn.end(); ++i)
    cout << *i << ' ';
  cout << endl;

  sim.reset();
}

TEST_F(SimulationTest, knn_pred)
{
  using namespace simulation::param;
  simulation::Simulation sim(_num_object=200);

  sim.run(200);
  sim.snapshot(100);

  int object = sim.random_object();

  boost::unordered_map<int, double> nn =
      sim.nearest_neighbors_pred(object, 5);
  for (auto i = nn.begin(); i != nn.end(); ++i)
    cout << i->first << ' '<< i->second << endl;

  sim.reset();
}

TEST_F(SimulationTest, range_query)
{
  using namespace simulation::param;
  simulation::Simulation sim(_num_object=200);

  sim.run(200);
  sim.snapshot(100);

  sim.random_window(0.01);
  boost::unordered_set<int> results = sim.range_query();
  for (auto i = results.begin(); i != results.end(); ++i)
    cout << *i << ' ';
  cout << endl;

  sim.reset();
}

int main(int argc, char** argv) {
  ::testing::GTEST_FLAG(filter) = "*particle_pos";
  // This allows the user to override the flag on the command line.
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
