#include <vector>
#include <iostream>
#include <fstream>
#include <utility>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/random/normal_distribution.hpp>

#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

#include "simulator.h"
#include "global.h"

using std::cout;
using std::endl;

static std::pair<int, int>
hit(const boost::unordered_set<int> &real,
    const boost::unordered_map<int, double> &fake)
{
  int a = 0, b = 0;
  for (auto it = fake.cbegin(); it != fake.cend(); ++it) {
    if (it->second >= 0.0) {
      ++b;
      if (real.end() != real.find(it->first))
        ++a;
    }
  }
  return std::make_pair(a, b);
}

static double
recall(const boost::unordered_set<int> &real,
       const boost::unordered_map<int, double> &fake)
{
  if (real.size() == 0) return 0.0;
  return 1.0 * hit(real, fake).first / real.size();
}

static double
precision(const boost::unordered_set<int> &real,
          const boost::unordered_map<int, double> &fake)
{
  if (real.size() == 0) return 0.0;
  std::pair<int, int> count = hit(real, fake);
  return 0 == count.second ? 0.0 : 1.0 * count.first / count.second;
}

double
f1score(const boost::unordered_set<int> &real,
        const boost::unordered_map<int, double> &fake)
{
  double a = recall(real, fake);
  double b = precision(real, fake);
  return (a + b < std::numeric_limits<double>::epsilon()) ? 0.0
      : 2.0 * b * a / (b + a);
}

int main()
{
  const double DURATION = 200.0;
  const int NUM_TIMESTAMP = 10;
  const double RATIO = 0.02;
  const double K = 5;

  using namespace simulation::param;

  simulation::Simulator sim(_num_object=200);

  sim.run(DURATION);

  // typedef boost::accumulators::accumulator_set<
  //   double, boost::accumulators::features<
  //             boost::accumulators::tag::mean,
  //             boost::accumulators::tag::variance> > accumulators;

  // using namespace simulation::param;
  // simulation::Simulation sim(_num_object=200);

  // sim.run(DURATION);
  // sim.detecting();

  // accumulators acc;
  // boost::random::uniform_real_distribution<> unifd(100.0, 200.0);
  // for (int ts = 0; ts < NUM_TIMESTAMP; ++ts) {

  //   sim.snapshot(unifd(gen));

  //   sim.random_window(0.02);

  //   boost::unordered_set<int> real = sim.range_query();
  //   boost::unordered_map<int, double> fake = sim.range_query_pred();

  //   acc(recall(real, fake));

  //   sim.reset();

  //   cout << "finishing..." << ts << endl;
  //   cout << boost::accumulators::mean(acc) << ' '
  //        << boost::accumulators::variance(acc) << endl;

  // }


  return 0;
}
