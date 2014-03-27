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

#include "global.h"
#include "simulator.h"
#include "range_query.h"
#include "nearest_neighbor.h"

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

void
test_range_query()
{
  const double DURATION = 200;
  const int NUM_TIMESTAMP = 100;

  typedef boost::accumulators::accumulator_set<
    double, boost::accumulators::features<
      boost::accumulators::tag::mean,
              boost::accumulators::tag::variance> > accumulators;

  using namespace simulation::param;
  simulation::Simulator sim(_num_object=200);

  sim.run(DURATION);

  simulation::RangeQuery rq(sim);

  accumulators acc;
  boost::random::uniform_real_distribution<>
      unifd(DURATION / 4.0, DURATION * 3.0 / 4.0);

  for (int ts = 0; ts < NUM_TIMESTAMP; ++ts) {

    while (true) {
      if (rq.random_window(0.01))
        break;
    }

    rq.prepare(unifd(gen));
    boost::unordered_set<int> real = rq.query();

    if (real.empty()) {
      --ts;
      continue;
    }
    boost::unordered_map<int, double> fake = rq.predict();

    acc(recall(real, fake));
    cout << ts << '/' << NUM_TIMESTAMP << endl;
  }

  cout << boost::accumulators::mean(acc) << ' '
       << boost::accumulators::variance(acc) << endl;
}

void
test_knn()
{
  const double DURATION = 200;
  const int NUM_TIMESTAMP = 100;
  const int MAX_K = 50;

  typedef boost::accumulators::accumulator_set<
    double, boost::accumulators::features<
      boost::accumulators::tag::mean,
              boost::accumulators::tag::variance> > accumulators;

  using namespace simulation::param;
  simulation::Simulator sim(_num_object=200);

  sim.run(DURATION);

  simulation::NearestNeighbor nn(sim);

  accumulators acc[MAX_K];
  boost::random::uniform_real_distribution<>
      unifd(DURATION / 4.0, DURATION * 3.0 / 4.0);

  for (int ts = 0; ts < NUM_TIMESTAMP; ++ts) {

    nn.random_object();

    nn.prepare(unifd(gen));

    for (int k = 1; k <= MAX_K; ++k) {
      boost::unordered_set<int> real = nn.query(k);
      boost::unordered_map<int, double> fake = nn.predict(k);
      acc[k-1](recall(real, fake));
    }
    cout << ts << '/' << NUM_TIMESTAMP << endl;
  }

  std::ofstream fout("data.txt");
  for (int i = 0; i < MAX_K; ++i)
    fout << boost::accumulators::mean(acc[i]) << ' '
         << boost::accumulators::variance(acc[i]) << endl;
  fout.close();
}

void
test_threshold()
{
  const double DURATION = 200;
  const int NUM_TIMESTAMP = 100;

  typedef boost::accumulators::accumulator_set<
    double, boost::accumulators::features<
      boost::accumulators::tag::mean,
              boost::accumulators::tag::variance> > accumulators;

  using namespace simulation::param;
  simulation::Simulator sim(_num_object=200);

  sim.run(DURATION);

  simulation::NearestNeighbor nn(sim);

  accumulators acc0, acc1;
  boost::random::uniform_real_distribution<>
      unifd(DURATION / 4.0, DURATION * 3.0 / 4.0);

  for (int ts = 0; ts < NUM_TIMESTAMP; ++ts) {

    nn.random_object();

    nn.prepare(unifd(gen));

    boost::unordered_set<int> real = nn.query(3);
    boost::unordered_map<int, double> fake = nn.predict(3);

    for (auto i = fake.begin(); i != fake.end(); ++i) {
      if (real.find(i->first) != real.end())
        acc0(i->second);
      else
        acc1(i->second);
    }

    cout << ts << '/' << NUM_TIMESTAMP << endl;
  }

  cout << boost::accumulators::mean(acc0) << ' '
       << boost::accumulators::variance(acc0) << endl;

  cout << boost::accumulators::mean(acc1) << ' '
       << boost::accumulators::variance(acc1) << endl;
}

void
test_range_query_cont()
{
  const double DURATION = 200;
  const int MAX_TIMESTAMP = 50;
  const int NUM_TEST = 100;

  typedef boost::accumulators::accumulator_set<
    double, boost::accumulators::features<
      boost::accumulators::tag::mean,
              boost::accumulators::tag::variance> > accumulators;

  using namespace simulation::param;
  simulation::Simulator sim(_num_object=200);

  sim.run(DURATION);

  simulation::RangeQuery rq(sim);

  accumulators acc[MAX_TIMESTAMP];
  double time = DURATION - MAX_TIMESTAMP;
  boost::random::uniform_real_distribution<>
      unifd(time / 2.0, time);

  while (true) {
    if (rq.random_window(0.01))
      break;
  }

  for (int test = 0; test < NUM_TEST; ++test) {

    double start = unifd(gen);
    boost::unordered_set<int> back;
    for (int ts = 0; ts < MAX_TIMESTAMP; ++ts) {

      rq.prepare(start + ts);
      boost::unordered_set<int> real = rq.query();

      if (real.empty() || real == back)
        continue;

      boost::unordered_map<int, double> fake = rq.predict();

      acc[ts](recall(real, fake));

      back = real;
    }
    cout << test << '/' << NUM_TEST << endl;
  }

  std::ofstream fout("data.txt");
  for (int i = 0; i < MAX_TIMESTAMP; ++i)
    fout << boost::accumulators::mean(acc[i]) << ' '
         << boost::accumulators::variance(acc[i]) << endl;
  fout.close();
}

void
test_knn_cont()
{
  const double DURATION = 200;
  const int MAX_TIMESTAMP = 50;
  const int NUM_TEST = 100;

  typedef boost::accumulators::accumulator_set<
    double, boost::accumulators::features<
      boost::accumulators::tag::mean,
              boost::accumulators::tag::variance> > accumulators;

  using namespace simulation::param;
  simulation::Simulator sim(_num_object=200);

  sim.run(DURATION);

  simulation::NearestNeighbor nn(sim);

  accumulators acc[NUM_TEST];
  double time = DURATION - MAX_TIMESTAMP;
  boost::random::uniform_real_distribution<>
      unifd(time / 2.0, time);

  for (int test = 0; test < NUM_TEST; ++test) {

    double start = unifd(gen);
    nn.random_object();
    boost::unordered_set<int> back;

    for (int ts = 0; ts < MAX_TIMESTAMP; ++ts) {

      nn.prepare(start + ts);

      boost::unordered_set<int> real = nn.query(5);
      boost::unordered_map<int, double> fake = nn.predict(5);

      if (real != back) {
        acc[ts](recall(real, fake));
        back = real;
      }
    }

    cout << test << '/' << NUM_TEST << endl;
  }

  std::ofstream fout("data.txt");
  for (int i = 0; i < MAX_TIMESTAMP; ++i)
    fout << boost::accumulators::mean(acc[i]) << ' '
         << boost::accumulators::variance(acc[i]) << endl;
  fout.close();
}

int main()
{
  test_knn_cont();
  return 0;
}
