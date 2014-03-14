#ifndef RANGE_QUERY_H_
#define RANGE_QUERY_H_

#pragma once

#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

#include "walkinggraph.h"

namespace simulation {

class RangeQuery
{
 public:
  RangeQuery(const WalkingGraph &g)
      : g_(g) { }

  void
  prepare(const std::vector<Point_2> &vec);

  void
  random_window(double ratio);

  boost::unordered_set<int>
  query();

  boost::unordered_map<int, double>
  predict(const boost::unordered_map<
          int, boost::unordered_map<int, double> > &probs);

 private:
  const WalkingGraph &g_;
};

}

#endif  // RANGE_QUERY_H_
