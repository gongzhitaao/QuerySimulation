#ifndef NEAREST_NEIGHBOR_H_
#define NEAREST_NEIGHBOR_H_

#pragma once

#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

#include "simulator.h"

class NearestNeighbor
{
 public:
  NearestNeighbor(Simulator &sim)
      : sim_(sim) { }

  void
  prepare(double t);

  boost::unordered_set<int>
  query();

  boost::unordered_map<int, double>
  predict();

 private:
  Simulator &sim_;
};

#endif  // NEAREST_NEIGHBOR_H_
