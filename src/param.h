#ifndef SRC_DEFS_H_
#define SRC_DEFS_H_

#pragma once

#include <vector>

#include <boost/random.hpp>
#include <boost/random/random_device.hpp>

extern boost::random::random_device rd;
extern boost::random::mt19937 gen;

extern const double SUCCESS_RATE;
extern const int DURATION;
extern const int NUM_OBJECT;
extern const int NUM_PARTICLE;
extern const int NUM_TIMESTAMP;
extern const int NUM_TEST_PER_TIMESTAMP;
extern const double RADIUS;
extern const double UNIT_LENGTH;
extern const double KNOCK_DOOR;
extern const double ENTER_ROOM;
extern const double THRESHOLD;

extern const std::vector<double> WINDOW_RATIOS;

#endif  // SRC_DEFS_H_
