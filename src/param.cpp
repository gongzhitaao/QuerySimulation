#include "param.h"

boost::random::random_device rd;
boost::random::mt19937 gen(rd());

const double HIT_RATE            = 0.99;  // probability of an object in detection range being detected
const int DURATION               = 100;   // length of simulation
const int NUM_OBJECT             = 100;   // number of moving objects under observation
const int NUM_PARTICLE           = 128;   // number of sub-particles each object Generate
const int NUM_TIMESTAMP          = 10;    // number of timestamps to test against
const int NUM_TEST_PER_TIMESTAMP = 10;    // number of tests to run for each set of parameters
const double RADIUS              = 120.0; // detection range of RFID readers
const double UNIT_LENGTH         = 20.0;  // distance between anchor points along each axis.
const double KNOCK_DOOR          = 0.1;   // probability of entering a room
const double ENTER_ROOM          = 0.1;   // probability of entering a room
const double THRESHOLD           = 0.4;   // objects with probability above this is considered

const std::vector<double> WINDOW_RATIOS = {
  0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
  0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9
};
