#include "defs.h"

boost::random::random_device rd;
boost::random::mt19937 gen(rd());

const double HIT_RATE            = 0.99;  // probability of an object in detection range being detected
const int DURATION               = 100;   // length of simulation
const int NUM_OBJECT             = 1;     // number of moving objects under observation
const int NUM_PARTICLE           = 128;   // number of sub-particles each object Generate
const int NUM_TIMESTAMP          = 1;     // number of timestamps to test against
const int NUM_TEST_PER_TIMESTAMP = 1;     // number of tests to run for each set of parameters
const double RADIUS              = 100.0; // detection range of RFID readers
const double UNIT_LENGTH         = 50.0;  // distance between anchor points along each axis.
