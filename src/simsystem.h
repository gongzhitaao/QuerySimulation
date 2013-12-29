#ifndef SRC_SIMSYSTEM_H_
#define SRC_SIMSYSTEM_H_

#pragma once

#include <fstream>

#include <vector>
#include <string>
#include <utility>

#include "defs.h"

namespace simsys {

class WalkingGraph;

class SimSystem
{
 public:
  SimSystem() {}
  ~SimSystem() {}

  void run(const WalkingGraph *g, double duration, int num_object);

  std::vector<Point_2> get_snapshot(double t) const {
    std::vector<Point_2> positions;
    for (unsigned i = 0; i < records_.size(); ++i)
      positions.push_back(get_pos_by_timestamp(records_[i], t));
    return positions;
  }

  std::vector<trace_t> get_snapshot(double start, double end, double step = 0.1) const {
    std::vector<trace_t> result;
    for (unsigned i = 0; i < records_.size(); ++i) {
      trace_t trace;
      for (double t = start; t <= end; t += step)
        trace.push_back(std::make_pair(t, get_pos_by_timestamp(records_[i], t)));
      result.push_back(trace);
    }

    return result;
  }

  std::vector<std::vector<int> > get_prior_readings(double t, int id = -1, int limit = 2) {
  }

  void load_trace(const std::string &fname) {}

  void save_trace(const std::string &fname) const {
    std::ofstream fout(fname);
    char msg[64];
    for (unsigned i = 0; i < records_.size(); ++i) {
      for (unsigned j = 0; j < records_[i].size(); ++j) {
        sprintf(msg, "%d %.6f %.6f %.6f",
                i, records_[i][j].first, records_[i][j].second.x(), records_[i][j].second.y());
        fout << msg << std::endl;
      }
    }
    fout.close();
  }

 private:
  Point_2 get_pos_by_timestamp(const trace_t &trace, double t) const {
    int low = 0, high = trace.size() - 1, mid;

    if (t >= trace[high].first) return trace[high].second;

    while (low < high) {
      mid = (low + high) / 2;
      if (trace[mid].first > t) high = mid;
      else low = mid + 1;
    }

    double ratio = (t - trace[low - 1].first) / (trace[low].first - trace[low - 1].first);
    return interpolate(trace[low - 1].second, trace[low].second, ratio);
  }

  Point_2 interpolate(const Point_2 &start, const Point_2 &target, double ratio) const {
    return start + (target - start) * ratio;
  }

  std::vector<trace_t> records_;
  std::vector<std::vector<int> > readings_;
};

}

#endif  // SRC_SIMSYSTEM_H_
