#ifndef SRC_UTILS_H_
#define SRC_UTILS_H_

#pragma once

#include <vector>
#include <utility>

#include "walkinggraph.h"
#include "particle.h"

void build_graph(simsys::WalkingGraph &g);

std::vector<std::vector<int> >
detect(simsys::WalkingGraph &g, const std::vector<simsys::Particle> &objects);

typedef std::map<int, std::map<int, double> > AnchorMap;

bool
predict(simsys::WalkingGraph &g, int id, const std::vector<int> &reading,
        double t, AnchorMap &anchors, int limit = 2);

enum Measurement { Precision, Recall, F1Score, MEASUREMENT };

const std::vector<std::pair<double, double> > &
measure(Measurement m);

void
range_query_windowsize(
    simsys::WalkingGraph &g,
    const std::vector<simsys::Particle> &objects,
    const std::vector<std::vector<int> > &readings,
    const std::vector<double> &winsizes);

#endif  // SRC_UTILS_H_
