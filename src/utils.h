#ifndef SRC_UTILS_H_
#define SRC_UTILS_H_

#pragma once

#include <vector>
#include <utility>

#include "walkinggraph.h"
#include "particle.h"

void build_graph(simsys::WalkingGraph &g);

std::vector<std::vector<int> >
detect(simsys::WalkingGraph &g, const std::vector<simsys::Particle> &particles);

typedef std::map<int, std::map<int, double> > AnchorMap;

bool
predict(simsys::WalkingGraph &g, int id, const std::vector<int> &reading,
        double t, AnchorMap &anchors, int limit = 2);

std::vector<double>
range_query_hitrate_vs_windowsize(
    simsys::WalkingGraph &g,
    const std::vector<simsys::Particle> &objects,
    const std::vector<std::vector<int> > &readings);

#endif  // SRC_UTILS_H_
