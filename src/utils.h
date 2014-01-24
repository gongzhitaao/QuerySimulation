#ifndef SRC_UTILS_H_
#define SRC_UTILS_H_

#pragma once

#include <vector>
#include <utility>

#include "walkinggraph.h"
#include "particle.h"

typedef std::map<int, std::map<int, double> > AnchorMap;

std::vector<std::vector<int> >
detect(simsys::WalkingGraph &g, const std::vector<simsys::Particle> &particles);

bool
predict(simsys::WalkingGraph &g, int id, const std::vector<int> &reading,
        double t, AnchorMap &anchors, int limit = 2);

simsys::IsoRect_2
random_window(double ratio, double xmax, double ymax);

std::vector<std::pair<simsys::Fuzzy_iso_box, double> >
intersect_room(const simsys::IsoRect_2 &win,
               const std::vector<simsys::IsoRect_2> &rooms);

std::vector<std::pair<simsys::Fuzzy_iso_box, double> >
intersect_hall(const simsys::IsoRect_2 &win,
               const std::vector<std::pair<simsys::IsoRect_2, int> > &halls);

std::vector<double>
range_query_hitrate_vs_windowsize(
    simsys::WalkingGraph &g,
    const std::vector<simsys::Particle> &objects,
    const std::vector<std::vector<int> > &readings,
    const std::vector<simsys::IsoRect_2> &rooms,
    const std::vector<std::pair<simsys::IsoRect_2, int> > &halls,
    double xmax, double ymax);

#endif  // SRC_UTILS_H_
