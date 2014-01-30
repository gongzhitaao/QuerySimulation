#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <utility>

#include "param.h"
#include "walkinggraph.h"
#include "particle.h"
#include "utils.h"

using std::cout;
using std::endl;

int main()
{
  simsys::WalkingGraph g;

  build_graph(g);

  // for (auto it = boost::vertices(g()); it.first != it.second; ++it.first) {
  //   cout << g.indices()[*(it.first)] << ": ";
  //   for (auto out = boost::out_edges(*(it.first), g()); out.first != out.second; ++out.first) {
  //     cout << "(" << g.indices()[source(*(out.first), g())]
  //          << "," << g.indices()[target(*(out.first), g())] << ") ";
  //   }
  //   cout << endl;
  // }

  // Initialize and run each particle along the graph for *DURATION*.
  std::vector<simsys::Particle> objects;
  for (int i = 0; i < NUM_OBJECT; ++i) {
    objects.push_back(simsys::Particle(g, i));
    objects[i].advance(g, DURATION);
    // objects[i].print(g);
  }

  // Generate RFID readings for each object.
  std::vector<std::vector<int> > readings = detect(g, objects);

  // for (size_t i = 0; i < readings.size(); ++i) {
  //   for (size_t j = 0; j < readings[i].size(); ++j)
  //     std::cout << j << " (" << objects[i].pos(g, j) << ") " << readings[i][j] << std::endl;
  //   std::cout << std::endl;
  // }

  range_query_windowsize(g, objects, readings, WINDOW_RATIOS);

  std::vector<std::pair<double, double> >
      precision = measure(Precision),
      recall = measure(Recall),
      f1score = measure(F1Score);

  {
    std::ofstream fout("precision.txt");
    char msg[64];
    for (size_t i = 0; i < precision.size(); ++i) {
      sprintf(msg, "%.02f %.6f %.6f",
              WINDOW_RATIOS[i], precision[i].first, precision[i].second);
      fout << msg << endl;
    }
    fout.close();
  }

  {
    std::ofstream fout("recall.txt");
    char msg[64];
    for (size_t i = 0; i < recall.size(); ++i) {
      sprintf(msg, "%.02f %.6f %.6f",
              WINDOW_RATIOS[i], recall[i].first, recall[i].second);
      fout << msg << endl;
    }
    fout.close();
  }

  {
    std::ofstream fout("f1score.txt");
    char msg[64];
    for (size_t i = 0; i < f1score.size(); ++i) {
      sprintf(msg, "%.02f %.6f %.6f",
              WINDOW_RATIOS[i], f1score[i].first, f1score[i].second);
      fout << msg << endl;
    }
    fout.close();
  }

  return 0;
}
