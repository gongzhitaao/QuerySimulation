#ifndef SRC_SIMSYSTEM_H_
#define SRC_SIMSYSTEM_H_

#pragma once

#include <fstream>

#include <vector>
#include <map>
#include <string>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

namespace simsys {

struct Vec2 { double x, y; };

struct Vec3 {
  double x, y, t;
  Vec3(const Vec2 &v2, double tt = 0) : x(v2.x), y(v2.y), t(tt) {}
  operator const Vec2() const { return {x, y}; }
};

typedef std::vector<Vec3> trace;

class SimSystem
{
 public:
  enum sim_vertex_color_enum { ROOM, DOOR, HALL, VERTEX_COLOR_ENUM };

 private:
  // vertex property
  struct vertex_coord { typedef boost::vertex_property_tag kind; };
  typedef boost::property<vertex_coord, Vec2> CoordProperty;
  typedef boost::property<boost::vertex_color_t, sim_vertex_color_enum, CoordProperty> VertexProperty;

  // edge property
  typedef boost::property<boost::edge_weight_t, double> EdgeProperty;

  // graph type
  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                                VertexProperty, EdgeProperty> UndirectedGraph;

  // vertex type
  typedef boost::graph_traits<UndirectedGraph>::vertex_descriptor Vertex;

  // edge type
  typedef boost::graph_traits<UndirectedGraph>::edge_descriptor Edge;

  // get property by descriptor
  typedef boost::property_map<UndirectedGraph, boost::vertex_color_t>::type ColorMap;
  typedef boost::property_map<UndirectedGraph, vertex_coord>::type CoordMap;
  typedef boost::property_map<UndirectedGraph, boost::edge_weight_t>::type WeightMap;

 public:
  SimSystem() {
    vertex_color_ = boost::get(boost::vertex_color, g_);
    vertex_coord_ = boost::get(vertex_coord(), g_);
    edge_weight_ = boost::get(boost::edge_weight, g_);
  }

  void add_node(int label, double x, double y, sim_vertex_color_enum color = HALL) {
    Vertex u = boost::add_vertex(g_);
    vertex_color_[u] = color;
    vertex_coord_[u] = {x, y};
    label2vertex_[label] = u;
  }

  void add_edge(int src, int des) {
    Vertex u = label2vertex_[src];
    Vertex v = label2vertex_[des];
    Edge e = (boost::add_edge(u, v, g_)).first;

    Vec2 beg = vertex_coord_[u];
    Vec2 end = vertex_coord_[v];
    edge_weight_[e] = sqrt((beg.x - end.y) * (beg.x - end.x) +
                           (beg.y - end.y) * (beg.y - end.y));
  }

  void run(double duration, int num_object, int reader_id = -1);

  std::vector<Vec3> get_snapshot(double t) const {
    std::vector<Vec3> position;
    for (unsigned i = 0; i < record_.size(); ++i) {
      position.push_back(get_pos_by_timestamp(record_[i], t));
    }
    return position;
  }

  std::vector<trace> get_snapshot(double start, double end, double step = 0.1) const {
    std::vector<trace> result;
    for (unsigned i = 0; i < record_.size(); ++i) {
      trace tr;
      for (double t = start; t <= end; t += step)
        tr.push_back(get_pos_by_timestamp(record_[i], t));
      result.push_back(tr);
    }
    return result;
  }

  void load_trace(const std::string &fname);

  void save_trace(const std::string &fname) const {
    std::ofstream fout(fname);
    char msg[64];
    for (unsigned i = 0; i < record_.size(); ++i) {
      for (unsigned j = 0; j < record_[i].size(); ++j) {
        sprintf(msg, "%d %.6f %.6f %.6f",
                i, record_[i][j].t, record_[i][j].x, record_[i][j].y);
        fout << msg << std::endl;
      }
    }
    fout.close();
  }

 private:
  Vec3 get_pos_by_timestamp(const trace &tr, double t) const {
    int low = 0, high = tr.size() - 1, mid;
    if (t >= tr[high].t) return tr[high];

    while (low < high) {
      mid = (low + high) / 2;
      if (tr[mid].t > t) high = mid;
      else low = mid + 1;
    }

    double ratio = (t - tr[low - 1].t) / (tr[low].t - tr[low - 1].t);
    return Vec3(interpolate(tr[low - 1], tr[low], ratio), t);
  }

  Vec2 interpolate(const Vec2 &start, const Vec2 &end, double ratio) const {
    return {start.x + (end.x - start.x) * ratio,
          start.y + (end.y - start.y) * ratio};
  }

  ColorMap vertex_color_;
  CoordMap vertex_coord_;
  WeightMap edge_weight_;

  std::map<int, Vertex> label2vertex_;

  // the walking graph
  UndirectedGraph g_;

  std::vector<trace> record_;
};

}

#endif  // SRC_SIMSYSTEM_H_
