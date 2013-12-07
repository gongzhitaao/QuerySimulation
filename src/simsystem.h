#ifndef SRC_SIMSYSTEM_H_
#define SRC_SIMSYSTEM_H_

#pragma once

#include <fstream>

#include <vector>
#include <map>
#include <string>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

class SimSystem
{
public:
  enum sim_vertex_color_enum { ROOM, DOOR, HALL, VERTEX_COLOR_ENUM };

private:
  struct Vec2 { double x, y; };

  struct Vec3 {
    double x, y, t;
    Vec3(const Vec2 &v2, double tt = 0) : x(v2.x), y(v2.y), t(tt) {}
  };
  typedef std::vector<Vec3> trace;

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

  // void add_reader(double x, double y) {
  //   readers_.push_back(std::make_pair(x, y));
  // }

  void run(double duration, int num_object, int reader_id = -1);

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
  Vec2 get_coord(const Vec2 &pa, const Vec2 &pb, double ratio) {
    return {pa.x + (pb.x - pa.x) * ratio,
	pa.y + (pb.y - pa.y) * ratio};
  }

  ColorMap vertex_color_;
  CoordMap vertex_coord_;
  WeightMap edge_weight_;

  std::map<int, Vertex> label2vertex_;

  // the walking graph
  UndirectedGraph g_;

  std::vector<trace> record_;
};

#endif  // SRC_SIMSYSTEM_H_
