#include <cmath>
#include <fstream>

#include <boost/property_map/property_map.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/graph_utility.hpp>

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#include <boost/iterator/zip_iterator.hpp>

#include <boost/algorithm/string.hpp>

#include "walkinggraph.h"
#include "global.h"

namespace simulation {

using std::cout;
using std::endl;

Point_2
linear_interpolate(const Point_2 &p0, const Point_2 &p1, double a)
{
  return p0 + (p1 - p0) * a;
}

WalkingGraph::WalkingGraph()
    : vnames_(boost::get(boost::vertex_name, wg_))
    , enames_(boost::get(boost::edge_name, wg_))
{
  initialize();
  insert_anchors();
}

void
WalkingGraph::initialize()
{
  const std::string Commenter("//");
  const std::string DataRoot("/home/gongzhitaao/Documents/"
                             "simulator/data/jiao/");

  // Read in vertices for walking graph.
  {
    std::ifstream fin(DataRoot + "node.txt");
    std::string line;
    while (std::getline(fin, line)) {
      boost::algorithm::trim_left(line);
      if (line.compare(0, Commenter.size(), Commenter) == 0) continue;

      std::stringstream ss(line);
      int name, type;
      double x, y;
      ss >> name >> x >> y >> type;
      vertices_[name] = boost::add_vertex(name, wg_);
      vertices_[name] = boost::add_vertex(name, ag_);
      vp_[name] = {Point_2(x, y), (vertex_color_enum) type};
    }
    fin.close();
  }

  // Read in edge for walking graph
  {
    std::ifstream fin(DataRoot + "edge.txt");
    std::string line;
    for (int id = 0; std::getline(fin, line); ++id) {
      boost::algorithm::trim_left(line);
      if (line.compare(0, Commenter.size(), Commenter) == 0) continue;

      std::stringstream ss(line);
      int n1, n2;
      ss >> n1 >> n2;
      boost::add_edge(vertices_.at(n1), vertices_.at(n2), id, wg_);
      boost::add_edge(vertices_.at(n1), vertices_.at(n2), id, ag_);
      ep_[id].weight = std::sqrt(
          CGAL::squared_distance(coord(n1), coord(n2)));
    }
    fin.close();
  }

  // Read in RFID readers.
  {
    std::vector<int> infos;
    std::vector<Point_2> points;

    std::ifstream fin(DataRoot + "rfid2.txt");
    std::string line;
    while (std::getline(fin, line)) {
      boost::algorithm::trim_left(line);
      if (line.compare(0, Commenter.size(), Commenter) == 0) continue;

      std::stringstream ss(line);
      int id, v1, v2;
      double x, y;
      ss >> id >> x >> y >> v1 >> v2;

      Point_2 p0(x, y), p1 = coord(v1), p2 = coord(v2);
      double p = std::sqrt(CGAL::squared_distance(p0, p1) /
                           CGAL::squared_distance(p1, p2));
      infos.push_back(id);
      points.push_back(Point_2(x, y));
      readermap_[id] = {v1, v2, p};
    }
    fin.close();

    readerset_.insert(
        boost::make_zip_iterator(boost::make_tuple(points.begin(),
                                                   infos.begin())),
        boost::make_zip_iterator(boost::make_tuple(points.end(),
                                                   infos.end())));
  }

  // Read in rooms configuration.
  {
    std::ifstream fin(DataRoot + "room.txt");
    std::string line;
    while (std::getline(fin, line)) {
      boost::algorithm::trim_left(line);
      if (line.compare(0, Commenter.size(), Commenter) == 0) continue;

      std::stringstream ss(line);
      double x0, y0, x1, y1;
      int id;
      ss >> id >> x0 >> y0 >> x1 >> y1;
      rooms_.push_back(IsoRect_2(x0, y0, x1, y1));

      if (x1 > xmax_) xmax_ = x1;
      if (y1 > ymax_) ymax_ = y1;
    }
    fin.close();
  }

  // Read in halls configuration
  {
    std::ifstream fin(DataRoot + "hall.txt");
    std::string line;
    while (std::getline(fin, line)) {
      boost::algorithm::trim_left(line);
      if (line.compare(0, Commenter.size(), Commenter) == 0) continue;
      std::stringstream ss(line);
      double x0, y0, x1, y1;
      int dir;
      ss >> x0 >> y0 >> x1 >> y1 >> dir;
      halls_.push_back(IsoRect_2(x0, y0, x1, y1));
      dirs_.push_back(dir);
    }
    fin.close();
  }
}

class InsertAnchor : public boost::default_bfs_visitor
{
 public:
  InsertAnchor(WalkingGraph &g, double unit = 20.0)
      : g_(g), name_(g_.ANCHORID), unit_(unit) { }

  void
  initialize_vertex(Vertex v, const UndirectedGraph &)
  { dist_[v] = 0.0; }

  void
  tree_edge(Edge e, const UndirectedGraph &g)
  {
    Vertex source = boost::source(e, g),
           target = boost::target(e, g);
    int id = boost::get(g_.enames_, e);

    double d = dist_[source],
           w = g_.ep_.at(id).weight;

    g_.ep_.at(id).anchors.push_back(
        std::make_pair(boost::get(g_.vnames_, source), 0.0));
    while (d < w) {
      g_.ep_.at(id).anchors.push_back(std::make_pair(name_++, d / w));
      d += unit_;
    }
    g_.ep_.at(id).anchors.push_back(
        std::make_pair(boost::get(g_.vnames_, target), 1.0));
    dist_[target] = d - w;
  }

 private:
  boost::unordered_map<Vertex, double> dist_;
  WalkingGraph &g_;
  int name_;
  double unit_;
};

void
WalkingGraph::insert_anchors(double unit)
{
  // calculate the anchor points in a bfs fashion
  InsertAnchor vis(*this, 20);
  boost::breadth_first_search(wg_, boost::random_vertex(wg_, gen),
                              boost::visitor(vis));

  int id = boost::num_edges(ag_);

  std::vector<int> infos;
  std::vector<Point_2> points;

  for (auto it = boost::edges(wg_); it.first != it.second;
       ++it.first, ++id) {

    const std::vector<std::pair<int, double> > &vec =
        ep_.at(boost::get(enames_, *it.first)).anchors;
    if (vec.size() <= 2) continue;

    Vertex source = vertices_[vec[0].first],
           target = vertices_[vec.back().first],
                u = source;

    double w = weight(vec[0].first, vec.back().first);
    int size = vec.size();

    boost::remove_edge(source, target, ag_);

    for (int i = 1; i < size - 1; ++i) {
      int from = boost::get(vnames_, source),
            to = boost::get(vnames_, target);
      Point_2 p = linear_interpolate(coord(from),
                                     coord(to), vec[i].second);

      Vertex v = boost::add_vertex(vec[i].first, ag_);
      vertex_color_enum c = HALL;
      if (color(from) == ROOM || color(to) == ROOM) c = ROOM;
      vertices_[vec[i].first] = v;

      vp_[vec[i].first] = {p, c};
      ap_[vec[i].first] = {from, to, vec[i].second};

      points.push_back(p);
      infos.push_back(vec[i].first);

      boost::add_edge(u, v, id, ag_);
      ep_[id].weight = w * (vec[i].second - vec[i-1].second);

      u = v;
    }

    boost::add_edge(u, vertices_[vec.back().first], id, ag_);
    ep_[id].weight = w * (1 - vec[size-2].second);
  }

  for (auto it = boost::vertices(wg_); it.first != it.second;
       ++it.first) {
    int name = boost::get(vnames_, *it.first);
    infos.push_back(name);
    points.push_back(coord(name));
  }

  anchorset_.insert(
      boost::make_zip_iterator(boost::make_tuple(points.begin(),
                                                 infos.begin())),
      boost::make_zip_iterator(boost::make_tuple(points.end(),
                                                 infos.end())));
}

// Randomly choose next vertex to advance to.  If u which is where the
// particle came from is present, then the randomly chosen vertex may
// not be u unless the out degree of v is only one in which we have no
// choice.  In this way, the particle preserves its direction.
int
WalkingGraph::random_next(int to, int from) const
{
  Vertex cur = vertices_.at(to);
  Vertex pre = from < 0 ? wg_.null_vertex() : vertices_.at(from);

  boost::unordered_set<Vertex> hall, door, room;

  auto pit = boost::out_edges(cur, wg_);
  for (auto it = pit.first; it != pit.second; ++it) {
    Vertex v = boost::target(*it, wg_);
    switch (color(boost::get(vnames_, v))) {
      case HALL: hall.insert(v); break;
      case DOOR: door.insert(v); break;
      case ROOM: room.insert(v); break;
      default: break;
    }
  }

  if (hall.size() == 0 || wg_.null_vertex() == pre) {
    if (hall.size() > 0) {
      boost::random::uniform_int_distribution<>
          unifi(0, hall.size() - 1);
      return boost::get(vnames_,
                        *boost::next(hall.begin(), unifi(gen)));
    }
    return boost::get(vnames_, *door.begin());
  }

  boost::random::uniform_real_distribution<> unifd(0, 1);

  if (room.size() > 0 && color(boost::get(vnames_, pre)) != ROOM &&
      unifd(gen) < enter_room_)
    return boost::get(vnames_, *room.begin());

  if (door.size() > 0 && (color(boost::get(vnames_, cur)) == ROOM ||
                          unifd(gen) < knock_door_))
    return boost::get(vnames_, *door.begin());

  if (hall.size() > 1) hall.erase(pre);

  boost::random::uniform_int_distribution<> unifi(0, hall.size() - 1);
  return boost::get(vnames_, *boost::next(hall.begin(), unifi(gen)));
}

landmark_t
WalkingGraph::random_pos() const
{
  Edge e = boost::random_edge(wg_, gen);
  boost::random::uniform_real_distribution<> unifd(0, 1);
  return boost::make_tuple(vnames_[boost::source(e, wg_)],
                           vnames_[boost::target(e, wg_)],
                           unifd(gen));
}

std::vector<std::pair<IsoRect_2, double> >
WalkingGraph::random_window(double ratio) const
{
  double r = 1 - std::sqrt(ratio);
  boost::random::uniform_real_distribution<>
      unifx(0, xmax_ * r), unify(0, ymax_ * r);

  double xmin = unifx(gen), ymin = unify(gen);

  IsoRect_2 win = IsoRect_2(xmin, ymin,
                            xmin + xmax_ * (1 - r),
                            ymin + ymax_ * (1 - r));

  std::vector<std::pair<IsoRect_2, double> > results;

  // Intersection with rooms
  for (size_t i = 0; i < rooms_.size(); ++i) {
    auto res = CGAL::intersection(win, rooms_[i]);
    // if the query intersects with a room, then the intersected part
    // extends to the whole room.
    if (res) {
      const IsoRect_2 tmp = *boost::get<IsoRect_2>(&*res);
      const IsoRect_2 room = rooms_[i];
      results.push_back(std::make_pair(tmp, tmp.area() / room.area()));
    }
  }

  // Intersection with halls.
  for (size_t i = 0; i < halls_.size(); ++i) {
    auto res = CGAL::intersection(win, halls_[i]);
    if (res) {
      const IsoRect_2 tmp = *boost::get<IsoRect_2>(&*res);
      const IsoRect_2 hall = halls_[i];
      if (0 == dirs_[i])
        results.push_back(std::make_pair(
            IsoRect_2(Point_2(tmp.xmin(), hall.ymin()),
                      Point_2(tmp.xmax(), hall.ymax())),
            (tmp.ymax() - tmp.ymin()) / (hall.ymax() - hall.ymax())));
      else
        results.push_back(std::make_pair(
            IsoRect_2(Point_2(hall.xmin(), tmp.ymin()),
                      Point_2(hall.xmax(), tmp.ymax())),
            (tmp.xmax() - tmp.xmin()) / (hall.xmax() - hall.xmin())));
    }
  }

  return results;
}

int
WalkingGraph::detected_by(const landmark_t &pos, double radius)
{
  Point_2 center = linear_interpolate(
      coord(pos.get<0>()), coord(pos.get<1>()), pos.get<2>());
  Circle_2 circle(center, radius);

  std::vector<vertex_handle> res;
  readerset_.range_search(circle, std::back_inserter(res));

  if (res.empty()) return -1;

  return (*res.begin())->info();
}

void
WalkingGraph::print(std::ostream &os) const
{
  // boost::graph_traits<UndirectedGraph>::vertex_iterator vi, end;
  // for (boost::tie(vi, end) = boost::vertices(wg_); vi != end; ++vi)
  //   os << boost::get(coords_, *vi) << std::endl;
  boost::graph_traits<UndirectedGraph>::edge_iterator ei, eend;
  for (boost::tie(ei, eend) = boost::edges(wg_); ei != eend; ++ei) {
    Vertex v = boost::source(*ei, wg_);
    Vertex u = boost::target(*ei, wg_);
    os << coord(boost::get(vnames_, v)) << ' '
       << coord(boost::get(vnames_, u)) << endl;
  }
}

int
WalkingGraph::align(const landmark_t &p)
{
  Point_2 pt = coord(p);
  return anchorset_.nearest_neighbor(pt)->info();
}

std::vector<int>
WalkingGraph::anchors_in_win(const IsoRect_2 &w)
{
  std::vector<vertex_handle> tmp;
  anchorset_.range_search(w.vertex(0), w.vertex(1),
                          w.vertex(2), w.vertex(3),
                          std::back_inserter(tmp));
  std::vector<int> results;
  std::transform(tmp.begin(), tmp.end(), std::back_inserter(results),
                 [] (const vertex_handle vh) {
                   return vh->info();
                 });
  return results;
}

}
