#include <cmath>
#include <fstream>

#include <boost/property_map/property_map.hpp>
#include <boost/graph/breadth_first_search.hpp>

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#include <boost/iterator/zip_iterator.hpp>

#include <boost/algorithm/string.hpp>

#include "walkinggraph.h"
#include "param.h"

namespace simulation {

using std::cout;
using std::endl;

WalkingGraph::WalkingGraph()
    : names_(boost::get(boost::vertex_name, g_))
    , coords_(boost::get(vertex_coord_t(), g_))
    , colors_(boost::get(boost::vertex_color, g_))
    , indices_(boost::get(boost::edge_index, g_))
    , weights_(boost::get(boost::edge_weight, g_))
    , fg_(make_filtered_graph(
        g_, boost::keep_all(),
        remove_anchor(boost::get(boost::vertex_name, g_))))
{
  initialize();
  insert_anchors();
}

void
WalkingGraph::initialize()
{
  const std::string Commenter("//");
  const std::string DataRoot("/home/gongzhitaao/Documents/simsystem/data/jiao/");

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
      vertices_[name] = boost::add_vertex(
          {name, {Point_2(x, y), {(vertex_color_enum) type}}}, g_);
    }
    fin.close();
  }

  // Read in edge for walking graph
  {
    std::ifstream fin(DataRoot + "edge.txt");
    std::string line;
    int id = 0;
    while (std::getline(fin, line)) {
      boost::algorithm::trim_left(line);
      if (line.compare(0, Commenter.size(), Commenter) == 0) continue;

      std::stringstream ss(line);
      int n1, n2;
      ss >> n1 >> n2;
      auto pair = boost::add_edge(
          vertices_[n1], vertices_[n2],
          {id, {std::sqrt(CGAL::squared_distance(
              coord(n1), coord(n2)))}}, g_);
      if (pair.second)
        edges_[id] = pair.first;
      id++;
    }
    fin.close();
  }

  // Read in RFID readers.
  {
    std::vector<std::pair<int, landmark_t> > infos;
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
      infos.push_back({id, {v1, v2, p}});
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

Point_2
linear_interpolate(const Point_2 &p0, const Point_2 &p1, double a)
{
  return p0 + (p1 - p0) * a;
}

class InsertAnchor : public boost::default_bfs_visitor
{
 public:
  InsertAnchor(UndirectedGraph &g, anchor_map_t &anchors, double unit)
      : names_(boost::get(boost::vertex_name, g))
      , indices_(boost::get(boost::edge_index, g))
      , weights_(boost::get(boost::edge_weight, g))
      , anchors_(anchors)
      , unit_(unit)
      , name_(-1) { }

  void
  initialize_vertex(Vertex v, const UndirectedGraph &g)
  { dist_[v] = 0.0; }

  void
  tree_edge(Edge e, const UndirectedGraph &g)
  {
    Vertex source = boost::source(e, g),
           target = boost::target(e, g);

    double d = dist_[source],
           w = weights_[e];

    int id = indices_[e];
    anchors_[id].push_back(std::make_pair(names_[source], 0.0));
    while (d <= w) {
      anchors_[id].push_back(std::make_pair(name_--, d / w));
      d += unit_;
    }
    anchors_[id].push_back(std::make_pair(names_[target], 1.0));
    dist_[target] = d - w;
  }

 private:
  boost::unordered_map<Vertex, double> dist_;
  name_map_t names_;
  index_map_t indices_;
  weight_map_t weights_;
  anchor_map_t &anchors_;
  double unit_;
  int name_;
};

void
WalkingGraph::insert_anchors(double unit)
{
  // calculate the anchor points in a bfs fashion
  InsertAnchor vis(g_, anchors_, 20);
  boost::breadth_first_search(g_, boost::random_vertex(g_, gen),
                              boost::visitor(vis));

  int id = boost::num_edges(g_);

  for (auto it = anchors_.begin(); it != anchors_.end(); ++it) {
    const std::vector<std::pair<int, double> > &vec = it->second;
    if (vec.size() <= 2) continue;

    Vertex source = vertices_[vec[0].first],
           target = vertices_[vec.back().first],
                u = source;

    double w = weight(vec[0].first, vec.back().first);

    for (size_t i = 1; i < vec.size(); ++i) {
      Point_2 p = linear_interpolate(coords_[source], coords_[target],
                                     vec[i].second);
      Vertex v = boost::add_vertex({vec[i].first, {p, {HALL}}}, g_);
      auto pair = boost::add_edge(
          u, v, {id, {w * (vec[i].second - vec[i-1].second)}}, g_);
      if (pair.second)
        edges_[id] = pair.first;
      ++id;
      u = v;
    }
  }
}

// Randomly choose next vertex to advance to.  If u which is where the
// particle came from is present, then the randomly chosen vertex may
// not be u unless the out degree of v is only one in which we have no
// choice.  In this way, the particle preserves its direction.
int
WalkingGraph::random_next(int to, int from) const
{
  Vertex cur = vertices_.at(to);
  Vertex pre = vertices_.at(from);

  std::set<Vertex> hall, door, room;
  auto pairit = boost::out_edges(cur, fg_);
  for (auto it = pairit.first; it != pairit.second; ++it) {
    Vertex v = boost::target(*it, fg_);
    switch (colors_(v)) {
      case HALL: hall.insert(v); break;
      case DOOR: door.insert(v); break;
      case ROOM: room.insert(v); break;
      default: break;
    }
  }

  if (fg_.null_vertex() == pre) {
    if (hall.size() > 0) {
      boost::random::uniform_int_distribution<>
          unifi(0, hall.size() - 1);
      return names_[*(boost::next(hall.begin(), unifi(gen)))];
    }
    return names_[*(door.begin())];
  }

  boost::random::uniform_real_distribution<> unifd(0, 1);

  if (room.size() > 0 && colors_[pre] != ROOM &&
      unifd(gen) < ENTER_ROOM)
    return names_[*(room.begin())];

  if (door.size() > 0
      && (colors_[cur] == ROOM ||
          unifd(gen) < KNOCK_DOOR))
      return names_[*(door.begin())];

  if (hall.size() > 1) hall.erase(pre);

  boost::random::uniform_int_distribution<> unifi(0, hall.size() - 1);
  return names_[*(boost::next(hall.begin(), unifi(gen)))];
}

landmark_t
WalkingGraph::random_pos() const
{
  Edge e = boost::random_edge(fg_, gen);
  boost::random::uniform_real_distribution<> unifd(0, 1);
  return boost::make_tuple(names_[boost::source(e, fg_)],
                           names_[boost::target(e, fg_)],
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

  std::vector<Vertex_handle> res;
  readerset_.range_search(circle, std::back_inserter(res));

  if (res.empty()) return -1;

  return (*res.begin())->info().first;
}

static int count_;
static std::vector<int> edge_with_objects_;

void
WalkingGraph::insert_objects(const std::vector<landmark_t> &objects)
{
  {
    int id = OBJECT_START;
    count_ = objects.size();
    for (size_t i = 0; i < objects.size(); ++i) {
      landmark_t pos = objects[i];
      Edge e = boost::edge(vertices_[pos.get<0>()],
                           vertices_[pos.get<1>()], fg_).first;

      int u = names_[boost::source(e, g_)];
      if (u == pos.get<1>()) {
        int t = pos.get<0>();
        pos.get<0>() = pos.get<1>();
        pos.get<1>() = t;
        pos.get<2>() = 1 - pos.get<2>();
      }

      int ind = indices_[e];

      if (objects_[ind].empty())
        edge_with_objects_.push_back(ind);

      objects_[ind].push_back(std::make_pair(id + i, pos.get<2>()));
    }
  }

  {
    int id = boost::num_edges(g_);
    for (auto it = edge_with_objects_.begin();
         it != edge_with_objects_.end(); ++it) {
      std::vector<std::pair<int, double> > &vec = objects_[*it];
      Edge e = edges_[*it];
      Vertex source = boost::source(e, g_);
      Vertex target = boost::target(e, g_);
      double w = weights_[e];

      sort(vec.begin(), vec.end(),
           [](const std::pair<int, double> &a,
              const std::pair<int, double> &b) {
             return a.second < b.second;
           });

      Vertex u = source;
      double pre = 0.0;
      for (auto vit = vec.begin(); vit != vec.end(); ++vit) {
        Point_2 p = linear_interpolate(
            coords_[source], coords_[target], vit->second);
        Vertex v = boost::add_vertex({vit->first, {p, {HALL}}}, g_);
        auto pair = boost::add_edge(u, v,
                                 {id, {w * vit->second - pre}}, g_);
        if (pair.second)
          edges_[id] = pair.first;
        ++id;

        pre = w * vit->second;
      }
    }
  }
}

void
WalkingGraph::clear_objects()
{
  for (int i = OBJECT_START; i < count_; ++i)
    boost::remove_vertex(vertices_[i], g_);

  for (auto it = edge_with_objects_.begin();
       it != edge_with_objects_.end(); ++it) {
    objects_[*it].clear();
  }
}

struct found_knn {};

struct knn_visitor : public boost::default_bfs_visitor
{
  // knn_visitor() { }
  knn_visitor(std::vector<int> &results, int k,
              name_map_t names)
      : results_(results)
      , count_(k)
      , names_(names){ }

  template <typename Graph>
  void
  discover_vertex(const Vertex v, const Graph &g)
  {
    int name = names_[v];
    if (name >= 1000) {
      --count_;

      results_.push_back(name);

      if (0 == count_)
        throw found_knn();
    }
  }

  std::vector<int> &results_;
  int count_;
  name_map_t names_;
};

std::vector<int>
WalkingGraph::nearest_neighbor(int object, int k)
{
  std::vector<int> results;
  knn_visitor vis(results, k, names_);
  try {
    boost::breadth_first_search(fg_, vertices_[OBJECT_START + object],
                                boost::visitor(vis));
  } catch (found_knn fk) { }

  return results;
}

void
WalkingGraph::print(std::ostream &os) const
{
  boost::graph_traits<UndirectedGraph>::vertex_iterator vi, end;
  for (boost::tie(vi, end) = boost::vertices(g_); vi != end; ++vi)
    os << coords_[*vi] << std::endl;
}

// int
// align(const Point_2 &p)
// {
//   K_neighbor_search search(anchortree_, p, 1);
//   Point_2 center = boost::get<0>(search.begin()->first);
//   return boost::get<1>(search.begin()->first);
// }

// std::vector<Point_and_int>
// WalkingGraph::anchors(const Fuzzy_iso_box &win)
// {
//   std::vector<Point_and_int> res;
//   anchortree_.search(std::back_inserter(res), win);
//   return res;
// }

// struct found_goal {};

// template <typename Vertex>
// class astar_goal_visitor : public boost::default_astar_visitor
// {
//  public:
//   astar_goal_visitor(Vertex goal) : goal_(goal) {}

//   template <typename Graph>
//   void examine_vertex(Vertex &u, const Graph &g) {
//     if (u == goal_) throw found_goal();
//   }
//  private:
//   Vertex goal_;
// };

// template<typename Graph, typename CoordMap>
// class heuristic : public boost::astar_heuristic<Graph, int>
// {
//  public:
//   typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;

//   heuristic(Vertex &goal, const CoordMap &coord)
//       :  goal_(goal), coord_(coord) {}

//   int operator()(Vertex u) {
//     return std::sqrt(CGAL::squared_distance(coord_[goal_], coord_[u]));
//   }

//  private:
//   Vertex goal_;
//   const CoordMap &coord_;
// };

// std::vector<Vertex>
// WalkingGraph::path(Vertex source, Vertex target) const
// {
//   std::vector<Vertex> p(boost::num_vertices(g_));
//   astar_goal_visitor<Vertex> visitor(source);
//   try {
//     boost::astar_search(
//         g_, target, heuristic<UndirectedGraph, CoordMap>(target, coords_),
//         boost::predecessor_map(boost::make_iterator_property_map(
//             p.begin(), boost::get(boost::vertex_index, g_))).
//         visitor(visitor));
//   } catch (found_goal fg) {
//     return p;
//   }

//   return std::vector<Vertex>();
// }

}
