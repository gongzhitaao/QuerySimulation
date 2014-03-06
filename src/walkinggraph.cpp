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
#include "param.h"
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
    : vnames_(boost::get(boost::vertex_name, g_))
    , coords_(boost::get(vertex_coord_t(), g_))
    , colors_(boost::get(boost::vertex_color, g_))
    , enames_(boost::get(boost::edge_name, g_))
    , weights_(boost::get(boost::edge_weight, g_))
    , fg_(make_filtered_graph(g_,
        positive_index<edge_name_t>(
            boost::get(boost::edge_name, g_)),
        positive_index<vertex_name_t>(
            boost::get(boost::vertex_name, g_))))
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
      edges_[id] = boost::add_edge(
          vertices_.at(n1), vertices_.at(n2),
          {id++, {std::sqrt(CGAL::squared_distance(coord(n1),
                                                   coord(n2)))}},
          g_).first;
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

class InsertAnchor : public boost::default_bfs_visitor
{
 public:
  InsertAnchor(UndirectedGraph &g, anchor_map_t &anchors, double unit)
      : vnames_(boost::get(boost::vertex_name, g))
      , enames_(boost::get(boost::edge_name, g))
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
           w = boost::get(weights_, e);

    int id = boost::get(enames_, e);
    anchors_[id].push_back(
        std::make_pair(boost::get(vnames_, source), 0.0));
    while (d < w) {
      anchors_[id].push_back(std::make_pair(name_--, d / w));
      d += unit_;
    }
    anchors_[id].push_back(
        std::make_pair(boost::get(vnames_, target), 1.0));
    dist_[target] = d - w;
  }

 private:
  boost::unordered_map<Vertex, double> dist_;
  vertex_name_t vnames_;
  edge_name_t enames_;
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

  int id = -1;

  for (auto it = anchors_.begin(); it != anchors_.end(); ++it) {
    const std::vector<std::pair<int, double> > &vec = it->second;
    if (vec.size() <= 2) continue;

    Vertex source = vertices_[vec[0].first],
           target = vertices_[vec.back().first],
                u = source;

    double w = weight(vec[0].first, vec.back().first);
    int size = vec.size();

    for (int i = 1; i < size - 1; ++i) {
      Point_2 p = linear_interpolate(boost::get(coords_, source),
                                     boost::get(coords_, target),
                                     vec[i].second);
      Vertex v = boost::add_vertex({vec[i].first, {p, {HALL}}}, g_);

      vertices_[vec[i].first] = v;

      Edge e = boost::add_edge(
          u, v, {id, {w * (vec[i].second - vec[i-1].second)}},
          g_).first;

      edges_[id] = e;

      --id;
      u = v;
    }

    Edge e = boost::add_edge(u, vertices_[vec.back().first],
                             {id, {w - w * vec[size-2].second}},
                             g_).first;
    edges_[id] = e;
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
  Vertex pre = from < 0 ? fg_.null_vertex() : vertices_.at(from);

  std::set<Vertex> hall, door, room;

  auto pairit = boost::out_edges(cur, fg_);
  for (auto it = pairit.first; it != pairit.second; ++it) {
    Vertex v = boost::target(*it, fg_);
    switch (boost::get(colors_, v)) {
      case HALL: hall.insert(v); break;
      case DOOR: door.insert(v); break;
      case ROOM: room.insert(v); break;
      default: break;
    }
  }

  if (hall.size() == 0 || fg_.null_vertex() == pre) {
    if (hall.size() > 0) {
      boost::random::uniform_int_distribution<>
          unifi(0, hall.size() - 1);
      return boost::get(vnames_,
                        *boost::next(hall.begin(), unifi(gen)));
    }
    return boost::get(vnames_, *door.begin());
  }

  boost::random::uniform_real_distribution<> unifd(0, 1);

  if (room.size() > 0 && boost::get(colors_, pre) != ROOM &&
      unifd(gen) < ENTER_ROOM)
    return boost::get(vnames_, *room.begin());

  if (door.size() > 0 && (boost::get(colors_, cur) == ROOM ||
                          unifd(gen) < KNOCK_DOOR))
    return boost::get(vnames_, *door.begin());

  if (hall.size() > 1) hall.erase(pre);

  boost::random::uniform_int_distribution<> unifi(0, hall.size() - 1);
  return boost::get(vnames_, *boost::next(hall.begin(), unifi(gen)));
}

landmark_t
WalkingGraph::random_pos() const
{
  Edge e = boost::random_edge(fg_, gen);
  while (!boost::edge(boost::source(e, fg_),
                   boost::target(e, fg_), fg_).second)
    e = boost::random_edge(fg_, gen);

  boost::random::uniform_real_distribution<> unifd(0, 1);
  return boost::make_tuple(vnames_[boost::source(e, g_)],
                           vnames_[boost::target(e, g_)],
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

      int u = boost::get(vnames_, boost::source(e, g_));
      if (u == pos.get<1>()) {
        int t = pos.get<0>();
        pos.get<0>() = pos.get<1>();
        pos.get<1>() = t;
        pos.get<2>() = 1 - pos.get<2>();
      }

      int ind = boost::get(enames_, e);

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
      double w = boost::get(weights_, e);

      sort(vec.begin(), vec.end(),
           [](const std::pair<int, double> &a,
              const std::pair<int, double> &b) {
             return a.second < b.second;
           });

      Vertex u = source;
      double pre = 0.0;
      for (auto vit = vec.begin(); vit != vec.end(); ++vit) {
        Point_2 p = linear_interpolate(
            boost::get(coords_, source), boost::get(coords_, target),
            vit->second);
        Vertex v = boost::add_vertex({vit->first, {p, {HALL}}}, g_);

        vertices_[vit->first] = v;

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
  edge_with_objects_.clear();
}

struct found_knn {};

template <typename VertexName>
struct knn_visitor : public boost::default_bfs_visitor
{
  // knn_visitor() { }
  knn_visitor(std::vector<int> &results, int k, VertexName vnames)
      : results_(results)
      , count_(k)
      , vnames_(vnames){ }

  template <typename Graph>
  void
  discover_vertex(const Vertex v, const Graph &g)
  {
    int name = boost::get(vnames_, v);
    if (name >= 1000) {
      --count_;

      results_.push_back(name - 1000);

      if (0 == count_)
        throw found_knn();
    }
  }

  std::vector<int> &results_;
  int count_;
  VertexName vnames_;
};

std::vector<int>
WalkingGraph::nearest_neighbors(int object, int k)
{
  std::vector<int> results;
  knn_visitor<vertex_name_t> vis(results, k, vnames_);

  try {
    boost::breadth_first_search(fg_, vertices_[OBJECT_START + object],
                                boost::visitor(vis));
  } catch (found_knn fk) { }

  return results;
}

void
WalkingGraph::print(std::ostream &os) const
{
  // boost::graph_traits<UndirectedGraph>::vertex_iterator vi, end;
  // for (boost::tie(vi, end) = boost::vertices(g_); vi != end; ++vi)
  //   os << boost::get(coords_, *vi) << std::endl;
  boost::graph_traits<UndirectedGraph>::edge_iterator ei, eend;
  for (boost::tie(ei, eend) = boost::edges(g_); ei != eend; ++ei) {
    Vertex v = boost::source(*ei, g_);
    Vertex u = boost::target(*ei, g_);
    os << boost::get(coords_, v) << ' '
       << boost::get(coords_, u) << endl;
  }
}

template <typename VertexName>
struct nearest_anchor : public boost::default_bfs_visitor
{
  nearest_anchor(int &result, VertexName vnames)
      : result_(result)
      , vnames_(vnames){ }

  template <typename Graph>
  void
  discover_vertex(const Vertex v, const Graph &g)
  {
    int name = boost::get(vnames_, v);
    if (name < 0) {
      result_ = boost::get(vnames_, v);
      throw found_knn();
    }
  }

  int &result_;
  VertexName vnames_;
};

int
WalkingGraph::align(const landmark_t &p)
{
  int result = 0;
  nearest_anchor<vertex_name_t> vis(result, vnames_);

  Edge e = boost::edge(vertices_.at(p.get<0>()),
                       vertices_.at(p.get<1>()), g_).first;

  std::vector<std::pair<int, double> > &vec = anchors_[enames_[e]];

  if (vec.empty()) {
    vec.push_back(std::make_pair(
        boost::get(vnames_, boost::source(e, g_)), 0.0));
    vec.push_back(std::make_pair(
        boost::get(vnames_, boost::target(e, g_)), 1.0));
  }

  auto pos = std::lower_bound(
      vec.begin(), vec.end(), p.get<2>(),
      [](const std::pair<int, double> &p, const double d) {
        return p.second < d;
      });

  assert(pos != vec.end());

  if (vec.begin() == pos)
    std::advance(pos, 1);

  assert(pos != vec.begin());

  int id = 5000;

  pos = vec.insert(pos, std::make_pair(id, p.get<2>()));

  Vertex curr = boost::add_vertex(id, g_);

  auto p0 = std::prev(pos),
       p1 = std::next(pos);

  double w = boost::get(weights_, e);

  // cout << p0->first << ' ' << p1->first << endl;

  Vertex prev = vertices_.at(p0->first);
  Vertex next = vertices_.at(p1->first);

  boost::add_edge(prev, curr,
                  {id, {(p.get<2>() - p0->second) * w}}, g_);
  boost::add_edge(curr, next,
                  {id + 1, {(p1->second - p.get<2>()) * w}}, g_);

  Edge e_old = boost::edge(prev, next, g_).first;
  int id_old = boost::get(enames_, e_old);
  double w_old = boost::get(weights_, e_old);
  boost::remove_edge(prev, next, g_);

  try {
    boost::breadth_first_search(g_, curr, boost::visitor(vis));
  } catch (found_knn fk) { }

  boost::add_edge(prev, next, {id_old, {w_old}}, g_);
  boost::remove_edge(prev, curr, g_);
  boost::remove_edge(curr, next, g_);
  vec.erase(pos);

  return result;
}

}
