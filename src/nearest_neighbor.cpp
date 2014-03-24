#include <algorithm>
#include <iterator>
#include <utility>

#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>

#include "nearest_neighbor.h"
#include "global.h"

namespace simulation {

using std::cout;
using std::endl;

template <typename W, typename ID>
class weight_map
{
 public:
  typedef typename boost::property_traits<ID>::key_type key_type;
  typedef typename W::mapped_type value_type;
  typedef typename W::mapped_type &reference;
  typedef boost::lvalue_property_map_tag category;

  weight_map()
      : id_(0), w_(0){ }

  weight_map(const ID *id, W *w = 0)
      : id_(id), w_(w) { }

  void
  clear()
  {
    w_->clear();
  }

  const ID *id_;
  W *w_;
};

template <typename W, typename ID>
weight_map<W, ID>
make_my_weight_map(const ID *id, W *w = 0)
{
  return weight_map<W, ID>(id, w);
}

template <typename W, typename ID>
typename W::mapped_type
get(const weight_map<W, ID> &w,
    typename boost::property_traits<ID>::key_type key)
{
  return w.w_->at((*w.id_)[key]);
}

template <typename W, typename ID>
void
put(weight_map<W, ID> &w,
    typename boost::property_traits<ID>::key_type key,
    const typename W::mapped_type &value)
{
  (*w.w_)[(*w.id_)[key]] = value;
}

template <typename W, typename ID>
typename W::mapped_type &
at(const weight_map<W, ID> &w,
   typename boost::property_traits<ID>::key_type key)
{
  return w.w_->at((*w.id_)[key]);
}

typedef weight_map<boost::unordered_map<int, double>,
                   edge_name_t> WeightMap;

static int object_;
static boost::unordered_map<
  int, std::vector<std::pair<int, double> > > landmarks_;

static UndirectedGraph g_copy_;
static UndirectedGraph a_copy_;
static edge_name_t enames_;
static boost::unordered_map<int, Edge> edges_;
static vertex_name_t vnames_;
static boost::unordered_map<int, Vertex> vertices_;
static WeightMap weights_;
static int idbase_;

NearestNeighbor::NearestNeighbor(Simulator &sim)
    : sim_(sim)
{
  boost::copy_graph(sim_.g_.ag_, a_copy_);
  // for (auto i = boost::edges(a_copy_); i.first != i.second;
  //      ++i.first)
}

void
NearestNeighbor::copy_graph()
{
  g_copy_.clear();

  edges_.clear();
  vertices_.clear();

  boost::copy_graph(sim_.g_.wg_, g_copy_);

  enames_ = boost::get(boost::edge_name, g_copy_);
  for (auto i = boost::edges(g_copy_); i.first != i.second;
       ++i.first)
    edges_[enames_[*(i.first)]] = *(i.first);

  vnames_ = boost::get(boost::vertex_name, g_copy_);
  for (auto i = boost::vertices(g_copy_); i.first != i.second;
       ++i.first)
    vertices_[vnames_[*(i.first)]] = *(i.first);

  weights_ = make_my_weight_map<boost::unordered_map<int, double>,
                                edge_name_t>(&enames_);

  for (auto i = sim_.g_.ep_.begin(); i != sim_.g_.ep_.end(); ++i) {
    cout << i->first << endl;
    put(weights_, edges_.at(i->first), i->second.weight);
  }

  cout << "helo" << endl;
}

void
NearestNeighbor::random_object()
{
  boost::random::uniform_int_distribution<> unifd(0, sim_.num_object_);
  object_ = unifd(gen);
}

void
NearestNeighbor::prepare(double t)
{
  copy_graph();

  std::vector<landmark_t> points = sim_.positions(t);

  for (size_t i = 0; i < points.size(); ++i) {
    landmark_t &p = points[i];
    Edge e = boost::edge(vertices_.at(p.get<0>()),
                         vertices_.at(p.get<1>()),
                         g_copy_).first;
    int s = boost::get(vnames_, boost::source(e, g_copy_));
    std::vector<std::pair<int, double> > &vec =
        landmarks_[boost::get(enames_, e)];

    if (vec.empty())
      vec.push_back({s, 0.0});

    vec.push_back({i, s == p.get<0>() ? p.get<2>() : 1 - p.get<2>()});
  }

  int eid = boost::num_edges(g_copy_);
  idbase_ = sim_.g_.OBJECTID;

  for (auto i = landmarks_.begin(); i != landmarks_.end(); ++i) {
    Edge e = edges_.at(i->first);
    double w = get(weights_, e);
    std::vector<std::pair<int, double> > &vec = i->second;

    std::stable_sort(vec.begin(), vec.end(),
                     [](const std::pair<int, double> &a,
                        const std::pair<int, double> &b) {
                       return a.second < b.second;
                     });

    boost::remove_edge(e, g_copy_);
    Vertex u = vertices_.at(vec.begin()->first);

    for (auto j = std::next(vec.begin()); j != vec.end(); ++j, eid++) {
      Vertex v = boost::add_vertex(idbase_ + j->first, g_copy_);
      Edge t = boost::add_edge(u, v, eid, g_copy_).first;
      put(weights_, t, w * (j->second - std::prev(j)->second));
      u = v;
    }

    // find the right `target'
    Vertex v = vertices_.at(vec.begin()->first);
    if (boost::target(e, sim_.g_()) == v)
      v = boost::source(e, sim_.g_());
    else v = boost::target(e, sim_.g_());

    Edge t = boost::add_edge(u, v, eid, g_copy_).first;
    put(weights_, t, w * (1 - vec.back().second));
  }
}

struct found_nn {};

class nearest_neighbor_visitor : public boost::default_dijkstra_visitor
{
 public:
  nearest_neighbor_visitor(int k, boost::unordered_set<int> &result)
      : k_(k), result_(result) { }

  template <typename V, typename G>
  void
  discover_vertex(V v, const G &g)
  {
    if (vnames_[v] >= idbase_) {
      result_.insert(vnames_[v]);
      if (--k_ == 0)
        throw found_nn();
    }
  }

 private:
  int k_;
  boost::unordered_set<int> &result_;
};

boost::unordered_set<int>
NearestNeighbor::query(int k)
{
  boost::unordered_set<int> result;
  nearest_neighbor_visitor vis(k, result);
  boost::dijkstra_shortest_paths(
      g_copy_, vertices_[object_],
      boost::visitor(vis).weight_map(weights_));
  return result;
}

boost::unordered_map<int, double>
NearestNeighbor::predict(int k)
{
  boost::unordered_map<int, double> result;
  return result;
}

}
