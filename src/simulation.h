#ifndef SIMULATION_H_
#define SIMULATION_H_

#include <boost/parameter/name.hpp>
#include <boost/parameter/preprocessor.hpp>

#include <boost/tuple/tuple.hpp>

#include <CGAL/Point_set_2.h>

#include "walkinggraph.h"
#include "particle.h"

namespace simulation {

BOOST_PARAMETER_NAME(num_object)
BOOST_PARAMETER_NAME(num_particle)
BOOST_PARAMETER_NAME(radius)
BOOST_PARAMETER_NAME(unit)
BOOST_PARAMETER_NAME(knock_door_prob)
BOOST_PARAMETER_NAME(enter_room_prob)
BOOST_PARAMETER_NAME(threshold)
BOOST_PARAMETER_NAME(success_rate)

class Simulation_impl_
{
 public:

  void
  run(double duration);

  void
  snapshot(double t);

  void
  random_window(double ratio) const;

  std::vector<int>
  range_query();

  std::map<int, double>
  range_query_pred();

  void
  random_object();

  std::vector<int>
  nearest_neighbors(int id, int k);

  std::map<int, double>
  nearest_neighbors_pred(int k);

 protected:

  Simulation_impl_() {}

  template <class ArgumentPack>
  Simulation_impl_(const ArgumentPack &args)
      : num_object_(args[_num_object])
      , num_particle_(args[_num_particle | 128])
      , radius_(args[_radius | 120.0])
      , unit_(args[_unit | 20])
      , knock_door_prob_(args[_knock_door_prob | 0.1])
      , enter_room_prob_(args[_enter_room_prob | 0.1])
      , threshold_(args[_threshold | 0.4])
      , success_rate_(args[_success_rate | 0.95])
  { }

  landmark_t
  random_inside_reader(int i) const;

  void
  detect();

 private:

  // system parameters
  int num_object_;
  int num_particle_;
  double radius_;
  double unit_;
  double knock_door_prob_;
  double enter_room_prob_;
  double threshold_;
  double success_rate_;

  WalkingGraph g_;

  // objects in the system
  std::vector<Particle> objects_;
};

class Simulation : public Simulation_impl_
{
 public:
  BOOST_PARAMETER_CONSTRUCTOR(
      Simulation, (Simulation_impl_), tag,
      (required
       (num_object, *))
      (optional
       (num_particle, *)
       (radius, *)
       (unit, *)
       (knock_door_prob, *)
       (enter_room_prob, *)
       (threshold, *)))
};

}

#endif  // SIMULATION_H_
