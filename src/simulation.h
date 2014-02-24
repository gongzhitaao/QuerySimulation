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
  nearest_neighbor(int k);

  std::map<int, double>
  nearest_neighbor_pred(int k);

 protected:

  template <class ArgumentPack>
  simulation_impl(const ArgumentPack &args)
      : num_object_(args[_num_object])
      , num_particle_(args[_num_particle])
      , radius_(args[_radius])
      , unit_(args[_unit])
      , knock_door_prob_(args[_knock_door_prob])
      , enter_room_prob_(args[_enter_room_prob])
      , threshold_(args[_threshold])
      , success_rate_(args[_success_rate])
  {
    initialize();
  }

  void
  initialize();

  void
  detect();

  landmark_t
  random_inside_reader(int i) const;

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
      myclass, (myclass_impl), tag,
      (required (num_object,*))
      (optional
       (num_particle, *, 128)
       (radius, *, 120.0)
       (unit, *, 20.0)
       (knock_door_prob, *, 0.1)
       (enter_room_prob, *, 0.1),
       (threshold, *, 0.4)))
};

}

#endif  // SIMULATION_H_
