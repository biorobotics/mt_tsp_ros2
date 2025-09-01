#include "mt_tsp_ros2/time_constrained_dubins_planning/custom_rrt.h"

CustomRRT::CustomRRT(const base::SpaceInformationPtr &si) : RRT(si) {
  rng_ = ompl::RNG(1);
}
