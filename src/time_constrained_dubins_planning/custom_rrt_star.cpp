#include "mt_tsp_ros2/time_constrained_dubins_planning/custom_rrt_star.h"

CustomRRTstar::CustomRRTstar(const base::SpaceInformationPtr &si) : RRTstar(si) {
  rng_ = ompl::RNG(1);
}
