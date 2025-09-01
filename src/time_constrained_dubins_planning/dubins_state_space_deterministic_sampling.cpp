#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_state_space_deterministic_sampling.h"

DubinsStateSpaceDeterministicSampling::DubinsStateSpaceDeterministicSampling(double rho) : ob::DubinsStateSpace(rho) {
}

ob::StateSamplerPtr DubinsStateSpaceDeterministicSampling::allocStateSampler() const {
  return std::make_shared<ob::SE2DeterministicStateSampler>(this);
}
