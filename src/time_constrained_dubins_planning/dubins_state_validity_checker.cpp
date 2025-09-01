#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_state_validity_checker.h"

DubinsStateValidityChecker::DubinsStateValidityChecker(const ob::SpaceInformationPtr &si, RowMatrixXbRef_const &occupancy, Vector2dRef_const map_lb, Vector2dRef_const map_ub) : ob::StateValidityChecker(si), occupancy(occupancy), map_lb(map_lb), map_ub(map_ub),
                                                                                                                                                     cell_size((map_ub(0) - map_lb(0))/occupancy.rows(), (map_ub(1) - map_lb(1))/occupancy.cols()) {
}

bool DubinsStateValidityChecker::isValid(const ob::State *state) const {
  double x = state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getX();
  double y = state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getY();
  if (x < map_lb(0) || x > map_ub(0) || y < map_lb(1) || y > map_ub(1)) {
    return false;
  }
  int xcell = (int)((x - map_lb(0))/cell_size(0));
  int ycell = (int)((y - map_lb(1))/cell_size(1));
  return !occupancy(xcell, ycell);
}
