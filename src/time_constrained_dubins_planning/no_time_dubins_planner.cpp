#include "mt_tsp_ros2/time_constrained_dubins_planning/no_time_dubins_planner.h"

NoTimeDubinsPlanner::NoTimeDubinsPlanner(double rho, RowMatrixXbRef_const occupancy, Vector2dRef_const map_lb, Vector2dRef_const map_ub) {
  space = std::make_shared<DubinsStateSpaceDeterministicSampling>(rho);
  ob::RealVectorBounds bounds(2);
  bounds.setLow(0, map_lb(0));
  bounds.setLow(1, map_lb(1));
  bounds.setHigh(0, map_ub(0));
  bounds.setHigh(1, map_ub(1));
  space->setBounds(bounds);

  si = std::make_shared<ob::SpaceInformation>(space);

  state_checker = std::make_shared<NoTimeDubinsStateValidityChecker>(si, occupancy, map_lb, map_ub);
  si->setStateValidityChecker(state_checker);

  motion_validator = std::make_shared<NoTimeDubinsMotionValidator>(si, rho);
  si->setMotionValidator(motion_validator);

  ss = std::make_shared<og::SimpleSetup>(si);
  // planner = std::make_shared<og::ABITstar>(si);
  // planner = std::make_shared<og::RRTstar>(si);
  planner = std::make_shared<CustomRRTstar>(si);
  // planner->setNearestNeighbors<ompl::NearestNeighborsSqrtApprox>();
  // planner = std::make_shared<og::RRTConnect>(si);
  ss->setPlanner(planner);
}

// My occupancy grid to rectangle code is in python, otherwise I wouldn't have a separate constructor here
NoTimeDubinsPlanner::NoTimeDubinsPlanner(double rho, RowMatrixXbRef_const occupancy, RowMatrixXdRef_const rects, Vector2dRef_const map_lb, Vector2dRef_const map_ub) {
  space = std::make_shared<DubinsStateSpaceDeterministicSampling>(rho);
  ob::RealVectorBounds bounds(2);
  bounds.setLow(0, map_lb(0));
  bounds.setLow(1, map_lb(1));
  bounds.setHigh(0, map_ub(0));
  bounds.setHigh(1, map_ub(1));
  space->setBounds(bounds);

  si = std::make_shared<ob::SpaceInformation>(space);

  state_checker = std::make_shared<NoTimeDubinsStateValidityChecker>(si, occupancy, map_lb, map_ub);
  si->setStateValidityChecker(state_checker);

  motion_validator = std::make_shared<NoTimeDubinsMotionValidatorRects>(si, rho, rects, map_lb, map_ub);
  si->setMotionValidator(motion_validator);

  ss = std::make_shared<og::SimpleSetup>(si);

  // planner = std::make_shared<og::ABITstar>(si);
  // planner = std::make_shared<og::RRTstar>(si);
  planner = std::make_shared<CustomRRTstar>(si);
  // planner->setNearestNeighbors<ompl::NearestNeighborsSqrtApprox>();
  // planner = std::make_shared<og::RRTConnect>(si);
  ss->setPlanner(planner);
}

RowMatrixXd NoTimeDubinsPlanner::plan(VectorXdRef_const start, VectorXdRef_const goal, double time_limit, int max_iter) {
  ompl::msg::setLogLevel(ompl::msg::LogLevel::LOG_NONE);
  ob::ScopedState<> start_state(space);
  ob::ScopedState<> goal_state(space);
  for (int state_idx = 0; state_idx < start.size(); ++state_idx) {
    start_state[state_idx] = start(state_idx);
    goal_state[state_idx] = goal(state_idx);
  }

  ss->clear();
  ss->setStartAndGoalStates(start_state, goal_state);

  ob::PlannerStatus status = ss->solve(ob::plannerOrTerminationCondition(ob::timedPlannerTerminationCondition(time_limit), ob::IterationTerminationCondition(max_iter)));

  if (status == ob::PlannerStatus::StatusType::EXACT_SOLUTION) {
    // std::cout << "Found solution:" << std::endl;
    og::PathGeometric &solutionPath = ss->getSolutionPath();
    int num_steps = solutionPath.getStateCount();
    RowMatrixXd ret(num_steps, 3);
    for (int step = 0; step < num_steps; ++step) {
      const ob::State* state = solutionPath.getState(step);
      ret(step, 0) = state->as<ob::SE2StateSpace::StateType>()->getX();
      ret(step, 1) = state->as<ob::SE2StateSpace::StateType>()->getY();
      ret(step, 2) = state->as<ob::SE2StateSpace::StateType>()->getYaw();
    }
    return ret;
    // ss->getSolutionPath().print(std::cout);
  } else {
    // std::cout << "No solution found" << std::endl;
    return std::numeric_limits<double>::infinity()*RowMatrixXd::Ones(1, 3);
  }
}

// Test functions
bool NoTimeDubinsPlanner::is_state_valid(VectorXdRef_const state_vec) {
  ob::ScopedState<> state(space);
  for (int state_idx = 0; state_idx < state_vec.size(); ++state_idx) {
    state[state_idx] = state_vec(state_idx);
  }
  return state_checker->isValid(state.get());
}

double NoTimeDubinsPlanner::checkMotion(VectorXdRef s_valid_vec, VectorXdRef_const s1_vec, VectorXdRef_const s2_vec) {
  ob::ScopedState<> s1(space);
  ob::ScopedState<> s2(space);
  ob::ScopedState<> s_valid(space);

  for (int state_idx = 0; state_idx < s1_vec.size(); ++state_idx) {
    s1[state_idx] = s1_vec(state_idx);
    s2[state_idx] = s2_vec(state_idx);
    s_valid[state_idx] = s2_vec(state_idx);
  }
  std::pair<ob::State*, double> lastValid;
  lastValid.first = s_valid.get();
  lastValid.second = 1.;
  bool valid = motion_validator->checkMotion(s1.get(), s2.get(), lastValid);
  for (int state_idx = 0; state_idx < s1_vec.size(); ++state_idx) {
    s_valid_vec[state_idx] = s_valid[state_idx];
  }
  return lastValid.second;
}

bool NoTimeDubinsPlanner::collision_free_get_intersection_point(double x, double y, double theta, double turn_dir, double turn_dist, double next_x, double next_y, double &x_collision, double &y_collision, double &theta_collision, double &collision_dist) {
  return std::static_pointer_cast<NoTimeDubinsMotionValidatorRects>(motion_validator)->collision_free_get_intersection_point(x, y, theta, turn_dir, turn_dist, next_x, next_y, x_collision, y_collision, theta_collision, collision_dist);
}

bool NoTimeDubinsPlanner::collision_free(double x, double y, double theta, double turn_dir, double turn_dist, double next_x, double next_y) {
  return std::static_pointer_cast<NoTimeDubinsMotionValidatorRects>(motion_validator)->collision_free(x, y, theta, turn_dir, turn_dist, next_x, next_y);
}
