#include "mt_tsp_ros2/time_constrained_dubins_planning/time_constrained_dubins_planner.h"

TimeConstrainedDubinsPlanner::TimeConstrainedDubinsPlanner(double vmax, double rho, RowMatrixXbRef_const occupancy, Vector2dRef_const map_lb, Vector2dRef_const map_ub, bool nn_sort_by_time, bool use_approx_dist) : vmax(vmax), rho(rho), do_rects(false), map_lb(map_lb), map_ub(map_ub), nn_sort_by_time(nn_sort_by_time), use_approx_dist(use_approx_dist) {
  std::shared_ptr<ob::SE2StateSpace> cspace = std::make_shared<ob::SE2StateSpace>();
  ob::RealVectorBounds bounds(2);
  bounds.setLow(0, map_lb(0));
  bounds.setLow(1, map_lb(1));
  bounds.setHigh(0, map_ub(0));
  bounds.setHigh(1, map_ub(1));
  cspace->setBounds(bounds);
  space = std::make_shared<DubinsTimeStateSpace>(cspace, vmax, rho);

  si = std::make_shared<ob::SpaceInformation>(space);

  state_checker = std::make_shared<DubinsStateValidityChecker>(si, occupancy, map_lb, map_ub);
  si->setStateValidityChecker(state_checker);

  motion_validator = std::make_shared<DubinsMotionValidator>(si, vmax, rho);
  si->setMotionValidator(motion_validator);

  ss = std::make_shared<og::SimpleSetup>(si);
  planner = std::make_shared<CustomRRTConnect>(si, nn_sort_by_time, use_approx_dist);
  // planner = std::make_shared<og::RRTConnect>(si);
  // planner = std::make_shared<og::RRT>(si);
  // planner = std::make_shared<CustomRRT>(si);
  ss->setPlanner(planner);

  nearest_neighbor_time = 0.;
  path_elongation_time = 0.;
  path_elongation_check_time = 0.;
  collision_check_time = 0.;
}

// My occupancy grid to rectangle code is in python, otherwise I wouldn't have a separate constructor here
TimeConstrainedDubinsPlanner::TimeConstrainedDubinsPlanner(double vmax, double rho, RowMatrixXbRef_const occupancy, RowMatrixXdRef_const rects, Vector2dRef_const map_lb, Vector2dRef_const map_ub, bool monte_carlo_prop, bool nn_sort_by_time, bool use_approx_dist) : vmax(vmax), rho(rho), rects(rects), do_rects(true), map_lb(map_lb), map_ub(map_ub), nn_sort_by_time(nn_sort_by_time), use_approx_dist(use_approx_dist) {
  std::shared_ptr<ob::SE2StateSpace> cspace = std::make_shared<ob::SE2StateSpace>();
  ob::RealVectorBounds bounds(2);
  bounds.setLow(0, map_lb(0));
  bounds.setLow(1, map_lb(1));
  bounds.setHigh(0, map_ub(0));
  bounds.setHigh(1, map_ub(1));
  cspace->setBounds(bounds);
  space = std::make_shared<DubinsTimeStateSpace>(cspace, vmax, rho);

  si = std::make_shared<ob::SpaceInformation>(space);

  motion_validator = std::make_shared<DubinsMotionValidatorRects>(si, vmax, rho, rects, map_lb, map_ub, VectorXd::Zero(4), VectorXd::Zero(4), monte_carlo_prop);
  si->setMotionValidator(motion_validator);

  ss = std::make_shared<og::SimpleSetup>(si);

  state_checker = std::make_shared<DubinsStateValidityChecker>(si, occupancy, map_lb, map_ub);
  si->setStateValidityChecker(state_checker);

  path_elongation_time = 0.;
  path_elongation_check_time = 0.;
  collision_check_time = 0.;
}

int TimeConstrainedDubinsPlanner::get_num_discarded_samples() {
  return motion_validator->get_num_discarded_samples();
}

int TimeConstrainedDubinsPlanner::get_num_samples() {
  return std::static_pointer_cast<CustomRRTConnect>(planner)->get_num_samples();
}

RowMatrixXd TimeConstrainedDubinsPlanner::plan(VectorXdRef_const start, VectorXdRef_const goal, double time_limit, int max_iter, std::vector<RowMatrixXd> &turns_chain) {
  // I'm avoiding reconstructing the motion validator for now so we don't have to reconstruct the AABB tree
  /*
  if (do_rects) {
    motion_validator = std::make_shared<DubinsMotionValidatorRects>(si, vmax, rho, rects, map_lb, map_ub, start, goal);
    si->setMotionValidator(motion_validator);
  }
  */

  space->reset_timing_info();
  motion_validator->reset_timing_info();
  motion_validator->reset_num_discarded_samples();

  ompl::msg::setLogLevel(ompl::msg::LogLevel::LOG_NONE);
  planner = std::make_shared<CustomRRTConnect>(si, nn_sort_by_time, use_approx_dist);
  ss->setPlanner(planner);
  space->set_start_and_goal(start(0), start(1), start(2), start(3), 
                            goal(0), goal(1), goal(2), goal(3));
  ob::ScopedState<> start_state(space);
  ob::ScopedState<> goal_state(space);
  for (int state_idx = 0; state_idx < start.size(); ++state_idx) {
    start_state[state_idx] = start(state_idx);
    goal_state[state_idx] = goal(state_idx);
  }

  ss->clear();
  ss->setStartAndGoalStates(start_state, goal_state);

  // std::cout << "Solve" << std::endl;
  ob::PlannerStatus status;
  // Just use a time limit
  if (max_iter == -1) {
     status = ss->solve(ob::timedPlannerTerminationCondition(time_limit));
  } else {
     status = ss->solve(ob::plannerOrTerminationCondition(ob::timedPlannerTerminationCondition(time_limit), ob::IterationTerminationCondition(max_iter)));
  }
  // std::cout << "Solved" << std::endl;

  nearest_neighbor_time = std::static_pointer_cast<CustomRRTConnect>(planner)->get_nearest_neighbor_time();
  path_elongation_time = motion_validator->get_path_elongation_time();
  path_elongation_check_time = space->get_path_elongation_check_time();
  collision_check_time = motion_validator->get_collision_check_time();
  sampling_time = std::static_pointer_cast<CustomRRTConnect>(planner)->get_sampling_time();
  add_to_tree_time = std::static_pointer_cast<CustomRRTConnect>(planner)->get_add_to_tree_time();

  // std::cout << "got profiling data" << std::endl;

  if (status == ob::PlannerStatus::StatusType::EXACT_SOLUTION) {
    // std::cout << "Got soln" << std::endl;
    // std::cout << "Found solution:" << std::endl;
    DubinsSegmentChain &solutionPath = static_cast<DubinsSegmentChain&>(ss->getSolutionPath());
    int num_steps = solutionPath.getStateCount();
    RowMatrixXd ret(num_steps, 4);
    for (int step = 0; step < num_steps; ++step) {
      const ob::State* state = solutionPath.getState(step);
      ret(step, 0) = state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getX();
      ret(step, 1) = state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getY();
      ret(step, 2) = state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getYaw();
      ret(step, 3) = state->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position;
    }
    turns_chain = solutionPath.turns_chain;
    // std::cout << "returning soln" << std::endl;
    return ret;
    // ss->getSolutionPath().print(std::cout);
  } else {
    // std::cout << "No solution found" << std::endl;
    // std::cout << "returning inf" << std::endl;
    return std::numeric_limits<double>::infinity()*RowMatrixXd::Ones(1, 4);
  }
}

double TimeConstrainedDubinsPlanner::get_num_tree_nodes() const {
  return std::static_pointer_cast<CustomRRTConnect>(planner)->get_num_tree_nodes();
}

double TimeConstrainedDubinsPlanner::get_nearest_neighbor_time() const {
  return nearest_neighbor_time;
}

double TimeConstrainedDubinsPlanner::get_path_elongation_time() const {
  return path_elongation_time;
}

double TimeConstrainedDubinsPlanner::get_path_elongation_check_time() const {
  return path_elongation_check_time;
}

double TimeConstrainedDubinsPlanner::get_collision_check_time() const {
  return collision_check_time;
}

double TimeConstrainedDubinsPlanner::get_sampling_time() const {
  return sampling_time;
}

double TimeConstrainedDubinsPlanner::get_add_to_tree_time() const {
  return add_to_tree_time;
}

// Test functions
bool TimeConstrainedDubinsPlanner::is_state_valid(VectorXdRef_const state_vec) {
  ob::ScopedState<> state(space);
  for (int state_idx = 0; state_idx < state_vec.size(); ++state_idx) {
    state[state_idx] = state_vec(state_idx);
  }
  return state_checker->isValid(state.get());
}

VectorXd TimeConstrainedDubinsPlanner::sample_random_state(VectorXdRef_const start, VectorXdRef_const goal) {
  space->set_start_and_goal(start(0), start(1), start(2), start(3), 
                            goal(0), goal(1), goal(2), goal(3));
  ob::StateSamplerPtr sampler = si->allocStateSampler();
  ob::ScopedState<> sampled_state(space);
  std::cout << "sampling state" << std::endl;
  sampler->sampleUniform(sampled_state.get());
  std::cout << "sampled state" << std::endl;
  VectorXd ret(4);
  for (int state_idx = 0; state_idx < ret.size(); ++state_idx) {
    ret(state_idx) = sampled_state[state_idx];
  }
  return ret;
}

double TimeConstrainedDubinsPlanner::checkMotion(VectorXdRef s_valid_vec, VectorXdRef_const s1_vec, VectorXdRef_const s2_vec) {
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

bool TimeConstrainedDubinsPlanner::checkMotionForwardBackward(VectorXdRef s_valid_vec, VectorXdRef_const s1_vec, VectorXdRef_const s2_vec, double maxDuration, bool forward) {
  ob::ScopedState<> s1(space);
  ob::ScopedState<> s2(space);
  ob::ScopedState<> s_valid(space);

  for (int state_idx = 0; state_idx < s1_vec.size(); ++state_idx) {
    s1[state_idx] = s1_vec(state_idx);
    s2[state_idx] = s2_vec(state_idx);
    s_valid[state_idx] = s2_vec(state_idx);
  }
  bool reach;
  bool validMotion;
  if (forward) {
    motion_validator->checkMotionForward(s1.get(), s2.get(), maxDuration, s_valid.get(), reach, validMotion);
    return validMotion;
  } else {
    motion_validator->checkMotionBackward(s1.get(), s2.get(), maxDuration, s_valid.get(), reach, validMotion);
    return validMotion;
  }
}
