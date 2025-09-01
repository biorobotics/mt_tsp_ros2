#include "mt_tsp_ros2/time_constrained_dubins_planning/time_constrained_dubins_planner_forward_prop.h"

// Dummy function, not used
void propagate(const ob::State *start, const oc::Control *control, const double duration, ob::State *result)
{
}

TimeConstrainedDubinsPlannerForwardProp::TimeConstrainedDubinsPlannerForwardProp(double vmax, double rho, RowMatrixXbRef_const occupancy, RowMatrixXdRef_const rects, Vector2dRef_const map_lb, Vector2dRef_const map_ub, bool monte_carlo_prop) : vmax(vmax), rho(rho), rects(rects), do_rects(true), map_lb(map_lb), map_ub(map_ub) {
  std::shared_ptr<ob::SE2StateSpace> cspace = std::make_shared<ob::SE2StateSpace>();
  ob::RealVectorBounds bounds(2);
  bounds.setLow(0, map_lb(0));
  bounds.setLow(1, map_lb(1));
  bounds.setHigh(0, map_ub(0));
  bounds.setHigh(1, map_ub(1));
  cspace->setBounds(bounds);
  space = std::make_shared<DubinsTimeStateSpace>(cspace, vmax, rho);

  // The first control is angular velocity, the second is duration.
  // Duration bounds will depend on the start and goal, so don't set them here
  control_dim = 2;
  control_space = std::make_shared<oc::RealVectorControlSpace>(space, control_dim);

  si = std::make_shared<oc::SpaceInformation>(space, control_space);

  state_checker = std::make_shared<DubinsStateValidityChecker>(si, occupancy, map_lb, map_ub);
  si->setStateValidityChecker(state_checker);
  si->setStatePropagator(propagate);

  motion_validator = std::make_shared<DubinsMotionValidatorRects>(si, vmax, rho, rects, map_lb, map_ub, VectorXd::Zero(4), VectorXd::Zero(4), monte_carlo_prop);

  ss = std::make_shared<oc::SimpleSetup>(si);
  planner = std::make_shared<CustomControlRRT>(si, motion_validator);
  ss->setPlanner(planner);

  path_elongation_time = 0.;
  path_elongation_check_time = 0.;
  collision_check_time = 0.;
}

RowMatrixXd TimeConstrainedDubinsPlannerForwardProp::plan(VectorXdRef_const start, VectorXdRef_const goal, double time_limit, int max_iter, std::vector<RowMatrixXd> &turns_chain) {
  space->reset_timing_info();
  motion_validator->reset_timing_info();

  // Angular velocity bound and time bound. The control space's
  // control bounds aren't used anymore
  /*
  ob::RealVectorBounds control_bounds(control_dim);
  control_bounds.setLow(0, -vmax/rho);
  control_bounds.setLow(1, 0.);
  control_bounds.setHigh(0, vmax/rho);
  control_bounds.setHigh(1, 0.2*(goal(3) - start(3)));
  si->getControlSpace()->setBounds(control_bounds);
  */

  si->setDirectedControlSamplerAllocator(DubinsControlSamplerAllocator(motion_validator, vmax/rho, 0.2*(goal(3) - start(3))));

  ompl::msg::setLogLevel(ompl::msg::LogLevel::LOG_NONE);
  planner = std::make_shared<CustomControlRRT>(si, motion_validator);
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
  ob::PlannerStatus status = ss->solve(ob::plannerOrTerminationCondition(ob::timedPlannerTerminationCondition(time_limit), ob::IterationTerminationCondition(max_iter)));
  // std::cout << "Solved" << std::endl;

  path_elongation_time = motion_validator->get_path_elongation_time();
  path_elongation_check_time = space->get_path_elongation_check_time();
  collision_check_time = std::static_pointer_cast<CustomControlRRT>(planner)->get_collision_check_time() + motion_validator->get_collision_check_time();

  if (status == ob::PlannerStatus::StatusType::EXACT_SOLUTION) {
    // std::cout << "Got soln" << std::endl;
    // std::cout << "Found solution:" << std::endl;
    PathControl &solutionPath = ss->getSolutionPath();
    int num_steps = solutionPath.getStateCount();
    RowMatrixXd ret(num_steps, 4);
    for (int step = 0; step < num_steps; ++step) {
      const ob::State* state = solutionPath.getState(step);
      ret(step, 0) = state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getX();
      ret(step, 1) = state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getY();
      ret(step, 2) = state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getYaw();
      ret(step, 3) = state->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position;
    }
    // std::cout << "returning soln" << std::endl;
    return ret;
    // ss->getSolutionPath().print(std::cout);
  } else {
    // std::cout << "No solution found" << std::endl;
    // std::cout << "returning inf" << std::endl;
    return std::numeric_limits<double>::infinity()*RowMatrixXd::Ones(1, 4);
  }
}

double TimeConstrainedDubinsPlannerForwardProp::get_path_elongation_time() const {
  return path_elongation_time;
}

double TimeConstrainedDubinsPlannerForwardProp::get_path_elongation_check_time() const {
  return path_elongation_check_time;
}

double TimeConstrainedDubinsPlannerForwardProp::get_collision_check_time() const {
  return collision_check_time;
}
