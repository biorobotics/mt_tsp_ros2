#pragma once
#include "mt_tsp_ros2/time_constrained_dubins_planning/no_time_dubins_motion_validator.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/no_time_dubins_state_validity_checker.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_state_space_deterministic_sampling.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/custom_rrt_star.h"
#include <ompl/base/spaces/RealVectorBounds.h>
#include <ompl/geometric/planners/rrt/RRTConnect.h>
#include <ompl/geometric/planners/rrt/RRTstar.h>
#include <ompl/geometric/planners/informedtrees/ABITstar.h>
#include <ompl/geometric/SimpleSetup.h>
#include <ompl/base/ScopedState.h>
#include <ompl/base/PlannerStatus.h>
#include <ompl/base/PlannerTerminationCondition.h>
#include <ompl/base/terminationconditions/IterationTerminationCondition.h>
#include <ompl/datastructures/NearestNeighborsSqrtApprox.h>

namespace ob = ompl::base;
namespace og = ompl::geometric;

typedef Ref<VectorXd> VectorXdRef;

class NoTimeDubinsPlanner {
  public:
    NoTimeDubinsPlanner(double rho, RowMatrixXbRef_const occupancy, Vector2dRef_const map_lb, Vector2dRef_const map_ub) {
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

    RowMatrixXd plan(VectorXdRef_const start, VectorXdRef_const goal, double time_limit, int max_iter) {
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
    bool is_state_valid(VectorXdRef_const state_vec) {
      ob::ScopedState<> state(space);
      for (int state_idx = 0; state_idx < state_vec.size(); ++state_idx) {
        state[state_idx] = state_vec(state_idx);
      }
      return state_checker->isValid(state.get());
    }

    double checkMotion(VectorXdRef s_valid_vec, VectorXdRef_const s1_vec, VectorXdRef_const s2_vec) {
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

  private:
    std::shared_ptr<DubinsStateSpaceDeterministicSampling> space;
    ob::SpaceInformationPtr si;
    std::shared_ptr<NoTimeDubinsStateValidityChecker> state_checker;
    std::shared_ptr<NoTimeDubinsMotionValidator> motion_validator;

    std::shared_ptr<og::SimpleSetup> ss;
    std::shared_ptr<CustomRRTstar> planner;
    // std::shared_ptr<og::RRTstar> planner;
};
