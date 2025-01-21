#pragma once
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_time_state_space.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_motion_validator.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_state_sampler.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_state_validity_checker.h"
#include <ompl/base/spaces/RealVectorBounds.h>
#include <ompl/geometric/planners/rrt/RRTConnect.h>
#include <ompl/geometric/SimpleSetup.h>
#include <ompl/base/ScopedState.h>
#include <ompl/base/PlannerStatus.h>
#include <ompl/base/PlannerTerminationCondition.h>
#include <ompl/base/terminationconditions/IterationTerminationCondition.h>

namespace ob = ompl::base;
namespace og = ompl::geometric;

typedef Ref<VectorXd> VectorXdRef;

class TimeConstrainedDubinsPlanner {
  public:
    TimeConstrainedDubinsPlanner(double vmax, double rho, RowMatrixXbRef_const occupancy, Vector2dRef_const map_lb, Vector2dRef_const map_ub) {
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
      rrt_connect = std::make_shared<og::RRTConnect>(si);
      ss->setPlanner(rrt_connect);
    }

    RowMatrixXd plan(VectorXdRef_const start, VectorXdRef_const goal, double time_limit, int max_iter) {
      space->set_start_and_goal(start(0), start(1), start(2), start(3), 
                                goal(0), goal(1), goal(2), goal(3));

      ob::ScopedState<> start_state(space);
      ob::ScopedState<> goal_state(space);
      for (int state_idx = 0; state_idx < start.size(); ++state_idx) {
        start_state[state_idx] = start(state_idx);
        goal_state[state_idx] = goal(state_idx);
      }

      ss->setStartAndGoalStates(start_state, goal_state);

      // attempt to solve the problem within one second of planning time
      ob::PlannerStatus status = ss->solve(ob::plannerOrTerminationCondition(ob::timedPlannerTerminationCondition(time_limit), ob::IterationTerminationCondition(max_iter)));

      if (status == ob::PlannerStatus::StatusType::EXACT_SOLUTION) {
        std::cout << "Found solution:" << std::endl;
        // print the path to screen
        ss->getSolutionPath().print(std::cout);
      } else {
        std::cout << "No solution found" << std::endl;
      }
      return RowMatrixXd::Zero(1, 1);
    }

    // Test functions
    bool is_state_valid(VectorXdRef_const state_vec) {
      ob::ScopedState<> state(space);
      for (int state_idx = 0; state_idx < state_vec.size(); ++state_idx) {
        state[state_idx] = state_vec(state_idx);
      }
      return state_checker->isValid(state.get());
    }

    VectorXd sample_random_state(VectorXdRef_const start, VectorXdRef_const goal) {
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
    std::shared_ptr<DubinsTimeStateSpace> space;
    ob::SpaceInformationPtr si;
    std::shared_ptr<DubinsStateValidityChecker> state_checker;
    std::shared_ptr<DubinsMotionValidator> motion_validator;

    std::shared_ptr<og::SimpleSetup> ss;
    std::shared_ptr<og::RRTConnect> rrt_connect;
};
