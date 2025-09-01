#pragma once
#include "mt_tsp_ros2/time_constrained_dubins_planning/no_time_dubins_motion_validator.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/no_time_dubins_motion_validator_rects.h"
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
    NoTimeDubinsPlanner(double rho, RowMatrixXbRef_const occupancy, Vector2dRef_const map_lb, Vector2dRef_const map_ub);

    // My occupancy grid to rectangle code is in python, otherwise I wouldn't have a separate constructor here
    NoTimeDubinsPlanner(double rho, RowMatrixXbRef_const occupancy, RowMatrixXdRef_const rects, Vector2dRef_const map_lb, Vector2dRef_const map_ub);

    RowMatrixXd plan(VectorXdRef_const start, VectorXdRef_const goal, double time_limit, int max_iter);

    // Test functions
    bool is_state_valid(VectorXdRef_const state_vec);

    double checkMotion(VectorXdRef s_valid_vec, VectorXdRef_const s1_vec, VectorXdRef_const s2_vec);

    bool collision_free_get_intersection_point(double x, double y, double theta, double turn_dir, double turn_dist, double next_x, double next_y, double &x_collision, double &y_collision, double &theta_collision, double &collision_dist);

  private:
    std::shared_ptr<DubinsStateSpaceDeterministicSampling> space;
    ob::SpaceInformationPtr si;
    std::shared_ptr<NoTimeDubinsStateValidityChecker> state_checker;
    std::shared_ptr<NoTimeDubinsMotionValidator> motion_validator;

    std::shared_ptr<og::SimpleSetup> ss;
    std::shared_ptr<CustomRRTstar> planner;
    // std::shared_ptr<og::RRTstar> planner;
};
