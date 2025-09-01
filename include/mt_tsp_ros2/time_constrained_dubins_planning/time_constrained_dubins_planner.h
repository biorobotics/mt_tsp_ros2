#pragma once
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_time_state_space.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_motion_validator.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_motion_validator_rects.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_state_sampler.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_state_validity_checker.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/custom_rrt_connect.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/custom_rrt.h"
#include <ompl/base/spaces/RealVectorBounds.h>
#include <ompl/geometric/planners/rrt/RRTConnect.h>
#include <ompl/geometric/planners/rrt/RRT.h>
#include <ompl/geometric/SimpleSetup.h>
#include <ompl/base/ScopedState.h>
#include <ompl/base/PlannerStatus.h>
#include <ompl/base/PlannerTerminationCondition.h>
#include <ompl/base/terminationconditions/IterationTerminationCondition.h>

namespace ob = ompl::base;
namespace og = ompl::geometric;

typedef Ref<VectorXd> VectorXdRef;
typedef const Ref<const RowMatrixXd>& RowMatrixXdRef_const;

class TimeConstrainedDubinsPlanner {
  public:
    TimeConstrainedDubinsPlanner(double vmax, double rho, RowMatrixXbRef_const occupancy, Vector2dRef_const map_lb, Vector2dRef_const map_ub, bool nn_sort_by_time, bool use_approx_dist);

    // My occupancy grid to rectangle code is in python, otherwise I wouldn't have a separate constructor here
    TimeConstrainedDubinsPlanner(double vmax, double rho, RowMatrixXbRef_const occupancy, RowMatrixXdRef_const rects, Vector2dRef_const map_lb, Vector2dRef_const map_ub, bool monte_carlo_prop, bool nn_sort_by_time, bool use_approx_dist);

    int get_num_discarded_samples();

    int get_num_samples();

    RowMatrixXd plan(VectorXdRef_const start, VectorXdRef_const goal, double time_limit, int max_iter, std::vector<RowMatrixXd> &turns_chain);

    double get_num_tree_nodes() const;

    double get_nearest_neighbor_time() const;

    double get_path_elongation_time() const;

    double get_path_elongation_check_time() const;

    double get_collision_check_time() const;

    double get_sampling_time() const;

    double get_add_to_tree_time() const;

    // Test functions
    bool is_state_valid(VectorXdRef_const state_vec);

    VectorXd sample_random_state(VectorXdRef_const start, VectorXdRef_const goal);

    double checkMotion(VectorXdRef s_valid_vec, VectorXdRef_const s1_vec, VectorXdRef_const s2_vec);

    bool checkMotionForwardBackward(VectorXdRef s_valid_vec, VectorXdRef_const s1_vec, VectorXdRef_const s2_vec, double maxDuration, bool forward);

  private:
    double vmax;
    double rho;
    RowMatrixXd rects;
    bool do_rects;
    Vector2d map_lb;
    Vector2d map_ub;
    bool nn_sort_by_time;
    bool use_approx_dist;

    std::shared_ptr<DubinsTimeStateSpace> space;
    ob::SpaceInformationPtr si;
    std::shared_ptr<DubinsStateValidityChecker> state_checker;
    std::shared_ptr<DubinsMotionValidator> motion_validator;

    std::shared_ptr<og::SimpleSetup> ss;
    std::shared_ptr<ob::Planner> planner;

    double nearest_neighbor_time;
    double path_elongation_time;
    double path_elongation_check_time;
    double collision_check_time;
    double sampling_time;
    double add_to_tree_time;
};
