#pragma once
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_time_state_space.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_motion_validator_rects.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_state_sampler.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_state_validity_checker.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/custom_control_rrt.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_control_sampler.h"
#include <ompl/base/spaces/RealVectorBounds.h>
#include <ompl/control/SimpleSetup.h>
#include <ompl/base/ScopedState.h>
#include <ompl/base/PlannerStatus.h>
#include <ompl/base/PlannerTerminationCondition.h>
#include <ompl/base/terminationconditions/IterationTerminationCondition.h>
#include <ompl/control/spaces/RealVectorControlSpace.h>

namespace ob = ompl::base;
namespace oc = ompl::control;

typedef Ref<VectorXd> VectorXdRef;
typedef const Ref<const RowMatrixXd>& RowMatrixXdRef_const;

class TimeConstrainedDubinsPlannerForwardProp {
  public:
    TimeConstrainedDubinsPlannerForwardProp(double vmax, double rho, RowMatrixXbRef_const occupancy, RowMatrixXdRef_const rects, Vector2dRef_const map_lb, Vector2dRef_const map_ub, bool monte_carlo_prop);

    RowMatrixXd plan(VectorXdRef_const start, VectorXdRef_const goal, double time_limit, int max_iter, std::vector<RowMatrixXd> &turns_chain);

    double get_path_elongation_time() const;

    double get_path_elongation_check_time() const;

    double get_collision_check_time() const;

  private:
    double vmax;
    double rho;
    RowMatrixXd rects;
    bool do_rects;
    Vector2d map_lb;
    Vector2d map_ub;

    int control_dim = 2;

    std::shared_ptr<DubinsTimeStateSpace> space;
    std::shared_ptr<ControlSpace> control_space;
    oc::SpaceInformationPtr si;
    std::shared_ptr<DubinsStateValidityChecker> state_checker;
    std::shared_ptr<DubinsMotionValidatorRects> motion_validator;

    std::shared_ptr<oc::SimpleSetup> ss;
    std::shared_ptr<ob::Planner> planner;

    double path_elongation_time;
    double path_elongation_check_time;
    double collision_check_time;
};
