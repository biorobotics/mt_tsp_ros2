#pragma once
 
#include <ompl/control/planners/rrt/RRT.h>
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_motion_validator_rects.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_control_sampler.h"
#include <ompl/control/spaces/RealVectorControlSpace.h>

using namespace ompl::control;
namespace ob = ompl::base;
 
class CustomControlRRT : public ompl::control::RRT
{
public:
    CustomControlRRT(const SpaceInformationPtr &si, std::shared_ptr<DubinsMotionValidatorRects> motion_validator);

    ob::PlannerStatus solve(const ob::PlannerTerminationCondition &ptc) override;

    double get_collision_check_time();

protected:
    using ompl::control::RRT::rng_;
    std::shared_ptr<DubinsMotionValidatorRects> motion_validator;
};
