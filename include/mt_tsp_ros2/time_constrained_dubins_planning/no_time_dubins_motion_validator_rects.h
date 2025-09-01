#pragma once
#include <ompl/base/SpaceInformation.h>
#include <ompl/base/MotionValidator.h>
#include "mt_tsp_ros2/elongate_dubins_path.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/angle_mod.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/no_time_dubins_motion_validator.h"
#include <chrono>

#include "mt_tsp_ros2/time_constrained_dubins_planning/geometry_utils.h"

namespace ob = ompl::base;

class NoTimeDubinsMotionValidatorRects : public NoTimeDubinsMotionValidator {
  public:
    explicit NoTimeDubinsMotionValidatorRects(const ob::SpaceInformationPtr si, double rho, RowMatrixXdRef_const &rects, Vector2dRef_const map_lb, Vector2dRef_const map_ub);

    const Tree &get_aabb_tree();

    bool collision_free(double x, double y, double theta, double turn_dir, double turn_dist, double next_x, double next_y) const;

    bool collision_free_get_intersection_point(double x, double y, double theta, double turn_dir, double turn_dist, double next_x, double next_y, double &x_collision, double &y_collision, double &theta_collision, double &collision_dist) const;

    bool checkMotion(const ob::State *s1, const ob::State *s2, std::pair<ob::State*, double> &lastValid) const override;

  protected:
    RowMatrixXd rects;
    Vector2d map_lb;
    Vector2d map_ub;

    // Need segments to persist in memory while the tree is in use
    std::list<Segment> segments;
    Tree aabb_tree;
    mutable ompl::RNG rng_;
};
