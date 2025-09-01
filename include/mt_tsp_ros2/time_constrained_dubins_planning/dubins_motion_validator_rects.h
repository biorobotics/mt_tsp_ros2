#pragma once
#include <ompl/base/SpaceInformation.h>
#include <ompl/base/MotionValidator.h>
#include "mt_tsp_ros2/elongate_dubins_path.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/angle_mod.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_motion_validator.h"
#include <chrono>

#include "mt_tsp_ros2/time_constrained_dubins_planning/geometry_utils.h"

namespace ob = ompl::base;

class DubinsMotionValidatorRects : public DubinsMotionValidator {
  public:
    explicit DubinsMotionValidatorRects(const ob::SpaceInformationPtr si, double vmax, double rho, RowMatrixXdRef_const &rects, Vector2dRef_const map_lb, Vector2dRef_const map_ub, const Ref<const VectorXd> &start, const Ref<const VectorXd> &goal, bool monte_carlo_prop);

    const Tree &get_aabb_tree();

    bool collision_free(double x, double y, double theta, double turn_dir, double turn_dist, double next_x, double next_y) const;

    // If we set validMotion = false, stopState isn't used so we don't have to populate it. Same deal with turns
    RowMatrixXd checkMotionForward(const ob::State *s1, const ob::State *s2, double maxDuration, ob::State *stopState, bool &reach, bool &validMotion, int num_checks = 1000) const override;

    RowMatrixXd checkMotionForward_internal(const ob::State *s1, const ob::State *s2, double maxDuration, ob::State *stopState, bool &reach, bool &validMotion) const;

    RowMatrixXd checkMotionBackward(const ob::State *s1, const ob::State *s2, double maxDuration, ob::State *stopState, bool &reach, bool &validMotion, int num_checks = 1000) const override;

    RowMatrixXd checkMotionBackward_internal(const ob::State *s1, const ob::State *s2, double maxDuration, ob::State *stopState, bool &reach, bool &validMotion) const;

    void checkMotionForward_connect(std::vector<RowMatrixXd> &turns_vec, std::vector<Vector4d> &states_vec, const ob::State *s1, const ob::State *s2, double maxDuration, bool &reach) const;

    void checkMotionBackward_connect(std::vector<RowMatrixXd> &turns_vec, std::vector<Vector4d> &states_vec, const ob::State *s1, const ob::State *s2, double maxDuration, bool &reach) const;

    void set_monte_carlo_prop(bool monte_carlo_prop);

  protected:
    RowMatrixXd rects;
    Vector2d map_lb;
    Vector2d map_ub;
    VectorXd start;
    VectorXd goal;
    bool monte_carlo_prop;

    // Need segments to persist in memory while the tree is in use
    std::list<Segment> segments;
    Tree aabb_tree;
    mutable ompl::RNG rng_;
};
