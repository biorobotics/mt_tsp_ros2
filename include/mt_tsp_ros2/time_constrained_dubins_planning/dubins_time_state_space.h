#pragma once
#include <ompl/base/spaces/SE2StateSpace.h>
#include <ompl/base/spaces/TimeStateSpace.h>
#include "mt_tsp_ros2/elongate_dubins_path.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_state_sampler.h"
#include <chrono>
#include "mt_tsp_ros2/time_constrained_dubins_planning/angle_mod.h"

namespace ob = ompl::base;

class DubinsTimeStateSpace : public ob::CompoundStateSpace {
  public:
    explicit DubinsTimeStateSpace(const std::shared_ptr<ob::SE2StateSpace>& spaceComponent, double vmax, double rho);

    void set_start_and_goal(double start_x, double start_y, double start_theta, double start_t,
                            double goal_x, double goal_y, double goal_theta, double goal_t);

    ob::StateSamplerPtr allocStateSampler() const override;

     bool isMetricSpace() const override;

     bool hasSymmetricDistance() const override;

     bool hasSymmetricInterpolate() const override;

    double distance(const ob::State *state1, const ob::State *state2) const override;

    double approx_distance(const ob::State *state1, const ob::State *state2) const;

    void interpolate(const ob::State *s1, const ob::State *s2, double fraction, ob::State *interp_state) const override;

    double get_path_elongation_check_time() const;

    void reset_timing_info();

  private:
    double vmax;
    double rho;
    double wmax;
    double start_x;
    double start_y;
    double start_theta;
    double start_t;
    double goal_x;
    double goal_y;
    double goal_theta;
    double goal_t;

    mutable double path_elongation_check_time;
    mutable double path_elongation_time;
};
