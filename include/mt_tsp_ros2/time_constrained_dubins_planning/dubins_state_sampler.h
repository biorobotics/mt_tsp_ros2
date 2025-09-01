#pragma once
#include <ompl/base/StateSampler.h>
#include <ompl/base/samplers/DeterministicStateSampler.h>
#include <ompl/util/Exception.h>
#include <ompl/base/samplers/deterministic/HaltonSequence.h>
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_time_state_space.h"
#include "mt_tsp_ros2/elongate_dubins_path.h"
#include <chrono>

namespace ob = ompl::base;
 
class DubinsTimeStateSampler : public ob::DeterministicStateSampler {
  public:
    DubinsTimeStateSampler(const ob::StateSpace *space, 
                           double start_x, double start_y, double start_theta, double start_t,
                           double goal_x, double goal_y, double goal_theta, double goal_t,
                           double vmax, double rho);

    void sampleUniform(ob::State *state) override;

    // We don't need this
    void sampleUniformNear(ob::State* state, const ob::State* near, const double distance) override;

    // We don't need this
    void sampleGaussian(ob::State* state, const ob::State* mean, const double stdDev) override;

    double get_path_elongation_intervals_time() const;

    void reset_timing_info();

  private:
    std::shared_ptr<ob::SE2DeterministicStateSampler> config_sampler;
    std::shared_ptr<ob::HaltonSequence1D> time_sequence_ptr;
    double start_x;
    double start_y;
    double start_theta;
    double start_t;
    double goal_x;
    double goal_y;
    double goal_theta;
    double goal_t;
    double vmax;
    double rho;
    mutable double path_elongation_intervals_time;
};
