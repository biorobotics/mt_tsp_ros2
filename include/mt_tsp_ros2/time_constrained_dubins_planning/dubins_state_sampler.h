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
                           double vmax, double rho) : ob::DeterministicStateSampler(space),
                                                      start_x(start_x), start_y(start_y), start_theta(start_theta), start_t(start_t),
                                                      goal_x(goal_x), goal_y(goal_y), goal_theta(goal_theta), goal_t(goal_t),
                                                      vmax(vmax), rho(rho) {

      config_sampler = std::make_shared<ob::SE2DeterministicStateSampler>(space->as<ob::CompoundStateSpace>()->as<ob::SE2StateSpace>(0));
      time_sequence_ptr = std::make_shared<ob::HaltonSequence1D>();
      /*
      std::cout << check_elongation_possible(start_x, start_y, start_theta, goal_x, goal_y, goal_theta, vmax*(goal_t - start_t), rho) << std::endl;;
      std::cout << turns_for_dubins_path(start_x, start_y, start_theta, goal_x, goal_y, goal_theta, rho).col(1).sum() << std::endl;
      std::cout << vmax*(goal_t - start_t) << std::endl;
      */

      path_elongation_intervals_time = 0.;
    }

    void sampleUniform(ob::State *state) override {
      // Sample a configuration q = (x, y, t) uniformly at random
      config_sampler->sampleUniform(state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0));
      double raw_sample = time_sequence_ptr->sample();
      state->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position = start_t + raw_sample*(goal_t - start_t);

      /*
      int max_iter = 1000;
      for (int i = 0; i < max_iter; ++i) {
        // Sample a configuration q = (x, y, t) uniformly at random
        config_sampler->sampleUniform(state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0));

        double x = state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getX();
        double y = state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getY();
        double theta = state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getYaw();

        auto timer_start = std::chrono::high_resolution_clock::now();
        Vector3d start_elongation_intervals = start_t*Vector3d::Ones() + get_elongation_intervals(start_x, start_y, start_theta, x, y, theta, rho)/vmax;
        Vector3d goal_elongation_intervals = goal_t*Vector3d::Ones() - get_elongation_intervals(x, y, theta, goal_x, goal_y, goal_theta, rho).reverse()/vmax;
        auto timer_stop = std::chrono::high_resolution_clock::now();
        auto nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
        path_elongation_intervals_time += ((double)nanos)/1e9;

        std::vector<Vector2d> valid_intervals;

        // Check intersection of early start interval and early goal interval
        Vector2d intersect1;
        if (std::isinf(goal_elongation_intervals(0))) {
          // There is no early goal interval, only a late one
          intersect1(0) = std::numeric_limits<double>::infinity();
          intersect1(1) = -std::numeric_limits<double>::infinity();
        } {
          // There is an early goal interval
          intersect1(0) = start_elongation_intervals(0);
          intersect1(1) = std::min(start_elongation_intervals(1), goal_elongation_intervals(0));
        }

        if (intersect1(1) >= intersect1(0)) {
          valid_intervals.push_back(intersect1);
        }

        // Check intersection of early start interval and late goal interval
        Vector2d intersect2 = Vector2d(std::max(start_elongation_intervals(0), goal_elongation_intervals(1)), std::min(start_elongation_intervals(1), goal_elongation_intervals(2)));

        if (intersect2(1) >= intersect2(0)) {
          valid_intervals.push_back(intersect2);
        }

        // Check intersection of late start interval and early goal interval
        Vector2d intersect3;
        if (std::isinf(goal_elongation_intervals(0))) {
          // There is no early goal interval, only a late one
          intersect3(0) = std::numeric_limits<double>::infinity();
          intersect3(1) = -std::numeric_limits<double>::infinity();
        } else if (std::isinf(start_elongation_intervals(1))) {
          // There is no late start interval, only a early one
          intersect3(0) = std::numeric_limits<double>::infinity();
          intersect3(1) = -std::numeric_limits<double>::infinity();
        } else {
          // There is a late start interval and early goal interval
          intersect3(0) = start_elongation_intervals(2);
          intersect3(1) = goal_elongation_intervals(0);
        }

        if (intersect3(1) >= intersect3(0)) {
          valid_intervals.push_back(intersect3);
        }

        // Check intersection of late start interval and late goal interval
        Vector2d intersect4;
        if (std::isinf(start_elongation_intervals(1))) {
          // There is no late start interval, only a early one
          intersect4(0) = std::numeric_limits<double>::infinity();
          intersect4(1) = -std::numeric_limits<double>::infinity();
        } else {
          // There is a late start interval
          intersect4(0) = std::max(start_elongation_intervals(2), goal_elongation_intervals(1));
          intersect4(1) = goal_elongation_intervals(2);
        }

        if (intersect4(1) >= intersect4(0)) {
          valid_intervals.push_back(intersect4);
        }

        if (valid_intervals.size() == 0) {
          continue;
        }

        double raw_sample = time_sequence_ptr->sample();
        int interval_idx = (int)(raw_sample*(double)valid_intervals.size());
        raw_sample = time_sequence_ptr->sample();
        double t = valid_intervals[interval_idx](0) + raw_sample*(valid_intervals[interval_idx](1) - valid_intervals[interval_idx](0));

        state->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position = t;

        return;
      }

      double raw_sample = time_sequence_ptr->sample();

      ob::State* start_state = space_->allocState();
      start_state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setX(start_x);
      start_state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setY(start_y);
      start_state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setYaw(start_theta);
      start_state->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position = start_t;

      ob::State* goal_state = space_->allocState();
      goal_state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setX(goal_x);
      goal_state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setY(goal_y);
      goal_state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setYaw(goal_theta);
      goal_state->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position = goal_t;

      space_->interpolate(start_state, goal_state, raw_sample, state);
      space_->freeState(start_state);
      space_->freeState(goal_state);
      */
    }

    // We don't need this
    void sampleUniformNear(ob::State* state, const ob::State* near, const double distance) override {
      throw ompl::Exception("DubinsTimeStateSampler::sampleUniformNear", "not implemented");
    }

    // We don't need this
    void sampleGaussian(ob::State* state, const ob::State* mean, const double stdDev) override {
      throw ompl::Exception("DubinsTimeStateSampler::sampleGaussian", "not implemented");
    }

    double get_path_elongation_intervals_time() const {
      return path_elongation_intervals_time;
    }

    void reset_timing_info() {
      path_elongation_intervals_time = 0.;
    }

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
