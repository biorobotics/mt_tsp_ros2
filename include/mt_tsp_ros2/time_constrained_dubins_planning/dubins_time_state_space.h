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
    explicit DubinsTimeStateSpace(const std::shared_ptr<ob::SE2StateSpace>& spaceComponent, double vmax, double rho) : vmax(vmax), rho(rho), wmax(vmax/rho) {
      addSubspace(spaceComponent, 0.5); // space component
      addSubspace(std::make_shared<ob::TimeStateSpace>(), 0.5); // time component
      // lock();

      path_elongation_check_time = 0.;
      path_elongation_time = 0.;
    }

    void set_start_and_goal(double start_x, double start_y, double start_theta, double start_t,
                            double goal_x, double goal_y, double goal_theta, double goal_t) {
      this->start_x = start_x;
      this->start_y = start_y;
      this->start_theta = start_theta;
      this->start_t = start_t;
      this->goal_x = goal_x;
      this->goal_y = goal_y;
      this->goal_theta = goal_theta;
      this->goal_t = goal_t;
      this->as<ob::TimeStateSpace>(1)->setBounds(start_t, goal_t);
    }

    ob::StateSamplerPtr allocStateSampler() const override {
      return std::make_shared<DubinsTimeStateSampler>(this,
                                                      start_x, start_y, start_theta, start_t,
                                                      goal_x, goal_y, goal_theta, goal_t,
                                                      vmax, rho);
    }

     bool isMetricSpace() const override {
         return false;
     }

     bool hasSymmetricDistance() const override {
         return false;
     }

     bool hasSymmetricInterpolate() const override {
         return false;
     }

    double distance(const ob::State *state1, const ob::State *state2) const override {
      /*
      double x1 = state1->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getX();
      double y1 = state1->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getY();
      double theta1 = state1->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getYaw();
      double t1 = state1->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position;

      double x2 = state2->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getX();
      double y2 = state2->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getY();
      double theta2 = state2->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getYaw();
      double t2 = state2->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position;
      double delta_x = x2 - x1;
      double delta_y = y2 - y1;
      double travel_time_lb = sqrt(delta_x*delta_x + delta_y*delta_y)/vmax;
      double delta_t = t2 - t1;
      if (delta_t < travel_time_lb) {
        return std::numeric_limits<double>::infinity();
      }
      return delta_t;
      */
      double x1 = state1->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getX();
      double y1 = state1->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getY();
      double theta1 = state1->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getYaw();
      double t1 = state1->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position;

      double x2 = state2->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getX();
      double y2 = state2->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getY();
      double theta2 = state2->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getYaw();
      double t2 = state2->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position;

      auto timer_start = std::chrono::high_resolution_clock::now();
      double delta_x = x2 - x1;
      double delta_y = y2 - y1;
      double delta_theta = angdiff(theta1, theta2);
      double travel_time_lb = std::max(sqrt(delta_x*delta_x + delta_y*delta_y)/vmax, std::abs(delta_theta)/wmax);
      double delta_t = t2 - t1;
      if (delta_t < travel_time_lb || !check_elongation_possible(x1, y1, theta1, x2, y2, theta2, vmax*(t2 - t1), rho)) {
        auto timer_stop = std::chrono::high_resolution_clock::now();
        auto nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
        path_elongation_check_time += ((double)nanos)/1e9;
        return std::numeric_limits<double>::infinity();
      }

      auto timer_stop = std::chrono::high_resolution_clock::now();
      auto nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
      path_elongation_check_time += ((double)nanos)/1e9;
      return t2 - t1;
    }

    void interpolate(const ob::State *s1, const ob::State *s2, double fraction, ob::State *interp_state) const override {
      double x1 = s1->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getX();
      double y1 = s1->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getY();
      double theta1 = s1->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getYaw();
      double t1 = s1->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position;

      double x2 = s2->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getX();
      double y2 = s2->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getY();
      double theta2 = s2->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getYaw();
      double t2 = s2->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position;

      auto timer_start = std::chrono::high_resolution_clock::now();
      RowMatrixXd turns = elongated_dubins_path(x1, y1, theta1, x2, y2, theta2, vmax*(t2 - t1), rho, false);
      auto timer_stop = std::chrono::high_resolution_clock::now();
      auto nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
      path_elongation_time += ((double)nanos)/1e9;
      if (std::isinf(turns(0, 0))) {
        interp_state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setX(x1);
        interp_state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setY(y1);
        interp_state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setYaw(theta1);
        interp_state->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position = t1;
        return;
      }

      double x = x1;
      double y = y1;
      double theta = theta1;
      double t = t1;

      double path_length = turns.col(1).sum();

      for (int turn_idx = 0; turn_idx < turns.rows(); ++turn_idx) {
        double turn_dir = turns(turn_idx, 0);
        double turn_dist = turns(turn_idx, 1);

        bool interp_point_on_this_segment = turns.col(1).head(turn_idx).sum() + turn_dist >= fraction*path_length;

        double dist;
        if (interp_point_on_this_segment) {
          // fraction*path_length = turns.col(1).head(turn_idx).sum() + dist
          dist = fraction*path_length - turns.col(1).head(turn_idx).sum();
        } else {
          dist = turn_dist;
        }

        double rho_times_turn_dir = rho*turn_dir;

        double ctheta = cos(theta);
        double stheta = sin(theta);

        double next_x;
        double next_y;
        double next_theta;
        double next_t;

        next_t = t + dist/vmax;
        if (turn_dir == 0) {
          // S segment
          next_theta = theta;
          next_x = x + dist*ctheta;
          next_y = y + dist*stheta;
        } else {
          // C segment
          next_theta = theta + dist/rho_times_turn_dir;
          next_x = x + rho_times_turn_dir*(-stheta + sin(next_theta));
          next_y = y + rho_times_turn_dir*(ctheta - cos(next_theta));
        }

        if (interp_point_on_this_segment) {
          interp_state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setX(next_x);
          interp_state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setY(next_y);
          interp_state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setYaw(next_theta);
          interp_state->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position = next_t;
          return;
        } else {
          x = next_x;
          y = next_y;
          theta = next_theta;
          t = next_t;
        }
      }

      throw ompl::Exception("DubinsTimeStateSpace::interpolate", "reached line of code that should be impossible to reach");
    }

    double get_path_elongation_check_time() const {
      return path_elongation_check_time;
    }

    void reset_timing_info() {
      path_elongation_check_time = 0.;
    }

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
