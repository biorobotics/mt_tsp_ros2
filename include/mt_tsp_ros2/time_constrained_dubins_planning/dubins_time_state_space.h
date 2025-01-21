#pragma once
#include <ompl/base/spaces/SE2StateSpace.h>
#include <ompl/base/spaces/TimeStateSpace.h>
#include "mt_tsp_ros2/elongate_dubins_path.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_state_sampler.h"

namespace ob = ompl::base;

class DubinsTimeStateSpace : public ob::CompoundStateSpace {
  public:
    explicit DubinsTimeStateSpace(const std::shared_ptr<ob::SE2StateSpace>& spaceComponent, double vmax, double rho) : vmax(vmax), rho(rho) {
      addSubspace(spaceComponent, 0.5); // space component
      addSubspace(std::make_shared<ob::TimeStateSpace>(), 0.5); // time component
      // lock();
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

  private:
    double vmax;
    double rho;
    double start_x;
    double start_y;
    double start_theta;
    double start_t;
    double goal_x;
    double goal_y;
    double goal_theta;
    double goal_t;
};
