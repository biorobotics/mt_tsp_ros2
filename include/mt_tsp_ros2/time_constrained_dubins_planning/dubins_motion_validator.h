#pragma once
#include <ompl/base/SpaceInformation.h>
#include <ompl/base/MotionValidator.h>
#include "mt_tsp_ros2/elongate_dubins_path.h"

namespace ob = ompl::base;

class DubinsMotionValidator : public ob::MotionValidator {
  public:
    explicit DubinsMotionValidator(const ob::SpaceInformationPtr si, double vmax, double rho) : ob::MotionValidator(si), vmax(vmax), rho(rho) {
    }

    bool checkMotion(const ob::State *s1, const ob::State *s2, std::pair<ob::State*, double> &lastValid) const override {
      // This function assumes s1 is valid, but not necessarily s2
      if (!si_->isValid(s2)) {
        return false;
      }

      double x1 = s1->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getX();
      double y1 = s1->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getY();
      double theta1 = s1->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getYaw();
      double t1 = s1->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position;

      double x2 = s2->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getX();
      double y2 = s2->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getY();
      double theta2 = s2->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getYaw();
      double t2 = s2->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position;

      RowMatrixXd turns = elongated_dubins_path(x1, y1, theta1, x2, y2, theta2, vmax*(t2 - t1), rho, false);
      if (std::isinf(turns(0, 0))) {
        if (lastValid.first != nullptr) {
          lastValid.first->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setX(x1);
          lastValid.first->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setY(y1);
          lastValid.first->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setYaw(theta1);
          lastValid.first->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position = t1;
        }
        lastValid.second = 0.;
        return false;
      }

      ob::State* next_s = si_->getStateSpace()->allocState();

      double x = x1;
      double y = y1;
      double theta = theta1;
      double t = t1;

      double valid_x = x1;
      double valid_y = y1;
      double valid_theta = theta1;
      double valid_t = t1;

      double valid_path_length = 0.;

      for (int turn_idx = 0; turn_idx < turns.rows(); ++turn_idx) {
        double turn_dir = turns(turn_idx, 0);
        double turn_dist = turns(turn_idx, 1);
        double rho_times_turn_dir = rho*turn_dir;

        double ctheta = cos(theta);
        double stheta = sin(theta);

        int num_checks = 100;
        double next_x;
        double next_y;
        double next_theta;
        double next_t;
        for (int check_idx = 0; check_idx < num_checks; ++check_idx) {
          double dist = (check_idx + 1)/(double)num_checks*turn_dist;
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

          next_s->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setX(next_x);
          next_s->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setY(next_y);
          next_s->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setYaw(next_theta);
          next_s->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position = next_t;

          if (!si_->isValid(next_s)) {
            if (lastValid.first != nullptr) {
              lastValid.first->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setX(valid_x);
              lastValid.first->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setY(valid_y);
              lastValid.first->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setYaw(valid_theta);
              lastValid.first->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position = valid_t;
            }
            double path_length = turns.col(1).sum();
            lastValid.second = valid_path_length/path_length;
            si_->getStateSpace()->freeState(next_s);
            return false;
          }

          valid_x = next_x;
          valid_y = next_y;
          valid_theta = next_theta;
          valid_t = next_t;

          valid_path_length = turns.col(1).head(turn_idx).sum() + dist;
        }
        x = next_x;
        y = next_y;
        theta = next_theta;
        t = next_t;
      }

      si_->getStateSpace()->freeState(next_s);
      return true;
    }

    bool checkMotion(const ob::State *s1, const ob::State *s2) const override {
      std::pair<ob::State*, double> lastValid;
      return checkMotion(s1, s2, lastValid);
    }

  private:
    double vmax;
    double rho;
};
