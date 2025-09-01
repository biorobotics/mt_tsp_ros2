#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_motion_validator.h"
#include <ompl/base/spaces/SE2StateSpace.h>
#include <ompl/base/spaces/TimeStateSpace.h>

DubinsMotionValidator::DubinsMotionValidator(const ob::SpaceInformationPtr si, double vmax, double rho) : ob::MotionValidator(si), vmax(vmax), rho(rho) {
  path_elongation_time = 0.;
  collision_check_time = 0.;

  num_discarded_samples = 0;
}

RowMatrixXd DubinsMotionValidator::checkMotionForward(const ob::State *s1, const ob::State *s2, double maxDuration, ob::State *stopState, bool &reach, bool &validMotion, int num_checks) const {
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
    reach = false;
    stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setX(x1);
    stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setY(y1);
    stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setYaw(theta1);
    stopState->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position = t1;
    validMotion = false;
    ++num_discarded_samples;
    return turns;
  }

  /*
  std::cout << "checking forward truncated" << std::endl;
  std::cout << x1 << " " << y1 << " " << theta1 << " " << t1 << std::endl;
  std::cout << x2 << " " << y2 << " " << theta2 << " " << t2 << std::endl;
  std::cout << turns << std::endl;
  */

  timer_start = std::chrono::high_resolution_clock::now();

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
        reach = false;
        stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setX(valid_x);
        stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setY(valid_y);
        stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setYaw(valid_theta);
        stopState->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position = valid_t;
        si_->getStateSpace()->freeState(next_s);
        validMotion = false;
        turns(turn_idx, 1) = dist;

        timer_stop = std::chrono::high_resolution_clock::now();
        nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
        collision_check_time += ((double)nanos)/1e9;
        return turns.topRows(turn_idx + 1);
      }

      valid_x = next_x;
      valid_y = next_y;
      valid_theta = next_theta;
      valid_t = next_t;

      /*
      if (turn_idx == 2) {
        std::cout << next_x << " " << next_y << " " << next_theta << " " << next_t << std::endl;
      }
      */

      valid_path_length = turns.col(1).head(turn_idx).sum() + dist;

      if (valid_path_length/vmax >= maxDuration) {
        reach = false;
        stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setX(valid_x);
        stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setY(valid_y);
        stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setYaw(valid_theta);
        stopState->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position = valid_t;
        si_->getStateSpace()->freeState(next_s);
        /*
        std::cout << "advance" << std::endl;
        std::cout << next_x << " " << next_y << " " << next_theta << " " << next_t << std::endl;
        */
        validMotion = true;
        turns(turn_idx, 1) = dist;

        timer_stop = std::chrono::high_resolution_clock::now();
        nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
        collision_check_time += ((double)nanos)/1e9;

        return turns.topRows(turn_idx + 1);
      }
    }
    x = next_x;
    y = next_y;
    theta = next_theta;
    t = next_t;
  }

  reach = true;
  stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setX(valid_x);
  stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setY(valid_y);
  stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setYaw(valid_theta);
  stopState->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position = valid_t;
  si_->getStateSpace()->freeState(next_s);
  validMotion = true;

  timer_stop = std::chrono::high_resolution_clock::now();
  nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
  collision_check_time += ((double)nanos)/1e9;

  return turns;
}

RowMatrixXd DubinsMotionValidator::checkMotionBackward(const ob::State *s1, const ob::State *s2, double maxDuration, ob::State *stopState, bool &reach, bool &validMotion, int num_checks) const {
  if (!si_->isValid(s2)) {
    reach = false;
    validMotion = false;
    return std::numeric_limits<double>::infinity()*RowMatrixXd::Ones(1, 2);
  }
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
  if (std::isinf(turns(0, 0))) {
    reach = false;
    stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setX(x2);
    stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setY(y2);
    stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setYaw(theta2);
    stopState->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position = t2;
    validMotion = false;
    ++num_discarded_samples;
    return turns;
  }

  timer_start = std::chrono::high_resolution_clock::now();

  ob::State* next_s = si_->getStateSpace()->allocState();

  double x = x2;
  double y = y2;
  double theta = theta2;
  double t = t2;

  double valid_x = x2;
  double valid_y = y2;
  double valid_theta = theta2;
  double valid_t = t2;

  double valid_path_length = 0.;

  // std::cout << "checking backwards truncated" << std::endl;
  for (int turn_idx = turns.rows() - 1; turn_idx >= 0; --turn_idx) {
    double turn_dir = turns(turn_idx, 0);
    double turn_dist = turns(turn_idx, 1);
    double rho_times_turn_dir = rho*turn_dir;

    double ctheta = cos(theta);
    double stheta = sin(theta);

    double next_x;
    double next_y;
    double next_theta;
    double next_t;
    for (int check_idx = 0; check_idx < num_checks; ++check_idx) {
      double dist = (check_idx + 1)/(double)num_checks*turn_dist;
      next_t = t - dist/vmax;
      if (turn_dir == 0) {
        // S segment
        next_theta = theta;
        next_x = x - dist*ctheta;
        next_y = y - dist*stheta;
      } else {
        // C segment
        next_theta = theta - dist/rho_times_turn_dir;
        /*
        next_x = x + rho_times_turn_dir*(-stheta + sin(next_theta));
        next_y = y + rho_times_turn_dir*(ctheta - cos(next_theta));
        */
        next_x = x - rho_times_turn_dir*(stheta - sin(next_theta));
        next_y = y - rho_times_turn_dir*(-ctheta + cos(next_theta));
      }

      next_s->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setX(next_x);
      next_s->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setY(next_y);
      next_s->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setYaw(next_theta);
      next_s->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position = next_t;

      if (!si_->isValid(next_s)) {
        reach = false;
        stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setX(valid_x);
        stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setY(valid_y);
        stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setYaw(valid_theta);
        stopState->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position = valid_t;
        si_->getStateSpace()->freeState(next_s);
        validMotion = false;
        turns(turn_idx, 1) = dist;

        timer_stop = std::chrono::high_resolution_clock::now();
        nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
        collision_check_time += ((double)nanos)/1e9;

        return turns.bottomRows(turns.rows() - turn_idx);
      }

      valid_x = next_x;
      valid_y = next_y;
      valid_theta = next_theta;
      valid_t = next_t;

      valid_path_length = turns.col(1).tail(turns.rows() - 1 - turn_idx).sum() + dist;

      if (valid_path_length/vmax >= maxDuration) {
        reach = false;
        stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setX(valid_x);
        stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setY(valid_y);
        stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setYaw(valid_theta);
        stopState->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position = valid_t;
        si_->getStateSpace()->freeState(next_s);
        validMotion = true;
        turns(turn_idx, 1) = dist;

        timer_stop = std::chrono::high_resolution_clock::now();
        nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
        collision_check_time += ((double)nanos)/1e9;

        return turns.bottomRows(turns.rows() - turn_idx);
      }
    }
    x = next_x;
    y = next_y;
    theta = next_theta;
    t = next_t;
  }

  reach = true;
  stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setX(valid_x);
  stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setY(valid_y);
  stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setYaw(valid_theta);
  stopState->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position = valid_t;
  si_->getStateSpace()->freeState(next_s);
  validMotion = true;

  timer_stop = std::chrono::high_resolution_clock::now();
  nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
  collision_check_time += ((double)nanos)/1e9;

  return turns;
}

bool DubinsMotionValidator::checkMotion(const ob::State *s1, const ob::State *s2, std::pair<ob::State*, double> &lastValid) const {
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

  /*
  std::cout << "checking forward" << std::endl;
  std::cout << x1 << " " << y1 << " " << theta1 << " " << t1 << std::endl;
  std::cout << x2 << " " << y2 << " " << theta2 << " " << t2 << std::endl;
  std::cout << turns << std::endl;
  */
  // std::cout << turns << std::endl;
  for (int turn_idx = 0; turn_idx < turns.rows(); ++turn_idx) {
    double turn_dir = turns(turn_idx, 0);
    double turn_dist = turns(turn_idx, 1);
    double rho_times_turn_dir = rho*turn_dir;

    double ctheta = cos(theta);
    double stheta = sin(theta);

    int num_checks = 1000;
    double next_x;
    double next_y;
    double next_theta;
    double next_t;
    for (int check_idx = 0; check_idx < num_checks; ++check_idx) {
      double dist = (check_idx + 1)/(double)num_checks*turn_dist; // This is like doing np.linspace with num_checks + 1 points and ignoring the first element because we already checked it
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

      /*
      if (turn_idx == 2) {
        std::cout << next_x << " " << next_y << " " << next_theta << " " << next_t << std::endl;
      }
      */
      // std::cout << next_x << " " << next_y << " " << next_theta << " " << next_t << " " << dist << std::endl;
      /*
      if (turn_idx == 3) {
        std::cout << next_x << " " << next_y << " " << next_theta << " " << next_t << std::endl;
      }
      */

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
        // std::cout << "failed at turn idx" << turn_idx << std::endl;
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

bool DubinsMotionValidator::checkMotion(const ob::State *s1, const ob::State *s2) const {
  std::pair<ob::State*, double> lastValid;
  return checkMotion(s1, s2, lastValid);
}

double DubinsMotionValidator::get_path_elongation_time() const {
  return path_elongation_time;
}

double DubinsMotionValidator::get_collision_check_time() const {
  return collision_check_time;
}

void DubinsMotionValidator::reset_timing_info() {
  path_elongation_time = 0.;
  collision_check_time = 0.;
}

double DubinsMotionValidator::get_vmax() const {
  return vmax;
}

int DubinsMotionValidator::get_num_discarded_samples() {
  return num_discarded_samples;
}

void DubinsMotionValidator::reset_num_discarded_samples() {
  num_discarded_samples = 0;
}
