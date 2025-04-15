#pragma once
#include <ompl/base/SpaceInformation.h>
#include <ompl/base/MotionValidator.h>
#include "mt_tsp_ros2/elongate_dubins_path.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/angle_mod.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_motion_validator.h"
#include <chrono>

namespace ob = ompl::base;

// https://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
double dist_point_to_line_segment(double x1, double y1, double x2, double y2, double x3, double y3) { // x3,y3 is the point
  double px = x2-x1;
  double py = y2-y1;

  double norm = px*px + py*py;

  if (norm == 0) {
    return sqrt((x3 - x2)*(x3 - x2) + (y3 - y2)*(y3 - y2));
  }

  double u =  ((x3 - x1) * px + (y3 - y1) * py) / norm;

  if (u > 1) {
    u = 1;
  } else if (u < 0) {
    u = 0;
  }

  double x = x1 + u * px;
  double y = y1 + u * py;

  double dx = x - x3;
  double dy = y - y3;

  // Note: If the actual distance does not matter,
  // if you only want to compare what this function
  // returns to other results of this function, you
  // can just return the squared distance instead
  // (i.e. remove the sqrt) to gain a little performance

  double dist = sqrt(dx*dx + dy*dy);

  return dist;
}

// Check if (x3, y3) is on the line segment from (x1, y1) to (x2, y2). Points don't have to be collinear
bool onSegment(double x1, double y1, double x2, double y2, double x3, double y3) {
  return dist_point_to_line_segment(x1, y1, x2, y2, x3, y3) == 0;
}

double lineSegmentsIntersect(double x1i, double y1i, double x1f, double y1f,
                             double x2i, double y2i, double x2f, double y2f) {
  double a = x1f - x1i;
  double b = x2i - x2f;
  double c = y1f - y1i;
  double d = y2i - y2f;

  double det = a*d - b*c;

  if (det == 0) {
    // Segments are parallel
    return onSegment(x1i, y1i, x2i, y2i, x1f, y1f) ||
           onSegment(x1i, y1i, x2f, y2f, x1f, y1f) ||
           onSegment(x2i, y2i, x1i, y1i, x2f, y2f) ||
           onSegment(x2i, y2i, x1f, y1f, x2f, y2f);
  }

  double lx = x2i - x1i;
  double ly = y2i - y1i;

  double t1 = (d*lx - b*ly)/det;
  double t2 = (-c*lx + a*ly)/det;
  return 0 <= t1 && t1 <= 1 && 0 <= t2 && t2 <= 1;
}

bool point_in_rect(double x, double y, \
                   double xir, double yir, double xfr, double yfr) {
  return xir <= x && x <= xfr &&
         yir <= y && y <= yfr;
}

bool line_segment_intersects_rect(double x1i, double y1i, double x1f, double y1f,
                                  double xir, double yir,
                                  double xfr, double yfr) {
  return point_in_rect(x1i, y1i,
                       xir, yir, xfr, yfr) ||
         point_in_rect(x1f, y1f, \
                       xir, yir, xfr, yfr) ||
         lineSegmentsIntersect(x1i, y1i, x1f, y1f,
                               xir, yir, xir, yfr) ||
         lineSegmentsIntersect(x1i, y1i, x1f, y1f,
                               xir, yfr, xfr, yfr) ||
         lineSegmentsIntersect(x1i, y1i, x1f, y1f,
                               xfr, yfr, xfr, yir) ||
         lineSegmentsIntersect(x1i, y1i, x1f, y1f,
                               xfr, yir, xir, yir);
}

// https://stackoverflow.com/questions/1073336/circle-line-segment-collision-detection-algorithm
Matrix2d line_segment_intersects_circle(double x1, double y1, double x2, double y2, double cx, double cy, double r) {
  Vector2d d(x2 - x1, y2 - y1);
  Vector2d f(x1 - cx, y1 - cy);
  double a = d.dot(d) ;
  double b = 2*f.dot(d);
  double c = f.dot(f) - r*r;

  Matrix2d ret = std::numeric_limits<double>::infinity()*Matrix2d::Ones();

  double discriminant = b*b-4*a*c;
  if (discriminant >= 0) {
    // line containing segment didn't totally miss sphere,
    // so there is a solution to
    // the equation.

    discriminant = sqrt(discriminant);

    // either solution may be on or off the segment so need to test both
    // t1 is always the smaller value, because BOTH discriminant and
    // a are nonnegative.
    double t1 = (-b - discriminant)/(2*a);
    double t2 = (-b + discriminant)/(2*a);

    // 3x HIT cases:
    //          -o->             --|-->  |            |  --|->
    // Impale(t1 hit,t2 hit), Poke(t1 hit,t2>1), ExitWound(t1<0, t2 hit),

    // 3x MISS cases:
    //       ->  o                     o ->              | -> |
    // FallShort (t1>1,t2>1), Past (t1<0,t2<0), CompletelyInside(t1<0, t2>1)

    if (t1 >= 0 && t1 <= 1)
    {
      // t1 is an intersection
      // Impale or Poke

      if(t2 >= 0 && t2 <= 1) {
        // Impale
        ret(0, 0) = x1 + d(0)*t1;
        ret(0, 1) = y1 + d(1)*t1;

        ret(1, 0) = x1 + d(0)*t2;
        ret(1, 1) = y1 + d(1)*t2;
      } else {
        // Poke
        ret(0, 0) = x1 + d(0)*t1;
        ret(0, 1) = y1 + d(1)*t1;
      }
    }

    // here t1 didn't intersect so we either started
    // inside the sphere or are completely past it
    if( t2 >= 0 && t2 <= 1 )
    {
      // ExitWound
      ret(0, 0) = x1 + d(0)*t2;
      ret(0, 1) = y1 + d(1)*t2;
    }

    // no intersection: FallShort, Past, CompletelyInside
  }
  return ret;
}


bool arc_intersects_line_segment(double xai, double yai, double thetaai, double turn_dir, double turn_angle_magnitude, double r,
                                 double x1, double y1,
                                 double x2, double y2) {
  double cai = cos(thetaai);
  double sai = sin(thetaai);

  Vector2d dirai(cai, sai);
  Vector2d perpai(-sai*turn_dir, cai*turn_dir);
  Vector2d centerai = Vector2d(xai, yai) + perpai*r;

  Matrix2d intersections = line_segment_intersects_circle(x1, y1, x2, y2, centerai(0), centerai(1), r);
  if (std::isinf(intersections(0, 0))) {
    return false;
  }

  double phiai = atan2(yai, xai);

  double phi1 = atan2(intersections(0, 1), intersections(0, 0));

  double diff1i = angdiff(phiai, phi1);
  if (diff1i < 0 && turn_dir == 1) {
    diff1i = diff1i + 2*M_PI;
  } else if (diff1i > 0 && turn_dir == -1) {
    diff1i = diff1i - 2*M_PI;
  }
  double abs_diff1i = std::abs(diff1i);
  if (abs_diff1i < turn_angle_magnitude) {
    return true;
  }

  if (std::isinf(intersections(1, 0))) {
    return false;
  }

  double phi2 = atan2(intersections(1, 1), intersections(1, 0));

  double diff2i = angdiff(phiai, phi2);
  if (diff2i < 0 && turn_dir == 1) {
    diff2i = diff2i + 2*M_PI;
  } else if (diff2i > 0 && turn_dir == -1) {
    diff2i = diff2i - 2*M_PI;
  }
  double abs_diff2i = std::abs(diff2i);
  if (abs_diff2i < turn_angle_magnitude) {
    return true;
  }

  return false;
}

bool arc_intersects_rect(double xai, double yai, double thetaai, double turn_dir, double turn_angle_magnitude, double r,
                         double xaf, double yaf,
                         double xir, double yir,
                         double xfr, double yfr) {
  return point_in_rect(xai, yai,
                       xir, yir, xfr, yfr) ||
         point_in_rect(xaf, yaf, \
                       xir, yir, xfr, yfr) ||
         arc_intersects_line_segment(xai, yai, thetaai, turn_dir, turn_angle_magnitude, r,
                                     xir, yir, xir, yfr) ||
         arc_intersects_line_segment(xai, yai, thetaai, turn_dir, turn_angle_magnitude, r,
                                     xir, yfr, xfr, yfr) ||
         arc_intersects_line_segment(xai, yai, thetaai, turn_dir, turn_angle_magnitude, r,
                                     xfr, yfr, xfr, yir) ||
         arc_intersects_line_segment(xai, yai, thetaai, turn_dir, turn_angle_magnitude, r,
                                     xfr, yir, xir, yir);
}

class DubinsMotionValidatorRects : public DubinsMotionValidator {
  public:
    explicit DubinsMotionValidatorRects(const ob::SpaceInformationPtr si, double vmax, double rho, RowMatrixXdRef_const &rects, Vector2dRef_const map_lb, Vector2dRef_const map_ub) : DubinsMotionValidator(si, vmax, rho), rects(rects), map_lb(map_lb), map_ub(map_ub) {
    }

    bool collision_free(double x, double y, double theta, double turn_dir, double turn_dist, double next_x, double next_y) const {
      // Line segment
      for (int rect_idx = 0; rect_idx < rects.rows(); ++rect_idx) {
        if (turn_dir == 0 && line_segment_intersects_rect(x, y, next_x, next_y, 
                                                          rects(rect_idx, 0), rects(rect_idx, 1),
                                                          rects(rect_idx, 2), rects(rect_idx, 3))) {
          return false;
        } else if (turn_dir == 1) {
          if (arc_intersects_rect(x, y, theta, turn_dir, turn_dist/vmax, rho, 
                                  next_x, next_y,
                                  rects(rect_idx, 0), rects(rect_idx, 1),
                                  rects(rect_idx, 2), rects(rect_idx, 3))) {
            return false;
          }
        }
      }

      // Check intersection with map boundary
      if (turn_dir == 0 && line_segment_intersects_rect(x, y, next_x, next_y, 
                                                        map_lb(0), map_lb(1),
                                                        map_ub(2), map_ub(3))) {
        return false;
      } else if (turn_dir == 1) {
        if (arc_intersects_rect(x, y, theta, turn_dir, turn_dist/vmax, rho, 
                                next_x, next_y,
                                map_lb(0), map_lb(1),
                                map_ub(2), map_ub(3))) {
          return false;
        }
      }

      return true;
    }

    // If we set validMotion = false, stopState isn't used so we don't have to populate it. Same deal with turns
    RowMatrixXd checkMotionForward(const ob::State *s1, const ob::State *s2, double maxDuration, ob::State *stopState, bool &reach, bool &validMotion) const override {
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
      auto micros = std::chrono::duration_cast<std::chrono::microseconds>(timer_stop - timer_start).count();
      path_elongation_time += ((double)micros)/1e6;
      if (std::isinf(turns(0, 0))) {
        reach = false;
        validMotion = false;
        return turns;
      }

      timer_start = std::chrono::high_resolution_clock::now();

      double x = x1;
      double y = y1;
      double theta = theta1;
      double t = t1;

      double next_x;
      double next_y;
      double next_theta;
      double next_t;

      double valid_path_length = 0.;

      for (int turn_idx = 0; turn_idx < turns.rows(); ++turn_idx) {
        double turn_dir = turns(turn_idx, 0);
        double turn_dist = turns(turn_idx, 1);
        bool stop_on_this_turn = valid_path_length + turn_dist >= maxDuration*vmax;
        if (stop_on_this_turn) {
          turn_dist = maxDuration*vmax - valid_path_length;
        }
        double rho_times_turn_dir = rho*turn_dir;

        double ctheta = cos(theta);
        double stheta = sin(theta);

        if (turn_dir == 0) {
          // S segment
          next_theta = theta;
          next_x = x + turn_dist*ctheta;
          next_y = y + turn_dist*stheta;
        } else {
          // C segment
          next_theta = theta + turn_dist/rho_times_turn_dir;
          next_x = x + rho_times_turn_dir*(-stheta + sin(next_theta));
          next_y = y + rho_times_turn_dir*(ctheta - cos(next_theta));
        }

        next_t = t + turn_dist/vmax;

        if (!collision_free(x, y, theta, turn_dir, turn_dist, next_x, next_y)) {
          reach = false;
          validMotion = false;
          turns(turn_idx, 1) = turn_dist; // Not needed, I think

          timer_stop = std::chrono::high_resolution_clock::now();
          micros = std::chrono::duration_cast<std::chrono::microseconds>(timer_stop - timer_start).count();
          collision_check_time += ((double)micros)/1e6;
          return turns.topRows(turn_idx + 1);
        }

        valid_path_length += turn_dist;

        if (stop_on_this_turn) {
          reach = false;
          stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setX(next_x);
          stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setY(next_y);
          stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setYaw(next_theta);
          stopState->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position = next_t;
          validMotion = true;
          turns(turn_idx, 1) = turn_dist;

          timer_stop = std::chrono::high_resolution_clock::now();
          micros = std::chrono::duration_cast<std::chrono::microseconds>(timer_stop - timer_start).count();
          collision_check_time += ((double)micros)/1e6;

          return turns.topRows(turn_idx + 1);
        }

        x = next_x;
        y = next_y;
        theta = next_theta;
        t = next_t;
      }

      reach = true;
      stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setX(next_x);
      stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setY(next_y);
      stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setYaw(next_theta);
      stopState->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position = next_t;
      validMotion = true;

      timer_stop = std::chrono::high_resolution_clock::now();
      micros = std::chrono::duration_cast<std::chrono::microseconds>(timer_stop - timer_start).count();
      collision_check_time += ((double)micros)/1e6;

      return turns;
    }

    RowMatrixXd checkMotionBackward(const ob::State *s1, const ob::State *s2, double maxDuration, ob::State *stopState, bool &reach, bool &validMotion) const override {
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
      auto micros = std::chrono::duration_cast<std::chrono::microseconds>(timer_stop - timer_start).count();
      if (std::isinf(turns(0, 0))) {
        reach = false;
        validMotion = false;
        return turns;
      }

      timer_start = std::chrono::high_resolution_clock::now();

      double x = x2;
      double y = y2;
      double theta = theta2;
      double t = t2;

      double next_x;
      double next_y;
      double next_theta;
      double next_t;

      double valid_path_length = 0.;

      // std::cout << "checking backwards truncated" << std::endl;
      for (int turn_idx = turns.rows() - 1; turn_idx >= 0; --turn_idx) {
        double turn_dir = turns(turn_idx, 0);
        double turn_dist = turns(turn_idx, 1);
        bool stop_on_this_turn = valid_path_length + turn_dist >= maxDuration*vmax;
        if (stop_on_this_turn) {
          turn_dist = maxDuration*vmax - valid_path_length;
        }
        double rho_times_turn_dir = rho*turn_dir;

        double ctheta = cos(theta);
        double stheta = sin(theta);

        next_t = t - turn_dist/vmax;
        if (turn_dir == 0) {
          // S segment
          next_theta = theta;
          next_x = x - turn_dist*ctheta;
          next_y = y - turn_dist*stheta;
        } else {
          // C segment
          next_theta = theta - turn_dist/rho_times_turn_dir;
          next_x = x - rho_times_turn_dir*(stheta - sin(next_theta));
          next_y = y - rho_times_turn_dir*(-ctheta + cos(next_theta));
        }

        // Use next_x etc because it's checkMotionBackward
        if (!collision_free(next_x, next_y, next_theta, turn_dir, turn_dist, x, y)) {
          reach = false;
          validMotion = false;
          turns(turn_idx, 1) = turn_dist; // Not needed, I think

          timer_stop = std::chrono::high_resolution_clock::now();
          micros = std::chrono::duration_cast<std::chrono::microseconds>(timer_stop - timer_start).count();
          collision_check_time += ((double)micros)/1e6;

          return turns.bottomRows(turns.rows() - turn_idx);
        }

        valid_path_length += turn_dist;

        if (stop_on_this_turn) {
          reach = false;
          stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setX(next_x);
          stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setY(next_y);
          stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setYaw(next_theta);
          stopState->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position = next_t;
          validMotion = true;
          turns(turn_idx, 1) = turn_dist;

          timer_stop = std::chrono::high_resolution_clock::now();
          micros = std::chrono::duration_cast<std::chrono::microseconds>(timer_stop - timer_start).count();
          collision_check_time += ((double)micros)/1e6;

          return turns.bottomRows(turns.rows() - turn_idx);
        }
        x = next_x;
        y = next_y;
        theta = next_theta;
        t = next_t;
      }

      reach = true;
      stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setX(next_x);
      stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setY(next_y);
      stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setYaw(next_theta);
      stopState->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position = next_t;
      validMotion = true;

      timer_stop = std::chrono::high_resolution_clock::now();
      micros = std::chrono::duration_cast<std::chrono::microseconds>(timer_stop - timer_start).count();
      collision_check_time += ((double)micros)/1e6;

      return turns;
    }

  protected:
    RowMatrixXd rects;
    Vector2d map_lb;
    Vector2d map_ub;
};
