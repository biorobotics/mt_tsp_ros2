#pragma once
#include <ompl/base/SpaceInformation.h>
#include <ompl/base/MotionValidator.h>
#include "mt_tsp_ros2/elongate_dubins_path.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/angle_mod.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_motion_validator.h"
#include <chrono>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Exact_circular_kernel_2.h>
#include <CGAL/AABB_traits_2.h>
#include <CGAL/AABB_segment_primitive_2.h>

typedef CGAL::Exact_circular_kernel_2 Circular_k;
typedef CGAL::Circular_arc_2<Circular_k> Circular_arc_2;
typedef CGAL::Circle_2<Circular_k> Circle_2;
typedef Circular_k::Point_2 Point_Circular_k;
typedef Circular_k::FT FT_Circular_k;

typedef CGAL::Simple_cartesian<double> K;

typedef K::Point_2 Point;
typedef K::Segment_2 Segment;

typedef std::list<Segment>::iterator Iterator;
typedef CGAL::AABB_segment_primitive_2<K, Iterator> Primitive;
typedef CGAL::AABB_traits_2<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Primitive_id Primitive_id;

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
    return onSegment(x1i, y1i, x1f, y1f, x2f, y2f) ||
           onSegment(x1i, y1i, x1f, y1f, x2i, y2i) ||
           onSegment(x2i, y2i, x2f, y2f, x1i, y1i) ||
           onSegment(x2i, y2i, x2f, y2f, x1f, y1f);
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
        return ret;
      } else {
        // Poke
        ret(0, 0) = x1 + d(0)*t1;
        ret(0, 1) = y1 + d(1)*t1;
        return ret;
      }
    }

    // here t1 didn't intersect so we either started
    // inside the sphere or are completely past it
    if( t2 >= 0 && t2 <= 1 )
    {
      // ExitWound
      ret(0, 0) = x1 + d(0)*t2;
      ret(0, 1) = y1 + d(1)*t2;
      return ret;
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
  double phiai = atan2(-perpai(1), -perpai(0));

  double phi1 = atan2(intersections(0, 1) - centerai(1), intersections(0, 0) - centerai(0));

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

  double phi2 = atan2(intersections(1, 1) - centerai(1), intersections(1, 0) - centerai(0));

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
    explicit DubinsMotionValidatorRects(const ob::SpaceInformationPtr si, double vmax, double rho, RowMatrixXdRef_const &rects, Vector2dRef_const map_lb, Vector2dRef_const map_ub, const Ref<const VectorXd> &start, const Ref<const VectorXd> &goal) : DubinsMotionValidator(si, vmax, rho), rects(rects), map_lb(map_lb), map_ub(map_ub), start(start), goal(goal) {
      for (int row = 0; row < rects.rows(); ++row) {
        Point point1(rects(row, 0), rects(row, 1)); // xlow, ylow
        Point point2(rects(row, 0), rects(row, 3)); // xlow, yhigh
        Point point3(rects(row, 2), rects(row, 3)); // xhigh, yhigh
        Point point4(rects(row, 2), rects(row, 1)); // xhigh, ylow
        segments.push_back(Segment(point1, point2));
        segments.push_back(Segment(point2, point3));
        segments.push_back(Segment(point3, point4));
        segments.push_back(Segment(point4, point1));
      }
      // Add the map boundaries
      Point point1(map_lb(0), map_lb(1)); // xlow, ylow
      Point point2(map_lb(0), map_ub(1)); // xlow, yhigh
      Point point3(map_ub(0), map_ub(1)); // xhigh, yhigh
      Point point4(map_ub(0), map_lb(1)); // xhigh, ylow
      segments.push_back(Segment(point1, point2));
      segments.push_back(Segment(point2, point3));
      segments.push_back(Segment(point3, point4));
      segments.push_back(Segment(point4, point1));

      aabb_tree = Tree(segments.begin(), segments.end()); 
    }

    bool collision_free(double x, double y, double theta, double turn_dir, double turn_dist, double next_x, double next_y) const {
      /*
      // Line segment
      for (int rect_idx = 0; rect_idx < rects.rows(); ++rect_idx) {
        if (turn_dir == 0 && line_segment_intersects_rect(x, y, next_x, next_y, 
                                                          rects(rect_idx, 0), rects(rect_idx, 1),
                                                          rects(rect_idx, 2), rects(rect_idx, 3))) {
          return false;
        } else if (turn_dir != 0) {
          if (arc_intersects_rect(x, y, theta, turn_dir, turn_dist/rho, rho, 
                                  next_x, next_y,
                                  rects(rect_idx, 0), rects(rect_idx, 1),
                                  rects(rect_idx, 2), rects(rect_idx, 3))) {
            return false;
          }
        }
      }

      // Check if we exit the map
      if (!point_in_rect(x, y, 
                         map_lb(0), map_lb(1),
                         map_ub(0), map_ub(1)) ||
          !point_in_rect(next_x, next_y, 
                         map_lb(0), map_lb(1),
                         map_ub(0), map_ub(1))) {
        return false;
      }

      if (turn_dir != 0) {
        double xir = map_lb(0);
        double xfr = map_ub(0);
        double yir = map_lb(1);
        double yfr = map_ub(1);
        if (arc_intersects_line_segment(x, y, theta, turn_dir, turn_dist/rho, rho,
                                        xir, yir, xir, yfr) ||
            arc_intersects_line_segment(x, y, theta, turn_dir, turn_dist/rho, rho,
                                        xir, yfr, xfr, yfr) ||
            arc_intersects_line_segment(x, y, theta, turn_dir, turn_dist/rho, rho,
                                        xfr, yfr, xfr, yir) ||
            arc_intersects_line_segment(x, y, theta, turn_dir, turn_dist/rho, rho,
                                        xfr, yir, xir, yir)) {
          return false;
        }
      }

      return true;
      */

      // S segment
      if (turn_dir == 0) {
        Point point1(x, y);
        Point point2(next_x, next_y);
        Segment segment(point1, point2);
        return !aabb_tree.do_intersect(segment);
      }

      // C segment

      double c = cos(theta);
      double s = sin(theta);

      Vector2d dir(c, s);
      Vector2d perp(-s*turn_dir, c*turn_dir);
      Vector2d center = Vector2d(x, y) + perp*rho;

      Point_Circular_k center_CGAL(center(0), center(1));
      Circle_2 circle_CGAL(center_CGAL, rho*rho);

      std::list<Primitive_id> primitives;
      aabb_tree.all_intersected_primitives(circle_CGAL.bbox(), std::back_inserter(primitives));

      for (Primitive_id primitive : primitives) {
        if (arc_intersects_line_segment(x, y, theta, turn_dir, turn_dist/rho, rho,
                                        primitive->source().x(), primitive->source().y(), primitive->target().x(), primitive->target().y())) {
          return false;
        }
      }
      return true;

      // Below code doesn't work. I think it's becauseCGAL always assumes the arc goes counterclockwise from point1 to point3
      /*
      if (turn_dist >= 2*M_PI*rho) {
        Vector2d dir(c, s);
        Vector2d perp(-s*turn_dir, c*turn_dir);
        Vector2d center = Vector2d(x, y) + perp*rho;

        Point_Circular_k center_CGAL(center(0), center(1));
        Circle_2 circle_CGAL(center_CGAL, rho*rho);

        std::list<Primitive_id> primitives;
        aabb_tree.all_intersected_primitives(circle_CGAL.bbox(), std::back_inserter(primitives));
        for (Primitive_id primitive : primitives) {
          Matrix2d intersections = line_segment_intersects_circle(primitive->source().x(), primitive->source().y(), primitive->target().x(), primitive->target().y(), center(0), center(1), rho);
          if (std::isfinite(intersections(0, 0))) {
            return false;
          }
        }
        return true;
      } else {
        double rho_times_turn_dir = rho*turn_dir;
        double theta_mid = theta + 0.5*turn_dist/rho_times_turn_dir;
        double x_mid = x + rho_times_turn_dir*(-s + sin(theta_mid));
        double y_mid = y + rho_times_turn_dir*(c - cos(theta_mid));

        Point_Circular_k point1_CGAL(x, y);
        Point_Circular_k point2_CGAL(x_mid, y_mid);
        Point_Circular_k point3_CGAL(next_x, next_y);

        Circular_arc_2 arc_CGAL(point1_CGAL, point2_CGAL, point3_CGAL);

        std::list<Primitive_id> primitives;
        aabb_tree.all_intersected_primitives(arc_CGAL.bbox(), std::back_inserter(primitives));
        for (Primitive_id primitive : primitives) {
          if (arc_intersects_line_segment(x, y, theta, turn_dir, turn_dist/rho, rho,
                                          primitive->source().x(), primitive->source().y(), primitive->target().x(), primitive->target().y())) {
            return false;
          }
        }
        return true;
      }
      */
    }

    // If we set validMotion = false, stopState isn't used so we don't have to populate it. Same deal with turns
    RowMatrixXd checkMotionForward(const ob::State *s1, const ob::State *s2, double maxDuration, ob::State *stopState, bool &reach, bool &validMotion, int num_checks = 1000) const override {
      return checkMotionForward_internal(s1, s2, maxDuration, stopState, reach, validMotion);
      /*
      RowMatrixXd turns1 = checkMotionForward_internal(s1, s2, maxDuration, stopState, reach, validMotion);
      bool reach2;
      bool valid2;
      RowMatrixXd turns2 = DubinsMotionValidator::checkMotionForward(s1, s2, maxDuration, stopState, reach2, valid2);
      if (validMotion != valid2) {
        std::cout << "mismatch, using finer collision check resolution" << std::endl;
        DubinsMotionValidator::checkMotionForward(s1, s2, maxDuration, stopState, reach2, valid2, 100000);
        if (validMotion != valid2) {
          throw std::runtime_error("Mismatch on valid forward");
        }
      }
      return turns1;
      */
    }

    RowMatrixXd checkMotionForward_internal(const ob::State *s1, const ob::State *s2, double maxDuration, ob::State *stopState, bool &reach, bool &validMotion) const {
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
        // std::cout << "failed on initial forward elongation check" << std::endl;
        return turns;
      }

      // Check if the state reached after maxDuration can get to the goal (assuming no obstacles)
      double x = x1;
      double y = y1;
      double theta = theta1;
      double t = t1;

      double next_x = x;
      double next_y = y;
      double next_theta = theta;
      double next_t = t;

      double valid_path_length = 0.;

      /*
      for (int turn_idx = 0; turn_idx < turns.rows(); ++turn_idx) {
        double turn_dir = turns(turn_idx, 0);
        double turn_dist = turns(turn_idx, 1);
        if (turn_dist == 0) {
          continue;
        }
        bool stop_on_this_turn = t2 - t1 > maxDuration + 1e-2 && valid_path_length + turn_dist >= maxDuration*vmax;
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

        valid_path_length += turn_dist;

        if (stop_on_this_turn) {
          if (!check_elongation_possible(next_x, next_y, next_theta, goal(0), goal(1), goal(2), vmax*(goal(3) - next_t), rho)) {
            reach = false;
            validMotion = false;
            // std::cout << "failed on forward elongation check" << std::endl;
            return turns;
          }
          break;
        }

        x = next_x;
        y = next_y;
        theta = next_theta;
        t = next_t;
      }
      */

      // Now perform collision-checks
      timer_start = std::chrono::high_resolution_clock::now();

      x = x1;
      y = y1;
      theta = theta1;
      t = t1;

      next_x = x;
      next_y = y;
      next_theta = theta;
      next_t = t;

      valid_path_length = 0.;

      for (int turn_idx = 0; turn_idx < turns.rows(); ++turn_idx) {
        double turn_dir = turns(turn_idx, 0);
        double turn_dist = turns(turn_idx, 1);
        if (turn_dist == 0) {
          continue;
        }
        bool stop_on_this_turn = t2 - t1 > maxDuration + 1e-2 && valid_path_length + turn_dist >= maxDuration*vmax;
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
          // std::cout << "forward collision check failed" << std::endl;
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

    RowMatrixXd checkMotionBackward(const ob::State *s1, const ob::State *s2, double maxDuration, ob::State *stopState, bool &reach, bool &validMotion, int num_checks = 1000) const override {
      return checkMotionBackward_internal(s1, s2, maxDuration, stopState, reach, validMotion);
      /*
      RowMatrixXd turns1 = checkMotionBackward_internal(s1, s2, maxDuration, stopState, reach, validMotion);
      bool reach2;
      bool valid2;
      RowMatrixXd turns2 = DubinsMotionValidator::checkMotionBackward(s1, s2, maxDuration, stopState, reach2, valid2);
      if (validMotion != valid2) {
        DubinsMotionValidator::checkMotionBackward(s1, s2, maxDuration, stopState, reach2, valid2, 100000);
        if (validMotion != valid2) {
          std::cout << validMotion << " " << valid2 << std::endl;
          throw std::runtime_error("Mismatch on valid backward");
        }
      }
      return turns1;
      */
    }

    RowMatrixXd checkMotionBackward_internal(const ob::State *s1, const ob::State *s2, double maxDuration, ob::State *stopState, bool &reach, bool &validMotion) const {
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
        // std::cout << "failed on initial backward elongation check" << std::endl;
        return turns;
      }

      // Check if the state reached after maxDuration can get to the start (assuming no obstacles)
      double x = x2;
      double y = y2;
      double theta = theta2;
      double t = t2;

      double next_x = x;
      double next_y = y;
      double next_theta = theta;
      double next_t = t;

      double valid_path_length = 0.;

      /*
      for (int turn_idx = turns.rows() - 1; turn_idx >= 0; --turn_idx) {
        double turn_dir = turns(turn_idx, 0);
        double turn_dist = turns(turn_idx, 1);
        if (turn_dist == 0) {
          continue;
        }
        bool stop_on_this_turn = t2 - t1 > maxDuration + 1e-2 && valid_path_length + turn_dist >= maxDuration*vmax;
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

        valid_path_length += turn_dist;

        if (stop_on_this_turn) {
          if (!check_elongation_possible(start(0), start(1), start(2), next_x, next_y, next_theta, vmax*(next_t - start(3)), rho)) {
            reach = false;
            validMotion = false;
            // std::cout << "failed on backward elongation check" << std::endl;
            return turns;
          }
        }
        x = next_x;
        y = next_y;
        theta = next_theta;
        t = next_t;
      }

      timer_start = std::chrono::high_resolution_clock::now();

      x = x2;
      y = y2;
      theta = theta2;
      t = t2;

      next_x = x;
      next_y = y;
      next_theta = theta;
      next_t = t;

      valid_path_length = 0.;
      */

      // std::cout << "checking backwards truncated" << std::endl;
      for (int turn_idx = turns.rows() - 1; turn_idx >= 0; --turn_idx) {
        double turn_dir = turns(turn_idx, 0);
        double turn_dist = turns(turn_idx, 1);
        if (turn_dist == 0) {
          continue;
        }
        bool stop_on_this_turn = t2 - t1 > maxDuration + 1e-2 && valid_path_length + turn_dist >= maxDuration*vmax;
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

          // std::cout << "backward collision check failed" << std::endl;
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
    VectorXd start;
    VectorXd goal;

    // Need segments to persist in memory while the tree is in use
    std::list<Segment> segments;
    Tree aabb_tree;
};
