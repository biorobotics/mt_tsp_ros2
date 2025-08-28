#pragma once
#include <ompl/base/SpaceInformation.h>
#include <ompl/base/MotionValidator.h>
#include "mt_tsp_ros2/elongate_dubins_path.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/angle_mod.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/no_time_dubins_motion_validator.h"
#include <chrono>

#include "mt_tsp_ros2/time_constrained_dubins_planning/geometry_utils.h"

namespace ob = ompl::base;

class NoTimeDubinsMotionValidatorRects : public NoTimeDubinsMotionValidator {
  public:
    explicit NoTimeDubinsMotionValidatorRects(const ob::SpaceInformationPtr si, double rho, RowMatrixXdRef_const &rects, Vector2dRef_const map_lb, Vector2dRef_const map_ub) : NoTimeDubinsMotionValidator(si, rho), rects(rects), map_lb(map_lb), map_ub(map_ub) {
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

      rng_ = ompl::RNG(1);
    }

    const Tree &get_aabb_tree() {
      return aabb_tree;
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

      CGAL::Bbox_2 bbox(center(0) - rho, center(1) - rho, center(0) + rho, center(1) + rho);
      std::list<Primitive_id> primitives;
      aabb_tree.all_intersected_primitives(bbox, std::back_inserter(primitives));

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

    bool collision_free_get_intersection_point(double x, double y, double theta, double turn_dir, double turn_dist, double next_x, double next_y, double &x_collision, double &y_collision, double &theta_collision, double &collision_dist) const {
      // S segment
      if (turn_dir == 0) {
        Point point1(x, y);
        Point point2(next_x, next_y);
        Segment segment(point1, point2);
        Segment_intersection intersection = aabb_tree.any_intersection(segment);
        if (!intersection) {
          return true;
        }
        const Point* p = std::get_if<Point>(&(intersection->first));
        if (!p) {
          throw std::runtime_error("Intersection is not a point");
        }
        x_collision = p->x();
        y_collision = p->y();
        theta_collision = theta;
        double delta_x = x_collision - x;
        double delta_y = y_collision - y;
        collision_dist = sqrt(delta_x*delta_x + delta_y*delta_y);
        return false;
      }

      // C segment
      double c = cos(theta);
      double s = sin(theta);

      Vector2d dir(c, s);
      Vector2d perp(-s*turn_dir, c*turn_dir);
      Vector2d center = Vector2d(x, y) + perp*rho;

      CGAL::Bbox_2 bbox(center(0) - rho, center(1) - rho, center(0) + rho, center(1) + rho);
      std::list<Primitive_id> primitives;
      aabb_tree.all_intersected_primitives(bbox, std::back_inserter(primitives));

      for (Primitive_id primitive : primitives) {
        if (arc_intersects_line_segment_get_intersection_point(x, y, theta, turn_dir, turn_dist/rho, rho,
                                                               primitive->source().x(), primitive->source().y(), primitive->target().x(), primitive->target().y(),
                                                               x_collision, y_collision, theta_collision, collision_dist)) {
          return false;
        }
      }
      return true;
    }

    bool checkMotion(const ob::State *s1, const ob::State *s2, std::pair<ob::State*, double> &lastValid) const override {
      double x1 = s1->as<ob::SE2StateSpace::StateType>()->getX();
      double y1 = s1->as<ob::SE2StateSpace::StateType>()->getY();
      double theta1 = s1->as<ob::SE2StateSpace::StateType>()->getYaw();

      double x2 = s2->as<ob::SE2StateSpace::StateType>()->getX();
      double y2 = s2->as<ob::SE2StateSpace::StateType>()->getY();
      double theta2 = s2->as<ob::SE2StateSpace::StateType>()->getYaw();

      RowMatrixXd turns = turns_for_dubins_path(x1, y1, theta1, x2, y2, theta2, rho);

      double ctheta1 = cos(theta1);
      double stheta1 = sin(theta1);

      double x = x1;
      double y = y1;
      double theta = theta1;

      double ctheta = ctheta1;
      double stheta = stheta1;

      double next_x = x;
      double next_y = y;
      double next_theta = theta;

      double next_ctheta = ctheta;
      double next_stheta = stheta;

      // Perform collision-checks
      for (int turn_idx = 0; turn_idx < turns.rows(); ++turn_idx) {
        double turn_dir = turns(turn_idx, 0);
        double turn_dist = turns(turn_idx, 1);
        if (turn_dist == 0) {
          continue;
        }
        double rho_times_turn_dir = rho*turn_dir;

        if (turn_dir == 0) {
          // S segment
          next_theta = theta;
          next_x = x + turn_dist*ctheta;
          next_y = y + turn_dist*stheta;
        } else {
          // C segment
          next_theta = theta + turn_dist/rho_times_turn_dir;
          next_ctheta = cos(next_theta);
          next_stheta = sin(next_theta);
          next_x = x + rho_times_turn_dir*(-stheta + next_stheta);
          next_y = y + rho_times_turn_dir*(ctheta - next_ctheta);
        }

        if (!collision_free(x, y, theta, turn_dir, turn_dist, next_x, next_y)) {
          return false;
        }

        x = next_x;
        y = next_y;
        theta = next_theta;

        ctheta = next_ctheta;
        stheta = next_stheta;
      }

      return true;
    }

  protected:
    RowMatrixXd rects;
    Vector2d map_lb;
    Vector2d map_ub;

    // Need segments to persist in memory while the tree is in use
    std::list<Segment> segments;
    Tree aabb_tree;
    mutable ompl::RNG rng_;
};
