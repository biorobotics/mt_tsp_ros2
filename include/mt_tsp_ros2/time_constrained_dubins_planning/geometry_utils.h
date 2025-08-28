#pragma once
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

typedef std::optional< Tree::Intersection_and_primitive_id<Segment>::Type > Segment_intersection;

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

bool arc_intersects_line_segment_get_intersection_point(double xai, double yai, double thetaai, double turn_dir, double turn_angle_magnitude, double r,
                                                        double x1, double y1,
                                                        double x2, double y2,
                                                        double &x_collision, double &y_collision, double &theta_collision, double &collision_dist) {
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
    x_collision = intersections(0, 0);
    y_collision = intersections(0, 1);
    theta_collision = thetaai + diff1i;
    collision_dist = r*abs_diff1i;
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
    x_collision = intersections(1, 0);
    y_collision = intersections(1, 1);
    theta_collision = thetaai + diff2i;
    collision_dist = r*abs_diff2i;
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
