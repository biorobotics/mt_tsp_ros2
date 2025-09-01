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

using namespace Eigen;

// https://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
double dist_point_to_line_segment(double x1, double y1, double x2, double y2, double x3, double y3); // x3,y3 is the point

// Check if (x3, y3) is on the line segment from (x1, y1) to (x2, y2). Points don't have to be collinear
bool onSegment(double x1, double y1, double x2, double y2, double x3, double y3);

double lineSegmentsIntersect(double x1i, double y1i, double x1f, double y1f,
                             double x2i, double y2i, double x2f, double y2f);

bool point_in_rect(double x, double y, \
                   double xir, double yir, double xfr, double yfr);

bool line_segment_intersects_rect(double x1i, double y1i, double x1f, double y1f,
                                  double xir, double yir,
                                  double xfr, double yfr);

// https://stackoverflow.com/questions/1073336/circle-line-segment-collision-detection-algorithm
Matrix2d line_segment_intersects_circle(double x1, double y1, double x2, double y2, double cx, double cy, double r);


bool arc_intersects_line_segment(double xai, double yai, double thetaai, double turn_dir, double turn_angle_magnitude, double r,
                                 double x1, double y1,
                                 double x2, double y2);

bool arc_intersects_line_segment_get_intersection_point(double xai, double yai, double thetaai, double turn_dir, double turn_angle_magnitude, double r,
                                                        double x1, double y1,
                                                        double x2, double y2,
                                                        double &x_collision, double &y_collision, double &theta_collision, double &collision_dist);

bool arc_intersects_rect(double xai, double yai, double thetaai, double turn_dir, double turn_angle_magnitude, double r,
                         double xaf, double yaf,
                         double xir, double yir,
                         double xfr, double yfr);
