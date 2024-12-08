#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include "gmdm.h"
#include "mt_tsp_ros2/dubins.h"
#include <stdexcept>

using namespace Eigen;
typedef Ref<Matrix<double, Dynamic, 1>> VectorXdRef;
typedef Matrix<double, Dynamic, Dynamic, RowMajor> RowMatrixXd;

typedef const Ref<const Vector2d>& Vector2dRef_const;

double root_find(double rho, double s, double l_m) {
  double theta = M_PI/4;
  double r = 4*rho*theta - (s - l_m + 4*rho*sin(theta));
  int max_iter = 10;
  for (int i = 0; i < max_iter; ++i) {
    if (abs(r) < 1e-4) {
      break;
    }
    double dr_dtheta = 4*rho - 4*rho*cos(theta);
    theta -= r/dr_dtheta;
    r = 4*rho*theta - (s - l_m + 4*rho*sin(theta));
  }
  return theta;
}

double arclength(const Vector2d &v1, const Vector2d &v2, bool left, double rho) {
  double theta = atan2(v2(1), v2(0)) - atan2(v1(1), v1(0));
  if (theta < 0 && left) {
    theta += 2*M_PI;
  } else if (theta > 0 && !left) {
    theta -= 2*M_PI;
  }
  return abs(theta*rho);
}

Vector2d line_segments_intersect_point(double x1i, double y1i, double x1f, double y1f,
                                       double x2i, double y2i, double x2f, double y2f) {
  double a = x1f - x1i;
  double b = x2i - x2f;
  double c = y1f - y1i;
  double d = y2i - y2f;

  double det = a*d - b*c;

  if (det == 0) {
    throw std::runtime_error("Line segments are parallel");
  }

  double lx = x2i - x1i;
  double ly = y2i - y1i;

  double t1 = (d*lx - b*ly)/det;
  double t2 = (-c*lx + a*ly)/det;
  if (!(0 <= t1 and t1 <= 1 and 0 <= t2 and t2 <= 1)) {
    return std::numeric_limits<double>::infinity()*Vector2d::Ones();
  }
  return Vector2d(x1i + t1*(x1f - x1i), y1i + t1*(y1f - y1i));
}

RowMatrixXd turns_for_dubins_path(double x_0, double y_0, double theta_0, double x_f, double y_f, double theta_f, double rho) {
  DubinsPath path;
  double q_0[3] = {x_0, y_0, theta_0};
  double q_f[3] = {x_f, y_f, theta_f};
  dubins_shortest_path(&path, q_0, q_f, rho);

  RowMatrixXd turns = RowMatrixXd::Zero(3, 2);
  if (path.type == DubinsPathType::LSL) {
    turns(0, 0) = 1.;
    turns(1, 0) = 0.;
    turns(2, 0) = 1.;
  } else if (path.type == DubinsPathType::LSR) {
    turns(0, 0) = 1.;
    turns(1, 0) = 0.;
    turns(2, 0) = -1.;
  } else if (path.type == DubinsPathType::RSR) {
    turns(0, 0) = -1.;
    turns(1, 0) = 0.;
    turns(2, 0) = -1.;
  } else if (path.type == DubinsPathType::RSL) {
    turns(0, 0) = -1.;
    turns(1, 0) = 0.;
    turns(2, 0) = 1.;
  } else if (path.type == DubinsPathType::RLR) {
    turns(0, 0) = -1.;
    turns(1, 0) = 1.;
    turns(2, 0) = -1.;
  } else if (path.type == DubinsPathType::LRL) {
    turns(0, 0) = 1.;
    turns(1, 0) = -1.;
    turns(2, 0) = 1.;
  }
  turns(0, 1) = path.param[0];
  turns(1, 1) = path.param[1];
  turns(2, 1) = path.param[2];

  return turns;
}

RowMatrixXd path_with_waypoint(double x_0, double y_0, double theta_0, double x_wp, double y_wp, double theta_wp, double x_f, double y_f, double theta_f, double rho) {
  RowMatrixXd turns(6, 2);
  turns.topRows(3) = turns_for_dubins_path(x_0, y_0, theta_0, x_wp, y_wp, theta_wp, rho);
  turns.bottomRows(3) = turns_for_dubins_path(x_wp, y_wp, theta_wp, x_f, y_f, theta_f, rho);
  return turns;
}

double path_length_with_waypoint(double x_0, double y_0, double theta_0, double x_wp, double y_wp, double theta_wp, double x_f, double y_f, double theta_f, double rho) {
  RowMatrixXd turns = path_with_waypoint(x_0, y_0, theta_0, x_wp, y_wp, theta_wp, x_f, y_f, theta_f, rho);
  return turns.col(1).sum();
}

bool check_elongation_possible(double x_0, double y_0, double theta_0, double x_f, double y_f, double theta_f, double s, double rho) {
  // LRL
  double LRL_dist_A = std::numeric_limits<double>::infinity();
  double LRL_dist_B = std::numeric_limits<double>::infinity();

  double c_0 = cos(theta_0);
  double s_0 = sin(theta_0);

  Vector2d dir_0(c_0, s_0);
  Vector2d perp_0(-s_0, c_0);
  Vector2d center1_L = Vector2d(x_0, y_0) + perp_0*rho;

  double c_f = cos(theta_f);
  double s_f = sin(theta_f);
  Vector2d dir_f(c_f, s_f);
  Vector2d perp_f(-s_f, c_f);
  Vector2d center3_L = Vector2d(x_f, y_f) + perp_f*rho;

  Vector2d V = center3_L - center1_L;
  double D = V.norm();
  if (D <= 4*rho) {
    double gamma = atan2(V(1), V(0));
    double theta = acos(D/(4*rho));

    double theta_A = gamma + theta;
    double c_A = cos(theta_A);
    double s_A = sin(theta_A);
    Vector2d vec_A = Vector2d(c_A, s_A);
    Vector2d center2 = center1_L + vec_A*rho*2;
    LRL_dist_A = arclength(-perp_0, vec_A, true, rho);
    LRL_dist_A += arclength(center1_L - center2, center3_L - center2, false, rho);
    LRL_dist_A += arclength(center2 - center3_L, -perp_f, true, rho);

    double theta_B = gamma - theta;
    double c_B = cos(theta_B);
    double s_B = sin(theta_B);
    Vector2d vec_B = Vector2d(c_B, s_B);
    center2 = center1_L + vec_B*rho*2;
    LRL_dist_B = arclength(-perp_0, vec_B, true, rho);
    LRL_dist_B += arclength(center1_L - center2, center3_L - center2, false, rho);
    LRL_dist_B += arclength(center2 - center3_L, -perp_f, true, rho);
  }

  // RLR
  double RLR_dist_A = std::numeric_limits<double>::infinity();
  double RLR_dist_B = std::numeric_limits<double>::infinity();

  perp_0 = Vector2d(s_0, -c_0);
  Vector2d center1_R = Vector2d(x_0, y_0) + perp_0*rho;

  perp_f = Vector2d(s_f, -c_f);
  Vector2d center3_R = Vector2d(x_f, y_f) + perp_f*rho;

  V = center3_R - center1_R;
  D = V.norm();
  if (D <= 4*rho) {
    double gamma = atan2(V(1), V(0));
    double theta = acos(D/(4*rho));

    double theta_A = gamma + theta;
    double c_A = cos(theta_A);
    double s_A = sin(theta_A);
    Vector2d vec_A = Vector2d(c_A, s_A);
    Vector2d center2 = center1_R + vec_A*rho*2;
    RLR_dist_A = arclength(-perp_0, vec_A, false, rho);
    RLR_dist_A += arclength(center1_R - center2, center3_R - center2, true, rho);
    RLR_dist_A += arclength(center2 - center3_R, -perp_f, false, rho);

    double theta_B = gamma - theta;
    double c_B = cos(theta_B);
    double s_B = sin(theta_B);
    Vector2d vec_B = Vector2d(c_B, s_B);
    center2 = center1_R + Vector2d(c_B, s_B)*rho*2;
    RLR_dist_B = arclength(-perp_0, vec_B, false, rho);
    RLR_dist_B += arclength(center1_R - center2, center3_R - center2, true, rho);
    RLR_dist_B += arclength(center2 - center3_R, -perp_f, false, rho);
  }

  double l_LRL_s = std::min(LRL_dist_A, LRL_dist_B);
  double l_RLR_s = std::min(RLR_dist_A, RLR_dist_B);

  double l_LRL_l = std::max(LRL_dist_A, LRL_dist_B);
  double l_RLR_l = std::max(RLR_dist_A, RLR_dist_B);

  double l_m = std::min(l_LRL_s, l_RLR_s);

  // LSL
  double l_LSL = std::numeric_limits<double>::infinity();
  VectorXd path_LSL = csc_inverse(x_0, y_0, theta_0, x_f, y_f, theta_f, rho, 1., rho, rho, 1.);
  if (!std::isinf(path_LSL(0))) {
    l_m = std::min(l_m, path_LSL(3));
    l_LSL = path_LSL(3);
  }

  // LSR
  double l_LSR = std::numeric_limits<double>::infinity();
  VectorXd path_LSR = csc_inverse(x_0, y_0, theta_0, x_f, y_f, theta_f, rho, 1., rho, rho, -1.);
  if (!std::isinf(path_LSR(0))) {
    l_m = std::min(l_m, path_LSR(3));
    l_LSR = path_LSR(3);
  }

  // RSR
  double l_RSR = std::numeric_limits<double>::infinity();
  VectorXd path_RSR = csc_inverse(x_0, y_0, theta_0, x_f, y_f, theta_f, rho, -1., rho, rho, -1.);
  if (!std::isinf(path_RSR(0))) {
    l_m = std::min(l_m, path_RSR(3));
    l_RSR = path_RSR(3);
  }

  // RSL
  double l_RSL = std::numeric_limits<double>::infinity();
  VectorXd path_RSL = csc_inverse(x_0, y_0, theta_0, x_f, y_f, theta_f, rho, -1., rho, rho, 1.);
  if (!std::isinf(path_RSL(0))) {
    l_m = std::min(l_m, path_RSL(3));
    l_RSL = path_RSL(3);
  }

  if (l_m > s) {
    return false; // We can't elongate a Dubins path to length s if the shortest Dubins path has length larger than s
  }

  if (l_m == l_LRL_s || l_m == l_RLR_s) {
    // If shortest path is CCC, we can always elongate it
    return true;
  }

  // Check if path is in O
  if (l_m == l_LSL || l_m == l_LSR || l_m == l_RSR || l_m == l_RSL) {
    double ldist = (center3_L - center1_L).norm();
    double rdist = (center3_R - center1_R).norm();
    // O4 and O5
    if (rdist >= 4*rho || ldist >= 4*rho) {
      return true;
    }
    // O1, O2, O3
    if ((l_m == l_LSL && (path_LSL(0) >= M_PI || path_LSL(2) >= M_PI || path_LSL(1) >= 4)) ||
        (l_m == l_LSR && (path_LSR(0) >= M_PI || path_LSR(2) >= M_PI || path_LSR(1) >= 4)) ||
        (l_m == l_RSR && (path_RSR(0) >= M_PI || path_RSR(2) >= M_PI || path_RSR(1) >= 4)) ||
        (l_m == l_RSL && (path_RSL(0) >= M_PI || path_RSL(2) >= M_PI || path_RSL(1) >= 4))) {
      return true;
    }
  }

  double l1 = std::max(l_LRL_s, l_RLR_s);
  double l2 = l_m + 2*M_PI*rho;
  if (l_LRL_l < l2) {
    l2 = l_LRL_l;
  }
  if (l_RLR_l < l2) {
    l2 = l_RLR_l;
  }
  if (l_RSR != l_m && l_RSR < l2) {
    l2 = l_RSR;
  }
  if (l_RSL != l_m && l_RSL < l2) {
    l2 = l_RSL;
  }
  if (l_LSR != l_m && l_LSR < l2) {
    l2 = l_LSR;
  }
  if (l_LSL != l_m && l_LSL < l2) {
    l2 = l_LSL;
  }

  return s <= l1 || s >= l2;
}

Vector3d get_elongation_intervals(double x_0, double y_0, double theta_0, double x_f, double y_f, double theta_f, double rho) {
  // LRL
  double LRL_dist_A = std::numeric_limits<double>::infinity();
  double LRL_dist_B = std::numeric_limits<double>::infinity();

  double c_0 = cos(theta_0);
  double s_0 = sin(theta_0);

  Vector2d dir_0(c_0, s_0);
  Vector2d perp_0(-s_0, c_0);
  Vector2d center1_L = Vector2d(x_0, y_0) + perp_0*rho;

  double c_f = cos(theta_f);
  double s_f = sin(theta_f);
  Vector2d dir_f(c_f, s_f);
  Vector2d perp_f(-s_f, c_f);
  Vector2d center3_L = Vector2d(x_f, y_f) + perp_f*rho;

  Vector2d V = center3_L - center1_L;
  double D = V.norm();
  if (D <= 4*rho) {
    double gamma = atan2(V(1), V(0));
    double theta = acos(D/(4*rho));

    double theta_A = gamma + theta;
    double c_A = cos(theta_A);
    double s_A = sin(theta_A);
    Vector2d vec_A = Vector2d(c_A, s_A);
    Vector2d center2 = center1_L + vec_A*rho*2;
    LRL_dist_A = arclength(-perp_0, vec_A, true, rho);
    LRL_dist_A += arclength(center1_L - center2, center3_L - center2, false, rho);
    LRL_dist_A += arclength(center2 - center3_L, -perp_f, true, rho);

    double theta_B = gamma - theta;
    double c_B = cos(theta_B);
    double s_B = sin(theta_B);
    Vector2d vec_B = Vector2d(c_B, s_B);
    center2 = center1_L + vec_B*rho*2;
    LRL_dist_B = arclength(-perp_0, vec_B, true, rho);
    LRL_dist_B += arclength(center1_L - center2, center3_L - center2, false, rho);
    LRL_dist_B += arclength(center2 - center3_L, -perp_f, true, rho);
  }

  // RLR
  double RLR_dist_A = std::numeric_limits<double>::infinity();
  double RLR_dist_B = std::numeric_limits<double>::infinity();

  perp_0 = Vector2d(s_0, -c_0);
  Vector2d center1_R = Vector2d(x_0, y_0) + perp_0*rho;

  perp_f = Vector2d(s_f, -c_f);
  Vector2d center3_R = Vector2d(x_f, y_f) + perp_f*rho;

  V = center3_R - center1_R;
  D = V.norm();
  if (D <= 4*rho) {
    double gamma = atan2(V(1), V(0));
    double theta = acos(D/(4*rho));

    double theta_A = gamma + theta;
    double c_A = cos(theta_A);
    double s_A = sin(theta_A);
    Vector2d vec_A = Vector2d(c_A, s_A);
    Vector2d center2 = center1_R + vec_A*rho*2;
    RLR_dist_A = arclength(-perp_0, vec_A, false, rho);
    RLR_dist_A += arclength(center1_R - center2, center3_R - center2, true, rho);
    RLR_dist_A += arclength(center2 - center3_R, -perp_f, false, rho);

    double theta_B = gamma - theta;
    double c_B = cos(theta_B);
    double s_B = sin(theta_B);
    Vector2d vec_B = Vector2d(c_B, s_B);
    center2 = center1_R + Vector2d(c_B, s_B)*rho*2;
    RLR_dist_B = arclength(-perp_0, vec_B, false, rho);
    RLR_dist_B += arclength(center1_R - center2, center3_R - center2, true, rho);
    RLR_dist_B += arclength(center2 - center3_R, -perp_f, false, rho);
  }

  double l_LRL_s = std::min(LRL_dist_A, LRL_dist_B);
  double l_RLR_s = std::min(RLR_dist_A, RLR_dist_B);

  double l_LRL_l = std::max(LRL_dist_A, LRL_dist_B);
  double l_RLR_l = std::max(RLR_dist_A, RLR_dist_B);

  double l_m = std::min(l_LRL_s, l_RLR_s);

  // LSL
  double l_LSL = std::numeric_limits<double>::infinity();
  VectorXd path_LSL = csc_inverse(x_0, y_0, theta_0, x_f, y_f, theta_f, rho, 1., rho, rho, 1.);
  if (!std::isinf(path_LSL(0))) {
    l_m = std::min(l_m, path_LSL(3));
    l_LSL = path_LSL(3);
  }

  // LSR
  double l_LSR = std::numeric_limits<double>::infinity();
  VectorXd path_LSR = csc_inverse(x_0, y_0, theta_0, x_f, y_f, theta_f, rho, 1., rho, rho, -1.);
  if (!std::isinf(path_LSR(0))) {
    l_m = std::min(l_m, path_LSR(3));
    l_LSR = path_LSR(3);
  }

  // RSR
  double l_RSR = std::numeric_limits<double>::infinity();
  VectorXd path_RSR = csc_inverse(x_0, y_0, theta_0, x_f, y_f, theta_f, rho, -1., rho, rho, -1.);
  if (!std::isinf(path_RSR(0))) {
    l_m = std::min(l_m, path_RSR(3));
    l_RSR = path_RSR(3);
  }

  // RSL
  double l_RSL = std::numeric_limits<double>::infinity();
  VectorXd path_RSL = csc_inverse(x_0, y_0, theta_0, x_f, y_f, theta_f, rho, -1., rho, rho, 1.);
  if (!std::isinf(path_RSL(0))) {
    l_m = std::min(l_m, path_RSL(3));
    l_RSL = path_RSL(3);
  }

  if (l_m == l_LRL_s || l_m == l_RLR_s) {
    // If shortest path is CCC, we can always elongate it
    return Vector3d(l_m, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
  }

  // Check if path is in O
  if (l_m == l_LSL || l_m == l_LSR || l_m == l_RSR || l_m == l_RSL) {
    double ldist = (center3_L - center1_L).norm();
    double rdist = (center3_R - center1_R).norm();
    // O4 and O5
    if (rdist >= 4*rho || ldist >= 4*rho) {
      return Vector3d(l_m, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
    }
    // O1, O2, O3
    if ((l_m == l_LSL && (path_LSL(0) >= M_PI || path_LSL(2) >= M_PI || path_LSL(1) >= 4)) ||
        (l_m == l_LSR && (path_LSR(0) >= M_PI || path_LSR(2) >= M_PI || path_LSR(1) >= 4)) ||
        (l_m == l_RSR && (path_RSR(0) >= M_PI || path_RSR(2) >= M_PI || path_RSR(1) >= 4)) ||
        (l_m == l_RSL && (path_RSL(0) >= M_PI || path_RSL(2) >= M_PI || path_RSL(1) >= 4))) {
      return Vector3d(l_m, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
    }
  }

  double l1 = std::max(l_LRL_s, l_RLR_s);
  double l2 = l_m + 2*M_PI*rho;
  if (l_LRL_l < l2) {
    l2 = l_LRL_l;
  }
  if (l_RLR_l < l2) {
    l2 = l_RLR_l;
  }
  if (l_RSR != l_m && l_RSR < l2) {
    l2 = l_RSR;
  }
  if (l_RSL != l_m && l_RSL < l2) {
    l2 = l_RSL;
  }
  if (l_LSR != l_m && l_LSR < l2) {
    l2 = l_LSR;
  }
  if (l_LSL != l_m && l_LSL < l2) {
    l2 = l_LSL;
  }

  return Vector3d(l_m, l1, l2);
}

// Returns a sequence of (turn direction, dist) pairs
RowMatrixXd elongated_dubins_path(double x_0, double y_0, double theta_0, double x_f, double y_f, double theta_f, double s, double rho, bool verbose) {
  // LRL
  double LRL_dist_A = std::numeric_limits<double>::infinity();
  double LRL_dist_B = std::numeric_limits<double>::infinity();

  double c_0 = cos(theta_0);
  double s_0 = sin(theta_0);

  Vector2d dir_0(c_0, s_0);
  Vector2d perp_0(-s_0, c_0);
  Vector2d center1_L = Vector2d(x_0, y_0) + perp_0*rho;

  double c_f = cos(theta_f);
  double s_f = sin(theta_f);
  Vector2d dir_f(c_f, s_f);
  Vector2d perp_f(-s_f, c_f);
  Vector2d center3_L = Vector2d(x_f, y_f) + perp_f*rho;

  Vector2d V = center3_L - center1_L;
  double D = V.norm();

  double LRL_l_dist1;
  double LRL_l_dist2;
  double LRL_l_dist3;
  if (D <= 4*rho) {
    double gamma = atan2(V(1), V(0));
    double theta = acos(D/(4*rho));

    double theta_A = gamma + theta;
    double c_A = cos(theta_A);
    double s_A = sin(theta_A);
    Vector2d vec_A = Vector2d(c_A, s_A);
    Vector2d center2 = center1_L + vec_A*rho*2;
    double dist1_A = arclength(-perp_0, vec_A, true, rho);
    double dist2_A = arclength(center1_L - center2, center3_L - center2, false, rho);
    double dist3_A = arclength(center2 - center3_L, -perp_f, true, rho);
    LRL_dist_A = dist1_A + dist2_A + dist3_A;

    double theta_B = gamma - theta;
    double c_B = cos(theta_B);
    double s_B = sin(theta_B);
    Vector2d vec_B = Vector2d(c_B, s_B);
    center2 = center1_L + vec_B*rho*2;
    double dist1_B = arclength(-perp_0, vec_B, true, rho);
    double dist2_B = arclength(center1_L - center2, center3_L - center2, false, rho);
    double dist3_B = arclength(center2 - center3_L, -perp_f, true, rho);
    LRL_dist_B = dist1_B + dist2_B + dist3_B;

    if (LRL_dist_A > LRL_dist_B) {
      LRL_l_dist1 = dist1_A;
      LRL_l_dist2 = dist2_A;
      LRL_l_dist3 = dist3_A;
    } else {
      LRL_l_dist1 = dist1_B;
      LRL_l_dist2 = dist2_B;
      LRL_l_dist3 = dist3_B;
    }
  }

  // RLR
  double RLR_dist_A = std::numeric_limits<double>::infinity();
  double RLR_dist_B = std::numeric_limits<double>::infinity();

  perp_0 = Vector2d(s_0, -c_0);
  Vector2d center1_R = Vector2d(x_0, y_0) + perp_0*rho;

  perp_f = Vector2d(s_f, -c_f);
  Vector2d center3_R = Vector2d(x_f, y_f) + perp_f*rho;

  V = center3_R - center1_R;
  D = V.norm();

  double RLR_l_dist1;
  double RLR_l_dist2;
  double RLR_l_dist3;
  if (D <= 4*rho) {
    double gamma = atan2(V(1), V(0));
    double theta = acos(D/(4*rho));

    double theta_A = gamma + theta;
    double c_A = cos(theta_A);
    double s_A = sin(theta_A);
    Vector2d vec_A = Vector2d(c_A, s_A);
    Vector2d center2 = center1_R + vec_A*rho*2;
    double dist1_A = arclength(-perp_0, vec_A, false, rho);
    double dist2_A = arclength(center1_R - center2, center3_R - center2, true, rho);
    double dist3_A = arclength(center2 - center3_R, -perp_f, false, rho);
    RLR_dist_A = dist1_A + dist2_A + dist3_A;

    double theta_B = gamma - theta;
    double c_B = cos(theta_B);
    double s_B = sin(theta_B);
    Vector2d vec_B = Vector2d(c_B, s_B);
    center2 = center1_R + Vector2d(c_B, s_B)*rho*2;
    double dist1_B = arclength(-perp_0, vec_B, false, rho);
    double dist2_B = arclength(center1_R - center2, center3_R - center2, true, rho);
    double dist3_B = arclength(center2 - center3_R, -perp_f, false, rho);
    RLR_dist_B = dist1_B + dist2_B + dist3_B;

    if (RLR_dist_A > RLR_dist_B) {
      RLR_l_dist1 = dist1_A;
      RLR_l_dist2 = dist2_A;
      RLR_l_dist3 = dist3_A;
    } else {
      RLR_l_dist1 = dist1_B;
      RLR_l_dist2 = dist2_B;
      RLR_l_dist3 = dist3_B;
    }
  }

  double l_LRL_s = std::min(LRL_dist_A, LRL_dist_B);
  double l_RLR_s = std::min(RLR_dist_A, RLR_dist_B);

  double l_LRL_l = std::max(LRL_dist_A, LRL_dist_B);
  double l_RLR_l = std::max(RLR_dist_A, RLR_dist_B);

  double l_m = std::min(l_LRL_s, l_RLR_s);

  // LSL
  double l_LSL = std::numeric_limits<double>::infinity();
  VectorXd path_LSL = csc_inverse(x_0, y_0, theta_0, x_f, y_f, theta_f, rho, 1., rho, rho, 1.);
  if (!std::isinf(path_LSL(0))) {
    l_m = std::min(l_m, path_LSL(3));
    l_LSL = path_LSL(3);
  }

  // LSR
  double l_LSR = std::numeric_limits<double>::infinity();
  VectorXd path_LSR = csc_inverse(x_0, y_0, theta_0, x_f, y_f, theta_f, rho, 1., rho, rho, -1.);
  if (!std::isinf(path_LSR(0))) {
    l_m = std::min(l_m, path_LSR(3));
    l_LSR = path_LSR(3);
  }

  // RSR
  double l_RSR = std::numeric_limits<double>::infinity();
  VectorXd path_RSR = csc_inverse(x_0, y_0, theta_0, x_f, y_f, theta_f, rho, -1., rho, rho, -1.);
  if (!std::isinf(path_RSR(0))) {
    l_m = std::min(l_m, path_RSR(3));
    l_RSR = path_RSR(3);
  }

  // RSL
  double l_RSL = std::numeric_limits<double>::infinity();
  VectorXd path_RSL = csc_inverse(x_0, y_0, theta_0, x_f, y_f, theta_f, rho, -1., rho, rho, 1.);
  if (!std::isinf(path_RSL(0))) {
    l_m = std::min(l_m, path_RSL(3));
    l_RSL = path_RSL(3);
  }

  if (l_m > s) {
    return std::numeric_limits<double>::infinity()*RowMatrixXd::Ones(1, 2); // We can't elongate a Dubins path to length s if the shortest Dubins path has length larger than s
  }

  // If shortest path is CCC, we can always elongate it
  if (l_m == l_LRL_s) { // Tested
    if (verbose) {
      std::cout << "Elongating LRL path" << std::endl;
    }
    RowMatrixXd turns = RowMatrixXd::Zero(6, 2);

    // LRL
    VectorXd path_LRL = ccc_inverse(x_0, y_0, theta_0, x_f, y_f, theta_f, rho, 1., rho, -1., rho, 1.);
    assert(!std::isinf(path_LRL(0)));

    // Left turn
    turns(0, 0) = 1.;
    turns(0, 1) = path_LRL(0)*rho;

    // Straight (this is a choice I made, we could also turn right for a bit)
    turns(1, 0) = 0.;
    turns(1, 1) = (s - l_m)/2;

    // Right turn
    turns(2, 0) = -1.;
    turns(2, 1) = rho*M_PI;

    // Straight
    turns(3, 0) = 0.;
    turns(3, 1) = (s - l_m)/2;

    // Right turn
    turns(4, 0) = -1.;
    turns(4, 1) = path_LRL(1)*rho - rho*M_PI; // Total distance for the original R segment minus distance we traveled in the first split-up R segment

    // Left turn
    turns(5, 0) = 1.;
    turns(5, 1) = path_LRL(2)*rho;

    return turns;
  }

  if (l_m == l_RLR_s) { // Tested
    if (verbose) {
      std::cout << "Elongating RLR path" << std::endl;
    }
    RowMatrixXd turns = RowMatrixXd::Zero(6, 2);

    // RLR
    VectorXd path_RLR = ccc_inverse(x_0, y_0, theta_0, x_f, y_f, theta_f, rho, -1., rho, 1., rho, -1.);
    assert(!std::isinf(path_RLR(0)));

    // Right turn
    turns(0, 0) = -1.;
    turns(0, 1) = path_RLR(0)*rho;

    // Straight (this is a choice I made, we could also turn left for a bit)
    turns(1, 0) = 0.;
    turns(1, 1) = (s - l_m)/2;

    // Left turn
    turns(2, 0) = 1.;
    turns(2, 1) = rho*M_PI;

    // Straight
    turns(3, 0) = 0.;
    turns(3, 1) = (s - l_m)/2;

    // Left turn
    turns(4, 0) = 1.;
    turns(4, 1) = path_RLR(1)*rho - rho*M_PI; // Total distance for the original R segment minus distance we traveled in the first split-up R segment

    // Right turn
    turns(5, 0) = -1.;
    turns(5, 1) = path_RLR(2)*rho;

    return turns;
  }

  // Check if path is in O
  if (l_m == l_LSL || l_m == l_LSR || l_m == l_RSR || l_m == l_RSL) {
    // O1, O2, O3
    double first_C_sign = 0.;
    double first_C_angle = 0.;
    double second_C_sign = 0.;
    double second_C_angle = 0.;
    double S_dist = 0.;
    if (l_m == l_LSL) {
      first_C_angle = path_LSL(0);
      first_C_sign = 1.;
      second_C_angle = path_LSL(2);
      second_C_sign = 1.;
      S_dist = path_LSL(1)*rho;
    } else if (l_m == l_LSR) {
      first_C_angle = path_LSR(0);
      first_C_sign = 1.;
      second_C_angle = path_LSR(2);
      second_C_sign = -1.;
      S_dist = path_LSR(1)*rho;
    } else if (l_m == l_RSR) {
      first_C_angle = path_RSR(0);
      first_C_sign = -1.;
      second_C_angle = path_RSR(2);
      second_C_sign = -1.;
      S_dist = path_RSR(1)*rho;
    } else {
      first_C_angle = path_RSL(0);
      first_C_sign = -1.;
      second_C_angle = path_RSL(2);
      second_C_sign = 1.;
      S_dist = path_RSL(1)*rho;
    }

    if (first_C_angle >= M_PI) { // Tested
      if (verbose) {
        std::cout << "Elongating first C segment of O1 path" << std::endl;
      }
      // O1. Elongate first C segment
      RowMatrixXd turns = RowMatrixXd::Zero(6, 2);

      // Straight (this is a choice I made, we could also turn for a bit)
      turns(0, 0) = 0.;
      turns(0, 1) = (s - l_m)/2;

      // Turn
      turns(1, 0) = first_C_sign;
      turns(1, 1) = M_PI*rho;

      // Straight
      turns(2, 0) = 0.;
      turns(2, 1) = (s - l_m)/2;

      // Turn
      turns(3, 0) = first_C_sign;
      turns(3, 1) = first_C_angle*rho - rho*M_PI;

      // Straight (from Dubins path)
      turns(4, 0) = 0.;
      turns(4, 1) = S_dist;

      // Turn (from Dubins path)
      turns(5, 0) = second_C_sign;
      turns(5, 1) = second_C_angle*rho;

      return turns;
    }

    if (second_C_angle >= M_PI) { // Tested
      if (verbose) {
        std::cout << "Elongating second C segment of O2 path" << std::endl;
      }
      // O2. Elongate second C segment
      RowMatrixXd turns = RowMatrixXd::Zero(6, 2);

      // Turn (from Dubins path)
      turns(0, 0) = first_C_sign;
      turns(0, 1) = first_C_angle*rho;

      // Straight (from Dubins path)
      turns(1, 0) = 0.;
      turns(1, 1) = S_dist;

      // Straight (this is a choice I made, we could also turn for a bit)
      turns(2, 0) = 0.;
      turns(2, 1) = (s - l_m)/2;

      // Turn
      turns(3, 0) = second_C_sign;
      turns(3, 1) = M_PI*rho;

      // Straight
      turns(4, 0) = 0.;
      turns(4, 1) = (s - l_m)/2;

      // Left turn
      turns(5, 0) = second_C_sign;
      turns(5, 1) = second_C_angle*rho - rho*M_PI;

      return turns;
    }

    if (S_dist >= 4*rho) { // Tested both cases
      if (verbose) {
        std::cout << "Elongating S segment of O3 path" << std::endl;
      }
      // O3. Elongate S segment

      // Elongate by turning left immediately after the first C segment (choice I made, we could also go straight for a bit, and/or elongate via right turn). First, check if we need an LRL or LSRSL
      if (s - l_m < 2*M_PI*rho - 4*rho) { // Tested
        if (verbose) {
          std::cout << "Short case" << std::endl;
        }
        // Elongate via LRL
        RowMatrixXd turns = RowMatrixXd::Zero(6, 2);

        // Turn (from Dubins path)
        turns(0, 0) = first_C_sign;
        turns(0, 1) = first_C_angle*rho;

        // LRL

        // All four of these arcs (left arc, right arc split in half, and left arc) will be the same arclength. Need the arclength to sum to s - l_m + 4*rho*sin(theta). So 4*rho*theta = s - l_m + 4*rho*sin(theta)
        double theta = root_find(rho, s, l_m);

        turns(1, 0) = 1.;
        turns(1, 1) = rho*theta;

        turns(2, 0) = -1.;
        turns(2, 1) = 2*rho*theta;

        turns(3, 0) = 1.;
        turns(3, 1) = rho*theta;

        // Straight (from Dubins path, but shorter)
        turns(4, 0) = 0.;
        turns(4, 1) = S_dist - 4*rho*sin(theta);

        // Turn (from Dubins path)
        turns(5, 0) = second_C_sign;
        turns(5, 1) = second_C_angle*rho;

        return turns;
      } else { // Tested
        if (verbose) {
          std::cout << "Long case" << std::endl;
        }
        // Elongate via LSRSL. In doing so, we subtract 4rho from the Dubins S segment. The added C segments travel 2*pi*rho. The added S segments must therefore travel a total distance of s - l_m + 4rho - 2*pi*rho
        RowMatrixXd turns = RowMatrixXd::Zero(8, 2);

        // Turn (from Dubins path)
        turns(0, 0) = first_C_sign;
        turns(0, 1) = first_C_angle*rho;

        // Left turn
        turns(1, 0) = 1.;
        turns(1, 1) = M_PI/2*rho;

        // Straight
        turns(2, 0) = 0.;
        turns(2, 1) = (s - l_m + 4*rho - 2*M_PI*rho)/2*rho;

        // Right turn
        turns(3, 0) = -1.;
        turns(3, 1) = M_PI*rho;

        // Straight
        turns(4, 0) = 0.;
        turns(4, 1) = (s - l_m + 4*rho - 2*M_PI*rho)/2*rho;

        // Left turn
        turns(5, 0) = 1.;
        turns(5, 1) = M_PI/2*rho;

        // Straight (from Dubins path, but shorter)
        turns(6, 0) = 0.;
        turns(6, 1) = S_dist - 4*rho;

        // Turn (from Dubins path)
        turns(7, 0) = second_C_sign;
        turns(7, 1) = second_C_angle*rho;

        return turns;
      }
    }

    double ldist = (center3_L - center1_L).norm();
    double rdist = (center3_R - center1_R).norm();
    // O4 and O5
    if (rdist >= 4*rho || ldist >= 4*rho) { // Tested both cases
      if (verbose) {
        if (rdist >= 4*rho) {
          std::cout << "Elongating O4 path" << std::endl;
        } else {
          std::cout << "Elongating O5 path" << std::endl;
        }
      }
      Vector2d midpoint;
      if (rdist >= 4*rho) {
        midpoint = (center1_R + center3_R)/2;
      } else {
        midpoint = (center1_L + center3_L)/2;
      }
      Vector2d direction = center3_R - center1_R;
      Vector2d perp(direction(1), -direction(0));

      // Pose at start of S segment
      double theta_S_start = theta_0 + first_C_sign*first_C_angle*rho;
      double c_S_start = cos(theta_S_start);
      double s_S_start = sin(theta_S_start);
      double x_S_start = x_0 + rho/first_C_sign*(-s_0 + s_S_start);
      double y_S_start = y_0 + rho/first_C_sign*(c_0 - c_S_start);

      // Pose at end of S segment
      double x_S_end = x_S_start + S_dist*c_S_start;
      double y_S_end = y_S_start + S_dist*s_S_start;

      Vector2d S_dir = Vector2d(x_S_end - x_S_start, y_S_end - y_S_start).normalized();
      std::cout << atan2(S_dir(1), S_dir(0)) << std::endl;
      Vector2d perp_to_S = Vector2d(S_dir(1), -S_dir(0)).normalized();

      // I'm assuming without proof for now that we have an intersection here
      Vector2d intersect_point = line_segments_intersect_point(midpoint(0), midpoint(1), midpoint(0) + perp(0)*10*rho, midpoint(1) + perp(1)*10*rho,
                                                               x_S_start, y_S_start, x_S_end, y_S_end);
      if (std::isinf(intersect_point(0))) {
        perp = -perp;
        intersect_point = line_segments_intersect_point(midpoint(0), midpoint(1), midpoint(0) + perp(0)*10*rho, midpoint(1) + perp(1)*10*rho,
                                                        x_S_start, y_S_start, x_S_end, y_S_end);
      }

      Vector2d disc_motion_dir = (midpoint - intersect_point).normalized();
      if (perp_to_S.dot(disc_motion_dir) < 0) {
        perp_to_S = -perp_to_S;
      }

      double s_angle = disc_motion_dir(0)*S_dir(1) - disc_motion_dir(1)*S_dir(0);;
      double additional_move_dist = rho/abs(s_angle);
  
      Vector2d disc_start_point = intersect_point + additional_move_dist*(intersect_point - midpoint).normalized();

      double disc_motion_dist_lb = 0.;
      double disc_motion_dist_ub = 1.;

      // Find upper bound
      int max_iter = 20;
      for (int i = 0; i < max_iter; ++i) {
        Vector2d disc_pos = disc_start_point + disc_motion_dir*disc_motion_dist_ub;
        Vector2d wp_pos = disc_pos + perp_to_S*rho;
        double l = path_length_with_waypoint(x_0, y_0, theta_0, wp_pos(0), wp_pos(1), theta_S_start, x_f, y_f, theta_f, rho);
        if (l >= s) {
          break;
        }
        disc_motion_dist_ub *= 2;
      }

      RowMatrixXd turns;
      for (int i = 0; i < max_iter; ++i) {
        double disc_motion_dist = (disc_motion_dist_ub + disc_motion_dist_lb)/2;
        Vector2d disc_pos = disc_start_point + disc_motion_dir*disc_motion_dist;
        Vector2d wp_pos = disc_pos + perp_to_S*rho;
        double l = path_length_with_waypoint(x_0, y_0, theta_0, wp_pos(0), wp_pos(1), theta_S_start, x_f, y_f, theta_f, rho);
        if (l >= s) {
          disc_motion_dist_ub = disc_motion_dist;
        } else {
          disc_motion_dist_lb = disc_motion_dist;
        }
      }

      Vector2d disc_pos = disc_start_point + disc_motion_dir*disc_motion_dist_lb;
      Vector2d wp_pos = disc_pos + perp_to_S*rho;
      return path_with_waypoint(x_0, y_0, theta_0, wp_pos(0), wp_pos(1), theta_S_start, x_f, y_f, theta_f, rho);
    }
  }

  double l1 = std::max(l_LRL_s, l_RLR_s);
  double l2 = l_m + 2*M_PI*rho;
  if (l_LRL_l < l2) {
    l2 = l_LRL_l;
  }
  if (l_RLR_l < l2) {
    l2 = l_RLR_l;
  }
  if (l_RSR != l_m && l_RSR < l2) {
    l2 = l_RSR;
  }
  if (l_RSL != l_m && l_RSL < l2) {
    l2 = l_RSL;
  }
  if (l_LSR != l_m && l_LSR < l2) {
    l2 = l_LSR;
  }
  if (l_LSL != l_m && l_LSL < l2) {
    l2 = l_LSL;
  }

  if (!(s <= l1 || s >= l2)) {
    return std::numeric_limits<double>::infinity()*RowMatrixXd::Ones(1, 2);
  }

  if (s <= l1) { // Tested
    if (verbose) {
      std::cout << "Elongating nabla O path, with s <= l1" << std::endl;
    }
    // Need to move a disc such that the path around the disc approaches the short LRL or RLR path
    double first_C_sign = 0.;
    double first_C_dist_low = 0.;
    double first_C_dist_high = 0.;
    if (s <= l_RLR_s) {
      first_C_sign = -1.;
      if (l_m == l_RSR) {
        first_C_dist_low = path_RSR(0);
      } else if (l_m == l_RSL) {
        first_C_dist_low = path_RSL(0);
      }

      VectorXd path_RLR = ccc_inverse(x_0, y_0, theta_0, x_f, y_f, theta_f, rho, -1., rho, 1., rho, -1.);
      first_C_dist_high = path_RLR(0);
    } else { // s <= l_LRL_s
      first_C_sign = 1.;
      if (l_m == l_LSL) {
        first_C_dist_low = path_LSL(0);
      } else if (l_m == l_LSR) {
        first_C_dist_low = path_LSR(0);
      }

      VectorXd path_LRL = ccc_inverse(x_0, y_0, theta_0, x_f, y_f, theta_f, rho, 1., rho, -1., rho, 1.);
      first_C_dist_high = path_LRL(0);
    }

    int max_iter = 10;
    double q_f[3] = {x_f, y_f, theta_f};
    DubinsPath path;
    for (int i = 0; i < max_iter; ++i) {
      double first_C_dist_mid = (first_C_dist_low + first_C_dist_high)/2;
      // Compute position and heading after traveling this first C segment
      double theta_mid = theta_0 + first_C_sign*first_C_dist_mid/rho;
      double c_mid = cos(theta_mid);
      double s_mid = sin(theta_mid);
      double x_mid = x_0 + rho/first_C_sign*(-s_0 + s_mid);
      double y_mid = y_0 + rho/first_C_sign*(c_0 - c_mid);

      // Compute Dubins path length
      double q_mid[3] = {x_mid, y_mid, theta_mid};
      dubins_shortest_path(&path, q_mid, q_f, rho);
      double l = dubins_path_length(&path) + first_C_dist_mid;
      if (l > s) {
        first_C_dist_high = first_C_dist_mid;
      } else {
        first_C_dist_low = first_C_dist_mid;
      }
    }

    RowMatrixXd turns = RowMatrixXd::Zero(4, 2);
    turns(0, 0) = first_C_sign;
    turns(0, 1) = first_C_dist_low;

    double theta_low = theta_0 + first_C_sign*first_C_dist_low/rho;
    double c_low = cos(theta_low);
    double s_low = sin(theta_low);
    double x_low = x_0 + rho/first_C_sign*(-s_0 + s_low);
    double y_low = y_0 + rho/first_C_sign*(c_0 - c_low);

    turns.bottomRows(3) = turns_for_dubins_path(x_low, y_low, theta_low, x_f, y_f, theta_f, rho);
    return turns;
  }

  if (l2 == l_m + 2*M_PI*rho) { // Tested
  // if (s >= l_m + 2*M_PI*rho) { // Uncomment this line and comment the above to force elongation via this method (for testing)
    if (verbose) {
      std::cout << "Elongating nabla O path, with l2 == l_m + 2*M_PI*rho" << std::endl;
    }
    // Add a loop to beginning of the path, then elongate via parallel tangents

    // Dubins path has to be CSC 
    double first_C_sign = 0.;
    double first_C_angle = 0.;
    double second_C_sign = 0.;
    double second_C_angle = 0.;
    double S_dist = 0.;
    if (l_m == l_LSL) {
      first_C_angle = path_LSL(0);
      first_C_sign = 1.;
      second_C_angle = path_LSL(2);
      second_C_sign = 1.;
      S_dist = path_LSL(1)*rho;
    } else if (l_m == l_LSR) {
      first_C_angle = path_LSR(0);
      first_C_sign = 1.;
      second_C_angle = path_LSR(2);
      second_C_sign = -1.;
      S_dist = path_LSR(1)*rho;
    } else if (l_m == l_RSR) {
      first_C_angle = path_RSR(0);
      first_C_sign = -1.;
      second_C_angle = path_RSR(2);
      second_C_sign = -1.;
      S_dist = path_RSR(1)*rho;
    } else {
      first_C_angle = path_RSL(0);
      first_C_sign = -1.;
      second_C_angle = path_RSL(2);
      second_C_sign = 1.;
      S_dist = path_RSL(1)*rho;
    }

    RowMatrixXd turns = RowMatrixXd::Zero(7, 2);

    // Straight (this is a choice I made, we could also turn for a bit)
    turns(0, 0) = 0.;
    turns(0, 1) = (s - 2*M_PI*rho - l_m)/2;

    turns(1, 0) = 1.;
    turns(1, 1) = M_PI*rho;

    // Straight
    turns(2, 0) = 0.;
    turns(2, 1) = (s - 2*M_PI*rho - l_m)/2;

    // Turn
    turns(3, 0) = 1.;
    turns(3, 1) = M_PI*rho;

    // Turn (from Dubins path)
    turns(4, 0) = first_C_sign;
    turns(4, 1) = first_C_angle*rho;

    // Straight (from Dubins path)
    turns(5, 0) = 0.;
    turns(5, 1) = S_dist;

    // Turn (from Dubins path)
    turns(6, 0) = second_C_sign;
    turns(6, 1) = second_C_angle*rho;

    return turns;
  }

  if (l2 == l_LRL_l) { // Tested
    if (verbose) {
      std::cout << "Elongating nabla O path, with l2 == long LRL path length" << std::endl;
    }

    RowMatrixXd turns = RowMatrixXd::Zero(6, 2);

    // Left turn
    turns(0, 0) = 1.;
    turns(0, 1) = LRL_l_dist1*rho;

    // Straight (this is a choice I made, we could also turn right for a bit)
    turns(1, 0) = 0.;
    turns(1, 1) = (s - l2)/2;

    // Right turn
    turns(2, 0) = -1.;
    turns(2, 1) = rho*M_PI;

    // Straight
    turns(3, 0) = 0.;
    turns(3, 1) = (s - l2)/2;

    // Right turn
    turns(4, 0) = -1.;
    turns(4, 1) = LRL_l_dist2*rho - rho*M_PI; // Total distance for the original R segment minus distance we traveled in the first split-up R segment

    // Left turn
    turns(5, 0) = 1.;
    turns(5, 1) = LRL_l_dist3*rho;

    return turns;
  }

  if (l2 == l_RLR_l) { // Tested
    if (verbose) {
      std::cout << "Elongating nabla O path, with l2 == long RLR path length" << std::endl;
    }
    RowMatrixXd turns = RowMatrixXd::Zero(6, 2);

    // Right turn
    turns(0, 0) = -1.;
    turns(0, 1) = RLR_l_dist1*rho;

    // Straight (this is a choice I made, we could also turn left for a bit)
    turns(1, 0) = 0.;
    turns(1, 1) = (s - l2)/2;

    // Left turn
    turns(2, 0) = 1.;
    turns(2, 1) = rho*M_PI;

    // Straight
    turns(3, 0) = 0.;
    turns(3, 1) = (s - l2)/2;

    // Left turn
    turns(4, 0) = 1.;
    turns(4, 1) = RLR_l_dist2*rho - rho*M_PI; // Total distance for the original R segment minus distance we traveled in the first split-up R segment

    // Right turn
    turns(5, 0) = -1.;
    turns(5, 1) = RLR_l_dist3*rho;

    return turns;
  }

  // Uncomment one of the below lines and comment the cases for 2pi-based elongation and CCC based elongation of nabla O paths to force this elongation method (for testing)
  // l2 = l_RSR; // To test second segment elongation
  // l2 = l_LSR; // To test first segment elongation
  if (l2 == l_RSR || l2 == l_RSL || l2 == l_LSR || l2 == l_LSL) { // Tested
    if (verbose) {
      std::cout << "Elongating nabla O path, with l2 == some CSC path length" << std::endl;
    }
    assert(l2 != l_m);

    double first_C_sign = 0.;
    double first_C_angle = 0.;
    double second_C_sign = 0.;
    double second_C_angle = 0.;
    double S_dist = 0.;

    if (l2 == l_RSR) {
      first_C_angle = path_RSR(0);
      first_C_sign = -1.;
      second_C_angle = path_RSR(2);
      second_C_sign = -1.;
      S_dist = path_RSR(1)*rho;
    } else if (l2 == l_RSL) {
      first_C_angle = path_RSL(0);
      first_C_sign = -1.;
      second_C_angle = path_RSL(2);
      second_C_sign = 1.;
      S_dist = path_RSL(1)*rho;
    } else if (l2 == l_LSR) {
      first_C_angle = path_LSR(0);
      first_C_sign = 1.;
      second_C_angle = path_LSR(2);
      second_C_sign = -1.;
      S_dist = path_LSR(1)*rho;
    } else if (l2 == l_LSL) {
      first_C_angle = path_LSL(0);
      first_C_sign = 1.;
      second_C_angle = path_LSL(2);
      second_C_sign = 1.;
      S_dist = path_LSL(1)*rho;
    }

    if (first_C_angle >= M_PI) {
      if (verbose) {
        std::cout << "Elongating first C segment" << std::endl;
      }
      // Elongate first C segment
      RowMatrixXd turns = RowMatrixXd::Zero(6, 2);

      // Straight (this is a choice I made, we could also turn for a bit)
      turns(0, 0) = 0.;
      turns(0, 1) = (s - l2)/2;

      // Turn
      turns(1, 0) = first_C_sign;
      turns(1, 1) = M_PI*rho;

      // Straight
      turns(2, 0) = 0.;
      turns(2, 1) = (s - l2)/2;

      // Turn
      turns(3, 0) = first_C_sign;
      turns(3, 1) = first_C_angle*rho - rho*M_PI;

      // Straight (from Dubins path)
      turns(4, 0) = 0.;
      turns(4, 1) = S_dist;

      // Turn (from Dubins path)
      turns(5, 0) = second_C_sign;
      turns(5, 1) = second_C_angle*rho;

      return turns;
    } else if (second_C_angle >= M_PI) {
      if (verbose) {
        std::cout << "Elongating second C segment" << std::endl;
      }
      // Elongate second C segment
      RowMatrixXd turns = RowMatrixXd::Zero(6, 2);

      // Turn (from Dubins path)
      turns(0, 0) = first_C_sign;
      turns(0, 1) = first_C_angle*rho;

      // Straight (from Dubins path)
      turns(1, 0) = 0.;
      turns(1, 1) = S_dist;

      // Straight (this is a choice I made, we could also turn for a bit)
      turns(2, 0) = 0.;
      turns(2, 1) = (s - l2)/2;

      // Turn
      turns(3, 0) = second_C_sign;
      turns(3, 1) = M_PI*rho;

      // Straight
      turns(4, 0) = 0.;
      turns(4, 1) = (s - l2)/2;

      // Left turn
      turns(5, 0) = second_C_sign;
      turns(5, 1) = second_C_angle*rho - rho*M_PI;

      return turns;
    } else {
      throw std::runtime_error("One of the C segments in a CSC path should have parallel tangents");
    }
  }

  throw std::runtime_error("Should not reach the end of elongated_dubins_path function");
  return std::numeric_limits<double>::infinity()*RowMatrixXd::Ones(1, 2); // We shouldn't ever get here
}
