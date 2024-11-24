#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include "gmdm.h"

using namespace Eigen;
typedef Ref<Matrix<double, Dynamic, 1>> VectorXdRef;

double arclength(const Vector2d &v1, const Vector2d &v2, bool left, double rho) {
  double theta = atan2(v2(1), v2(0)) - atan2(v1(1), v1(0));
  if (theta < 0 && left) {
    theta += 2*M_PI;
  } else if (theta > 0 && !left) {
    theta -= 2*M_PI;
  }
  return abs(theta*rho);
}

bool check_elongation_possible(double x_0, double y_0, double theta_0, double x_f, double y_f, double theta_f, double s, double rho) {
  // LRL
  double LRL_dist_A = std::numeric_limits<double>::infinity();
  double LRL_dist_B = std::numeric_limits<double>::infinity();

  double c_0 = cos(theta_0);
  double s_0 = sin(theta_0);

  Vector2d dir_0(c_0, s_0);
  Vector2d perp_0(-s_0, c_0);
  Vector2d center1 = Vector2d(x_0, y_0) + perp_0*rho;

  double c_f = cos(theta_f);
  double s_f = sin(theta_f);
  Vector2d dir_f(c_f, s_f);
  Vector2d perp_f(-s_f, c_f);
  Vector2d center3 = Vector2d(x_f, y_f) + perp_f*rho;

  Vector2d V = center3 - center1;
  double D = V.norm();
  if (D <= 4*rho) {
    double gamma = atan2(V(1), V(0));
    double theta = acos(D/(4*rho));

    double theta_A = gamma + theta;
    double c_A = cos(theta_A);
    double s_A = sin(theta_A);
    Vector2d vec_A = Vector2d(c_A, s_A);
    Vector2d center2 = center1 + vec_A*rho*2;
    LRL_dist_A = arclength(-perp_0, vec_A, true, rho);
    LRL_dist_A += arclength(center1 - center2, center3 - center2, false, rho);
    LRL_dist_A += arclength(center2 - center3, -perp_f, true, rho);

    double theta_B = gamma - theta;
    double c_B = cos(theta_B);
    double s_B = sin(theta_B);
    Vector2d vec_B = Vector2d(c_B, s_B);
    center2 = center1 + vec_B*rho*2;
    LRL_dist_B = arclength(-perp_0, vec_B, true, rho);
    LRL_dist_B += arclength(center1 - center2, center3 - center2, false, rho);
    LRL_dist_B += arclength(center2 - center3, -perp_f, true, rho);
  }

  // RLR
  double RLR_dist_A = std::numeric_limits<double>::infinity();
  double RLR_dist_B = std::numeric_limits<double>::infinity();

  perp_0 = Vector2d(s_0, -c_0);
  center1 = Vector2d(x_0, y_0) + perp_0*rho;

  perp_f = Vector2d(s_f, -c_f);
  center3 = Vector2d(x_f, y_f) + perp_f*rho;

  V = center3 - center1;
  D = V.norm();
  if (D <= 4*rho) {
    double gamma = atan2(V(1), V(0));
    double theta = acos(D/(4*rho));

    double theta_A = gamma + theta;
    double c_A = cos(theta_A);
    double s_A = sin(theta_A);
    Vector2d vec_A = Vector2d(c_A, s_A);
    Vector2d center2 = center1 + vec_A*rho*2;
    RLR_dist_A = arclength(-perp_0, vec_A, false, rho);
    RLR_dist_A += arclength(center1 - center2, center3 - center2, true, rho);
    RLR_dist_A += arclength(center2 - center3, -perp_f, false, rho);

    double theta_B = gamma - theta;
    double c_B = cos(theta_B);
    double s_B = sin(theta_B);
    Vector2d vec_B = Vector2d(c_B, s_B);
    center2 = center1 + Vector2d(c_B, s_B)*rho*2;
    RLR_dist_B = arclength(-perp_0, vec_B, false, rho);
    RLR_dist_B += arclength(center1 - center2, center3 - center2, true, rho);
    RLR_dist_B += arclength(center2 - center3, -perp_f, false, rho);
  }

  double l_LRL_s = std::min(LRL_dist_A, LRL_dist_B);
  double l_RLR_s = std::min(RLR_dist_A, RLR_dist_B);

  double l_LRL_l = std::max(LRL_dist_A, LRL_dist_B);
  double l_RLR_l = std::max(RLR_dist_A, RLR_dist_B);

  double l_m = std::min(l_LRL_s, l_LRL_l);

  // LSL
  double l_LSL = std::numeric_limits<double>::infinity();
  VectorXd ret = csc_inverse(x_0, y_0, theta_0, x_f, y_f, theta_f, rho, 1., rho, rho, 1.);
  if (!std::isinf(ret(0))) {
    l_m = std::min(l_m, ret(3));
    l_LSL = ret(3);
  }

  // LSR
  double l_LSR = std::numeric_limits<double>::infinity();
  ret = csc_inverse(x_0, y_0, theta_0, x_f, y_f, theta_f, rho, 1., rho, rho, -1.);
  if (!std::isinf(ret(0))) {
    l_m = std::min(l_m, ret(3));
    l_LSR = ret(3);
  }

  // RSR
  double l_RSR = std::numeric_limits<double>::infinity();
  ret = csc_inverse(x_0, y_0, theta_0, x_f, y_f, theta_f, rho, -1., rho, rho, -1.);
  if (!std::isinf(ret(0))) {
    l_m = std::min(l_m, ret(3));
    l_RSR = ret(3);
  }

  // RSL
  double l_RSL = std::numeric_limits<double>::infinity();
  ret = csc_inverse(x_0, y_0, theta_0, x_f, y_f, theta_f, rho, -1., rho, rho, 1.);
  if (!std::isinf(ret(0))) {
    l_m = std::min(l_m, ret(3));
    l_RSL = ret(3);
  }

  double l1 = std::max(l_LRL_s, l_RLR_s);
  double l2 = l_m + 2*M_PI;
  if (l_LRL_l < l2) {
    l2 = l_LRL_l;
  }
  if (l_RLR_l < l2) {
    l2 = l_RLR_l;
  }
  if (l_RSR < l2) {
    l2 = l_RSR;
  }
  if (l_RSL < l2) {
    l2 = l_RSL;
  }
  if (l_LSR < l2) {
    l2 = l_LSR;
  }
  if (l_LSL < l2) {
    l2 = l_LSL;
  }

  return s <= l1 || s >= l2;
}
