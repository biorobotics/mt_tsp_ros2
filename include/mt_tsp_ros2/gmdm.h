#pragma once
#include <Eigen/Dense>
#include <random>
#include <cmath>
#include <iostream>
#include <product.hpp>

using namespace Eigen;
typedef Ref<Matrix<double, Dynamic, 1>> VectorXdRef;

double mod(double a, double m) {
  return a - m*std::floor(a/m);
}

// https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

double theta_ij(double theta_i, double theta_j, double w_i) {
  return mod(theta_i - theta_j, 2*M_PI*sgn<double>(w_i));
}

VectorXd csc_inverse(double x_0, double y_0, double theta_0, double x_f, double y_f, double theta_f, double v_1, double w_1, double v_2, double v_3, double w_3) {
  // Turning radius per segment
  double r_1 = v_1/w_1;
  double r_3 = v_3/w_3;

  double a = x_f - x_0 + r_1*sin(theta_0) - r_3*sin(theta_f);
  double b = y_f - y_0 - r_1*cos(theta_0) + r_3*cos(theta_f);

  double r_31 = r_3 - r_1;

  // Check reachable
  double a2 = a*a;
  double b2 = b*b;
  if (a2 + b2 < r_31*r_31) {
    return std::numeric_limits<double>::infinity()*VectorXd::Ones(1);
  }

  double theta_1 = asin(-r_31/sqrt(a2 + b2)) - atan2(-b, a);
  double theta_10 = theta_ij(theta_1, theta_0, w_1);
  double tau_1 = theta_10/w_1;

  // delta_i is the distance traveled in segment i
  double delta_2 = sqrt(a2 + b2 - r_31*r_31);
  double tau_2 = delta_2/v_2;

  double theta_31 = theta_ij(theta_f, theta_1, w_3);
  double tau_3 = theta_31/w_3;

  double dist = v_1*tau_1 + delta_2 + v_3*tau_3;

  assert(dist >= 0);
  VectorXd ret(4);
  ret(0) = tau_1;
  ret(1) = tau_2;
  ret(2) = tau_3;
  ret(3) = dist;
  return ret;
}

VectorXd ccc_inverse(double x_0, double y_0, double theta_0, double x_f, double y_f, double theta_f, double v_1, double w_1, double v_2, double w_2, double v_3, double w_3) {
  // Turning radius per segment
  double r_1 = v_1/w_1;
  double r_2 = v_2/w_2;
  double r_3 = v_3/w_3;
  double r_12 = r_1 - r_2;
  double r_23 = r_2 - r_3;
  double a = x_f - x_0 + r_1*sin(theta_0) - r_3*sin(theta_f);
  double b = y_f - y_0 - r_1*cos(theta_0) + r_3*cos(theta_f);

  // Check reachable
  double r_13 = r_1 - r_3;
  double a2 = a*a;
  double b2 = b*b;
  if (r_13*r_13 > a2 + b2 || a2 + b2 > (r_12 - r_23)*(r_12 - r_23)) {
    return std::numeric_limits<double>::infinity()*VectorXd::Ones(1);
  }

  double theta_1 = M_PI - asin((a2 + b2 + r_12*r_12 - r_23*r_23)/(2*r_12*sqrt(a2 + b2))) - atan2(-b, a);
  double theta_2 = M_PI - asin((a2 + b2 + r_23*r_23 - r_12*r_12)/(2*r_23*sqrt(a2 + b2))) - atan2(-b, a);

  double theta_10 = theta_ij(theta_1, theta_0, w_1);
  double theta_21 = theta_ij(theta_2, theta_1, w_2);
  double theta_32 = theta_ij(theta_f, theta_2, w_3);

  double tau_1 = theta_10/w_1;
  double tau_2 = theta_21/w_2;
  double tau_3 = theta_32/w_3;

  double dist = v_1*tau_1 + v_2*tau_2 + v_3*tau_3;
  assert(dist >= 0);
  VectorXd ret(4);
  ret(0) = tau_1;
  ret(1) = tau_2;
  ret(2) = tau_3;
  ret(3) = dist;
  return ret;
}

VectorXd csc_inverse_free_v2(double x_0, double y_0, double theta_0, double x_f, double y_f, double theta_f, double v_1, double w_1, double v_3, double w_3, double T, double v_min, double v_max) {
  // Turning radius per segment
  double r_1 = v_1/w_1;
  double r_3 = v_3/w_3;

  double a = x_f - x_0 + r_1*sin(theta_0) - r_3*sin(theta_f);
  double b = y_f - y_0 - r_1*cos(theta_0) + r_3*cos(theta_f);

  double r_31 = r_3 - r_1;

  // Check reachable
  double a2 = a*a;
  double b2 = b*b;
  if (a2 + b2 < r_31*r_31) {
    return std::numeric_limits<double>::infinity()*VectorXd::Ones(1);
  }

  double theta_1 = asin(-r_31/sqrt(a2 + b2)) - atan2(-b, a);
  double theta_10 = theta_ij(theta_1, theta_0, w_1);
  double tau_1 = theta_10/w_1;

  // delta_i is the distance traveled in segment i
  double delta_2 = sqrt(a2 + b2 - r_31*r_31);

  double theta_31 = theta_ij(theta_f, theta_1, w_3);
  double tau_3 = theta_31/w_3;

  double v_2;
  double tau_2;
  if (std::isinf(T)) {
    v_2 = v_max;
    tau_2 = delta_2/v_2;
  } else {
    tau_2 = T - tau_1 - tau_3;
    if (v_min*tau_2 <= delta_2 && delta_2 <= v_max*tau_2) {
      assert(tau_2 >= 0);
      v_2 = delta_2/tau_2;
    } else {
      return std::numeric_limits<double>::infinity()*VectorXd::Ones(1);
    }
  }

  double dist = v_1*tau_1 + delta_2 + v_3*tau_3;

  assert(dist >= 0);
  VectorXd ret(5);
  ret(0) = tau_1;
  ret(1) = tau_2;
  ret(2) = tau_3;
  ret(3) = v_2;
  ret(4) = dist;
  return ret;
}

VectorXd ccc_inverse_free_v2(double x_0, double y_0, double theta_0, double x_f, double y_f, double theta_f, double v_1, double w_1, double r_2, double sign_w_2, double v_3, double w_3, double T, double v_min, double v_max, double w_max) {
  // Turning radius per segment
  double r_1 = v_1/w_1;
  double r_3 = v_3/w_3;
  double r_12 = r_1 - r_2;
  double r_23 = r_2 - r_3;
  double a = x_f - x_0 + r_1*sin(theta_0) - r_3*sin(theta_f);
  double b = y_f - y_0 - r_1*cos(theta_0) + r_3*cos(theta_f);

  // Check reachable
  double r_13 = r_1 - r_3;
  double a2 = a*a;
  double b2 = b*b;
  if (r_13*r_13 > a2 + b2 || a2 + b2 > (r_12 - r_23)*(r_12 - r_23)) {
    return std::numeric_limits<double>::infinity()*VectorXd::Ones(1);
  }

  double theta_1 = M_PI - asin((a2 + b2 + r_12*r_12 - r_23*r_23)/(2*r_12*sqrt(a2 + b2))) - atan2(-b, a);
  double theta_2 = M_PI - asin((a2 + b2 + r_23*r_23 - r_12*r_12)/(2*r_23*sqrt(a2 + b2))) - atan2(-b, a);

  double theta_10 = theta_ij(theta_1, theta_0, w_1);
  double theta_21 = theta_ij(theta_2, theta_1, sign_w_2);
  double theta_32 = theta_ij(theta_f, theta_2, w_3);

  double tau_1 = theta_10/w_1;
  double tau_3 = theta_32/w_3;

  double w_2;
  double tau_2;
  if (std::isinf(T)) {
    w_2 = sign_w_2*w_max;
    tau_2 = theta_21/w_2;
  } else {
    tau_2 = T - tau_1 - tau_3;
    if (abs(theta_21) <= w_max*tau_2) {
      w_2 = theta_21/tau_2;
    } else {
      return std::numeric_limits<double>::infinity()*VectorXd::Ones(1);
    }
  }

  double v_2 = r_2*w_2;
  if (v_min > v_2 || v_2 > v_max) {
    return std::numeric_limits<double>::infinity()*VectorXd::Ones(1);
  }

  double dist = v_1*tau_1 + v_2*tau_2 + v_3*tau_3;
  assert(dist >= 0);
  VectorXd ret(6);
  ret(0) = tau_1;
  ret(1) = tau_2;
  ret(2) = tau_3;
  ret(3) = v_2;
  ret(4) = w_2;
  ret(5) = dist;
  return ret;
}

double random_gmdm_planner(VectorXdRef best_trj, double t_0, double x_0, double y_0, double theta_0, 
                           double t_f, double x_f, double y_f, double theta_f, 
                           double w_max, double v_min, double v_max,
                           int num_sample, int random_seed) {
  assert(best_trj.size() == 9);
  double T = t_f - t_0;

  double min_dist = std::numeric_limits<double>::infinity();

  std::mt19937 rng;
  rng.seed(random_seed);

  std::uniform_int_distribution<> binary_dist(0, 1);
  std::uniform_real_distribution<> v_dist(v_min, v_max);
  std::uniform_real_distribution<> w_dist(-w_max, w_max);
  std::uniform_real_distribution<> r_2_dist(v_min/w_max, v_max/w_max);

  for (int sample_idx = 0; sample_idx < num_sample; ++sample_idx) {
    bool do_CSC = binary_dist(rng);
    if (do_CSC) {
      double v_1 = v_dist(rng);
      // double w_1 = w_dist(rng);
      double w_1 = (binary_dist(rng) ? 1. : -1.)*w_max;
      double v_3 = v_dist(rng);
      // double w_3 = w_dist(rng);
      double w_3 = (binary_dist(rng) ? 1. : -1.)*w_max;
      VectorXd res = csc_inverse_free_v2(x_0, y_0, theta_0, x_f, y_f, theta_f, v_1, w_1, v_3, w_3, T, v_min, v_max);
      if (std::isinf(res(0))) {
        continue;
      }
      assert(res.size() == 5);
      double tau_1 = res(0);
      double tau_2 = res(1);
      double tau_3 = res(2);
      double v_2 = res(3);
      double dist = res(4);
      if (dist < min_dist) {
        min_dist = dist;
        best_trj(0) = v_1;
        best_trj(1) = w_1;
        best_trj(2) = v_2;
        best_trj(3) = 0.;
        best_trj(4) = v_3;
        best_trj(5) = w_3;
        best_trj(6) = tau_1;
        best_trj(7) = tau_2;
        best_trj(8) = tau_3;
      }
    } else {
      double v_1 = v_dist(rng);
      // double w_1 = w_dist(rng);
      double w_1 = (binary_dist(rng) ? 1. : -1.)*w_max;
      double r_2 = r_2_dist(rng);
      double sign_w_2 = binary_dist(rng) ? 1. : -1.;
      double v_3 = v_dist(rng);
      // double w_3 = w_dist(rng);
      double w_3 = (binary_dist(rng) ? 1. : -1.)*w_max;
      VectorXd res = ccc_inverse_free_v2(x_0, y_0, theta_0, x_f, y_f, theta_f, v_1, w_1, r_2, sign_w_2, v_3, w_3, T, v_min, v_max, w_max);
      if (std::isinf(res(0))) {
        continue;
      }
      assert(res.size()  == 6);
      double tau_1 = res(0);
      double tau_2 = res(1);
      double tau_3 = res(2);
      double v_2 = res(3);
      double w_2 = res(4);
      double dist = res(5);
      if (dist < min_dist) {
        min_dist = dist;
        best_trj(0) = v_1;
        best_trj(1) = w_1;
        best_trj(2) = v_2;
        best_trj(3) = w_2;
        best_trj(4) = v_3;
        best_trj(5) = w_3;
        best_trj(6) = tau_1;
        best_trj(7) = tau_2;
        best_trj(8) = tau_3;
      }
    }
  }

  return min_dist;
}

double gmdm_planner_primitives(VectorXdRef best_trj, double t_0, double x_0, double y_0, double theta_0, 
                               double t_f, double x_f, double y_f, double theta_f, 
                               double w_max, double v_min, double v_max,
                               int num_sample) {
  assert(best_trj.size() == 9);
  double T = t_f - t_0;

  double min_dist = std::numeric_limits<double>::infinity();

  double w_1_vals[2] = {-w_max, w_max};
  double w_2_signs[3] = {-1., 0., 1.};
  double w_3_vals[2] = {-w_max, w_max};

  double v_1_vals[2] = {v_min, v_max};
  double r_2_vals[2] = {v_min/w_max, v_max/w_max};
  double v_3_vals[2] = {v_min, v_max};

  // CSC
  for (auto&& [v_1, w_1, v_3, w_3] : iter::product(v_1_vals, w_1_vals, v_3_vals, w_3_vals)) {
    VectorXd res = csc_inverse_free_v2(x_0, y_0, theta_0, x_f, y_f, theta_f, v_1, w_1, v_3, w_3, T, v_min, v_max);
    if (std::isinf(res(0))) {
      continue;
    }
    assert(res.size() == 5);
    double tau_1 = res(0);
    double tau_2 = res(1);
    double tau_3 = res(2);
    double v_2 = res(3);
    double dist = res(4);
    if (dist < min_dist) {
      min_dist = dist;
      best_trj(0) = v_1;
      best_trj(1) = w_1;
      best_trj(2) = v_2;
      best_trj(3) = 0.;
      best_trj(4) = v_3;
      best_trj(5) = w_3;
      best_trj(6) = tau_1;
      best_trj(7) = tau_2;
      best_trj(8) = tau_3;
    }
  }
  // CCC
  for (auto&& [v_1, w_1, r_2, sign_w_2, v_3, w_3] : iter::product(v_1_vals, w_1_vals, r_2_vals, w_2_signs, v_3_vals, w_3_vals)) {
    VectorXd res = ccc_inverse_free_v2(x_0, y_0, theta_0, x_f, y_f, theta_f, v_1, w_1, r_2, sign_w_2, v_3, w_3, T, v_min, v_max, w_max);
    if (std::isinf(res(0))) {
      continue;
    }
    assert(res.size()  == 6);
    double tau_1 = res(0);
    double tau_2 = res(1);
    double tau_3 = res(2);
    double v_2 = res(3);
    double w_2 = res(4);
    double dist = res(5);
    if (dist < min_dist) {
      min_dist = dist;
      best_trj(0) = v_1;
      best_trj(1) = w_1;
      best_trj(2) = v_2;
      best_trj(3) = w_2;
      best_trj(4) = v_3;
      best_trj(5) = w_3;
      best_trj(6) = tau_1;
      best_trj(7) = tau_2;
      best_trj(8) = tau_3;
    }
  }

  return min_dist;
}
