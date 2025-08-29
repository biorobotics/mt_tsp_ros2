#pragma once
#include "mt_tsp_ros2/elongate_dubins_path.h"

// Returns a sequence of (turn direction, dist) pairs.
// turns should have two or more rows.
// amount is how much to shorten the segment in the sliding direction.
// direction = true means slide forward, false means slide backward.
// Performs sliding on segment turns.rows() - 2 (zero-based indexing) if forward, and on segment 1 if backward
RowMatrixXd ic_sliding(RowMatrixXdRef_const turns, double x_0, double y_0, double theta_0, double x_f, double y_f, double theta_f, double rho, double amount, bool direction) {
  int slide_idx = direction ? turns.rows() - 2 : 1;
  if (turns(slide_idx, 0) == 0) {
    throw std::runtime_error("Cannot slide an S segment");
  }

  int turn_to_shorten = direction ? slide_idx + 1 : slide_idx - 1;
  // Compute configuration after/before shortened segment (after if sliding backward, before if sliding forward)
  double x_mid;
  double y_mid;
  double theta_mid;
  double dist = turns(turn_to_shorten, 1) - amount;
  if (turns(turn_to_shorten, 0) == 0) {
    // S segment
    if (direction) {
      // Sliding forward
      double ctheta_f = cos(theta_f);
      double stheta_f = sin(theta_f);
      x_mid = x_f - ctheta_f*dist;
      y_mid = y_f - stheta_f*dist;
      theta_mid = theta_f;
    } else {
      // Sliding backward
      double ctheta_0 = cos(theta_0);
      double stheta_0 = sin(theta_0);
      x_mid = x_0 + ctheta_0*dist;
      y_mid = y_0 + stheta_0*dist;
      theta_mid = theta_0;
    }
  } else {
    // C segment
    double C_sign = turns(turn_to_shorten, 0);
    if (direction) {
      // Sliding forward
      double c_f = cos(theta_f);
      double s_f = sin(theta_f);
      theta_mid = theta_f - C_sign*dist/rho;
      double c_mid = cos(theta_mid);
      double s_mid = sin(theta_mid);
      x_mid = x_f - rho/C_sign*(-s_mid + s_f);
      y_mid = y_f - rho/C_sign*(c_mid - c_f);
    } else {
      // Sliding backward
      double c_0 = cos(theta_0);
      double s_0 = sin(theta_0);
      theta_mid = theta_0 + C_sign*dist/rho;
      double c_mid = cos(theta_mid);
      double s_mid = sin(theta_mid);
      x_mid = x_0 + rho/C_sign*(-s_0 + s_mid);
      y_mid = y_0 + rho/C_sign*(c_0 - c_mid);
    }
  }

  // Compute Dubins path
  RowMatrixXd ret(4, 2);
  if (direction) {
    // Sliding forward
    ret(3, 0) = turns(2, 0);
    ret(3, 1) = dist;
    ret.topRows<3>() = turns_for_dubins_path(x_0, y_0, theta_0, x_mid, y_mid, theta_mid, rho);
  } else {
    // Sliding backward
    ret(0, 0) = turns(0, 0);
    ret(0, 1) = dist;
    ret.bottomRows<3>() = turns_for_dubins_path(x_mid, y_mid, theta_mid, x_f, y_f, theta_f, rho);
  }
  return ret;
}
