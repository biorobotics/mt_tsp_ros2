#pragma once
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "mt_tsp_ros2/time_constrained_dubins_planning/angle_mod.h"
#include "mt_tsp_ros2/elongate_dubins_path.h"
#include <iomanip>

using namespace Eigen;
typedef Matrix<double, Dynamic, Dynamic, RowMajor> RowMatrixXd;

const int max_bisection_iter = 100;

RowMatrixXd turns_for_CS_path(double x_0, double y_0, double theta_0, double x_f, double y_f, double rho, bool left_turn) {
  if (true) {
    double c_0 = cos(theta_0);
    double s_0 = sin(theta_0);

    Vector2d p_0(x_0, y_0);
    Vector2d dir_0(c_0, s_0);
    Vector2d perp_0(-s_0, c_0);
    Vector2d center;
    if (left_turn) {
      center = p_0 + rho*perp_0;
    } else {
      center = p_0 - rho*perp_0;
    }
    Vector2d P(x_f, y_f); 

    double dist_from_center = (P - center).norm();
    if (dist_from_center < rho) {
      return std::numeric_limits<double>::infinity()*RowMatrixXd::Ones(1, 2);
    }

    double dist_from_tangent_point = sqrt(dist_from_center*dist_from_center - rho*rho);

    double intersect_point_angle = acos(rho/dist_from_center);
    double Pangle = atan2(P(1) - center(1), P(0) - center(0));

    // There could be two tangent points. Here's the first candidate
    double total_angle = Pangle - intersect_point_angle;
    Vector2d center_to_tangent_point(rho*cos(total_angle), rho*sin(total_angle));
    Vector2d tangent_point = center + center_to_tangent_point;
    Vector2d tangent_point_to_P = P - tangent_point;
    Vector2d dir_at_tangent_point;
    if (left_turn) {
      dir_at_tangent_point(0) = -center_to_tangent_point(1);
      dir_at_tangent_point(1) = center_to_tangent_point(0);
    } else {
      dir_at_tangent_point(0) = center_to_tangent_point(1);
      dir_at_tangent_point(1) = -center_to_tangent_point(0);
    }

    if (tangent_point_to_P.dot(dir_at_tangent_point) < 0) {
      // Need to use the other candidate
      total_angle = Pangle + intersect_point_angle;
    }

    Vector2d rel_p_0 = p_0 - center;
    double rel_p_0_angle = atan2(rel_p_0(1), rel_p_0(0));

    double diff = angdiff(rel_p_0_angle, total_angle);
    if (diff < 0 && left_turn) {
      diff = diff + 2*M_PI;
    } else if (diff > 0 && !left_turn) {
      diff = diff - 2*M_PI;
    }

    if (!left_turn) {
      diff = -diff;
    }

    double C_dist = diff*rho;
    RowMatrixXd turns(2, 2);
    turns(0, 0) = left_turn ? rho : -rho;
    turns(0, 1) = C_dist;
    turns(1, 0) = 0;
    turns(1, 1) = dist_from_tangent_point;
    return turns;
  } else {
    // From the GDIP code
    double c_0 = cos(theta_0);
    double s_0 = sin(theta_0);

    Vector2d p_0(x_0, y_0);
    Vector2d dir_0(c_0, s_0);
    Vector2d perp_0(-s_0, c_0);

    Vector2d center;
    if (left_turn) {
      center = p_0 + rho*perp_0;
    } else {
      center = p_0 - rho*perp_0;
    }
    Vector2d P(x_f, y_f); 

    Vector2d center_to_P = P - center;
    double l = center_to_P.norm();
    double alpha = asin(rho/l);
    double center_to_P_angle = atan2(center_to_P(1), center_to_P(0));
    double tangent_direction = left_turn ? center_to_P_angle + alpha : center_to_P_angle - alpha;
    double diff = atan2(sin(tangent_direction), cos(tangent_direction)) - atan2(s_0, c_0);
    if (diff < 0 && left_turn) {
      diff += 2*M_PI;
    } else if (diff > 0 && !left_turn) {
      diff -= 2*M_PI;
    }

    RowMatrixXd turns(2, 2);
    turns(0, 0) = left_turn ? rho : -rho;
    turns(0, 1) = std::abs(diff)*rho;
    turns(1, 0) = 0;
    turns(1, 1) = l*cos(alpha);
    return turns;
  }
}

RowMatrixXd turns_for_LR_path(double x_0, double y_0, double theta_0, double x_f, double y_f, double rho1, double rho2, double xi, bool gt_xi) {
  RowMatrixXd turns(2, 2);

  double c_0 = cos(theta_0);
  double s_0 = sin(theta_0);

  Vector2d p_0(x_0, y_0);
  Vector2d dir_0(c_0, s_0);
  Vector2d perp_0(-s_0, c_0);

  double angle_p0 = atan2(-perp_0(1), -perp_0(0));

  Vector2d center_L = p_0 + rho1*perp_0;
  Vector2d P(x_f, y_f);

  double normL = (P - center_L).norm();

  double sum_rho = rho1 + rho2;
  double alpha = acos((normL*normL + sum_rho*sum_rho - rho2*rho2)/(2*normL*sum_rho));

  Vector2d P_wrt_center_L = P - center_L;

  double beta = atan2(P_wrt_center_L(1), P_wrt_center_L(0));

  double gamma = alpha + beta;

  // Intersection between first L and second R turning circles
  Vector2d intersect_dir(cos(gamma), sin(gamma));

  double angle_intersect_point_wrt_center_L = atan2(intersect_dir(1), intersect_dir(0));

  // Since atan2 returns values in [-pi, pi], diff1 is in range [-2pi, 2pi]
  double diff1 = angle_intersect_point_wrt_center_L - angle_p0;
  if (diff1 < 0) {
    diff1 += 2*M_PI; // Since we are turning left
  }

  if ((gt_xi && diff1 <= xi) || (!gt_xi && diff1 > xi)) {
    gamma = -alpha + beta;

    // Intersection between first R and second L turning circles
    intersect_dir = Vector2d(cos(gamma), sin(gamma));

    angle_intersect_point_wrt_center_L = atan2(intersect_dir(1), intersect_dir(0));

    // Since atan2 returns values in [-pi, pi], diff1 is in range [-2pi, 2pi]
    diff1 = angle_intersect_point_wrt_center_L - angle_p0;
    if (diff1 < 0) {
      diff1 += 2*M_PI; // Since we are turning left
    }

    if ((gt_xi && diff1 <= xi) || (!gt_xi && diff1 > xi)) {
      throw std::runtime_error("LR path computation failed during elongation");
    }
  }

  Vector2d center_R2 = center_L + sum_rho*intersect_dir;
  Vector2d intersect_point_wrt_center_R2 = -rho2*intersect_dir;
  Vector2d P_wrt_center_R2 = P - center_R2;
  double angleP = atan2(P_wrt_center_R2(1), P_wrt_center_R2(0));
  double angle_intersect_point_wrt_center_R2 = atan2(intersect_point_wrt_center_R2(1), intersect_point_wrt_center_R2(0));
  // Since atan2 returns values in [-pi, pi], diff2 is in range [-2pi, 2pi]
  double diff2 = angleP - angle_intersect_point_wrt_center_R2;
  if (diff2 > 0) { 
    diff2 -= 2*M_PI; // Since we are turning right
  }

  turns(0, 0) = rho1;
  turns(0, 1) = rho1*diff1;
  turns(1, 0) = -rho2;
  turns(1, 1) = -rho2*diff2;    
  return turns;
}

RowMatrixXd turns_for_one_sided_dubins_path(double x_0, double y_0, double theta_0, double x_f, double y_f, double rho) {
  RowMatrixXd turns(2, 2);

  double c_0 = cos(theta_0);
  double s_0 = sin(theta_0);

  Vector2d p_0(x_0, y_0);
  Vector2d dir_0(c_0, s_0);
  Vector2d perp_0(-s_0, c_0);
  Vector2d center_L = p_0 + rho*perp_0;
  Vector2d center_R = p_0 - rho*perp_0;
  Vector2d P(x_f, y_f);

  double normL = (P - center_L).norm();
  double normR = (P - center_R).norm();

  if (normL < rho) {
    // RL case
    double alpha = acos((normR*normR + 3*rho*rho)/(4*rho*normR));

    Vector2d P_wrt_center_R = P - center_R;

    double beta = atan2(P_wrt_center_R(1), P_wrt_center_R(0));

    double gamma = alpha + beta;

    // Intersection between first R and second L turning circles
    Vector2d intersect_point_wrt_center_R(rho*cos(gamma), rho*sin(gamma));
    Vector2d center_L2 = center_R + 2*intersect_point_wrt_center_R;
    Vector2d intersect_point_wrt_center_L2 = -intersect_point_wrt_center_R;
    Vector2d P_wrt_center_L2 = P - center_L2;
    double angleP = atan2(P_wrt_center_L2(1), P_wrt_center_L2(0));
    double angle_intersect_point = atan2(intersect_point_wrt_center_L2(1), intersect_point_wrt_center_L2(0));
    // Since atan2 returns values in [-pi, pi], diff2 is in range [-2pi, 2pi]
    double diff2 = angleP - angle_intersect_point;
    if (diff2 < 0) { 
      diff2 += 2*M_PI; // Since we are turning left
    }
    if (diff2 < M_PI) {
      gamma = -alpha + beta;

      // Intersection between first R and second L turning circles
      intersect_point_wrt_center_R = Vector2d(rho*cos(gamma), rho*sin(gamma));
      center_L2 = center_R + 2*intersect_point_wrt_center_R;
      intersect_point_wrt_center_L2 = -intersect_point_wrt_center_R;
      P_wrt_center_L2 = P - center_L2;
      angleP = atan2(P_wrt_center_L2(1), P_wrt_center_L2(0));
      angle_intersect_point = atan2(intersect_point_wrt_center_L2(1), intersect_point_wrt_center_L2(0));
      // Since atan2 returns values in [-pi, pi], diff2 is in range [-2pi, 2pi]
      diff2 = angleP - angle_intersect_point;

      if (diff2 < 0) {
        diff2 += 2*M_PI; // Since we are turning left
      }

      if (diff2 < M_PI) {
        throw std::runtime_error("Error in computing RL path");
      }
    }

    angle_intersect_point = atan2(intersect_point_wrt_center_R(1), intersect_point_wrt_center_R(0));
    double angle_p0 = atan2(perp_0(1), perp_0(0));

    // Since atan2 returns values in [-pi, pi], diff1 is in range [-2pi, 2pi]
    double diff1 = angle_intersect_point - angle_p0;
    if (diff1 > 0) {
      diff1 -= 2*M_PI; // Since we are turning right
    }

    turns(0, 0) = -rho;
    turns(0, 1) = -rho*diff1;
    turns(1, 0) = rho;
    turns(1, 1) = rho*diff2;
    return turns;
  } else if (normR < rho) {
    // LR case
    double alpha = acos((normL*normL + 3*rho*rho)/(4*rho*normL));

    Vector2d P_wrt_center_L = P - center_L;

    double beta = atan2(P_wrt_center_L(1), P_wrt_center_L(0));

    double gamma = alpha + beta;

    // Intersection between first L and second R turning circles
    Vector2d intersect_point_wrt_center_L(rho*cos(gamma), rho*sin(gamma));
    Vector2d center_R2 = center_L + 2*intersect_point_wrt_center_L;
    Vector2d intersect_point_wrt_center_R2 = -intersect_point_wrt_center_L;
    Vector2d P_wrt_center_R2 = P - center_R2;
    double angleP = atan2(P_wrt_center_R2(1), P_wrt_center_R2(0));
    double angle_intersect_point = atan2(intersect_point_wrt_center_R2(1), intersect_point_wrt_center_R2(0));
    // Since atan2 returns values in [-pi, pi], diff2 is in range [-2pi, 2pi]
    double diff2 = angleP - angle_intersect_point;
    if (diff2 > 0) { 
      diff2 -= 2*M_PI; // Since we are turning right
    }
    if (diff2 > -M_PI) {
      gamma = -alpha + beta;

      // Intersection between first R and second L turning circles
      intersect_point_wrt_center_L = Vector2d(rho*cos(gamma), rho*sin(gamma));
      center_R2 = center_L + 2*intersect_point_wrt_center_L;
      intersect_point_wrt_center_R2 = -intersect_point_wrt_center_L;
      P_wrt_center_R2 = P - center_R2;
      angleP = atan2(P_wrt_center_R2(1), P_wrt_center_R2(0));
      angle_intersect_point = atan2(intersect_point_wrt_center_R2(1), intersect_point_wrt_center_R2(0));
      // Since atan2 returns values in [-pi, pi], diff2 is in range [-2pi, 2pi]
      diff2 = angleP - angle_intersect_point;

      if (diff2 > 0) {
        diff2 -= 2*M_PI; // Since we are turning right
      }

      if (diff2 > -M_PI) {
        throw std::runtime_error("Error in computing LR path");
      }
    }

    angle_intersect_point = atan2(intersect_point_wrt_center_L(1), intersect_point_wrt_center_L(0));
    double angle_p0 = atan2(-perp_0(1), -perp_0(0));

    // Since atan2 returns values in [-pi, pi], diff1 is in range [-2pi, 2pi]
    double diff1 = angle_intersect_point - angle_p0;
    if (diff1 < 0) {
      diff1 += 2*M_PI; // Since we are turning left
    }

    turns(0, 0) = rho;
    turns(0, 1) = rho*diff1;
    turns(1, 0) = -rho;
    turns(1, 1) = -rho*diff2;    
    return turns;
  } else if (perp_0.dot(P - p_0) >= 0) {
    turns = turns_for_CS_path(x_0, y_0, theta_0, x_f, y_f, rho, true);
  } else {
    turns = turns_for_CS_path(x_0, y_0, theta_0, x_f, y_f, rho, false);
  }

  return turns;
}

RowMatrixXd elongated_dubins_path_one_sided(double x_0, double y_0, double theta_0, double x_f, double y_f, double s, double rho, double tol, bool verbose) {
  RowMatrixXd shortest_turns = turns_for_one_sided_dubins_path(x_0, y_0, theta_0, x_f, y_f, rho);
  double length = shortest_turns.col(1).sum();

  if (std::abs(length - s) < tol) {
    double theta_f = theta_0;
    for (int row = 0; row < shortest_turns.rows(); ++row) {
      if (shortest_turns(row, 0) != 0) {
        theta_f += shortest_turns(row, 1)/shortest_turns(row, 0);
      }
    }
    if (check_elongation_possible(x_0, y_0, theta_0, x_f, y_f, theta_f, s, rho)) {
      return shortest_turns;
    }
  }

  if (length > s) {
    return std::numeric_limits<double>::infinity()*RowMatrixXd::Ones(1, 2);
  }

  // Check if endpoint is in D_I (interior of union of left and right turning circles). If so, we can elongate to arbitrary length
  double c_0 = cos(theta_0);
  double s_0 = sin(theta_0);

  Vector2d p_0(x_0, y_0);
  Vector2d dir_0(c_0, s_0);
  Vector2d perp_0(-s_0, c_0);
  Vector2d center_L = p_0 + rho*perp_0;
  Vector2d P(x_f, y_f);

  // Left case
  if ((P - center_L).norm() < rho) {
    // Apply parallel tangents
    RowMatrixXd turns(5, 2);
    if (shortest_turns(0, 0) != -rho || shortest_turns(1, 0) != rho) {
      throw std::runtime_error("Did not get expected path type for D_I left");
    }
    turns(0, 0) = -rho; // Turn right
    turns(0, 1) = shortest_turns(0, 1);

    double extra_length = s - shortest_turns.col(1).sum();

    turns(1, 0) = 0; // Straight
    turns(1, 1) = extra_length/2;

    turns(2, 0) = rho; // Turn left
    turns(2, 1) = M_PI*rho;

    turns(3, 0) = 0; // Straight
    turns(3, 1) = turns(1, 1);

    turns(4, 0) = rho; // Turn left
    turns(4, 1) = shortest_turns(1, 1) - turns(2, 1);

    double theta_f = theta_0;
    for (int row = 0; row < turns.rows(); ++row) {
      if (turns(row, 0) != 0) {
        theta_f += turns(row, 1)/turns(row, 0);
      }
    }
    if (!check_elongation_possible(x_0, y_0, theta_0, x_f, y_f, theta_f, s, rho)) {
      throw std::runtime_error("RL elongation inconsistent");
    }
    return turns;
  }

  Vector2d center_R = p_0 - rho*perp_0;

  // Right case
  if ((P - center_R).norm() < rho) {
    // Apply parallel tangents
    RowMatrixXd turns(5, 2);
    if (shortest_turns(0, 0) != rho || shortest_turns(1, 0) != -rho) {
      throw std::runtime_error("Did not get expected path type for D_I right");
    }
    turns(0, 0) = rho; // Turn left
    turns(0, 1) = shortest_turns(0, 1);

    double extra_length = s - shortest_turns.col(1).sum();

    turns(1, 0) = 0; // Straight
    turns(1, 1) = extra_length/2;

    turns(2, 0) = -rho; // Turn right
    turns(2, 1) = M_PI*rho;

    turns(3, 0) = 0; // Straight
    turns(3, 1) = turns(1, 1);

    turns(4, 0) = -rho; // Turn right
    turns(4, 1) = shortest_turns(1, 1) - turns(2, 1);

    double theta_f = theta_0;
    for (int row = 0; row < turns.rows(); ++row) {
      if (turns(row, 0) != 0) {
        theta_f += turns(row, 1)/turns(row, 0);
      }
    }
    if (!check_elongation_possible(x_0, y_0, theta_0, x_f, y_f, theta_f, s, rho)) {
      throw std::runtime_error("LR elongation inconsistent");
    }
    return turns;
  }

  // Check if endpoint is in D_II = Omega - int(Omega_+ union Omega_L union Omega_R)
  // We already know the endpoint isn't in the interior of Omega_L union Omega_R (D_I).
  // int(Omega_+ union Omega_L union Omega_R contains points on the boundary of Omega_L union Omega_R, but these points
  // are in the interior of Omega_+ so we don't need to check them separately from checking the interior of Omega_+
  // Omega_+ = Omega_EL intersect Omega_ER intersect {y : y >= 0}
  Vector2d P_wrt_p0 = P - p_0;
  bool left_turn = perp_0.dot(P_wrt_p0) >= 0;

  if (!(Vector2d(c_0, s_0).dot(P_wrt_p0) >= 0 && (P - center_L).norm() < 3*rho && (P - center_R).norm() < 3*rho)) {
    if (left_turn && (shortest_turns(0, 0) != rho || shortest_turns(1, 0) != 0)) {
      throw std::runtime_error("Did not get expected path type for D_II left");
    }

    if (!left_turn && (shortest_turns(0, 0) != -rho || shortest_turns(1, 0) != 0)) {
      throw std::runtime_error("Did not get expected path type for D_II right");
    }

    // Compute pose after turning opposite direction from first segment for angle pi
    double opposite_turn_dir = left_turn ? -1. : 1;
    double theta_after_opposite_turn = left_turn ? theta_0 - M_PI : theta_0 + M_PI;
    Vector2d pos_after_opposite_turn = p_0 + rho*opposite_turn_dir*Vector2d(-s_0 + sin(theta_after_opposite_turn),
                                                                            c_0 - cos(theta_after_opposite_turn));
    RowMatrixXd remaining_turns = turns_for_CS_path(pos_after_opposite_turn(0),
                                                    pos_after_opposite_turn(1),
                                                    theta_after_opposite_turn, x_f, y_f, rho, left_turn);
    double d_I = M_PI*rho + remaining_turns.col(1).sum();
    if (s <= d_I) {
      // Elongate by performing bisection on the opposite turn angle
      double theta_opposite_min = 0;
      double theta_opposite_max = M_PI;

      for (int bisection_iter = 0; bisection_iter < max_bisection_iter; ++bisection_iter) {
        double theta_opposite_mid = 0.5*(theta_opposite_min + theta_opposite_max);

        double theta_after_opposite_turn = left_turn ? theta_0 - theta_opposite_mid : theta_0 + theta_opposite_mid;
        Vector2d pos_after_opposite_turn = p_0 + rho*opposite_turn_dir*Vector2d(-s_0 + sin(theta_after_opposite_turn),
                                                                                c_0 - cos(theta_after_opposite_turn));

        RowMatrixXd remaining_turns = turns_for_CS_path(pos_after_opposite_turn(0),
                                                        pos_after_opposite_turn(1),
                                                        theta_after_opposite_turn, x_f, y_f, rho, left_turn);
        double dist = theta_opposite_mid*rho + remaining_turns.col(1).sum();
        if (std::abs(dist - s) < tol) {
          RowMatrixXd turns(3, 2);
          turns(0, 0) = left_turn ? -rho : rho; // If left turn, opposite turn is right turn, and vice versa
          turns(0, 1) = rho*theta_opposite_mid;
          turns.bottomRows<2>() = remaining_turns.bottomRows<2>();

          double theta_f = theta_0;
          for (int row = 0; row < turns.rows(); ++row) {
            if (turns(row, 0) != 0) {
              theta_f += turns(row, 1)/turns(row, 0);
            }
          }
          if (check_elongation_possible(x_0, y_0, theta_0, x_f, y_f, theta_f, s, rho)) {
            return turns;
          }
        }
        if (dist > s) {
          theta_opposite_max = theta_opposite_mid;
        } else {
          theta_opposite_min = theta_opposite_mid;
        }
      }
      throw std::runtime_error("Bisection s <= d_I failed");
    } else {
      // Elongate by performing bisection on the middle turn radius
      double rho_min = rho;
      double rho_max = rho*2;
      bool found_ub = false;
      for (int bisection_iter = 0; bisection_iter < max_bisection_iter; ++bisection_iter) {
        RowMatrixXd remaining_turns = turns_for_CS_path(pos_after_opposite_turn(0),
                                                        pos_after_opposite_turn(1),
                                                        theta_after_opposite_turn, x_f, y_f, rho_max, left_turn);
        if (M_PI*rho + remaining_turns.col(1).sum() > s) {
          found_ub = true;
          break;
        }
        rho_max *= 2;
      }
      if (!found_ub) {
        throw std::runtime_error("Bisection s > d_I upper bound search failed");
      }

      for (int bisection_iter = 0; bisection_iter < max_bisection_iter; ++bisection_iter) {
        double rho_mid = 0.5*(rho_min + rho_max);
        RowMatrixXd remaining_turns = turns_for_CS_path(pos_after_opposite_turn(0),
                                                        pos_after_opposite_turn(1),
                                                        theta_after_opposite_turn, x_f, y_f, rho_mid, left_turn);
        double dist = M_PI*rho + remaining_turns.col(1).sum();
        if (std::abs(dist - s) < tol) {
          RowMatrixXd turns(3, 2);
          turns(0, 0) = left_turn ? -rho : rho; // If left turn, opposite turn is right turn, and vice versa
          turns(0, 1) = rho*M_PI;
          turns.bottomRows<2>() = remaining_turns.bottomRows<2>();

          double theta_f = theta_0;
          for (int row = 0; row < turns.rows(); ++row) {
            if (turns(row, 0) != 0) {
              theta_f += turns(row, 1)/turns(row, 0);
            }
          }
          if (check_elongation_possible(x_0, y_0, theta_0, x_f, y_f, theta_f, s, rho)) {
            return turns;
          }
        }
        if (dist > s) {
          rho_max = rho_mid;
        } else {
          rho_min = rho_mid;
        }
      }
      throw std::runtime_error("Bisection s > d_I failed");
    }
  }

  // Now we're in D_III. Handle LS paths via mirroring
  Vector2d P_original = P;
  if (left_turn) {
    P = P - 2*perp_0*perp_0.dot(P - p_0);
    x_f = P(0);
    y_f = P(1);
  }

  double delta_x = x_f - x_0;
  double delta_y = y_f - y_0;
  double r_M = (delta_x*delta_x + delta_y*delta_y)/(2*std::abs(s_0*delta_x - c_0*delta_y));

  // Find length of R^{rM} path
  Vector2d center_R_rM = p_0 - r_M*perp_0;
  Vector2d p_0_wrt_center = p_0 - center_R_rM;

  double p_0_angle = atan2(p_0_wrt_center(1), p_0_wrt_center(0));

  Vector2d P_wrt_center = P - center_R_rM;

  double P_angle = atan2(P_wrt_center(1), P_wrt_center(0));

  // Since atan2 returns values in [-pi, pi], diff is in range [-2pi, 2pi]
  double diff = P_angle - p_0_angle;
  if (diff > 0) {
    diff = diff - 2*M_PI; // Since we're turning right
  }
  diff = -diff; // Since we're turning right
  double L_rM = diff*r_M;
  if (s <= L_rM) {
    // Elongate by performing bisection on the turn radius
    double rho_min = rho;
    double rho_max = r_M;
    for (int bisection_iter = 0; bisection_iter < max_bisection_iter; ++bisection_iter) {
      double rho_mid = 0.5*(rho_min + rho_max);
      RowMatrixXd turns = turns_for_CS_path(x_0,
                                            y_0,
                                            theta_0, x_f, y_f, rho_mid, false);
      double dist = turns.col(1).sum();
      if (std::abs(dist - s) < tol) {
        if (left_turn) {
          turns(0, 0) = -turns(0, 0);
        }

        double theta_f = theta_0;
        for (int row = 0; row < turns.rows(); ++row) {
          if (turns(row, 0) != 0) {
            theta_f += turns(row, 1)/turns(row, 0);
          }
        }
        if (check_elongation_possible(x_0, y_0, theta_0, P_original(0), P_original(1), theta_f, s, rho)) {
          return turns;
        }
      }
      if (dist > s) {
        rho_max = rho_mid;
      } else {
        rho_min = rho_mid;
      }
    }
    throw std::runtime_error("Bisection s <= L_rM failed");
  }

  Vector2d P_wrt_L = P - center_L;
  double xi = acos(P_wrt_L.dot(-perp_0)/P_wrt_L.norm());

  // Ding 2019 paper defines rho as distance from turning circle center and P
  double paper_rho = (P - center_L).norm();
  double paper_rho2 = paper_rho*paper_rho;
  double rho2 = rho*rho;

  double acos1 = acos((paper_rho2 + 3*rho2)/(4*paper_rho*rho));
  double acos2 = acos((5*rho2 - paper_rho2)/(4*rho2));

  double length_lambda_minus = rho*(xi - acos1 + acos2);
  double length_beta_minus = rho*(xi + acos1 + 2*M_PI - acos2);

  if (length_lambda_minus < s && s < length_beta_minus) {
    // Elongation impossible
    return std::numeric_limits<double>::infinity()*RowMatrixXd::Ones(1, 2);
  }

  if (s <= length_lambda_minus) {
    // Elongate by performing bisection on the turn radius
    double rho_min = rho;
    double rho_max = r_M;
    RowMatrixXd turns(2, 2);
    for (int bisection_iter = 0; bisection_iter < max_bisection_iter; ++bisection_iter) {
      double rho_mid = 0.5*(rho_min + rho_max);
      RowMatrixXd turns = turns_for_LR_path(x_0,
                                            y_0,
                                            theta_0, x_f, y_f, rho, rho_mid, xi, false);
      double dist = turns.col(1).sum();
      if (std::abs(dist - s) < tol) {
        if (left_turn) {
          turns.col(0) = -turns.col(0);
        }

        double theta_f = theta_0;
        for (int row = 0; row < turns.rows(); ++row) {
          if (turns(row, 0) != 0) {
            theta_f += turns(row, 1)/turns(row, 0);
          }
        }
        if (check_elongation_possible(x_0, y_0, theta_0, P_original(0), P_original(1), theta_f, s, rho)) {
          return turns;
        }
      }
      // This is reversed because decreasing rho increases path length
      if (dist > s) {
        rho_min = rho_mid;
      } else {
        rho_max = rho_mid;
      }
    }
    throw std::runtime_error("Bisection s <= length_lambda_minus failed");
  }

  // s >= length_beta_minus
  // Elongate by performing bisection on the turn radius
  double rho_min = rho;
  double rho_max = rho*2;
  bool found_ub = false;
  for (int bisection_iter = 0; bisection_iter < max_bisection_iter; ++bisection_iter) {
    RowMatrixXd turns = turns_for_LR_path(x_0,
                                          y_0,
                                          theta_0, x_f, y_f, rho, rho_max, xi, true);
    if (turns.col(1).sum() > s) {
      found_ub = true;
      break;
    }
    rho_max *= 2;
  }
  if (!found_ub) {
    throw std::runtime_error("Bisection upper bound search failed");
  }

  for (int bisection_iter = 0; bisection_iter < max_bisection_iter; ++bisection_iter) {
    double rho_mid = 0.5*(rho_min + rho_max);
    RowMatrixXd turns = turns_for_LR_path(x_0,
                                          y_0,
                                          theta_0, x_f, y_f, rho, rho_mid, xi, true);
    double dist = turns.col(1).sum();
    if (std::abs(dist - s) < tol) {
      if (left_turn) {
        turns.col(0) = -turns.col(0);
      }
      double theta_f = theta_0;
      // double x = x_0;
      // double y = y_0;
      for (int row = 0; row < turns.rows(); ++row) {
        if (turns(row, 0) != 0) {
          // double theta_f_before = theta_f;
          theta_f += turns(row, 1)/turns(row, 0);
          /*
          x += turns(row, 0)*(-sin(theta_f_before) + sin(theta_f));
          y += turns(row, 0)*(cos(theta_f_before) - cos(theta_f));
          */
        } else {
          throw std::runtime_error("Supposed to be CC path");
        }
      }
      if (check_elongation_possible(x_0, y_0, theta_0, P_original(0), P_original(1), theta_f, s, rho)) {
        return turns;
      }
      /*
      Dubins path(AngleInterval(Point(x_0, y_0), theta_0, 0), AngleInterval(Point(P_original(0), P_original(1)), theta_f, 1e-4), rho);
      // std::cout << turns_for_dubins_path(x_0, y_0, theta_0, x_f, y_f, theta_f, rho).col(1).sum() << " " << s << " " << dist << std::endl;
      std::cout << path.getLength() << " " << s << " " << dist << " " << x << " " << y << " " << P_original(0) << " " << P_original(1) << " " << left_turn << std::endl;
      */
    }
    if (dist > s) {
      rho_max = rho_mid;
    } else {
      rho_min = rho_mid;
    }
  }
  
  // std::cout << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << x_0 << " " << y_0 << " " << theta_0 << " " << x_f << " " << y_f << " " << s << std::endl;

  throw std::runtime_error("Bisection s >= length_beta_minus failed");
}

RowMatrixXd elongated_dubins_path_one_sided(double x_0, double y_0, double theta_0, double x_f, double y_f, double s, double rho) {
  return elongated_dubins_path_one_sided(x_0, y_0, theta_0, x_f, y_f, s, rho, 1e-4, false);
}

RowMatrixXd elongated_dubins_path_one_sided(double x_0, double y_0, double theta_0, double x_f, double y_f, double s, double rho, double tol) {
  return elongated_dubins_path_one_sided(x_0, y_0, theta_0, x_f, y_f, s, rho, tol, false);
}
