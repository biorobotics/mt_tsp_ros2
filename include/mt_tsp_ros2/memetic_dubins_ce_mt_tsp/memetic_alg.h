#pragma once
#include <chrono>
#include <omp.h>
#include <Eigen/Dense>
#include <random>
#include <set>
#include "mt_tsp_ros2/cpp_spline.h"
#include "mt_tsp_ros2/extended_cpp_spline.h"
#include "mt_tsp_ros2/elongate_one_sided_dubins_path.h"
#include "mt_tsp_ros2/memetic_dubins_ce_mt_tsp/memetic_alg_params.h"
#include <iomanip>

using namespace Eigen;

typedef Matrix<double, Dynamic, Dynamic, RowMajor> RowMatrixXd;

typedef Matrix<bool, Dynamic, 1> VectorXb;

typedef Matrix<double, 1, 1> Vector1d;

const int gene_size = 3; // target index, theta, and delta t

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values
  std::stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

bool repair_chromosome(Ref<MatrixXd> X, const Ref<const RowMatrixXd> &tw_per_target, const Ref<const VectorXd> &target_radii, const std::vector<ExtendedCppSpline> &q_trj_per_target, const Ref<const Vector2d> &p0, double heading0, double vmax, double &cost, bool dubins, double rho, int &max_newton_iter_for_success_repair, const MemeticAlgParams &params, bool no_tw, double &final_heading, double t0) {
  int num_targets = tw_per_target.rows();

  bool repair_failed = false;
  double t = t0;
  Vector2d pos = p0;
  double heading = heading0;
  Vector2d next_rel_pos;
  Vector2d next_pos;
  cost = 0;
  for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
    int target_idx = X(seq_idx, 0);
    double theta = X(seq_idx, 1);
    double delta_t = X(seq_idx, 2);
    double next_t = t + delta_t;

    if (!no_tw) {
      if (next_t > tw_per_target(target_idx, 1)) {
        next_t = tw_per_target(target_idx, 1);
        delta_t = next_t - t;
        if (delta_t < 0) {
          repair_failed = true;
          break;
        }
        X(seq_idx, 2) = delta_t;
      } else if (next_t < tw_per_target(target_idx, 0)) {
        next_t = tw_per_target(target_idx, 0);
        delta_t = next_t - t;
        X(seq_idx, 2) = delta_t;
      }
    }

    next_rel_pos = target_radii[target_idx]*Vector2d(cos(theta), sin(theta));

    // Check if interception is feasible
    if (dubins) {
      // Run Newton to restore feasibility if needed
      int max_newton_iter = 2776; // On one of the 10 target instances I ran, we needed at most 1388 iterations for successful transformation so I'm using 2x that number to declare failure
      double delta_delta_t_finite_diff = 1e-4;
      double repair_tol = 1e-4;
      bool newton_succeeded = false;
      // std::cout << "starting newton" << std::endl;
      for (int newton_iter = 0; newton_iter < max_newton_iter; ++newton_iter) {
        next_pos = q_trj_per_target[target_idx](next_t) + next_rel_pos;

        // std::cout << delta_t << " desired length " << vmax*delta_t << " shortest path length " << turns_for_one_sided_dubins_path(pos(0), pos(1), heading, next_pos(0), next_pos(1), rho).col(1).sum() << std::endl;
        RowMatrixXd turns = elongated_dubins_path_one_sided(pos(0), pos(1), heading, next_pos(0), next_pos(1), vmax*delta_t, rho, 1e-2);
        if (std::isfinite(turns(0, 0))) {
          max_newton_iter_for_success_repair = std::max(newton_iter, max_newton_iter_for_success_repair);
          newton_succeeded = true;
          for (int row = 0; row < turns.rows(); ++row) {
            if (turns(row, 0) != 0) {
              heading += turns(row, 1)/turns(row, 0);
            }
          }
          break;
        }

        RowMatrixXd shortest_path_turns = turns_for_one_sided_dubins_path(pos(0), pos(1), heading, next_pos(0), next_pos(1), rho);
        double shortest_path_dist = shortest_path_turns.col(1).sum();
        double c = shortest_path_dist - vmax*delta_t;

        // Finite-diff
        double delta_t_plus = delta_t + delta_delta_t_finite_diff;

        double next_t_plus = t + delta_t_plus;

        Vector2d next_pos_plus = q_trj_per_target[target_idx](next_t_plus) + next_rel_pos;

        RowMatrixXd shortest_path_turns_plus = turns_for_one_sided_dubins_path(pos(0), pos(1), heading, next_pos_plus(0), next_pos_plus(1), rho);
        double shortest_path_dist_plus = shortest_path_turns_plus.col(1).sum();
        double c_plus = shortest_path_dist_plus - vmax*delta_t_plus;

        double deriv = (c_plus - c)/delta_delta_t_finite_diff;

        delta_t -= params.repair_step_size*c/deriv;

        next_t = t + delta_t;

        if (no_tw) {
          if (next_t < t) {
            next_t = t;
            delta_t = 0.;
          }
        } else {
          if (next_t > tw_per_target(target_idx, 1)) {
            next_t = tw_per_target(target_idx, 1);
            delta_t = next_t - t;
            if (delta_t < 0) {
              throw std::runtime_error("Projection of delta t failed in repair");
            }
          } else if (next_t < tw_per_target(target_idx, 0)) {
            next_t = tw_per_target(target_idx, 0);
            delta_t = next_t - t;
          }
        }
        X(seq_idx, 2) = delta_t;
      }

      if (!newton_succeeded) {
        repair_failed = true;
        break;
      }

      if (no_tw && seq_idx == num_targets - 1) {
        cost = next_t;
      } else {
        cost += next_t - tw_per_target(target_idx, 0);
      }
    } else {
      next_pos = q_trj_per_target[target_idx](next_t) + next_rel_pos;
      double dist = (next_pos - pos).norm();
      if (dist <= vmax*delta_t) {
        // Travel is feasible, no need to repair
        cost += dist;
        t = next_t;
        pos = next_pos;
        continue;
      }

      double t_high;
      double t_high_dist;

      if (no_tw) {
        bool found_ub = false;
        t_high = 2*t;
        delta_t = t_high - t;
        for (int i = 0; i < 100; ++i) {
          // Check if travel is feasible to next_pos
          next_pos = q_trj_per_target[target_idx](t_high) + next_rel_pos;
          dist = (next_pos - pos).norm();
          if (dist <= vmax*delta_t) {
            found_ub = true;
            t_high_dist = dist;
            break;
          }
          t_high *= 2;
          delta_t = t_high - t;
        }
        if (!found_ub) {
          throw std::runtime_error("Did not find upper bound for bisection");
        }
      } else {
        // Check if travel is feasible to next_rel_pos at end of time window
        next_pos = q_trj_per_target[target_idx](tw_per_target(target_idx, 1)) + next_rel_pos;
        delta_t = tw_per_target(target_idx, 1) - t;
        dist = (next_pos - pos).norm();
        if (dist > vmax*delta_t) {
          // Travel is infeasible even to end of time window
          repair_failed = true;
          break;
        }
        t_high = tw_per_target(target_idx, 1);
        t_high_dist = dist;
      }

      // Run bisection to find earliest time such that interception is feasible
      double t_low = next_t;
      int num_bisection_iter = 10;
      for (int bisection_iter = 0; bisection_iter < num_bisection_iter; ++bisection_iter) {
        double t_mid = 0.5*(t_low + t_high);
        delta_t = t_mid - t;
        Vector2d pos_mid = q_trj_per_target[target_idx](t_mid) + next_rel_pos;
        dist = (pos_mid - pos).norm();
        if (dist > vmax*delta_t) {
          // Travel is infeasible
          t_low = t_mid;
        } else {
          // Travel is feasible
          t_high = t_mid;
          t_high_dist = dist;
          next_pos = pos_mid;
        }
      }
      next_t = t_high;

      delta_t = t_high - t;
      X(seq_idx, 2) = delta_t;
      cost += t_high_dist;
    }

    // Update cost, time, and position
    t = next_t;
    pos = next_pos;
  }
  if (repair_failed) {
    if (no_tw && !dubins) {
      throw std::runtime_error("Repair failed and no time windows or min turning radius");
    }
    cost = std::numeric_limits<double>::infinity();
  }
  final_heading = heading;

  /*
  if (!dubins && !repair_failed) {
    std::cout << "Checking after repair" << std::endl;
    double t = t0;
    Vector2d pos = p0;
    for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
      int target_idx = X(seq_idx, 0);
      double theta = X(seq_idx, 1);
      double delta_t = X(seq_idx, 2);
      t += delta_t;
      Vector2d next_pos = q_trj_per_target[target_idx](t) + target_radii[target_idx]*Vector2d(cos(theta), sin(theta));
      std::cout << t << " " << next_pos.transpose() << std::endl;
      if ((next_pos - pos).norm() > vmax*delta_t) {
        std::cout << (next_pos - pos).norm() - vmax*delta_t << std::endl;
        throw std::runtime_error("Speed constraint violated after repair");
      }
      pos = next_pos;

      if (!no_tw && (t < tw_per_target(target_idx, 0) - 1e-4 || t > tw_per_target(target_idx, 1) + 1e-4)) {
        throw std::runtime_error("t out of window");
      }
    }
  }
  */

  return repair_failed;
}

bool repair_chromosome(Ref<MatrixXd> X, const Ref<const RowMatrixXd> &tw_per_target, const Ref<const VectorXd> &target_radii, const std::vector<ExtendedCppSpline> &q_trj_per_target, const Ref<const Vector2d> &p0, double heading0, double vmax, Ref<Vector1d> cost, bool dubins, double rho, int &max_newton_iter_for_success_repair, const MemeticAlgParams &params, bool no_tw) {
  double tmp_cost = 0.;
  double final_heading;
  bool repair_failed = repair_chromosome(X, tw_per_target, target_radii, q_trj_per_target, p0, heading0, vmax, tmp_cost, dubins, rho, max_newton_iter_for_success_repair, params, no_tw, final_heading, 0.);
  cost(0) = tmp_cost;
  return repair_failed;
}

// next_heading is Ref<Vector1d> rather than double& so I can test in python
double find_earliest_arrival_time_dubins(Ref<Vector2d> next_pos, Ref<Vector1d> next_heading, const ExtendedCppSpline &q_trj, const Ref<const Vector2d> &pos, double heading, double tw_start, double tw_end, double radius, double theta, double t, double vmax_agent, double rho, double vmax_target, int &max_bisection_iter_for_success_transformation, bool no_tw) {
  Vector2d next_rel_pos = radius*Vector2d(cos(theta), sin(theta));

  Vector2d pos_t = q_trj(t) + next_rel_pos;

  RowMatrixXd turns = turns_for_one_sided_dubins_path(pos(0), pos(1), heading, pos_t(0), pos_t(1), rho);
  double length = turns.col(1).sum();

  double t_min = t + length/(vmax_agent + vmax_target);
  double t_max = t + length/(vmax_agent - vmax_target);
  // double t_min = tw_start;
  // double t_max = tw_end;

  if (!no_tw) {
    t_min = std::max(t_min, tw_start);
    t_max = std::min(t_max, tw_end);
  }

  Vector2d pos_min = q_trj(t_min) + next_rel_pos;
  turns = turns_for_one_sided_dubins_path(pos(0), pos(1), heading, pos_min(0), pos_min(1), rho);
  bool CC_detected = turns(1, 0) != 0;
  length = turns.col(1).sum();
  double delta_min = length - vmax_agent*(t_min - t);
  bool min_CC = turns(1, 0) != 0;
  double length_min = length;

  Vector2d pos_max = q_trj(t_max) + next_rel_pos;
  turns = turns_for_one_sided_dubins_path(pos(0), pos(1), heading, pos_max(0), pos_max(1), rho);
  CC_detected |= turns(1, 0) != 0;
  length = turns.col(1).sum();
  double delta_max = length - vmax_agent*(t_max - t);
  bool max_CC = turns(1, 0) != 0;
  double length_max = length;
  if (delta_min*delta_max > 0) {
    // throw std::runtime_error("No bracket"); 
    return std::numeric_limits<double>::infinity();
  }

  int max_bisection_iter = no_tw ? 100 : 34; // On one of the 10 target instances I ran, we needed at most 17 iterations for successful transformation so I'm using 2x that number to declare failure. I just put the 100 in here just in case for the no_tw option
  double bisection_tol = 1e-4;
  for (int bisection_iter = 0; bisection_iter < max_bisection_iter; ++bisection_iter) {
    double t_mid = 0.5*(t_min + t_max);
    next_pos = q_trj(t_mid) + next_rel_pos;
    turns = turns_for_one_sided_dubins_path(pos(0), pos(1), heading, next_pos(0), next_pos(1), rho);
    CC_detected |= turns(1, 0) != 0;
    double length = turns.col(1).sum();
    double delta = length - vmax_agent*(t_mid - t);
    if (std::abs(delta) < bisection_tol) {
      next_heading(0) = heading;
      for (int row = 0; row < turns.rows(); ++row) {
        if (turns(row, 0) != 0) {
          next_heading(0) += turns(row, 1)/turns(row, 0);
        }
      }
      max_bisection_iter_for_success_transformation = std::max(max_bisection_iter_for_success_transformation, bisection_iter);
      return t_mid;
    }
    if (delta*delta_min > 0) {
      t_min = t_mid;
      min_CC = turns(1, 0) != 0;
      length_min = length;
    } else {
      t_max = t_mid;
      max_CC = turns(1, 0) != 0;
      length_max = length;
    }
  }
  // throw std::runtime_error("Transformation method bisection ran out of iterations"); 
  if (!CC_detected) {
    throw std::runtime_error("Transformation bisection failed and no CC path detected");
  }
  return std::numeric_limits<double>::infinity();
}

bool transform_chromosome(Ref<MatrixXd> X, const Ref<const RowMatrixXd> &tw_per_target, const Ref<const VectorXd> &target_radii, const std::vector<ExtendedCppSpline> &q_trj_per_target, const Ref<const Vector2d> &p0, double heading0, double vmax, double rho, double &cost, const Ref<const VectorXd> &speed_upper_bounds, int &max_newton_iter_for_success_repair, int &max_bisection_iter_for_success_transformation, const MemeticAlgParams &params, bool no_tw) {
  int num_targets = tw_per_target.rows();

  double t = 0;
  Vector2d pos = p0;
  double heading = heading0;
  Vector2d next_rel_pos;
  Vector2d next_pos;
  Vector1d next_heading;
  cost = 0;
  for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
    int target_idx = X(seq_idx, 0);

    // Changing previous delta_t values may make next interception infeasible. If so, repair
    double tmp_cost;
    std::vector<ExtendedCppSpline> tmp_q_trj;
    tmp_q_trj.push_back(q_trj_per_target[target_idx]);
    MatrixXd tmp_chromosome = X.block(seq_idx, 0, 1, gene_size);
    tmp_chromosome(0, 0) = 0;
    double tmp_next_heading;
    bool repair_failed = repair_chromosome(tmp_chromosome, tw_per_target.block(target_idx, 0, 1, 2), target_radii.segment(target_idx, 1), tmp_q_trj, pos, heading, vmax, tmp_cost, true, rho, max_newton_iter_for_success_repair, params, no_tw, tmp_next_heading, t);
    if (repair_failed) {
      cost = std::numeric_limits<double>::infinity();
      return false;
    }
    next_heading(0) = tmp_next_heading;
    X.block(seq_idx, 0, 1, gene_size) = tmp_chromosome;
    X(seq_idx, 0) = target_idx;

    double theta = X(seq_idx, 1);
    double delta_t = X(seq_idx, 2);
    double next_t = t + delta_t;

    if (seq_idx != num_targets - 1) {
      // Check the encounter pattern to the target at seq_idx + 1. If catchup, then optimize the arrival time to target at seq_idx.
      // If meeting, don't optimize the arrival time. The following method of determining the encounter pattern was obtained
      // from correspondence with the authors
      int following_target_idx = X(seq_idx + 1, 0);
      double following_theta = X(seq_idx + 1, 1);
      double following_delta_t = X(seq_idx + 1, 2);
      Vector2d following_rel_pos = target_radii[following_target_idx]*Vector2d(cos(following_theta), sin(following_theta));
      Vector2d following_pos = q_trj_per_target[following_target_idx](next_t + following_delta_t) + following_rel_pos;

      RowMatrixXd turns = turns_for_one_sided_dubins_path(p0(0), p0(1), heading0, following_pos(0), following_pos(1), rho);

      double length = turns.col(1).sum();

      Vector2d following_pos_plus = q_trj_per_target[following_target_idx](next_t + following_delta_t + 0.1) + following_rel_pos;
      turns = turns_for_one_sided_dubins_path(p0(0), p0(1), heading0, following_pos_plus(0), following_pos_plus(1), rho);

      double length_plus = turns.col(1).sum();

      if (length_plus < length) {
        next_rel_pos = target_radii[target_idx]*Vector2d(cos(theta), sin(theta));
        next_pos = q_trj_per_target[target_idx](next_t);

        // Meeting pattern
        if (!no_tw) {
          cost += next_t - tw_per_target(target_idx, 0);
        }
        t = next_t;
        pos = next_pos;
        heading = next_heading(0);

        continue;
      }
      // Catch-up pattern
    }

    double next_t_tmp = find_earliest_arrival_time_dubins(next_pos, next_heading, q_trj_per_target[target_idx], pos, heading, tw_per_target(target_idx, 0), tw_per_target(target_idx, 1), target_radii(target_idx), theta, t, vmax, rho, speed_upper_bounds(target_idx), max_bisection_iter_for_success_transformation, no_tw);
    if (std::isfinite(next_t_tmp)) {
      next_t = next_t_tmp;
      delta_t = next_t - t;
    } else {
      cost = std::numeric_limits<double>::infinity();
      return false;
      // Comment the above two lines if we want to continue the transformation on subsequent targets even if bisection failed.
      // However, I tried this on a 50 target instance and it brought the final cost from 522.868936 to 789.149723
      next_pos = q_trj_per_target[target_idx](next_t);
      next_heading(0) = tmp_next_heading;
    }
    X(seq_idx, 2) = delta_t;

    // Update cost, time, and position
    // Assume dubins
    if (no_tw && seq_idx == num_targets - 1) {
      cost = next_t;
    } else {
      cost += next_t - tw_per_target(target_idx, 0);
    }
    t = next_t;
    pos = next_pos;
    heading = next_heading(0);
  }
  return true;
}

bool transform_chromosome_no_dubins(Ref<MatrixXd> X, const Ref<const RowMatrixXd> &tw_per_target, const Ref<const VectorXd> &target_radii, const std::vector<ExtendedCppSpline> &q_trj_per_target, const Ref<const Vector2d> &p0, double vmax, double &cost, int &max_newton_iter_for_success_repair, int &max_bisection_iter_for_success_transformation, const MemeticAlgParams &params, bool no_tw) {
  int num_targets = tw_per_target.rows();

  double t = 0;
  Vector2d pos = p0;
  Vector2d next_rel_pos;
  Vector2d next_pos;
  cost = 0;
  bool made_change = false;
  for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
    int target_idx = X(seq_idx, 0);

    if (seq_idx != 0) {
      // Changing previous delta_t values may make next interception infeasible. If so, repair
      double tmp_cost;
      std::vector<ExtendedCppSpline> tmp_q_trj;
      tmp_q_trj.push_back(q_trj_per_target[target_idx]);
      MatrixXd tmp_chromosome = X.block(seq_idx, 0, 1, gene_size);
      tmp_chromosome(0, 0) = 0;
      double dummy_next_heading;
      bool repair_failed = repair_chromosome(tmp_chromosome, tw_per_target.block(target_idx, 0, 1, 2), target_radii.segment(target_idx, 1), tmp_q_trj, pos, 0., vmax, tmp_cost, false, 0., max_newton_iter_for_success_repair, params, no_tw, dummy_next_heading, t);
      if (repair_failed) {
        cost = std::numeric_limits<double>::infinity();
        return false;
      }
      X.block(seq_idx, 0, 1, gene_size) = tmp_chromosome;
      X(seq_idx, 0) = target_idx;
    }

    double theta = X(seq_idx, 1);
    double delta_t = X(seq_idx, 2);

    Vector2d next_rel_pos = target_radii[target_idx]*Vector2d(cos(theta), sin(theta));

    int max_gd_iter = 10;
    int max_backtrack_iter = 3;
    double next_t = t + delta_t;
    for (int gd_iter = 0; gd_iter < max_gd_iter; ++gd_iter) {
      next_pos = q_trj_per_target[target_idx](next_t) + next_rel_pos;
      double dist = (next_pos - pos).norm();
      double deriv = (next_pos - pos).dot(q_trj_per_target[target_idx].derivatives(next_t));
      // Limit step size to avoid going outside time window
      double step;
      if (no_tw) {
        step = std::min(std::max(-deriv, tw_per_target(target_idx, 0) - next_t), tw_per_target(target_idx, 1) - next_t);
      } else {
        step = -deriv;
      }
      double step_size = step/(-deriv);

      double b = 0.01; // From Zac's class (he mentioned to set b between 1e-4 and 0.1)
      double c = 0.5; // From Zac's class
      bool reduction = false;
      for (int backtrack_iter = 0; backtrack_iter < max_backtrack_iter; ++backtrack_iter) {
        double next_t_cand = next_t - step_size*deriv;
        Vector2d next_pos_cand = q_trj_per_target[target_idx](next_t_cand) + next_rel_pos;
        double dist_cand = (next_pos_cand - pos).norm();
        // Armijo rule.
        // Mutliply change in next_t by derivative to get expected change in cost.
        // We're checking if the actual cost reduction is at least b times the expected
        if (dist_cand - dist <= b*deriv*step) {
          next_t = next_t_cand;
          next_pos = next_pos_cand;
          reduction = true;
          break;
        }
        step_size *= c;
      }
      if (!reduction) {
        break;
      }
    }

    delta_t = next_t - t;
    X(seq_idx, 2) = delta_t;

    // Update cost, time, and position
    // Assume dubins
    cost += (next_pos - pos).norm();
    t = next_t;
    pos = next_pos;
  }
  return true;
}

void check_chromosome_feasible(const Ref<const MatrixXd> &chromosome, const Ref<const RowMatrixXd> &tw_per_target, const Ref<const Vector2d> &p0, double vmax, const std::vector<ExtendedCppSpline> &q_trj_per_target, const Ref<const VectorXd> &target_radii) {
  int num_targets = tw_per_target.rows();
  double t = 0.;
  Vector2d pos = p0;
  for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
    int target_idx = chromosome(seq_idx, 0);
    double theta = chromosome(seq_idx, 1);
    double delta_t = chromosome(seq_idx, 2);
    t += delta_t;

    if (t < tw_per_target(target_idx, 0) || t > tw_per_target(target_idx, 1)) {
      throw std::runtime_error("chromosome infeasible");
    }

    Vector2d next_rel_pos = target_radii[target_idx]*Vector2d(cos(theta), sin(theta));
    Vector2d next_pos = q_trj_per_target[target_idx](t) + next_rel_pos;
    double dist = (next_pos - pos).norm();
    if (dist > vmax*delta_t + 1e-4) {
      throw std::runtime_error("chromosome infeasible");
    }
    pos = next_pos;
  }
}

void get_speed_upper_bounds(Ref<VectorXd> upper_bounds, const std::vector<ExtendedCppSpline> &q_trj_per_target, const Ref<const RowMatrixXd> &tw_per_target) {
  upper_bounds.setConstant(2.);
  /*
  #pragma omp parallel for
  for (int target_idx = 0; target_idx < tw_per_target.rows(); ++target_idx) {
    int num_sample = 100;
    for (int sample_idx = 0; sample_idx < num_sample; ++sample_idx) {
      double alpha = sample_idx/(double)num_sample;
      double t = alpha*tw_per_target(target_idx, 0) + (1 - alpha)*tw_per_target(target_idx, 1);
      Vector2d deriv = q_trj_per_target[target_idx].derivatives(t);
      if (deriv.norm() > upper_bounds(target_idx)) {
        throw std::runtime_error("Upper bound incorrect");
      }
    }
  }
  */
}

// initial_population should have number of rows equal to num_targets*pop_size, and gene_size columns
RowMatrixXd memetic_alg(Ref<RowMatrixXd> selected_pts_per_target, const Ref<const RowMatrixXd> &initial_population, const Ref<const VectorXd> &initial_costs, const std::vector<ExtendedCppSpline> &q_trj_per_target_python, const Ref<const RowMatrixXd> &tw_per_target, const Ref<const VectorXd> &target_radii, double time_limit, double rho, double vmax, const Ref<const Vector2d> &p0, double heading0, int num_openmp_threads, std::vector<RowMatrixXd> &turns_chain, const MemeticAlgParams &params, Ref<Matrix<long, 1, 1>> num_feas_final, bool no_tw) {
  std::vector<std::pair<double, double>> cost_vs_time;
  
  auto timer_start = std::chrono::high_resolution_clock::now();

  double min_cost_record_time = 0.;

  omp_set_num_threads(num_openmp_threads);

  int num_targets = tw_per_target.rows();

  // I'm doing this regardless of whether there are time windows because I'm worried about the GIL
  std::vector<ExtendedCppSpline> q_trj_per_target;
  for (int target_idx = 0; target_idx < num_targets; ++target_idx) {
    q_trj_per_target.push_back(ExtendedCppSpline(q_trj_per_target_python[target_idx].get_knots(), q_trj_per_target_python[target_idx].get_coeffs(), tw_per_target(target_idx, 0), tw_per_target(target_idx, 1)));
  }

  bool dubins = rho != 0;

  VectorXd speed_upper_bounds(1);
  if (dubins) {
    speed_upper_bounds = VectorXd(num_targets);
    // Whether or not time windows are a constraint in the problem, the way the instances are generated,
    // we have targets that move along a straight line, and then a cubic B-spline, then a straight line,
    // and the time windows denote where they have each motion pattern
    get_speed_upper_bounds(speed_upper_bounds, q_trj_per_target, tw_per_target);

    if (speed_upper_bounds.maxCoeff() > vmax) {
      throw std::runtime_error("Max speed of some target is larger than agent speed");
    }
  }

  int pop_size = initial_costs.size();
  if (pop_size != initial_population.rows()/num_targets) {
    throw std::runtime_error("Population size does not match the number of provided cost values");
  }

  std::vector<MatrixXd> population1(pop_size);
  std::vector<double> population_costs1(pop_size);
  Map<VectorXd>(population_costs1.data(), pop_size) = initial_costs;

  std::vector<unsigned int> seeds{4223100027, 870236586, 1683737518, 3430707182, 1613429085, 1714341085, 853110547, 1988005940, 2629786018, 1139192408};
  std::vector<std::mt19937> rngs_per_thread;
  if (num_openmp_threads > seeds.size()) {
    throw std::runtime_error("Too many threads, not enough stored random seeds");
  }
  for (int thread_idx = 0; thread_idx < num_openmp_threads; ++thread_idx) {
    rngs_per_thread.push_back(std::mt19937(seeds[thread_idx]));
  }
  std::uniform_int_distribution<int> parent_distribution(0, pop_size - 1);
  std::uniform_int_distribution<int> crossover_distribution(0, 1);
  std::uniform_real_distribution<double> mutation_distribution(0, 1);

  std::uniform_int_distribution<int> mutation_operator1_distribution(0, num_targets - 1);

  std::uniform_int_distribution<int> mutation_operator2_seq_idx_distribution(0, num_targets - 1);
  std::uniform_real_distribution<double> mutation_operator2_angle_distribution(0, 2*M_PI);

  std::uniform_int_distribution<int> mutation_operator3_seq_idx_distribution(0, num_targets - 1);
  std::uniform_real_distribution<double> mutation_operator3_delta_t_distribution(0, 1);

  std::uniform_int_distribution<int> local_search_elite_distribution(0, pop_size/2);
  std::uniform_int_distribution<int> local_search_gene_idx_distribution(0, num_targets - 1);
  std::uniform_int_distribution<int> local_search_grad_vs_sampling_distribution(0, 1);

  std::uniform_real_distribution<double> local_search_sampling_theta_distribution(0, 2*M_PI);

  std::vector<int> max_newton_iter_for_success_repair_per_thread(num_openmp_threads, 0);
  std::vector<int> max_bisection_iter_for_success_transformation_per_thread(num_openmp_threads, 0);

  // Initialize population
  #pragma omp parallel for
  for (int i = 0; i < pop_size; ++i) {
    population1[i] = initial_population.block(num_targets*i, 0, num_targets, gene_size);
    // check_chromosome_feasible(population1[i], tw_per_target, p0, vmax, q_trj_per_target, target_radii);
  }
  cost_vs_time.push_back(std::pair<double, double>(0., Map<VectorXd>(population_costs1.data(), population_costs1.size()).minCoeff()));
  int num_finite_cost = 0;
  for (auto cost : population_costs1) {
    if (std::isfinite(cost)){ 
      ++num_finite_cost;
    }
  }
  std::cout << "Initialized population of " << num_finite_cost << " finite-cost solutions" << std::endl;

  std::vector<MatrixXd> population2 = population1;
  std::vector<double> population_costs2 = population_costs1;

  int gen_idx = 0;

  std::vector<MatrixXd> *population = &population1;
  std::vector<double> *population_costs = &population_costs1;

  std::vector<MatrixXd> *updated_population = &population1;
  std::vector<double> *updated_population_costs = &population_costs1;

  while (true) {
    auto timer_stop = std::chrono::high_resolution_clock::now();
    auto nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
    if (((double)nanos)/1e9 > time_limit) {
      break;
    }

    // std::cout << "Beginning generation " << gen_idx << std::endl;

    if (gen_idx%2) {
      population = &population2;
      population_costs = &population_costs2;

      updated_population = &population1;
      updated_population_costs = &population_costs1;
    } else {
      population = &population1;
      population_costs = &population_costs1;

      updated_population = &population2;
      updated_population_costs = &population_costs2;
    }

    #pragma omp parallel for
    for (int chromosome_idx = 0; chromosome_idx < pop_size; ++chromosome_idx) {
      int parent1_idx = chromosome_idx;
      int parent2_idx = parent_distribution(rngs_per_thread[omp_get_thread_num()]);
      // Crossover
      VectorXb inserted_targets = VectorXb::Zero(num_targets);
      MatrixXd Xnew(num_targets, gene_size);
      int parent1_counter = 0;
      int parent2_counter = 0;

      for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
        if (crossover_distribution(rngs_per_thread[omp_get_thread_num()]) == 0) {
          while (inserted_targets((int)((*population)[parent1_idx](parent1_counter, 0)))) {
            ++parent1_counter;
          }
          if (parent1_counter >= num_targets) {
            throw std::runtime_error("parent 1 out of bounds");
          }
          int target_idx = (*population)[parent1_idx](parent1_counter, 0);
          Xnew.row(seq_idx) = (*population)[parent1_idx].row(parent1_counter);
          inserted_targets(target_idx) = true;
        } else {
          while (inserted_targets((int)((*population)[parent2_idx](parent2_counter, 0)))) {
            ++parent2_counter;
          }
          if (parent2_counter >= num_targets) {
            throw std::runtime_error("parent 2 out of bounds");
          }
          int target_idx = (*population)[parent2_idx](parent2_counter, 0);
          Xnew.row(seq_idx) = (*population)[parent2_idx].row(parent2_counter);
          inserted_targets(target_idx) = true;
        }
      }

      double mutation_sample = mutation_distribution(rngs_per_thread[omp_get_thread_num()]);
      if (mutation_sample < params.mutation_prob) {
        // Mutation
        if (mutation_sample < params.mutation_prob/3) {
          int seq_idx1 = mutation_operator1_distribution(rngs_per_thread[omp_get_thread_num()]);
          int seq_idx2 = mutation_operator1_distribution(rngs_per_thread[omp_get_thread_num()]);
          RowVectorXd tmp = Xnew.row(seq_idx1);
          Xnew.row(seq_idx1) = Xnew.row(seq_idx2);
          Xnew.row(seq_idx2) = tmp;
        } else if (mutation_sample < 2*params.mutation_prob/3) {
          int seq_idx = mutation_operator2_seq_idx_distribution(rngs_per_thread[omp_get_thread_num()]);
          double theta = mutation_operator2_angle_distribution(rngs_per_thread[omp_get_thread_num()]);
          Xnew(seq_idx, 1) = theta;
        } else {
          int seq_idx = mutation_operator3_seq_idx_distribution(rngs_per_thread[omp_get_thread_num()]);
          double raw_sample = mutation_operator3_delta_t_distribution(rngs_per_thread[omp_get_thread_num()]);
          int target_idx = Xnew(seq_idx, 0);
          double delta_t;
          if (no_tw) {
            double max_delta_t = 108; // TODO: if Yulong answers my question, replace this
            delta_t = max_delta_t*raw_sample;
          } else {
            double min_t = tw_per_target(target_idx, 0);
            double max_t = tw_per_target(target_idx, 1);
            delta_t = (max_t - min_t)*raw_sample;
          }
          Xnew(seq_idx, 2) = delta_t;
        }
      }

      // Repair to restore feasibility
      double cost = 0.;
      double final_heading;
      bool repair_failed = repair_chromosome(Xnew, tw_per_target, target_radii, q_trj_per_target, p0, heading0, vmax, cost, dubins, rho, max_newton_iter_for_success_repair_per_thread[omp_get_thread_num()], params, no_tw, final_heading, 0.);
      if (repair_failed || cost >= (*population_costs)[chromosome_idx]) {
        (*updated_population)[chromosome_idx] = (*population)[chromosome_idx];
        (*updated_population_costs)[chromosome_idx] = (*population_costs)[chromosome_idx];
        continue;
      }

      // check_chromosome_feasible((*updated_population)[chromosome_idx], tw_per_target, p0, vmax, q_trj_per_target, target_radii);

      (*updated_population)[chromosome_idx] = Xnew;
      (*updated_population_costs)[chromosome_idx] = cost;

      // Transformation to reduce cost
      bool success;
      if (dubins) {
        success = transform_chromosome(Xnew, tw_per_target, target_radii, q_trj_per_target, p0, heading0, vmax, rho, cost, speed_upper_bounds, max_newton_iter_for_success_repair_per_thread[omp_get_thread_num()], max_bisection_iter_for_success_transformation_per_thread[omp_get_thread_num()], params, no_tw);
      } else {
        if (no_tw) {
          success = transform_chromosome_no_dubins(Xnew, tw_per_target, target_radii, q_trj_per_target, p0, vmax, cost, max_newton_iter_for_success_repair_per_thread[omp_get_thread_num()], max_bisection_iter_for_success_transformation_per_thread[omp_get_thread_num()], params, no_tw);
        } else {
          success = false;
        }
      }
      if (success && cost < (*updated_population_costs)[chromosome_idx]) {
        (*updated_population)[chromosome_idx] = Xnew;
        (*updated_population_costs)[chromosome_idx] = cost;
        // std::cout << "transformation reduced cost" << std::endl;
      }
    }

    /*
    for (int chromosome_idx = 0; chromosome_idx < pop_size; ++chromosome_idx) {
      check_chromosome_feasible((*updated_population)[chromosome_idx], tw_per_target);
    }
    */
    
    ++gen_idx;
    // -1 means don't run local search
    if (params.Tlp != -1 && gen_idx%params.Tlp == 0) {
      std::vector<size_t> sort_idx = sort_indexes(*updated_population_costs);
      for (int j = 0; j < gen_idx/params.Tlp; ++j) {
        // Get individual from top 50% and run local search
        int chromosome_idx = sort_idx[local_search_elite_distribution(rngs_per_thread[omp_get_thread_num()])];
        int gene_idx = local_search_gene_idx_distribution(rngs_per_thread[omp_get_thread_num()]);
        double theta = (*updated_population)[chromosome_idx](gene_idx, 1);
        // Only run gradient-based local search on feasible solutions (only relevant if I do Dubins close-enough)
        if (std::isfinite((*updated_population_costs)[chromosome_idx]) && local_search_grad_vs_sampling_distribution(rngs_per_thread[omp_get_thread_num()]) == 0) {
          // Gradient-based local search
          bool improvement = true;
          while (improvement) {
            int max_backtrack_fd = 3;
            double delta_theta = 0.01; // To compute approximate gradient
            double new_cost = std::numeric_limits<double>::infinity();
            for (int backtrack_iter = 0; backtrack_iter < max_backtrack_fd; ++backtrack_iter) {
              double new_theta = theta + delta_theta;

              MatrixXd local_modification = (*updated_population)[chromosome_idx];
              local_modification(gene_idx, 1) = new_theta;

              double final_heading;
              bool repair_failed = repair_chromosome(local_modification, tw_per_target, target_radii, q_trj_per_target, p0, heading0, vmax, new_cost, dubins, rho, max_newton_iter_for_success_repair_per_thread[omp_get_thread_num()], params, no_tw, final_heading, 0.);
              if (!repair_failed) {
                double tmp_cost;
                bool transformation_succeeded;
                if (dubins) {
                  transformation_succeeded = transform_chromosome(local_modification, tw_per_target, target_radii, q_trj_per_target, p0, heading0, vmax, rho, tmp_cost, speed_upper_bounds, max_newton_iter_for_success_repair_per_thread[omp_get_thread_num()], max_bisection_iter_for_success_transformation_per_thread[omp_get_thread_num()], params, no_tw);
                } else {
                  transformation_succeeded = false; // transform_chromosome_no_dubins(local_modification, tw_per_target, target_radii, q_trj_per_target, p0, vmax, tmp_cost, max_newton_iter_for_success_repair_per_thread[omp_get_thread_num()], max_bisection_iter_for_success_transformation_per_thread[omp_get_thread_num()], params);
                }
                if (transformation_succeeded && tmp_cost < new_cost) {
                  new_cost = tmp_cost;
                }

                break;
              }
              delta_theta *= 0.1;
            }
            if (std::isinf(new_cost)) {
              break;
            }
            double gradient = (new_cost - (*updated_population_costs)[chromosome_idx])/delta_theta;

            int max_backtrack_gd = 3;
            double gd_step_size = params.local_search_gd_step_size;
            for (int backtrack_iter = 0; backtrack_iter < max_backtrack_gd; ++backtrack_iter) {
              double new_theta = theta - gd_step_size*gradient;

              MatrixXd local_modification = (*updated_population)[chromosome_idx];
              local_modification(gene_idx, 1) = new_theta;

              double final_heading;
              bool repair_failed = repair_chromosome(local_modification, tw_per_target, target_radii, q_trj_per_target, p0, heading0, vmax, new_cost, dubins, rho, max_newton_iter_for_success_repair_per_thread[omp_get_thread_num()], params, no_tw, final_heading, 0.);
              if (!repair_failed && new_cost < (*updated_population_costs)[chromosome_idx]) {
                double tmp_cost;
                bool transformation_succeeded;
                MatrixXd tmp_local_modification = local_modification;
                if (dubins) {
                  transformation_succeeded = transform_chromosome(tmp_local_modification, tw_per_target, target_radii, q_trj_per_target, p0, heading0, vmax, rho, tmp_cost, speed_upper_bounds, max_newton_iter_for_success_repair_per_thread[omp_get_thread_num()], max_bisection_iter_for_success_transformation_per_thread[omp_get_thread_num()], params, no_tw);
                } else {
                  if (no_tw) {
                    transformation_succeeded = transform_chromosome_no_dubins(tmp_local_modification, tw_per_target, target_radii, q_trj_per_target, p0, vmax, tmp_cost, max_newton_iter_for_success_repair_per_thread[omp_get_thread_num()], max_bisection_iter_for_success_transformation_per_thread[omp_get_thread_num()], params, no_tw);
                  } else {
                    transformation_succeeded = false;
                  }
                }
                if (transformation_succeeded && tmp_cost < new_cost) {
                  (*updated_population)[chromosome_idx] = tmp_local_modification;
                  (*updated_population_costs)[chromosome_idx] = tmp_cost;
                } else {
                  (*updated_population)[chromosome_idx] = local_modification;
                  (*updated_population_costs)[chromosome_idx] = new_cost;
                }
                break;
              }
              gd_step_size *= 0.1;
            }

            improvement = new_cost < (*updated_population_costs)[chromosome_idx] - 1e-4;
          }
        } else {
          double theta = (*updated_population)[chromosome_idx](gene_idx, 1);

          // Do sampling-based local search
          std::vector<MatrixXd> local_modifications(params.local_search_num_samples);
          std::vector<double> local_modification_costs(params.local_search_num_samples);
          #pragma omp parallel for
          for (int sample_idx = 0; sample_idx < params.local_search_num_samples; ++sample_idx) {
            // double new_theta = local_search_sampling_theta_distribution(rngs_per_thread[omp_get_thread_num()]);
            double new_theta = (2*M_PI*sample_idx)/params.local_search_num_samples;

            double new_cost;
            local_modifications[sample_idx] = (*updated_population)[chromosome_idx];
            local_modifications[sample_idx](gene_idx, 1) = new_theta;

            double final_heading;
            bool repair_failed = repair_chromosome(local_modifications[sample_idx], tw_per_target, target_radii, q_trj_per_target, p0, heading0, vmax, new_cost, dubins, rho, max_newton_iter_for_success_repair_per_thread[omp_get_thread_num()], params, no_tw, final_heading, 0.);

            if (!(repair_failed || new_cost >= (*updated_population_costs)[chromosome_idx])) {
              if (dubins) {
                double tmp_cost;
                MatrixXd tmp_local_modification = local_modifications[sample_idx];
                bool transformation_succeeded = transform_chromosome(tmp_local_modification, tw_per_target, target_radii, q_trj_per_target, p0, heading0, vmax, rho, tmp_cost, speed_upper_bounds, max_newton_iter_for_success_repair_per_thread[omp_get_thread_num()], max_bisection_iter_for_success_transformation_per_thread[omp_get_thread_num()], params, no_tw);
                if (transformation_succeeded && tmp_cost < new_cost) {
                  local_modification_costs[sample_idx] = tmp_cost;
                } else {
                  local_modification_costs[sample_idx] = new_cost;
                  local_modifications[sample_idx] = tmp_local_modification;
                }
              } else {
                local_modification_costs[sample_idx] = new_cost;
              }
              // check_chromosome_feasible(local_modifications[sample_idx], tw_per_target);
            } else {
              local_modification_costs[sample_idx] = std::numeric_limits<double>::infinity();
            }
          } // Sampling-based local search iterations

          auto it = std::min_element(local_modification_costs.begin(), local_modification_costs.end());
          int min_idx = it - local_modification_costs.begin();
          if (local_modification_costs[min_idx] < (*updated_population_costs)[chromosome_idx]) {
            (*updated_population)[chromosome_idx] = local_modifications[min_idx];
            (*updated_population_costs)[chromosome_idx] = local_modification_costs[min_idx];
            // check_chromosome_feasible((*updated_population)[chromosome_idx], tw_per_target, p0, vmax, q_trj_per_target, target_radii);
          }
        } // Sampling-based local search instead of gradient
      } // Pick a chromosome for local search
    } // Check if local search condition has been met

    auto timer_start2 = std::chrono::high_resolution_clock::now();
    auto it2 = std::min_element((*updated_population_costs).begin(), (*updated_population_costs).end());
    auto timer_stop2 = std::chrono::high_resolution_clock::now();
    auto nanos2 = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop2 - timer_start2).count();
    min_cost_record_time += ((double)nanos2)/1e9;

    timer_stop = std::chrono::high_resolution_clock::now();
    nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
    if (cost_vs_time.size() && *it2 > cost_vs_time.back().second) {
      throw std::runtime_error("Cost increased after memetic alg iteration");
    }
    if (cost_vs_time.size() == 0 || *it2 < cost_vs_time.back().second) {
      cost_vs_time.push_back(std::pair<double, double>(((double)nanos)/1e9, *it2));
    }
  } // Overall loop

  num_feas_final(0) = 0;
  for (auto cost : (*updated_population_costs)) {
    if (std::isfinite(cost)){ 
      ++num_feas_final(0);
    }
  }
  std::cout << "Final population contains " << num_feas_final(0) << " feas solns" << std::endl;

  auto it = std::min_element((*updated_population_costs).begin(), (*updated_population_costs).end());
  int min_idx = it - (*updated_population_costs).begin();
  const Ref<const MatrixXd> &Xbest = (*updated_population)[min_idx];

  if (std::isinf(*it)) {
  } else if (dubins) {
    double t = 0;
    Vector2d pos = p0;
    Vector2d next_pos;
    double heading = heading0;
    for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
      int target_idx = Xbest(seq_idx, 0);
      double theta = Xbest(seq_idx, 1);
      double delta_t = Xbest(seq_idx, 2);
      t += delta_t;
      selected_pts_per_target(target_idx, 0) = t;
      next_pos = q_trj_per_target[target_idx](t) + target_radii[target_idx]*Vector2d(cos(theta), sin(theta));

      RowMatrixXd turns = elongated_dubins_path_one_sided(pos(0), pos(1), heading, next_pos(0), next_pos(1), vmax*delta_t, rho, 1e-2);
      if (std::isinf(turns(0, 0))) {
        throw std::runtime_error("elongated path generation failed");
      }

      pos = next_pos;
      for (int row = 0; row < turns.rows(); ++row) {
        if (turns(row, 0) != 0) {
          heading += turns(row, 1)/turns(row, 0);
        }
      }

      turns_chain.push_back(turns);

      selected_pts_per_target(target_idx, 1) = pos(0);
      selected_pts_per_target(target_idx, 2) = pos(1);
      selected_pts_per_target(target_idx, 3) = heading;

      if (!no_tw && (t < tw_per_target(target_idx, 0) - 1e-4 || t > tw_per_target(target_idx, 1) + 1e-4)) {
        throw std::runtime_error("t out of window");
      }
    }
  } else {
    double t = 0;
    Vector2d pos;
    for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
      int target_idx = Xbest(seq_idx, 0);
      double theta = Xbest(seq_idx, 1);
      double delta_t = Xbest(seq_idx, 2);
      t += delta_t;
      selected_pts_per_target(target_idx, 0) = t;
      pos = q_trj_per_target[target_idx](t) + target_radii[target_idx]*Vector2d(cos(theta), sin(theta));
      selected_pts_per_target(target_idx, 1) = pos(0);
      selected_pts_per_target(target_idx, 2) = pos(1);

      if (!no_tw && (t < tw_per_target(target_idx, 0) - 1e-4 || t > tw_per_target(target_idx, 1) + 1e-4)) {
        throw std::runtime_error("t out of window");
      }
    }
  }

  RowMatrixXd cost_vs_time_mat(cost_vs_time.size(), 2);
  for (int i = 0; i < cost_vs_time.size(); ++i) {
    cost_vs_time_mat(i, 0) = cost_vs_time[i].first;
    cost_vs_time_mat(i, 1) = cost_vs_time[i].second;
  }

  std::cout << "Max newton iterations for successful repair: " << *(std::min_element(max_newton_iter_for_success_repair_per_thread.begin(), max_newton_iter_for_success_repair_per_thread.end())) << std::endl;
  std::cout << "Max bisection iterations for successful transformation: " << *(std::min_element(max_bisection_iter_for_success_transformation_per_thread.begin(), max_bisection_iter_for_success_transformation_per_thread.end())) << std::endl;

  std::cout << "Spent " << min_cost_record_time << " s tracking what the min cost was after each iteration (just making sure this is not too large)" << std::endl;
  return cost_vs_time_mat;
}

void get_selected_pts(Ref<RowMatrixXd> selected_pts_per_target, const std::vector<ExtendedCppSpline> q_trj_per_target, const Ref<const Vector2d> &p0, double heading0, double vmax, double rho, const Ref<const MatrixXd> &X, const Ref<const RowMatrixXd> &tw_per_target, const Ref<const VectorXd> &target_radii, bool no_tw) {
  bool dubins = rho != 0.;
  int num_targets = tw_per_target.rows();
  if (dubins) {
    double t = 0;
    Vector2d pos = p0;
    Vector2d next_pos;
    double heading = heading0;
    for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
      int target_idx = X(seq_idx, 0);
      double theta = X(seq_idx, 1);
      double delta_t = X(seq_idx, 2);
      t += delta_t;
      selected_pts_per_target(target_idx, 0) = t;
      next_pos = q_trj_per_target[target_idx](t) + target_radii[target_idx]*Vector2d(cos(theta), sin(theta));

      RowMatrixXd turns = elongated_dubins_path_one_sided(pos(0), pos(1), heading, next_pos(0), next_pos(1), vmax*delta_t, rho, 1e-2);
      if (std::isinf(turns(0, 0))) {
        throw std::runtime_error("elongated path generation failed");
      }

      pos = next_pos;
      for (int row = 0; row < turns.rows(); ++row) {
        if (turns(row, 0) != 0) {
          heading += turns(row, 1)/turns(row, 0);
        }
      }

      selected_pts_per_target(target_idx, 1) = pos(0);
      selected_pts_per_target(target_idx, 2) = pos(1);
      selected_pts_per_target(target_idx, 3) = heading;

      if (!no_tw && (t < tw_per_target(target_idx, 0) - 1e-4 || t > tw_per_target(target_idx, 1) + 1e-4)) {
        throw std::runtime_error("t out of window");
      }
    }
  } else {
    double t = 0;
    Vector2d pos = p0;
    std::cout << "Checking feasibility in get selected pts" << std::endl;
    for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
      int target_idx = X(seq_idx, 0);
      double theta = X(seq_idx, 1);
      double delta_t = X(seq_idx, 2);
      t += delta_t;
      selected_pts_per_target(target_idx, 0) = t;
      Vector2d next_pos = q_trj_per_target[target_idx](t) + target_radii[target_idx]*Vector2d(cos(theta), sin(theta));
      if ((next_pos - pos).norm() > vmax*delta_t) {
        throw std::runtime_error("Speed constraint violated in get selected pts");
      }
      pos = next_pos;
      selected_pts_per_target(target_idx, 1) = pos(0);
      selected_pts_per_target(target_idx, 2) = pos(1);

      if (!no_tw && (t < tw_per_target(target_idx, 0) - 1e-4 || t > tw_per_target(target_idx, 1) + 1e-4)) {
        throw std::runtime_error("t out of window");
      }
    }
  }
}
