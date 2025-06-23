#pragma once
#include <chrono>
#include <omp.h>
#include <Eigen/Dense>
#include <random>
#include <set>
#include "mt_tsp_ros2/cpp_spline.h"
#include "mt_tsp_ros2/extended_cpp_spline.h"
#include <iomanip>

using namespace Eigen;

typedef Matrix<double, Dynamic, Dynamic, RowMajor> RowMatrixXd;

typedef Matrix<bool, Dynamic, 1> VectorXb;

typedef Matrix<double, 1, 1> Vector1d;

typedef const Ref<const Matrix<long, Dynamic, 1>> &VectorXlRef_const;
typedef const Ref<const Matrix<double, Dynamic, 1>> &VectorXdRef_const;
typedef const Ref<const RowMatrixXd> &RowMatrixXdRef_const;

const double root_finding_tol = 1e-2;

class EFAT {
  public:
    EFAT(const Ref<const RowMatrixXd> &tw_per_target, const std::vector<ExtendedCppSpline> &q_trj_per_target, const Ref<const Vector2d> &p0, double vmax, bool no_tw, double t0) : tw_per_target(tw_per_target), q_trj_per_target(q_trj_per_target), p0(p0), vmax(vmax), no_tw(no_tw), t0(t0) {
    }

    bool efat_chain(Ref<Vector1d> cost, Ref<RowMatrixXd> selected_pts_per_target, VectorXlRef_const target_seq, bool feasible_times) {
      int num_targets = tw_per_target.rows();

      if (selected_pts_per_target.rows() != num_targets || selected_pts_per_target.cols() != 3) {
        throw std::runtime_error("selected_pts_per_target has incorrect size");
      }

      double t = t0;
      Vector2d pos = p0;
      Vector2d next_rel_pos;
      Vector2d next_pos;
      cost(0) = 0;
      for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
        int target_idx = target_seq(seq_idx);

        next_pos = selected_pts_per_target.block<1, 2>(target_idx, 1).transpose();
        next_rel_pos = next_pos - q_trj_per_target[target_idx](selected_pts_per_target(target_idx, 0));

        // Newton
        // Check if travel is feasible to next_rel_pos at start of time window
        double next_t = tw_per_target(target_idx, 0);
        Vector2d next_pos = q_trj_per_target[target_idx](next_t) + next_rel_pos;
        double delta_t = next_t - t;
        double dist = (next_pos - pos).norm();
        if (dist <= vmax*delta_t) {
          // Travel is feasible to start of time window
          selected_pts_per_target(target_idx, 0) = next_t;
          selected_pts_per_target(target_idx, 1) = next_pos(0);
          selected_pts_per_target(target_idx, 2) = next_pos(1);
          cost(0) += next_t - tw_per_target(target_idx, 0);
          t = next_t;
          pos = next_pos;
          continue;
        }

        next_t = feasible_times ? selected_pts_per_target(target_idx, 0) : tw_per_target(target_idx, 1);
        delta_t = next_t - t;
        next_pos = q_trj_per_target[target_idx](next_t) + next_rel_pos;
        dist = (next_pos - pos).norm();

        // Check if travel is feasible to next_rel_pos at end of time window
        if (!feasible_times && dist > vmax*delta_t) {
          // Travel is infeasible to end of time window
          return false;
        }

        int max_newton_iter = 10;
        double feas_next_t = next_t;
        Vector2d feas_next_pos = next_pos;
        for (int newton_iter = 0; newton_iter < max_newton_iter; ++newton_iter) {
          // Find root of dist - vmax*delta_t
          double resid = dist - vmax*delta_t + 1e-4; // The 1e-4 is so we actually get to a feasible solution
          // delta_vs_iterations.push_back(resid);
          double deriv = 1/(2*dist)*(next_pos - pos).dot(q_trj_per_target[target_idx].derivatives(next_t)) - vmax;
          delta_t -= resid/deriv;
          next_t = t + delta_t;
          next_pos = q_trj_per_target[target_idx](next_t) + next_rel_pos;
          dist = (next_pos - pos).norm();
          if (dist <= vmax*delta_t) {
            if (next_t < feas_next_t) {
              feas_next_t = next_t;
              feas_next_pos = next_pos;
            }

            if (std::abs(resid) < root_finding_tol) {
              break;
            }
          }
        }

        next_t = feas_next_t;
        next_pos = feas_next_pos;

        selected_pts_per_target(target_idx, 0) = next_t;
        selected_pts_per_target(target_idx, 1) = next_pos(0);
        selected_pts_per_target(target_idx, 2) = next_pos(1);

        // Update cost, time, and position
        if (no_tw) {
          cost(0) += next_t;
        } else {
          cost(0) += next_t - tw_per_target(target_idx, 0);
        }

        t = next_t;
        pos = next_pos;

        // Bisection version
        /*
        double t_low;

        double t_high = feasible_times ? selected_pts_per_target(target_idx, 0) : tw_per_target(target_idx, 1);

        if (no_tw) {
          if (!feasible_times) {
            throw std::runtime_error("Did not implement efat chain for no time window case without given upper bounds");
          }
          t_low = t;
        } else {
          if (!feasible_times) {
            // Check if travel is feasible to next_rel_pos at end of time window
            double next_t = tw_per_target(target_idx, 1);
            Vector2d next_pos = q_trj_per_target[target_idx](next_t) + next_rel_pos;
            double delta_t = next_t - t;
            double dist = (next_pos - pos).norm();
            if (dist > vmax*delta_t) {
              // Travel is infeasible to end of time window
              return false;
            }
          }

          // Check if travel is feasible to next_rel_pos at start of time window
          double next_t = tw_per_target(target_idx, 0);
          Vector2d next_pos = q_trj_per_target[target_idx](next_t) + next_rel_pos;
          double delta_t = next_t - t;
          double dist = (next_pos - pos).norm();
          if (dist <= vmax*delta_t) {
            // Travel is feasible to start of time window
            selected_pts_per_target(target_idx, 0) = next_t;
            selected_pts_per_target(target_idx, 1) = next_pos(0);
            selected_pts_per_target(target_idx, 2) = next_pos(1);
            cost(0) += next_t - tw_per_target(target_idx, 0);
            t = next_t;
            pos = next_pos;
            continue;
          }

          t_low = std::max(t, tw_per_target(target_idx, 0));
        }

        // Run bisection to find earliest time such that interception is feasible
        int num_bisection_iter = 10;
        for (int bisection_iter = 0; bisection_iter < num_bisection_iter; ++bisection_iter) {
          double t_mid = 0.5*(t_low + t_high);
          double delta_t = t_mid - t;
          Vector2d pos_mid = q_trj_per_target[target_idx](t_mid) + next_rel_pos;
          double dist = (pos_mid - pos).norm();
          if (dist > vmax*delta_t) {
            // Travel is infeasible
            t_low = t_mid;
          } else {
            // Travel is feasible
            t_high = t_mid;
            next_pos = pos_mid;
          }
        }

        selected_pts_per_target(target_idx, 0) = t_high;
        selected_pts_per_target(target_idx, 1) = next_pos(0);
        selected_pts_per_target(target_idx, 2) = next_pos(1);

        // Update cost, time, and position
        if (no_tw) {
          cost(0) += t_high;
        } else {
          cost(0) += t_high - tw_per_target(target_idx, 0);
        }

        t = t_high;
        pos = next_pos;
        */
      }
      return true;
    }

  private:
    RowMatrixXd tw_per_target;
    std::vector<ExtendedCppSpline> q_trj_per_target;
    Vector2d p0;
    double vmax;
    bool no_tw;
    double t0;
};
