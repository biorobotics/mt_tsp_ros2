#pragma once
#include <chrono>
#include <omp.h>
#include <Eigen/Dense>
#include <random>
#include <set>
#include "mt_tsp_ros2/cpp_spline.h"
#include "mt_tsp_ros2/extended_cpp_spline.h"
#include "mt_tsp_ros2/elongate_one_sided_dubins_path.h"
#include <iomanip>

using namespace Eigen;

typedef Matrix<double, Dynamic, Dynamic, RowMajor> RowMatrixXd;

typedef Matrix<bool, Dynamic, 1> VectorXb;

typedef Matrix<double, 1, 1> Vector1d;

typedef const Ref<const Matrix<long, Dynamic, 1>> &VectorXlRef_const;
typedef const Ref<const Matrix<double, Dynamic, 1>> &VectorXdRef_const;
typedef const Ref<const RowMatrixXd> &RowMatrixXdRef_const;

class DubinsTrjThroughSeqOfTargets {
  public:
    DubinsTrjThroughSeqOfTargets(const Ref<const RowMatrixXd> &tw_per_target, const std::vector<ExtendedCppSpline> &q_trj_per_target, const Ref<const Vector2d> &p0, double heading0, double vmax, double rho, bool no_tw, double t0, bool min_latency, bool min_time, double newton_step_size) : tw_per_target(tw_per_target), q_trj_per_target(q_trj_per_target), p0(p0), heading0(heading0), vmax(vmax), rho(rho), no_tw(no_tw), t0(t0), min_latency(min_latency), min_time(min_time), newton_step_size(newton_step_size) {
    }

    bool optimize_trj(Ref<Vector1d> cost, Ref<RowMatrixXd> selected_pts_per_target, VectorXlRef_const target_seq) {
      cost(0) = 0.;
      int num_targets = tw_per_target.rows();

      if (selected_pts_per_target.rows() != num_targets || selected_pts_per_target.cols() != 4) {
        throw std::runtime_error("selected_pts_per_target has incorrect size");
      }

      double t = t0;
      Vector2d pos = p0;
      double heading = heading0;
      Vector2d next_pos;
      for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
        int target_idx = target_seq(seq_idx);

        double next_t;
        if (no_tw) {
          throw std::runtime_error("Did not implement no-time-window case");
        } else {
          // Check start of time window
          next_t = tw_per_target(target_idx, 0);
          next_pos = q_trj_per_target[target_idx](next_t);
          double delta_t = next_t - t;
          RowMatrixXd turns = elongated_dubins_path_one_sided(pos(0), pos(1), heading, next_pos(0), next_pos(1), vmax*delta_t, rho, 1e-4);
          if (std::isfinite(turns(0, 0))) {
            selected_pts_per_target(target_idx, 0) = next_t;
            selected_pts_per_target(target_idx, 1) = next_pos(0);
            selected_pts_per_target(target_idx, 2) = next_pos(1);
            double next_heading = heading;
            for (int row = 0; row < turns.rows(); ++row) {
              if (turns(row, 0) != 0) {
                next_heading += turns(row, 1)/turns(row, 0);
              }
            }
            selected_pts_per_target(target_idx, 3) = next_heading;
            t = next_t;
            pos = next_pos;
            heading = next_heading;

            if (no_tw) {
              if (min_latency) {
                cost(0) += next_t;
              } else if (seq_idx == num_targets - 1) {
                if (min_time) {
                  cost(0) += next_t;
                } else {
                  cost(0) += vmax*next_t;
                }
              }
            } else {
              if (min_latency) {
                cost(0) += next_t - tw_per_target(target_idx, 0);
              } else if (seq_idx == num_targets - 1) {
                if (min_time) {
                  cost(0) += next_t;
                } else {
                  cost(0) += vmax*next_t;
                }
              }
            }

            continue;
          }

          // Check end of time window
          next_t = tw_per_target(target_idx, 1);
          next_pos = q_trj_per_target[target_idx](next_t);
          delta_t = next_t - t;
          turns = elongated_dubins_path_one_sided(pos(0), pos(1), heading, next_pos(0), next_pos(1), vmax*delta_t, rho, 1e-4);
          if (std::isinf(turns(0, 0))) {
            return false;
          }
        }

        // Run Newton to restore feasibility if needed
        int max_newton_iter = 2000;
        double delta_delta_t_finite_diff = 1e-4;
        bool newton_succeeded = false;
        // std::cout << "starting newton" << std::endl;
        double feas_next_t = std::numeric_limits<double>::infinity();
        Vector2d feas_next_pos = std::numeric_limits<double>::infinity()*Vector2d::Ones();
        double feas_next_heading = std::numeric_limits<double>::infinity();
        double delta_t = next_t - t;
        for (int newton_iter = 0; newton_iter < max_newton_iter; ++newton_iter) {
          next_pos = q_trj_per_target[target_idx](next_t);

          RowMatrixXd turns = elongated_dubins_path_one_sided(pos(0), pos(1), heading, next_pos(0), next_pos(1), vmax*delta_t, rho, 1e-4);
          if (std::isfinite(turns(0, 0)) && next_t < feas_next_t) {
            newton_succeeded = true;
            feas_next_t = next_t;
            feas_next_pos = next_pos;
            feas_next_heading = heading;
            for (int row = 0; row < turns.rows(); ++row) {
              if (turns(row, 0) != 0) {
                feas_next_heading += turns(row, 1)/turns(row, 0);
              }
            }
          }

          RowMatrixXd shortest_path_turns = turns_for_one_sided_dubins_path(pos(0), pos(1), heading, next_pos(0), next_pos(1), rho);
          double shortest_path_dist = shortest_path_turns.col(1).sum();

          double c = shortest_path_dist - vmax*delta_t;

          if (std::abs(c) < 1e-4) {
            if (std::isinf(turns(0, 0))) {
              throw std::runtime_error("Path elongation should have just returned the shortest dubins path but the elongation actually failed");
            }
            break;
          }

          // Finite-diff
          double delta_t_plus = delta_t + delta_delta_t_finite_diff;

          double next_t_plus = t + delta_t_plus;

          Vector2d next_pos_plus = q_trj_per_target[target_idx](next_t_plus);

          RowMatrixXd shortest_path_turns_plus = turns_for_one_sided_dubins_path(pos(0), pos(1), heading, next_pos_plus(0), next_pos_plus(1), rho);
          double shortest_path_dist_plus = shortest_path_turns_plus.col(1).sum();
          double c_plus = shortest_path_dist_plus - vmax*delta_t_plus;

          double deriv = (c_plus - c)/delta_delta_t_finite_diff;

          delta_t -= newton_step_size*c/deriv;

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
        }

        if (!newton_succeeded) {
          std::cout << "newton failed at seq_idx" << seq_idx << std::endl;
          return false;
        }

        next_t = feas_next_t;
        next_pos = feas_next_pos;
        heading = feas_next_heading;
        selected_pts_per_target(target_idx, 0) = next_t;
        selected_pts_per_target(target_idx, 1) = next_pos(0);
        selected_pts_per_target(target_idx, 2) = next_pos(1);
        selected_pts_per_target(target_idx, 3) = heading;

        if (no_tw) {
          if (min_latency) {
            cost(0) += next_t;
          } else if (seq_idx == num_targets - 1) {
            if (min_time) {
              cost(0) += next_t;
            } else {
              cost(0) += vmax*next_t;
            }
          }
        } else {
          if (min_latency) {
            cost(0) += next_t - tw_per_target(target_idx, 0);
          } else if (seq_idx == num_targets - 1) {
            if (min_time) {
              cost(0) += next_t;
            } else {
              cost(0) += vmax*next_t;
            }
          }
        }

        t = next_t;
        pos = next_pos;
      }
      return true;
    }

  private:
    RowMatrixXd tw_per_target;
    std::vector<ExtendedCppSpline> q_trj_per_target;
    Vector2d p0;
    double heading0;
    double vmax;
    double rho;
    bool no_tw;
    double t0;
    bool min_latency;
    bool min_time;
    double newton_step_size;
};
