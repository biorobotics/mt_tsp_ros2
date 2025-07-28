#pragma once
#include <chrono>
#include <omp.h>
#include <Eigen/Dense>
#include <random>
#include <set>
#include "mt_tsp_ros2/cpp_spline.h"
#include "mt_tsp_ros2/extended_cpp_spline.h"
#include "mt_tsp_ros2/memetic_dubins_ce_mt_tsp/memetic_alg_vs_dubins.h"
#include <iomanip>

using namespace Eigen;

typedef Matrix<double, Dynamic, Dynamic, RowMajor> RowMatrixXd;

typedef Matrix<bool, Dynamic, 1> VectorXb;

typedef Matrix<double, 1, 1> Vector1d;

typedef const Ref<const Matrix<long, Dynamic, 1>> &VectorXlRef_const;
typedef const Ref<const Matrix<double, Dynamic, 1>> &VectorXdRef_const;
typedef const Ref<const RowMatrixXd> &RowMatrixXdRef_const;

class VSDubinsTrjThroughSeqOfTargets {
  public:
    VSDubinsTrjThroughSeqOfTargets(const Ref<const RowMatrixXd> &tw_per_target, const std::vector<ExtendedCppSpline> &q_trj_per_target, const Ref<const Vector2d> &p0, double heading0, const Ref<const VectorXd> &speed_options, double wmax, double t0, bool min_latency, bool min_time, bool tree_search) : tw_per_target(tw_per_target), q_trj_per_target(q_trj_per_target), p0(p0), heading0(heading0), speed_options(speed_options), wmax(wmax), t0(t0), params(0., 0., 0, -1, 1., min_latency, min_time), tree_search(tree_search), num_feas_chromosomes_where_transformation_improved_cost(1), num_feas_chromosomes_generated(1), num_feas_chromosomes_where_transformation_was_feasible(1) {
    }

    bool optimize_trj(Ref<Vector1d> cost, Ref<RowMatrixXd> selected_pts_per_target, VectorXlRef_const target_seq, double time_limit) {
      int num_targets = tw_per_target.rows();
      MatrixXd chromosome = MatrixXd::Zero(num_targets, gene_size);
      chromosome.col(0) = target_seq.cast<double>();
      double tmp_cost;
      VectorXl num_newton_successes_when_tw_end_check_failed(1);
      VectorXl num_newton_solves(1);
      VectorXl num_newton_successes(1);
      repair_chromosome(chromosome, tw_per_target, q_trj_per_target, p0, heading0, speed_options, tmp_cost, wmax, params, t0, time_limit, selected_pts_per_target, true, tree_search, num_feas_chromosomes_where_transformation_improved_cost, num_feas_chromosomes_generated, num_feas_chromosomes_where_transformation_was_feasible, 0, num_newton_successes_when_tw_end_check_failed, num_newton_solves, num_newton_successes);
      cost(0) = tmp_cost;
      return std::isfinite(tmp_cost);
    }

  private:
    RowMatrixXd tw_per_target;
    std::vector<ExtendedCppSpline> q_trj_per_target;
    Vector2d p0;
    double heading0;
    VectorXd speed_options;
    double wmax;
    double t0;
    MemeticAlgParams params;
    bool tree_search;
    VectorXl num_feas_chromosomes_where_transformation_improved_cost;
    VectorXl num_feas_chromosomes_generated;
    VectorXl num_feas_chromosomes_where_transformation_was_feasible;
};
