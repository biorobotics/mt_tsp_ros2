#pragma once
#include "mt_tsp_ros2/gcs_node.h"
#include <memory>

struct TrajoptThroughConvexSets {
  TrajoptThroughConvexSets(int dim_q, const Ref<const VectorXd> &set_lb, const Ref<const VectorXd> &set_ub, const Ref<const VectorXd> &start_pos_and_time, double vmax_agent);

  void add_convex_set(const GCSNode &convex_set);
  double solve_trajopt(std::vector<VectorXd> &pos_seq, std::vector<double> &time_seq, bool min_time, bool verbose);

  int dim_q;
  int set_dim;
  VectorXd set_lb;
  VectorXd set_ub;
  VectorXd start_pos_and_time;

  int z_idx;
  int l_idx;

  int vars_per_step;
  int steps;

  VectorXi zero_cone_rows_per_step;
  VectorXi zero_cone_cols_per_step;
  VectorXd zero_cone_vals_per_step;
  VectorXd zero_cone_b_per_step;
  int num_zero_cone_per_step;

  int num_dist_cons_rows;

  VectorXi linear_cone_rows_per_step;
  VectorXi linear_cone_cols_per_step;
  VectorXd linear_cone_vals_per_step;
  VectorXd linear_cone_b_per_step;
  int num_linear_cone_per_step;

  VectorXi soc_rows_per_step;
  VectorXi soc_cols_per_step;
  VectorXd soc_vals_per_step;
  VectorXd soc_b_per_step;
  int num_soc_per_step;

  std::vector<VectorXi> zero_cone_rows;
  std::vector<VectorXi> zero_cone_cols;
  std::vector<VectorXd> zero_cone_vals;
  std::vector<VectorXd> zero_cone_b;

  std::vector<VectorXi> linear_cone_rows;
  std::vector<VectorXi> linear_cone_cols;
  std::vector<VectorXd> linear_cone_vals;
  std::vector<VectorXd> linear_cone_b;

  std::vector<VectorXi> soc_rows;
  std::vector<VectorXi> soc_cols;
  std::vector<VectorXd> soc_vals;
  std::vector<VectorXd> soc_b;
  std::vector<int> soc_sizes;

  int num_zero_cone;
  int num_linear_cone;
  int num_soc;
};
