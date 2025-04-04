#pragma once
#include <memory>
#include "mt_tsp_ros2/gcs_node.h"
#include "mt_tsp_ros2/trajopt_through_convex_sets.h"

struct GCSPath {
  GCSPath(const std::vector<int> &node_seq, 
          int seq_list_idx,
          double f_val,
          double cost_lb,
          double g_val,
          double h_val,
          int tie_break_val,
          const TrajoptThroughConvexSets &trajopt) : node_seq(node_seq),
                                                     seq_list_idx(seq_list_idx),
                                                     f_val(f_val),
                                                     cost_lb(cost_lb),
                                                     g_val(g_val),
                                                     h_val(h_val),
                                                     tie_break_val(tie_break_val),
                                                     trajopt(trajopt) {
  }
  GCSPath(const TrajoptThroughConvexSets &trajopt) : trajopt(trajopt) {
  }
  std::vector<int> node_seq;
  int seq_list_idx;
  double f_val;
  double cost_lb;
  double g_val;
  double h_val;
  int tie_break_val;
  TrajoptThroughConvexSets trajopt;
};

class MTIxGStar {
  public:
    MTIxGStar(const std::vector<GCSNode> &gcs_nodes, const Ref<const Matrix<bool, Dynamic, Dynamic, RowMajor>> &gcs_adj_mat, int dim_q, const Ref<const VectorXd> &set_lb, const Ref<const VectorXd> &set_ub, double vmax_agent);

    MatrixXd solve(int start_idx, int goal_idx, double w, VectorXd start_pos_and_time, double max_final_time, const Ref<const VectorXl> &ignore_set_vec, bool verbose, double time_limit, double tw_start);

  private:
    int dim_q;
    VectorXd set_lb;
    VectorXd set_ub;
    double vmax_agent;
    std::vector<GCSNode> gcs_nodes;
    Matrix<bool, Dynamic, Dynamic> gcs_adj_mat;

    std::shared_ptr<GCSPath> extend_path(std::shared_ptr<GCSPath> path, int next_node_idx, int goal_node_idx, double w, double max_final_time, std::vector<MatrixXd> &pos_seq_list, std::vector<VectorXd> &time_seq_list);
};
