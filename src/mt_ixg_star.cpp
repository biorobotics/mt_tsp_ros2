#include "mt_tsp_ros2/mt_ixg_star.h"
#include <chrono>
#include <limits>
#include <unordered_set>
#include <queue>

using namespace std::chrono;

double nan_num() {
  return std::numeric_limits<double>::quiet_NaN();
}

struct compare_gcs_paths {
  bool operator() (std::shared_ptr<GCSPath> &path1,
                   std::shared_ptr<GCSPath> &path2) {
    if (path1->f_val > path2->f_val) {
      return true;
    }
    if (path1->f_val < path2->f_val) {
      return false;
    }
    return path1->tie_break_val > path2->tie_break_val;
  }
};


MTIxGStar::MTIxGStar(const std::vector<GCSNode> &gcs_nodes, const Ref<const Matrix<bool, Dynamic, Dynamic, RowMajor>> &gcs_adj_mat, int dim_q, const Ref<const VectorXd> &set_lb, const Ref<const VectorXd> &set_ub, double vmax_agent) : gcs_nodes(gcs_nodes),
                                                                                                                                                                                                                                          gcs_adj_mat(gcs_adj_mat),
                                                                                                                                                                                                                                          dim_q(dim_q),
                                                                                                                                                                                                                                          set_lb(set_lb), 
                                                                                                                                                                                                                                          set_ub(set_ub),
                                                                                                                                                                                                                                          vmax_agent(vmax_agent) {
}

MatrixXd MTIxGStar::solve(int start_idx, int goal_idx, double w, VectorXd start_pos_and_time, double max_final_time, const Ref<const VectorXl> &ignore_set_vec, bool verbose, double time_limit, double tw_start) {
  auto ixg_start_time = std::chrono::high_resolution_clock::now();

  std::priority_queue<std::shared_ptr<GCSPath>, std::vector<std::shared_ptr<GCSPath>>, compare_gcs_paths> open_list;
  TrajoptThroughConvexSets trajopt(dim_q, set_lb, set_ub, start_pos_and_time, vmax_agent);
  trajopt.add_convex_set(gcs_nodes[start_idx]);
  std::vector<int> node_seq;
  node_seq.push_back(start_idx);
  int tie_break_val = 0;
  open_list.push(std::make_shared<GCSPath>(node_seq, -1, 0, 0, 0, 0, tie_break_val, trajopt));
  tie_break_val += 1;
  std::vector<MatrixXd> pos_seq_list;
  std::vector<VectorXd> time_seq_list;
  std::unordered_set<int> ignore_set;
  for (int i = 0; i < ignore_set_vec.size(); ++i) {
    if (ignore_set_vec(0) == -1) {
      break;
    }
    ignore_set.insert(ignore_set_vec(i));
  }

  while (open_list.size() != 0) {
    auto cur_time = std::chrono::high_resolution_clock::now();
    auto millis = std::chrono::duration_cast<std::chrono::milliseconds>(cur_time - ixg_start_time).count();
    if (millis >= time_limit*1000) {
      // Stop planning!
      break;
    }

    if (time_seq_list.size() && time_seq_list.back().tail<1>()(0) <= w*open_list.top()->f_val) {
      MatrixXd ret(pos_seq_list.back().rows(), dim_q + 1);
      ret.leftCols(dim_q) = pos_seq_list.back();
      ret.col(dim_q) = time_seq_list.back();
      return ret;
    }

    std::shared_ptr<GCSPath> path = open_list.top();
    open_list.pop();
    if (verbose) {
      for (auto node_idx : path->node_seq) {
        std::cout << node_idx << " ";
      }
      std::cout << " cost_lb: " << path->cost_lb << " g " << path->g_val << " h " << path->h_val;
      std::cout << std::endl;
    }

    for (int next_node_idx = 0; next_node_idx < gcs_nodes.size(); ++next_node_idx) {
      if (!gcs_adj_mat(path->node_seq.back(), next_node_idx) || 
          ignore_set.find(next_node_idx) != ignore_set.end() ||
          std::find(path->node_seq.begin(), path->node_seq.end(), next_node_idx) != path->node_seq.end()) {
        continue;
      }

      std::shared_ptr<GCSPath> neighbor = extend_path(path, next_node_idx, goal_idx, w, max_final_time, pos_seq_list, time_seq_list);

      if (neighbor == nullptr) {
        continue;
      }

      neighbor->tie_break_val = tie_break_val;
      tie_break_val += 1;

      open_list.push(neighbor);

      if (next_node_idx == goal_idx) {
        max_final_time = time_seq_list.back().tail<1>()[0];

        if (verbose) {
          std::cout << "Goal is successor, cost = " << max_final_time << std::endl;
        }

        /*
        if (max_final_time < tw_start + 1e-8) {
          MatrixXd ret(pos_seq_list.back().rows(), dim_q + 1);
          ret.leftCols(dim_q) = pos_seq_list.back();
          ret.col(dim_q) = time_seq_list.back();
          return ret;
        }
        */
      }
    }
  }

  if (verbose) {
    std::cout << "open empty" << std::endl;
  }
  return nan_num()*MatrixXd::Ones(1, 1);
}

std::shared_ptr<GCSPath> MTIxGStar::extend_path(std::shared_ptr<GCSPath> path, int next_node_idx, int goal_node_idx, double w, double max_final_time, std::vector<MatrixXd> &pos_seq_list, std::vector<VectorXd> &time_seq_list) {
  std::shared_ptr<GCSPath> ret = std::make_shared<GCSPath>(path->trajopt);
  ret->node_seq = path->node_seq;
  ret->node_seq.push_back(next_node_idx);
  ret->trajopt.add_convex_set(gcs_nodes[next_node_idx]);

  MatrixXd A_eq;
  VectorXd b_eq;
  MatrixXd A_ineq;
  VectorXd b_ineq;
  if (next_node_idx != goal_node_idx) {
    if (std::isnan(gcs_nodes[next_node_idx].A_eq_pt(0, 0)) && std::isnan(gcs_nodes[goal_node_idx].A_eq_pt(0, 0))) {
      A_eq = nan_num()*MatrixXd::Ones(1, 1);
      b_eq = nan_num()*VectorXd::Ones(1);
    } else if (std::isnan(gcs_nodes[next_node_idx].A_eq_pt(0, 0))) {
      A_eq = MatrixXd(gcs_nodes[goal_node_idx].A_eq_pt.rows(), gcs_nodes[goal_node_idx].A_eq_pt.cols()*2);
      A_eq.leftCols(gcs_nodes[goal_node_idx].A_eq_pt.cols()).setZero();
      A_eq.rightCols(gcs_nodes[goal_node_idx].A_eq_pt.cols()) = gcs_nodes[goal_node_idx].A_eq_pt;
      b_eq = gcs_nodes[goal_node_idx].b_eq_pt;
    } else if (std::isnan(gcs_nodes[goal_node_idx].A_eq_pt(0, 0))) {
      A_eq = MatrixXd(gcs_nodes[next_node_idx].A_eq_pt.rows(), gcs_nodes[next_node_idx].A_eq_pt.cols()*2);
      A_eq.leftCols(gcs_nodes[next_node_idx].A_eq_pt.cols()) = gcs_nodes[next_node_idx].A_eq_pt;
      A_eq.rightCols(gcs_nodes[next_node_idx].A_eq_pt.cols()).setZero();
      b_eq = gcs_nodes[next_node_idx].b_eq_pt;
    } else {
      A_eq = MatrixXd(gcs_nodes[next_node_idx].A_eq_pt.rows() + gcs_nodes[goal_node_idx].A_eq_pt.rows(), gcs_nodes[next_node_idx].A_eq_pt.cols()*2);
      A_eq.topRightCorner(gcs_nodes[next_node_idx].A_eq_pt.rows(), gcs_nodes[next_node_idx].A_eq_pt.cols()).setZero();
      A_eq.bottomLeftCorner(gcs_nodes[next_node_idx].A_eq_pt.rows(), gcs_nodes[next_node_idx].A_eq_pt.cols()).setZero();
      A_eq.topLeftCorner(gcs_nodes[next_node_idx].A_eq_pt.rows(), gcs_nodes[next_node_idx].A_eq_pt.cols()) = gcs_nodes[next_node_idx].A_eq_pt;
      A_eq.bottomRightCorner(gcs_nodes[goal_node_idx].A_eq_pt.rows(), gcs_nodes[goal_node_idx].A_eq_pt.cols()) = gcs_nodes[goal_node_idx].A_eq_pt;
      b_eq = VectorXd(gcs_nodes[next_node_idx].b_eq_pt.rows() + gcs_nodes[goal_node_idx].b_eq_pt.rows());
      b_eq.head(gcs_nodes[next_node_idx].b_eq_pt.rows()) = gcs_nodes[next_node_idx].b_eq_pt;
      b_eq.tail(gcs_nodes[goal_node_idx].b_eq_pt.rows()) = gcs_nodes[goal_node_idx].b_eq_pt;
    }

    if (std::isnan(gcs_nodes[next_node_idx].A_ineq_pt(0, 0)) && std::isnan(gcs_nodes[goal_node_idx].A_ineq_pt(0, 0))) {
      A_ineq = nan_num()*MatrixXd::Ones(1, 1);
      b_ineq = nan_num()*VectorXd::Ones(1);
    } else if (std::isnan(gcs_nodes[next_node_idx].A_ineq_pt(0, 0))) {
      A_ineq = MatrixXd(gcs_nodes[goal_node_idx].A_ineq_pt.rows(), gcs_nodes[goal_node_idx].A_ineq_pt.cols()*2);
      A_ineq.leftCols(gcs_nodes[goal_node_idx].A_ineq_pt.cols()).setZero();
      A_ineq.rightCols(gcs_nodes[goal_node_idx].A_ineq_pt.cols()) = gcs_nodes[goal_node_idx].A_ineq_pt;
      b_ineq = gcs_nodes[goal_node_idx].b_ineq_pt;
    } else if (std::isnan(gcs_nodes[goal_node_idx].A_ineq_pt(0, 0))) {
      A_ineq = MatrixXd(gcs_nodes[next_node_idx].A_ineq_pt.rows(), gcs_nodes[next_node_idx].A_ineq_pt.cols()*2);
      A_ineq.leftCols(gcs_nodes[next_node_idx].A_ineq_pt.cols()) = gcs_nodes[next_node_idx].A_ineq_pt;
      A_ineq.rightCols(gcs_nodes[next_node_idx].A_ineq_pt.cols()).setZero();
      b_ineq = gcs_nodes[next_node_idx].b_ineq_pt;
    } else {
      A_ineq = MatrixXd(gcs_nodes[next_node_idx].A_ineq_pt.rows() + gcs_nodes[goal_node_idx].A_ineq_pt.rows(), gcs_nodes[next_node_idx].A_ineq_pt.cols()*2);
      A_ineq.topRightCorner(gcs_nodes[next_node_idx].A_ineq_pt.rows(), gcs_nodes[next_node_idx].A_ineq_pt.cols()).setZero();
      A_ineq.bottomLeftCorner(gcs_nodes[next_node_idx].A_ineq_pt.rows(), gcs_nodes[next_node_idx].A_ineq_pt.cols()).setZero();
      A_ineq.topLeftCorner(gcs_nodes[next_node_idx].A_ineq_pt.rows(), gcs_nodes[next_node_idx].A_ineq_pt.cols()) = gcs_nodes[next_node_idx].A_ineq_pt;
      A_ineq.bottomRightCorner(gcs_nodes[goal_node_idx].A_ineq_pt.rows(), gcs_nodes[goal_node_idx].A_ineq_pt.cols()) = gcs_nodes[goal_node_idx].A_ineq_pt;
      b_ineq = VectorXd(gcs_nodes[next_node_idx].b_ineq_pt.rows() + gcs_nodes[goal_node_idx].b_ineq_pt.rows());
      b_ineq.head(gcs_nodes[next_node_idx].b_ineq_pt.rows()) = gcs_nodes[next_node_idx].b_ineq_pt;
      b_ineq.tail(gcs_nodes[goal_node_idx].b_ineq_pt.rows()) = gcs_nodes[goal_node_idx].b_ineq_pt;
    }

    TrajoptThroughConvexSets trajopt_f(ret->trajopt);
    GCSNode final_convex_set(A_eq, b_eq, A_ineq, b_ineq, nan_num()*MatrixXd::Ones(1, 1), nan_num()*VectorXd::Ones(1), nan_num()*MatrixXd::Ones(1, 1), nan_num()*VectorXd::Ones(1), true);
    trajopt_f.add_convex_set(final_convex_set);

    if (!std::isinf(max_final_time)) {
      if (std::isnan(A_ineq(0, 0))) {
        trajopt_f.linear_cone_b.back()[trajopt_f.num_dist_cons_rows + dim_q + 1 + dim_q] = max_final_time;
      } else {
        trajopt_f.linear_cone_b[trajopt_f.linear_cone_b.size() - 2][trajopt_f.num_dist_cons_rows + dim_q + 1 + dim_q] = max_final_time;
      }
    }
    std::vector<VectorXd> extended_pos_seq;
    std::vector<double> extended_time_seq;
    double cost = trajopt_f.solve_trajopt(extended_pos_seq, extended_time_seq, true, false);

    if (std::isinf(cost)) {
      return nullptr;
    }
    ret->g_val = extended_time_seq[extended_time_seq.size() - 2];
    ret->h_val = extended_time_seq.back() - ret->g_val;
    ret->cost_lb = extended_time_seq.back();
    ret->f_val = ret->g_val + w*(ret->h_val);
    return ret;
  } else {
    if (!std::isinf(max_final_time)) {
      if (std::isnan(gcs_nodes[next_node_idx].A_ineq(0, 0))) {
        ret->trajopt.linear_cone_b.back()[ret->trajopt.num_dist_cons_rows + dim_q + 1 + dim_q] = max_final_time;
      } else {
        ret->trajopt.linear_cone_b[ret->trajopt.linear_cone_b.size() - 2][ret->trajopt.num_dist_cons_rows + dim_q + 1 + dim_q] = max_final_time;
      }
    }
    std::vector<VectorXd> extended_pos_seq;
    std::vector<double> extended_time_seq;
    double cost = ret->trajopt.solve_trajopt(extended_pos_seq, extended_time_seq, true, false);
    if (std::isinf(cost)) {
      return nullptr;
    }
    MatrixXd extended_pos_seq_mat(extended_pos_seq.size(), dim_q);
    for (int i = 0; i < extended_pos_seq.size(); ++i) {
      extended_pos_seq_mat.row(i) = extended_pos_seq[i];
    }
    pos_seq_list.clear();
    time_seq_list.clear();
    pos_seq_list.push_back(extended_pos_seq_mat);
    time_seq_list.push_back(Map<VectorXd>(extended_time_seq.data(), extended_time_seq.size()));
    ret->g_val = extended_time_seq.back();
    ret->h_val = 0.;
    ret->cost_lb = extended_time_seq.back();
    ret->f_val = ret->g_val;
    ret->seq_list_idx = pos_seq_list.size() - 1;
    return ret;
  }
}
