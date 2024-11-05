#pragma once
#include "mt_tsp_ros2/astar_problem.h"
#include <memory>
#include "mt_tsp_ros2/arastar.h"

typedef Matrix<long, Dynamic, 1> VectorXl;
typedef Matrix<double, Dynamic, Dynamic, RowMajor> RowMatrixXd;

class DTPAStarProblem : public AStarProblem {
  public:
    DTPAStarProblem(const Ref<const VectorXl> &start_indices,
                    const Ref<const VectorXl> &goal_indices,
                    const Ref<const RowMatrixXd> &masked_cost_mat,
                    const Ref<const RowMatrixXd> &cost_mat) : start_indices(start_indices), 
                                                              goal_indices(goal_indices), 
                                                              masked_cost_mat(masked_cost_mat),
                                                              cost_mat(cost_mat) {
      
    }

    virtual void generate_successors(std::vector<AStarCell> &succ, std::vector<double> &transition_costs, const AStarNodePtr &node) const override {
      succ.clear();
      transition_costs.clear();

      int node_idx = node->get_cell()(0);
      if (is_goal_node_idx(node_idx)) {
        return;
      }

      VectorXi next_cell(1);

      if (node_idx == -1) {
        if (start_indices.size() == 1 && start_indices(0) == 0) {
          // Solving full DTP where we're returning to start q
          for (int next_node_idx = 0; next_node_idx < cost_mat.cols(); ++next_node_idx) {
            if (std::isinf(masked_cost_mat(0, next_node_idx))) {
              continue;
            }
            next_cell(0) = next_node_idx;
            succ.push_back(next_cell);
            transition_costs.push_back(masked_cost_mat(0, next_node_idx));
          }
        } else {
          for (auto next_node_idx : start_indices) {
            next_cell(0) = next_node_idx;
            succ.push_back(next_cell);
            transition_costs.push_back(0.);
          }
        }
      } else {
        for (int next_node_idx = 0; next_node_idx < cost_mat.cols(); ++next_node_idx) {
          if (std::isinf(masked_cost_mat(node_idx, next_node_idx))) {
            continue;
          }
          next_cell(0) = next_node_idx;
          succ.push_back(next_cell);
          transition_costs.push_back(masked_cost_mat(node_idx, next_node_idx));
        }
      }
    }

    bool is_goal_node_idx(int node_idx) const {
      return std::find(goal_indices.begin(), goal_indices.end(), node_idx) != goal_indices.end();
    }

    virtual bool is_goal(const AStarNodePtr &node) const override {
      return is_goal_node_idx(node->get_cell()(0));
    }

    virtual double heuristic(const AStarCell &cell) const override {
      int node_idx = cell(0);
      if (is_goal_node_idx(node_idx)) {
        return 0;
      }
      
      double h = 0;
      if (node_idx != -1) {
        for (auto goal_idx : goal_indices) {
          h = std::min(h, cost_mat(node_idx, goal_idx));
        }
      }
      return h;
    }

  private:
    const Ref<const VectorXl> &start_indices;
    const Ref<const VectorXl> &goal_indices;
    const Ref<const RowMatrixXd> &masked_cost_mat;
    const Ref<const RowMatrixXd> &cost_mat;
};

VectorXl solve_dtp_astar_problem(Ref<Matrix<double, 1, 1>> cost,
                                 const Ref<const VectorXl> &start_indices,
                                 const Ref<const VectorXl> &goal_indices,
                                 const Ref<const RowMatrixXd> &masked_cost_mat,
                                 const Ref<const RowMatrixXd> &cost_mat) {
  std::shared_ptr<DTPAStarProblem> problem =  std::make_shared<DTPAStarProblem>(start_indices,
                                                                                goal_indices, 
                                                                                masked_cost_mat, 
                                                                                cost_mat);

  AStarCell start_cell = VectorXi(1);
  start_cell(0) = -1;
  std::shared_ptr<ARAStarData> astar_data = std::make_shared<ARAStarData>(problem, start_cell);
  AStarPath path;
  if (arastar(path, astar_data, 1.0, 10, 10000)) {
    VectorXl node_seq(path.size() - 1);
    for (int i = 0; i < path.size() - 1; ++i) {
      node_seq(i) = path[i](0);
    }
    cost(0) = astar_data->get_f_goal();
    return node_seq;
  } else {
    // std::cout << "No solution. Returning path where agent stays at start config" << std::endl;
    VectorXl node_seq(1);
    node_seq(0) = -1;
    cost(0) = std::numeric_limits<double>::infinity();
    return node_seq;
  }
}
