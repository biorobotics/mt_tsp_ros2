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
      
      double h;
      if (node_idx == -1) {
        h = 0;
      } else {
        h = std::numeric_limits<double>::infinity();
        for (auto goal_idx : goal_indices) {
          h = std::min(h, cost_mat(node_idx, goal_idx));
        }
      }
      return h;
    }

  private:
    VectorXl start_indices;
    VectorXl goal_indices;
    RowMatrixXd masked_cost_mat;
    RowMatrixXd cost_mat;
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
    VectorXl node_seq(path.size());
    for (int i = 0; i < path.size(); ++i) {
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

typedef VectorXi TargetSeq;

// Taken from https://jimmy-shen.medium.com/stl-map-unordered-map-with-a-vector-for-the-key-f30e5f670bae#:~:text=unordered_map%20uses%20vector%20as%20the%20key&text=You%20can%20use%20the%20following,make%20the%20best%20of%20STL.&text=%7D%3B,so%20that%20collisions%20are%20minimized
struct TargetSeqHash {
  int operator()(const TargetSeq &V) const {
    int hash = V.size();
    for(int i = 0; i < V.size(); ++i) {
      hash ^= V[i] + 0x9e3779b9 + (hash << 6) + (hash >> 2);
    }
    return hash;
  }
};

class DTPAStarSolver {
  public:
    DTPAStarSolver(const Ref<const RowMatrixXd> &cost_mat) : cost_mat(cost_mat) {
    }

    VectorXl solve(Ref<Matrix<double, 1, 1>> cost,
                   const Ref<const VectorXl> &start_indices,
                   const Ref<const VectorXl> &goal_indices,
                   const Ref<const RowMatrixXd> &masked_cost_mat,
                   const Ref<const VectorXl> &target_seq) {
      std::shared_ptr<DTPAStarProblem> problem =  std::make_shared<DTPAStarProblem>(start_indices,
                                                                                    goal_indices, 
                                                                                    masked_cost_mat, 
                                                                                    cost_mat);
      std::shared_ptr<ARAStarData> astar_data;
      VectorXi target_seq_int = target_seq.cast<int>();
      // assert(stored_data.find(target_seq_int) == stored_data.end());
      std::unordered_map<TargetSeq, std::shared_ptr<ARAStarData>, TargetSeqHash>::iterator it;
      for (int size = target_seq_int.size() - 1; size >= 1; --size) {
        it = stored_data.find(target_seq_int.head(target_seq_int.size() - 1));
        if (it != stored_data.end()) {
          break;
        }
      }
      if (it == stored_data.end()) {
        AStarCell start_cell = VectorXi(1);
        start_cell(0) = -1;
        astar_data = std::make_shared<ARAStarData>(problem, start_cell);
      } else {
        astar_data = it->second->copy();
        astar_data->update_problem(problem);
      }

      stored_data[target_seq_int] = astar_data;

      AStarPath path;
      if (arastar(path, astar_data, 1.0, 10, 10000)) {
        VectorXl node_seq(path.size());
        for (int i = 0; i < path.size(); ++i) {
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
  private:
    RowMatrixXd cost_mat;
    std::unordered_map<TargetSeq, std::shared_ptr<ARAStarData>, TargetSeqHash> stored_data;
};
