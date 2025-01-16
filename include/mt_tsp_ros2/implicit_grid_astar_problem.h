#pragma once
#include "mt_tsp_ros2/astar_problem.h"
#include <memory>
#include "mt_tsp_ros2/arastar.h"
#include <pybind11/pybind11.h>

namespace py = pybind11;

typedef Matrix<long, Dynamic, 1> VectorXl;
typedef Matrix<long, Dynamic, Dynamic, RowMajor> RowMatrixXl;

class ImplicitGridAStarProblem : public AStarProblem {
  public:
    ImplicitGridAStarProblem(const Ref<const VectorXi> &goal_cell,
                             py::object collision_checker,
                             const Ref<const VectorXi> &num_cells,
                             const Ref<const VectorXd> &cell_size,
                             double cost_bound) : goal_cell(goal_cell), collision_checker(collision_checker), num_cells(num_cells), cell_size(cell_size), cost_bound(cost_bound) {
      
    }

    virtual void generate_successors(std::vector<AStarCell> &succ, std::vector<double> &transition_costs, const AStarNodePtr &node) const override {
      succ.clear();
      transition_costs.clear();

      const Ref<const VectorXi> &cur_cell = node->get_cell();

      if (is_goal(node)) {
        return;
      }

      int dim_q = goal_cell.size();
      VectorXi next_cell = cur_cell;

      for (int q_idx = 0; q_idx < dim_q; ++q_idx) {
        for (int direction = -1; direction <= 1; direction += 2) {
          if ((cur_cell(q_idx) == 0 && direction == -1) || (cur_cell(q_idx) == num_cells(q_idx) && direction == 1)) {
            continue;
          }

          next_cell(q_idx) = cur_cell(q_idx) + direction;

          if (collision_checker.attr("CheckConfigCollisionFree")(next_cell) && node->get_g() + cell_size(q_idx) + heuristic(next_cell) <= cost_bound) {
            succ.push_back(next_cell);
            transition_costs.push_back(cell_size(q_idx));
          }

          next_cell(q_idx) = cur_cell(q_idx);
        }
      }
    }

    virtual bool is_goal(const AStarNodePtr &node) const override {
      return (goal_cell - node->get_cell()).lpNorm<1>() == 0;
    }

    virtual double heuristic(const AStarCell &cell) const override {
      return (goal_cell - cell).cast<double>().cwiseProduct(cell_size).lpNorm<1>();
    }

  private:
    VectorXi goal_cell;
    py::object collision_checker;
    VectorXi num_cells;
    VectorXd cell_size;
    double cost_bound;
};

RowMatrixXl solve_implicit_grid_astar_problem(Ref<Matrix<double, 1, 1>> cost,
                                              const Ref<const VectorXl> &start_cell,
                                              const Ref<const VectorXl> &goal_cell,
                                              py::object collision_checker,
                                              const Ref<const VectorXl> &num_cells,
                                              const Ref<const VectorXd> &cell_size,
                                              double cost_bound,
                                              double time_limit,
                                              double w) {
  std::shared_ptr<ImplicitGridAStarProblem> problem =  std::make_shared<ImplicitGridAStarProblem>(goal_cell.cast<int>(),
                                                                                                  collision_checker, 
                                                                                                  num_cells.cast<int>(), 
                                                                                                  cell_size, 
                                                                                                  cost_bound);

  std::shared_ptr<ARAStarData> astar_data = std::make_shared<ARAStarData>(problem, start_cell.cast<int>());
  AStarPath path;
  if (arastar(path, astar_data, w, time_limit*1000, time_limit*1000, w)) {
    RowMatrixXl cell_seq(path.size(), num_cells.size());
    for (int i = 0; i < path.size(); ++i) {
      cell_seq.row(i) = path[i].transpose().cast<long>();
    }
    cost(0) = astar_data->get_f_goal();
    return cell_seq;
  } else {
    // std::cout << "No solution. Returning path where agent stays at start config" << std::endl;
    RowMatrixXl cell_seq(1, num_cells.size());
    cell_seq.row(0) = start_cell.transpose();
    cost(0) = std::numeric_limits<double>::infinity();
    return cell_seq;
  }
}
