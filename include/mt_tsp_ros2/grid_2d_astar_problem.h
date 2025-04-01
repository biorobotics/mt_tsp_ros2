#pragma once
#include "mt_tsp_ros2/astar_problem.h"
#include <memory>
#include "mt_tsp_ros2/arastar.h"

class Grid2dAStarProblem : public AStarProblem {
  public:
    Grid2dAStarProblem(const Ref<const Matrix<bool, Dynamic, Dynamic, RowMajor>> &occupancy, int gx, int gy) : occupancy(occupancy), gx(gx), gy(gy) {
    }

    virtual void generate_successors(std::vector<AStarCell> &succ, std::vector<double> &transition_costs, const AStarNodePtr &node) const override {
      succ.clear();
      transition_costs.clear();

      const VectorXi &coords = node->get_cell();
      VectorXi next_coords(2);

      // +x
      if (coords(0) + 1 < occupancy.rows() && !occupancy(coords(0) + 1, coords(1))) {
        next_coords(0) = coords(0) + 1;
        next_coords(1) = coords(1);
        succ.push_back(next_coords);
        transition_costs.push_back(1);
      }
      // -x
      if (coords(0) - 1 >= 0 && !occupancy(coords(0) - 1, coords(1))) {
        next_coords(0) = coords(0) - 1;
        next_coords(1) = coords(1);
        succ.push_back(next_coords);
        transition_costs.push_back(1);
      }
      // +y
      if (coords(1) + 1 < occupancy.cols() && !occupancy(coords(0), coords(1) + 1)) {
        next_coords(0) = coords(0);
        next_coords(1) = coords(1) + 1;
        succ.push_back(next_coords);
        transition_costs.push_back(1);
      }
      // -y
      if (coords(1) - 1 >= 0 && !occupancy(coords(0), coords(1) - 1)) {
        next_coords(0) = coords(0);
        next_coords(1) = coords(1) - 1;
        succ.push_back(next_coords);
        transition_costs.push_back(1);
      }
    }

    virtual bool is_goal(const AStarNodePtr &node) const override {
      const VectorXi &coords = node->get_cell();
      return coords(0) == gx && coords(1) == gy;
    }

    void get_grid_path(MatrixXi &grid_path, const AStarPath &path) {
      grid_path = MatrixXi::Zero(path.size(), 2);
      for (int i = 0; i < path.size(); ++i) {
        grid_path.row(i) = path[i].transpose();
      }
    }

    virtual double heuristic(const AStarCell &cell) const override {
      return abs(cell(0) - gx) + abs(cell(1) - gy);
    }

  private:
    Matrix<bool, Dynamic, Dynamic, RowMajor> occupancy;
    int gx;
    int gy;
};

MatrixXi solve_grid_2d_astar_problem(const Ref<const Matrix<bool, Dynamic, Dynamic, RowMajor>> &occupancy, int sx, int sy, int gx, int gy) {
  std::shared_ptr<Grid2dAStarProblem> problem = std::make_shared<Grid2dAStarProblem>(occupancy, gx, gy);
  AStarCell start_cell = Vector2i(sx, sy);
  std::shared_ptr<ARAStarData> data = std::make_shared<ARAStarData>(problem, start_cell);
  AStarPath path;
  if (arastar(path, data, 1.0, 10, 10000)) {
    MatrixXi grid_path;
    problem->get_grid_path(grid_path, path);
    return grid_path;
  } else {
    std::cout << "No solution. Returning path where agent stays at start config" << std::endl;
    return RowVector2i(sx, sy);
  }
}
