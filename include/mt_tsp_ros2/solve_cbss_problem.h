#pragma once

#include <pybind11/pybind11.h>
#include "mt_tsp_ros2/cbss_data.h"
#include "mt_tsp_ros2/cbss.h"

namespace py = pybind11;

MatrixXi solve_cbss_problem(const Ref<const Matrix<bool, Dynamic, Dynamic, RowMajor>> &occupancy, int num_agents, py::object k_best_mamttsp_solver, int timeout_millis, int low_level_planner_timeout_millis,
                            Ref<RowMatrixXl> high_level_cell_seqs,
                            Ref<VectorXl> high_level_time_seqs,
                            Ref<VectorXl> agent_cell_seq_start_ptr) {
  std::cout << "making problem" << std::endl;
  std::shared_ptr<CBSSProblem> problem = std::make_shared<CBSSProblem>(occupancy, num_agents, low_level_planner_timeout_millis);
  std::cout << "making data" << std::endl;
  std::shared_ptr<CBSSData> data = std::make_shared<CBSSData>(problem, k_best_mamttsp_solver);
  std::cout << "running cbss" << std::endl;
  AStarPath path;
  if (cbss(path, data, 10, timeout_millis, high_level_cell_seqs, high_level_time_seqs, agent_cell_seq_start_ptr)) {
    MatrixXi grid_path;
    problem->get_grid_path(grid_path, path);
    return grid_path;
  } else {
    std::cout << "No solution. Returning empty path" << std::endl;
    MatrixXi grid_path(0, problem->get_num_agents()*2);
    return grid_path;
  }
}
