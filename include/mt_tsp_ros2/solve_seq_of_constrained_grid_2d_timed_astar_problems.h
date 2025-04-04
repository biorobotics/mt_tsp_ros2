#include "mt_tsp_ros2/constrained_grid_2d_timed_astar_problem.h"
#include "mt_tsp_ros2/cbs_node.h"
#include "mt_tsp_ros2/arastar.h"
#include <omp.h>

bool solve_seq_of_constrained_grid_2d_timed_astar_problems(GridPathPtr path, const Matrix<bool, Dynamic, Dynamic> &occupancy, const AStarConstraintList &constraints, Ref<RowMatrixXl> cell_seq, Ref<VectorXl> time_seq, double low_level_planner_timeout_millis, long T) {
  int dim_q = 2;

  int num_partial_paths = std::max((int)(time_seq.size() - 1), 1);
  std::vector<GridPathPtr> partial_paths(num_partial_paths);

  volatile bool infeas = false;

  #pragma omp parallel for
  for (int i = 0; i < num_partial_paths; ++i) {
    if (infeas) {
      continue;
    }
    // Instantiate astar problem with constraints
    std::shared_ptr<ConstrainedGrid2dTimedAStarProblem> problem;
    if (time_seq.size() == 1) {
      problem = std::make_shared<ConstrainedGrid2dTimedAStarProblem>(occupancy, constraints, cell_seq(i, 0), cell_seq(i, 1), T);
    } else {
      problem = std::make_shared<ConstrainedGrid2dTimedAStarProblem>(occupancy, constraints, cell_seq(i + 1, 0), cell_seq(i + 1, 1), time_seq(i + 1));
    }
    VectorXi time_aug_start_cell(dim_q + 1);
    time_aug_start_cell.head(dim_q) = cell_seq.row(i).transpose().cast<int>();
    time_aug_start_cell(dim_q) = time_seq(i);
    std::shared_ptr<ARAStarData> data = std::make_shared<ARAStarData>(problem, time_aug_start_cell);

    // Run astar
    AStarPath astar_path;
    if (!arastar(astar_path, data, 1.0, 10, low_level_planner_timeout_millis)) {
      infeas = true;
      continue;
    }
    partial_paths[i] = std::make_shared<GridPath>();
    problem->get_grid_path(partial_paths[i]->path, astar_path);
  }
  if (infeas) {
    return false;
  }

  // The -1's are important so we don't repeat time steps
  int total_path_length = 1;
  for (auto partial_path : partial_paths) {
    total_path_length += partial_path->path.rows() - 1;
  }

  path->cost = 0;
  path->path = MatrixXi::Zero(total_path_length, dim_q);
  path->path.topRows(1) = partial_paths[0]->path.topRows(1);
  int length_so_far = 1;
  for (auto partial_path : partial_paths) {
    path->path.block(length_so_far, 0, partial_path->path.rows() - 1, dim_q) = partial_path->path.bottomRows(partial_path->path.rows() - 1);
    length_so_far += partial_path->path.rows() - 1;
    path->cost += partial_path->cost;
  }
  return true;
}
