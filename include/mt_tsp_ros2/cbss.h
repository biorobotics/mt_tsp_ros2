#pragma once
#include <Eigen/Dense>
#include <random>
#include <memory>
#include <chrono>
#include <limits>
#include <iostream>
#include "mt_tsp_ros2/cbss_data.h"
#include <fstream>

using namespace Eigen;
using namespace std::chrono;

// We assume the successor generator only generates valid transitions
// We also assume the start node is valid.
// If we take more milliseconds that acceptable_millis, we stop as
// soon as we find a feasible solution, even if it's suboptimal.
// If we take more milliseconds than timeout_millis, or if we exhaust the
// open list, we stop whether we've found a feasible solution or not, and
// return an empty path
bool cbss(AStarPath &path,
          std::shared_ptr<CBSSData> data,
          int acceptable_millis,
          int timeout_millis,
          Ref<RowMatrixXl> high_level_cell_seqs,
          Ref<VectorXl> high_level_time_seqs,
          Ref<VectorXl> agent_cell_seq_start_ptr) {
  auto start_time = std::chrono::high_resolution_clock::now();

  int num_exp = 0;

  CBSSNodePtr soln = nullptr;
  while (!data->open_list_empty()) {
    CBSSNodePtr pop = data->expand_next();
    CBSConstraintList new_constraints = data->get_new_constraints_from_conflicts(pop);
    if (new_constraints.size() == 0) {
      soln = pop;
      break;
    }

    for (auto new_constraint : new_constraints) {
      data->add_constraint(pop, new_constraint);
    }

    ++num_exp;

    auto stop_time = std::chrono::high_resolution_clock::now();
    auto millis = std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time).count();

    if (millis >= timeout_millis) {
      std::cout << "Timed out. Returning empty path" << std::endl;
      std::cout << "Number of expansions: " << num_exp << std::endl;
      return false;
    }
  } 

  if (soln == nullptr) {
    std::cout << "No path. Returning empty path" << std::endl;
    std::cout << "Number of expansions: " << num_exp << std::endl;
    return false;
  }

  std::cout << "getting path" << std::endl;
  soln->get_path_through_joint_state_space(path);
  std::cout << "got path" << std::endl;
  std::cout << "Number of expansions: " << num_exp << std::endl;

  high_level_cell_seqs = soln->get_ma_mt_tsp_soln()->get_cell_seqs();
  high_level_time_seqs = soln->get_ma_mt_tsp_soln()->get_time_seqs();
  agent_cell_seq_start_ptr = soln->get_ma_mt_tsp_soln()->get_agent_cell_seq_start_ptr();

  return true;
}
