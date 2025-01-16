#pragma once
#include <Eigen/Dense>
#include <random>
#include <memory>
#include <chrono>
#include <limits>
#include <iostream>
#include "mt_tsp_ros2/arastar_data.h"

using namespace Eigen;
using namespace std::chrono;

// We assume the successor generator only generates valid transitions
// We also assume the start node is valid.
// If we take more milliseconds that acceptable_millis, we stop as
// soon as we find a feasible solution, even if it's suboptimal.
// If we take more milliseconds than timeout_millis, or if we exhaust the
// open list, we stop whether we've found a feasible solution or not, and return
// the partial solution ending with the best h-value.
bool arastar(AStarPath &path, std::shared_ptr<ARAStarData> data,
             double start_eps,
             int acceptable_millis, int timeout_millis, double end_eps = 1.0) {
  auto start_time = std::chrono::high_resolution_clock::now();

  for (double eps = start_eps; eps >= end_eps; eps = (eps == end_eps) ? 0 : std::max(end_eps, eps / 2)) {
    data->reset(eps);
    while (!data->open_list_empty() &&
           data->get_f_goal() >
           data->get_next()->get_f()) {
      AStarNodePtr pop = data->expand_next();
      double g = pop->get_g();

      std::vector<AStarCell> succ;
      std::vector<double> transition_costs;
      data->generate_successors(succ, transition_costs, pop);
      for (int i = 0; i < succ.size(); ++i) {
        data->update_path_to_node(pop, succ[i], transition_costs[i]);
      }

      auto stop_time = std::chrono::high_resolution_clock::now();
      auto millis = std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time).count();
      if (eps != start_eps && millis >= acceptable_millis) {
        // Stop planning!
        return true;
      }

      if (millis >= timeout_millis) {
        // std::cout << "A* timed out. Returning partial solution with the best terminal h-value" << std::endl;
        data->get_path(path);
        return false;
      }
    } 

    if (data->open_list_empty()) {
      // std::cout << "No path. Returning partial solution with the best terminal h-value" << std::endl;
      data->get_path(path);
      return false;
    }

    data->get_path(path);

    auto stop_time = std::chrono::high_resolution_clock::now();
    auto millis = std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time).count();
    if (millis >= acceptable_millis) {
      // Stop planning!
      return true;
    }
  }
  return true;
}
