#pragma once
#include <memory>
#include <unordered_set>
#include <set>
#include <vector>
#include <chrono>
#include <Eigen/Dense>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

using namespace Eigen;
namespace py = pybind11;

typedef const Ref<const Matrix<long, Dynamic, 1>> &VectorXlRef_const;
typedef const Ref<const Matrix<double, Dynamic, Dynamic, RowMajor>> &RowMatrixXdRef_const;
typedef const Ref<const Matrix<bool, Dynamic, Dynamic, RowMajor>> &RowMatrixXbRef_const;
typedef Matrix<bool, Dynamic, Dynamic> MatrixXb;
typedef Matrix<bool, Dynamic, 1> VectorXb;
typedef Matrix<long, Dynamic, 1> VectorXl;

// bool reopt_gtsp_tour(Ref<VectorXl> tour, RowMatrixXdRef_const gtsp_cost_mat, VectorXlRef_const pt_to_target_ptr, const std::vector<py::array_t<long>> &target_to_pt_ptr, const Ref<const VectorXl> &target_seq, Ref<Matrix<long, 1, 1>> longest_feasible_prefix) {
bool reopt_gtsp_tour(Ref<VectorXl> tour, RowMatrixXdRef_const gtsp_cost_mat, VectorXlRef_const pt_to_target_ptr, const std::vector<py::array_t<long>> &target_to_pt_ptr, const Ref<const VectorXl> &target_seq) {
  int num_targets = target_seq.size();
  int num_nodes = gtsp_cost_mat.rows();
  if (tour.size() != num_targets + 2) {
    throw std::runtime_error("Tour not correct size");
  }
  tour(0) = 0;
  tour(num_targets + 1) = 0;
  VectorXd g_vals = std::numeric_limits<double>::infinity()*VectorXd::Ones(num_nodes);
  g_vals(0) = 0;
  VectorXi backpointers = VectorXi::Zero(num_nodes);
  int target_idx = num_targets;
  for (int tour_idx = 1; tour_idx < num_targets + 1; ++tour_idx) {
    int next_target_idx = target_seq(tour_idx - 1);
    auto ptr = target_to_pt_ptr[target_idx].unchecked<1>();
    auto next_ptr = target_to_pt_ptr[next_target_idx].unchecked<1>();
    for (int ptr_idx = 0; ptr_idx < ptr.size(); ++ptr_idx) {
      int node_idx = ptr(ptr_idx);
      for (int next_ptr_idx = 0; next_ptr_idx < next_ptr.size(); ++next_ptr_idx) {
        int next_node_idx = next_ptr(next_ptr_idx);
        double g_cand = g_vals(node_idx) + gtsp_cost_mat(node_idx, next_node_idx);
        /*
        if (std::isfinite(g_cand)) {
          longest_feasible_prefix(0) = tour_idx;
        }
        */
        if (g_cand < g_vals(next_node_idx)) {
          g_vals(next_node_idx) = g_cand;
          backpointers(next_node_idx) = node_idx;
        }
      }
    }
    target_idx = next_target_idx;
  }

  target_idx = target_seq(num_targets - 1);
  auto ptr = target_to_pt_ptr[target_idx].unchecked<1>();
  int best_final_node_idx = 0;
  double best_g_val = std::numeric_limits<double>::infinity();
  for (int ptr_idx = 0; ptr_idx < ptr.size(); ++ptr_idx) {
    int node_idx = ptr(ptr_idx);
    if (g_vals(node_idx) < best_g_val) {
      best_final_node_idx = node_idx;
      best_g_val = g_vals(node_idx);
    }
  }

  if (std::isinf(best_g_val)) {
    return false;
  }

  int node_idx = best_final_node_idx;
  for (int tour_idx = num_targets; tour_idx > 0; --tour_idx) {
    tour(tour_idx) = node_idx;
    node_idx = backpointers(node_idx);
  }

  return true;
}

double reopt_gtsp_tour_with_list_of_cost_matrices(Ref<VectorXl> tour, const std::vector<py::array_t<double>> &flattened_cost_mat_list, const Ref<const VectorXl> &num_nodes_per_tour_idx) {
  if (tour.size() - 1 != flattened_cost_mat_list.size()) {
    throw std::runtime_error("Number of cost matrices should be tour length - 1");
  }

  if (tour.size() != num_nodes_per_tour_idx.size()) {
    throw std::runtime_error("Tour size should match num_nods_per_tour_idx size");
  }

  tour(0) = 0;
  VectorXd prev_g_vals = VectorXd::Zero(1);
  std::vector<VectorXi> backpointers_per_layer;
  for (int tour_idx = 1; tour_idx < tour.size(); ++tour_idx) {
    int num_nodes1 = num_nodes_per_tour_idx(tour_idx - 1);
    int num_nodes2 = num_nodes_per_tour_idx(tour_idx);
    VectorXd g_vals = std::numeric_limits<double>::infinity()*VectorXd::Ones(num_nodes2);
    backpointers_per_layer.push_back(VectorXi::Zero(num_nodes2));

    auto cost_mat = flattened_cost_mat_list[tour_idx - 1].unchecked<1>();
    for (int node_idx1 = 0; node_idx1 < num_nodes1; ++node_idx1) {
      for (int node_idx2 = 0; node_idx2 < num_nodes2; ++node_idx2) {
        double g_cand = prev_g_vals(node_idx1) + cost_mat(node_idx1*num_nodes2 + node_idx2);
        if (g_cand < g_vals(node_idx2)) {
          g_vals(node_idx2) = g_cand;
          backpointers_per_layer[tour_idx - 1](node_idx2) = node_idx1;
        }
      }
    }
    prev_g_vals = g_vals;
  }

  int best_final_node_idx = 0;
  double best_g_val = std::numeric_limits<double>::infinity();
  for (int node_idx2 = 0; node_idx2 < num_nodes_per_tour_idx(tour.size() - 1); ++node_idx2) {
    if (prev_g_vals(node_idx2) < best_g_val) {
      best_final_node_idx = node_idx2;
      best_g_val = prev_g_vals(node_idx2);
    }
  }

  if (std::isinf(best_g_val)) {
    return best_g_val;
  }

  int node_idx = best_final_node_idx;
  for (int tour_idx = tour.size() - 1; tour_idx > 0; --tour_idx) {
    tour(tour_idx) = node_idx;
    node_idx = backpointers_per_layer[tour_idx - 1](node_idx);
  }

  return best_g_val;
}
