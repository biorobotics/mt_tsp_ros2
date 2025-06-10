#pragma once
#include <memory>
#include <unordered_set>
#include <set>
#include <vector>
#include <chrono>
#include <Eigen/Dense>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "mt_tsp_ros2/dag_dfs.h"
#include <random>
#include <tuple>

typedef const Ref<const Matrix<long, Dynamic, Dynamic, RowMajor>> &RowMatrixXlRef_const;

class LifelongDAGDFSPlanner {
  public:
    LifelongDAGDFSPlanner(RowMatrixXdRef_const gtsp_cost_mat, VectorXlRef_const pt_to_target_ptr, const std::vector<py::array_t<long>> &target_to_pt_ptr, bool do_prune, bool do_sort, RowMatrixXdRef_const all_pts) {
      before_time = 0.;
      if (do_prune) {
        auto timer_start = std::chrono::high_resolution_clock::now();
        int num_nodes = gtsp_cost_mat.rows();
        int num_targets = target_to_pt_ptr.size() - 1; // -1 because we have a dummy target associated with the depot
        before = MatrixXb::Ones(num_nodes, num_targets);
        for (int target_idx = 0; target_idx < num_targets; ++target_idx) {
          auto ptr = target_to_pt_ptr[target_idx].unchecked<1>();
          for (int node_idx = 0; node_idx < num_nodes; ++node_idx) {
            for (int ptr_idx = 0; ptr_idx < ptr.size(); ++ptr_idx) {
              int node_idx2 = ptr(ptr_idx);
              if (std::isfinite(gtsp_cost_mat(node_idx, node_idx2))) {
                before(node_idx, target_idx) = false;
                break;
              }
            }
          }
        }

        for (int node_idx = 1; node_idx < num_nodes; ++node_idx) {
          before(node_idx, pt_to_target_ptr(node_idx)) = false;
        }
        auto timer_stop = std::chrono::high_resolution_clock::now();
        auto nanos = std::chrono::duration_cast<std::chrono::microseconds>(timer_stop - timer_start).count();
        before_time = ((double)nanos)/1e9;
      }
    }

    VectorXd plan_biased(Ref<VectorXl> tour, RowMatrixXdRef_const gtsp_cost_mat, VectorXlRef_const pt_to_target_ptr, const std::vector<py::array_t<long>> &target_to_pt_ptr, bool do_prune, bool do_sort, double time_limit, py::object other_tour_queue, RowMatrixXdRef_const all_pts, RowMatrixXlRef_const deleted_edges, const Ref<const VectorXl> &bias_tour) {
      auto timer_start = std::chrono::high_resolution_clock::now();
      VectorXd profiling_data = VectorXd::Zero(4);
      int num_targets = target_to_pt_ptr.size() - 1; // -1 because we have a dummy target associated with the depot

      std::mt19937 rng;

      if (do_prune) {
        // Update the before sets
        auto tmp_timer_start = std::chrono::high_resolution_clock::now();
        for (int edge_idx = 0; edge_idx < deleted_edges.rows(); ++edge_idx) {
          int node_idx1 = deleted_edges(edge_idx, 0);
          int node_idx2 = deleted_edges(edge_idx, 1);

          int target_idx = pt_to_target_ptr(node_idx2);

          auto ptr = target_to_pt_ptr[target_idx].unchecked<1>();
          before(node_idx1, target_idx) = true; // Assume we can't get from node_idx1 to target_idx anymore...
          for (int ptr_idx = 0; ptr_idx < ptr.size(); ++ptr_idx) {
            int node_idx2 = ptr(ptr_idx);
            if (std::isfinite(gtsp_cost_mat(node_idx1, node_idx2))) {
              // ... unless we can still get to target_idx via some other node
              before(node_idx1, target_idx) = false;
              break;
            }
          }
        }
        auto tmp_timer_stop = std::chrono::high_resolution_clock::now();
        auto tmp_nanos = std::chrono::duration_cast<std::chrono::microseconds>(tmp_timer_stop - tmp_timer_start).count();
        profiling_data(0) = ((double)tmp_nanos)/1e9;
        before_time += profiling_data(0);
      }

      std::vector<std::tuple<VectorXb, int, int>> bias_edges;
      VectorXb visited_targets(num_targets);
      for (int i = 0; i < bias_tour.size() - 1; ++i) {
        bias_edges.push_back(std::tuple<VectorXb, int, int>(visited_targets, bias_tour(i), bias_tour(i + 1)));
        if (bias_tour(i) != 0) {
          visited_targets(pt_to_target_ptr(bias_tour(i))) = true;
        }
      }

      std::unordered_set<VectorXi, key_hash> closed_list;
      std::vector<DFSNodePtr> stack;
      stack.push_back(std::make_shared<DFSNode>(nullptr, 0, -1, num_targets));

      int num_nodes = gtsp_cost_mat.rows();

      while (stack.size()) {
        auto timer_stop = std::chrono::high_resolution_clock::now();
        auto nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
        if (((double)nanos)/1e9 > time_limit) {
          break;
        }

        if (other_tour_queue.attr("qsize")().cast<long>()) {
          break;
        }

        DFSNodePtr pop = stack.back();
        stack.pop_back();

        if (closed_list.find(pop->key) != closed_list.end()) {
          continue;
        }

        closed_list.insert(pop->key);

        std::vector<int> neighbors;
        std::vector<double> neighbor_times;
        int num_visited_targets = pop->visited_targets.cast<int>().sum();
        if (num_visited_targets == num_targets) {
          std::vector<long> tour_vec;
          tour_vec.push_back(pop->final_pt_idx);
          DFSNodePtr node = pop->parent;
          while (node != nullptr) {
            tour_vec.push_back(node->final_pt_idx);
            node = node->parent;
          }
          std::reverse(tour_vec.begin(), tour_vec.end());
          // Assume open tsp
          tour_vec.push_back(0);
          tour = Map<VectorXl>(tour_vec.data(), tour_vec.size());
          return profiling_data;
        } else {
          for (int target_idx = 0; target_idx < num_targets; ++target_idx) {
            if (pop->visited_targets(target_idx)) {
              continue;
            }
            auto ptr = target_to_pt_ptr[target_idx].unchecked<1>();
            for (int ptr_idx = 0; ptr_idx < ptr.size(); ++ptr_idx) {
              int node_idx = ptr(ptr_idx);
              if (std::isfinite(gtsp_cost_mat(pop->final_pt_idx, node_idx))) {
                neighbors.push_back(node_idx);
                neighbor_times.push_back(all_pts(node_idx, 0));
              }
            }
          }

          if (neighbors.size() == 0) {
            continue;
          }
        }

        std::vector<size_t> sort_idx;
        if (do_sort) {
          // tmp_timer_start = std::chrono::high_resolution_clock::now();
          sort_idx = sort_indexes(neighbor_times);
          // tmp_timer_stop = std::chrono::high_resolution_clock::now();
          // tmp_nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(tmp_timer_stop - tmp_timer_start).count();
          // profiling_data(1) += ((double)tmp_nanos)/1e9; // sort time
          std::reverse(sort_idx.begin(), sort_idx.end());
        } else {
          sort_idx.resize(neighbor_times.size());
          for (int node_idx = 0; node_idx < neighbor_times.size(); ++node_idx) {
            sort_idx[node_idx] = node_idx;
          }
          std::shuffle(sort_idx.begin(), sort_idx.end(), rng);
          // throw std::runtime_error("Random successor ordering not implemented");
        }

        if (bias_tour.size() && std::get<1>(bias_edges[num_visited_targets]) == pop->final_pt_idx) {
          bool subset_same = true;
          for (int target_idx = 0; target_idx < num_targets; ++target_idx) {
            bool subset_same = true;
            if (pop->visited_targets(target_idx) != std::get<0>(bias_edges[num_visited_targets])(target_idx)) {
              subset_same = false;
              break;
            }
          }
          if (subset_same) {
            for (int i = 0; i < sort_idx.size(); ++i) {
              int j = sort_idx[i];
              if (std::get<2>(bias_edges[num_visited_targets]) == neighbors[j]) {
                sort_idx.erase(sort_idx.begin() + i);
                sort_idx.push_back(j);
                break;
              }
            }
          }
        }

        for (int neighbor_idx : sort_idx) {
          int pt_idx = neighbors[neighbor_idx];
          if (do_prune) {
            // tmp_timer_start = std::chrono::high_resolution_clock::now();
            bool prune = false;
            for (int target_idx = 0; target_idx < num_targets; ++target_idx) {
              if (target_idx != pt_to_target_ptr(pt_idx) && before(pt_idx, target_idx) && !pop->visited_targets(target_idx)) {
                prune = true;
                break;
              }
            }

            // tmp_timer_stop = std::chrono::high_resolution_clock::now();
            // tmp_nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(tmp_timer_stop - tmp_timer_start).count();
            // profiling_data(2) += ((double)tmp_nanos)/1e9; // prune time

            if (prune) {
              continue;
            }
          }

          // tmp_timer_start = std::chrono::high_resolution_clock::now();

          DFSNodePtr neighbor_node = std::make_shared<DFSNode>(pop, pt_idx, pt_to_target_ptr(pt_idx));

          // tmp_timer_stop = std::chrono::high_resolution_clock::now();
          // tmp_nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(tmp_timer_stop - tmp_timer_start).count();
          // profiling_data(3) += ((double)tmp_nanos)/1e9; // node gen time

          if (closed_list.find(neighbor_node->key) != closed_list.end()) {
            continue;
          }

          stack.push_back(neighbor_node);
        }
      }
      tour(0) = -1;
      return profiling_data;
    }

    VectorXd plan(Ref<VectorXl> tour, RowMatrixXdRef_const gtsp_cost_mat, VectorXlRef_const pt_to_target_ptr, const std::vector<py::array_t<long>> &target_to_pt_ptr, bool do_prune, bool do_sort, double time_limit, py::object other_tour_queue, RowMatrixXdRef_const all_pts, RowMatrixXlRef_const deleted_edges) {
      VectorXl dummy_bias_tour(0);
      return plan_biased(tour, gtsp_cost_mat, pt_to_target_ptr, target_to_pt_ptr, do_prune, do_sort, time_limit, other_tour_queue, all_pts, deleted_edges, dummy_bias_tour);
    }

    double get_before_time() {
      return before_time;
    }
  private:
    MatrixXb before;
    double before_time;
};

class LifelongDAGDFSPlannerUnorderedSet {
  public:
    LifelongDAGDFSPlannerUnorderedSet(RowMatrixXdRef_const gtsp_cost_mat, VectorXlRef_const pt_to_target_ptr, const std::vector<py::array_t<long>> &target_to_pt_ptr, bool do_prune, bool do_sort, RowMatrixXdRef_const all_pts) {
      before_time = 0.;
      if (do_prune) {
        auto timer_start = std::chrono::high_resolution_clock::now();
        int num_nodes = gtsp_cost_mat.rows();
        int num_targets = target_to_pt_ptr.size() - 1; // -1 because we have a dummy target associated with the depot
        before = MatrixXb::Ones(num_nodes, num_targets);
        for (int target_idx = 0; target_idx < num_targets; ++target_idx) {
          auto ptr = target_to_pt_ptr[target_idx].unchecked<1>();
          for (int node_idx = 0; node_idx < num_nodes; ++node_idx) {
            for (int ptr_idx = 0; ptr_idx < ptr.size(); ++ptr_idx) {
              int node_idx2 = ptr(ptr_idx);
              if (std::isfinite(gtsp_cost_mat(node_idx, node_idx2))) {
                before(node_idx, target_idx) = false;
                break;
              }
            }
          }
        }

        for (int node_idx = 1; node_idx < num_nodes; ++node_idx) {
          before(node_idx, pt_to_target_ptr(node_idx)) = false;
        }
        auto timer_stop = std::chrono::high_resolution_clock::now();
        auto nanos = std::chrono::duration_cast<std::chrono::microseconds>(timer_stop - timer_start).count();
        before_time = ((double)nanos)/1e9;
      }
    }

    VectorXd plan(Ref<VectorXl> tour, RowMatrixXdRef_const gtsp_cost_mat, VectorXlRef_const pt_to_target_ptr, const std::vector<py::array_t<long>> &target_to_pt_ptr, bool do_prune, bool do_sort, double time_limit, py::object other_tour_queue, RowMatrixXdRef_const all_pts, RowMatrixXlRef_const deleted_edges) {
      auto timer_start = std::chrono::high_resolution_clock::now();
      VectorXd profiling_data = VectorXd::Zero(4);
      int num_targets = target_to_pt_ptr.size() - 1; // -1 because we have a dummy target associated with the depot

      if (do_prune) {
        // Update the before sets
        auto tmp_timer_start = std::chrono::high_resolution_clock::now();
        for (int edge_idx = 0; edge_idx < deleted_edges.rows(); ++edge_idx) {
          int node_idx1 = deleted_edges(edge_idx, 0);
          int node_idx2 = deleted_edges(edge_idx, 1);

          int target_idx = pt_to_target_ptr(node_idx2);

          auto ptr = target_to_pt_ptr[target_idx].unchecked<1>();
          before(node_idx1, target_idx) = true; // Assume we can't get from node_idx1 to target_idx anymore...
          for (int ptr_idx = 0; ptr_idx < ptr.size(); ++ptr_idx) {
            int node_idx2 = ptr(ptr_idx);
            if (std::isfinite(gtsp_cost_mat(node_idx1, node_idx2))) {
              // ... unless we can still get to target_idx via some other node
              before(node_idx1, target_idx) = false;
              break;
            }
          }
        }
        auto tmp_timer_stop = std::chrono::high_resolution_clock::now();
        auto tmp_nanos = std::chrono::duration_cast<std::chrono::microseconds>(tmp_timer_stop - tmp_timer_start).count();
        profiling_data(0) = ((double)tmp_nanos)/1e9;
        before_time += profiling_data(0);
      }

      std::unordered_set<VectorXi, key_hash> closed_list;
      std::vector<DFSNodeUnorderedSetPtr> stack;
      stack.push_back(std::make_shared<DFSNodeUnorderedSet>(nullptr, 0, -1));

      int num_nodes = gtsp_cost_mat.rows();

      while (stack.size()) {
        auto timer_stop = std::chrono::high_resolution_clock::now();
        auto nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
        if (((double)nanos)/1e9 > time_limit) {
          break;
        }

        if (other_tour_queue.attr("qsize")().cast<long>()) {
          break;
        }

        DFSNodeUnorderedSetPtr pop = stack.back();
        stack.pop_back();

        if (closed_list.find(pop->key) != closed_list.end()) {
          continue;
        }

        closed_list.insert(pop->key);

        std::vector<int> neighbors;
        std::vector<double> neighbor_times;
        if (pop->visited_targets.size() == num_targets) {
          std::vector<long> tour_vec;
          tour_vec.push_back(pop->final_pt_idx);
          DFSNodeUnorderedSetPtr node = pop->parent;
          while (node != nullptr) {
            tour_vec.push_back(node->final_pt_idx);
            node = node->parent;
          }
          std::reverse(tour_vec.begin(), tour_vec.end());
          // Assume open tsp
          tour_vec.push_back(0);
          tour = Map<VectorXl>(tour_vec.data(), tour_vec.size());
          return profiling_data;
        } else {
          for (int target_idx = 0; target_idx < num_targets; ++target_idx) {
            if (pop->visited_targets.find(target_idx) != pop->visited_targets.end()) {
              continue;
            }
            auto ptr = target_to_pt_ptr[target_idx].unchecked<1>();
            for (int ptr_idx = 0; ptr_idx < ptr.size(); ++ptr_idx) {
              int node_idx = ptr(ptr_idx);
              if (std::isfinite(gtsp_cost_mat(pop->final_pt_idx, node_idx))) {
                neighbors.push_back(node_idx);
                neighbor_times.push_back(all_pts(node_idx, 0));
              }
            }
          }

          if (neighbors.size() == 0) {
            continue;
          }
        }

        std::vector<size_t> sort_idx;
        if (do_sort) {
          // tmp_timer_start = std::chrono::high_resolution_clock::now();
          sort_idx = sort_indexes(neighbor_times);
          // tmp_timer_stop = std::chrono::high_resolution_clock::now();
          // tmp_nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(tmp_timer_stop - tmp_timer_start).count();
          // profiling_data(1) += ((double)tmp_nanos)/1e9; // sort time
          std::reverse(sort_idx.begin(), sort_idx.end());
        } else {
          throw std::runtime_error("Random successor ordering not implemented");
        }
        for (int neighbor_idx : sort_idx) {
          int pt_idx = neighbors[neighbor_idx];
          if (do_prune) {
            // tmp_timer_start = std::chrono::high_resolution_clock::now();
            bool prune = false;
            for (int target_idx = 0; target_idx < num_targets; ++target_idx) {
              if (target_idx != pt_to_target_ptr(pt_idx) && before(pt_idx, target_idx) && pop->visited_targets.find(target_idx) == pop->visited_targets.end()) {
                prune = true;
                break;
              }
            }

            // tmp_timer_stop = std::chrono::high_resolution_clock::now();
            // tmp_nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(tmp_timer_stop - tmp_timer_start).count();
            // profiling_data(2) += ((double)tmp_nanos)/1e9; // prune time

            if (prune) {
              continue;
            }
          }

          // tmp_timer_start = std::chrono::high_resolution_clock::now();

          DFSNodeUnorderedSetPtr neighbor_node = std::make_shared<DFSNodeUnorderedSet>(pop, pt_idx, pt_to_target_ptr(pt_idx));

          // tmp_timer_stop = std::chrono::high_resolution_clock::now();
          // tmp_nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(tmp_timer_stop - tmp_timer_start).count();
          // profiling_data(3) += ((double)tmp_nanos)/1e9; // node gen time

          if (closed_list.find(neighbor_node->key) != closed_list.end()) {
            continue;
          }

          stack.push_back(neighbor_node);
        }
      }
      tour(0) = -1;
      return profiling_data;
    }

    double get_before_time() {
      return before_time;
    }
  private:
    MatrixXb before;
    double before_time;
};
