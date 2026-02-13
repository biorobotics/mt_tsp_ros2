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
#include <omp.h>

typedef const Ref<const Matrix<long, Dynamic, Dynamic, RowMajor>> &RowMatrixXlRef_const;
typedef const Ref<const Matrix<bool, Dynamic, Dynamic, RowMajor>> &RowMatrixXbRef_const;

class LifelongDAGDFSPlanner {
  public:
    LifelongDAGDFSPlanner(RowMatrixXdRef_const gtsp_cost_mat, VectorXlRef_const pt_to_target_ptr, const std::vector<py::array_t<long>> &target_to_pt_ptr, bool do_prune, bool do_sort, RowMatrixXdRef_const all_pts, int num_threads) : num_threads(num_threads) {
      before_time = 0.;
      num_targets = target_to_pt_ptr.size() - 1; // -1 because we have a dummy target associated with the depot
      if (do_prune) {
        auto timer_start = std::chrono::high_resolution_clock::now();
        int num_nodes = gtsp_cost_mat.rows();
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

      max_observed_stack_size = 0;
      max_observed_closed_list_size = 0;
    }

    MatrixXb compute_before_target_to_target(const std::vector<py::array_t<long>> &target_to_pt_ptr) {
      before_target_to_target = MatrixXb::Ones(num_targets, num_targets);
      for (int target_idx1 = 0; target_idx1 < num_targets; ++target_idx1) {
        for (int target_idx2 = 0; target_idx2 < num_targets; ++target_idx2) {
          if (target_idx1 == target_idx2) {
            continue;
          }
          auto ptr = target_to_pt_ptr[target_idx1].unchecked<1>();
          for (int ptr_idx = 0; ptr_idx < ptr.size(); ++ptr_idx) {
            int node_idx1 = ptr(ptr_idx);
            if (!before(node_idx1, target_idx2)) {
              before_target_to_target(target_idx1, target_idx2) = false;
              break;
            }
          }
        }
      }
      return before_target_to_target;
    }

    void adjust_before_using_tour(RowMatrixXdRef_const gtsp_cost_mat, VectorXlRef_const pt_to_target_ptr, VectorXlRef_const tour) {
      for (int tour_idx = 1; tour_idx < tour.size() - 1; ++tour_idx) {
        int node_idx = tour(tour_idx);
        for (int next_tour_idx = tour_idx + 1; next_tour_idx < tour.size() - 1; ++next_tour_idx) {
          int next_node_idx = tour(next_tour_idx);
          before(node_idx, pt_to_target_ptr(next_node_idx)) = false;
        }
      }
    }

    VectorXd plan_considering_evaluations(Ref<VectorXl> tour, RowMatrixXdRef_const gtsp_cost_mat, VectorXlRef_const pt_to_target_ptr, const std::vector<py::array_t<long>> &target_to_pt_ptr, bool do_prune, bool do_sort, double time_limit, py::object other_tour_queue, RowMatrixXdRef_const all_pts, RowMatrixXlRef_const deleted_edges, const Ref<const VectorXl> &bias_tour, bool sort_by_time, RowMatrixXbRef_const unevaluated_edges) {
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
      VectorXb visited_targets = VectorXb::Zero(num_targets);
      for (int i = 0; i < bias_tour.size() - 1; ++i) {
        bias_edges.push_back(std::tuple<VectorXb, int, int>(visited_targets, bias_tour(i), bias_tour(i + 1)));
        if (bias_tour(i + 1) != 0) {
          visited_targets(pt_to_target_ptr(bias_tour(i + 1))) = true;
        }
      }

      std::unordered_set<VectorXi, key_hash> closed_list;
      std::vector<DFSNodePtr> stack;
      stack.push_back(std::make_shared<DFSNode>(nullptr, 0, -1, num_targets));

      int num_nodes = gtsp_cost_mat.rows();

      omp_set_num_threads(num_threads);

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

        int num_visited_targets = pop->visited_targets.cast<int>().sum();
        /*
        if (num_visited_targets != num_targets && bias_tour.size() &&
            std::isfinite(gtsp_cost_mat(pop->final_pt_idx, std::get<2>(bias_edges[num_visited_targets])))) {
          bool subset_same = true;
          for (int target_idx = 0; target_idx < num_targets; ++target_idx) {
            if (pop->visited_targets(target_idx) != std::get<0>(bias_edges[num_visited_targets])(target_idx)) {
              subset_same = false;
              break;
            }
          }
          if (subset_same) {
            int pt_idx = std::get<2>(bias_edges[num_visited_targets]);
            DFSNodePtr neighbor_node = std::make_shared<DFSNode>(pop, pt_idx, pt_to_target_ptr(pt_idx));

            if (closed_list.find(neighbor_node->key) == closed_list.end()) {
              stack.push_back(pop);
              stack.push_back(neighbor_node);
              continue;
            }
          }
        }
        */
        bool pushed = false;
        for (auto edge : bias_edges) {
          int pt_idx1 = std::get<1>(edge);
          int pt_idx2 = std::get<2>(edge);
          if (pt_idx2 == 0) {
            continue;
          }
          if (pop->final_pt_idx == pt_idx1 && std::isfinite(gtsp_cost_mat(pop->final_pt_idx, pt_idx2)) && !pop->visited_targets(pt_to_target_ptr(pt_idx2))) {
            DFSNodePtr neighbor_node = std::make_shared<DFSNode>(pop, pt_idx2, pt_to_target_ptr(pt_idx2));
            if (closed_list.find(neighbor_node->key) == closed_list.end()) {
              stack.push_back(pop);
              stack.push_back(neighbor_node);
              pushed = true;
              break;
            }
          }
        }
        if (pushed) {
          continue;
        }

        closed_list.insert(pop->key);

        std::vector<std::vector<DFSNodePtr>> neighbors_per_thread(num_threads);
        // First pair element is whether the edge is unevaluated, second pair is either cost or time depending on whether sort_by_time is true
        std::vector<std::vector<std::pair<bool, double>>> neighbor_sort_vals_per_thread(num_threads);
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
          #pragma omp parallel for
          for (int node_idx = 1; node_idx < num_nodes; ++node_idx) {
            int target_idx = pt_to_target_ptr(node_idx);
            if (pop->visited_targets(target_idx)) {
              continue;
            }
            if (std::isfinite(gtsp_cost_mat(pop->final_pt_idx, node_idx))) {
              if (do_prune) {
                // tmp_timer_start = std::chrono::high_resolution_clock::now();
                bool prune = false;
                for (int target_idx2 = 0; target_idx2 < num_targets; ++target_idx2) {
                  if (target_idx2 != target_idx && before(node_idx, target_idx2) && !pop->visited_targets(target_idx2)) {
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

              DFSNodePtr neighbor_node = std::make_shared<DFSNode>(pop, node_idx, target_idx);

              // tmp_timer_stop = std::chrono::high_resolution_clock::now();
              // tmp_nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(tmp_timer_stop - tmp_timer_start).count();
              // profiling_data(3) += ((double)tmp_nanos)/1e9; // node gen time

              if (closed_list.find(neighbor_node->key) != closed_list.end()) {
                continue;
              }

              neighbors_per_thread[omp_get_thread_num()].push_back(neighbor_node);
              if (sort_by_time) {
                neighbor_sort_vals_per_thread[omp_get_thread_num()].push_back(std::pair<bool, double>(unevaluated_edges(pop->final_pt_idx, node_idx), all_pts(node_idx, 0)));
              } else {
                // Sort by cost
                neighbor_sort_vals_per_thread[omp_get_thread_num()].push_back(std::pair<bool, double>(unevaluated_edges(pop->final_pt_idx, node_idx), gtsp_cost_mat(pop->final_pt_idx, node_idx)));
              }
            }
          }
        }

        int num_neighbors = 0;
        for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
          num_neighbors += neighbors_per_thread[thread_idx].size();
        }
        std::vector<std::pair<bool, double>> neighbor_sort_vals(num_neighbors);
        std::vector<DFSNodePtr> neighbors(num_neighbors);

        int i = 0;
        for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
          if (neighbor_sort_vals_per_thread[thread_idx].size() != neighbors_per_thread[thread_idx].size()) {
            throw std::runtime_error("Size error");
          }
          for (int neighbor_idx = 0; neighbor_idx < neighbor_sort_vals_per_thread[thread_idx].size(); ++neighbor_idx) {
            if (i >= num_neighbors) {
              throw std::runtime_error("Index error");
            }
            neighbor_sort_vals[i] = neighbor_sort_vals_per_thread[thread_idx][neighbor_idx];
            neighbors[i] = neighbors_per_thread[thread_idx][neighbor_idx];
            ++i;
          }
        }

        std::vector<size_t> sort_idx;
        if (do_sort) {
          // tmp_timer_start = std::chrono::high_resolution_clock::now();
          sort_idx = sort_indexes(neighbor_sort_vals);
          // tmp_timer_stop = std::chrono::high_resolution_clock::now();
          // tmp_nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(tmp_timer_stop - tmp_timer_start).count();
          // profiling_data(1) += ((double)tmp_nanos)/1e9; // sort time
          std::reverse(sort_idx.begin(), sort_idx.end());
        } else {
          sort_idx.resize(neighbor_sort_vals.size());
          for (int node_idx = 0; node_idx < neighbor_sort_vals.size(); ++node_idx) {
            sort_idx[node_idx] = node_idx;
          }
          std::shuffle(sort_idx.begin(), sort_idx.end(), rng);
          // throw std::runtime_error("Random successor ordering not implemented");
        }

        for (int neighbor_idx : sort_idx) {
          stack.push_back(neighbors[neighbor_idx]);
        }
      }
      tour(0) = -1;
      return profiling_data;
    }

    VectorXd plan_biased(Ref<VectorXl> tour, RowMatrixXdRef_const gtsp_cost_mat, VectorXlRef_const pt_to_target_ptr, const std::vector<py::array_t<long>> &target_to_pt_ptr, bool do_prune, bool do_sort, double time_limit, py::object other_tour_queue, RowMatrixXdRef_const all_pts, RowMatrixXlRef_const deleted_edges, const Ref<const VectorXl> &bias_tour, bool sort_by_time) {
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
      VectorXb visited_targets = VectorXb::Zero(num_targets);
      for (int i = 0; i < bias_tour.size() - 1; ++i) {
        bias_edges.push_back(std::tuple<VectorXb, int, int>(visited_targets, bias_tour(i), bias_tour(i + 1)));
        if (bias_tour(i + 1) != 0) {
          visited_targets(pt_to_target_ptr(bias_tour(i + 1))) = true;
        }
      }

      std::unordered_set<VectorXi, key_hash> closed_list;
      std::vector<DFSNodePtr> stack;
      stack.push_back(std::make_shared<DFSNode>(nullptr, 0, -1, num_targets));

      int num_nodes = gtsp_cost_mat.rows();

      omp_set_num_threads(num_threads);

      while (stack.size()) {
        max_observed_stack_size = std::max(max_observed_stack_size, (int)(stack.size()));
        max_observed_closed_list_size = std::max(max_observed_closed_list_size, (int)(closed_list.size()));
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

        int num_visited_targets = pop->visited_targets.cast<int>().sum();
        if (num_visited_targets != num_targets && bias_tour.size() &&
            std::isfinite(gtsp_cost_mat(pop->final_pt_idx, std::get<2>(bias_edges[num_visited_targets])))) {
          bool subset_same = true;
          for (int target_idx = 0; target_idx < num_targets; ++target_idx) {
            if (pop->visited_targets(target_idx) != std::get<0>(bias_edges[num_visited_targets])(target_idx)) {
              subset_same = false;
              break;
            }
          }
          if (subset_same) {
            int pt_idx = std::get<2>(bias_edges[num_visited_targets]);
            DFSNodePtr neighbor_node = std::make_shared<DFSNode>(pop, pt_idx, pt_to_target_ptr(pt_idx));

            if (closed_list.find(neighbor_node->key) == closed_list.end()) {
              stack.push_back(pop);
              stack.push_back(neighbor_node);
              continue;
            }
          }
        }

        closed_list.insert(pop->key);

        std::vector<std::vector<DFSNodePtr>> neighbors_per_thread(num_threads);
        std::vector<std::vector<double>> neighbor_sort_vals_per_thread(num_threads);
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
          #pragma omp parallel for
          for (int node_idx = 1; node_idx < num_nodes; ++node_idx) {
            int target_idx = pt_to_target_ptr(node_idx);
            if (pop->visited_targets(target_idx)) {
              continue;
            }
            if (std::isfinite(gtsp_cost_mat(pop->final_pt_idx, node_idx))) {
              if (do_prune) {
                // tmp_timer_start = std::chrono::high_resolution_clock::now();
                bool prune = false;
                for (int target_idx2 = 0; target_idx2 < num_targets; ++target_idx2) {
                  if (target_idx2 != target_idx && before(node_idx, target_idx2) && !pop->visited_targets(target_idx2)) {
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

              DFSNodePtr neighbor_node = std::make_shared<DFSNode>(pop, node_idx, target_idx);

              // tmp_timer_stop = std::chrono::high_resolution_clock::now();
              // tmp_nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(tmp_timer_stop - tmp_timer_start).count();
              // profiling_data(3) += ((double)tmp_nanos)/1e9; // node gen time

              if (closed_list.find(neighbor_node->key) != closed_list.end()) {
                continue;
              }

              neighbors_per_thread[omp_get_thread_num()].push_back(neighbor_node);
              if (sort_by_time) {
                neighbor_sort_vals_per_thread[omp_get_thread_num()].push_back(all_pts(node_idx, 0));
              } else {
                // Sort by cost
                neighbor_sort_vals_per_thread[omp_get_thread_num()].push_back(gtsp_cost_mat(pop->final_pt_idx, node_idx));
              }
            }
          }
        }

        int num_neighbors = 0;
        for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
          num_neighbors += neighbors_per_thread[thread_idx].size();
        }
        std::vector<double> neighbor_sort_vals(num_neighbors);
        std::vector<DFSNodePtr> neighbors(num_neighbors);

        int i = 0;
        for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
          if (neighbor_sort_vals_per_thread[thread_idx].size() != neighbors_per_thread[thread_idx].size()) {
            throw std::runtime_error("Size error");
          }
          for (int neighbor_idx = 0; neighbor_idx < neighbor_sort_vals_per_thread[thread_idx].size(); ++neighbor_idx) {
            if (i >= num_neighbors) {
              throw std::runtime_error("Index error");
            }
            neighbor_sort_vals[i] = neighbor_sort_vals_per_thread[thread_idx][neighbor_idx];
            neighbors[i] = neighbors_per_thread[thread_idx][neighbor_idx];
            ++i;
          }
        }

        std::vector<size_t> sort_idx;
        if (do_sort) {
          // tmp_timer_start = std::chrono::high_resolution_clock::now();
          sort_idx = sort_indexes(neighbor_sort_vals);
          // tmp_timer_stop = std::chrono::high_resolution_clock::now();
          // tmp_nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(tmp_timer_stop - tmp_timer_start).count();
          // profiling_data(1) += ((double)tmp_nanos)/1e9; // sort time
          std::reverse(sort_idx.begin(), sort_idx.end());
        } else {
          sort_idx.resize(neighbor_sort_vals.size());
          for (int node_idx = 0; node_idx < neighbor_sort_vals.size(); ++node_idx) {
            sort_idx[node_idx] = node_idx;
          }
          std::shuffle(sort_idx.begin(), sort_idx.end(), rng);
          // throw std::runtime_error("Random successor ordering not implemented");
        }

        for (int neighbor_idx : sort_idx) {
          stack.push_back(neighbors[neighbor_idx]);
        }
      }
      tour(0) = -1;
      return profiling_data;
    }

    VectorXd plan(Ref<VectorXl> tour, RowMatrixXdRef_const gtsp_cost_mat, VectorXlRef_const pt_to_target_ptr, const std::vector<py::array_t<long>> &target_to_pt_ptr, bool do_prune, bool do_sort, double time_limit, py::object other_tour_queue, RowMatrixXdRef_const all_pts, RowMatrixXlRef_const deleted_edges, bool sort_by_time) {
      // VectorXl dummy_bias_tour(0);
      // return plan_biased(tour, gtsp_cost_mat, pt_to_target_ptr, target_to_pt_ptr, do_prune, do_sort, time_limit, other_tour_queue, all_pts, deleted_edges, dummy_bias_tour, sort_by_time);

      VectorXl dummy_bias_tour(0);
      return plan_biased(tour, gtsp_cost_mat, pt_to_target_ptr, target_to_pt_ptr, do_prune, do_sort, time_limit, other_tour_queue, all_pts, deleted_edges, dummy_bias_tour, sort_by_time);
    }

    double get_before_time() {
      return before_time;
    }

    int get_max_observed_stack_size() {
      return max_observed_stack_size;
    }

    int get_max_observed_closed_list_size() {
      return max_observed_closed_list_size;
    }

  private:
    MatrixXb before;
    double before_time;
    int num_threads;
    MatrixXb before_target_to_target;
    int num_targets;
    int max_observed_stack_size;
    int max_observed_closed_list_size;
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
