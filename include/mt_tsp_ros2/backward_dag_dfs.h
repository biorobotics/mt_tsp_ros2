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

struct BackwardDFSNode { 
  BackwardDFSNode(std::shared_ptr<BackwardDFSNode> parent, 
                  int final_pt_idx,
                  int removed_target_idx,
                  int num_targets = 0) : parent(parent), 
                                         final_pt_idx(final_pt_idx) {
    if (parent != nullptr) {
      visited_targets = parent->visited_targets;
      if (removed_target_idx != -1) {
        visited_targets(removed_target_idx) = false;
      }
    } else {
      if (num_targets == 0) {
        throw std::runtime_error("num_targets can't be zero if parent == nullptr");
      }
      visited_targets = VectorXb::Ones(num_targets);
    }
    key = VectorXi::Zero(visited_targets.size() + 1);
    for (int i = 0; i < visited_targets.size(); ++i) {
      key(i) = visited_targets(i);
    }
    key(visited_targets.size()) = final_pt_idx;
  }

  std::shared_ptr<BackwardDFSNode> parent;
  VectorXb visited_targets;
  int final_pt_idx;
  VectorXi key;
};

typedef std::shared_ptr<BackwardDFSNode> BackwardDFSNodePtr;

class BackwardDAGDFSPlanner {
  public:
    BackwardDAGDFSPlanner(RowMatrixXdRef_const gtsp_cost_mat, VectorXlRef_const pt_to_target_ptr, const std::vector<py::array_t<long>> &target_to_pt_ptr, bool do_prune, bool do_sort, RowMatrixXdRef_const all_pts, int num_threads) : num_threads(num_threads) {
      num_targets = target_to_pt_ptr.size() - 1; // -1 because we have a dummy target associated with the depot
      if (do_prune) {
        auto timer_start = std::chrono::high_resolution_clock::now();
        int num_nodes = gtsp_cost_mat.rows();
        after = MatrixXb::Ones(num_nodes, num_targets);
        for (int target_idx = 0; target_idx < num_targets; ++target_idx) {
          auto ptr = target_to_pt_ptr[target_idx].unchecked<1>();
          for (int node_idx = 0; node_idx < num_nodes; ++node_idx) {
            for (int ptr_idx = 0; ptr_idx < ptr.size(); ++ptr_idx) {
              int node_idx1 = ptr(ptr_idx);
              if (std::isfinite(gtsp_cost_mat(node_idx1, node_idx))) {
                after(node_idx, target_idx) = false;
                break;
              }
            }
          }
        }

        for (int node_idx = 1; node_idx < num_nodes; ++node_idx) {
          after(node_idx, pt_to_target_ptr(node_idx)) = false;
        }
      }
    }

    VectorXd plan(Ref<VectorXl> tour, RowMatrixXdRef_const gtsp_cost_mat, VectorXlRef_const pt_to_target_ptr, const std::vector<py::array_t<long>> &target_to_pt_ptr, bool do_prune, bool do_sort, double time_limit, py::object other_tour_queue, RowMatrixXdRef_const all_pts, RowMatrixXlRef_const deleted_edges, const Ref<const VectorXl> &bias_tour) {
      auto timer_start = std::chrono::high_resolution_clock::now();
      VectorXd profiling_data = VectorXd::Zero(4);
      int num_targets = target_to_pt_ptr.size() - 1; // -1 because we have a dummy target associated with the depot

      std::mt19937 rng;

      if (do_prune) {
        // Update the after sets
        for (int edge_idx = 0; edge_idx < deleted_edges.rows(); ++edge_idx) {
          int node_idx1 = deleted_edges(edge_idx, 0);
          int node_idx2 = deleted_edges(edge_idx, 1);

          int target_idx = pt_to_target_ptr(node_idx1);

          auto ptr = target_to_pt_ptr[target_idx].unchecked<1>();
          after(node_idx2, target_idx) = true; // Assume we can't get from target_idx to node_idx2 anymore...
          for (int ptr_idx = 0; ptr_idx < ptr.size(); ++ptr_idx) {
            int node_idx1 = ptr(ptr_idx);
            if (std::isfinite(gtsp_cost_mat(node_idx1, node_idx2))) {
              // ... unless we can still get to node_idx2 via some other node
              after(node_idx2, target_idx) = false;
              break;
            }
          }
        }
      }

      std::vector<std::pair<VectorXb, int>> bias_dag_nodes;
      VectorXb visited_targets = VectorXb::Zero(num_targets);
      for (int i = 0; i < bias_tour.size() - 1; ++i) {
        bias_dag_nodes.push_back(std::pair<VectorXb, int>(visited_targets, bias_tour(i)));
        if (bias_tour(i + 1) != 0) {
          visited_targets(pt_to_target_ptr(bias_tour(i + 1))) = true;
        }
      }

      std::unordered_set<VectorXi, key_hash> closed_list;
      std::vector<BackwardDFSNodePtr> stack;
      stack.push_back(std::make_shared<BackwardDFSNode>(nullptr, 0, -1, num_targets));

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

        BackwardDFSNodePtr pop = stack.back();
        stack.pop_back();

        if (closed_list.find(pop->key) != closed_list.end()) {
          continue;
        }

        int num_visited_targets = pop->visited_targets.cast<int>().sum();
        if (pop->final_pt_idx != 0 && bias_tour.size() &&
            std::isfinite(gtsp_cost_mat(bias_dag_nodes[num_visited_targets - 1].second, pop->final_pt_idx))) {
          bool subset_same = true;
          for (int target_idx = 0; target_idx < num_targets; ++target_idx) {
            if (target_idx == pt_to_target_ptr(pop->final_pt_idx)) {
              if (bias_dag_nodes[num_visited_targets - 1].first(target_idx)) {
                subset_same = false;
                break;
              }
            } else if (pop->visited_targets(target_idx) != bias_dag_nodes[num_visited_targets - 1].first(target_idx)) {
              subset_same = false;
              break;
            }
          }
          if (subset_same) {
            int pt_idx = bias_dag_nodes[num_visited_targets - 1].second;
            BackwardDFSNodePtr neighbor_node = std::make_shared<BackwardDFSNode>(pop, pt_idx, pt_to_target_ptr(pop->final_pt_idx));

            if (closed_list.find(neighbor_node->key) == closed_list.end()) {
              stack.push_back(pop);
              stack.push_back(neighbor_node);
              continue;
            }
          }
        }

        closed_list.insert(pop->key);

        std::vector<std::vector<BackwardDFSNodePtr>> neighbors_per_thread(num_threads);
        std::vector<std::vector<double>> neighbor_sort_vals_per_thread(num_threads);
        if (num_visited_targets == 0) {
          std::vector<long> tour_vec;
          tour_vec.push_back(pop->final_pt_idx);
          BackwardDFSNodePtr node = pop->parent;
          while (node != nullptr) {
            tour_vec.push_back(node->final_pt_idx);
            node = node->parent;
          }
          tour = Map<VectorXl>(tour_vec.data(), tour_vec.size());
          if (tour(0) != 0 || tour(num_targets + 1) != 0) {
            throw std::runtime_error("Tour does not start and end with node 0");
          }
          return profiling_data;
        } else {
          // Generate depot successor
          if (num_visited_targets == 1) {
            if (std::isfinite(gtsp_cost_mat(0, pop->final_pt_idx))) {
              BackwardDFSNodePtr neighbor_node = std::make_shared<BackwardDFSNode>(pop, 0, pt_to_target_ptr(pop->final_pt_idx));
              neighbors_per_thread[0].push_back(neighbor_node);
              neighbor_sort_vals_per_thread[0].push_back(0);
            }
          } else {
            #pragma omp parallel for
            for (int node_idx = 1; node_idx < num_nodes; ++node_idx) {
              int target_idx = pt_to_target_ptr(node_idx);
              if (!pop->visited_targets(target_idx)) {
                continue;
              }
              if (std::isfinite(gtsp_cost_mat(node_idx, pop->final_pt_idx))) {
                if (do_prune) {
                  // tmp_timer_start = std::chrono::high_resolution_clock::now();
                  bool prune = false;
                  for (int target_idx2 = 0; target_idx2 < num_targets; ++target_idx2) {
                    if (target_idx2 != target_idx && after(node_idx, target_idx2) && pop->visited_targets(target_idx2)) {
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

                BackwardDFSNodePtr neighbor_node = std::make_shared<BackwardDFSNode>(pop, node_idx, pop->final_pt_idx == 0 ? -1 : pt_to_target_ptr(pop->final_pt_idx));

                // tmp_timer_stop = std::chrono::high_resolution_clock::now();
                // tmp_nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(tmp_timer_stop - tmp_timer_start).count();
                // profiling_data(3) += ((double)tmp_nanos)/1e9; // node gen time

                if (closed_list.find(neighbor_node->key) != closed_list.end()) {
                  continue;
                }

                neighbors_per_thread[omp_get_thread_num()].push_back(neighbor_node);
                if (pop->final_pt_idx == 0) {
                  // Expand earliest successors first when at root node
                  neighbor_sort_vals_per_thread[omp_get_thread_num()].push_back(-all_pts(node_idx, 0));
                } else {
                  // Expand latest successors first
                  neighbor_sort_vals_per_thread[omp_get_thread_num()].push_back(all_pts(node_idx, 0));
                }
              }
            }
          }
        }

        int num_neighbors = 0;
        for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
          num_neighbors += neighbors_per_thread[thread_idx].size();
        }
        std::vector<double> neighbor_sort_vals(num_neighbors);
        std::vector<BackwardDFSNodePtr> neighbors(num_neighbors);

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
        } else {
          sort_idx.resize(neighbor_sort_vals.size());
          for (int node_idx = 0; node_idx < neighbor_sort_vals.size(); ++node_idx) {
            sort_idx[node_idx] = node_idx;
          }
          std::shuffle(sort_idx.begin(), sort_idx.end(), rng);
          // throw std::runtime_error("Random successor ordering not implemented");
        }

        for (int sort_idx_idx = 0; sort_idx_idx < sort_idx.size(); ++sort_idx_idx) {
          int neighbor_idx = sort_idx[sort_idx_idx];
          // Need to compute costs between points associated with the same target for this pruning method to work
          /*
          if (pop->final_pt_idx != 0 && num_visited_targets != 1) {
            bool prune_this_neighbor = false;
            for (int sort_idx_idx2 = sort_idx_idx + 1; sort_idx_idx2 < sort_idx.size(); ++sort_idx_idx2) {
              int neighbor_idx2 = sort_idx[sort_idx_idx2];
              // Due to the sort, neighbor_idx2 must be later than neighbor_idx
              // If same target and neighbor_idx can get to neighbor_idx2, prune neighbor_idx
              if (pt_to_target_ptr(neighbors[neighbor_idx]->final_pt_idx) == pt_to_target_ptr(neighbors[neighbor_idx2]->final_pt_idx) &&
                  std::isfinite(gtsp_cost_mat(neighbors[neighbor_idx]->final_pt_idx, neighbors[neighbor_idx2]->final_pt_idx))) {
                prune_this_neighbor = true;
                break;
              }
            }
            if (prune_this_neighbor) {
              continue;
            }
          }
          */
          stack.push_back(neighbors[neighbor_idx]);
        }
      }
      tour(0) = -1;
      return profiling_data;
    }
  private:
    MatrixXb after;
    int num_threads;
    int num_targets;
};
