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

// THIS IS UNFINISHED

typedef const Ref<const Matrix<long, Dynamic, Dynamic, RowMajor>> &RowMatrixXlRef_const;
typedef const Ref<const Matrix<bool, Dynamic, Dynamic, RowMajor>> &RowMatrixXbRef_const;

struct MultiAgentDFSNode { 
  MultiAgentDFSNode(std::shared_ptr<MultiAgentDFSNode> parent, 
                    int final_pt_idx,
                    int final_target_idx,
                    int final_agent_idx,
                    int num_targets = 0) : parent(parent), 
                                           final_agent_idx(final_agent_idx),
                                           final_pt_idx(final_pt_idx) {
    if (parent != nullptr) {
      visited_targets = parent->visited_targets;
      num_targets = visited_targets.size();
      if (final_target_idx != num_targets) {
        visited_targets(final_target_idx) = true;
      }
    } else {
      if (num_targets == 0) {
        throw std::runtime_error("num_targets can't be zero if parent == nullptr");
      }
      visited_targets = VectorXb::Zero(num_targets);
    }
    num_targets = visited_targets.size();
    key = VectorXi::Zero(num_targets + 2);
    for (int i = 0; i < num_targets; ++i) {
      key(i) = visited_targets(i);
    }
    key(num_targets) = final_pt_idx;
    key(num_targets + 1) = final_agent_idx;
  }

  std::shared_ptr<MultiAgentDFSNode> parent;
  VectorXb visited_targets;
  int final_pt_idx;
  int final_agent_idx;
  VectorXi key;
};

typedef std::shared_ptr<MultiAgentDFSNode> MultiAgentDFSNodePtr;

class LifelongMultiAgentDAGDFSPlanner {
  public:
    LifelongMultiAgentDAGDFSPlanner(RowMatrixXdRef_const gtsp_cost_mat, VectorXlRef_const pt_to_target_ptr, const std::vector<py::array_t<long>> &target_to_pt_ptr, bool do_prune, bool do_sort, RowMatrixXdRef_const all_pts, int num_threads, int num_agents) : num_threads(num_threads), num_agents(num_agents) {
      before_time = 0.;
      num_targets = target_to_pt_ptr.size() - 1; // -1 because we have a dummy target associated with the depots
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

        for (int node_idx = num_agents; node_idx < num_nodes; ++node_idx) {
          before(node_idx, pt_to_target_ptr(node_idx)) = false;
        }
        auto timer_stop = std::chrono::high_resolution_clock::now();
        auto nanos = std::chrono::duration_cast<std::chrono::microseconds>(timer_stop - timer_start).count();
        before_time = ((double)nanos)/1e9;

        // By default, assume all agents can visit all targets
        backwards_cumulative_visitable_targets_per_agent = MatrixXb::Ones(num_agents, num_targets);
        for (int agent_idx = num_agents - 1; agent_idx >= 0; --agent_idx) {
          for (int target_idx = 0; target_idx < num_targets; ++target_idx) {
            if (!before(agent_idx, target_idx)) {
              // If agent_idx can visit target_idx, we're done
              continue;
            }

            int next_agent_idx = agent_idx + 1;
            if (next_agent_idx != num_agents && backwards_cumulative_visitable_targets_per_agent(next_agent_idx, target_idx)) {
              // If some other agent with larger index can visit target_idx, we're done
              continue;
            }

            backwards_cumulative_visitable_targets_per_agent(agent_idx, target_idx) = false;
          }
        }
      }
    }

    VectorXd plan(Ref<VectorXl> tour, RowMatrixXdRef_const gtsp_cost_mat, VectorXlRef_const pt_to_target_ptr, const std::vector<py::array_t<long>> &target_to_pt_ptr, bool do_prune, bool do_sort, double time_limit, py::object other_tour_queue, RowMatrixXdRef_const all_pts, RowMatrixXlRef_const deleted_edges, bool sort_by_time) {
      auto timer_start = std::chrono::high_resolution_clock::now();
      VectorXd profiling_data = VectorXd::Zero(4);
      int num_targets = target_to_pt_ptr.size() - 1; // -1 because we have a dummy target associated with the depots

      std::mt19937 rng;

      if (do_prune) {
        // Update the before sets
        auto tmp_timer_start = std::chrono::high_resolution_clock::now();
        for (int edge_idx = 0; edge_idx < deleted_edges.rows(); ++edge_idx) {
          int node_idx1 = deleted_edges(edge_idx, 0);
          int node_idx2 = deleted_edges(edge_idx, 1);

          if (node_idx2 < num_agents) {
            throw std::runtime_error("We deleted an edge entering a depot");
          }

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

        if (deleted_edges.rows()) {
          // By default, assume all agents can visit all targets
          backwards_cumulative_visitable_targets_per_agent = MatrixXb::Ones(num_agents, num_targets);
          for (int agent_idx = num_agents - 1; agent_idx >= 0; --agent_idx) {
            for (int target_idx = 0; target_idx < num_targets; ++target_idx) {
              if (!before(agent_idx, target_idx)) {
                // If agent_idx can visit target_idx, we're done
                continue;
              }

              int next_agent_idx = agent_idx + 1;
              if (next_agent_idx != num_agents && backwards_cumulative_visitable_targets_per_agent(next_agent_idx, target_idx)) {
                // If some other agent with larger index can visit target_idx, we're done
                continue;
              }

              backwards_cumulative_visitable_targets_per_agent(agent_idx, target_idx) = false;
            }
          }
        }

        auto tmp_timer_stop = std::chrono::high_resolution_clock::now();
        auto tmp_nanos = std::chrono::duration_cast<std::chrono::microseconds>(tmp_timer_stop - tmp_timer_start).count();
        profiling_data(0) = ((double)tmp_nanos)/1e9;
        before_time += profiling_data(0);
      }

      std::unordered_set<VectorXi, key_hash> closed_list;
      std::vector<MultiAgentDFSNodePtr> stack;
      stack.push_back(std::make_shared<MultiAgentDFSNode>(nullptr, 0, -1, 0, num_targets));

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

        MultiAgentDFSNodePtr pop = stack.back();
        stack.pop_back();

        if (closed_list.find(pop->key) != closed_list.end()) {
          continue;
        }

        int num_visited_targets = pop->visited_targets.cast<int>().sum();

        closed_list.insert(pop->key);

        std::vector<std::vector<MultiAgentDFSNodePtr>> neighbors_per_thread(num_threads);
        std::vector<std::vector<double>> neighbor_sort_vals_per_thread(num_threads);
        if (num_visited_targets == num_targets) {
          std::vector<long> tour_vec;
          tour_vec.push_back(pop->final_pt_idx);
          MultiAgentDFSNodePtr node = pop->parent;
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
          // Add successors associated with nodes of unvisited targets
          #pragma omp parallel for
          for (int node_idx = num_agents; node_idx < num_nodes; ++node_idx) {
            int target_idx = pt_to_target_ptr(node_idx);
            if (pop->visited_targets(target_idx)) {
              continue;
            }
            if (std::isfinite(gtsp_cost_mat(pop->final_pt_idx, node_idx))) {
              if (do_prune) {
                int next_agent_idx = pop->final_agent_idx + 1;

                // tmp_timer_start = std::chrono::high_resolution_clock::now();
                bool prune = false;
                for (int target_idx2 = 0; target_idx2 < num_targets; ++target_idx2) {
                  if (target_idx2 != target_idx && before(node_idx, target_idx2) && !pop->visited_targets(target_idx2) && (next_agent_idx == num_agents || !backwards_cumulative_visitable_targets_per_agent(next_agent_idx, target_idx2))) {
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

              MultiAgentDFSNodePtr neighbor_node = std::make_shared<MultiAgentDFSNode>(pop, node_idx, target_idx, pop->final_agent_idx);

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
        std::vector<MultiAgentDFSNodePtr> neighbors(num_neighbors);

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

        // Prioritize going to next depot last
        if (pop->final_agent_idx != num_agents - 1) {
          MultiAgentDFSNodePtr neighbor_node = std::make_shared<MultiAgentDFSNode>(pop, pop->final_agent_idx + 1, num_targets, pop->final_agent_idx + 1);
          if (closed_list.find(neighbor_node->key) == closed_list.end()) {
            stack.push_back(neighbor_node);
          }
        }

        for (int neighbor_idx : sort_idx) {
          stack.push_back(neighbors[neighbor_idx]);
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
    MatrixXb backwards_cumulative_visitable_targets_per_agent;
    double before_time;
    int num_threads;
    int num_targets;
    int num_agents;
};
