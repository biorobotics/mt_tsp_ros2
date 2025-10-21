#pragma once
#include <memory>
#include <unordered_set>
#include <set>
#include <vector>
#include <chrono>
#include <random>
#include <Eigen/Dense>
#include "mt_tsp_ros2/cpp_spline.h"

using namespace Eigen;
namespace py = pybind11;

typedef const Ref<const Matrix<long, Dynamic, 1>> &VectorXlRef_const;
typedef const Ref<const Matrix<double, Dynamic, Dynamic, RowMajor>> &RowMatrixXdRef_const;
typedef const Ref<const Matrix<bool, Dynamic, Dynamic, RowMajor>> &RowMatrixXbRef_const;
typedef Matrix<bool, Dynamic, Dynamic> MatrixXb;
typedef Matrix<bool, Dynamic, 1> VectorXb;
typedef Matrix<long, Dynamic, 1> VectorXl;

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values
  std::stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

// Taken from https://jimmy-shen.medium.com/stl-map-unordered-map-with-a-vector-for-the-key-f30e5f670bae#:~:text=unordered_map%20uses%20vector%20as%20the%20key&text=You%20can%20use%20the%20following,make%20the%20best%20of%20STL.&text=%7D%3B,so%20that%20collisions%20are%20minimized
struct key_hash {
  int operator()(const VectorXi &k) const {
    int hash = k.size();
    for(int i = 0; i < k.size(); ++i) {
      hash ^= k(i) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
    }
    return hash;
  }
};

struct DFSNode { 
  DFSNode(std::shared_ptr<DFSNode> parent, 
          double final_time,
          const Ref<const Vector2d> &final_pos,
          int final_target_idx,
          int num_targets = 0) : parent(parent), 
                                 final_time(final_time),
                                 final_pos(final_pos),
                                 final_target_idx(final_target_idx) {
    if (parent != nullptr) {
      visited_targets = parent->visited_targets;
      visited_targets(final_target_idx) = true;
    } else {
      if (num_targets == 0) {
        throw std::runtime_error("num_targets can't be zero if parent == nullptr");
      }
      visited_targets = VectorXb::Zero(num_targets);
    }
    key = VectorXi::Zero(visited_targets.size() + 1);
    key.head(visited_targets.size()) = visited_targets.cast<int>();
    key(num_targets) = final_target_idx;
  }

  std::shared_ptr<DFSNode> parent;
  VectorXb visited_targets;
  double final_time;
  Vector2d final_pos;
  int final_target_idx;
  VectorXi key;
};


double sft(Ref<Vector2d> next_pos, const Ref<const Vector2d> &pos, double t, double tw_start, double tw_end, double vmax, const CppSpline& trj, bool negate_trj_input) {
  double time_multiplier = negate_trj_input ? -1 : 1;
  double next_t = tw_end;
  double delta_t = next_t - t;
  next_pos = trj(time_multiplier*next_t);
  double dist = (next_pos - pos).norm();

  if (dist > vmax*delta_t) {
    return std::numeric_limits<double>::infinity();
  }

  // Check if we can intercept at start of time window
  double start_t = tw_start;
  Vector2d start_pos = trj(time_multiplier*start_t);
  double start_dist = (start_pos - pos).norm();
  double start_delta_t = start_t - t;
  if (start_dist <= vmax*start_delta_t) {
    next_pos = start_pos;
    return start_delta_t;
  }

  int max_newton_iter = 100;
  double feas_next_t = next_t;
  Vector2d feas_next_pos = next_pos;
  double feas_dist = dist;
  bool newton_success = false;
  for (int newton_iter = 0; newton_iter < max_newton_iter; ++newton_iter) {
    // Find root of dist - vmax*delta_t
    double resid = dist - vmax*delta_t + 1e-4; // The 1e-4 is so we actually get to a feasible solution
    double deriv = 1/dist*(next_pos - pos).dot(time_multiplier*trj.derivatives(time_multiplier*next_t)) - vmax;
    delta_t -= resid/deriv;
    next_t = t + delta_t;
    next_pos = trj(time_multiplier*next_t);
    dist = (next_pos - pos).norm();
    if (dist <= vmax*delta_t) {
      feas_next_t = next_t;
      feas_next_pos = next_pos;
      feas_dist = dist;

      if (std::abs(resid) < 1e-4) {
        newton_success = true;
        break;
      }
    }
  }

  if (!newton_success) {
    throw std::runtime_error("Newton did not converge. Increase the max number of iterations");
  }
  next_pos = feas_next_pos;
  return feas_next_t - t;
}

typedef std::shared_ptr<DFSNode> DFSNodePtr;

VectorXd continuous_dag_dfs_no_obstacles(Ref<VectorXl> target_seq, Ref<VectorXd> time_seq, const std::vector<CppSpline> &q_trj_per_target_python, const Ref<const RowMatrixXd> &tw_per_target, double vmax, const Ref<const Vector2d> &p0, double t0, int num_threads, bool do_prune, bool do_sort, double time_limit) {
  auto timer_start = std::chrono::high_resolution_clock::now();
  VectorXd profiling_data = VectorXd::Zero(4);
  int num_targets = tw_per_target.rows();

  // I'm doing this regardless of whether there are time windows because I'm worried about the GIL
  std::vector<CppSpline> q_trj_per_target;
  for (int target_idx = 0; target_idx < num_targets; ++target_idx) {
    q_trj_per_target.push_back(CppSpline(q_trj_per_target_python[target_idx].get_knots(), q_trj_per_target_python[target_idx].get_coeffs()));
  }

  std::mt19937 rng;

  MatrixXd lfdt_mat(1, 1);

  if (do_prune) {
    lfdt_mat = MatrixXd(num_targets + 1, num_targets);
    // Compute LFDT matrix
    #pragma omp parallel for
    for (int target_idx1 = 0; target_idx1 < num_targets; ++target_idx1) {
      #pragma omp parallel for
      for (int target_idx2 = 0; target_idx2 < num_targets; ++target_idx2) {
        double t = -tw_per_target(target_idx2, 1);
        Vector2d pos = q_trj_per_target[target_idx2](tw_per_target(target_idx2, 1));
        Vector2d next_pos;
        double delta_t = sft(next_pos, pos, t, -tw_per_target(target_idx1, 1), -tw_per_target(target_idx1, 0), vmax, q_trj_per_target[target_idx1], true);
        lfdt_mat(target_idx1, target_idx2) = tw_per_target(target_idx2, 1) - delta_t;
      }
    }
    lfdt_mat.row(num_targets).setConstant(t0); // Leave p0 immediately
  }

  std::unordered_map<VectorXi, DFSNodePtr, key_hash> earliest_node_per_key;
  std::vector<DFSNodePtr> stack;
  stack.push_back(std::make_shared<DFSNode>(nullptr, t0, p0, num_targets, num_targets));

  omp_set_num_threads(num_threads);

  while (stack.size()) {
    auto timer_stop = std::chrono::high_resolution_clock::now();
    auto nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
    if (((double)nanos)/1e9 > time_limit) {
      break;
    }

    DFSNodePtr pop = stack.back();
    stack.pop_back();

    auto it = earliest_node_per_key.find(pop->key);
    if (it != earliest_node_per_key.end() && it->second != pop && it->second->final_time <= pop->final_time) {
      continue;
    }

    int num_visited_targets = pop->visited_targets.cast<int>().sum();

    std::vector<std::vector<DFSNodePtr>> neighbors_per_thread(num_threads);
    std::vector<std::vector<double>> neighbor_sort_vals_per_thread(num_threads);
    if (num_visited_targets == num_targets) {
      std::vector<double> time_seq_vec;
      std::vector<long> target_seq_vec;
      time_seq_vec.push_back(pop->final_time);
      target_seq_vec.push_back(pop->final_target_idx);
      DFSNodePtr node = pop->parent;
      std::cout << "continuous-time DAG-DFS found solution" << std::endl;
      while (true) {
        time_seq_vec.push_back(node->final_time);
        if (node->final_target_idx == num_targets) {
          // Back to depot
          break;
        }
        target_seq_vec.push_back(node->final_target_idx);
        node = node->parent;
      }
      std::reverse(time_seq_vec.begin(), time_seq_vec.end());
      std::reverse(target_seq_vec.begin(), target_seq_vec.end());
      time_seq = Map<VectorXd>(time_seq_vec.data(), time_seq_vec.size());
      target_seq = Map<VectorXl>(target_seq_vec.data(), target_seq_vec.size());
      return profiling_data;
    } else {
      #pragma omp parallel for
      for (int target_idx = 0; target_idx < num_targets; ++target_idx) {
        if (pop->visited_targets(target_idx)) {
          continue;
        }

        Vector2d next_pos;
        double next_t = pop->final_time + sft(next_pos, pop->final_pos, pop->final_time, tw_per_target(target_idx, 0), tw_per_target(target_idx, 1), vmax, q_trj_per_target[target_idx], false);

        if (std::isinf(next_t) && do_prune) {
          std::cout << "pop->final_time = " << pop->final_time << std::endl;
          std::cout << "lfdt_mat(pop->final_target_idx, target_idx) = " << lfdt_mat(pop->final_target_idx, target_idx) << std::endl;
          throw std::runtime_error("We are starting no later than the LFDT to this target but sft says travel is infeasible");
        }

        if (do_prune) {
          // tmp_timer_start = std::chrono::high_resolution_clock::now();
          bool prune = false;
          for (int target_idx2 = 0; target_idx2 < num_targets; ++target_idx2) {
            if (target_idx2 != target_idx && next_t > lfdt_mat(target_idx, target_idx2) && !pop->visited_targets(target_idx2)) {
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

        DFSNodePtr neighbor_node = std::make_shared<DFSNode>(pop, next_t, next_pos, target_idx);

        // tmp_timer_stop = std::chrono::high_resolution_clock::now();
        // tmp_nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(tmp_timer_stop - tmp_timer_start).count();
        // profiling_data(3) += ((double)tmp_nanos)/1e9; // node gen time

        auto it = earliest_node_per_key.find(neighbor_node->key);
        if (it != earliest_node_per_key.end() && it->second->final_time <= pop->final_time) {
          continue;
        }

        neighbors_per_thread[omp_get_thread_num()].push_back(neighbor_node);
        neighbor_sort_vals_per_thread[omp_get_thread_num()].push_back(next_t);
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
      earliest_node_per_key[neighbors[neighbor_idx]->key] = neighbors[neighbor_idx];
    }
  }
  time_seq(0) = std::numeric_limits<double>::infinity();
  target_seq(0) = -1;
  return profiling_data;
}
