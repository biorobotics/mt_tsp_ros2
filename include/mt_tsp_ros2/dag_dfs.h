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
          int final_pt_idx,
          int final_target_idx) : parent(parent), 
                                  final_pt_idx(final_pt_idx) {
    if (parent != nullptr) {
      visited_targets = parent->visited_targets;
      visited_targets.insert(final_target_idx);
    }
    std::vector<int> sorted_visited_targets(visited_targets.begin(), visited_targets.end());
    std::sort(sorted_visited_targets.begin(), sorted_visited_targets.end());
    key = VectorXi::Zero(visited_targets.size() + 1);
    for (int i = 0; i < visited_targets.size(); ++i) {
      key(i) = sorted_visited_targets[i];
    }
    key(visited_targets.size()) = final_pt_idx;
  }

  std::shared_ptr<DFSNode> parent;
  std::unordered_set<int> visited_targets;
  int final_pt_idx;
  VectorXi key;
};

typedef std::shared_ptr<DFSNode> DFSNodePtr;

VectorXl dag_dfs(RowMatrixXdRef_const gtsp_cost_mat, VectorXlRef_const pt_to_target_ptr, const std::vector<py::array_t<long>> &target_to_pt_ptr, bool do_prune, bool do_sort, double time_limit, py::object other_tour_queue, RowMatrixXdRef_const all_pts) {
  auto timer_start = std::chrono::high_resolution_clock::now();

  int num_targets = target_to_pt_ptr.size() - 1; // -1 because we have a dummy target associated with the depot

  std::unordered_set<VectorXi, key_hash> closed_list;
  std::vector<DFSNodePtr> stack;
  stack.push_back(std::make_shared<DFSNode>(nullptr, 0, -1));

  int num_nodes = gtsp_cost_mat.rows();

  MatrixXb before(1, 1);
  if (do_prune) {
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
  }

  while (stack.size()) {
    auto timer_stop = std::chrono::high_resolution_clock::now();
    auto micros = std::chrono::duration_cast<std::chrono::microseconds>(timer_stop - timer_start).count();
    if (((double)micros)/1e6 > time_limit) {
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
    if (pop->visited_targets.size() == num_targets) {
      std::vector<long> tour;
      tour.push_back(pop->final_pt_idx);
      DFSNodePtr node = pop->parent;
      while (node != nullptr) {
        tour.push_back(node->final_pt_idx);
        node = node->parent;
      }
      std::reverse(tour.begin(), tour.end());
      // Assume open tsp
      tour.push_back(0);
      return Map<VectorXl>(tour.data(), tour.size());
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
      sort_idx = sort_indexes(neighbor_times);
      std::reverse(sort_idx.begin(), sort_idx.end());
    } else {
      throw std::runtime_error("Random successor ordering not implemented");
    }
    for (int neighbor_idx : sort_idx) {
      int pt_idx = neighbors[neighbor_idx];
      if (do_prune) {
        bool prune = false;
        for (int target_idx = 0; target_idx < num_targets; ++target_idx) {
          if (target_idx != pt_to_target_ptr(pt_idx) && before(pt_idx, target_idx) && pop->visited_targets.find(target_idx) == pop->visited_targets.end()) {
            prune = true;
            break;
          }
        }

        if (prune) {
          continue;
        }
      }

      DFSNodePtr neighbor_node = std::make_shared<DFSNode>(pop, pt_idx, pt_to_target_ptr(pt_idx));

      if (closed_list.find(neighbor_node->key) != closed_list.end()) {
        continue;
      }

      stack.push_back(neighbor_node);
    }
  }
  return -1*VectorXl::Ones(1);
}
