#pragma once
#include <unordered_map>
#include <vector>
#include <chrono>
#include <Eigen/Dense>
#include <omp.h>

using namespace Eigen;

typedef const Ref<const Matrix<long, Dynamic, Dynamic, RowMajor>> &RowMatrixXlRef_const;

class AdjacencyList {
  public:
    void add_edge(int source, int dest) {
      if (adjacency_list.find(source) == adjacency_list.end()) {
        adjacency_list[source] = std::vector<int>();
      }
      adjacency_list[source].push_back(dest);
    }

    VectorXi dfs(int source, int dest) {
      std::vector<int> stack;
      std::unordered_map<int, int> backpointers;
      stack.push_back(source);
      while (stack.size()) {
        int pop = stack.back();
        stack.pop_back();
        if (pop == dest) {
          // Backpointer traversal
          std::vector<int> path;
          int node = pop;
          while (node != source) {
            path.push_back(node);
            node = backpointers[node];
          }
          path.push_back(node);
          std::reverse(path.begin(), path.end());
          return Map<VectorXi>(path.data(), path.size());
        }

        if (adjacency_list.find(pop) == adjacency_list.end()) {
          continue;
        }
        for (int succ : adjacency_list[pop]) {
          if (backpointers.find(succ) != backpointers.end()) {
            continue;
          }
          backpointers[succ] = pop;
          stack.push_back(succ);
        }
      }
      return VectorXi::Zero(0);
    }

    void batch_dfs(RowMatrixXlRef_const source_dest_pairs, int num_threads, std::vector<VectorXi> &paths) {
      paths.resize(source_dest_pairs.rows());
      omp_set_num_threads(num_threads);
      #pragma omp parallel for
      for (int pair_idx = 0; pair_idx < source_dest_pairs.rows(); ++pair_idx) {
        int source = source_dest_pairs(pair_idx, 0);
        int dest = source_dest_pairs(pair_idx, 1);
        paths[pair_idx] = dfs(source, dest);
      }
    }

  private:
    std::unordered_map<int, std::vector<int>> adjacency_list;
};
