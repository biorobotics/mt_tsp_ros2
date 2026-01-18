#pragma once
#include <chrono>
#include <omp.h>
#include <Eigen/Dense>
#include <random>
#include <set>
#include <iomanip>
#include "mt_tsp_ros2/reopt_gtsp_tour.h"

using namespace Eigen;

typedef Matrix<double, Dynamic, Dynamic, RowMajor> RowMatrixXd;

typedef Matrix<bool, Dynamic, 1> VectorXb;
typedef Matrix<long, Dynamic, 1> VectorXl;

typedef Matrix<double, 1, 1> Vector1d;

bool repair_chromosome(const Ref<const VectorXl> &target_seq, const py::array_t<long> &gtsp_cost_mat_rounded_and_scaled_flat, const py::array_t<double> &gtsp_cost_mat_flat, const std::vector<py::array_t<long>> &target_to_pt_ptr, const py::array_t<double> &all_pts, double &cost, long inf_val, int num_nodes, Ref<VectorXl> tour) {
  bool repair_succeeded = reopt_gtsp_tour(tour, gtsp_cost_mat_rounded_and_scaled_flat, target_to_pt_ptr, target_seq, inf_val, num_nodes);
  if (!repair_succeeded) {
    cost = std::numeric_limits<double>::infinity();
    return true;
  }
  cost = 0.;
  int prev_node_idx = 0;

  auto gtsp_cost_mat_rounded_and_scaled_flat_unchecked = gtsp_cost_mat_rounded_and_scaled_flat.unchecked<1>();
  auto gtsp_cost_mat_flat_unchecked = gtsp_cost_mat_flat.unchecked<1>();
  auto all_pts_unchecked = all_pts.unchecked<2>();

  int pt_dim = all_pts_unchecked.shape(1);

  for (int tour_idx = 1; tour_idx < tour.size() - 1; ++tour_idx) {
    int node_idx = tour(tour_idx);
    cost += gtsp_cost_mat_flat_unchecked(prev_node_idx*num_nodes + node_idx);
    // std::cout << prev_node_idx << " " << node_idx << " " << num_nodes << " " << gtsp_cost_mat_flat_unchecked(prev_node_idx*num_nodes + node_idx) << " " << gtsp_cost_mat_rounded_and_scaled_flat_unchecked(prev_node_idx*num_nodes + node_idx) << std::endl;
    prev_node_idx = node_idx;
  }
  if (std::isinf(cost)) {
    throw std::runtime_error("reopt_gtsp_tour returned infinite cost solution but reported success");
  }
  return false;
}

class MemeticPCGUtilsSampling {
  public:
    MemeticPCGUtilsSampling(int num_openmp_threads, double mutation_prob, int num_targets, int gene_size) : num_openmp_threads(num_openmp_threads), mutation_prob(mutation_prob) {
      std::vector<unsigned int> seeds{4223100027, 870236586, 1683737518, 3430707182, 1613429085, 1714341085, 853110547, 1988005940, 2629786018, 1139192408};

      if (num_openmp_threads > seeds.size()) {
        throw std::runtime_error("Too many threads, not enough stored random seeds");
      }
      for (int thread_idx = 0; thread_idx < num_openmp_threads; ++thread_idx) {
        rngs_per_thread.push_back(std::mt19937(seeds[thread_idx]));
        // Preallocate
        selected_pts_per_target_per_thread.push_back(RowMatrixXd::Zero(num_targets, gene_size - 1));
        target_seq_per_thread.push_back(VectorXl::Zero(num_targets));
        inserted_targets_per_thread.push_back(VectorXb::Zero(num_targets));
      }
    }

    void crossover(Ref<RowMatrixXd> population, Ref<VectorXd> population_costs, const py::array_t<long> &gtsp_cost_mat_rounded_and_scaled_flat, const py::array_t<double> &gtsp_cost_mat_flat, const std::vector<py::array_t<long>> &target_to_pt_ptr, const py::array_t<double> &all_pts, long inf_val, int num_nodes) {
      omp_set_num_threads(num_openmp_threads);
      int pop_size = population_costs.size();
      int num_targets = population.rows()/pop_size;
      int gene_size = population.cols();

      std::uniform_int_distribution<int> parent_distribution(0, pop_size - 1);
      std::uniform_int_distribution<int> crossover_distribution(0, 1);
      std::uniform_real_distribution<double> mutation_distribution(0, 1);
      std::uniform_int_distribution<int> mutation_operator1_distribution(0, num_targets - 1);

      RowMatrixXd updated_population = population;
      VectorXd updated_population_costs = population_costs;

      // #pragma omp parallel for
      for (int chromosome_idx = 0; chromosome_idx < pop_size; ++chromosome_idx) {
        int parent1_idx = chromosome_idx;
        int parent2_idx = parent_distribution(rngs_per_thread[omp_get_thread_num()]);
        // Crossover

        int parent1_start_row = parent1_idx*num_targets;
        int parent2_start_row = parent2_idx*num_targets;

        int parent1_counter = parent1_start_row;
        int parent2_counter = parent2_start_row;

        int updated_start_row = parent1_start_row;

        int thread_idx = omp_get_thread_num();
        inserted_targets_per_thread[thread_idx].setZero();

        for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
          if (crossover_distribution(rngs_per_thread[omp_get_thread_num()]) == 0) {
            while (inserted_targets_per_thread[thread_idx]((int)(population(parent1_counter, 0)))) {
              ++parent1_counter;
            }
            if (parent1_counter >= parent1_start_row + num_targets) {
              throw std::runtime_error("parent 1 out of bounds");
            }
            int target_idx = population(parent1_counter, 0);
            selected_pts_per_target_per_thread[thread_idx].row(target_idx) = population.block(parent1_counter, 1, 1, gene_size - 1);
            target_seq_per_thread[thread_idx](seq_idx) = target_idx;
            inserted_targets_per_thread[thread_idx](target_idx) = true;
          } else {
            while (inserted_targets_per_thread[thread_idx]((int)(population(parent2_counter, 0)))) {
              ++parent2_counter;
            }
            if (parent2_counter >= parent2_start_row + num_targets) {
              throw std::runtime_error("parent 2 out of bounds");
            }
            int target_idx = population(parent2_counter, 0);
            selected_pts_per_target_per_thread[thread_idx].row(target_idx) = population.block(parent2_counter, 1, 1, gene_size - 1);
            target_seq_per_thread[thread_idx](seq_idx) = target_idx;
            inserted_targets_per_thread[thread_idx](target_idx) = true;
          }
        }

        double mutation_sample = mutation_distribution(rngs_per_thread[omp_get_thread_num()]);
        if (mutation_sample < mutation_prob) {
          int seq_idx1 = mutation_operator1_distribution(rngs_per_thread[omp_get_thread_num()]);
          int seq_idx2 = mutation_operator1_distribution(rngs_per_thread[omp_get_thread_num()]);
          int target_idx1 = target_seq_per_thread[thread_idx](seq_idx1);
          int target_idx2 = target_seq_per_thread[thread_idx](seq_idx2);
          target_seq_per_thread[thread_idx](seq_idx1) = target_idx2;
          target_seq_per_thread[thread_idx](seq_idx2) = target_idx1;
        }

        double cost = -1.;
        VectorXl tour = -VectorXl::Ones(num_targets + 2);
        bool repair_failed = repair_chromosome(target_seq_per_thread[thread_idx], gtsp_cost_mat_rounded_and_scaled_flat, gtsp_cost_mat_flat, target_to_pt_ptr, all_pts, cost, inf_val, num_nodes, tour);
        bool success = !repair_failed;
        if (cost == -1.) {
          throw std::runtime_error("Did not update cost");
        }
        if (std::isfinite(cost) && tour.minCoeff() == -1) {
          throw std::runtime_error("Did not populate some element of tour");
        }

        if (!success || cost >= population_costs(chromosome_idx)) {
          continue;
        }

        for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
          int target_idx = target_seq_per_thread[thread_idx](seq_idx);
          updated_population.block(updated_start_row + seq_idx, 1, 1, gene_size - 1) = selected_pts_per_target_per_thread[thread_idx].row(target_idx);
          updated_population(updated_start_row + seq_idx, 0) = target_idx;
        }
        updated_population_costs(chromosome_idx) = cost;
      }

      population = updated_population;
      population_costs = updated_population_costs;
    }
  private:
    std::vector<std::mt19937> rngs_per_thread;
    std::vector<RowMatrixXd> selected_pts_per_target_per_thread;
    std::vector<VectorXl> target_seq_per_thread;
    std::vector<VectorXb> inserted_targets_per_thread;
    int num_openmp_threads;
    double mutation_prob;
};
