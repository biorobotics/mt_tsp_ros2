#pragma once
#include <chrono>
#include <omp.h>
#include <Eigen/Dense>
#include <random>
#include <set>
#include <iomanip>
#include "mt_tsp_ros2/efat.h"
#include "mt_tsp_ros2/dubins_trj_through_seq_of_targets.h"

using namespace Eigen;

typedef Matrix<double, Dynamic, Dynamic, RowMajor> RowMatrixXd;

typedef Matrix<bool, Dynamic, 1> VectorXb;
typedef Matrix<long, Dynamic, 1> VectorXl;

typedef Matrix<double, 1, 1> Vector1d;

class MemeticPCGUtils {
  public:
    MemeticPCGUtils(const Ref<const RowMatrixXd> &tw_per_target, const std::vector<ExtendedCppSpline> &q_trj_per_target, const Ref<const Vector2d> &p0, double heading0, double vmax, double rho, bool min_latency, bool min_time, double newton_step_size, int num_openmp_threads) {
      int num_targets = tw_per_target.rows();
      int gene_size;
      if (rho == 0) {
        gene_size = 4; // target_idx, time, pos
      } else {
        gene_size = 5; // target_idx, time, pos, heading
      }
      if (rho == 0) {
        if (!min_latency) {
          throw std::runtime_error("Only accepting min-latency for close-enough MT-TSP currently");
        }
        efat_obj = std::make_shared<EFAT>(tw_per_target, q_trj_per_target, p0, vmax, false, 0.);
      } else {
        dubins_obj = std::make_shared<DubinsTrjThroughSeqOfTargets>(tw_per_target, q_trj_per_target, p0, heading0, vmax, rho, false, 0., min_latency, min_time, newton_step_size);
      }

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
      omp_set_num_threads(num_openmp_threads);
    }

    void crossover(Ref<RowMatrixXd> population, Ref<VectorXd> population_costs) {
      int pop_size = population_costs.size();
      int num_targets = population.rows()/pop_size;
      int gene_size = population.cols();

      std::uniform_int_distribution<int> parent_distribution(0, pop_size - 1);
      std::uniform_int_distribution<int> crossover_distribution(0, 1);

      RowMatrixXd updated_population = population;
      VectorXd updated_population_costs = population_costs;

      #pragma omp parallel for
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

        bool success = false;
        Vector1d cost;
        if (efat_obj != nullptr) {
          success = efat_obj->efat_chain(cost, selected_pts_per_target_per_thread[thread_idx], target_seq_per_thread[thread_idx], false);
        } else if (dubins_obj != nullptr) {
          success = dubins_obj->optimize_trj(cost, selected_pts_per_target_per_thread[thread_idx], target_seq_per_thread[thread_idx]);
        } else {
          throw std::runtime_error("Repair objects not initialized");
        }

        if (!success || cost(0) >= population_costs(chromosome_idx)) {
          continue;
        }

        for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
          int target_idx = target_seq_per_thread[thread_idx](seq_idx);
          updated_population.block(updated_start_row + seq_idx, 1, 1, gene_size - 1) = selected_pts_per_target_per_thread[thread_idx].row(target_idx);
          updated_population(updated_start_row + seq_idx, 0) = target_idx;
        }
        updated_population_costs(chromosome_idx) = cost(0);
      }

      population = updated_population;
      population_costs = updated_population_costs;
    }
  private:
    std::shared_ptr<EFAT> efat_obj;
    std::shared_ptr<DubinsTrjThroughSeqOfTargets> dubins_obj;
    std::vector<std::mt19937> rngs_per_thread;
    std::vector<RowMatrixXd> selected_pts_per_target_per_thread;
    std::vector<VectorXl> target_seq_per_thread;
    std::vector<VectorXb> inserted_targets_per_thread;
};
