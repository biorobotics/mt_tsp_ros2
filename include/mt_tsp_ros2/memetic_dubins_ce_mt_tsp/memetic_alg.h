#pragma once
#include <chrono>
#include <omp.h>
#include <Eigen/Dense>
#include <pybind11/pybind11.h>
#include <random>
#include <set>

using namespace Eigen;
namespace py = pybind11;

typedef Matrix<double, Dynamic, Dynamic, RowMajor> RowMatrixXd;

typedef Matrix<bool, Dynamic, 1> VectorXb;

// initial_population should have number of rows equal to num_targets*pop_size, and gene_size columns
double memetic_alg(Ref<RowMatrixXd> selected_pts_per_target, const Ref<const RowMatrixXd> &initial_population, const Ref<const VectorXd> &initial_costs, const std::vector<py::object> &q_trj_per_target, const Ref<const RowMatrixXd> &tw_per_target, const Ref<const VectorXd> &target_radii, double time_limit, int pop_size, double rho, double vmax, const Ref<const Vector2d> &p0, double theta0) {
  auto timer_start = std::chrono::high_resolution_clock::now();

  bool dubins = rho != 0;

  int num_targets = q_trj_per_target.size();

  std::vector<MatrixXd> population(pop_size);
  std::vector<double> population_costs(pop_size);
  Map<VectorXd>(population_costs.data(), pop_size) = initial_costs;

  std::mt19937 rng;
  std::uniform_int_distribution<int> parent_distribution(0, pop_size - 1);
  std::uniform_int_distribution<int> crossover_distribution(0, 1);
  std::uniform_real_distribution<double> mutation_distribution(0, 1);
  std::uniform_int_distribution<int> mutation_operator_distribution(0, 2);

  std::uniform_int_distribution<int> mutation_operator1_distribution(0, num_targets - 1);

  std::uniform_int_distribution<int> mutation_operator2_seq_idx_distribution(0, num_targets - 1);
  std::uniform_real_distribution<double> mutation_operator2_angle_distribution(0, 2*M_PI);

  std::uniform_int_distribution<int> mutation_operator3_seq_idx_distribution(0, num_targets - 1);
  std::uniform_real_distribution<double> mutation_operator3_delta_t_distribution(0, 1);

  double mutation_prob = 0.5;
  int Tlp = 1; // TODO, how often do we run local search

  const int gene_size = 3; // target index, theta, and delta t

  // Initialize population
  for (int i = 0; i < pop_size; ++i) {
    population[i] = initial_population.block(num_targets*i, 0, num_targets, gene_size);
  }

  int gen_idx = 0;
  while (true) {
    auto timer_stop = std::chrono::high_resolution_clock::now();
    auto nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
    if (((double)nanos)/1e9 > time_limit) {
      break;
    }

    for (int chromosome_idx = 0; chromosome_idx < pop_size; ++chromosome_idx) {
      int parent1_idx = chromosome_idx;
      int parent2_idx = parent_distribution(rng);
      // Crossover
      VectorXb inserted_targets = VectorXb::Zero(num_targets);
      MatrixXd Xnew(num_targets, gene_size);
      int parent1_counter = 0;
      int parent2_counter = 0;
      for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
        if (crossover_distribution(rng) == 0) {
          while (inserted_targets((int)(population[parent1_idx](parent1_counter, 0)))) {
            ++parent1_counter;
          }
          if (parent1_counter >= num_targets) {
            throw std::runtime_error("parent 1 out of bounds");
          }
          int target_idx = population[parent1_idx](parent1_counter, 0);
          Xnew.row(seq_idx) = population[parent1_idx].row(parent1_counter);
          inserted_targets(target_idx) = true;
        } else {
          while (inserted_targets((int)(population[parent2_idx](parent2_counter, 0)))) {
            ++parent2_counter;
          }
          if (parent2_counter >= num_targets) {
            throw std::runtime_error("parent 2 out of bounds");
          }
          int target_idx = population[parent2_idx](parent2_counter, 0);
          Xnew.row(seq_idx) = population[parent2_idx].row(parent2_counter);
          inserted_targets(target_idx) = true;
        }
      }

      if (mutation_distribution(rng) < mutation_prob) {
        // Mutation
        int mutation_operator = mutation_operator_distribution(rng);
        if (mutation_operator == 0) {
          int seq_idx1 = mutation_operator1_distribution(rng);
          int seq_idx2 = mutation_operator1_distribution(rng);
          RowVectorXd tmp = Xnew.row(seq_idx1);
          Xnew.row(seq_idx1) = Xnew.row(seq_idx2);
          Xnew.row(seq_idx2) = tmp;
        } else if (mutation_operator == 1) {
          int seq_idx = mutation_operator2_seq_idx_distribution(rng);
          double theta = mutation_operator2_angle_distribution(rng);
          Xnew(seq_idx, 1) = theta;
        } else {
          int seq_idx = mutation_operator3_seq_idx_distribution(rng);
          double raw_sample = mutation_operator3_delta_t_distribution(rng);
          int target_idx = Xnew(seq_idx, 0);
          double min_t = tw_per_target(target_idx, 0);
          double max_t = tw_per_target(target_idx, 1);
          double delta_t = (max_t - min_t)*raw_sample;
          Xnew(seq_idx, 2) = delta_t;
        }
      }

      // Repair to restore feasibility
      bool repair_failed = false;
      double t = 0;
      Vector2d pos = p0;
      Vector2d next_rel_pos;
      Vector2d next_pos;
      double cost = 0;
      for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
        int target_idx = Xnew(seq_idx, 0);
        double theta = Xnew(seq_idx, 1);
        double delta_t = Xnew(seq_idx, 2);
        if (dubins) {
          throw std::runtime_error("Did not implement Dubins yet");
        } else {
          double next_t = t + delta_t;
          if (next_t > tw_per_target(target_idx, 1)) {
            next_t = tw_per_target(target_idx, 1);
            delta_t = next_t - t;
            if (delta_t < 0) {
              repair_failed = true;
              break;
            }
            Xnew(seq_idx, 2) = delta_t;
          } else if (next_t < tw_per_target(target_idx, 0)) {
            next_t = tw_per_target(target_idx, 0);
            delta_t = next_t - t;
            Xnew(seq_idx, 2) = delta_t;
          }

          // Check if interception is feasible
          next_rel_pos = target_radii[target_idx]*Vector2d(cos(theta), sin(theta));
          q_trj_per_target[target_idx].attr("__call__")(next_t, std::ref(next_pos));
          next_pos = next_pos + next_rel_pos;
          double dist = (next_pos - pos).norm();
          if (dist <= vmax*delta_t) {
            // Travel is feasible, no need to repair
            cost += dist;
            t = next_t;
            pos = next_pos;

            continue;
          }

          // Check if travel is feasible to next_rel_pos at end of time window
          q_trj_per_target[target_idx].attr("__call__")(tw_per_target(target_idx, 1), std::ref(next_pos));
          next_pos = next_pos + next_rel_pos;
          dist = (next_pos - pos).norm();
          if (dist > vmax*delta_t) {
            // Travel is infeasible even to end of time window
            repair_failed = true;
            break;
          }

          // Run bisection to find earliest time such that interception is feasible
          double t_low = next_t;
          double t_high = tw_per_target(target_idx, 1);
          int num_bisection_iter = 10;
          double t_high_dist = dist;
          for (int bisection_iter = 0; bisection_iter < num_bisection_iter; ++bisection_iter) {
            double t_mid = 0.5*(t_low + t_high);
            delta_t = t_mid - t;
            q_trj_per_target[target_idx].attr("__call__")(t_mid, std::ref(next_pos));
            next_pos = next_pos + next_rel_pos;
            dist = (next_pos - pos).norm();
            if (dist > vmax*delta_t) {
              // Travel is infeasible
              t_low = t_mid;
            } else {
              // Travel is feasible
              t_high = t_mid;
              t_high_dist = dist;
            }
          }

          delta_t = t_high - t;
          Xnew(seq_idx, 2) = delta_t;

          // Update cost, time, and position
          cost += t_high_dist;
          t = t_high;
          pos = next_pos;
        }
      }
      if (repair_failed) {
        continue;
      }

      if (cost < population_costs[chromosome_idx]) {
        population[chromosome_idx] = Xnew;
        population_costs[chromosome_idx] = cost;
      }

      // TODO: transformation to reduce cost (only for Dubins)
    }

    ++gen_idx;
    /*
    if (gen_idx%Tlp == 0) {
      for (int j = 0; j < gen_idx/Tlp; ++j) {
        // TODO: Get individual from top 50% and run local search
      }
    }
    */
  }

  auto it = std::min_element(population_costs.begin(), population_costs.end());
  int min_idx = it - population_costs.begin();
  const Ref<const MatrixXd> &Xbest = population[min_idx];

  int dim_q = dubins ? 3 : 2;
  if (dubins) {
    throw std::runtime_error("Did not implement Dubins yet");
  } else {
    double t = 0;
    Vector2d pos;
    for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
      int target_idx = Xbest(seq_idx, 0);
      double theta = Xbest(seq_idx, 1);
      double delta_t = Xbest(seq_idx, 2);
      t += delta_t;
      selected_pts_per_target(target_idx, 0) = t;
      q_trj_per_target[target_idx].attr("__call__")(t, std::ref(pos));
      pos = pos + target_radii[target_idx]*Vector2d(cos(theta), sin(theta));
      selected_pts_per_target(target_idx, 1) = pos(0);
      selected_pts_per_target(target_idx, 2) = pos(1);

      if (t < tw_per_target(target_idx, 0) || t > tw_per_target(target_idx, 1)) {
        throw std::runtime_error("error");
      }
    }
  }

  return population_costs[min_idx];
}
