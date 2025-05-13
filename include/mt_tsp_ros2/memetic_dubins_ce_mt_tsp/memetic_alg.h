#pragma once
#include <chrono>
#include <omp.h>
#include <Eigen/Dense>
#include <random>
#include <set>
#include "mt_tsp_ros2/cpp_spline.h"
#include <iomanip>

using namespace Eigen;

typedef Matrix<double, Dynamic, Dynamic, RowMajor> RowMatrixXd;

typedef Matrix<bool, Dynamic, 1> VectorXb;

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

bool repair_chromosome(Ref<MatrixXd> X, const Ref<const RowMatrixXd> &tw_per_target, const Ref<const VectorXd> &target_radii, const std::vector<CppSpline> &q_trj_per_target, const Ref<const Vector2d> &p0, double vmax, double &cost, bool dubins) {
  int num_targets = tw_per_target.rows();

  bool repair_failed = false;
  double t = 0;
  Vector2d pos = p0;
  Vector2d next_rel_pos;
  Vector2d next_pos;
  cost = 0;
  for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
    int target_idx = X(seq_idx, 0);
    double theta = X(seq_idx, 1);
    double delta_t = X(seq_idx, 2);
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
        X(seq_idx, 2) = delta_t;
      } else if (next_t < tw_per_target(target_idx, 0)) {
        next_t = tw_per_target(target_idx, 0);
        delta_t = next_t - t;
        X(seq_idx, 2) = delta_t;
      }

      // Check if interception is feasible
      next_rel_pos = target_radii[target_idx]*Vector2d(cos(theta), sin(theta));
      next_pos = q_trj_per_target[target_idx](next_t) + next_rel_pos;
      double dist = (next_pos - pos).norm();
      if (dist <= vmax*delta_t) {
        // Travel is feasible, no need to repair
        cost += dist;
        t = next_t;
        pos = next_pos;
        continue;
      }

      // Check if travel is feasible to next_rel_pos at end of time window
      next_pos = q_trj_per_target[target_idx](tw_per_target(target_idx, 1)) + next_rel_pos;
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
        next_pos = q_trj_per_target[target_idx](t_mid) + next_rel_pos;
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
      X(seq_idx, 2) = delta_t;

      // Update cost, time, and position
      cost += t_high_dist;
      t = t_high;
      pos = next_pos;
    }
  }
  if (repair_failed) {
    cost = std::numeric_limits<double>::infinity();
  }
  return repair_failed;
}

void check_chromosome_feasible(const Ref<const MatrixXd> &chromosome, const Ref<const RowMatrixXd> &tw_per_target) {
  int num_targets = tw_per_target.rows();
  double t = 0.;
  for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
    int target_idx = chromosome(seq_idx, 0);
    double theta = chromosome(seq_idx, 1);
    double delta_t = chromosome(seq_idx, 2);
    t += delta_t;

    if (t < tw_per_target(target_idx, 0) || t > tw_per_target(target_idx, 1)) {
      throw std::runtime_error("chromosome infeasible");
    }
  }
}

// initial_population should have number of rows equal to num_targets*pop_size, and gene_size columns
RowMatrixXd memetic_alg(Ref<RowMatrixXd> selected_pts_per_target, const Ref<const RowMatrixXd> &initial_population, const Ref<const VectorXd> &initial_costs, const std::vector<CppSpline> &q_trj_per_target_python, const Ref<const RowMatrixXd> &tw_per_target, const Ref<const VectorXd> &target_radii, double time_limit, int pop_size, double rho, double vmax, const Ref<const Vector2d> &p0, double theta0, int num_openmp_threads) {
  std::vector<std::pair<double, double>> cost_vs_time;
  cost_vs_time.push_back(std::pair<double, double>(0., initial_costs.minCoeff()));

  auto timer_start = std::chrono::high_resolution_clock::now();

  double min_cost_record_time = 0.;

  omp_set_num_threads(num_openmp_threads);

  // I'm doing this because I'm worried about the GIL
  std::vector<CppSpline> q_trj_per_target;
  for (auto q_trj : q_trj_per_target_python) {
    q_trj_per_target.push_back(CppSpline(q_trj));
  }

  bool dubins = rho != 0;

  int num_targets = tw_per_target.rows();

  std::vector<MatrixXd> population1(pop_size);
  std::vector<double> population_costs1(pop_size);
  Map<VectorXd>(population_costs1.data(), pop_size) = initial_costs;

  std::vector<unsigned int> seeds{4223100027, 870236586, 1683737518, 3430707182, 1613429085, 1714341085, 853110547, 1988005940, 2629786018, 1139192408};
  std::vector<std::mt19937> rngs_per_thread;
  if (num_openmp_threads > seeds.size()) {
    throw std::runtime_error("Too many threads, not enough stored random seeds");
  }
  for (int thread_idx = 0; thread_idx < num_openmp_threads; ++thread_idx) {
    rngs_per_thread.push_back(std::mt19937(seeds[thread_idx]));
  }
  std::uniform_int_distribution<int> parent_distribution(0, pop_size - 1);
  std::uniform_int_distribution<int> crossover_distribution(0, 1);
  std::uniform_real_distribution<double> mutation_distribution(0, 1);

  std::uniform_int_distribution<int> mutation_operator1_distribution(0, num_targets - 1);

  std::uniform_int_distribution<int> mutation_operator2_seq_idx_distribution(0, num_targets - 1);
  std::uniform_real_distribution<double> mutation_operator2_angle_distribution(0, 2*M_PI);

  std::uniform_int_distribution<int> mutation_operator3_seq_idx_distribution(0, num_targets - 1);
  std::uniform_real_distribution<double> mutation_operator3_delta_t_distribution(0, 1);

  std::uniform_int_distribution<int> local_search_elite_distribution(0, pop_size/2);
  std::uniform_int_distribution<int> local_search_gene_idx_distribution(0, num_targets - 1);
  std::uniform_int_distribution<int> local_search_grad_vs_sampling_distribution(0, 1);

  std::uniform_real_distribution<double> local_search_sampling_theta_distribution(0, 2*M_PI);

  int Tlp = 2; // From paper

  const int gene_size = 3; // target index, theta, and delta t

  const double mutation_prob = 0.1; // From paper

  // Initialize population
  for (int i = 0; i < pop_size; ++i) {
    population1[i] = initial_population.block(num_targets*i, 0, num_targets, gene_size);
  }

  std::vector<MatrixXd> population2 = population1;
  std::vector<double> population_costs2 = population_costs1;

  int gen_idx = 0;

  std::vector<MatrixXd> *population = &population1;
  std::vector<double> *population_costs = &population_costs1;

  std::vector<MatrixXd> *updated_population = &population1;
  std::vector<double> *updated_population_costs = &population_costs1;

  while (true) {
    auto timer_stop = std::chrono::high_resolution_clock::now();
    auto nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
    if (((double)nanos)/1e9 > time_limit) {
      break;
    }

    if (gen_idx%2) {
      population = &population2;
      population_costs = &population_costs2;

      updated_population = &population1;
      updated_population_costs = &population_costs1;
    } else {
      population = &population1;
      population_costs = &population_costs1;

      updated_population = &population2;
      updated_population_costs = &population_costs2;
    }

    #pragma omp parallel for
    for (int chromosome_idx = 0; chromosome_idx < pop_size; ++chromosome_idx) {
      int parent1_idx = chromosome_idx;
      int parent2_idx = parent_distribution(rngs_per_thread[omp_get_thread_num()]);
      // Crossover
      VectorXb inserted_targets = VectorXb::Zero(num_targets);
      MatrixXd Xnew(num_targets, gene_size);
      int parent1_counter = 0;
      int parent2_counter = 0;

      // check_chromosome_feasible((*population)[chromosome_idx], tw_per_target);

      for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
        if (crossover_distribution(rngs_per_thread[omp_get_thread_num()]) == 0) {
          while (inserted_targets((int)((*population)[parent1_idx](parent1_counter, 0)))) {
            ++parent1_counter;
          }
          if (parent1_counter >= num_targets) {
            throw std::runtime_error("parent 1 out of bounds");
          }
          int target_idx = (*population)[parent1_idx](parent1_counter, 0);
          Xnew.row(seq_idx) = (*population)[parent1_idx].row(parent1_counter);
          inserted_targets(target_idx) = true;
        } else {
          while (inserted_targets((int)((*population)[parent2_idx](parent2_counter, 0)))) {
            ++parent2_counter;
          }
          if (parent2_counter >= num_targets) {
            throw std::runtime_error("parent 2 out of bounds");
          }
          int target_idx = (*population)[parent2_idx](parent2_counter, 0);
          Xnew.row(seq_idx) = (*population)[parent2_idx].row(parent2_counter);
          inserted_targets(target_idx) = true;
        }
      }

      double mutation_sample = mutation_distribution(rngs_per_thread[omp_get_thread_num()]);
      if (mutation_sample < mutation_prob) {
        // Mutation
        if (mutation_sample < mutation_prob/3) {
          int seq_idx1 = mutation_operator1_distribution(rngs_per_thread[omp_get_thread_num()]);
          int seq_idx2 = mutation_operator1_distribution(rngs_per_thread[omp_get_thread_num()]);
          RowVectorXd tmp = Xnew.row(seq_idx1);
          Xnew.row(seq_idx1) = Xnew.row(seq_idx2);
          Xnew.row(seq_idx2) = tmp;
        } else if (mutation_sample < 2*mutation_prob/3) {
          int seq_idx = mutation_operator2_seq_idx_distribution(rngs_per_thread[omp_get_thread_num()]);
          double theta = mutation_operator2_angle_distribution(rngs_per_thread[omp_get_thread_num()]);
          Xnew(seq_idx, 1) = theta;
        } else {
          int seq_idx = mutation_operator3_seq_idx_distribution(rngs_per_thread[omp_get_thread_num()]);
          double raw_sample = mutation_operator3_delta_t_distribution(rngs_per_thread[omp_get_thread_num()]);
          int target_idx = Xnew(seq_idx, 0);
          double min_t = tw_per_target(target_idx, 0);
          double max_t = tw_per_target(target_idx, 1);
          double delta_t = (max_t - min_t)*raw_sample;
          Xnew(seq_idx, 2) = delta_t;
        }
      }

      // Repair to restore feasibility
      double cost = 0.;
      bool repair_failed = repair_chromosome(Xnew, tw_per_target, target_radii, q_trj_per_target, p0, vmax, cost, dubins);
      if (repair_failed) {
        (*updated_population)[chromosome_idx] = (*population)[chromosome_idx];
        (*updated_population_costs)[chromosome_idx] = (*population_costs)[chromosome_idx];
        continue;
      }

      if (cost < (*population_costs)[chromosome_idx]) {
        (*updated_population)[chromosome_idx] = Xnew;
        (*updated_population_costs)[chromosome_idx] = cost;
        // check_chromosome_feasible((*updated_population)[chromosome_idx], tw_per_target);
      }

      // TODO: transformation to reduce cost (only for Dubins)
    }

    /*
    for (int chromosome_idx = 0; chromosome_idx < pop_size; ++chromosome_idx) {
      check_chromosome_feasible((*updated_population)[chromosome_idx], tw_per_target);
    }
    */
    
    ++gen_idx;
    if (gen_idx%Tlp == 0) {
      std::vector<size_t> sort_idx = sort_indexes(*updated_population_costs);
      for (int j = 0; j < gen_idx/Tlp; ++j) {
        // Get individual from top 50% and run local search
        int chromosome_idx = sort_idx[local_search_elite_distribution(rngs_per_thread[omp_get_thread_num()])];
        int gene_idx = local_search_gene_idx_distribution(rngs_per_thread[omp_get_thread_num()]);
        double theta = (*updated_population)[chromosome_idx](gene_idx, 1);
        if (local_search_grad_vs_sampling_distribution(rngs_per_thread[omp_get_thread_num()]) == 0) {
          // Gradient-based local search
          bool improvement = true;
          while (improvement) {
            int max_backtrack_fd = 3;
            double delta_theta = 0.01; // To compute approximate gradient
            double new_cost = std::numeric_limits<double>::infinity();
            for (int backtrack_iter = 0; backtrack_iter < max_backtrack_fd; ++backtrack_iter) {
              double new_theta = theta + delta_theta;

              MatrixXd local_modification = (*updated_population)[chromosome_idx];
              local_modification(gene_idx, 1) = new_theta;

              bool repair_failed = repair_chromosome(local_modification, tw_per_target, target_radii, q_trj_per_target, p0, vmax, new_cost, dubins);
              if (!repair_failed) {
                break;
              }
              delta_theta *= 0.1;
            }
            if (std::isinf(new_cost)) {
              break;
            }
            double gradient = (new_cost - (*updated_population_costs)[chromosome_idx])/delta_theta;

            int max_backtrack_gd = 3;
            double gd_step_size = 0.01;
            for (int backtrack_iter = 0; backtrack_iter < max_backtrack_gd; ++backtrack_iter) {
              double new_theta = theta - gd_step_size*gradient;

              MatrixXd local_modification = (*updated_population)[chromosome_idx];
              local_modification(gene_idx, 1) = new_theta;

              bool repair_failed = repair_chromosome(local_modification, tw_per_target, target_radii, q_trj_per_target, p0, vmax, new_cost, dubins);
              if (!repair_failed) {
                if (new_cost < (*updated_population_costs)[chromosome_idx]) {
                  (*updated_population)[chromosome_idx] = local_modification;
                  (*updated_population_costs)[chromosome_idx] = new_cost;
                }

                break;
              }
              gd_step_size *= 0.1;
            }

            improvement = new_cost < (*updated_population_costs)[chromosome_idx] - 1e-4;
          }
        } else {
          if (dubins) {
            throw std::runtime_error("Did not implement Dubins yet");
          } else {
            int target_idx;
            double t = 0;
            double theta;
            double delta_t;
            double prev_t = 0.;
            for (int seq_idx = 0; seq_idx <= gene_idx; ++seq_idx) {
              target_idx = (*updated_population)[chromosome_idx](seq_idx, 0);
              theta = (*updated_population)[chromosome_idx](seq_idx, 1);
              delta_t = (*updated_population)[chromosome_idx](seq_idx, 2);
              prev_t = t;
              t += delta_t;
            }
            // check_chromosome_feasible((*updated_population)[chromosome_idx], tw_per_target);

            // Do sampling-based local search
            const int num_samples = 20; // From paper
            std::vector<MatrixXd> local_modifications(num_samples);
            std::vector<double> local_modification_costs(num_samples);
            #pragma omp parallel for
            for (int sample_idx = 0; sample_idx < num_samples; ++sample_idx) {
              // double new_theta = local_search_sampling_theta_distribution(rngs_per_thread[omp_get_thread_num()]);
              double new_theta = (2*M_PI*sample_idx)/num_samples;

              double new_cost;
              local_modifications[sample_idx] = (*updated_population)[chromosome_idx];
              local_modifications[sample_idx](gene_idx, 1) = new_theta;

              bool repair_failed = repair_chromosome(local_modifications[sample_idx], tw_per_target, target_radii, q_trj_per_target, p0, vmax, new_cost, dubins);

              if (!(repair_failed || new_cost >= (*updated_population_costs)[chromosome_idx])) {
                local_modification_costs[sample_idx] = new_cost;
                // check_chromosome_feasible(local_modifications[sample_idx], tw_per_target);
              } else {
                local_modification_costs[sample_idx] = std::numeric_limits<double>::infinity();
              }
            } // Sampling-based local search iterations

            auto it = std::min_element(local_modification_costs.begin(), local_modification_costs.end());
            int min_idx = it - local_modification_costs.begin();
            if (local_modification_costs[min_idx] < (*updated_population_costs)[chromosome_idx]) {
              (*updated_population)[chromosome_idx] = local_modifications[min_idx];
              (*updated_population_costs)[chromosome_idx] = local_modification_costs[min_idx];
            }
          } // Not Dubins
        } // Sampling-based local search instead of gradient
      } // Pick a chromosome for local search
    } // Check if local search condition has been met

    auto timer_start2 = std::chrono::high_resolution_clock::now();
    auto it2 = std::min_element((*updated_population_costs).begin(), (*updated_population_costs).end());
    auto timer_stop2 = std::chrono::high_resolution_clock::now();
    auto nanos2 = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop2 - timer_start2).count();
    min_cost_record_time += ((double)nanos2)/1e9;

    timer_stop = std::chrono::high_resolution_clock::now();
    nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
    if (cost_vs_time.size() && *it2 > cost_vs_time.back().second) {
      throw std::runtime_error("Cost increased after memetic alg iteration");
    }
    cost_vs_time.push_back(std::pair<double, double>(((double)nanos)/1e9, *it2));
  } // Overall loop

  auto it = std::min_element((*updated_population_costs).begin(), (*updated_population_costs).end());
  int min_idx = it - (*updated_population_costs).begin();
  const Ref<const MatrixXd> &Xbest = (*updated_population)[min_idx];

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
      pos = q_trj_per_target[target_idx](t) + target_radii[target_idx]*Vector2d(cos(theta), sin(theta));
      selected_pts_per_target(target_idx, 1) = pos(0);
      selected_pts_per_target(target_idx, 2) = pos(1);

      if (t < tw_per_target(target_idx, 0) || t > tw_per_target(target_idx, 1)) {
        throw std::runtime_error("t out of window");
      }
    }
  }

  RowMatrixXd cost_vs_time_mat(cost_vs_time.size(), 2);
  for (int i = 0; i < cost_vs_time.size(); ++i) {
    cost_vs_time_mat(i, 0) = cost_vs_time[i].first;
    cost_vs_time_mat(i, 1) = cost_vs_time[i].second;
  }

  std::cout << "Spent " << min_cost_record_time << " s tracking what the min cost was after each iteration (just making sure this is not too large)" << std::endl;
  return cost_vs_time_mat;
}
