#pragma once
#include <chrono>
#include <omp.h>
#include <Eigen/Dense>
#include <random>
#include <set>
#include "mt_tsp_ros2/reopt_gtsp_tour.h"
#include "mt_tsp_ros2/memetic_dubins_ce_mt_tsp/memetic_alg_params.h"
#include <iomanip>
#include <fstream>
#include <queue>

using namespace Eigen;

typedef Matrix<double, Dynamic, Dynamic, RowMajor> RowMatrixXd;

typedef Matrix<bool, Dynamic, 1> VectorXb;

typedef Matrix<double, 1, 1> Vector1d;

const double root_finding_tol = 1e-2;

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

bool repair_chromosome(Ref<MatrixXd> X, const py::array_t<long> &gtsp_cost_mat_rounded_and_scaled_flat, const py::array_t<double> &gtsp_cost_mat_flat, const std::vector<py::array_t<long>> &target_to_pt_ptr, const py::array_t<double> &all_pts, double &cost, long inf_val, int num_nodes, Ref<VectorXl> tour) {
  VectorXl target_seq = X.col(0).cast<long>();
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

  int pt_dim = X.cols() - 1;

  for (int tour_idx = 1; tour_idx < tour.size() - 1; ++tour_idx) {
    int node_idx = tour(tour_idx);
    for (int i = 0; i < pt_dim; ++i) {
      X(tour_idx - 1, 1 + i) = all_pts_unchecked(node_idx, i);
    }
    cost += gtsp_cost_mat_flat_unchecked(prev_node_idx*num_nodes + node_idx);
    // std::cout << prev_node_idx << " " << node_idx << " " << num_nodes << " " << gtsp_cost_mat_flat_unchecked(prev_node_idx*num_nodes + node_idx) << " " << gtsp_cost_mat_rounded_and_scaled_flat_unchecked(prev_node_idx*num_nodes + node_idx) << std::endl;
    prev_node_idx = node_idx;
  }
  if (std::isinf(cost)) {
    throw std::runtime_error("reopt_gtsp_tour returned infinite cost solution but reported success");
  }
  return false;
}

RowMatrixXd memetic_alg_sampling(Ref<RowMatrixXd> selected_pts_per_target, const Ref<const RowMatrixXd> &initial_population, const Ref<const VectorXd> &initial_costs, double time_limit, int num_openmp_threads, const MemeticAlgParams &params, Ref<Matrix<long, 1, 1>> num_feas_final, std::vector<double> &cost_vs_iterations, int max_generations, std::vector<int> &best_target_seq_change_per_iteration, py::object sample_point_graph_generator, int num_random_points_per_target, int nonimproving_iterations_between_resample, long inf_val, Ref<VectorXd> timing_info) {
  std::vector<std::pair<double, double>> cost_vs_time;

  auto timer_start = std::chrono::high_resolution_clock::now();

  double min_cost_record_time = 0.;

  omp_set_num_threads(num_openmp_threads);

  int num_targets = selected_pts_per_target.rows();

  int pop_size = initial_costs.size();
  if (pop_size != initial_population.rows()/num_targets) {
    throw std::runtime_error("Population size does not match the number of provided cost values");
  }

  int pt_dim = selected_pts_per_target.cols();
  int gene_size = 1 + pt_dim;

  // Max because some chromosomes may have duplicate points
  int max_interception_points = 1 + (pop_size + num_random_points_per_target)*num_targets;
  py::array_t<double> gtsp_cost_mat_flat({max_interception_points*max_interception_points}, {sizeof(double)});
  py::array_t<long> gtsp_cost_mat_rounded_and_scaled_flat({max_interception_points*max_interception_points}, {sizeof(double)});;
  py::array_t<double> all_pts({max_interception_points, pt_dim}, {sizeof(double)*pt_dim, sizeof(double)});
  std::vector<py::array_t<long>> target_to_pt_ptr(num_targets);

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

  // Initialize population
  #pragma omp parallel for
  for (int i = 0; i < pop_size; ++i) {
    population1[i] = initial_population.block(num_targets*i, 0, num_targets, gene_size);
  }
  int min_cost_idx = -1;
  cost_vs_time.push_back(std::pair<double, double>(0., Map<VectorXd>(population_costs1.data(), population_costs1.size()).minCoeff(&min_cost_idx)));
  cost_vs_iterations.push_back(cost_vs_time.back().second);
  VectorXi best_target_seq = population1[min_cost_idx].col(0).cast<int>();
  int num_finite_cost = 0;
  for (auto cost : population_costs1) {
    if (std::isfinite(cost)){ 
      ++num_finite_cost;
    }
  }
  std::cout << "Initialized population of " << num_finite_cost << " finite-cost solutions" << std::endl;

  std::vector<MatrixXd> population2 = population1;
  std::vector<double> population_costs2 = population_costs1;

  int gen_idx = 0;

  std::vector<MatrixXd> *population = &population1;
  std::vector<double> *population_costs = &population_costs1;

  std::vector<MatrixXd> *updated_population = &population1;
  std::vector<double> *updated_population_costs = &population_costs1;

  int nonimproving_iterations_since_resample = 0;

  int num_nodes;

  double resampling_time = 0.;
  double crossover_mutation_repair_time = 0.;
  if (timing_info.size() != 2) {
    throw std::runtime_error("timing_info should be size 2");
  }

  std::vector<py::array_t<long>> tour_per_chromosome(pop_size);

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

    if (gen_idx == 0 || nonimproving_iterations_since_resample == nonimproving_iterations_between_resample) {
      // std::cout << "resampling" << std::endl;
      auto timer_start = std::chrono::high_resolution_clock::now();
      num_nodes = sample_point_graph_generator.attr("__call__")(std::ref(*population), 
                                                                gtsp_cost_mat_flat, 
                                                                gtsp_cost_mat_rounded_and_scaled_flat,
                                                                all_pts,
                                                                std::ref(target_to_pt_ptr), 
                                                                num_random_points_per_target,
                                                                std::ref(tour_per_chromosome)).cast<int>();

      auto timer_stop = std::chrono::high_resolution_clock::now();
      auto nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
      resampling_time += ((double)nanos)/1e9;
      /*
      std::cout << "flat cost mat" << std::endl;
      auto gtsp_cost_mat_flat_unchecked = gtsp_cost_mat_flat.unchecked<1>();
      for (int i = 0; i < num_nodes*num_nodes; ++i) {
        std::cout << gtsp_cost_mat_flat_unchecked(i) << " ";
      }
      std::cout << std::endl;
      std::cout << std::endl;

      std::cout << "rounded flat cost mat" << std::endl;
      auto gtsp_cost_mat_rounded_and_scaled_flat_unchecked = gtsp_cost_mat_rounded_and_scaled_flat.unchecked<1>();
      for (int i = 0; i < num_nodes*num_nodes; ++i) {
        std::cout << gtsp_cost_mat_rounded_and_scaled_flat_unchecked(i) << " ";
      }
      std::cout << std::endl;
      std::cout << std::endl;

      std::cout << "all pts" << std::endl;
      auto all_pts_unchecked = all_pts.unchecked<2>();
      for (int row = 0; row < num_nodes; ++row) {
        for (int i = 0; i < pt_dim; ++i) {
          std::cout << all_pts_unchecked(row, i) << " ";
        }
        std::cout << std::endl;
        std::cout << std::endl;
      }
      std::cout << "target to pt ptr" << std::endl;
      for (int target_idx = 0; target_idx < num_targets; ++target_idx) {
        auto ptr = target_to_pt_ptr[target_idx].unchecked<1>();
        std::cout << ptr.size() << std::endl;
        for (int ptr_idx = 0; ptr_idx < ptr.size(); ++ptr_idx) {
          std::cout << ptr(ptr_idx) << " ";
        }
        std::cout << std::endl << std::endl;
      }

      std::cout << "there" << std::endl;
      throw std::runtime_error("Done");
      */
      nonimproving_iterations_since_resample = 0;
    }

    // std::cout << "Beginning generation " << gen_idx << std::endl;
  
    auto crossover_mutation_repair_timer_start = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for
    for (int chromosome_idx = 0; chromosome_idx < pop_size; ++chromosome_idx) {
      int parent1_idx = chromosome_idx;
      int parent2_idx = parent_distribution(rngs_per_thread[omp_get_thread_num()]);
      // Crossover
      VectorXb inserted_targets = VectorXb::Zero(num_targets);
      MatrixXd Xnew(num_targets, gene_size);
      int parent1_counter = 0;
      int parent2_counter = 0;

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
      if (mutation_sample < params.mutation_prob) {
        // Mutation
        int seq_idx1 = mutation_operator1_distribution(rngs_per_thread[omp_get_thread_num()]);
        int seq_idx2 = mutation_operator1_distribution(rngs_per_thread[omp_get_thread_num()]);
        RowVectorXd tmp = Xnew.row(seq_idx1);
        Xnew.row(seq_idx1) = Xnew.row(seq_idx2);
        Xnew.row(seq_idx2) = tmp;
      }
      // Xnew = (*population)[chromosome_idx];

      // Repair to restore feasibility
      double cost = 0.;

      auto timer_stop = std::chrono::high_resolution_clock::now();
      auto nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
      double repair_time_limit = time_limit - ((double)nanos)/1e9 ;

      VectorXl tour(num_targets + 2);
      bool repair_failed = repair_chromosome(Xnew, gtsp_cost_mat_rounded_and_scaled_flat, gtsp_cost_mat_flat, target_to_pt_ptr, all_pts, cost, inf_val, num_nodes, tour);
      if (repair_failed || cost >= (*population_costs)[chromosome_idx]) {
        (*updated_population)[chromosome_idx] = (*population)[chromosome_idx];
        (*updated_population_costs)[chromosome_idx] = (*population_costs)[chromosome_idx];
        continue;
      }

      auto chromosome_tour_unchecked = tour_per_chromosome[chromosome_idx].mutable_unchecked<1>();
      for (int tour_idx = 1; tour_idx < tour.size() - 1; ++tour_idx) {
        chromosome_tour_unchecked(tour_idx - 1) = tour(tour_idx);
      }

      // std::cout << "repair improved cost" << std::endl;

      (*updated_population)[chromosome_idx] = Xnew;
      (*updated_population_costs)[chromosome_idx] = cost;
    }
    
    auto crossover_mutation_repair_timer_stop = std::chrono::high_resolution_clock::now();
    auto crossover_mutation_repair_nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(crossover_mutation_repair_timer_stop - crossover_mutation_repair_timer_start).count();
    crossover_mutation_repair_time += ((double)crossover_mutation_repair_nanos)/1e9;

    ++gen_idx;

    auto timer_start2 = std::chrono::high_resolution_clock::now();
    auto it2 = std::min_element((*updated_population_costs).begin(), (*updated_population_costs).end());

    int best_target_seq_change = (best_target_seq -  (*updated_population)[it2 - updated_population_costs->begin()].col(0).cast<int>()).cast<bool>().cast<int>().sum();
    best_target_seq_change_per_iteration.push_back(best_target_seq_change);
    best_target_seq = (*updated_population)[it2 - updated_population_costs->begin()].col(0).cast<int>();

    auto timer_stop2 = std::chrono::high_resolution_clock::now();
    auto nanos2 = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop2 - timer_start2).count();
    min_cost_record_time += ((double)nanos2)/1e9;

    timer_stop = std::chrono::high_resolution_clock::now();
    nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
    if (cost_vs_time.size() && *it2 > cost_vs_time.back().second) {
      throw std::runtime_error("Cost increased after memetic alg iteration");
    }
    if (cost_vs_time.size() == 0 || *it2 < cost_vs_time.back().second) {
      std::cout << "New best cost: " << *it2 << std::endl;
      cost_vs_time.push_back(std::pair<double, double>(((double)nanos)/1e9, *it2));
    } else {
      nonimproving_iterations_since_resample += 1;
    }
    cost_vs_iterations.push_back(*it2);

    if (max_generations != -1 && gen_idx >= max_generations) {
      break;
    }
  } // Overall loop

  num_feas_final(0) = 0;
  for (auto cost : (*updated_population_costs)) {
    if (std::isfinite(cost)){ 
      ++num_feas_final(0);
    }
  }
  std::cout << "Final population contains " << num_feas_final(0) << " feas solns" << std::endl;

  auto it = std::min_element((*updated_population_costs).begin(), (*updated_population_costs).end());
  int min_idx = it - (*updated_population_costs).begin();
  const Ref<const MatrixXd> &Xbest = (*updated_population)[min_idx];

  if (std::isinf(*it)) {
  } else {
    for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
      int target_idx = Xbest(seq_idx, 0);

      selected_pts_per_target.row(target_idx) = Xbest.block(seq_idx, 1, 1, gene_size - 1);
    }
  }

  RowMatrixXd cost_vs_time_mat(cost_vs_time.size(), 2);
  for (int i = 0; i < cost_vs_time.size(); ++i) {
    cost_vs_time_mat(i, 0) = cost_vs_time[i].first;
    cost_vs_time_mat(i, 1) = cost_vs_time[i].second;
  }

  std::cout << "Spent " << min_cost_record_time << " s tracking what the min cost was after each iteration (just making sure this is not too large)" << std::endl;

  timing_info(0) = resampling_time;
  timing_info(1) = crossover_mutation_repair_time;

  return cost_vs_time_mat;
}
