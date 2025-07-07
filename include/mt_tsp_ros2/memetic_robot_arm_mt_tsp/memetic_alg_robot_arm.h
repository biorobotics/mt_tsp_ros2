#pragma once
#include <chrono>
#include <omp.h>
#include <Eigen/Dense>
#include <random>
#include <set>
#include "mt_tsp_ros2/SE3_spline.h"
#include "mt_tsp_ros2/memetic_dubins_ce_mt_tsp/memetic_alg_params.h"
#include <iomanip>
#include "mt_tsp_ros2/memetic_robot_arm_mt_tsp/kuka_ik.h"
#include "mt_tsp_ros2/memetic_robot_arm_mt_tsp/robot_arm_nlp.h"

using namespace Eigen;

typedef Matrix<double, Dynamic, Dynamic, RowMajor> RowMatrixXd;

typedef Matrix<bool, Dynamic, 1> VectorXb;

typedef Matrix<double, 1, 1> Vector1d;

const int gene_size = 2 + dim_q; // target index, delta t, and configuration

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

// Return true if repair failed
bool repair_chromosome(Ref<MatrixXd> X, double &cost, std::shared_ptr<RobotArmNLPInfo> nlp_info, Ref<RowMatrixXd> trajectory, Ref<VectorXd> warm_start, int ipopt_print_level, Ref<VectorXb> restoration_info) {
  VectorXl target_seq = X.col(0).cast<long>();
  int num_targets = target_seq.size();
  RobotArmNLPSolver nlp_solver(nlp_info, target_seq, 500, ipopt_print_level); // TODO: tune max number of iterations

  double t = 0.;
  for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
    int target_idx = X(seq_idx, 0);
    double delta_t = X(seq_idx, 1);
    t += delta_t;
    if (t < nlp_info->tw_per_target(target_idx, 0)) {
      t = nlp_info->tw_per_target(target_idx, 0);
    } else if (t > nlp_info->tw_per_target(target_idx, 1)) {
      t = nlp_info->tw_per_target(target_idx, 1);
    }
    warm_start(seq_idx*(1 + dim_q)) = t;
    warm_start.segment(seq_idx*(1 + dim_q) + 1, dim_q) = X.block(seq_idx, 2, 1, dim_q).transpose();
  }
  nlp_solver.set_warm_start(warm_start);
  bool restoration_invoked;
  bool restoration_failed;
  double new_cost = nlp_solver.solve(trajectory, true, restoration_invoked, restoration_failed);
  if (restoration_info.size() == 2) {
    restoration_info(0) = restoration_invoked;
    restoration_info(1) = restoration_failed;
  } else {
    throw std::runtime_error("Restoration info not the correct size");
  }
  if (std::isinf(new_cost)) {
    return true;
  }
  if (new_cost < cost) {
    cost = new_cost;
    X.block(0, 2, num_targets, dim_q) = trajectory.rightCols(dim_q);
    t = 0.;
    for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
      double delta_t = trajectory(seq_idx, 0) - t;
      X(seq_idx, 1) = delta_t;
      t = trajectory(seq_idx, 0);
    }
  }
  return false;
}

bool repair_chromosome(Ref<MatrixXd> X, Ref<Vector1d> cost, const Ref<const RowMatrixXd> &tw_per_target, const std::vector<SE3Spline> &q_trj_per_target, const Ref<const VectorXd> &q0, const Ref<const VectorXd> &joint_limits, const Ref<const VectorXd> &vmax, const Ref<const Vector4d> &l, const Ref<const RowMatrixXd> &dh, Ref<RowMatrixXd> trajectory, Ref<VectorXd> warm_start, int ipopt_print_level, Ref<VectorXb> restoration_info) {
  std::shared_ptr<RobotArmNLPInfo> nlp_info = std::make_shared<RobotArmNLPInfo>(tw_per_target, q_trj_per_target, q0, joint_limits, vmax, l, dh);
  double dummy_cost;
  bool repair_failed = repair_chromosome(X, dummy_cost, nlp_info, trajectory, warm_start, ipopt_print_level, restoration_info);
  cost(0) = dummy_cost;
  return repair_failed;
}

// initial_population should have number of rows equal to num_targets*pop_size, and gene_size columns
RowMatrixXd memetic_alg(Ref<RowMatrixXd> selected_pts_per_target, const Ref<const RowMatrixXd> &initial_population, const Ref<const VectorXd> &initial_costs, const std::vector<SE3Spline> &q_trj_per_target_python, const Ref<const RowMatrixXd> &tw_per_target, double time_limit, const Ref<const VectorXd> &vmax, const Ref<const VectorXd> &joint_limits, const Ref<const Vector4d> &l, const Ref<const RowMatrixXd> &dh, const Ref<const VectorXd> &q0, int num_openmp_threads, const MemeticAlgParams &params, Ref<Matrix<long, 1, 1>> num_feas_final, std::vector<double> &cost_vs_iterations, int max_generations, std::vector<int> &best_target_seq_change_per_iteration) {
  std::vector<std::pair<double, double>> cost_vs_time;
  
  auto timer_start = std::chrono::high_resolution_clock::now();

  double min_cost_record_time = 0.;

  omp_set_num_threads(num_openmp_threads);

  int num_targets = tw_per_target.rows();

  std::shared_ptr<RobotArmNLPInfo> nlp_info = std::make_shared<RobotArmNLPInfo>(tw_per_target, q_trj_per_target_python, q0, joint_limits, vmax, l, dh);

  int pop_size = initial_costs.size();
  if (pop_size != initial_population.rows()/num_targets) {
    throw std::runtime_error("Population size does not match the number of provided cost values");
  }

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
  std::uniform_int_distribution<int> mutation_operator2_qidx_distribution(0, dim_q - 1);
  std::uniform_real_distribution<double> mutation_operator2_joint_angle_distribution(0, 1);

  std::uniform_int_distribution<int> mutation_operator3_seq_idx_distribution(0, num_targets - 1);
  std::uniform_real_distribution<double> mutation_operator3_delta_t_distribution(0, 1);

  std::vector<RowMatrixXd> trajectory_per_thread;
  std::vector<VectorXd> warm_start_per_thread;
  for (int thread_idx = 0; thread_idx < num_openmp_threads; ++thread_idx) {
    trajectory_per_thread.push_back(RowMatrixXd(num_targets, 1 + dim_q));
    warm_start_per_thread.push_back(VectorXd(num_targets*(1 + dim_q)));
  }

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
        if (mutation_sample < params.mutation_prob/3) {
          // Swap two targets
          int seq_idx1 = mutation_operator1_distribution(rngs_per_thread[omp_get_thread_num()]);
          int seq_idx2 = mutation_operator1_distribution(rngs_per_thread[omp_get_thread_num()]);
          RowVectorXd tmp = Xnew.row(seq_idx1);
          Xnew.row(seq_idx1) = Xnew.row(seq_idx2);
          Xnew.row(seq_idx2) = tmp;
        } else if (mutation_sample < 2*params.mutation_prob/3) {
          // Perturb a joint angle
          int seq_idx = mutation_operator2_seq_idx_distribution(rngs_per_thread[omp_get_thread_num()]);
          int q_idx = mutation_operator2_qidx_distribution(rngs_per_thread[omp_get_thread_num()]);
          double raw_sample = mutation_operator2_joint_angle_distribution(rngs_per_thread[omp_get_thread_num()]);
          Xnew(seq_idx, 2 + q_idx) = -joint_limits(q_idx) + 2*joint_limits(q_idx)*raw_sample;
        } else {
          // Perturb delta_t
          int seq_idx = mutation_operator3_seq_idx_distribution(rngs_per_thread[omp_get_thread_num()]);
          double raw_sample = mutation_operator3_delta_t_distribution(rngs_per_thread[omp_get_thread_num()]);
          int target_idx = Xnew(seq_idx, 0);
          double delta_t;
          double min_t = tw_per_target(target_idx, 0);
          double max_t = tw_per_target(target_idx, 1);
          delta_t = (max_t - min_t)*raw_sample;
          Xnew(seq_idx, 1) = delta_t;
        }
      }

      // Repair to restore feasibility
      double cost = (*population_costs)[chromosome_idx];
      VectorXb restoration_info(2);
      bool repair_failed = repair_chromosome(Xnew, cost, nlp_info, trajectory_per_thread[omp_get_thread_num()], warm_start_per_thread[omp_get_thread_num()], 0, restoration_info);
      if (repair_failed) {
        (*updated_population)[chromosome_idx] = (*population)[chromosome_idx];
        (*updated_population_costs)[chromosome_idx] = (*population_costs)[chromosome_idx];
      } else {
        (*updated_population)[chromosome_idx] = Xnew;
        (*updated_population_costs)[chromosome_idx] = cost;
        std::cout << "improved cost of chromosome" << std::endl;
      }
    }
 
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

  if (std::isfinite(*it)) {
    double t = 0;
    VectorXd q = q0;
    for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
      int target_idx = Xbest(seq_idx, 0);
      VectorXd next_q = Xbest.block(seq_idx, 2, 1, dim_q).transpose();
      double delta_t = Xbest(seq_idx, 1);
      t += delta_t;
      selected_pts_per_target(target_idx, 0) = t;
      double travel_time = ((next_q - q).array()/vmax.array()).matrix().lpNorm<Infinity>();
      if (travel_time > delta_t + 1e-4) {
        throw std::runtime_error("Speed constraint violated when getting trajectory associated with best chromosome");
      }
      q = next_q;
      selected_pts_per_target.block(target_idx, 1, 1, dim_q) = q.transpose();

      if ((t < tw_per_target(target_idx, 0) - 1e-4 || t > tw_per_target(target_idx, 1) + 1e-4)) {
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
