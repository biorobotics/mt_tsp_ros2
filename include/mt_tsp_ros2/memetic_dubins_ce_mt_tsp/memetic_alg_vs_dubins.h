#pragma once
#include <chrono>
#include <omp.h>
#include <Eigen/Dense>
#include <random>
#include <set>
#include "mt_tsp_ros2/cpp_spline.h"
#include "mt_tsp_ros2/extended_cpp_spline.h"
#include "mt_tsp_ros2/elongate_one_sided_dubins_path.h"
#include "mt_tsp_ros2/memetic_dubins_ce_mt_tsp/memetic_alg_params.h"
#include <iomanip>
#include <fstream>
#include <queue>

using namespace Eigen;

typedef Matrix<double, Dynamic, Dynamic, RowMajor> RowMatrixXd;

typedef Matrix<bool, Dynamic, 1> VectorXb;

typedef Matrix<double, 1, 1> Vector1d;

const int gene_size = 2; // target index, and delta t

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

struct RepairTreeNode {
  double t;
  Vector2d pos;
  double heading;
  int seq_idx;
  std::shared_ptr<RepairTreeNode> parent;
  double g;
  double h;
  double f;

  RepairTreeNode(double t, const Ref<const Vector2d> &pos, double heading, int seq_idx, std::shared_ptr<RepairTreeNode> parent, double g, double h) : t(t), pos(pos), heading(heading), seq_idx(seq_idx), parent(parent), g(g), h(h), f(g + h) {
  }
};

typedef std::shared_ptr<RepairTreeNode> RepairTreeNodePtr;

struct compare_repair_tree_nodes {
  bool operator() (RepairTreeNodePtr &elem1,
                   RepairTreeNodePtr &elem2) {
    return elem1->f > elem2->f;
  }
};

void find_next_interception_point(int target_idx, const Ref<const RowMatrixXd> &tw_per_target, const std::vector<ExtendedCppSpline> &q_trj_per_target, double v, double rho, double t, const Ref<const Vector2d> &pos, double heading, double &next_t, Ref<Vector2d> next_pos, double &next_heading, const MemeticAlgParams& params, double &transition_cost) {
  int max_newton_iter = 1000;

  int num_targets = tw_per_target.rows();

  double feas_next_t = std::numeric_limits<double>::infinity();
  Vector2d feas_next_pos = std::numeric_limits<double>::infinity()*Vector2d::Ones();
  double feas_next_heading = std::numeric_limits<double>::infinity();

  // Check start of time window
  next_t = tw_per_target(target_idx, 0);
  next_pos = q_trj_per_target[target_idx](next_t);
  double delta_t = next_t - t;
  RowMatrixXd turns = elongated_dubins_path_one_sided(pos(0), pos(1), heading, next_pos(0), next_pos(1), v*delta_t, rho, root_finding_tol);
  if (std::isfinite(turns(0, 0))) {
    feas_next_t = next_t;
    feas_next_pos = next_pos;
    feas_next_heading = heading;
    for (int row = 0; row < turns.rows(); ++row) {
      if (turns(row, 0) != 0) {
        feas_next_heading += turns(row, 1)/turns(row, 0);
      }
    }
    max_newton_iter = 0; // No need to run Newton because we can get to the start of the time window
  } else {
    // Check end of time window. If infeasible, return, and if feasible, use as starting point for Newton
    next_t = tw_per_target(target_idx, 1);
    next_pos = q_trj_per_target[target_idx](next_t);
    delta_t = next_t - t;
    turns = elongated_dubins_path_one_sided(pos(0), pos(1), heading, next_pos(0), next_pos(1), v*delta_t, rho, root_finding_tol);
    if (std::isinf(turns(0, 0))) {
      transition_cost = std::numeric_limits<double>::infinity();
      return;
    }
    feas_next_t = next_t;
    feas_next_pos = next_pos;
    feas_next_heading = heading;
    for (int row = 0; row < turns.rows(); ++row) {
      if (turns(row, 0) != 0) {
        feas_next_heading += turns(row, 1)/turns(row, 0);
      }
    }
  }

  // Run Newton to restore feasibility if needed
  double delta_delta_t_finite_diff = 1e-4;
  // std::cout << "starting newton" << std::endl;
  for (int newton_iter = 0; newton_iter < max_newton_iter; ++newton_iter) {
    RowMatrixXd shortest_path_turns = turns_for_one_sided_dubins_path(pos(0), pos(1), heading, next_pos(0), next_pos(1), rho);
    double shortest_path_dist = shortest_path_turns.col(1).sum();
    double c = shortest_path_dist - v*delta_t;
    // delta_vs_iterations.push_back(c);

    if (std::abs(c) < root_finding_tol) {
      break;
    }

    // Finite-diff
    double delta_t_plus = delta_t + delta_delta_t_finite_diff;

    double next_t_plus = t + delta_t_plus;

    Vector2d next_pos_plus = q_trj_per_target[target_idx](next_t_plus);

    RowMatrixXd shortest_path_turns_plus = turns_for_one_sided_dubins_path(pos(0), pos(1), heading, next_pos_plus(0), next_pos_plus(1), rho);
    double shortest_path_dist_plus = shortest_path_turns_plus.col(1).sum();
    double c_plus = shortest_path_dist_plus - v*delta_t_plus;

    double deriv = (c_plus - c)/delta_delta_t_finite_diff;

    delta_t -= params.repair_step_size*c/deriv;

    next_t = t + delta_t;

    if (next_t > tw_per_target(target_idx, 1)) {
      next_t = tw_per_target(target_idx, 1);
      delta_t = next_t - t;
      if (delta_t < 0) {
        throw std::runtime_error("Projection of delta t failed in repair");
      }
    } else if (next_t < tw_per_target(target_idx, 0)) {
      next_t = tw_per_target(target_idx, 0);
      delta_t = next_t - t;
    }

    next_pos = q_trj_per_target[target_idx](next_t);

    RowMatrixXd turns = elongated_dubins_path_one_sided(pos(0), pos(1), heading, next_pos(0), next_pos(1), v*delta_t, rho, root_finding_tol);
    if (std::isfinite(turns(0, 0)) && next_t < feas_next_t) {
      feas_next_t = next_t;
      feas_next_pos = next_pos;
      feas_next_heading = heading;
      for (int row = 0; row < turns.rows(); ++row) {
        if (turns(row, 0) != 0) {
          feas_next_heading += turns(row, 1)/turns(row, 0);
        }
      }
    }
  }

  next_t = feas_next_t;
  next_pos = feas_next_pos;
  next_heading = feas_next_heading;
  delta_t = next_t - t;

  if (params.min_latency) {
    transition_cost = next_t - tw_per_target(target_idx, 0);
  } else {
    if (params.min_time) {
      transition_cost = next_t - t;
    } else {
      transition_cost = v*(next_t - t);
    }
  }
}

// bool repair_chromosome(Ref<MatrixXd> X, const Ref<const RowMatrixXd> &tw_per_target, const std::vector<ExtendedCppSpline> &q_trj_per_target, const Ref<const Vector2d> &p0, double heading0, const Ref<const VectorXd> &speed_options, double &cost, double wmax, const MemeticAlgParams &params, double t0, std::vector<double> &delta_vs_iterations, Ref<RowMatrixXd> selected_pts_per_target, bool populate_selected_pts) {
bool repair_chromosome(Ref<MatrixXd> X, const Ref<const RowMatrixXd> &tw_per_target, const std::vector<ExtendedCppSpline> &q_trj_per_target, const Ref<const Vector2d> &p0, double heading0, const Ref<const VectorXd> &speed_options, double &cost, double wmax, const MemeticAlgParams &params, double t0, double time_limit, Ref<RowMatrixXd> selected_pts_per_target, bool populate_selected_pts) {
  auto timer_start = std::chrono::high_resolution_clock::now();

  int num_targets = tw_per_target.rows();

  /*
  std::priority_queue<RepairTreeNodePtr, std::vector<RepairTreeNodePtr>, compare_repair_tree_nodes> open_list;
  open_list.push(std::make_shared<RepairTreeNode>(t0, p0, heading0, -1, nullptr, 0., 0.));
  */

  std::vector<RepairTreeNodePtr> stack;
  stack.push_back(std::make_shared<RepairTreeNode>(t0, p0, heading0, -1, nullptr, 0., 0.));

  int expansion_limit = num_targets*10;
  int num_expansions = 0;

  RepairTreeNodePtr goal = nullptr;
  // while (open_list.size()) {
  while (stack.size()) {
    auto timer_stop = std::chrono::high_resolution_clock::now();
    auto nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
    if (((double)nanos)/1e9 > time_limit || num_expansions >= expansion_limit) {
      break;
    }

    /*
    RepairTreeNodePtr pop = open_list.top();
    open_list.pop();
    */

    RepairTreeNodePtr pop = stack.back();
    stack.pop_back();

    if (pop->seq_idx == num_targets - 1) {
      goal = pop;
      break;
    }

    ++num_expansions;

    int next_target_idx = X(pop->seq_idx + 1, 0);
    double next_t;
    Vector2d next_pos;
    double next_heading;
    double transition_cost;

    /*
    for (int speed_idx = 0; speed_idx < speed_options.size(); ++speed_idx) {
      double v = speed_options(speed_idx);
      find_next_interception_point(next_target_idx, tw_per_target, q_trj_per_target, v, v/wmax, pop->t, pop->pos, pop->heading, next_t, next_pos, next_heading, params, transition_cost);

      if (std::isinf(transition_cost)) {
        continue;
      }

      RepairTreeNodePtr successor = std::make_shared<RepairTreeNode>(next_t, next_pos, next_heading, pop->seq_idx + 1, pop, pop->g + transition_cost, 0.);
      open_list.push(successor);
    }
    */

    std::vector<double> sort_vals_per_successor;
    std::vector<RepairTreeNodePtr> successors;

    for (int speed_idx = 0; speed_idx < speed_options.size(); ++speed_idx) {
      double v = speed_options(speed_idx);
      find_next_interception_point(next_target_idx, tw_per_target, q_trj_per_target, v, v/wmax, pop->t, pop->pos, pop->heading, next_t, next_pos, next_heading, params, transition_cost);

      if (std::isinf(transition_cost)) {
        continue;
      }

      RepairTreeNodePtr successor = std::make_shared<RepairTreeNode>(next_t, next_pos, next_heading, pop->seq_idx + 1, pop, pop->g + transition_cost, 0.);
      successors.push_back(successor);
      sort_vals_per_successor.push_back(transition_cost);
      // sort_vals_per_successor.push_back(next_t);
    }

    std::vector<size_t> sort_idx = sort_indexes(sort_vals_per_successor);
    std::reverse(sort_idx.begin(), sort_idx.end());

    for (int neighbor_idx : sort_idx) {
      stack.push_back(successors[neighbor_idx]);
    }
  }

  if (goal == nullptr) {
    cost = std::numeric_limits<double>::infinity();
    return true;
  }
  // std::cout << num_expansions << std::endl;

  RepairTreeNodePtr node = goal;
  cost = goal->g;
  for (int seq_idx = num_targets - 1; seq_idx >= 0; --seq_idx) {
    if (node->parent == nullptr) {
      throw std::runtime_error("Parent is null where it should not be null during backpointer traversal");
    }
    if (node->seq_idx != seq_idx) {
      throw std::runtime_error("Node seq_idx doesnt match seq_idx during backpointer traversal");
    }
    X(seq_idx, 1) = node->t - node->parent->t;

    if (populate_selected_pts) {
      int target_idx = X(seq_idx, 0);
      selected_pts_per_target(target_idx, 0) = node->t;
      selected_pts_per_target(target_idx, 1) = node->pos(0);
      selected_pts_per_target(target_idx, 2) = node->pos(1);
      selected_pts_per_target(target_idx, 3) = node->heading;
    }

    node = node->parent;
  }
  if (node->parent != nullptr) {
    throw std::runtime_error("Parent is not null where it should be null during backpointer traversal");
  }

  return false;
}

// bool repair_chromosome(Ref<MatrixXd> X, const Ref<const RowMatrixXd> &tw_per_target, const std::vector<ExtendedCppSpline> &q_trj_per_target, const Ref<const Vector2d> &p0, double heading0, const Ref<const VectorXd> &speed_options, Ref<Vector1d> cost, double wmax, const MemeticAlgParams &params, std::vector<double> &delta_vs_iterations, double time_limit) {
bool repair_chromosome(Ref<MatrixXd> X, const Ref<const RowMatrixXd> &tw_per_target, const std::vector<ExtendedCppSpline> &q_trj_per_target, const Ref<const Vector2d> &p0, double heading0, const Ref<const VectorXd> &speed_options, Ref<Vector1d> cost, double wmax, const MemeticAlgParams &params, double time_limit) {
  double tmp_cost = 0.;
  RowMatrixXd selected_pts_per_target(0, 0);
  // bool repair_failed = repair_chromosome(X, tw_per_target, q_trj_per_target, p0, heading0, speed_options, tmp_cost, wmax, params, t0, delta_vs_iterations, selected_pts_per_target, false);
  bool repair_failed = repair_chromosome(X, tw_per_target, q_trj_per_target, p0, heading0, speed_options, tmp_cost, wmax, params, 0., time_limit, selected_pts_per_target, false);
  cost(0) = tmp_cost;
  return repair_failed;
}

bool repair_chromosome(Ref<MatrixXd> X, const Ref<const RowMatrixXd> &tw_per_target, const std::vector<ExtendedCppSpline> &q_trj_per_target, const Ref<const Vector2d> &p0, double heading0, const Ref<const VectorXd> &speed_options, Ref<Vector1d> cost, double wmax, const MemeticAlgParams &params, double time_limit, Ref<RowMatrixXd> selected_pts_per_target) {
  double tmp_cost = 0.;
  bool repair_failed = repair_chromosome(X, tw_per_target, q_trj_per_target, p0, heading0, speed_options, tmp_cost, wmax, params, 0., time_limit, selected_pts_per_target, true);
  cost(0) = tmp_cost;
  return repair_failed;
}

// initial_population should have number of rows equal to num_targets*pop_size, and gene_size columns
RowMatrixXd memetic_alg(Ref<RowMatrixXd> selected_pts_per_target, const Ref<const RowMatrixXd> &initial_population, const Ref<const VectorXd> &initial_costs, const std::vector<ExtendedCppSpline> &q_trj_per_target_python, const Ref<const RowMatrixXd> &tw_per_target, double time_limit, double wmax, const Ref<const VectorXd> &speed_options, const Ref<const Vector2d> &p0, double heading0, int num_openmp_threads, const MemeticAlgParams &params, Ref<Matrix<long, 1, 1>> num_feas_final, std::vector<double> &cost_vs_iterations, int max_generations, std::vector<int> &best_target_seq_change_per_iteration) {
  std::vector<std::pair<double, double>> cost_vs_time;
  
  auto timer_start = std::chrono::high_resolution_clock::now();

  double min_cost_record_time = 0.;

  omp_set_num_threads(num_openmp_threads);

  int num_targets = tw_per_target.rows();

  // I'm doing this regardless of whether there are time windows because I'm worried about the GIL
  std::vector<ExtendedCppSpline> q_trj_per_target;
  for (int target_idx = 0; target_idx < num_targets; ++target_idx) {
    q_trj_per_target.push_back(ExtendedCppSpline(q_trj_per_target_python[target_idx].get_knots(), q_trj_per_target_python[target_idx].get_coeffs(), tw_per_target(target_idx, 0), tw_per_target(target_idx, 1)));
  }

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

  while (true) {
    auto timer_stop = std::chrono::high_resolution_clock::now();
    auto nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
    if (((double)nanos)/1e9 > time_limit) {
      break;
    }

    // std::cout << "Beginning generation " << gen_idx << std::endl;

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
        int seq_idx1 = mutation_operator1_distribution(rngs_per_thread[omp_get_thread_num()]);
        int seq_idx2 = mutation_operator1_distribution(rngs_per_thread[omp_get_thread_num()]);
        RowVectorXd tmp = Xnew.row(seq_idx1);
        Xnew.row(seq_idx1) = Xnew.row(seq_idx2);
        Xnew.row(seq_idx2) = tmp;
      }

      // Repair to restore feasibility
      double cost = 0.;

      auto timer_stop = std::chrono::high_resolution_clock::now();
      auto nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
      double repair_time_limit = time_limit - ((double)nanos)/1e9 ;

      bool repair_failed = repair_chromosome(Xnew, tw_per_target, q_trj_per_target, p0, heading0, speed_options, cost, wmax, params, 0., repair_time_limit, selected_pts_per_target, false);
      if (repair_failed || cost >= (*population_costs)[chromosome_idx]) {
        (*updated_population)[chromosome_idx] = (*population)[chromosome_idx];
        (*updated_population_costs)[chromosome_idx] = (*population_costs)[chromosome_idx];
        continue;
      }

      (*updated_population)[chromosome_idx] = Xnew;
      (*updated_population_costs)[chromosome_idx] = cost;
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
  Ref<MatrixXd> Xbest = (*updated_population)[min_idx];

  if (std::isinf(*it)) {
  } else {
    double cost;
    bool repair_failed = repair_chromosome(Xbest, tw_per_target, q_trj_per_target, p0, heading0, speed_options, cost, wmax, params, 0., 60., selected_pts_per_target, true);
    if (repair_failed) {
      throw std::runtime_error("Ran repair on best chromosome to get interception points, and repair failed, even though best chromosome has finite cost");
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
