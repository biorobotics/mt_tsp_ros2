#pragma once
#include <chrono>
#include <omp.h>
#include <Eigen/Dense>
#include <random>
#include <set>
#include "mt_tsp_ros2/cpp_spline.h"
#include "mt_tsp_ros2/elongate_one_sided_dubins_path.h"
#include <iomanip>

using namespace Eigen;

typedef Matrix<double, Dynamic, Dynamic, RowMajor> RowMatrixXd;

typedef Matrix<bool, Dynamic, 1> VectorXb;

typedef Matrix<double, 1, 1> Vector1d;

const int gene_size = 3; // target index, theta, and delta t

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

bool repair_chromosome(Ref<MatrixXd> X, const Ref<const RowMatrixXd> &tw_per_target, const Ref<const VectorXd> &target_radii, const std::vector<CppSpline> &q_trj_per_target, const Ref<const Vector2d> &p0, double heading0, double vmax, double &cost, bool dubins, double rho, int &max_newton_iter_for_success_repair) {
  int num_targets = tw_per_target.rows();

  bool repair_failed = false;
  double t = 0;
  Vector2d pos = p0;
  double heading = heading0;
  Vector2d next_rel_pos;
  Vector2d next_pos;
  cost = 0;
  for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
    int target_idx = X(seq_idx, 0);
    double theta = X(seq_idx, 1);
    double delta_t = X(seq_idx, 2);
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

    next_rel_pos = target_radii[target_idx]*Vector2d(cos(theta), sin(theta));

    // Check if interception is feasible
    if (dubins) {
      // Run Newton to restore feasibility if needed
      int max_newton_iter = 2776; // On one of the 10 target instances I ran, we needed at most 1388 iterations for successful transformation so I'm using 2x that number to declare failure
      double delta_delta_t_finite_diff = 1e-4;
      double repair_tol = 1e-4;
      bool newton_succeeded = false;
      // std::cout << "starting newton" << std::endl;
      for (int newton_iter = 0; newton_iter < max_newton_iter; ++newton_iter) {
        next_pos = q_trj_per_target[target_idx](next_t) + next_rel_pos;

        // std::cout << delta_t << " desired length " << vmax*delta_t << " shortest path length " << turns_for_one_sided_dubins_path(pos(0), pos(1), heading, next_pos(0), next_pos(1), rho).col(1).sum() << std::endl;
        RowMatrixXd turns = elongated_dubins_path_one_sided(pos(0), pos(1), heading, next_pos(0), next_pos(1), vmax*delta_t, rho);
        if (std::isfinite(turns(0, 0))) {
          max_newton_iter_for_success_repair = std::max(newton_iter, max_newton_iter_for_success_repair);
          newton_succeeded = true;
          for (int row = 0; row < turns.rows(); ++row) {
            if (turns(row, 0) != 0) {
              heading += turns(row, 1)/turns(row, 0);
            }
          }
          break;
        }

        RowMatrixXd shortest_path_turns = turns_for_one_sided_dubins_path(pos(0), pos(1), heading, next_pos(0), next_pos(1), rho);
        double shortest_path_dist = shortest_path_turns.col(1).sum();
        double c = shortest_path_dist - vmax*delta_t;

        // Finite-diff
        double delta_t_plus = delta_t + delta_delta_t_finite_diff;

        double next_t_plus = t + delta_t_plus;

        Vector2d next_pos_plus = q_trj_per_target[target_idx](next_t_plus) + next_rel_pos;

        RowMatrixXd shortest_path_turns_plus = turns_for_one_sided_dubins_path(pos(0), pos(1), heading, next_pos_plus(0), next_pos_plus(1), rho);
        double shortest_path_dist_plus = shortest_path_turns_plus.col(1).sum();
        double c_plus = shortest_path_dist_plus - vmax*delta_t_plus;

        double deriv = (c_plus - c)/delta_delta_t_finite_diff;

        double step_size = 0.01;
        delta_t -= step_size*c/deriv;

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
        X(seq_idx, 2) = delta_t;
      }

      if (!newton_succeeded) {
        repair_failed = true;
        break;
      }

      cost += next_t - tw_per_target(target_idx, 0);
    } else {
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

      next_t = t_high;

      delta_t = t_high - t;
      X(seq_idx, 2) = delta_t;
      cost += t_high_dist;
    }

    // Update cost, time, and position
    t = next_t;
    pos = next_pos;
  }
  if (repair_failed) {
    cost = std::numeric_limits<double>::infinity();
  }
  return repair_failed;
}

// next_heading is Ref<Vector1d> rather than double& so I can test in python
double find_earliest_arrival_time_dubins(Ref<Vector2d> next_pos, Ref<Vector1d> next_heading, const CppSpline &q_trj, const Ref<const Vector2d> &pos, double heading, double tw_start, double tw_end, double radius, double theta, double t, double vmax_agent, double rho, double vmax_target, int &max_bisection_iter_for_success_transformation) {
  Vector2d next_rel_pos = radius*Vector2d(cos(theta), sin(theta));

  Vector2d pos_t = q_trj(t) + next_rel_pos;

  RowMatrixXd turns = turns_for_one_sided_dubins_path(pos(0), pos(1), heading, pos_t(0), pos_t(1), rho);
  double length = turns.col(1).sum();

  double length_min = length;

  double t_min = t + length/(vmax_agent + vmax_target);
  double t_max = t + length/(vmax_agent - vmax_target);
  // double t_min = tw_start;
  // double t_max = tw_end;

  t_min = std::max(t_min, tw_start);
  t_max = std::min(t_max, tw_end);

  Vector2d pos_min = q_trj(t_min) + next_rel_pos;
  turns = turns_for_one_sided_dubins_path(pos(0), pos(1), heading, pos_min(0), pos_min(1), rho);
  bool CC_detected = turns(1, 0) != 0;
  length = turns.col(1).sum();
  double delta_min = length - vmax_agent*(t_min - t);

  Vector2d pos_max = q_trj(t_max) + next_rel_pos;
  turns = turns_for_one_sided_dubins_path(pos(0), pos(1), heading, pos_max(0), pos_max(1), rho);
  CC_detected |= turns(1, 0) != 0;
  length = turns.col(1).sum();
  double delta_max = length - vmax_agent*(t_max - t);
  if (delta_min*delta_max > 0) {
    // throw std::runtime_error("No bracket"); 
    return std::numeric_limits<double>::infinity();
  }

  int max_bisection_iter = 34; // On one of the 10 target instances I ran, we needed at most 17 iterations for successful transformation so I'm using 2x that number to declare failure
  double bisection_tol = 1e-4;
  for (int bisection_iter = 0; bisection_iter < max_bisection_iter; ++bisection_iter) {
    double t_mid = 0.5*(t_min + t_max);
    next_pos = q_trj(t_mid) + next_rel_pos;
    turns = turns_for_one_sided_dubins_path(pos(0), pos(1), heading, next_pos(0), next_pos(1), rho);
    CC_detected |= turns(1, 0) != 0;
    double length = turns.col(1).sum();
    double delta = length - vmax_agent*(t_mid - t);
    if (std::abs(delta) < bisection_tol) {
      next_heading(0) = heading + turns.col(1).sum()/rho; // I know the first column is all -rho, 0, or rho because we weren't doing any elongation
      max_bisection_iter_for_success_transformation = std::max(max_bisection_iter_for_success_transformation, bisection_iter);
      return t_mid;
    }
    if (delta*delta_min > 0) {
      t_min = t_mid;
    } else {
      t_max = t_mid;
    }
  }
  // throw std::runtime_error("Transformation method bisection ran out of iterations"); 
  if (!CC_detected) {
    throw std::runtime_error("Transformation bisection failed and no CC path detected");
  }
  return std::numeric_limits<double>::infinity();
}

bool transform_chromosome(Ref<MatrixXd> X, const Ref<const RowMatrixXd> &tw_per_target, const Ref<const VectorXd> &target_radii, const std::vector<CppSpline> &q_trj_per_target, const Ref<const Vector2d> &p0, double heading0, double vmax, double rho, double &cost, const Ref<const VectorXd> &speed_upper_bounds, int &max_newton_iter_for_success_repair, int &max_bisection_iter_for_success_transformation) {
  int num_targets = tw_per_target.rows();

  double t = 0;
  Vector2d pos = p0;
  double heading = heading0;
  Vector2d next_rel_pos;
  Vector2d next_pos;
  Vector1d next_heading;
  cost = 0;
  bool made_change = false;
  for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
    int target_idx = X(seq_idx, 0);

    if (made_change) {
      // Changing previous delta_t values may make next interception infeasible. If so, repair
      double tmp_cost;
      std::vector<CppSpline> tmp_q_trj;
      tmp_q_trj.push_back(q_trj_per_target[target_idx]);
      MatrixXd tmp_chromosome = X.block(seq_idx, 0, 1, gene_size);
      tmp_chromosome(0, 0) = 0;
      bool repair_failed = repair_chromosome(tmp_chromosome, tw_per_target.block(target_idx, 0, 1, 2), target_radii.segment(target_idx, 1), tmp_q_trj, pos, heading, vmax, tmp_cost, false, rho, max_newton_iter_for_success_repair);
      if (repair_failed) {
        cost = std::numeric_limits<double>::infinity();
        return false;
      }
      X.block(seq_idx, 0, 1, gene_size) = tmp_chromosome;
      X(seq_idx, 0) = target_idx;
    }

    double theta = X(seq_idx, 1);
    double delta_t = X(seq_idx, 2);

    if (seq_idx != num_targets - 1) {
      double next_t = t + delta_t;
      // Check the encounter pattern to the target at seq_idx + 1. If catchup, then optimize the arrival time to target at seq_idx.
      // If meeting, don't optimize the arrival time. The following method of determining the encounter pattern was obtained
      // from correspondence with the authors
      int following_target_idx = X(seq_idx + 1, 0);
      double following_theta = X(seq_idx + 1, 1);
      double following_delta_t = X(seq_idx + 1, 2);
      Vector2d following_rel_pos = target_radii[following_target_idx]*Vector2d(cos(following_theta), sin(following_theta));
      Vector2d following_pos = q_trj_per_target[following_target_idx](next_t + following_delta_t);

      RowMatrixXd turns = turns_for_one_sided_dubins_path(p0(0), p0(1), heading0, following_pos(0), following_pos(1), rho);

      double length = turns.col(1).sum();

      Vector2d following_pos_plus = q_trj_per_target[following_target_idx](next_t + following_delta_t + 0.1);
      turns = turns_for_one_sided_dubins_path(p0(0), p0(1), heading0, following_pos_plus(0), following_pos_plus(1), rho);

      double length_plus = turns.col(1).sum();

      if (length_plus < length) {
        next_rel_pos = target_radii[target_idx]*Vector2d(cos(theta), sin(theta));
        next_pos = q_trj_per_target[target_idx](next_t);

        // Meeting pattern
        cost += next_t - tw_per_target(target_idx, 0);
        t = next_t;
        pos = next_pos;
        heading = next_heading(0);

        continue;
      }
      // Catch-up pattern
    }

    double next_t = find_earliest_arrival_time_dubins(next_pos, next_heading, q_trj_per_target[target_idx], pos, heading, tw_per_target(target_idx, 0), tw_per_target(target_idx, 1), target_radii(target_idx), theta, t, vmax, rho, speed_upper_bounds(target_idx), max_bisection_iter_for_success_transformation);
    if (std::isinf(next_t)) {
      // std::cout << "Transformation failed on seq_idx " << seq_idx << std::endl;
      return false;
    }
    made_change = true;
    delta_t = next_t - t;
    X(seq_idx, 2) = delta_t;

    // Update cost, time, and position
    // Assume dubins
    cost += next_t - tw_per_target(target_idx, 0);
    t = next_t;
    pos = next_pos;
    heading = next_heading(0);
  }
  return true;
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

void get_speed_upper_bounds(Ref<VectorXd> upper_bounds, const std::vector<CppSpline> &q_trj_per_target, const Ref<const RowMatrixXd> &tw_per_target) {
  #pragma omp parallel for
  for (int target_idx = 0; target_idx < tw_per_target.rows(); ++target_idx) {
    // Need to evaluate at three points to fit a quadratic
    double t1 = tw_per_target(target_idx, 0);
    double t2 = tw_per_target(target_idx, 1);
    double t3 = 0.5*t1 + 0.5*t2;

    Matrix3d lhs;
    lhs(0, 0) = t1*t1;
    lhs(0, 1) = t1;
    lhs(0, 2) = 1;

    lhs(1, 0) = t2*t2;
    lhs(1, 1) = t2;
    lhs(1, 2) = 1;

    lhs(2, 0) = t3*t3;
    lhs(2, 1) = t3;
    lhs(2, 2) = 1;

    Matrix<double, 3, 2> rhs;
    rhs.row(0) = q_trj_per_target[target_idx].derivatives(t1).transpose();
    rhs.row(1) = q_trj_per_target[target_idx].derivatives(t2).transpose();
    rhs.row(2) = q_trj_per_target[target_idx].derivatives(t3).transpose();

    Matrix<double, 3, 2> coeffs = lhs.colPivHouseholderQr().solve(rhs);

    // x dimension
    double a = coeffs(0);
    double b = coeffs(1);
    double c = coeffs(2);
    double t_crit = -b/(2*a);

    double max_speed_x = std::max(std::abs(rhs(0, 0)), std::abs(rhs(1, 0)));
    max_speed_x = std::max(max_speed_x, std::abs(q_trj_per_target[target_idx].derivatives(t_crit)(0))); // Unnecessary computation of derivative in all dimensions here, fix if this overall procedure is slow

    // y dimension
    a = coeffs(0);
    b = coeffs(1);
    c = coeffs(2);
    t_crit = -b/(2*a);

    double max_speed_y = std::max(std::abs(rhs(0, 1)), std::abs(rhs(1, 1)));
    max_speed_y = std::max(max_speed_y, std::abs(q_trj_per_target[target_idx].derivatives(t_crit)(1))); // Unnecessary computation of derivative in all dimensions here, fix if this overall procedure is slow

    upper_bounds(target_idx) = sqrt(max_speed_x*max_speed_x + max_speed_y*max_speed_y);

    /*
    int num_sample = 100;
    for (int sample_idx = 0; sample_idx < num_sample; ++sample_idx) {
      double alpha = sample_idx/(double)num_sample;
      double t = alpha*t2 + (1 - alpha)*t1;
      Vector2d deriv = q_trj_per_target[target_idx].derivatives(t);
      if (deriv.norm() > upper_bounds(target_idx)) {
        std::cout << "Upper bound incorrect" << std::endl;
      }
    }
    */
  }
}

// initial_population should have number of rows equal to num_targets*pop_size, and gene_size columns
RowMatrixXd memetic_alg(Ref<RowMatrixXd> selected_pts_per_target, const Ref<const RowMatrixXd> &initial_population, const Ref<const VectorXd> &initial_costs, const std::vector<CppSpline> &q_trj_per_target_python, const Ref<const RowMatrixXd> &tw_per_target, const Ref<const VectorXd> &target_radii, double time_limit, double rho, double vmax, const Ref<const Vector2d> &p0, double heading0, int num_openmp_threads, std::vector<RowMatrixXd> &turns_chain) {
  std::vector<std::pair<double, double>> cost_vs_time;
  
  auto timer_start = std::chrono::high_resolution_clock::now();

  double min_cost_record_time = 0.;

  omp_set_num_threads(num_openmp_threads);

  int num_targets = tw_per_target.rows();

  // I'm doing this because I'm worried about the GIL
  std::vector<CppSpline> q_trj_per_target;
  for (int target_idx = 0; target_idx < num_targets; ++target_idx) {
    q_trj_per_target.push_back(CppSpline(q_trj_per_target_python[target_idx]));
  }

  bool dubins = rho != 0;

  VectorXd speed_upper_bounds(1);
  if (dubins) {
    speed_upper_bounds = VectorXd(num_targets);
    get_speed_upper_bounds(speed_upper_bounds, q_trj_per_target, tw_per_target);
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

  std::uniform_int_distribution<int> local_search_elite_distribution(0, pop_size/2);
  std::uniform_int_distribution<int> local_search_gene_idx_distribution(0, num_targets - 1);
  std::uniform_int_distribution<int> local_search_grad_vs_sampling_distribution(0, 1);

  std::uniform_real_distribution<double> local_search_sampling_theta_distribution(0, 2*M_PI);

  int Tlp = 2; // From paper

  const double mutation_prob = 0.1; // From paper

  std::vector<int> max_newton_iter_for_success_repair_per_thread(num_openmp_threads, 0);
  std::vector<int> max_bisection_iter_for_success_transformation_per_thread(num_openmp_threads, 0);

  // Initialize population
  #pragma omp parallel for
  for (int i = 0; i < pop_size; ++i) {
    population1[i] = initial_population.block(num_targets*i, 0, num_targets, gene_size);

    if (dubins) {
      double cost;
      if (repair_chromosome(population1[i], tw_per_target, target_radii, q_trj_per_target, p0, heading0, vmax, cost, dubins, rho, max_newton_iter_for_success_repair_per_thread[omp_get_thread_num()])) {
        population_costs1[i] = std::numeric_limits<double>::infinity();
      } else {
        population_costs1[i] = cost;
      }
    }
  }
  cost_vs_time.push_back(std::pair<double, double>(0., Map<VectorXd>(population_costs1.data(), population_costs1.size()).minCoeff()));
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
      bool repair_failed = repair_chromosome(Xnew, tw_per_target, target_radii, q_trj_per_target, p0, heading0, vmax, cost, dubins, rho, max_newton_iter_for_success_repair_per_thread[omp_get_thread_num()]);
      if (repair_failed || cost >= (*population_costs)[chromosome_idx]) {
        (*updated_population)[chromosome_idx] = (*population)[chromosome_idx];
        (*updated_population_costs)[chromosome_idx] = (*population_costs)[chromosome_idx];
        continue;
      }

      (*updated_population)[chromosome_idx] = Xnew;
      (*updated_population_costs)[chromosome_idx] = cost;

      // Transformation to reduce cost (only for Dubins)
      if (dubins) {
        bool success = transform_chromosome(Xnew, tw_per_target, target_radii, q_trj_per_target, p0, heading0, vmax, rho, cost, speed_upper_bounds, max_newton_iter_for_success_repair_per_thread[omp_get_thread_num()], max_bisection_iter_for_success_transformation_per_thread[omp_get_thread_num()]);
        if (success && cost < (*updated_population_costs)[chromosome_idx]) {
          (*updated_population)[chromosome_idx] = Xnew;
          (*updated_population_costs)[chromosome_idx] = cost;
          std::cout << "transformation reduced cost" << std::endl;
        }
      }
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

              bool repair_failed = repair_chromosome(local_modification, tw_per_target, target_radii, q_trj_per_target, p0, heading0, vmax, new_cost, dubins, rho, max_newton_iter_for_success_repair_per_thread[omp_get_thread_num()]);
              if (!repair_failed) {
                if (dubins) {
                  double tmp_cost;
                  bool transformation_succeeded = transform_chromosome(local_modification, tw_per_target, target_radii, q_trj_per_target, p0, heading0, vmax, rho, tmp_cost, speed_upper_bounds, max_newton_iter_for_success_repair_per_thread[omp_get_thread_num()], max_bisection_iter_for_success_transformation_per_thread[omp_get_thread_num()]);
                  if (transformation_succeeded && tmp_cost < new_cost) {
                    new_cost = tmp_cost;
                  }
                }

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

              bool repair_failed = repair_chromosome(local_modification, tw_per_target, target_radii, q_trj_per_target, p0, heading0, vmax, new_cost, dubins, rho, max_newton_iter_for_success_repair_per_thread[omp_get_thread_num()]);
              if (!repair_failed && new_cost < (*updated_population_costs)[chromosome_idx]) {
                if (dubins) {
                  double tmp_cost;
                  MatrixXd tmp_local_modification = local_modification;
                  bool transformation_succeeded = transform_chromosome(tmp_local_modification, tw_per_target, target_radii, q_trj_per_target, p0, heading0, vmax, rho, tmp_cost, speed_upper_bounds, max_newton_iter_for_success_repair_per_thread[omp_get_thread_num()], max_bisection_iter_for_success_transformation_per_thread[omp_get_thread_num()]);
                  if (transformation_succeeded && tmp_cost < new_cost) {
                    (*updated_population)[chromosome_idx] = tmp_local_modification;
                    (*updated_population_costs)[chromosome_idx] = tmp_cost;
                  } else {
                    (*updated_population)[chromosome_idx] = local_modification;
                    (*updated_population_costs)[chromosome_idx] = new_cost;
                  }
                } else {
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
          double theta = (*updated_population)[chromosome_idx](gene_idx, 1);

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

            bool repair_failed = repair_chromosome(local_modifications[sample_idx], tw_per_target, target_radii, q_trj_per_target, p0, heading0, vmax, new_cost, dubins, rho, max_newton_iter_for_success_repair_per_thread[omp_get_thread_num()]);

            if (!(repair_failed || new_cost >= (*updated_population_costs)[chromosome_idx])) {
              if (dubins) {
                double tmp_cost;
                MatrixXd tmp_local_modification = local_modifications[sample_idx];
                bool transformation_succeeded = transform_chromosome(tmp_local_modification, tw_per_target, target_radii, q_trj_per_target, p0, heading0, vmax, rho, tmp_cost, speed_upper_bounds, max_newton_iter_for_success_repair_per_thread[omp_get_thread_num()], max_bisection_iter_for_success_transformation_per_thread[omp_get_thread_num()]);
                if (transformation_succeeded && tmp_cost < new_cost) {
                  local_modification_costs[sample_idx] = tmp_cost;
                } else {
                  local_modification_costs[sample_idx] = new_cost;
                  local_modifications[sample_idx] = tmp_local_modification;
                }
              } else {
                local_modification_costs[sample_idx] = new_cost;
              }
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

  if (std::isinf(*it)) {
  } else if (dubins) {
    double t = 0;
    Vector2d pos = p0;
    Vector2d next_pos;
    double heading = heading0;
    for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
      int target_idx = Xbest(seq_idx, 0);
      double theta = Xbest(seq_idx, 1);
      double delta_t = Xbest(seq_idx, 2);
      t += delta_t;
      selected_pts_per_target(target_idx, 0) = t;
      next_pos = q_trj_per_target[target_idx](t) + target_radii[target_idx]*Vector2d(cos(theta), sin(theta));

      RowMatrixXd turns = elongated_dubins_path_one_sided(pos(0), pos(1), heading, next_pos(0), next_pos(1), vmax*delta_t, rho);
      if (std::isinf(turns(0, 0))) {
        throw std::runtime_error("elongated path generation failed");
      }

      pos = next_pos;
      for (int row = 0; row < turns.rows(); ++row) {
        if (turns(row, 0) != 0) {
          heading += turns(row, 1)/turns(row, 0);
        }
      }

      turns_chain.push_back(turns);

      selected_pts_per_target(target_idx, 1) = pos(0);
      selected_pts_per_target(target_idx, 2) = pos(1);
      selected_pts_per_target(target_idx, 3) = heading;

      if (t < tw_per_target(target_idx, 0) - 1e-4 || t > tw_per_target(target_idx, 1) + 1e-4) {
        throw std::runtime_error("t out of window");
      }
    }
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

  std::cout << "Max newton iterations for successful repair: " << *(std::min_element(max_newton_iter_for_success_repair_per_thread.begin(), max_newton_iter_for_success_repair_per_thread.end())) << std::endl;
  std::cout << "Max bisection iterations for successful transformation: " << *(std::min_element(max_bisection_iter_for_success_transformation_per_thread.begin(), max_bisection_iter_for_success_transformation_per_thread.end())) << std::endl;

  std::cout << "Spent " << min_cost_record_time << " s tracking what the min cost was after each iteration (just making sure this is not too large)" << std::endl;
  return cost_vs_time_mat;
}
