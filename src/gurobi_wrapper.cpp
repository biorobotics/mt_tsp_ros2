#include "mt_tsp_ros2/gurobi_wrapper.h"
#include "gurobi_c.h"
#include "stdlib.h"
#include "stdio.h"
#include <chrono>
#include <iostream>

using namespace std::chrono;

void quit(GRBenv *env, GRBmodel* model) {
  printf("ERROR: %s\n", GRBgeterrormsg(env));

  /* Free model */

  GRBfreemodel(model);

  /* Free environment */

  GRBfreeenv(env);

  exit(1);
}

void quit_no_error(GRBenv *env, GRBmodel* model) {
  /* Free model */

  GRBfreemodel(model);

  /* Free environment */

  GRBfreeenv(env);
}

VectorXd solve_socp(VectorXdRef soln, VectorXlRef_const A_csr_indptr, VectorXlRef_const A_csr_indices, VectorXdRef_const A_csr_data, VectorXdRef_const b, VectorXlRef_const Prows, VectorXlRef_const Pcols, VectorXdRef_const Pvals, VectorXdRef_const gradient, int num_zero_cone, int num_linear_cone, bool add_soc, int num_ctrl_pts, int li_dim, int l_idx, int dim_q, int vars_per_step, int steps, int num_decision_vars, bool verbose, double time_limit) {
  auto start_time = std::chrono::high_resolution_clock::now();

  GRBenv   *env   = NULL;
  GRBmodel *model = NULL;

  int error = 0;
  error = GRBemptyenv(&env);
  if (error) quit(env, model);

  if (verbose) {
    error = GRBsetintparam(env, "OutputFlag", 1);
  } else {
    error = GRBsetintparam(env, "OutputFlag", 0);
  }

  error = GRBstartenv(env);
  if (error) quit(env, model);

  if (error) quit(env, model);
  error = GRBsetdblparam(env, "TimeLimit", time_limit);
  if (error) quit(env, model);

  // Note that once we create a model, the model has its own copy of the env. So if we want to set params later,
  // we need to use GRBgetenv (see https://www.gurobi.com/documentation/current/refman/c_setdblparam.html)

  error = GRBnewmodel(env, &model, "socp", 0, NULL, NULL, NULL, NULL, NULL);
  if (error) quit(env, model);

  for (int i = 0; i < num_decision_vars; ++i) {
    error = GRBaddvar(model, 0, NULL, NULL, gradient(i), -GRB_INFINITY, GRB_INFINITY, GRB_CONTINUOUS, ("x" + std::to_string(i)).c_str());
    if (error) quit(env, model);
  }

  if (Prows.size()) {
    std::vector<int> Prows_vec(Prows.size());
    std::vector<int> Pcols_vec(Prows.size());
    std::vector<double> Pvals_vec(Prows.size());
    VectorXi::Map(Prows_vec.data(), Prows.size()) = Prows.cast<int>();
    VectorXi::Map(Pcols_vec.data(), Prows.size()) = Pcols.cast<int>();
    VectorXd::Map(Pvals_vec.data(), Prows.size()) = 0.5*Pvals;
    GRBaddqpterms(model, Prows.size(), Prows_vec.data(), Pcols_vec.data(), Pvals_vec.data());
  }

  std::vector<int> ind;
  std::vector<double> val;

  for (int row = 0; row < num_zero_cone; ++row) { 
    int constr_start_idx = ind.size();
    for (int col_idx = A_csr_indptr(row); col_idx < A_csr_indptr(row + 1); ++col_idx) {
      ind.push_back(A_csr_indices(col_idx));
      val.push_back(A_csr_data(col_idx));
    }
    int num_constr = A_csr_indptr(row + 1) - A_csr_indptr(row);
    error = GRBaddconstr(model, num_constr, &ind[constr_start_idx], &val[constr_start_idx], GRB_EQUAL, b(row), ("cons" + std::to_string(row)).c_str());
    if (error) quit(env, model);
  }

  for (int row = num_zero_cone; row < num_zero_cone + num_linear_cone; ++row) { 
    int constr_start_idx = ind.size();
    for (int col_idx = A_csr_indptr(row); col_idx < A_csr_indptr(row + 1); ++col_idx) {
      ind.push_back(A_csr_indices(col_idx));
      val.push_back(A_csr_data(col_idx));
    }
    int num_constr = A_csr_indptr(row + 1) - A_csr_indptr(row);
    error = GRBaddconstr(model, num_constr, &ind[constr_start_idx], &val[constr_start_idx], GRB_LESS_EQUAL, b(row), ("cons" + std::to_string(row)).c_str());
    if (error) quit(env, model);
  }

  if (add_soc) {
    for (int li_idx = 0; li_idx < num_ctrl_pts - 1; ++li_idx) {
      for (int i = l_idx + li_dim*li_idx + dim_q; i < vars_per_step*steps; i += vars_per_step) {
        error = GRBsetdblattrelement(model, "LB", i, 0.);
        if (error) quit(env, model);
      }
    }
    for (int step = 0; step < steps; ++step) {
      for (int ctrl_pt_idx = 0; ctrl_pt_idx < num_ctrl_pts - 1; ++ctrl_pt_idx) {
        int constr_start_idx = ind.size();
        for (int i = vars_per_step*step + l_idx + li_dim*ctrl_pt_idx; i < vars_per_step*step + l_idx + li_dim*ctrl_pt_idx + li_dim - 1; ++i) {
          ind.push_back(i);
          val.push_back(1.);
        }
        ind.push_back(vars_per_step*step + l_idx + li_dim*ctrl_pt_idx + li_dim - 1);
        val.push_back(-1.);
        error = GRBaddqconstr(model, 0, NULL, NULL, li_dim, &ind[constr_start_idx], &ind[constr_start_idx], &val[constr_start_idx], GRB_LESS_EQUAL, 0., ("soc" + std::to_string(step) + "_" + std::to_string(ctrl_pt_idx)).c_str());
        if (error) quit(env, model);
      }
    }
  }

  auto stop_time = std::chrono::high_resolution_clock::now();
  auto setup_time = (double)(std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time).count())/1000;
  start_time = std::chrono::high_resolution_clock::now();
  error = GRBoptimize(model);
  if (error) quit(env, model);
  stop_time = std::chrono::high_resolution_clock::now();
  auto solve_time = (double)(std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time).count())/1000;

  VectorXd ret_vec = std::numeric_limits<double>::infinity()*VectorXd::Ones(3);
  ret_vec(0) = setup_time;
  ret_vec(1) = solve_time;

  
  /*
  start_time = std::chrono::high_resolution_clock::now();
  GRBwrite(model, "model.mps");
  stop_time = std::chrono::high_resolution_clock::now();
  auto write_time = (double)(std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time).count())/1000;
  std::cout << "Write time: " << write_time << std::endl;
  */

  int solcount = 0;
  error = GRBgetintattr(model, GRB_INT_ATTR_SOLCOUNT, &solcount);
  if (error) quit(env, model);
  if (solcount == 0) {
    quit_no_error(env, model);
    return ret_vec;
  }

  error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &ret_vec(2));
  if (error) quit(env, model);

  assert(soln.size() == num_decision_vars);
  for (int i = 0; i < num_decision_vars; ++i) {
    error = GRBgetdblattrelement(model, GRB_DBL_ATTR_X, i, &soln(i));
    if (error) quit(env, model);
  }

  quit_no_error(env, model);
  return ret_vec;
}
