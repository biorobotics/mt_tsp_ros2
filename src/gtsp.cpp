#include "mt_tsp_ros2/gtsp.h"
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

struct EdgeEvalData {
  py::object edge_evaluator;
  std::vector<double> sol;
  int num_groups;
  RowMatrixXi selected;
  VectorXd updated_costs;
  std::vector<std::pair<int, int>> sol_to_edge_ptr;
  double callback_time;
  RowMatrixXdRef_const lb_cost_mat;

  std::vector<int> updated_cost_sol_idx;
  std::vector<double> updated_cost_constr_vals;

  EdgeEvalData(py::object edge_evaluator, int sol_size, int num_groups, const std::vector<std::pair<int, int>> &sol_to_edge_ptr, RowMatrixXdRef_const lb_cost_mat) : edge_evaluator(edge_evaluator), sol(sol_size, 0.), num_groups(num_groups), sol_to_edge_ptr(sol_to_edge_ptr), lb_cost_mat(lb_cost_mat) {
    selected = RowMatrixXi(num_groups + 1, 2);
    updated_costs = VectorXd(num_groups + 1);
    callback_time = 0.;
  }
};

int __stdcall lazy_edge_eval_cb(GRBmodel *model,
                                void     *cbdata,
                                int       where,
                                void     *usrdata) {
  int error = 0;
  if (where == GRB_CB_MIPSOL) {
    auto start_time = std::chrono::high_resolution_clock::now();
    struct EdgeEvalData *edge_eval_data = (struct EdgeEvalData *)usrdata;
    GRBcbget(cbdata, where, GRB_CB_MIPSOL_SOL, edge_eval_data->sol.data());
    int selected_idx = 0;
    std::vector<int> selected_sol_indices;
    // Stop for loop at sol size - 1 because sol size - 1 is the index of theta, which isn't an edge indicator var
    for (int sol_idx = 0; sol_idx < edge_eval_data->sol.size() - 1; ++sol_idx) {
      if (edge_eval_data->sol[sol_idx] > 0.5 && edge_eval_data->sol_to_edge_ptr[sol_idx].first != edge_eval_data->sol_to_edge_ptr[sol_idx].second) {
        edge_eval_data->selected(selected_idx, 0) = edge_eval_data->sol_to_edge_ptr[sol_idx].first;
        edge_eval_data->selected(selected_idx, 1) = edge_eval_data->sol_to_edge_ptr[sol_idx].second;
        selected_sol_indices.push_back(sol_idx);
        selected_idx += 1;
      }
    }

    edge_eval_data->edge_evaluator.attr("__call__")(std::ref(edge_eval_data->updated_costs), edge_eval_data->selected);
    std::vector<int> inf_sol_idx;
    std::vector<double> inf_sol_constr_vals;
    bool did_update_costs = false;
    for (int selected_idx = 0; selected_idx < edge_eval_data->num_groups + 1; ++selected_idx) {
      if (std::isnan(edge_eval_data->updated_costs(selected_idx))) {
        continue;
      }
      int sol_idx = selected_sol_indices[selected_idx];
      if (std::isinf(edge_eval_data->updated_costs(selected_idx))) {
        inf_sol_idx.push_back(sol_idx);
        inf_sol_constr_vals.push_back(1.);
        continue;
      }
      edge_eval_data->updated_cost_sol_idx.push_back(sol_idx);
      int node_idx1 = edge_eval_data->selected(selected_idx, 0);
      int node_idx2 = edge_eval_data->selected(selected_idx, 1);
      edge_eval_data->updated_cost_constr_vals.push_back(edge_eval_data->updated_costs(selected_idx) - edge_eval_data->lb_cost_mat(node_idx1, node_idx2));
      did_update_costs = true;
    }

    // Handle infeasible edges
    if (!error && inf_sol_idx.size()) {
      error = GRBcblazy(cbdata, inf_sol_idx.size(), inf_sol_idx.data(), inf_sol_constr_vals.data(), GRB_EQUAL, 0.);
      if (error) {
        printf("ERROR at callback edge prohibition: %s\n", GRBgeterrormsg(GRBgetenv(model)));
      }
    }

    // Handle updated costs
    if (!error && edge_eval_data->updated_cost_sol_idx.size() && (did_update_costs or inf_sol_idx.size() == 0)) {
      std::vector<int> updated_cost_sol_idx = edge_eval_data->updated_cost_sol_idx;
      std::vector<double> updated_cost_constr_vals = edge_eval_data->updated_cost_constr_vals;
      updated_cost_sol_idx.push_back(edge_eval_data->sol.size() - 1);
      updated_cost_constr_vals.push_back(-1.);
      error = GRBcblazy(cbdata, updated_cost_sol_idx.size(), updated_cost_sol_idx.data(), updated_cost_constr_vals.data(), GRB_LESS_EQUAL, 0.);
      if (error) {
        printf("ERROR at callback cost update: %s\n", GRBgeterrormsg(GRBgetenv(model)));
      }
    }
    auto stop_time = std::chrono::high_resolution_clock::now();
    edge_eval_data->callback_time += (double)(std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time).count())/1000;
  }
  return error;
}

VectorXd solve_gtsp_no_gsec(VectorXlRef node_seq, RowMatrixXdRef_const cost_mat, VectorXlRef_const group_start_idx, bool verbose, double time_limit, std::string save_path, bool take_first_feas_soln, double mipgap, bool solve_relaxed) {
  auto start_time = std::chrono::high_resolution_clock::now();
  int num_nodes = cost_mat.rows();
  int num_groups = group_start_idx.size() - 1;
  MatrixXi edges = -MatrixXi::Ones(num_nodes, num_nodes);

  int nx = 1;

  int vars_per_edge = nx;

  double xlb = 0.;
  double xub = 1.;

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
  if (error) quit(env, model);

  if (save_path.length()) {
    error = GRBsetstrparam(env, "LogFile", (save_path + "/log.txt").c_str());
    if (error) quit(env, model);
  }

  error = GRBstartenv(env);
  if (error) quit(env, model);

  error = GRBsetdblparam(env, "TimeLimit", time_limit);
  if (error) quit(env, model);
  error = GRBsetdblparam(env, "MIPGap", mipgap);
  if (error) quit(env, model);
  if (take_first_feas_soln) {
    error = GRBsetintparam(env, "SolutionLimit", 1);
    if (error) quit(env, model);
  }

  // Note that once we create a model, the model has its own copy of the env. So if we want to set params later,
  // we need to use GRBgetenv (see https://www.gurobi.com/documentation/current/refman/c_setdblparam.html)

  error = GRBnewmodel(env, &model, "ip", 0, NULL, NULL, NULL, NULL, NULL);
  if (error) quit(env, model);

  int edge_start_idx = 0;
  for (int node_idx1 = 0; node_idx1 < num_nodes; ++node_idx1) {
    for (int node_idx2 = 0; node_idx2 < num_nodes; ++node_idx2) {
      // This is different from the rest of my code, in the fact that
      // I'm allowing x_ii, following Laporte 1987. x_ii = 0
      // if node i is part of the tour, and 1 otherwise
      if (node_idx1 != node_idx2 && std::isinf(cost_mat(node_idx1, node_idx2))) {
        continue;
      }

      if (node_idx1 == node_idx2 && node_idx1 == 0) {
        continue;
      }

      edges(node_idx1, node_idx2) = edge_start_idx;

      error = GRBaddvar(model, 0, NULL, NULL, (node_idx1 == node_idx2 ? 0. : cost_mat(node_idx1, node_idx2)), xlb, xub, solve_relaxed ? GRB_CONTINUOUS : GRB_BINARY, ("x" + std::to_string(node_idx1) + "_" + std::to_string(node_idx2)).c_str());
      if (error) quit(env, model);

      edge_start_idx += vars_per_edge;
    }
  }

  int num_decision_vars = edge_start_idx;

  std::vector<int> ind;
  std::vector<double> val;

  // Conservation of flow
  for (int node_idx = 0; node_idx < num_nodes; ++node_idx) {
    int constr_start_idx = ind.size();
    int num_constr = 0;
    for (int node_idx1 = 0; node_idx1 < num_nodes; ++node_idx1) {
      if (node_idx1 != node_idx && edges(node_idx1, node_idx) != -1) {
        ind.push_back(edges(node_idx1, node_idx));
        val.push_back(1.);
        num_constr += 1;
      }
    }
    for (int node_idx2 = 0; node_idx2 < num_nodes; ++node_idx2) {
      if (node_idx != node_idx2 && edges(node_idx, node_idx2) != -1) {
        ind.push_back(edges(node_idx, node_idx2));
        val.push_back(-1.);
        num_constr += 1;
      }
    }
    if (num_constr) {
      error = GRBaddconstr(model, num_constr, &ind[constr_start_idx], &val[constr_start_idx], GRB_EQUAL, 0.0, ("flow_cons" + std::to_string(node_idx)).c_str());
      if (error) quit(env, model);
    }
  }

  // Sum of edges in = 1
  for (int node_idx = 0; node_idx < num_nodes; ++node_idx) {
    int constr_start_idx = ind.size();
    int num_constr = 0;
    for (int node_idx1 = 0; node_idx1 < num_nodes; ++node_idx1) {
      if (edges(node_idx1, node_idx) != -1) {
        ind.push_back(edges(node_idx1, node_idx));
        val.push_back(1.);
        num_constr += 1;
      }
    }

    error = GRBaddconstr(model, num_constr, &ind[constr_start_idx], &val[constr_start_idx], GRB_EQUAL, 1.0, ("enter_once" + std::to_string(node_idx)).c_str());
    if (error) quit(env, model);
  }

  for (int group_idx = 0; group_idx < group_start_idx.size() - 1; ++group_idx) {
    int constr_start_idx = ind.size();
    int num_constr = 0;
    int group_size = group_start_idx(group_idx + 1) - group_start_idx(group_idx);
    assert(group_size != 0);
    for (int node_idx = group_start_idx(group_idx); node_idx < group_start_idx(group_idx + 1); ++node_idx) {
      ind.push_back(edges(node_idx, node_idx));
      val.push_back(1.);
      num_constr += 1;
    }
    error = GRBaddconstr(model, num_constr, &ind[constr_start_idx], &val[constr_start_idx], GRB_EQUAL, group_size - 1, ("visit_group" + std::to_string(group_idx)).c_str());
    if (error) quit(env, model);
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

  std::vector<double> soln(num_decision_vars);
  for (int i = 0; i < num_decision_vars; ++i) {
    error = GRBgetdblattrelement(model, GRB_DBL_ATTR_X, i, &soln[i]);
    if (error) quit(env, model);
  }

  if (solve_relaxed) {
    assert(node_seq.size() == num_groups);
    // Populate node_seq with the nodes with min exclusion value per group (ith element is for group i + 1)
    for (int group_idx = 0; group_idx < num_groups; ++group_idx) {
      double min_val = std::numeric_limits<double>::infinity();
      for (int node_idx = group_start_idx[group_idx]; node_idx < group_start_idx[group_idx + 1]; ++node_idx) {
        if (soln[edges(node_idx, node_idx)] < min_val) {
          node_seq(group_idx) = node_idx;
          min_val = soln[edges(node_idx, node_idx)];
        }
      }
      assert(!std::isinf(min_val));
    }
  } else {
    assert(node_seq.size() == num_groups + 2);
    node_seq(0) = 0;
    for (int i = 1; i < num_groups + 2; ++i) {
      int node_idx1 = node_seq(i - 1);
      bool found = false;
      for (int node_idx2 = 0; node_idx2 < num_nodes; ++node_idx2) {
        if (edges(node_idx1, node_idx2) != -1 && soln[edges(node_idx1, node_idx2)] > 0.5) {
          node_seq(i) = node_idx2;
          found = true;
          break;
        }
      }
      assert(found);
    }
  }

  quit_no_error(env, model);
  return ret_vec;
}

struct PCGData {
  py::object send_callback;
  py::object recv_callback;
  std::vector<double> sol;
  int num_groups;
  RowMatrixXi selected;
  std::vector<std::pair<int, int>> sol_to_edge_ptr;
  double callback_time;

  PCGData(py::object send_callback, py::object recv_callback, int sol_size, int num_groups, const std::vector<std::pair<int, int>> &sol_to_edge_ptr) : send_callback(send_callback), recv_callback(recv_callback), sol(sol_size, 0.), num_groups(num_groups), sol_to_edge_ptr(sol_to_edge_ptr) {
    selected = RowMatrixXi(num_groups + 1, 2);
    callback_time = 0.;
  }
};

int __stdcall pcg_cb(GRBmodel *model,
                     void     *cbdata,
                     int       where,
                     void     *usrdata) {
  int error = 0;

  // Receive tours from other procs
  if (where == GRB_CB_MIP) {
    auto start_time = std::chrono::high_resolution_clock::now();
    struct PCGData *pcg_data = (struct PCGData *)usrdata;
    double lb = 0.;
    GRBcbget(cbdata, where, GRB_CB_MIP_OBJBND, &lb);
    bool terminate = pcg_data->recv_callback.attr("__call__")(lb).cast<bool>();
    if (terminate) {
      GRBterminate(model);
    }
    auto stop_time = std::chrono::high_resolution_clock::now();
    pcg_data->callback_time += (double)(std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time).count())/1000;
  }

  // Send tours to other procs
  if (where == GRB_CB_MIPSOL) {
    auto start_time = std::chrono::high_resolution_clock::now();
    struct PCGData *pcg_data = (struct PCGData *)usrdata;
    GRBcbget(cbdata, where, GRB_CB_MIPSOL_SOL, pcg_data->sol.data());
    int selected_idx = 0;
    std::vector<int> selected_sol_indices;
    // Stop for loop at sol size - 1 because sol size - 1 is the index of theta, which isn't an edge indicator var
    for (int sol_idx = 0; sol_idx < pcg_data->sol.size() - 1; ++sol_idx) {
      if (pcg_data->sol[sol_idx] > 0.5 && pcg_data->sol_to_edge_ptr[sol_idx].first != pcg_data->sol_to_edge_ptr[sol_idx].second) {
        pcg_data->selected(selected_idx, 0) = pcg_data->sol_to_edge_ptr[sol_idx].first;
        pcg_data->selected(selected_idx, 1) = pcg_data->sol_to_edge_ptr[sol_idx].second;
        selected_sol_indices.push_back(sol_idx);
        selected_idx += 1;
      }
    }

    double cost;
    GRBcbget(cbdata, where, GRB_CB_MIPSOL_OBJ, &cost);
    pcg_data->send_callback.attr("__call__")(pcg_data->selected, cost);
    auto stop_time = std::chrono::high_resolution_clock::now();
    pcg_data->callback_time += (double)(std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time).count())/1000;
  }
  return error;
}

VectorXd pcg_gtsp(VectorXlRef node_seq, RowMatrixXdRef soln_mat, py::object send_callback, py::object recv_callback, RowMatrixXdRef_const cost_mat, VectorXlRef_const group_start_idx, bool verbose, double time_limit, std::string save_path, double mipgap, VectorXlRef_const known_feas_tour, bool solve_relaxed, int proc_idx) {
  auto start_time = std::chrono::high_resolution_clock::now();
  int num_nodes = cost_mat.rows();
  int num_groups = group_start_idx.size() - 1;
  MatrixXi edges = -MatrixXi::Ones(num_nodes, num_nodes);

  int nx = 1;

  int vars_per_edge = nx;

  double xlb = 0.;
  double xub = 1.;

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
  if (error) quit(env, model);

  // error = GRBsetintparam(env, "Threads", 1);
  // if (error) quit(env, model);

  if (save_path.length()) {
    error = GRBsetstrparam(env, "LogFile", (save_path + "/log.txt").c_str());
    if (error) quit(env, model);
  }

  error = GRBstartenv(env);
  if (error) quit(env, model);

  error = GRBsetdblparam(env, "TimeLimit", time_limit);
  if (error) quit(env, model);
  error = GRBsetdblparam(env, "MIPGap", mipgap);
  if (error) quit(env, model);

  // Note that once we create a model, the model has its own copy of the env. So if we want to set params later,
  // we need to use GRBgetenv (see https://www.gurobi.com/documentation/current/refman/c_setdblparam.html)

  error = GRBnewmodel(env, &model, "ip", 0, NULL, NULL, NULL, NULL, NULL);
  if (error) quit(env, model);

  int edge_start_idx = 0;
  std::vector<std::pair<int, int>> sol_to_edge_ptr;
  for (int node_idx1 = 0; node_idx1 < num_nodes; ++node_idx1) {
    for (int node_idx2 = 0; node_idx2 < num_nodes; ++node_idx2) {
      // This is different from the rest of my code, in the fact that
      // I'm allowing x_ii, following Laporte 1987. x_ii = 0
      // if node i is part of the tour, and 1 otherwise
      if (node_idx1 != node_idx2 && std::isinf(cost_mat(node_idx1, node_idx2))) {
        continue;
      }

      if (node_idx1 == node_idx2 && node_idx1 == 0) {
        continue;
      }

      edges(node_idx1, node_idx2) = edge_start_idx;

      sol_to_edge_ptr.push_back(std::pair<int, int>(node_idx1, node_idx2));

      error = GRBaddvar(model, 0, NULL, NULL, (node_idx1 == node_idx2 ? 0. : cost_mat(node_idx1, node_idx2)), xlb, xub, solve_relaxed ? GRB_CONTINUOUS : GRB_BINARY, ("x" + std::to_string(node_idx1) + "_" + std::to_string(node_idx2)).c_str());
      if (error) quit(env, model);

      edge_start_idx += vars_per_edge;
    }
  }

  int num_decision_vars = edge_start_idx;

  std::vector<int> ind;
  std::vector<double> val;

  // Conservation of flow
  for (int node_idx = 0; node_idx < num_nodes; ++node_idx) {
    int constr_start_idx = ind.size();
    int num_constr = 0;
    for (int node_idx1 = 0; node_idx1 < num_nodes; ++node_idx1) {
      if (node_idx1 != node_idx && edges(node_idx1, node_idx) != -1) {
        ind.push_back(edges(node_idx1, node_idx));
        val.push_back(1.);
        num_constr += 1;
      }
    }
    for (int node_idx2 = 0; node_idx2 < num_nodes; ++node_idx2) {
      if (node_idx != node_idx2 && edges(node_idx, node_idx2) != -1) {
        ind.push_back(edges(node_idx, node_idx2));
        val.push_back(-1.);
        num_constr += 1;
      }
    }
    if (num_constr) {
      error = GRBaddconstr(model, num_constr, &ind[constr_start_idx], &val[constr_start_idx], GRB_EQUAL, 0.0, ("flow_cons" + std::to_string(node_idx)).c_str());
      if (error) quit(env, model);
    }
  }

  // Sum of edges in = 1
  for (int node_idx = 0; node_idx < num_nodes; ++node_idx) {
    int constr_start_idx = ind.size();
    int num_constr = 0;
    for (int node_idx1 = 0; node_idx1 < num_nodes; ++node_idx1) {
      if (edges(node_idx1, node_idx) != -1) {
        ind.push_back(edges(node_idx1, node_idx));
        val.push_back(1.);
        num_constr += 1;
      }
    }

    error = GRBaddconstr(model, num_constr, &ind[constr_start_idx], &val[constr_start_idx], GRB_EQUAL, 1.0, ("enter_once" + std::to_string(node_idx)).c_str());
    if (error) quit(env, model);
  }

  for (int group_idx = 0; group_idx < group_start_idx.size() - 1; ++group_idx) {
    int constr_start_idx = ind.size();
    int num_constr = 0;
    int group_size = group_start_idx(group_idx + 1) - group_start_idx(group_idx);
    assert(group_size != 0);
    for (int node_idx = group_start_idx(group_idx); node_idx < group_start_idx(group_idx + 1); ++node_idx) {
      ind.push_back(edges(node_idx, node_idx));
      val.push_back(1.);
      num_constr += 1;
    }
    error = GRBaddconstr(model, num_constr, &ind[constr_start_idx], &val[constr_start_idx], GRB_EQUAL, group_size - 1, ("visit_group" + std::to_string(group_idx)).c_str());
    if (error) quit(env, model);
  }

  if (known_feas_tour(0) != -1) {
    assert(known_feas_tour(0) == 0);
    assert(known_feas_tour(-1) == 0);

    std::unordered_set<int> used_self_edge_set;
    for (int node_idx = 1; node_idx < num_nodes; ++node_idx) {
      // Will erase unused elements next
      used_self_edge_set.insert(edges(node_idx, node_idx));
    }

    std::unordered_set<int> unused_self_edge_set;
    for (int seq_idx = 1; seq_idx < known_feas_tour.size() - 1; ++seq_idx) {
      unused_self_edge_set.insert(edges(known_feas_tour(seq_idx), known_feas_tour(seq_idx))); // Self-edges that are not used in the known feas tour. Again, exclude depot
      used_self_edge_set.erase(edges(known_feas_tour(seq_idx), known_feas_tour(seq_idx)));
    }

    std::unordered_set<int> unused_nonself_edge_set;
    for (int node_idx1 = 0; node_idx1 < num_nodes; ++node_idx1) {
      for (int node_idx2 = 0; node_idx2 < num_nodes; ++node_idx2) {
        if (edges(node_idx1, node_idx2) == -1 || node_idx1 == node_idx2) {
          continue;
        }
        // Will erase used elements next
        unused_nonself_edge_set.insert(edges(node_idx1, node_idx2));
      }
    }

    std::unordered_set<int> used_nonself_edge_set;
    for (int seq_idx = 0; seq_idx < known_feas_tour.size() - 1; ++seq_idx) {
      used_nonself_edge_set.insert(edges(known_feas_tour(seq_idx), known_feas_tour(seq_idx + 1)));
      unused_nonself_edge_set.erase(edges(known_feas_tour(seq_idx), known_feas_tour(seq_idx + 1)));
    }

    assert(used_self_edge_set.size() + used_nonself_edge_set.size() + unused_self_edge_set.size() + unused_nonself_edge_set.size() == num_decision_vars);

    for (auto edge : used_self_edge_set) {
      assert(edge != -1);
      GRBsetdblattrelement(model, GRB_DBL_ATTR_START, edge, 1);
    }

    for (auto edge : used_nonself_edge_set) {
      assert(edge != -1);
      GRBsetdblattrelement(model, GRB_DBL_ATTR_START, edge, 1);
    }

    for (auto edge : unused_self_edge_set) {
      assert(edge != -1);
      GRBsetdblattrelement(model, GRB_DBL_ATTR_START, edge, 0);
    }

    for (auto edge : unused_nonself_edge_set) {
      assert(edge != -1);
      GRBsetdblattrelement(model, GRB_DBL_ATTR_START, edge, 0);
    }
  }


  std::shared_ptr<struct PCGData> pcg_data_ptr;
  if (!solve_relaxed) {
    pcg_data_ptr = std::make_shared<struct PCGData>(send_callback, recv_callback, num_decision_vars, num_groups, sol_to_edge_ptr);
    error = GRBsetcallbackfunc(model, pcg_cb, (void *) pcg_data_ptr.get());
  }
  if (error) quit(env, model);

  auto stop_time = std::chrono::high_resolution_clock::now();
  auto setup_time = (double)(std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time).count())/1000;
  start_time = std::chrono::high_resolution_clock::now();
  error = GRBoptimize(model);
  if (error) quit(env, model);
  stop_time = std::chrono::high_resolution_clock::now();
  auto solve_time = (double)(std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time).count())/1000;

  VectorXd ret_vec = std::numeric_limits<double>::infinity()*VectorXd::Ones(4);
  ret_vec(0) = setup_time;
  ret_vec(1) = solve_time;
  if (solve_relaxed) {
    ret_vec(2) = 0.;
  } else {
    ret_vec(2) = pcg_data_ptr->callback_time;
  }

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

  error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &ret_vec(3));
  if (error) quit(env, model);

  std::vector<double> soln(num_decision_vars);
  for (int i = 0; i < num_decision_vars; ++i) {
    error = GRBgetdblattrelement(model, GRB_DBL_ATTR_X, i, &soln[i]);
    if (error) quit(env, model);
  }

  
  if (soln_mat.rows() == num_nodes and soln_mat.cols() == num_nodes) {
    for (int node_idx1 = 0; node_idx1 < num_nodes; ++node_idx1) {
      for (int node_idx2 = 0; node_idx2 < num_nodes; ++node_idx2) {
        if (edges(node_idx1, node_idx2) == -1) {
          continue; // Assume already contained infinity
        } else {
          soln_mat(node_idx1, node_idx2) = soln[edges(node_idx1, node_idx2)];
        }
      }
    }
  }

  if (solve_relaxed) {
    assert(node_seq.size() == num_groups);
    // Populate node_seq with the nodes with min exclusion value per group (ith element is for group i + 1)
    for (int group_idx = 0; group_idx < num_groups; ++group_idx) {
      double min_val = std::numeric_limits<double>::infinity();
      for (int node_idx = group_start_idx[group_idx]; node_idx < group_start_idx[group_idx + 1]; ++node_idx) {
        if (soln[edges(node_idx, node_idx)] < min_val) {
          node_seq(group_idx) = node_idx;
          min_val = soln[edges(node_idx, node_idx)];
        }
      }
      assert(!std::isinf(min_val));
    }
  } else {
    assert(node_seq.size() == num_groups + 2);
    node_seq(0) = 0;
    for (int i = 1; i < num_groups + 2; ++i) {
      int node_idx1 = node_seq(i - 1);
      bool found = false;
      for (int node_idx2 = 0; node_idx2 < num_nodes; ++node_idx2) {
        if (edges(node_idx1, node_idx2) != -1 && soln[edges(node_idx1, node_idx2)] > 0.5) {
          node_seq(i) = node_idx2;
          found = true;
          break;
        }
      }
      assert(found);
    }
  }
  quit_no_error(env, model);
  return ret_vec;
}

VectorXd solve_gtsp_no_gsec_lazy_edge_eval(VectorXlRef node_seq, py::object edge_evaluator, RowMatrixXdRef_const lb_cost_mat, VectorXlRef_const group_start_idx, bool verbose, double time_limit, std::string save_path, bool take_first_feas_soln, double mipgap, VectorXlRef_const known_feas_tour) {
  auto start_time = std::chrono::high_resolution_clock::now();
  int num_nodes = lb_cost_mat.rows();
  int num_groups = group_start_idx.size() - 1;
  MatrixXi edges = -MatrixXi::Ones(num_nodes, num_nodes);

  int nx = 1;

  int vars_per_edge = nx;

  double xlb = 0.;
  double xub = 1.;

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
  if (error) quit(env, model);

  if (save_path.length()) {
    error = GRBsetstrparam(env, "LogFile", (save_path + "/log.txt").c_str());
    if (error) quit(env, model);
  }

  error = GRBstartenv(env);
  if (error) quit(env, model);

  error = GRBsetdblparam(env, "TimeLimit", time_limit);
  if (error) quit(env, model);
  error = GRBsetdblparam(env, "MIPGap", mipgap);
  if (error) quit(env, model);
  if (take_first_feas_soln) {
    error = GRBsetintparam(env, "SolutionLimit", 1);
    if (error) quit(env, model);
  }
  error = GRBsetintparam(env, "LazyConstraints", 1);

  // Note that once we create a model, the model has its own copy of the env. So if we want to set params later,
  // we need to use GRBgetenv (see https://www.gurobi.com/documentation/current/refman/c_setdblparam.html)

  error = GRBnewmodel(env, &model, "ip", 0, NULL, NULL, NULL, NULL, NULL);
  if (error) quit(env, model);

  int edge_start_idx = 0;
  std::vector<std::pair<int, int>> sol_to_edge_ptr;
  for (int node_idx1 = 0; node_idx1 < num_nodes; ++node_idx1) {
    for (int node_idx2 = 0; node_idx2 < num_nodes; ++node_idx2) {
      // This is different from the rest of my code, in the fact that
      // I'm allowing x_ii, following Laporte 1987. x_ii = 0
      // if node i is part of the tour, and 1 otherwise
      if (node_idx1 != node_idx2 && std::isinf(lb_cost_mat(node_idx1, node_idx2))) {
        continue;
      }

      if (node_idx1 == node_idx2 && node_idx1 == 0) {
        continue;
      }

      edges(node_idx1, node_idx2) = edge_start_idx;
      sol_to_edge_ptr.push_back(std::pair<int, int>(node_idx1, node_idx2));

      error = GRBaddvar(model, 0, NULL, NULL, (node_idx1 == node_idx2 ? 0. : lb_cost_mat(node_idx1, node_idx2)), xlb, xub, GRB_BINARY, ("x" + std::to_string(node_idx1) + "_" + std::to_string(node_idx2)).c_str());
      if (error) quit(env, model);

      edge_start_idx += vars_per_edge;
    }
  }

  error = GRBaddvar(model, 0, NULL, NULL, 1., 0., GRB_INFINITY, GRB_CONTINUOUS, "theta");
  if (error) quit(env, model);

  int num_decision_vars = edge_start_idx + 1; // +1 because of theta

  std::vector<int> ind;
  std::vector<double> val;

  // Conservation of flow
  for (int node_idx = 0; node_idx < num_nodes; ++node_idx) {
    int constr_start_idx = ind.size();
    int num_constr = 0;
    for (int node_idx1 = 0; node_idx1 < num_nodes; ++node_idx1) {
      if (node_idx1 != node_idx && edges(node_idx1, node_idx) != -1) {
        ind.push_back(edges(node_idx1, node_idx));
        val.push_back(1.);
        num_constr += 1;
      }
    }
    for (int node_idx2 = 0; node_idx2 < num_nodes; ++node_idx2) {
      if (node_idx != node_idx2 && edges(node_idx, node_idx2) != -1) {
        ind.push_back(edges(node_idx, node_idx2));
        val.push_back(-1.);
        num_constr += 1;
      }
    }
    if (num_constr) {
      error = GRBaddconstr(model, num_constr, &ind[constr_start_idx], &val[constr_start_idx], GRB_EQUAL, 0.0, ("flow_cons" + std::to_string(node_idx)).c_str());
      if (error) quit(env, model);
    }
  }

  // Sum of edges in = 1
  for (int node_idx = 0; node_idx < num_nodes; ++node_idx) {
    int constr_start_idx = ind.size();
    int num_constr = 0;
    for (int node_idx1 = 0; node_idx1 < num_nodes; ++node_idx1) {
      if (edges(node_idx1, node_idx) != -1) {
        ind.push_back(edges(node_idx1, node_idx));
        val.push_back(1.);
        num_constr += 1;
      }
    }

    error = GRBaddconstr(model, num_constr, &ind[constr_start_idx], &val[constr_start_idx], GRB_EQUAL, 1.0, ("enter_once" + std::to_string(node_idx)).c_str());
    if (error) quit(env, model);
  }

  for (int group_idx = 0; group_idx < group_start_idx.size() - 1; ++group_idx) {
    int constr_start_idx = ind.size();
    int num_constr = 0;
    int group_size = group_start_idx(group_idx + 1) - group_start_idx(group_idx);
    assert(group_size != 0);
    for (int node_idx = group_start_idx(group_idx); node_idx < group_start_idx(group_idx + 1); ++node_idx) {
      ind.push_back(edges(node_idx, node_idx));
      val.push_back(1.);
      num_constr += 1;
    }
    error = GRBaddconstr(model, num_constr, &ind[constr_start_idx], &val[constr_start_idx], GRB_EQUAL, group_size - 1, ("visit_group" + std::to_string(group_idx)).c_str());
    if (error) quit(env, model);
  }

  struct EdgeEvalData edge_eval_data(edge_evaluator, num_decision_vars, num_groups, sol_to_edge_ptr, lb_cost_mat);
  error = GRBsetcallbackfunc(model, lazy_edge_eval_cb, (void *) &edge_eval_data);
  if (error) quit(env, model);

  if (known_feas_tour(0) != -1) {
    assert(known_feas_tour(0) == 0);
    assert(known_feas_tour(-1) == 0);

    std::unordered_set<int> used_self_edge_set;
    for (int node_idx = 1; node_idx < num_nodes; ++node_idx) {
      // Will erase unused elements next
      used_self_edge_set.insert(edges(node_idx, node_idx));
    }

    std::unordered_set<int> unused_self_edge_set;
    for (int seq_idx = 1; seq_idx < known_feas_tour.size() - 1; ++seq_idx) {
      unused_self_edge_set.insert(edges(known_feas_tour(seq_idx), known_feas_tour(seq_idx))); // Self-edges that are not used in the known feas tour. Again, exclude depot
      used_self_edge_set.erase(edges(known_feas_tour(seq_idx), known_feas_tour(seq_idx)));
    }

    std::unordered_set<int> unused_nonself_edge_set;
    for (int node_idx1 = 0; node_idx1 < num_nodes; ++node_idx1) {
      for (int node_idx2 = 0; node_idx2 < num_nodes; ++node_idx2) {
        if (edges(node_idx1, node_idx2) == -1 || node_idx1 == node_idx2) {
          continue;
        }
        // Will erase used elements next
        unused_nonself_edge_set.insert(edges(node_idx1, node_idx2));
      }
    }

    std::unordered_set<int> used_nonself_edge_set;
    for (int seq_idx = 0; seq_idx < known_feas_tour.size() - 1; ++seq_idx) {
      used_nonself_edge_set.insert(edges(known_feas_tour(seq_idx), known_feas_tour(seq_idx + 1)));
      unused_nonself_edge_set.erase(edges(known_feas_tour(seq_idx), known_feas_tour(seq_idx + 1)));
    }

    assert(used_self_edge_set.size() + used_nonself_edge_set.size() + unused_self_edge_set.size() + unused_nonself_edge_set.size() == num_decision_vars);

    for (auto edge : used_self_edge_set) {
      assert(edge != -1);
      GRBsetdblattrelement(model, GRB_DBL_ATTR_START, edge, 1);
    }

    for (auto edge : used_nonself_edge_set) {
      assert(edge != -1);
      GRBsetdblattrelement(model, GRB_DBL_ATTR_START, edge, 1);
    }

    for (auto edge : unused_self_edge_set) {
      assert(edge != -1);
      GRBsetdblattrelement(model, GRB_DBL_ATTR_START, edge, 0);
    }

    for (auto edge : unused_nonself_edge_set) {
      assert(edge != -1);
      GRBsetdblattrelement(model, GRB_DBL_ATTR_START, edge, 0);
    }

    GRBsetdblattrelement(model, GRB_DBL_ATTR_START, num_decision_vars - 1, 0.); // Theta
  }

  auto stop_time = std::chrono::high_resolution_clock::now();
  auto setup_time = (double)(std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time).count())/1000;
  start_time = std::chrono::high_resolution_clock::now();
  error = GRBoptimize(model);
  if (error) quit(env, model);
  stop_time = std::chrono::high_resolution_clock::now();
  auto solve_time = (double)(std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time).count())/1000;

  VectorXd ret_vec = std::numeric_limits<double>::infinity()*VectorXd::Ones(5);
  ret_vec(0) = setup_time;
  ret_vec(1) = solve_time;
  ret_vec(2) = edge_eval_data.callback_time;

  
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

  error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &ret_vec(3));
  if (error) quit(env, model);
  error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJBOUND, &ret_vec(4));
  if (error) quit(env, model);

  std::vector<double> soln(num_decision_vars);
  for (int i = 0; i < num_decision_vars; ++i) {
    error = GRBgetdblattrelement(model, GRB_DBL_ATTR_X, i, &soln[i]);
    if (error) quit(env, model);
  }

  assert(node_seq.size() == num_groups + 2);

  node_seq(0) = 0;
  for (int i = 1; i < num_groups + 2; ++i) {
    int node_idx1 = node_seq(i - 1);
    bool found = false;
    for (int node_idx2 = 0; node_idx2 < num_nodes; ++node_idx2) {
      if (edges(node_idx1, node_idx2) != -1 && soln[edges(node_idx1, node_idx2)] > 0.5) {
        node_seq(i) = node_idx2;
        found = true;
        break;
      }
    }
    assert(found);
  }
  quit_no_error(env, model);
  return ret_vec;
}
