#pragma once
#include "mt_tsp_ros2/cbss_problem.h"
#include "mt_tsp_ros2/ma_mt_tsp_soln.h"
#include <algorithm>
#include <fstream>

namespace py = pybind11;

struct compare_cbss_nodes {
  bool operator() (const CBSSNodePtr &elem1,
                   const CBSSNodePtr &elem2) {
    return elem1->get_cost() > elem2->get_cost();
  }
};

class CBSSData {
  protected:
    std::priority_queue<CBSSNodePtr, std::vector<CBSSNodePtr>, compare_cbss_nodes> open_list;
    std::shared_ptr<CBSSProblem> problem;
    py::object k_best_ma_mt_tsp_solver;
    double most_recent_ma_mt_tsp_soln_cost;
    int num_targets;

  public:
    // Adds start node to the open list
    CBSSData(const std::shared_ptr<CBSSProblem> &problem, py::object k_best_ma_mt_tsp_solver) {
      this->problem = problem;
      this->k_best_ma_mt_tsp_solver = k_best_ma_mt_tsp_solver;
      num_targets = k_best_ma_mt_tsp_solver.attr("get_num_targets")().cast<int>();
      CBSConstraint dummy_constraint(0, 
                                     -1*VectorXi::Ones(problem->get_dim_q()),
                                     -1*VectorXi::Ones(problem->get_dim_q()),
                                     -1,
                                     false);
      std::shared_ptr<MAMTTSPSoln> ma_mt_tsp_soln = std::make_shared<MAMTTSPSoln>(problem->get_num_agents(), num_targets);
      CBSSNodePtr start_node = nullptr;
      while (start_node == nullptr) {
        most_recent_ma_mt_tsp_soln_cost = k_best_ma_mt_tsp_solver.attr("solve")(std::ref(ma_mt_tsp_soln->get_cell_seqs()), std::ref(ma_mt_tsp_soln->get_time_seqs()), std::ref(ma_mt_tsp_soln->get_agent_cell_seq_start_ptr())).cast<double>();
        if (std::isinf(most_recent_ma_mt_tsp_soln_cost)) {
          std::cout << "Problem infeasible. Could not generate individual agent paths for any MAMTTSP solution" << std::endl;
          exit(1);
        }
        start_node = problem->generate_successor(nullptr, dummy_constraint, ma_mt_tsp_soln);
      }
      open(start_node);
    }

    void add_constraint(const CBSSNodePtr &node, const CBSConstraint &new_constraint) {
      CBSSNodePtr successor = problem->generate_successor(node, new_constraint, node->get_ma_mt_tsp_soln());
      if (successor != nullptr) {
        open(successor);
      }
    }

    CBSConstraintList get_new_constraints_from_conflicts(const CBSSNodePtr &node) {
      return problem->get_new_constraints_from_conflicts(node);
    }

    /*
     * open: pushes the node onto the open list
     * ARGUMENTS
     * node: node
     */
    void open(const CBSSNodePtr &node) {
      open_list.push(node);
    }

    /*
     * expand_next: pops off the node at the top of the open list
     * and returns it
     * RETURN: node previously at the top of the open list
     * REQUIRES: !open_list_empty()
     */
    CBSSNodePtr expand_next() {
      CBSSNodePtr ret = open_list.top();

      // First, check if we need to generate a new root
      if (ret->get_cost() > most_recent_ma_mt_tsp_soln_cost) {
        std::shared_ptr<MAMTTSPSoln> ma_mt_tsp_soln = std::make_shared<MAMTTSPSoln>(problem->get_num_agents(), num_targets);
        most_recent_ma_mt_tsp_soln_cost = k_best_ma_mt_tsp_solver.attr("solve")(std::ref(ma_mt_tsp_soln->get_cell_seqs()), std::ref(ma_mt_tsp_soln->get_time_seqs()), std::ref(ma_mt_tsp_soln->get_agent_cell_seq_start_ptr())).cast<double>();
        if (!std::isinf(most_recent_ma_mt_tsp_soln_cost)) {
          CBSConstraint dummy_constraint(0, 
                                         -1*VectorXi::Ones(problem->get_dim_q()),
                                         -1*VectorXi::Ones(problem->get_dim_q()),
                                         -1,
                                         false);
          CBSSNodePtr new_root_node = problem->generate_successor(nullptr, dummy_constraint, ma_mt_tsp_soln);
          open(new_root_node);
          ret = open_list.top();
        }
      }

      open_list.pop();
      return ret;
    }

    /*
     * open_list_empty: determines whether the open list is empty
     * RETURN: true if the open list is empty, false otherwise
     */
    bool open_list_empty() {
      return open_list.size() == 0;
    }
};
