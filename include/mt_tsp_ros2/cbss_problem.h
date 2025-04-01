#include <Eigen/Dense>
#include "astar_node.h"
#include "mt_tsp_ros2/cbss_node.h"
#include "mt_tsp_ros2/ma_mt_tsp_soln.h"
#include "mt_tsp_ros2/solve_seq_of_constrained_grid_2d_timed_astar_problems.h"
#include <pybind11/pybind11.h>

class CBSSProblem {
  public:
    CBSSProblem(const Ref<const Matrix<bool, Dynamic, Dynamic, RowMajor>> &occupancy,
                int num_agents,
                int low_level_planner_timeout_millis) : occupancy(occupancy),
                                                        num_agents(num_agents),
                                                        low_level_planner_timeout_millis(low_level_planner_timeout_millis) {
    }

    const Matrix<bool, Dynamic, Dynamic> &get_occupancy() const {
      return occupancy;
    }

    int get_num_agents() const {
      return num_agents;
    }

    int get_dim_q() const {
      return dim_q;
    }

    CBSSNodePtr generate_successor(const CBSSNodePtr &parent, const CBSConstraint new_constraint, std::shared_ptr<MAMTTSPSoln> ma_mt_tsp_soln) const {
      if (parent == nullptr) {
        double cost = 0;
        std::vector<GridPathPtr> paths(num_agents);

        // If there's no parent node, then we have to plan paths for all agents.
        // This should only be true for the start node
        AStarConstraintList individual_constraints;
        long T = 0;
        for (int agent = 0; agent < num_agents; ++agent) {
          T = std::max(T, ma_mt_tsp_soln->get_time_seq(agent).tail<1>()(0));
        }
        for (int agent = 0; agent < num_agents; ++agent) {
          paths[agent] = std::make_shared<GridPath>();
          // Needs to populate the path and its cost into paths[agent]
          bool solved = solve_seq_of_constrained_grid_2d_timed_astar_problems(paths[agent], occupancy, individual_constraints, ma_mt_tsp_soln->get_cell_seq(agent), ma_mt_tsp_soln->get_time_seq(agent), low_level_planner_timeout_millis, T);
          if (!solved) {
            return nullptr;
          }
          cost += paths[agent]->cost;
        }
        return std::make_shared<CBSSNode>(CBSConstraintList(), paths, cost, ma_mt_tsp_soln);
      } else {
        double cost = parent->get_cost();
        CBSConstraintList input_constraints;

        long T = 0;
        for (int agent = 0; agent < num_agents; ++agent) {
          T = std::max(T, ma_mt_tsp_soln->get_time_seq(agent).tail<1>()(0));
        }

        int agent = new_constraint.get_agent();

        input_constraints.insert(input_constraints.end(), parent->get_input_constraints().begin(), parent->get_input_constraints().end());
        AStarConstraintList individual_constraints;
        for (auto constraint : parent->get_input_constraints()) {
          if (constraint.get_agent() == agent) {
            individual_constraints.push_back(AStarConstraint(constraint.get_prev_cell(),
                                                             constraint.get_cell(),
                                                             constraint.get_time_step(),
                                                             constraint.is_edge_constraint()));
          }
        }

        individual_constraints.push_back(AStarConstraint(new_constraint.get_prev_cell(),
                                                         new_constraint.get_cell(),
                                                         new_constraint.get_time_step(),
                                                         new_constraint.is_edge_constraint()));
        input_constraints.push_back(new_constraint);

        std::vector<GridPathPtr> paths = parent->get_paths();
        double cost_before = paths[agent]->cost;
        paths[agent] = std::make_shared<GridPath>();

        bool solved = solve_seq_of_constrained_grid_2d_timed_astar_problems(paths[agent], occupancy, individual_constraints, ma_mt_tsp_soln->get_cell_seq(agent), ma_mt_tsp_soln->get_time_seq(agent), low_level_planner_timeout_millis, T);
        if (!solved) {
          return nullptr;
        }
        double cost_after = paths[agent]->cost;
        cost += cost_after - cost_before;

        return std::make_shared<CBSSNode>(input_constraints, paths, cost, ma_mt_tsp_soln);
      }
    }

    CBSConstraintList get_new_constraints_from_conflicts(const CBSSNodePtr &node) const {
      CBSConstraintList ret;
      const std::vector<GridPathPtr> &paths = node->get_paths();
      int max_steps = 0;
      for (int i = 0; i < num_agents; ++i) {
        max_steps = std::max((int)(paths[i]->path.rows()), max_steps);
      }
      // We assume there's no vertex conflict at step 0 because the problem would be infeasible otherwise
      for (int step = 1; step < max_steps; ++step) {
        for (int agent1 = 0; agent1 < num_agents; ++agent1) {
          for (int agent2 = agent1 + 1; agent2 < num_agents; ++agent2) {
            int path1_step = std::min((int)(paths[agent1]->path.rows()) - 1, step);
            int path2_step = std::min((int)(paths[agent2]->path.rows()) - 1, step);

            Eigen::DenseBase<Eigen::MatrixXi>::RowXpr coords1 = paths[agent1]->path.row(path1_step);
            Eigen::DenseBase<Eigen::MatrixXi>::RowXpr coords2 = paths[agent2]->path.row(path2_step);

            Eigen::DenseBase<Eigen::MatrixXi>::RowXpr prev_coords1 = paths[agent1]->path.row(path1_step - 1);
            Eigen::DenseBase<Eigen::MatrixXi>::RowXpr prev_coords2 = paths[agent2]->path.row(path2_step - 1);

            // Vertex conflict
            if (!(coords1 - coords2).any()) {
              ret.push_back(CBSConstraint(agent1,
                                          prev_coords1.transpose(),
                                          coords1.transpose(),
                                          step,
                                          false));
              ret.push_back(CBSConstraint(agent2,
                                          prev_coords2.transpose(),
                                          coords2.transpose(),
                                          step,
                                          false));
            }

            // Edge conflict
            if (!(coords1 - prev_coords2).any() && !(coords2 - prev_coords1).any()) {
              ret.push_back(CBSConstraint(agent1,
                                          prev_coords1.transpose(),
                                          coords1.transpose(),
                                          step,
                                          true));
              ret.push_back(CBSConstraint(agent2,
                                          prev_coords2.transpose(),
                                          coords2.transpose(),
                                          step,
                                          true));
            }
          }
        }
      }

      return ret;
    }

    void get_grid_path(MatrixXi &grid_path, const AStarPath &path) {
      grid_path = MatrixXi::Zero(path.size(), 2*num_agents);
      for (int i = 0; i < path.size(); ++i) {
        grid_path.row(i) = path[i].transpose();
      }
    }

  private:
    Matrix<bool, Dynamic, Dynamic> occupancy;
    int num_agents;
    int low_level_planner_timeout_millis;
    int dim_q = 2;
};
