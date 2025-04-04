#pragma once
#include <Eigen/Dense>
#include <memory>
#include "mt_tsp_ros2/cbs_node.h"
#include "mt_tsp_ros2/ma_mt_tsp_soln.h"
#include <unordered_set>

class CBSSNode;
typedef std::shared_ptr<CBSSNode> CBSSNodePtr;

class CBSSNode {
  private:
    CBSConstraintList input_constraints;
    std::vector<GridPathPtr> paths; // Should have length num_agents
    double cost;
    std::shared_ptr<MAMTTSPSoln> ma_mt_tsp_soln;

  public:
    CBSSNode(const CBSConstraintList &input_constraints, const std::vector<GridPathPtr> &paths, double cost, std::shared_ptr<MAMTTSPSoln> ma_mt_tsp_soln) : input_constraints(input_constraints),
                                                                                                              paths(paths),
                                                                                                              cost(cost),
                                                                                                              ma_mt_tsp_soln(ma_mt_tsp_soln) {
    }

    const CBSConstraintList &get_input_constraints() const {
      return input_constraints;
    }

    const std::vector<GridPathPtr> &get_paths() {
      return paths;
    }

    double get_cost() const {
      return cost;
    }

    std::shared_ptr<MAMTTSPSoln> get_ma_mt_tsp_soln() const {
      return ma_mt_tsp_soln;
    }

    void get_path_through_joint_state_space(AStarPath &path) {
      path.clear();
      int max_steps = 0;
      int num_agents = paths.size();
      for (int i = 0; i < num_agents; ++i) {
        max_steps = std::max((int)(paths[i]->path.rows()), max_steps);
      }
      path.resize(max_steps);
      int dim_q = paths[0]->path.cols();
      for (int step = 0; step < max_steps; ++step) {
        path[step] = AStarCell(dim_q*num_agents);
        for (int agent = 0; agent < num_agents; ++agent) {
          int path_step = std::min((int)(paths[agent]->path.rows()) - 1, step);
          for (int i = 0; i < dim_q; ++i) {
            path[step](num_agents*i + agent) = paths[agent]->path(path_step, i);
          }
        }
      }
    }
};
