#pragma once
#include <Eigen/Dense>
#include <memory>
#include "mt_tsp_ros2/astar_node.h"
#include <unordered_set>

class CBSConstraint {
  public:
    CBSConstraint(int agent, 
                  const AStarCell &prev_cell,
                  const AStarCell &cell,
                  int time_step,
                  bool edge_constraint) : agent(agent),
                                          prev_cell(prev_cell),
                                          cell(cell),
                                          time_step(time_step),
                                          edge_constraint(edge_constraint) {}
    int get_agent() const {
      return agent;
    }
    const AStarCell &get_prev_cell() const {
      return prev_cell;
    }
    const AStarCell &get_cell() const {
      return cell;
    }
    int get_time_step() const {
      return time_step;
    }
    bool is_edge_constraint() const {
      return edge_constraint;
    }
  private:
    int agent;
    AStarCell prev_cell;
    AStarCell cell;
    int time_step; 
    bool edge_constraint;
};

struct GridPath {
  MatrixXi path;
  double cost;
};

typedef std::vector<CBSConstraint> CBSConstraintList;
typedef std::shared_ptr<GridPath> GridPathPtr;

class CBSNode;
typedef std::shared_ptr<CBSNode> CBSNodePtr;

class CBSNode {
  private:
    CBSConstraintList input_constraints;
    std::vector<GridPathPtr> paths; // Should have length num_agents
    double cost;

  public:
    CBSNode(const CBSConstraintList &input_constraints, const std::vector<GridPathPtr> &paths, double cost) : input_constraints(input_constraints),
                                                                                                              paths(paths),
                                                                                                              cost(cost) {
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
