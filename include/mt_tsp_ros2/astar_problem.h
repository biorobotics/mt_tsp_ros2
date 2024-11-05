#pragma once
#include <Eigen/Dense>
#include "mt_tsp_ros2/astar_node.h"

using namespace Eigen;

class AStarProblem {
  public:
    virtual void generate_successors(std::vector<AStarCell> &succ, std::vector<double> &transition_costs, const AStarNodePtr &node) const {
    }

    virtual bool is_goal(const AStarNodePtr &node) const {
      return true;
    }

    virtual double heuristic(const AStarCell &cell) const {
      return 0;
    }
};
