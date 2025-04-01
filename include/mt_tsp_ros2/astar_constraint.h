#pragma once
#include <Eigen/Dense>
#include "mt_tsp_ros2/astar_node.h"
using namespace Eigen;

class AStarConstraint {
  public:
    AStarConstraint(const AStarCell &prev_cell,
                    const AStarCell &cell,
                    int time_step,
                    bool edge_constraint) : prev_cell(prev_cell),
                                            cell(cell),
                                            time_step(time_step),
                                            edge_constraint(edge_constraint) {}
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
    AStarCell prev_cell;
    AStarCell cell;
    int time_step; 
    bool edge_constraint;
};
