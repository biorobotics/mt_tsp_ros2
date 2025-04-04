#pragma once
#include "mt_tsp_ros2/astar_problem.h"
#include "mt_tsp_ros2/astar_constraint.h"
#include <memory>

typedef std::vector<AStarConstraint> AStarConstraintList;

// In this problem, each AStarCell contains the agent coordinates, and then the time step
class ConstrainedGrid2dTimedAStarProblem : public AStarProblem {
  public:
    ConstrainedGrid2dTimedAStarProblem(const Matrix<bool, Dynamic, Dynamic> &occupancy, const AStarConstraintList &constraints, int gx, int gy, int gt) : occupancy(occupancy), constraints(constraints), gx(gx), gy(gy), gt(gt) {
    }

    virtual void generate_successors(std::vector<AStarCell> &succ, std::vector<double> &transition_costs, const AStarNodePtr &node) const override {
      succ.clear();
      transition_costs.clear();

      const VectorXi &coords = node->get_cell();

      if (coords(2) == gt) {
        return;
      }

      VectorXi next_coords(3);
      next_coords(2) = coords(2) + 1;

      // No motion
      next_coords(0) = coords(0);
      next_coords(1) = coords(1);
      if (satisfies_constraints(coords, next_coords)) {
        succ.push_back(next_coords);
        transition_costs.push_back(1);
      }

      // +x
      if (coords(0) + 1 < occupancy.rows() && !occupancy(coords(0) + 1, coords(1))) {
        next_coords(0) = coords(0) + 1;
        next_coords(1) = coords(1);
        if (satisfies_constraints(coords, next_coords)) {
          succ.push_back(next_coords);
          transition_costs.push_back(1);
        }
      }
      // -x
      if (coords(0) - 1 >= 0 && !occupancy(coords(0) - 1, coords(1))) {
        next_coords(0) = coords(0) - 1;
        next_coords(1) = coords(1);
        if (satisfies_constraints(coords, next_coords)) {
          succ.push_back(next_coords);
          transition_costs.push_back(1);
        }
      }
      // +y
      if (coords(1) + 1 < occupancy.cols() && !occupancy(coords(0), coords(1) + 1)) {
        next_coords(0) = coords(0);
        next_coords(1) = coords(1) + 1;
        if (satisfies_constraints(coords, next_coords)) {
          succ.push_back(next_coords);
          transition_costs.push_back(1);
        }
      }
      // -y
      if (coords(1) - 1 >= 0 && !occupancy(coords(0), coords(1) - 1)) {
        next_coords(0) = coords(0);
        next_coords(1) = coords(1) - 1;
        if (satisfies_constraints(coords, next_coords)) {
          succ.push_back(next_coords);
          transition_costs.push_back(1);
        }
      }
    }

    virtual bool is_goal(const AStarNodePtr &node) const override {
      const VectorXi &coords = node->get_cell();
      return coords(0) == gx && coords(1) == gy && coords(2) == gt;
    }

    void get_grid_path(MatrixXi &grid_path, const AStarPath &path) {
      grid_path = MatrixXi::Zero(path.size(), 2);
      for (int i = 0; i < path.size(); ++i) {
        grid_path.row(i) = path[i].head(2).transpose();
      }
    }

    virtual double heuristic(const AStarCell &cell) const override {
      return abs(cell(0) - gx) + abs(cell(1) - gy);
    }

  private:
    Matrix<bool, Dynamic, Dynamic> occupancy;
    AStarConstraintList constraints;
    int gx;
    int gy;
    int gt;

    bool satisfies_constraints(const VectorXi &coords, const VectorXi &next_coords) const {
      for (auto constraint : constraints) {
        if (next_coords(2) == constraint.get_time_step() &&
            !((next_coords.head(2) - constraint.get_cell()).any() ||
              (constraint.is_edge_constraint() &&
               (coords.head(2) - constraint.get_prev_cell()).any()))) {
          return false;
        }
      }
      return true;
    }
};
