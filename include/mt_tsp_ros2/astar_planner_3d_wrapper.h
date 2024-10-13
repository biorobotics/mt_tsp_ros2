#pragma once
#include "astar_planner_3d.h"
#include <memory>

class AStarPlanner3DWrapper {
  public:
    AStarPlanner3DWrapper(const Ref<const VectorXd> &occupancy_flat, const Ref<const Vector3d> &map_lb, const Ref<const Vector3d> &map_ub, const Ref<const Vector3i> &ncells);

    MatrixXd plan(const Ref<const Vector3d>& start_pos, const Ref<const Vector3d>& goal_pos);

  private:
    std::unique_ptr<AStarPlanner3D> planner;
    std::vector<Vector3d> path;

    void unravel_index(VectorXi &unraveled_index, int raveled_index, const VectorXi &shape) {
      unraveled_index = VectorXi::Zero(shape.size());
      int prod = 1;
      for (int dim = shape.size() - 1; dim >= 0; --dim) {
        unraveled_index(dim) = (raveled_index/prod)%shape(dim);
        prod *= shape(dim);
      }
    }
};

