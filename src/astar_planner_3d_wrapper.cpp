#include "mt_tsp_ros2/astar_planner_3d_wrapper.h"
#include <iostream>

AStarPlanner3DWrapper::AStarPlanner3DWrapper(const Ref<const VectorXd> &occupancy_flat, const Ref<const Vector3d> &map_lb, const Ref<const Vector3d> &map_ub, const Ref<const Vector3i> &ncells){
  Tensor<bool, 3> occupancy(ncells(0), ncells(1), ncells(2));
  occupancy.setZero();
  VectorXi unraveled_index(3);
  for (int i = 0; i < occupancy_flat.size(); ++i) {
    if (!occupancy_flat(i)) {
      continue;
    }
    unravel_index(unraveled_index, i, ncells);
    occupancy(unraveled_index(0), unraveled_index(1), unraveled_index(2)) = true;
  }

  planner = std::make_unique<AStarPlanner3D>(occupancy, map_lb, map_ub, ncells);
  // planner->set_connected_26(false);
}

MatrixXd AStarPlanner3DWrapper::plan(const Ref<const Vector3d>& start_pos, const Ref<const Vector3d>& goal_pos) {
  path.clear();
  bool success = planner->plan_r3_path(path, start_pos, goal_pos);
  if (success) {
    MatrixXd ret = MatrixXd::Zero(path.size(), 3);
    for (int i = 0; i < path.size(); ++i) {
      ret.row(i) = path[i].transpose();
    }
    return ret;
  } else {
    MatrixXd ret = MatrixXd::Zero(1, 3);
    ret.row(0) = start_pos;
    return ret;
  }
}
