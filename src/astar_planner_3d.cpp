#include "mt_tsp_ros2/astar_planner_3d.h"
#include <iostream>
#include <limits>
#include "mt_tsp_ros2/astar_data_3d.h"
#include <chrono>

using namespace std::chrono;

bool AStarPlanner3D::is_cell_valid(int x, int y, int z) {
  return x >= 0 && x < ncells(0) && 
         y >= 0 && y < ncells(1) &&
         z >= 0 && z < ncells(2) &&
         !occupancy(x, y, z);
}

bool AStarPlanner3D::is_cell_valid(const Vector3i &cell) {
  return is_cell_valid(cell(0), cell(1), cell(2));
}

Vector3i AStarPlanner3D::pos_to_cell(const Ref<const Vector3d>& pos) {
  return ((pos - map_lb).array()/cell_size.array()).matrix().cast<int>();
}

Vector3d AStarPlanner3D::cell_to_pos(const Ref<const Vector3i>& cell) {
  return Vector3d(cell(0), cell(1), cell(2)).cwiseProduct(cell_size) + cell_size/2 + map_lb;
}

AStarPlanner3D::AStarPlanner3D(const Tensor<bool, 3> &occupancy,
                           const Vector3d &map_lb,
                           const Vector3d &map_ub,
                           const Vector3i &ncells) : map_lb(map_lb), ncells(ncells), occupancy(occupancy)
                           {
  cell_size = (map_ub - map_lb).array()/ncells.cast<double>().array();
  max_successors = 26;
}

bool AStarPlanner3D::plan_r3_path(std::vector<Vector3d>& path, 
                                const Ref<const Vector3d>& start_pos, const Ref<const Vector3d>& goal_pos) {
  // Return if start position is occupied. 
  Vector3i start = pos_to_cell(start_pos);
  if (!is_cell_valid(start)) {
    std::cout << "Start position is in collision" << std::endl;
    path.resize(1);
    path[0] = start_pos;
    return false;
  }

  // Return if goal position is occupied. 
  Vector3i goal = pos_to_cell(goal_pos);
  if (!is_cell_valid(goal)) {
    std::cout << "Goal position is in collision" << std::endl;
    path.resize(1);
    path[0] = start_pos;
    return false;
  }

  auto start_time = std::chrono::high_resolution_clock::now();

  int dX[26] = {-1, -1, -1,  0, 0,  1, 1, 1, -1, -1, -1,  0,  0,  1,  1,  1, -1, -1, -1,  0,  0,  1, 1, 1,  0, 0};    
  int dY[26] = {-1,  0,  1, -1, 1, -1, 0, 1, -1,  0,  1, -1,  1, -1,  0,  1, -1,  0,  1, -1,  1, -1, 0, 1,  0, 0};
  int dZ[26] = {0,   0,  0,  0, 0,  0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1, 1, 1, -1, 1};

  if (max_successors == 6) {
    dX[0] = -1;
    dX[1] = 1;
    for (int i = 2; i < 6; ++i) {
      dX[i] = 0;
    }

    dY[2] = -1;
    dY[3] = 1;
    for (int i = 0; i < 2; ++i) {
      dY[i] = 0;
    }
    for (int i = 4; i < 6; ++i) {
      dY[i] = 0;
    }

    dZ[4] = -1;
    dZ[5] = 1;
    for (int i = 0; i < 4; ++i) {
      dZ[i] = 0;
    }
  }

  int max_millis = 1000;

  float start_eps = 1.;

  AStarData3D astar(start(0), start(1), start(2),
                    goal(0), goal(1), goal(2), start_eps);

  vector<std::array<int, 3>> int_path;
  for (float eps = start_eps; eps >= 1.0f; eps = (eps == 1.0f) ? 0 : std::max(1.0f, eps / 2)) {
    astar.reset(eps);
    while (!astar.open_list_empty() &&
           astar.get_f(goal(0), goal(1), goal(2)) >
           astar.get_next().get_f()) {
      Node pop = astar.expand_next();
      NodePos pos = pop.get_pos();
      int x = pos.x;
      int y = pos.y;
      int z = pos.z;
      float g = pop.get_g();

      for (int i = 0; i < max_successors; ++i) {
        int newx = x + dX[i];
        int newy = y + dY[i];
        int newz = z + dZ[i];

        bool succ_valid = is_cell_valid(newx, newy, newz);
        bool adj_valid = (is_cell_valid(newx, y, z) &&
                          (is_cell_valid(newx, y, newz) || is_cell_valid(newx, newy, z))) ||
                         (is_cell_valid(x, newy, z) &&
                          (is_cell_valid(x, newy, newz) || is_cell_valid(newx, newy, z))) ||
                         (is_cell_valid(x, y, newz) &&
                          (is_cell_valid(x, newy, newz) || is_cell_valid(newx, y, newz)));
                         

        if (succ_valid && adj_valid) {
          float cost = sqrt(dX[i]*dX[i] + dY[i]*dY[i] + dZ[i]*dZ[i]);
          if (g + cost < astar.get_g(newx, newy, newz)) {
            astar.set_predecessor(newx, newy, newz, x, y, z, g + cost);
            if (astar.is_closed(newx, newy, newz)) {
              astar.make_incons(newx, newy, newz);
            } else {
              astar.open(newx, newy, newz);
            }
          }
        }
      }

      auto stop_time = std::chrono::high_resolution_clock::now();
      auto millis = std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time).count();
      if (eps != start_eps && millis >= max_millis) {
        // Stop planning!
        return true;
      }
    } 

    if (astar.open_list_empty()) {
      std::cout << "NO PATH" << std::endl;
      return false;
    }

    astar.get_path(int_path);
    path.resize(int_path.size());
    for (int i = 0; i < path.size(); ++i) {
      path[i] = cell_to_pos(Map<Vector3i>(&int_path[i][0]));
    }

    auto stop_time = std::chrono::high_resolution_clock::now();
    auto millis = std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time).count();
    if (millis >= max_millis) {
      // Stop planning!
      return true;
    }
  }
  return true;
}
