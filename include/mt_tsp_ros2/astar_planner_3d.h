#pragma once

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <memory>

using namespace Eigen;
using namespace std;

class AStarPlanner3D {
  public:
    AStarPlanner3D(const Tensor<bool, 3> &occupancy, const Vector3d &map_lb, const Vector3d &map_ub, const Vector3i &ncells);

    // Each Eigen::VectorXd in the boxes vector should have length 6, containing the center of the box, then the size
    bool plan_r3_path(std::vector<Vector3d>& path, const Ref<const Vector3d>& start_pos, const Ref<const Vector3d>& goal_pos);

    void set_connected_26(bool connected26) {
      if (connected26) {
        max_successors = 26;
      } else {
        max_successors = 6;
      }
    }

  private:
    Vector3d map_lb;
    Vector3i ncells;
    Vector3d cell_size;
    Tensor<bool, 3> occupancy;

    bool is_cell_valid(int x, int y, int z);
    bool is_cell_valid(const Vector3i &cell);
    Vector3i pos_to_cell(const Ref<const Vector3d>& pos);
    Vector3d cell_to_pos(const Ref<const Vector3i>& cell);
    int max_successors;
};
