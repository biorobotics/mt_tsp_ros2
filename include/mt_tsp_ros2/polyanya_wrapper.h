#pragma once

#include <Eigen/Dense>
#include <memory>
#include "searchinstance.h"
#include "mesh.h"

using namespace Eigen;

class PolyanyaWrapper {
  public:
    PolyanyaWrapper(const Ref<const Matrix<bool, Dynamic, Dynamic, RowMajor>> &occupancy, const std::string &polyanya_path, const std::string &map_save_path);

    MatrixXd shortest_path(const Ref<const Vector2d>& start, const Ref<const Vector2d>& goal);
  private:
    std::shared_ptr<polyanya::Mesh> mesh;
    std::shared_ptr<polyanya::SearchInstance> si;
};

