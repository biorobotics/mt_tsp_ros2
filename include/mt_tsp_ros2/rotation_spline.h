#pragma once
#include "mt_tsp_ros2/cpp_ppoly.h"

class RotationSpline {
  public:
    RotationSpline(const CppPPoly &ppoly, const py::array_t<double> &rotations) : ppoly(ppoly) {
      auto rotations_unchecked = rotations.unchecked<3>();
      int num_rotations = rotations_unchecked.shape(0);
      if (rotations_unchecked.shape(1) != 3 || rotations_unchecked.shape(2) != 3) {
        throw std::runtime_error("Rotations not correct size, should be (n, 3, 3)");
      }

      for (int rotation_idx = 0; rotation_idx < num_rotations; ++rotation_idx) {
        this->rotations.push_back(Matrix3d());
        for (int row = 0; row < 3; ++row) {
          for (int col = 0; col < 3; ++col) {
            this->rotations.back()(row, col) = rotations_unchecked(rotation_idx, row, col);
          }
        }
      }
    }

    virtual Matrix3d operator()(double t) {
      int interval;
      Vector3d rotvec = ppoly(t, interval);
      double angle = rotvec.norm();
      Vector3d axis(0, 0, 0);
      if (angle != 0) {
        axis = rotvec/angle;
      }
      return rotations[interval]*AngleAxisd(angle, axis).toRotationMatrix();
    }

  private:
    CppPPoly ppoly;
    std::vector<Matrix3d> rotations;
};
