#pragma once

#include "mt_tsp_ros2/rotation_spline.h"

class ConstantVelocityRotationTrajectory {
  public:
    ConstantVelocityRotationTrajectory() : R0(Matrix3d::Identity()), axis(0, 0, 1), w(0.) {
    }

    ConstantVelocityRotationTrajectory(const Ref<const Matrix3d> &R0, const Ref<const Vector3d> &axis, double w) : R0(R0), axis(axis), w(w) {
    }

    ConstantVelocityRotationTrajectory(const RotationSpline &spline, double t) {
      Matrix3d R = spline(t);
      double eps = 1e-4;
      Matrix3d Rplus = spline(t + eps);
      AngleAxisd Rdiff(Rplus*R.transpose());
      Rdiff.angle() /= eps;
      axis = Rdiff.axis();
      w = Rdiff.angle()/eps;
      R0 = AngleAxisd(w*t, axis).inverse()*R;

      Matrix3d R_test = this->operator()(t);
      Matrix3d Rdiff_test = R_test.transpose()*R;
      if (std::abs(Rdiff_test(0, 0) - 1) > 1e-10 ||
          std::abs(Rdiff_test(1, 1) - 1) > 1e-10 ||
          std::abs(Rdiff_test(2, 2) - 1) > 1e-10) {
        throw std::runtime_error("Linearization not correct");
      }
    }

    virtual Matrix3d operator()(double t) const {
      return (AngleAxisd(w*t, axis)*R0).toRotationMatrix();
    }
  private:
    AngleAxisd R0;
    Vector3d axis;
    double w;
};
