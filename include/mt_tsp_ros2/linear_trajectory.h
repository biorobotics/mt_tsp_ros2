#pragma once
#include "mt_tsp_ros2/cpp_spline_3d.h"

class LinearTrajectory {
  public:
    LinearTrajectory(const Ref<const VectorXd> &p0, const Ref<const VectorXd> &v) : p0(p0), v(v) {
    }

    // Linearize the spline at time t
    LinearTrajectory(const CppSpline3D &spline, double t) {
      VectorXd p = spline(t);
      v = spline.derivatives(t);
      p0 = p - v*t;
    }

    virtual VectorXd operator()(double t) const {
      return p0 + v*t;
    }

    virtual VectorXd derivatives(double t) const {
      return v;
    }

  protected:
    VectorXd p0;
    VectorXd v;
};
