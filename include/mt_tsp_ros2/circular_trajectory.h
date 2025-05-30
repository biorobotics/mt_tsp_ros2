#pragma once
#include <Eigen/Dense>
#include <memory>

using namespace Eigen;

class CircularTrajectory {
  public:
    CircularTrajectory(const Ref<const Vector2d> &center, double rad, double omega, double theta0) : center(center), rad(rad), omega(omega), theta0(theta0) {
    }

    virtual Vector2d operator()(double t) const {
      double theta = theta0 + omega*t;
      return Vector2d(center(0) + rad*cos(theta), center(1) + rad*sin(theta));
    }

    virtual Vector2d derivatives(double t) const {
      double theta = theta0 + omega*t;
      return Vector2d(-rad*sin(theta)*omega, rad*cos(theta)*omega);
    }

    Vector2d get_center() {
      return center;
    }

    double get_rad() {
      return rad;
    }

    double get_omega() {
      return omega;
    }

    double get_theta0() {
      return theta0;
    }

  protected:
    Vector2d center;
    double rad;
    double omega;
    double theta0;
};
