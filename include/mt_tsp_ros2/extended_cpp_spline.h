#pragma once
#include <Eigen/Dense>
#include <memory>
#include <unsupported/Eigen/Splines>
#include "mt_tsp_ros2/cpp_spline.h"

using namespace Eigen;
typedef Matrix<double, Dynamic, Dynamic, RowMajor> RowMatrixXd;

namespace py = pybind11;

class ExtendedCppSpline : public CppSpline {
  public:
    ExtendedCppSpline(const Ref<const VectorXd> &knots, const Ref<const RowMatrixXd> &coeffs, double tw_start, double tw_end) : CppSpline(knots, coeffs), tw_start(tw_start), tw_end(tw_end) {
      pos_at_start = CppSpline::operator()(tw_start);
      pos_at_end = CppSpline::operator()(tw_end);
      vel_at_start = CppSpline::derivatives(tw_start);
      vel_at_end = CppSpline::derivatives(tw_end);
    }

    virtual Vector2d operator()(double t) const override {
      if (t < tw_start) {
        return pos_at_start + vel_at_start*(t - tw_start);
      }
      if (t > tw_end) {
        return pos_at_end + vel_at_end*(t - tw_end);
      }
      return (*mySpline)((t - min_knot)/knot_range);
    }

    virtual Vector2d derivatives(double t) const override {
      if (t < tw_start) {
        return vel_at_start;
      }
      if (t > tw_end) {
        return vel_at_end;
      }
      return mySpline->derivatives((t - min_knot)/knot_range, 1).col(1)/knot_range;
    }

  protected:
    double tw_start;
    double tw_end;
    Vector2d pos_at_start;
    Vector2d vel_at_start;
    Vector2d pos_at_end;
    Vector2d vel_at_end;
};
