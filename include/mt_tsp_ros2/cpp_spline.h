#pragma once
#include <Eigen/Dense>
#include <memory>
#include <unsupported/Eigen/Splines>

using namespace Eigen;
typedef Matrix<double, Dynamic, Dynamic, RowMajor> RowMatrixXd;

namespace py = pybind11;

class CppSpline {
  public:
    CppSpline(const Ref<const VectorXd> &knots, const Ref<const RowMatrixXd> &coeffs) {
      min_knot = knots(0);
      max_knot = knots(knots.size() - 1);
      knot_range = max_knot - min_knot;
      mySpline = std::make_shared<Spline<double, 2, 3>>((knots - min_knot*VectorXd::Ones(knots.size()))/knot_range, coeffs.transpose());
    }

    Vector2d operator()(double t) const {
      return (*mySpline)((t - min_knot)/knot_range);
    }

    Vector2d derivatives(double t) const {
      return mySpline->derivatives((t - min_knot)/knot_range, 1)/knot_range;
    }

  private:
    double min_knot;
    double max_knot;
    double knot_range;
    std::shared_ptr<Spline<double, 2, 3>> mySpline;
};
