#pragma once
#include <Eigen/Dense>
#include <memory>
#include <unsupported/Eigen/Splines>

using namespace Eigen;
typedef Matrix<double, Dynamic, Dynamic, RowMajor> RowMatrixXd;

class CppSpline3D {
  public:
    CppSpline3D(const Ref<const VectorXd> &knots, const Ref<const RowMatrixXd> &coeffs) : knots(knots), coeffs(coeffs) {
      min_knot = knots(0);
      max_knot = knots(knots.size() - 1);
      knot_range = max_knot - min_knot;
      mySpline = std::make_shared<Spline<double, 3, 3>>((knots - min_knot*VectorXd::Ones(knots.size()))/knot_range, coeffs.transpose());
    }

    virtual Vector3d operator()(double t) const {
      return (*mySpline)((t - min_knot)/knot_range);
    }

    virtual Vector3d derivatives(double t) const {
      return mySpline->derivatives((t - min_knot)/knot_range, 1).col(1)/knot_range;
    }

    const VectorXd &get_knots() const {
      return knots;
    }

    const RowMatrixXd &get_coeffs() const {
      return coeffs;
    }

  protected:
    double min_knot;
    double max_knot;
    double knot_range;
    std::shared_ptr<Spline<double, 3, 3>> mySpline;
    VectorXd knots;
    RowMatrixXd coeffs;
};
