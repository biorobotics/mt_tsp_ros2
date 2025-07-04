#pragma once
#include <Eigen/Dense>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

using namespace Eigen;
namespace py = pybind11;

class CppPPoly {
  public:
    CppPPoly(const py::array_t<double> &c, const Ref<const VectorXd> &x) : x(x) {
      auto c_unchecked = c.unchecked<3>();
      degree = c_unchecked.shape(0) - 1;
      num_intervals = c_unchecked.shape(1);
      output_dim = c_unchecked.shape(2);
      num_breakpoints = x.size();
      if (num_breakpoints != num_intervals + 1) {
        throw std::runtime_error("number of breakpoints should be number of intervals + 1");
      }
      // m is the degree of the coefficient
      for (int m = 0; m < degree + 1; ++m) {
        this->c.push_back(RowMatrixXd(num_intervals, output_dim));
        // i is the interval
        for (int i = 0; i < num_intervals; ++i) {
          // j is the component of the output
          for (int j = 0; j < output_dim; ++j) {
            this->c.back()(i, j) = c_unchecked(m, i, j);
          } 
        }
      }
    }

    virtual VectorXd operator()(double xp, int &i) const {
      VectorXd ret = VectorXd::Zero(output_dim);
      // Iterate over intervals
      bool found_interval = false;
      if (xp < x(0)) {
        i = 0;
        found_interval = true;
      } else if (xp > num_breakpoints) {
        i = num_intervals - 1;
        found_interval = true;
      } else {
        for (i = 0; i < num_intervals; ++i) {
          if (x(i) <= xp && xp <= x(i + 1)) {
            found_interval = true;
            break;
          }
        }
      }
      if (!found_interval) {
        throw std::runtime_error("Did not find interval");
      }

      for (int m = 0; m < degree + 1; ++m) {
        ret += c[m].row(i)*pow(xp - x(i), degree - m);
      }
      return ret;
    }

    virtual VectorXd operator()(double xp) const {
      int i;
      return this->operator()(xp, i);
    }

  private:
    std::vector<RowMatrixXd> c;
    VectorXd x;
    int degree;
    int num_intervals;
    int output_dim;
    int num_breakpoints;
};
