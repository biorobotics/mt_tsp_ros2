#pragma once
#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;

typedef Matrix<long, Dynamic, 1> VectorXl;
typedef Matrix<double, Dynamic, Dynamic, RowMajor> MatrixXdRM;

// Represents convex set Ax + b <= 0
struct GCSNode {
  GCSNode(const Ref<const MatrixXd> &A_eq, 
          const Ref<const VectorXd> &b_eq,
          const Ref<const MatrixXd> &A_ineq, 
          const Ref<const VectorXd> &b_ineq,
          const Ref<const MatrixXd> &A_eq_pt, 
          const Ref<const VectorXd> &b_eq_pt,
          const Ref<const MatrixXd> &A_ineq_pt, 
          const Ref<const VectorXd> &b_ineq_pt,
          bool in_cpp) : A_eq(A_eq),
                         b_eq(b_eq),
                         A_ineq(A_ineq),
                         b_ineq(b_ineq),
                         A_eq_pt(A_eq_pt),
                         b_eq_pt(b_eq_pt),
                         A_ineq_pt(A_ineq_pt),
                         b_ineq_pt(b_ineq_pt) {
  }

  GCSNode(const Ref<const MatrixXdRM> &A_eq, 
          const Ref<const VectorXd> &b_eq,
          const Ref<const MatrixXdRM> &A_ineq, 
          const Ref<const VectorXd> &b_ineq,
          const Ref<const MatrixXdRM> &A_eq_pt, 
          const Ref<const VectorXd> &b_eq_pt,
          const Ref<const MatrixXdRM> &A_ineq_pt, 
          const Ref<const VectorXd> &b_ineq_pt) : A_eq(A_eq),
                                                  b_eq(b_eq),
                                                  A_ineq(A_ineq),
                                                  b_ineq(b_ineq),
                                                  A_eq_pt(A_eq_pt),
                                                  b_eq_pt(b_eq_pt),
                                                  A_ineq_pt(A_ineq_pt),
                                                  b_ineq_pt(b_ineq_pt) {
  }

  bool contains(const Ref<const VectorXd> &x) {
    if (!std::isnan(A_eq(0, 0)) && (A_eq*x + b_eq).lpNorm<Infinity>() > 1e-4) {
      return false;
    }
    if (!std::isnan(A_ineq(0, 0)) && (A_ineq*x + b_ineq).maxCoeff() > 1e-4) {
      return false;
    }
    return true;
  }

  bool contains_pt(const Ref<const VectorXd> &x) {
    if (!std::isnan(A_eq_pt(0, 0)) && (A_eq_pt*x + b_eq_pt).lpNorm<Infinity>() > 1e-4) {
      return false;
    }
    if (!std::isnan(A_ineq_pt(0, 0)) && (A_ineq_pt*x + b_ineq_pt).maxCoeff() > 1e-4) {
      return false;
    }
    return true;
  }

  MatrixXd A_eq;
  VectorXd b_eq;
  MatrixXd A_ineq;
  VectorXd b_ineq;
  MatrixXd A_eq_pt;
  VectorXd b_eq_pt;
  MatrixXd A_ineq_pt;
  VectorXd b_ineq_pt;
};
