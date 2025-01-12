#include <Eigen/Dense>
#include <string>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

using namespace Eigen;
namespace py = pybind11;

typedef Ref<Matrix<long, Dynamic, 1>> VectorXlRef;
typedef Ref<Matrix<double, Dynamic, 1>> VectorXdRef;
typedef const Ref<const Matrix<long, Dynamic, 1>> VectorXlRef_const;
typedef const Ref<const Matrix<double, Dynamic, 1>> VectorXdRef_const;

VectorXd solve_socp(VectorXdRef soln, VectorXlRef_const A_csr_indptr, VectorXlRef_const A_csr_indices, VectorXdRef_const A_csr_data, VectorXdRef_const b, VectorXlRef_const Prows, VectorXlRef_const Pcols, VectorXdRef_const Pvals, VectorXdRef_const gradient, int num_zero_cone, int num_linear_cone, bool add_soc, int num_ctrl_pts, int li_dim, int l_idx, int dim_q, int vars_per_step, int steps, int num_decision_vars, bool verbose, double time_limit);
