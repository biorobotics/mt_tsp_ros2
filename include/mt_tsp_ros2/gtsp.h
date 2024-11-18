#include <Eigen/Dense>
#include <string>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

using namespace Eigen;
namespace py = pybind11;

typedef Ref<Matrix<long, Dynamic, 1>> VectorXlRef;
typedef const Ref<const Matrix<long, Dynamic, 1>> &VectorXlRef_const;
typedef const Ref<const Matrix<double, Dynamic, Dynamic, RowMajor>> &RowMatrixXdRef_const;
typedef Matrix<double, Dynamic, Dynamic, RowMajor> RowMatrixXd;
typedef Ref<Matrix<double, Dynamic, Dynamic, RowMajor>> RowMatrixXdRef;
typedef Matrix<int, Dynamic, Dynamic, RowMajor> RowMatrixXi;

VectorXd pcg_gtsp(VectorXlRef node_seq, py::object send_callback, py::object recv_callback, RowMatrixXdRef_const cost_mat, VectorXlRef_const group_start_idx, bool verbose, double time_limit, std::string save_path, double mipgap, VectorXlRef_const known_feas_tour, bool solve_relaxed, int proc_idx);

VectorXd solve_gtsp_no_gsec(VectorXlRef node_seq, RowMatrixXdRef_const cost_mat, VectorXlRef_const group_start_idx, bool verbose, double time_limit, std::string save_path, bool take_first_feas_soln, double mipgap, bool solve_relaxed);

VectorXd solve_gtsp_no_gsec_lazy_edge_eval(VectorXlRef node_seq, py::object edge_evaluator, RowMatrixXdRef_const lb_cost_mat, VectorXlRef_const group_start_idx, bool verbose, double time_limit, std::string save_path, bool take_first_feas_soln, double mipgap, VectorXlRef_const known_feas_tour);
