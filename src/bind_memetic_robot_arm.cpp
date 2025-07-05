#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/memetic_robot_arm_mt_tsp/robot_arm_nlp.h"
#include "mt_tsp_ros2/memetic_robot_arm_mt_tsp/memetic_alg_robot_arm.h"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<SE3Spline>)
PYBIND11_MAKE_OPAQUE(std::vector<int>)
PYBIND11_MAKE_OPAQUE(std::vector<double>)

PYBIND11_MODULE(memetic_robot_arm, m) {
  py::class_<RobotArmNLPSolver>(m, "RobotArmNLPSolver")
    .def(py::init<const Ref<const RowMatrixXd> &, const std::vector<SE3Spline> &, const Ref<const VectorXd> &, const Ref<const VectorXl> &, int, const Ref<const VectorXd> &, const Ref<const VectorXd> &, const Ref<const Vector4d> &, const Ref<const MatrixXd> &, int>())
    .def("solve", &RobotArmNLPSolver::solve)
    .def("get_Jrows", &RobotArmNLPSolver::get_Jrows)
    .def("get_Jcols", &RobotArmNLPSolver::get_Jcols)
    .def("get_Hrows", &RobotArmNLPSolver::get_Hrows)
    .def("get_Hcols", &RobotArmNLPSolver::get_Hcols)
    .def("get_num_constraints", &RobotArmNLPSolver::get_num_constraints)
    .def("get_num_decision_vars", &RobotArmNLPSolver::get_num_decision_vars)
    .def("set_warm_start", &RobotArmNLPSolver::set_warm_start)
    .def("get_dense_jacobian", &RobotArmNLPSolver::get_dense_jacobian)
    ;

  py::bind_vector<std::vector<SE3Spline>>(m, "VectorOfSE3Splines");
  py::bind_vector<std::vector<int>>(m, "VectorOfInts");
  py::bind_vector<std::vector<double>>(m, "VectorOfDoubles");

  m.def("memetic_alg", &memetic_alg);
  m.def("repair_chromosome", py::overload_cast<Ref<MatrixXd>, Ref<Vector1d>, const Ref<const RowMatrixXd> &, const std::vector<SE3Spline> &, const Ref<const VectorXd> &, const Ref<const VectorXd> &, const Ref<const VectorXd> &, const Ref<const Vector4d> &, const Ref<const RowMatrixXd> &, Ref<RowMatrixXd>, Ref<VectorXd>, int>(&repair_chromosome));
}
