#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/memetic_dubins_ce_mt_tsp/memetic_alg.h"
#include "mt_tsp_ros2/extended_cpp_spline.h"
#include "mt_tsp_ros2/circular_trajectory.h"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<py::object>)
PYBIND11_MAKE_OPAQUE(std::vector<CppSpline>)
PYBIND11_MAKE_OPAQUE(std::vector<ExtendedCppSpline>)
PYBIND11_MAKE_OPAQUE(std::vector<CircularTrajectory>)
PYBIND11_MAKE_OPAQUE(std::vector<RowMatrixXd>)

PYBIND11_MODULE(memetic_dubins_ce_mt_tsp, m) {
  m.def("memetic_alg", &memetic_alg);
  m.def("get_selected_pts", &get_selected_pts);
  m.def("repair_chromosome", py::overload_cast<Ref<MatrixXd>, const Ref<const RowMatrixXd>&, const Ref<const VectorXd>&, const std::vector<CircularTrajectory>&, const Ref<const Vector2d>&, double, double, Ref<Vector1d>, bool, double, int&, const MemeticAlgParams&, bool>(&repair_chromosome));

  py::class_<MemeticAlgParams>(m, "MemeticAlgParams")
    .def(py::init<>())
    .def(py::init<double, double, int, int, double>())
    .def("get_mutation_prob", &MemeticAlgParams::get_mutation_prob)
    .def("get_local_search_gd_step_size", &MemeticAlgParams::get_local_search_gd_step_size)
    .def("get_local_search_num_samples", &MemeticAlgParams::get_local_search_num_samples)
    .def("get_Tlp", &MemeticAlgParams::get_Tlp)
    .def("get_repair_step_size", &MemeticAlgParams::get_repair_step_size)
    ;

  py::class_<CppSpline>(m, "CppSpline")
    .def(py::init<const Ref<const VectorXd>&, const Ref<const RowMatrixXd>&>())
    .def("__call__", &CppSpline::operator())
    .def("derivatives", &CppSpline::derivatives)
    ;

  py::class_<ExtendedCppSpline>(m, "ExtendedCppSpline")
    .def(py::init<const Ref<const VectorXd>&, const Ref<const RowMatrixXd>&, double, double>())
    .def("__call__", &ExtendedCppSpline::operator())
    .def("derivatives", &ExtendedCppSpline::derivatives)
    ;

  py::class_<CircularTrajectory>(m, "CircularTrajectory")
    .def(py::init<const Ref<const Vector2d>&, double, double, double>())
    .def("__call__", &CircularTrajectory::operator())
    .def("derivatives", &CircularTrajectory::derivatives)
    .def("get_center", &CircularTrajectory::get_center)
    .def("get_rad", &CircularTrajectory::get_rad)
    .def("get_omega", &CircularTrajectory::get_omega)
    .def("get_theta0", &CircularTrajectory::get_theta0)
    ;

  py::bind_vector<std::vector<py::object>>(m, "VectorOfPyObjects");
  py::bind_vector<std::vector<CppSpline>>(m, "VectorOfCppSplines");
  py::bind_vector<std::vector<ExtendedCppSpline>>(m, "VectorOfExtendedCppSplines");
  py::bind_vector<std::vector<CircularTrajectory>>(m, "VectorOfCircularTrajectories");
  py::bind_vector<std::vector<RowMatrixXd>>(m, "VectorOfMatrices");
}
