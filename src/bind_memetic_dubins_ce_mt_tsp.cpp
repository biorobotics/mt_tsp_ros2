#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/memetic_dubins_ce_mt_tsp/memetic_alg.h"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<py::object>)
PYBIND11_MAKE_OPAQUE(std::vector<CppSpline>)
PYBIND11_MAKE_OPAQUE(std::vector<RowMatrixXd>)

PYBIND11_MODULE(memetic_dubins_ce_mt_tsp, m) {
  m.def("memetic_alg", &memetic_alg);
  m.def("repair_chromosome", py::overload_cast<Ref<MatrixXd>, const Ref<const RowMatrixXd>&, const Ref<const VectorXd>&, const std::vector<CppSpline>&, const Ref<const Vector2d>&, double, double, Ref<Vector1d>, bool, double, int&, const MemeticAlgParams&>(&repair_chromosome));

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
    ;

  py::bind_vector<std::vector<py::object>>(m, "VectorOfPyObjects");
  py::bind_vector<std::vector<CppSpline>>(m, "VectorOfCppSplines");
  py::bind_vector<std::vector<RowMatrixXd>>(m, "VectorOfMatrices");
}
