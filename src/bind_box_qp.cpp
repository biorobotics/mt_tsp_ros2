#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/box-qp.hpp"

namespace py = pybind11;

PYBIND11_MODULE(box_qp, m) {
  py::class_<BoxQP>(m, "BoxQP")
    .def(py::init<const std::size_t>())
    .def("solve", &BoxQP::solve)
    ;

  py::class_<BoxQPSolution>(m, "BoxQPSolution")
    .def_readwrite("x", &BoxQPSolution::x)
    ;
}
