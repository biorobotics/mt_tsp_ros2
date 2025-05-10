#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/cpp_spline.h"

namespace py = pybind11;

PYBIND11_MODULE(cpp_spline, m) {
  py::class_<CppSpline>(m, "CppSpline")
    .def(py::init<const Ref<const VectorXd>&, const Ref<const RowMatrixXd>&>())
    .def("__call__", &CppSpline::operator())
    ;
}
