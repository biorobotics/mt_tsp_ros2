#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/cpp_spline.h"
#include "mt_tsp_ros2/cpp_ppoly.h"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<py::array_t<double>>)

PYBIND11_MODULE(cpp_spline, m) {
  py::class_<CppSpline>(m, "CppSpline")
    .def(py::init<const Ref<const VectorXd>&, const Ref<const RowMatrixXd>&>())
    .def("__call__", &CppSpline::operator())
    ;
    
  py::class_<CppPPoly>(m, "CppPPoly")
    .def(py::init<const py::array_t<double> &, const Ref<const VectorXd> &>())
    .def("__call__", &CppPPoly::operator())
    ;

  py::bind_vector<std::vector<py::array_t<double>>>(m, "VectorOfDoubles");
}
