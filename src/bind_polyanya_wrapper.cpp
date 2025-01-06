#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/polyanya_wrapper.h"

namespace py = pybind11;

PYBIND11_MODULE(polyanya_wrapper, m) {
  py::class_<PolyanyaWrapper>(m, "PolyanyaWrapper")
    .def(py::init<const Ref<const Matrix<bool, Dynamic, Dynamic, RowMajor>>&, const std::string&, const std::string&>())
    .def("shortest_path", &PolyanyaWrapper::shortest_path)
    ;
}
