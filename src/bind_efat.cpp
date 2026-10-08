#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "mt_tsp_ros2/efat.h"

namespace py = pybind11;

PYBIND11_MODULE(efat, m) {
  py::class_<EFAT>(m, "EFAT")
    .def(py::init<const Ref<const RowMatrixXd>&, const std::vector<ExtendedCppSpline>&, const Ref<const Vector2d>&, double, bool, double, bool>())
    .def("efat_chain", &EFAT::efat_chain)
    ;
}
