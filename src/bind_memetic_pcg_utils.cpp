#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/memetic_dubins_ce_mt_tsp/memetic_pcg_utils.h"

namespace py = pybind11;

PYBIND11_MODULE(memetic_pcg_utils, m) {
  py::class_<MemeticPCGUtils>(m, "MemeticPCGUtils")
    .def(py::init<const Ref<const RowMatrixXd> &, const std::vector<ExtendedCppSpline> &, const Ref<const Vector2d> &, double, double, double, bool, bool, double, int>())
    .def("crossover", &MemeticPCGUtils::crossover)
    ;
}
