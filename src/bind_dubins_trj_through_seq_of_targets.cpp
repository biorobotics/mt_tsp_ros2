#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "mt_tsp_ros2/dubins_trj_through_seq_of_targets.h"

namespace py = pybind11;

PYBIND11_MODULE(dubins_trj_through_seq_of_targets, m) {
  py::class_<DubinsTrjThroughSeqOfTargets>(m, "DubinsTrjThroughSeqOfTargets")
    .def(py::init<const Ref<const RowMatrixXd>&, const std::vector<ExtendedCppSpline>&, const Ref<const Vector2d>&, double, double, double, bool, double, bool, bool, double>())
    .def("optimize_trj", &DubinsTrjThroughSeqOfTargets::optimize_trj)
    ;
}
