#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "mt_tsp_ros2/vs_dubins_trj_through_seq_of_targets.h"

namespace py = pybind11;

PYBIND11_MODULE(vs_dubins_trj_through_seq_of_targets, m) {
  py::class_<VSDubinsTrjThroughSeqOfTargets>(m, "VSDubinsTrjThroughSeqOfTargets")
    .def(py::init<const Ref<const RowMatrixXd>&, const std::vector<ExtendedCppSpline>&, const Ref<const Vector2d>&, double, const Ref<const VectorXd>&, double, double, bool, bool>())
    .def("optimize_trj", &VSDubinsTrjThroughSeqOfTargets::optimize_trj)
    ;
}
