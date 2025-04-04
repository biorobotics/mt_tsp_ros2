#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/astar_planner_3d_wrapper.h"

namespace py = pybind11;

PYBIND11_MODULE(astar_planner_3d_wrapper, m) {
  py::class_<AStarPlanner3DWrapper>(m, "AStarPlanner3DWrapper")
    .def(py::init<const Ref<const VectorXd>&, const Ref<const Vector3d>&, const Ref<const Vector3d>&, const Ref<const Vector3i>&>())
    .def("plan", &AStarPlanner3DWrapper::plan)
    .def("set_connected_26", &AStarPlanner3DWrapper::set_connected_26)
    ;
}
