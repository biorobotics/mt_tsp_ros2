#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/time_constrained_dubins_planning/time_constrained_dubins_planner.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/no_time_dubins_planner.h"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<RowMatrixXd>)

PYBIND11_MODULE(time_constrained_dubins_planning, m) {
  py::class_<TimeConstrainedDubinsPlanner>(m, "TimeConstrainedDubinsPlanner")
    .def(py::init<double, double, RowMatrixXbRef_const, Vector2dRef_const, Vector2dRef_const>())
    .def("plan", &TimeConstrainedDubinsPlanner::plan)
    .def("is_state_valid", &TimeConstrainedDubinsPlanner::is_state_valid)
    .def("sample_random_state", &TimeConstrainedDubinsPlanner::sample_random_state)
    .def("checkMotion", &TimeConstrainedDubinsPlanner::checkMotion)
    ;

  py::class_<NoTimeDubinsPlanner>(m, "NoTimeDubinsPlanner")
    .def(py::init<double, RowMatrixXbRef_const, Vector2dRef_const, Vector2dRef_const>())
    .def("plan", &NoTimeDubinsPlanner::plan)
    .def("is_state_valid", &NoTimeDubinsPlanner::is_state_valid)
    .def("checkMotion", &NoTimeDubinsPlanner::checkMotion)
    ;

  py::bind_vector<std::vector<RowMatrixXd>>(m, "VectorOfMatrices");
}
