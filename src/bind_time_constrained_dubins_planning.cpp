#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/time_constrained_dubins_planning/time_constrained_dubins_planner.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/no_time_dubins_planner.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/time_constrained_dubins_planner_forward_prop.h"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<RowMatrixXd>)

PYBIND11_MODULE(time_constrained_dubins_planning, m) {
  py::class_<TimeConstrainedDubinsPlanner>(m, "TimeConstrainedDubinsPlanner")
    .def(py::init<double, double, RowMatrixXbRef_const, Vector2dRef_const, Vector2dRef_const, bool, bool>())
    .def(py::init<double, double, RowMatrixXbRef_const, RowMatrixXdRef_const, Vector2dRef_const, Vector2dRef_const, bool, bool, bool>())
    .def("plan", &TimeConstrainedDubinsPlanner::plan)
    .def("is_state_valid", &TimeConstrainedDubinsPlanner::is_state_valid)
    .def("sample_random_state", &TimeConstrainedDubinsPlanner::sample_random_state)
    .def("checkMotion", &TimeConstrainedDubinsPlanner::checkMotion)
    .def("checkMotionForwardBackward", &TimeConstrainedDubinsPlanner::checkMotionForwardBackward)
    .def("get_nearest_neighbor_time", &TimeConstrainedDubinsPlanner::get_nearest_neighbor_time)
    .def("get_path_elongation_time", &TimeConstrainedDubinsPlanner::get_path_elongation_time)
    .def("get_path_elongation_check_time", &TimeConstrainedDubinsPlanner::get_path_elongation_check_time)
    .def("get_collision_check_time", &TimeConstrainedDubinsPlanner::get_collision_check_time)
    .def("get_sampling_time", &TimeConstrainedDubinsPlanner::get_sampling_time)
    .def("get_add_to_tree_time", &TimeConstrainedDubinsPlanner::get_add_to_tree_time)
    .def("get_num_tree_nodes", &TimeConstrainedDubinsPlanner::get_num_tree_nodes)
    .def("get_num_discarded_samples", &TimeConstrainedDubinsPlanner::get_num_discarded_samples)
    .def("get_num_samples", &TimeConstrainedDubinsPlanner::get_num_samples)
    ;

  py::class_<NoTimeDubinsPlanner>(m, "NoTimeDubinsPlanner")
    .def(py::init<double, RowMatrixXbRef_const, Vector2dRef_const, Vector2dRef_const>())
    .def("plan", &NoTimeDubinsPlanner::plan)
    .def("is_state_valid", &NoTimeDubinsPlanner::is_state_valid)
    .def("checkMotion", &NoTimeDubinsPlanner::checkMotion)
    .def("checkMotion", &NoTimeDubinsPlanner::checkMotion)
    ;

  py::bind_vector<std::vector<RowMatrixXd>>(m, "VectorOfMatrices");

  py::class_<TimeConstrainedDubinsPlannerForwardProp>(m, "TimeConstrainedDubinsPlannerForwardProp")
    .def(py::init<double, double, RowMatrixXbRef_const, RowMatrixXdRef_const, Vector2dRef_const, Vector2dRef_const, bool>())
    .def("plan", &TimeConstrainedDubinsPlannerForwardProp::plan)
    .def("get_path_elongation_time", &TimeConstrainedDubinsPlannerForwardProp::get_path_elongation_time)
    .def("get_path_elongation_check_time", &TimeConstrainedDubinsPlannerForwardProp::get_path_elongation_check_time)
    .def("get_collision_check_time", &TimeConstrainedDubinsPlannerForwardProp::get_collision_check_time)
    ;
}
