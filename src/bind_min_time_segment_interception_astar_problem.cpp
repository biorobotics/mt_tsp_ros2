#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/min_time_segment_interception_astar_problem.h"

namespace py = pybind11;

PYBIND11_MODULE(mtvg_problem, m) {
  m.def("solve_min_time_segment_interception_astar_problem", &solve_min_time_segment_interception_astar_problem);
}
