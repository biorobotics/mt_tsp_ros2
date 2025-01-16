#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/implicit_grid_astar_problem.h"

namespace py = pybind11;

PYBIND11_MODULE(implicit_grid_astar, m) {
  m.def("solve_implicit_grid_astar_problem", &solve_implicit_grid_astar_problem);
}
