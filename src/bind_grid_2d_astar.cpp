#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/grid_2d_astar_problem.h"

namespace py = pybind11;

PYBIND11_MODULE(grid_2d_astar, m) {
  m.def("solve_grid_2d_astar_problem", &solve_grid_2d_astar_problem);
}
