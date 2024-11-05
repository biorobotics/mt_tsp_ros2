#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/gurobi_wrapper.h"

namespace py = pybind11;

PYBIND11_MODULE(gurobi_wrapper, m) {
  m.def("solve_socp", &solve_socp);
}
