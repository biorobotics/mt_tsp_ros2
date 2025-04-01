#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/cbss_data.h"
#include "mt_tsp_ros2/cbss.h"
#include "mt_tsp_ros2/solve_cbss_problem.h"

namespace py = pybind11;

PYBIND11_MODULE(mt_cbss, m) {
  m.def("solve_cbss_problem", &solve_cbss_problem);
}
