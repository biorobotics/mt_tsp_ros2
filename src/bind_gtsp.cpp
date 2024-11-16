#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/gtsp.h"

namespace py = pybind11;

PYBIND11_MODULE(gtsp_wrapper, m) {
  m.def("solve_gtsp_no_gsec", &solve_gtsp_no_gsec);
  m.def("pcg_gtsp", &pcg_gtsp);
  m.def("solve_gtsp_no_gsec_lazy_edge_eval", &solve_gtsp_no_gsec_lazy_edge_eval);
}
