#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/gmdm.h"
#include "mt_tsp_ros2/elongate_dubins_path.h"

namespace py = pybind11;

PYBIND11_MODULE(gmdm_wrapper, m) {
  m.def("random_gmdm_planner", &random_gmdm_planner);
  m.def("gmdm_planner_primitives", &gmdm_planner_primitives);
  m.def("csc_inverse_free_v2", &csc_inverse_free_v2);
  m.def("ccc_inverse_free_v2", &ccc_inverse_free_v2);
  m.def("csc_inverse", &csc_inverse);
  m.def("ccc_inverse", &ccc_inverse);
  m.def("check_elongation_possible", &check_elongation_possible);
  m.def("get_elongation_intervals", &get_elongation_intervals);
}
