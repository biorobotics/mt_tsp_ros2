#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/elongate_one_sided_dubins_path.h"

namespace py = pybind11;

PYBIND11_MODULE(one_sided_dubins_elongation, m) {
  m.def("turns_for_CS_path", &turns_for_CS_path);
  m.def("turns_for_LR_path", &turns_for_LR_path);
  m.def("turns_for_one_sided_dubins_path", &turns_for_one_sided_dubins_path);
  m.def("elongated_dubins_path_one_sided", &elongated_dubins_path_one_sided);
}
