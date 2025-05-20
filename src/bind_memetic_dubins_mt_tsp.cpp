#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/memetic_dubins_ce_mt_tsp/memetic_alg_no_close_enough.h"

namespace py = pybind11;

PYBIND11_MODULE(memetic_dubins_mt_tsp, m) {
  m.def("memetic_alg_no_ce", &memetic_alg);
}
