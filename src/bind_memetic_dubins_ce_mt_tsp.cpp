#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/memetic_dubins_ce_mt_tsp/memetic_alg.h"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<py::object>)

PYBIND11_MODULE(memetic_dubins_ce_mt_tsp, m) {
  m.def("memetic_alg", &memetic_alg);

  py::bind_vector<std::vector<py::object>>(m, "VectorOfPyObjects");
}
