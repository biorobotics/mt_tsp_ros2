#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/memetic_dubins_ce_mt_tsp/memetic_alg_sampling.h"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<double>)
PYBIND11_MAKE_OPAQUE(std::vector<int>)
PYBIND11_MAKE_OPAQUE(std::vector<MatrixXd>)
PYBIND11_MAKE_OPAQUE(std::vector<py::array_t<double>>)

PYBIND11_MODULE(memetic_sampling, m) {
  m.def("memetic_alg_sampling", &memetic_alg_sampling);

  py::bind_vector<std::vector<double>>(m, "VectorOfDoubles");
  py::bind_vector<std::vector<int>>(m, "VectorOfInts");
  py::bind_vector<std::vector<MatrixXd>>(m, "VectorOfMatrices");
  py::bind_vector<std::vector<py::array_t<long>>>(m, "VectorOfLongArrays");
}
