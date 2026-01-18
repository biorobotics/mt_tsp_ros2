#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/memetic_dubins_ce_mt_tsp/memetic_pcg_utils_sampling.h"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<py::array_t<long>>)

PYBIND11_MODULE(memetic_pcg_utils_sampling, m) {
  py::bind_vector<std::vector<py::array_t<long>>>(m, "VectorOfLongArrays");

  py::class_<MemeticPCGUtilsSampling>(m, "MemeticPCGUtilsSampling")
    .def(py::init<int, double, int, int>())
    .def("crossover", &MemeticPCGUtilsSampling::crossover)
    ;
}
