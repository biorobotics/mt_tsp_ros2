#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/memetic_dubins_ce_mt_tsp/memetic_pcg_utils_sampling.h"

namespace py = pybind11;

PYBIND11_MODULE(memetic_pcg_utils_sampling, m) {
  py::class_<MemeticPCGUtils>(m, "MemeticPCGUtilsSampling")
    .def(py::init<int, double>())
    .def("crossover", &MemeticPCGUtilsSampling::crossover)
    ;
}
