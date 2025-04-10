#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "mt_tsp_ros2/nearest_neighbor_edge_posterior.h"

namespace py = pybind11;

PYBIND11_MODULE(nearest_neighbor_edge_posterior, m) {
  py::class_<NearestNeighborEdgePosterior>(m, "NearestNeighborEdgePosterior")
    .def(py::init<double, double, bool>())
    .def("sample", &NearestNeighborEdgePosterior::sample)
    .def("add_edge", &NearestNeighborEdgePosterior::add_edge)
    ;
}
