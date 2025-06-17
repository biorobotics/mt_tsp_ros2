#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/adjacency_list.h"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<VectorXi>)

PYBIND11_MODULE(adjacency_list, m) {
  py::class_<AdjacencyList>(m, "AdjacencyList")
    .def(py::init<>())
    .def("add_edge", &AdjacencyList::add_edge)
    .def("dfs", &AdjacencyList::dfs)
    .def("batch_dfs", &AdjacencyList::batch_dfs)
    ;

  py::bind_vector<std::vector<VectorXi>>(m, "VectorOfVectorXi");

  py::class_<VectorXi>(m, "VectorXi")
    .def(py::init<int>())
  ;
}
