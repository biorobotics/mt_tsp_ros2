#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/reopt_gtsp_tour.h"

namespace py = pybind11;

PYBIND11_MODULE(reopt_gtsp_tour, m) {
  m.def("reopt_gtsp_tour", &reopt_gtsp_tour);
}
