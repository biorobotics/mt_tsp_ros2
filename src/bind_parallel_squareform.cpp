#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "mt_tsp_ros2/parallel_squareform.h"

namespace py = pybind11;

PYBIND11_MODULE(parallel_squareform, m) {
  m.def("parallel_squareform", &parallel_squareform);
}
