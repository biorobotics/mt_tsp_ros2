#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/continuous_dag_dfs_no_obstacles.h"

namespace py = pybind11;

PYBIND11_MODULE(continuous_dag_dfs_no_obstacles, m) {
  m.def("continuous_dag_dfs_no_obstacles", &continuous_dag_dfs_no_obstacles);
}
