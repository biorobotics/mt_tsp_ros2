#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/dag_dfs.h"
#include "mt_tsp_ros2/lifelong_dag_dfs.h"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<py::array_t<long>>)

PYBIND11_MODULE(dag_dfs, m) {
  m.def("dag_dfs", &dag_dfs);
  py::bind_vector<std::vector<py::array_t<long>>>(m, "VectorOfLongArrays");

  py::class_<LifelongDAGDFSPlanner>(m, "LifelongDAGDFSPlanner")
    .def(py::init<RowMatrixXdRef_const, VectorXlRef_const, const std::vector<py::array_t<long>> &, bool, bool, RowMatrixXdRef_const>())
    .def("plan", &LifelongDAGDFSPlanner::plan)
    .def("get_before_time", &LifelongDAGDFSPlanner::get_before_time)
    ;
}
