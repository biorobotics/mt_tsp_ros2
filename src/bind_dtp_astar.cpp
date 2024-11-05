#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/dtp_astar_problem.h"

namespace py = pybind11;

PYBIND11_MODULE(dtp_astar, m) {
  m.def("solve_dtp_astar_problem", &solve_dtp_astar_problem);

  py::class_<DTPAStarSolver>(m, "DTPAStarSolver")
    .def(py::init<const Ref<const RowMatrixXd>&>())
    .def("solve", &DTPAStarSolver::solve)
    ;
}
