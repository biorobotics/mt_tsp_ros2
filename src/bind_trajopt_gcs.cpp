#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/mt_ixg_star.h"

namespace py = pybind11;

PYBIND11_MODULE(trajopt_gcs, m) {
  py::class_<GCSNode>(m, "GCSNode")
    .def(py::init<const Ref<const MatrixXdRM> &, const Ref<const MatrixXdRM> &, const Ref<const MatrixXdRM> &, const Ref<const MatrixXdRM> &, const Ref<const MatrixXdRM> &, const Ref<const MatrixXdRM> &, const Ref<const MatrixXdRM> &, const Ref<const MatrixXdRM> &>())
    .def("contains", &GCSNode::contains)
    .def("contains_pt", &GCSNode::contains_pt)
    ;

  py::class_<TrajoptThroughConvexSets>(m, "TrajoptThroughConvexSets")
    .def(py::init<int, const Ref<const VectorXd> &, const Ref<const VectorXd> &, const Ref<const VectorXd> &, double>())
    .def("add_convex_set", &TrajoptThroughConvexSets::add_convex_set)
    .def("solve_trajopt", &TrajoptThroughConvexSets::solve_trajopt)
    ;

  py::class_<MTIxGStar>(m, "MTIxGStar")
    .def(py::init<const std::vector<GCSNode> &, const Ref<const Matrix<bool, Dynamic, Dynamic, RowMajor>> &, int, const Ref<const VectorXd> &, const Ref<const VectorXd> &, double>())
    .def("solve", &MTIxGStar::solve)
    ;

  py::bind_vector<std::vector<VectorXd>>(m, "VectorOfVectors");
  py::bind_vector<std::vector<double>>(m, "VectorOfDoubles");
  py::bind_vector<std::vector<GCSNode>>(m, "VectorOfGCSNodes");
}
