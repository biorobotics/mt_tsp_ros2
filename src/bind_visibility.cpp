#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/cgal_visibility_example.h"
#include "mt_tsp_ros2/visibility_wrapper.h"

namespace py = pybind11;

PYBIND11_MODULE(visibility, m) {
  m.def("cgal_visibility_example", &cgal_visibility_example);
  py::class_<VisibilityWrapper>(m, "VisibilityWrapper")
    .def(py::init<Ref<Matrix<double, Dynamic, Dynamic, RowMajor>>, Ref<Matrix<double, Dynamic, Dynamic, RowMajor>>, Ref<VectorXl>>())
    .def("visibility_polygon", &VisibilityWrapper::visibility_polygon)
    .def("visibility_polygon_interior", &VisibilityWrapper::visibility_polygon_interior)
    ;
}
