#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/memetic_dubins_ce_mt_tsp/memetic_alg_no_close_enough.h"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<double>)

PYBIND11_MODULE(memetic_dubins_mt_tsp, m) {
  m.def("memetic_alg_no_ce", &memetic_alg);
  m.def("get_selected_pts_no_ce", &get_selected_pts);
  m.def("repair_chromosome_no_ce", py::overload_cast<Ref<MatrixXd>, const Ref<const RowMatrixXd>&, const std::vector<ExtendedCppSpline>&, const Ref<const Vector2d>&, double, double, Ref<Vector1d>, bool, double, int&, const MemeticAlgParams&, bool, double>(&repair_chromosome));
  // m.def("repair_chromosome_no_ce", py::overload_cast<Ref<MatrixXd>, const Ref<const RowMatrixXd>&, const std::vector<ExtendedCppSpline>&, const Ref<const Vector2d>&, double, double, Ref<Vector1d>, bool, double, int&, const MemeticAlgParams&, bool, double, std::vector<double>&>(&repair_chromosome));

  py::bind_vector<std::vector<double>>(m, "VectorOfDoubles");
}
