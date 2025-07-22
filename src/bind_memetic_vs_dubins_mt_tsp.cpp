#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>
#include "mt_tsp_ros2/memetic_dubins_ce_mt_tsp/memetic_alg_vs_dubins.h"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<double>)
PYBIND11_MAKE_OPAQUE(std::vector<int>)

PYBIND11_MODULE(memetic_vs_dubins_mt_tsp, m) {
  m.def("memetic_alg_vs_dubins_mt_tsp", &memetic_alg);
  m.def("repair_chromosome_vs_dubins_mt_tsp", py::overload_cast<Ref<MatrixXd>, const Ref<const RowMatrixXd>&, const std::vector<ExtendedCppSpline>&, const Ref<const Vector2d>&, double, const Ref<const VectorXd> &, Ref<Vector1d>, double, const MemeticAlgParams&, double, bool, Ref<VectorXl>, Ref<VectorXl>, Ref<VectorXl>>(&repair_chromosome));
  m.def("repair_chromosome_vs_dubins_mt_tsp", py::overload_cast<Ref<MatrixXd>, const Ref<const RowMatrixXd>&, const std::vector<ExtendedCppSpline>&, const Ref<const Vector2d>&, double, const Ref<const VectorXd> &, Ref<Vector1d>, double, const MemeticAlgParams&, double, Ref<RowMatrixXd>, bool, Ref<VectorXl>, Ref<VectorXl>, Ref<VectorXl>>(&repair_chromosome));

  py::bind_vector<std::vector<double>>(m, "VectorOfDoubles");
  py::bind_vector<std::vector<int>>(m, "VectorOfInts");
}
