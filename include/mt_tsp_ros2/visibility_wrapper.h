#pragma once

#include <Eigen/Dense>
#include <memory>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Triangular_expansion_visibility_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel       Kernel;
typedef Kernel::Point_2                                         Point_2;
typedef Kernel::Segment_2                                       Segment_2;
typedef CGAL::Arr_segment_traits_2<Kernel>                      Traits_2;
typedef CGAL::Arrangement_2<Traits_2>                           Arrangement_2;
typedef Arrangement_2::Halfedge_const_handle                    Halfedge_const_handle;
typedef Arrangement_2::Face_handle                              Face_handle;
typedef CGAL::Triangular_expansion_visibility_2<Arrangement_2>  TEV;

using namespace Eigen;

typedef Matrix<long, Dynamic, 1> VectorXl;

class VisibilityWrapper {
  public:
    VisibilityWrapper(Ref<Matrix<double, Dynamic, Dynamic, RowMajor>> boundary_vertices, Ref<Matrix<double, Dynamic, Dynamic, RowMajor>> obstacle_vertices, Ref<VectorXl> obstacle_start_indices);
    MatrixXd visibility_polygon(long obstacle_vertex_idx);
    MatrixXd visibility_polygon_interior(Ref<Vector2d> pt);
  private:
    Arrangement_2 env;
    std::unique_ptr<TEV> tev;
    MatrixXd obstacle_vertices;
    VectorXl obstacle_start_indices;
    std::vector<int> vertex_to_obstacle_ptr;
};
