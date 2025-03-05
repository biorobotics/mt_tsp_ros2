#include "mt_tsp_ros2/visibility_wrapper.h"

VisibilityWrapper::VisibilityWrapper(Ref<Matrix<double, Dynamic, Dynamic, RowMajor>> boundary_vertices, Ref<Matrix<double, Dynamic, Dynamic, RowMajor>> obstacle_vertices, Ref<VectorXl> obstacle_start_indices) : obstacle_vertices(obstacle_vertices), obstacle_start_indices(obstacle_start_indices) {
  std::vector<Segment_2> segments;
  for (int i = 0; i < boundary_vertices.rows() - 1; ++i) {
    double x1 = boundary_vertices(i, 0);
    double y1 = boundary_vertices(i, 1);
    double x2 = boundary_vertices(i + 1, 0);
    double y2 = boundary_vertices(i + 1, 1);
    Point_2 p1(x1, y1), p2(x2, y2);
    segments.push_back(Segment_2(p1, p2));
  }

  for (int i = 0; i < obstacle_start_indices.rows(); ++i) {
    int start_idx = obstacle_start_indices[i];
    int end_idx;
    if (i == obstacle_start_indices.rows() - 1) {
      end_idx = obstacle_vertices.rows();
    } else {
      end_idx = obstacle_start_indices[i + 1];
    }
    for (int j = start_idx; j < end_idx - 1; ++j) {
      vertex_to_obstacle_ptr.push_back(i);
      double x1 = obstacle_vertices(j, 0);
      double y1 = obstacle_vertices(j, 1);
      double x2 = obstacle_vertices(j + 1, 0);
      double y2 = obstacle_vertices(j + 1, 1);
      Point_2 p1(x1, y1), p2(x2, y2);
      segments.push_back(Segment_2(p1, p2));
    }
    vertex_to_obstacle_ptr.push_back(i);
  }
  CGAL::insert_non_intersecting_curves(env,segments.begin(),segments.end());
  tev = std::make_unique<TEV>(env);
}

MatrixXd VisibilityWrapper::visibility_polygon(long obstacle_vertex_idx) {
  double x1 = obstacle_vertices(obstacle_vertex_idx, 0);
  double y1 = obstacle_vertices(obstacle_vertex_idx, 1);

  int i = vertex_to_obstacle_ptr[obstacle_vertex_idx];
  int start_idx = obstacle_start_indices[i];
  int end_idx;

  if (obstacle_vertex_idx == start_idx) {
    std::cout << "Do not compute visibility from first vertex of obstacle polygon" << std::endl;
    exit(1);
  }

  double x2 = obstacle_vertices(obstacle_vertex_idx - 1, 0);
  double y2 = obstacle_vertices(obstacle_vertex_idx - 1, 1);

  Point_2 p1(x1, y1), p2(x2, y2);

  Halfedge_const_handle he = env.halfedges_begin();
  // while (he->source()->point() != p1 || he->target()->point() != p2) {
  while (he != env.halfedges_end() && (he->target()->point() != p1 || he->source()->point() != p2)) {
    he++;
  }
  if (he == env.halfedges_end()) {
    std::cout << "Query point does not correspond to any halfedge" << std::endl;
    exit(1);
  }

  Arrangement_2 output_arr;
  Face_handle fh = tev->compute_visibility(p1, he, output_arr);
  Arrangement_2::Ccb_halfedge_circulator curr = fh->outer_ccb();
  std::vector<double> edge_coords;
  while (++curr != fh->outer_ccb()) {
    edge_coords.push_back(CGAL::to_double(curr->source()->point().x()));
    edge_coords.push_back(CGAL::to_double(curr->source()->point().y()));
    edge_coords.push_back(CGAL::to_double(curr->target()->point().x()));
    edge_coords.push_back(CGAL::to_double(curr->target()->point().y()));
  }
  double lastx = edge_coords[edge_coords.size() - 2];
  double lasty = edge_coords[edge_coords.size() - 1];
  edge_coords.push_back(lastx);
  edge_coords.push_back(lasty);
  edge_coords.push_back(edge_coords[0]);
  edge_coords.push_back(edge_coords[1]);

  return Map<Matrix<double, Dynamic, Dynamic, RowMajor>>(edge_coords.data(), edge_coords.size()/4, 4);
}

MatrixXd VisibilityWrapper::visibility_polygon_interior(Ref<Vector2d> pt) {
  Point_2 query_pt(pt(0), pt(1));

  // From https://github.com/d-krupke/pyvispoly/blob/main/src/pyvispoly/_cgal_bindings.cpp
  auto face = env.unbounded_face();
  auto hole_it = face->holes_begin();
  assert(hole_it != face->holes_end());
  auto f = (*hole_it)->twin()->face();
  if (f->is_unbounded()) {
    throw std::runtime_error("Bad arrangement. Face should not be unbounded.");
  }
  auto interior_face = f;

  Arrangement_2 output_arr;
  Face_handle fh = tev->compute_visibility(query_pt, interior_face, output_arr);
  Arrangement_2::Ccb_halfedge_circulator curr = fh->outer_ccb();
  std::vector<double> edge_coords;
  while (++curr != fh->outer_ccb()) {
    edge_coords.push_back(CGAL::to_double(curr->source()->point().x()));
    edge_coords.push_back(CGAL::to_double(curr->source()->point().y()));
    edge_coords.push_back(CGAL::to_double(curr->target()->point().x()));
    edge_coords.push_back(CGAL::to_double(curr->target()->point().y()));
  }
  double lastx = edge_coords[edge_coords.size() - 2];
  double lasty = edge_coords[edge_coords.size() - 1];
  edge_coords.push_back(lastx);
  edge_coords.push_back(lasty);
  edge_coords.push_back(edge_coords[0]);
  edge_coords.push_back(edge_coords[1]);
  return Map<Matrix<double, Dynamic, Dynamic, RowMajor>>(edge_coords.data(), edge_coords.size()/4, 4);
}
