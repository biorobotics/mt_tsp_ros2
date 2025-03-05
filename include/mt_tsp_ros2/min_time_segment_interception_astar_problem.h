#pragma once
#include "mt_tsp_ros2/astar_problem.h"
#include <memory>
#include "mt_tsp_ros2/arastar.h"

typedef Matrix<long, Dynamic, 1> VectorXl;

class MinTimeSegmentInterceptionAStarProblem : public AStarProblem {
  public:
    MinTimeSegmentInterceptionAStarProblem(const Ref<const VectorXl> &indptr,
                                           const Ref<const VectorXl> &indices, 
                                           const Ref<const VectorXd> &data, 
                                           const Ref<const VectorXl> &visible_vertices,
                                           const Ref<const Matrix<double, Dynamic, 2, RowMajor>> &visible_intervals,
                                           const Ref<const Matrix<double, Dynamic, 2, RowMajor>> &vertex_locs,
                                           const Ref<const Vector2d> &target_p0,
                                           const Ref<const Vector2d> &target_vel,
                                           double segment_t0, double segment_t1,
                                           double vmax_agent,
                                           int goal_node_idx,
                                           double start_time) : indptr(indptr), 
                                                                indices(indices), 
                                                                data(data),
                                                                visible_intervals(visible_intervals), 
                                                                vertex_locs(vertex_locs), 
                                                                target_p0(target_p0), 
                                                                target_vel(target_vel),
                                                                vmax_agent(vmax_agent),
                                                                goal_node_idx(goal_node_idx), 
                                                                start_time(start_time) {

      target_pos_at_start = target_p0 + target_vel*segment_t0;
      target_pos_at_end = target_p0 + target_vel*segment_t1;
      target_vel_normalized = target_vel.normalized();
      target_speed = target_vel.norm();

      if (visible_vertices(0) != -1) {
        for (int i = 0; i < visible_vertices.size(); ++i) {
          if (visible_vertices_map.find(visible_vertices[i]) == visible_vertices_map.end()) {
            visible_vertices_map[visible_vertices[i]] = std::vector<int>();
          }
          visible_vertices_map[visible_vertices[i]].push_back(i);
        }
      }
    }

    virtual void generate_successors(std::vector<AStarCell> &succ, std::vector<double> &transition_costs, const AStarNodePtr &node) const override {
      succ.clear();
      transition_costs.clear();

      int node_idx = node->get_cell()(0);
      if (node_idx == goal_node_idx) {
        return;
      }

      VectorXi next_cell(1);

      // If this node is visible from the target segment, find min-time to target segment
      auto it = visible_vertices_map.find(node_idx);
      if (it != visible_vertices_map.end()) {
        double t = start_time + node->get_g();

        const std::vector<int> &visible_vertex_indices = it->second;
        for (auto visible_vertex_idx : visible_vertex_indices) {
          // Feasibility check
          Vector2d pos_at_end = target_p0 + target_vel*visible_intervals(visible_vertex_idx, 1);
          double travel_time = (pos_at_end - vertex_locs.row(node_idx).transpose()).norm()/vmax_agent;
          double time_to_end = t + travel_time;
          if (time_to_end > visible_intervals(visible_vertex_idx, 1)) {
            // Infeasible
            continue;
          }

          // Special case check
          Vector2d pos_at_start = target_p0 + target_vel*visible_intervals(visible_vertex_idx, 0);
          travel_time = (pos_at_start - vertex_locs.row(node_idx).transpose()).norm()/vmax_agent;
          double time_to_start = t + travel_time;
          if (time_to_start <= visible_intervals(visible_vertex_idx, 0)) {
            // Min-time is at the beginning of the segment
            next_cell(0) = goal_node_idx;
            succ.push_back(next_cell);
            transition_costs.push_back(visible_intervals(visible_vertex_idx, 0) - t);
            continue;
          }

          // Quadratic formula from C*
          double Bp = vmax_agent*vmax_agent;
          Vector2d Cip = target_p0 - vertex_locs.row(node_idx).transpose();
          double Cp = Cip.dot(target_vel);
          double Ap = -vmax_agent*vmax_agent;
          double Dp = 0.;
          double Ep = Cip.squaredNorm();

          double A = target_vel.squaredNorm() - vmax_agent*vmax_agent;
          double B = 2*Bp*t + 2*Cp;
          double C = Ap*t*t - Dp*t + Ep;
          double tnext;
          if (A == 0) {
            // Special case where target moves at agent's max speed
            tnext = -C/B;
          } else {
            double discrim = sqrt(B*B - 4*A*C);
            tnext = (-B + discrim)/(2*A);
            if (tnext < t) {
              tnext = (-B - discrim)/(2*A);
            }
          }

          next_cell(0) = goal_node_idx;
          succ.push_back(next_cell);
          transition_costs.push_back(tnext - t);
        }
      }

      for (int neighbor_idx = indptr[node_idx]; neighbor_idx < indptr[node_idx + 1]; ++neighbor_idx) {
        int next_node_idx = indices[neighbor_idx];
        next_cell(0) = next_node_idx;
        succ.push_back(next_cell);
        transition_costs.push_back(data[neighbor_idx]);
      }
    }

    virtual bool is_goal(const AStarNodePtr &node) const override {
      return node->get_cell()(0) == goal_node_idx;
    }

    virtual double heuristic(const AStarCell &cell) const override {
      int node_idx = cell(0);
      if (node_idx == goal_node_idx) {
        return 0;
      }

      Vector2d rel_pos = vertex_locs.row(node_idx).transpose() - target_pos_at_start;
      double dot_product = target_vel_normalized.dot(rel_pos);
      if (target_speed == 0) {
        return rel_pos.norm()/vmax_agent;
      } else {
        return (rel_pos - target_vel_normalized*dot_product).norm()/vmax_agent;
      }
    }

  private:
    const Ref<const VectorXl> &indptr;
    const Ref<const VectorXl> &indices;
    const Ref<const VectorXd> &data;
    std::unordered_map<int, std::vector<int>> visible_vertices_map;
    const Ref<const Matrix<double, Dynamic, 2, RowMajor>> &visible_intervals;

    const Ref<const Matrix<double, Dynamic, 2, RowMajor>> &vertex_locs;
    const Ref<const Vector2d> &target_p0;
    const Ref<const Vector2d> &target_vel;
    double segment_t0;
    double segment_t1;
    double vmax_agent;
    int goal_node_idx;
    double start_time;

    Vector2d target_pos_at_start;
    Vector2d target_pos_at_end;
    Vector2d target_vel_normalized;
    double target_speed;
};

VectorXl solve_min_time_segment_interception_astar_problem(Ref<Matrix<double, 1, 1>> cost,
                                                           int start_node_idx,
                                                           const Ref<const VectorXl> &indptr,
                                                           const Ref<const VectorXl> &indices, 
                                                           const Ref<const VectorXd> &data, 
                                                           const Ref<const VectorXl> &visible_vertices,
                                                           const Ref<const Matrix<double, Dynamic, 2, RowMajor>> &visible_intervals,
                                                           const Ref<const Matrix<double, Dynamic, 2, RowMajor>> &vertex_locs,
                                                           const Ref<const Vector2d> &target_p0,
                                                           const Ref<const Vector2d> &target_vel,
                                                           double segment_t0, double segment_t1,
                                                           double vmax_agent,
                                                           int goal_node_idx,
                                                           double start_time) {
  std::shared_ptr<MinTimeSegmentInterceptionAStarProblem> problem =  std::make_shared<MinTimeSegmentInterceptionAStarProblem>(indptr,
                                                                                                                              indices, 
                                                                                                                              data, 
                                                                                                                              visible_vertices,
                                                                                                                              visible_intervals,
                                                                                                                              vertex_locs,
                                                                                                                              target_p0,
                                                                                                                              target_vel,
                                                                                                                              segment_t0, segment_t1,
                                                                                                                              vmax_agent,
                                                                                                                              goal_node_idx,
                                                                                                                              start_time);

  AStarCell start_cell = VectorXi(1);
  start_cell(0) = start_node_idx;
  std::shared_ptr<ARAStarData> astar_data = std::make_shared<ARAStarData>(problem, start_cell);
  AStarPath path;
  if (arastar(path, astar_data, 1.0, 10, 10000)) {
    VectorXl node_seq(path.size() - 1);
    for (int i = 0; i < path.size() - 1; ++i) {
      node_seq(i) = path[i](0);
    }
    cost(0) = astar_data->get_f_goal();
    return node_seq;
  } else {
    // std::cout << "No solution. Returning path where agent stays at start config" << std::endl;
    VectorXl node_seq(1);
    node_seq(0) = start_node_idx;
    cost(0) = std::numeric_limits<double>::infinity();
    return node_seq;
  }
}
