#pragma once
#include "mt_tsp_ros2/time_constrained_dubins_planning/no_time_dubins_planner.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/yao2020/ic_sliding.h"
#include <fstream>

const int num_turns_after_slide = 4;

class Yao2020Planner {
  public:
    Yao2020Planner(double vmax, double rho, RowMatrixXbRef_const occupancy, RowMatrixXdRef_const rects, Vector2dRef_const map_lb, Vector2dRef_const map_ub);

    bool check_path_collides(double &x_collision, double &y_collision, double &theta_collision,
                             int &collision_turn_idx, double &collision_dist,
                             double x_0, double y_0, double theta_0,
                             double x_f, double y_f, double theta_f,
                             RowMatrixXdRef_const turns);

    // path_seq is the sequence of paths (turn sequences) we got via sliding
    RowMatrixXd plan(VectorXdRef_const start, VectorXdRef_const goal, double time_limit, int max_iter, std::vector<RowMatrixXd> &path_seq, bool just_return_rrt_trj, bool ignore_obstacles_in_normalization);

    bool is_state_valid(VectorXdRef_const state_vec);
  private:
    NoTimeDubinsPlanner spatial_planner;
    double vmax;
    double rho;
    RowMatrixXb occupancy;
    RowMatrixXd rects;
    Vector2d map_lb;
    Vector2d map_ub;
    // bool debug = true;
    bool debug = false;
};
