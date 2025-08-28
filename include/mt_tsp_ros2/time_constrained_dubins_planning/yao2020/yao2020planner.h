#pragma once
#include "mt_tsp_ros2/time_constrained_dubins_planning/no_time_dubins_planner.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/yao2020/ic_sliding.h"

const int num_turns_after_slide = 4;

class Yao2020Planner {
  public:
    Yao2020Planner(double vmax, double rho, RowMatrixXbRef_const occupancy, RowMatrixXdRef_const rects, Vector2dRef_const map_lb, Vector2dRef_const map_ub) : spatial_planner(rho, occupancy, rects, map_lb, map_ub), vmax(vmax), rho(rho) {
    }

    bool check_path_collides(double &x_collision, double &y_collision, double &theta_collision,
                             int &collision_turn_idx, double &collision_dist,
                             double x_0, double y_0, double theta_0,
                             double x_f, double y_f, double theta_f,
                             RowMatrixXdRef_const turns) {
      double x = x_0;
      double y = y_0;
      double theta = theta_0;
      double ctheta = cos(theta_0);
      double stheta = sin(theta_0);

      double next_x;
      double next_y;
      double next_theta;
      double next_ctheta;
      double next_stheta;
      for (int turn_idx = 0; turn_idx < turns.rows(); ++turn_idx) {
        double turn_dir = turns(turn_idx, 0);
        double turn_dist = turns(turn_idx, 1);
        double rho_times_turn_dir = rho*turn_dir;

        if (turn_dir == 0) {
          // S segment
          next_theta = theta;
          next_ctheta = ctheta;
          next_stheta = stheta;
          next_x = x + turn_dist*ctheta;
          next_y = y + turn_dist*stheta;
        } else {
          // C segment
          next_theta = theta + turn_dist/rho_times_turn_dir;
          next_ctheta = cos(next_theta);
          next_stheta = sin(next_theta);
          next_x = x + rho_times_turn_dir*(-stheta + next_stheta);
          next_y = y + rho_times_turn_dir*(ctheta - next_ctheta);
        }

        bool collision = !spatial_planner.collision_free_get_intersection_point(x, y, theta, turn_dir, turn_dist, next_x, next_y, x_collision, y_collision, theta_collision, collision_dist);
        if (collision) {
          collision_turn_idx = turn_idx;
          return true;
        }
        x = next_x;
        y = next_y;
        theta = next_theta;
        ctheta = next_ctheta;
        stheta = next_stheta;
      }
      return false;
    }

    RowMatrixXd plan(VectorXdRef_const start, VectorXdRef_const goal, double time_limit, int max_iter, std::vector<RowMatrixXd> &turns_chain) {
      RowMatrixXd pose_seq = spatial_planner.plan(start.head<3>(), goal.head<3>(), time_limit, max_iter);
      std::vector<std::vector<RowVector2d>> subpaths;
      std::vector<RowVector2d> subpath;
      subpath.push_back(RowVector2d::Zero()); // Initial sequence of L circles
      subpath.push_back(RowVector2d::Zero()); // Initial sequence of R circles
      // TODO: take out below 2 lines
      RowMatrixXd pose_and_time_seq(pose_seq.rows(), pose_seq.cols() + 1);
      pose_and_time_seq.leftCols(3) = pose_seq;
      double t = 0.;
      pose_and_time_seq(0, 3) = t;
      for (int pose_idx = 1; pose_idx < pose_seq.rows(); ++pose_idx) {
        double x_0 = pose_seq(pose_idx - 1, 0);
        double y_0 = pose_seq(pose_idx - 1, 1);
        double theta_0 = pose_seq(pose_idx - 1, 2);

        double x_f = pose_seq(pose_idx, 0);
        double y_f = pose_seq(pose_idx, 1);
        double theta_f = pose_seq(pose_idx, 2);

        RowMatrixXd turns = turns_for_dubins_path(x_0, y_0, theta_0, x_f, y_f, theta_f, rho);
        for (int turn_idx = 0; turn_idx < turns.rows(); ++turn_idx) {
          // Since we use an RRT variant, we may have multiple of the same kind of segment in a row (e.g. L, R, or S)
          if (subpath.size() && subpath.back()(0) == turns(turn_idx, 0)) {
            subpath.back()(1) += turns(turn_idx, 1);
          } else {
            subpath.push_back(turns.row(turn_idx));
          }
        }

        // TODO: take out below 3 lines
        turns_chain.push_back(turns);
        double next_t = t + turns.col(1).sum()/vmax;
        pose_and_time_seq(pose_idx, 3) = next_t;
        t = next_t;
      }
      // TODO: take out
      return pose_and_time_seq;

      subpaths.push_back(subpath);

      double x_0 = start(0);
      double y_0 = start(1);
      double theta_0 = start(2);

      double x_f = goal(0);
      double y_f = goal(1);
      double theta_f = goal(2);

      double x_f_tmp = x_f;
      double y_f_tmp = y_f;
      double theta_f_tmp = theta_f;

      double c_f_tmp_prev = cos(theta_f_tmp);
      double s_f_tmp_prev = sin(theta_f_tmp);

      // Normalize
      for (int subpath_idx = subpaths.size() - 1; subpath_idx >= 0; --subpath_idx) {
        // We do -4 here because we're skipping over the last 3 segments, and we have to do -1 to account for zero-based indexing
        int final_C_idx = subpaths[subpath_idx].size() - 4;
        if (subpaths[subpath_idx][final_C_idx](0) == 0) {
          ++final_C_idx;
        }
        // turn_idx is the turn we're sliding
        for (int turn_idx = final_C_idx; turn_idx > 0; --turn_idx) {
          const RowVector2d &prev_turn = subpaths[subpath_idx][turn_idx - 1];
          const RowVector2d &cur_turn = subpaths[subpath_idx][turn_idx];
          const RowVector2d &next_turn = subpaths[subpath_idx][turn_idx + 1];
          if (cur_turn(0) == 0) {
            throw std::runtime_error("Encountered S segment");
          }

          // Merge turns of the same type
          if (cur_turn(0) == prev_turn(0)) {
            double combined_dist = prev_turn(1) + cur_turn(1);
            subpaths[subpath_idx][turn_idx - 1](1) = combined_dist;
            subpaths[subpath_idx].erase(subpaths[subpath_idx].begin() + turn_idx);

            double twopirho = 2*M_PI*rho;
            if (subpaths[subpath_idx][turn_idx - 1](1) >= twopirho) {
              double combined_dist_mod_twopirho = fmod(combined_dist, twopirho);
              if (cur_turn(0) == 0) {
                subpaths[subpath_idx][0](1) = combined_dist - combined_dist_mod_twopirho;
              } else {
                // cur_turn(0) == 1
                subpaths[subpath_idx][1](1) = combined_dist - combined_dist_mod_twopirho;
              }
              subpaths[subpath_idx][turn_idx - 1](1) = combined_dist_mod_twopirho;
            }
            continue;
          }

          // Move backwards along the next turn
          if (next_turn(0) == 0) {
            // S segment
            x_f_tmp -= next_turn(1)*cos(theta_f_tmp);
            y_f_tmp -= next_turn(1)*sin(theta_f_tmp);
            // theta_f_tmp stays the same
          } else {
            // C segment
            double next_C_sign = next_turn(0);
            double next_C_dist = next_turn(1);
            theta_f_tmp -= next_C_sign*next_C_dist/rho;
            double c_f_tmp = cos(theta_f_tmp);
            double s_f_tmp = sin(theta_f_tmp);
            x_f_tmp -= rho/next_C_sign*(-s_f_tmp + s_f_tmp_prev);
            y_f_tmp -= rho/next_C_sign*(c_f_tmp - c_f_tmp_prev);
            c_f_tmp_prev = c_f_tmp;
            s_f_tmp_prev = s_f_tmp;
          }

          double x_0_tmp = x_f_tmp;
          double y_0_tmp = y_f_tmp;
          double theta_0_tmp = theta_f_tmp;

          double c_0_tmp_prev = c_f_tmp_prev;
          double s_0_tmp_prev = s_f_tmp_prev;

          // Move backwards along the current turn
          if (cur_turn(0) == 0) {
            // S segment
            x_0_tmp -= cur_turn(1)*cos(theta_0_tmp);
            y_0_tmp -= cur_turn(1)*sin(theta_0_tmp);
            // theta_0_tmp stays the same
          } else {
            // C segment
            double cur_C_sign = cur_turn(0);
            double cur_C_dist = cur_turn(1);
            theta_0_tmp -= cur_C_sign*cur_C_dist/rho;
            double c_0_tmp = cos(theta_0_tmp);
            double s_0_tmp = sin(theta_0_tmp);
            x_0_tmp -= rho/cur_C_sign*(-s_0_tmp + s_0_tmp_prev);
            y_0_tmp -= rho/cur_C_sign*(c_0_tmp - c_0_tmp_prev);
            c_0_tmp_prev = c_0_tmp;
            s_0_tmp_prev = s_0_tmp;
          }

          // Move backwards along the previous turn (we're moving back along two turns but decrementing turn_idx by 1, but
          // this is ok because we're gonna delete the previous turn)
          if (prev_turn(0) == 0) {
            // S segment
            x_0_tmp -= prev_turn(1)*cos(theta_0_tmp);
            y_0_tmp -= prev_turn(1)*sin(theta_0_tmp);
            // theta_0_tmp stays the same
          } else {
            // C segment
            double prev_C_sign = prev_turn(0);
            double prev_C_dist = prev_turn(1);
            theta_0_tmp -= prev_C_sign*prev_C_dist/rho;
            double c_0_tmp = cos(theta_0_tmp);
            double s_0_tmp = sin(theta_0_tmp);
            x_0_tmp -= rho/prev_C_sign*(-s_0_tmp + s_0_tmp_prev);
            y_0_tmp -= rho/prev_C_sign*(c_0_tmp - c_0_tmp_prev);
            c_0_tmp_prev = c_0_tmp;
            s_0_tmp_prev = s_0_tmp;
          }

          double step_size = 0.1;
          double slide_amount = 0.;
          RowMatrixXd turns_before_slide(3, 2);
          turns_before_slide.row(0) = prev_turn;
          turns_before_slide.row(1) = cur_turn;
          turns_before_slide.row(2) = next_turn;

          RowMatrixXd turns_after_prev_slide(4, 2);
          turns_after_prev_slide.topRows(3) = turns_before_slide;
          turns_after_prev_slide.row(3).setZero();

          bool collision = false;

          while (slide_amount < prev_turn(1)) {
            slide_amount += step_size;
            if (slide_amount > prev_turn(1)) {
              slide_amount = prev_turn(1);
            }

            RowMatrixXd turns_after_slide = ic_sliding(turns_before_slide, x_0_tmp, y_0_tmp, theta_0_tmp, x_f_tmp, y_f_tmp, theta_f_tmp, rho, slide_amount, false);

            double x_collision;
            double y_collision;
            double theta_collision;
            int local_collision_turn_idx;
            double collision_dist;
            collision = check_path_collides(x_collision, y_collision, theta_collision,
                                            local_collision_turn_idx, collision_dist,
                                            x_0_tmp, y_0_tmp, theta_0_tmp,
                                            x_f_tmp, y_f_tmp, theta_f_tmp,
                                            turns_after_slide);

            if (local_collision_turn_idx == 0) {
              throw std::runtime_error("local_collision_turn_idx should not equal 0, since we are simply truncating local turn 0 (i.e. turn_idx - 1)");
            }
            // local_collision_turn_idx indexes into turns_after_prev_slide
            if (collision) {
              // Check if turns_after_slide has the same first column as turns_after_prev_slide. Otherwise, probably need to reduce step size
              bool reduce_step_size = false;
              for (int local_turn_idx = 0; local_turn_idx < 4; ++local_turn_idx) {
                if (turns_after_slide(local_turn_idx, 0) != turns_after_prev_slide(local_turn_idx, 0)) {
                  reduce_step_size = true;
                }
              }
              if (reduce_step_size) {
                step_size /= 2;
                slide_amount -= step_size;
                continue;
              }

              // Split subpath into two subpaths

              // New subpath before collision
              subpaths.insert(subpaths.begin() + subpath_idx, std::vector<RowVector2d>(turn_idx + local_collision_turn_idx));
              // Really, the above is turn_idx - 1 + local_collision_turn_idx + 1. We want everything up to and EXCLUDING turn_idx - 1 (the unmodified turns), then
              // everything from turn_idx - 1 up to and INCLUDING turn_idx - 1 + local_collision_turn_idx. Since the last one is INCLUDING, we need to add 1 to the size
              for (int new_subpath_turn_idx = 0; new_subpath_turn_idx < turn_idx - 1; ++new_subpath_turn_idx) {
                subpaths[subpath_idx][new_subpath_turn_idx] = subpaths[subpath_idx + 1][new_subpath_turn_idx];
              }
              for (int new_subpath_turn_idx = turn_idx - 1; new_subpath_turn_idx < turn_idx - 1 + local_collision_turn_idx; ++new_subpath_turn_idx) {
                int local_turn_idx = new_subpath_turn_idx - (turn_idx - 1); // Indexes into turns_after_prev_slide
                subpaths[subpath_idx][new_subpath_turn_idx] = turns_after_prev_slide.row(local_turn_idx);
              }
              subpaths[subpath_idx][turn_idx - 1 + local_collision_turn_idx](0) = turns_after_slide(local_collision_turn_idx, 0);
              subpaths[subpath_idx][turn_idx - 1 + local_collision_turn_idx](1) = collision_dist;

              // New subpath after collision
              // Erase everything up to the last turn we're modifying (i.e. up to and INCLUDING turn_idx + 1)
              subpaths[subpath_idx + 1].erase(subpaths[subpath_idx + 1].begin(), subpaths[subpath_idx + 1].begin() + turn_idx + 1);
              // Insert the turns in turns_after_slide that occur after the collision point
              for (int local_turn_idx = local_collision_turn_idx; local_turn_idx < num_turns_after_slide; ++local_turn_idx) {
                subpaths[subpath_idx + 1].insert(subpaths[subpath_idx + 1].begin() + local_turn_idx - local_collision_turn_idx, turns_after_prev_slide.row(local_turn_idx));
              }
              subpaths[subpath_idx + 1][0](1) -= collision_dist; // Colliding turn

              // Since local_collision_turn_idx != 0, the second new subpath consists of the latter three or fewer segments of prev_turns_after_slide,
              // so it's a shortest Dubins path and thereby a quintet path. Thus we move on to the first new subpath
              turn_idx = 1; // So we exit the loop over turn_idx

              x_f_tmp = x_collision;
              y_f_tmp = y_collision;
              theta_f_tmp = theta_collision;
              break;
            }

            turns_after_prev_slide = turns_after_slide;
          }

          if (!collision) {
            if (slide_amount != prev_turn(1) || turns_after_prev_slide(0, 1) != 0) {
              throw std::runtime_error("slide_amount != prev_turn(1) || turns_after_prev_slide(0, 1) != 0");
            }
            // Replace cur_turn and next_turn with three turns of the shortest Dubins path, and delete prev_turn
            subpaths[subpath_idx][turn_idx] = turns_after_prev_slide.row(1);
            subpaths[subpath_idx][turn_idx + 1] = turns_after_prev_slide.row(2);
            subpaths[subpath_idx].insert(subpaths[subpath_idx].begin() + turn_idx + 2, turns_after_prev_slide.row(3));
            subpaths[subpath_idx].erase(subpaths[subpath_idx].begin() + turn_idx - 1);

            // Now we've gotta update the f_tmp variables
            // The next turn we're gonna slide is turn_idx - 1, which corresponds to turns_after_prev_slide.row(1).
            // So we need to move back along turns_after_prev_slide.row(3). And then at the next iteration
            // we'll move back along turns_after_prev_slide.row(2)
            if (turns_after_prev_slide(3, 0) == 0) {
              // S segment
              x_f_tmp -= turns_after_prev_slide(3, 1)*cos(theta_f_tmp);
              y_f_tmp -= turns_after_prev_slide(3, 1)*sin(theta_f_tmp);
              // theta_f_tmp stays the same
            } else {
              // C segment
              double next_C_sign = turns_after_prev_slide(3, 0);
              double next_C_dist = turns_after_prev_slide(3, 1);
              theta_f_tmp -= next_C_sign*next_C_dist/rho;
              double c_f_tmp = cos(theta_f_tmp);
              double s_f_tmp = sin(theta_f_tmp);
              x_f_tmp -= rho/next_C_sign*(-s_f_tmp + s_f_tmp_prev);
              y_f_tmp -= rho/next_C_sign*(c_f_tmp - c_f_tmp_prev);
              c_f_tmp_prev = c_f_tmp;
              s_f_tmp_prev = s_f_tmp;
            }
          }
        }
      }
    }

    bool is_state_valid(VectorXdRef_const state_vec) {
      return spatial_planner.is_state_valid(state_vec);
    }
  private:
    NoTimeDubinsPlanner spatial_planner;
    double vmax;
    double rho;
};
