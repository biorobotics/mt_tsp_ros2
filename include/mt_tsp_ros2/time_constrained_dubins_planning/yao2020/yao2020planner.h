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
      double twopirho = 2*M_PI*rho;

      RowMatrixXd pose_seq = spatial_planner.plan(start.head<3>(), goal.head<3>(), time_limit, max_iter);
      std::vector<std::vector<RowVector2d>> subpaths;
      std::vector<RowVector2d> subpath;
      subpath.push_back(RowVector2d::Zero()); // Initial sequence of L circles
      subpath.push_back(RowVector2d::Zero()); // Initial sequence of R circles
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

          double turn_dir = subpath.back()(0);
          double turn_dist = subpath.back()(1);
          if (turn_dir != 0 && turn_dist >= twopirho) {
            double turn_dist_mod_twopirho = fmod(turn_dist, twopirho);
            if (turn_dir == 1) {
              // L
              subpath[0](1) += turn_dist - turn_dist_mod_twopirho;
            } else {
              // turn_dir == -1
              // R
              subpath[1](1) += turn_dist - turn_dist_mod_twopirho;
            }
            subpath.back()(1) = turn_dist_mod_twopirho;
          }
        }
      }
      
      // We start by sliding the final C segment in the path.
      // When we do so, we will replace the final C segment and the subsequent
      // two segments by the shortest Dubins path, then delete the segment
      // before the final C segment.
      // If the final segment is type S, we add an additional dummy segment in front,
      // so the S and the dummy will get replaced. If the final segment is type C,
      // we will be sliding it, so we add two dummy segments in front to replace.
      if (subpath.back()(0) == 0) {
        // Final segment is type S
        subpath.push_back(RowVector2d::Zero());
      } else {
        // Final segment is type C
        subpath.push_back(RowVector2d::Zero());
        subpath.push_back(RowVector2d::Zero());
      }

      subpaths.push_back(subpath);

      double x_0 = start(0);
      double y_0 = start(1);
      double theta_0 = start(2);

      double x_f = goal(0);
      double y_f = goal(1);
      double theta_f = goal(2);

      // This will always be the pose after traversing backward from the final pose in the subpath
      // along the final four segments
      double x_0_tmp = x_f;
      double y_0_tmp = y_f;
      double theta_0_tmp = theta_f;
      double c_0_tmp_prev = cos(theta_0_tmp);
      double s_0_tmp_prev = sin(theta_0_tmp);

      double x_f_subpath = x_f;
      double y_f_subpath = y_f;
      double theta_f_subpath = theta_f;

      // Normalize
      for (int subpath_idx = subpaths.size() - 1; subpath_idx >= 0; --subpath_idx) {
        std::cout << "subpath " << subpath_idx << std::endl;
        int final_C_idx = subpaths[subpath_idx].size() - 3; // We appended dummies such that this is true
        // turn_idx is the turn we're sliding
        // loop guard is turn_idx > 2 because turn_idx 0 and 1 are the initial L and R circles and we don't want to slide those,
        // and we don't want to slide along those circles either so we don't want to slide turn_idx 2

        for (int turn_idx = subpaths[subpath_idx].size() - 1; turn_idx >= final_C_idx; --turn_idx) {
          const RowVector2d &cur_turn = subpaths[subpath_idx][turn_idx];

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
        }

        for (int turn_idx = final_C_idx; turn_idx > 2; --turn_idx) {
          std::cout << "turn " << turn_idx << std::endl;
          std::cout << "printing current subpath" << std::endl;
          for (int turn_idx2 = 0; turn_idx2 < subpaths[subpath_idx].size(); ++turn_idx2) {
            std::cout << "turn " << turn_idx2 << ": " << subpaths[subpath_idx][turn_idx2](0) << " " << subpaths[subpath_idx][turn_idx2](1) << std::endl;
          }
          std::cout << std::endl;

          const RowVector2d &prev_turn = subpaths[subpath_idx][turn_idx - 1];
          const RowVector2d &cur_turn = subpaths[subpath_idx][turn_idx];
          if (cur_turn(0) == 0) {
            throw std::runtime_error("Cannot IC-slide S segment");
          }

          // Move backwards along the previous turn
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

          // Merge turns of the same type
          if (cur_turn(0) == prev_turn(0)) {
            double combined_dist = prev_turn(1) + cur_turn(1);
            subpaths[subpath_idx][turn_idx - 1](1) = combined_dist;
            subpaths[subpath_idx].erase(subpaths[subpath_idx].begin() + turn_idx);

            if (subpaths[subpath_idx][turn_idx - 1](0) != 0 && subpaths[subpath_idx][turn_idx - 1](1) >= twopirho) {
              double combined_dist_mod_twopirho = fmod(combined_dist, twopirho);
              if (cur_turn(0) == 1) {
                // L
                subpaths[subpath_idx][0](1) += combined_dist - combined_dist_mod_twopirho;
              } else {
                // cur_turn(0) == -1
                // R
                subpaths[subpath_idx][1](1) += combined_dist - combined_dist_mod_twopirho;
              }
              subpaths[subpath_idx][turn_idx - 1](1) = combined_dist_mod_twopirho;
            }
            continue;
          }

          // TODO: take out
          bool do_check = true;
          if (do_check) {
            double x = x_0;
            double y = y_0;
            double theta = theta_0;
            double ctheta = cos(theta);
            double stheta = sin(theta);

            double next_x;
            double next_y;
            double next_theta;
            double next_ctheta;
            double next_stheta;
            for (int subpath_idx2 = 0; subpath_idx2 < subpaths.size(); ++subpath_idx2) {
              for (int turn_idx2 = 0; turn_idx2 < subpaths[subpath_idx2].size(); ++turn_idx2) {
                double turn_dir = subpaths[subpath_idx2][turn_idx2](0);
                double turn_dist = subpaths[subpath_idx2][turn_idx2](1);
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

                x = next_x;
                y = next_y;
                theta = next_theta;
                ctheta = next_ctheta;
                stheta = next_stheta;
              }
            }
            if (std::abs(x - x_f) > 1e-10 ||
                std::abs(y - y_f) > 1e-10 || 
                std::abs(angdiff(theta_f, theta)) > 1e-10) {
              std::cout << x << " " << y << " " << theta << std::endl;
              std::cout << x_f << " " << y_f << " " << theta_f << std::endl;
              throw std::runtime_error("Before sliding, integrating sequence of segments doesnt reach goal state");
            }
          }

          double step_size = 0.1;
          double slide_amount = 0.;
          RowMatrixXd turns_before_slide(num_turns_after_slide, 2);
          for (int local_turn_idx = 0; local_turn_idx < num_turns_after_slide; ++local_turn_idx) {
            turns_before_slide.row(local_turn_idx) = subpaths[subpath_idx][turn_idx - 1 + local_turn_idx];
          }

          RowMatrixXd turns_after_prev_slide = turns_before_slide;

          bool collision = false;

          while (slide_amount < prev_turn(1)) {
            slide_amount += step_size;
            if (slide_amount > prev_turn(1)) {
              slide_amount = prev_turn(1);
            }

            RowMatrixXd turns_after_slide = ic_sliding(turns_before_slide, x_0_tmp, y_0_tmp, theta_0_tmp, x_f_subpath, y_f_subpath, theta_f_subpath, rho, slide_amount, false);

            double x_collision;
            double y_collision;
            double theta_collision;
            int local_collision_turn_idx;
            double collision_dist;
            collision = check_path_collides(x_collision, y_collision, theta_collision,
                                            local_collision_turn_idx, collision_dist,
                                            x_0_tmp, y_0_tmp, theta_0_tmp,
                                            x_f_subpath, y_f_subpath, theta_f_subpath,
                                            turns_after_slide);

            // local_collision_turn_idx indexes into turns_after_prev_slide
            if (collision) {
              if (local_collision_turn_idx == 0) {
                throw std::runtime_error("local_collision_turn_idx should not equal 0, since we are simply truncating local turn 0 (i.e. turn_idx - 1)");
              }

              // Check if turns_after_slide has the same first column as turns_after_prev_slide. Otherwise, reduce step size
              bool reduce_step_size = false;
              for (int local_turn_idx = 0; local_turn_idx < num_turns_after_slide; ++local_turn_idx) {
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
              // everything from turn_idx - 1 up to and INCLUDING turn_idx - 1 + local_collision_turn_idx. Since the last one is INCLUDING, we need to add 1
              for (int new_subpath_turn_idx = 0; new_subpath_turn_idx < turn_idx - 1; ++new_subpath_turn_idx) {
                subpaths[subpath_idx][new_subpath_turn_idx] = subpaths[subpath_idx + 1][new_subpath_turn_idx]; // The subpath we were just working on is now subpath_idx + 1
              }
              for (int new_subpath_turn_idx = turn_idx - 1; new_subpath_turn_idx < turn_idx - 1 + local_collision_turn_idx; ++new_subpath_turn_idx) {
                int local_turn_idx = new_subpath_turn_idx - (turn_idx - 1); // Indexes into turns_after_prev_slide
                subpaths[subpath_idx][new_subpath_turn_idx] = turns_after_prev_slide.row(local_turn_idx);
              }
              subpaths[subpath_idx][turn_idx - 1 + local_collision_turn_idx](0) = turns_after_slide(local_collision_turn_idx, 0);
              subpaths[subpath_idx][turn_idx - 1 + local_collision_turn_idx](1) = collision_dist;

              // If final segment of new subpath is type S, add 1 dummy segment. Otherwise, 2
              if (subpaths[subpath_idx].back()(0) == 0) {
                subpaths[subpath_idx].push_back(RowVector2d::Zero());
              } else {
                subpaths[subpath_idx].push_back(RowVector2d::Zero());
                subpaths[subpath_idx].push_back(RowVector2d::Zero());
              }

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
              turn_idx = 0; // So we exit the loop over turn_idx

              x_f_subpath = x_collision;
              y_f_subpath = y_collision;
              theta_f_subpath = theta_collision;

              x_0_tmp = x_f_subpath;
              y_0_tmp = y_f_subpath;
              theta_0_tmp = theta_f_subpath;
              c_0_tmp_prev = cos(theta_0_tmp);
              s_0_tmp_prev = sin(theta_0_tmp);
              break;
            }

            turns_after_prev_slide = turns_after_slide;
          }

          if (!collision) {
            if (slide_amount != prev_turn(1) || turns_after_prev_slide(0, 1) != 0) {
              throw std::runtime_error("slide_amount != prev_turn(1) || turns_after_prev_slide(0, 1) != 0");
            }
            // Replace cur_turn and next_turn with three turns of the shortest Dubins path, and delete prev_turn
            subpaths[subpath_idx].erase(subpaths[subpath_idx].begin() + turn_idx - 1); // We reduced this segment to zero length
            std::cout << "editing" << std::endl;
            std::cout << subpaths[subpath_idx][turn_idx - 1] << std::endl;
            std::cout << subpaths[subpath_idx][turn_idx] << std::endl;
            std::cout << subpaths[subpath_idx][turn_idx + 1] << std::endl;
            subpaths[subpath_idx][turn_idx - 1] = turns_after_prev_slide.row(1);
            subpaths[subpath_idx][turn_idx] = turns_after_prev_slide.row(2);
            subpaths[subpath_idx][turn_idx + 1] = turns_after_prev_slide.row(3);
            std::cout << "done editing" << std::endl;
            std::cout << subpaths[subpath_idx][turn_idx - 1] << std::endl;
            std::cout << subpaths[subpath_idx][turn_idx] << std::endl;
            std::cout << subpaths[subpath_idx][turn_idx + 1] << std::endl;
            std::cout << std::endl;

            // TODO: take out
            bool do_check = true;
            if (do_check) {
              double x = x_0;
              double y = y_0;
              double theta = theta_0;
              double ctheta = cos(theta);
              double stheta = sin(theta);

              double next_x;
              double next_y;
              double next_theta;
              double next_ctheta;
              double next_stheta;
              for (int subpath_idx2 = 0; subpath_idx2 < subpaths.size(); ++subpath_idx2) {
                for (int turn_idx2 = 0; turn_idx2 < subpaths[subpath_idx2].size(); ++turn_idx2) {
                  double turn_dir = subpaths[subpath_idx2][turn_idx2](0);
                  double turn_dist = subpaths[subpath_idx2][turn_idx2](1);
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

                  x = next_x;
                  y = next_y;
                  theta = next_theta;
                  ctheta = next_ctheta;
                  stheta = next_stheta;
                }
              }
              if (std::abs(x - x_f) > 1e-10 ||
                  std::abs(y - y_f) > 1e-10 || 
                  std::abs(angdiff(theta_f, theta)) > 1e-10) {
                std::cout << x << " " << y << " " << theta << std::endl;
                std::cout << x_f << " " << y_f << " " << theta_f << std::endl;
                throw std::runtime_error("After sliding, integrating sequence of segments doesnt reach goal state");
              }
            }
          }
        }
      }

      double x = start(0);
      double y = start(1);
      double theta = start(2);
      double ctheta = cos(theta);
      double stheta = sin(theta);
      double t = 0.;

      double next_x;
      double next_y;
      double next_theta;
      double next_ctheta;
      double next_stheta;

      if (turns_chain.size()) {
        throw std::runtime_error("This function should not be called with a nonempty turns_chain");
      }

      RowMatrixXd pose_and_time_seq(subpaths.size() + 1, 4);
      int seq_idx = 0;
      for (auto subpath : subpaths) {
        turns_chain.push_back(RowMatrixXd(subpath.size(), 2));

        pose_and_time_seq(seq_idx, 0) = x;
        pose_and_time_seq(seq_idx, 1) = y;
        pose_and_time_seq(seq_idx, 2) = theta;
        pose_and_time_seq(seq_idx, 3) = t;
        ++seq_idx;

        for (int turn_idx = 0; turn_idx < subpath.size(); ++turn_idx) {
          turns_chain.back().row(turn_idx) = subpath[turn_idx];

          double turn_dir = subpath[turn_idx](0);
          double turn_dist = subpath[turn_idx](1);
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

          x = next_x;
          y = next_y;
          theta = next_theta;
          ctheta = next_ctheta;
          stheta = next_stheta;
        }
        t += turns_chain.back().col(1).sum()/vmax;
      }

      pose_and_time_seq(seq_idx, 0) = x;
      pose_and_time_seq(seq_idx, 1) = y;
      pose_and_time_seq(seq_idx, 2) = theta;
      pose_and_time_seq(seq_idx, 3) = t;

      return pose_and_time_seq;
    }

    bool is_state_valid(VectorXdRef_const state_vec) {
      return spatial_planner.is_state_valid(state_vec);
    }
  private:
    NoTimeDubinsPlanner spatial_planner;
    double vmax;
    double rho;
};
