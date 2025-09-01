#include "mt_tsp_ros2/time_constrained_dubins_planning/yao2020/yao2020planner.h"

Yao2020Planner::Yao2020Planner(double vmax, double rho, RowMatrixXbRef_const occupancy, RowMatrixXdRef_const rects, Vector2dRef_const map_lb, Vector2dRef_const map_ub) : spatial_planner(rho, occupancy, rects, map_lb, map_ub), vmax(vmax), rho(rho),
                                                                                                                                                          occupancy(occupancy), rects(rects), map_lb(map_lb), map_ub(map_ub) {
}

bool Yao2020Planner::check_path_collides(double &x_collision, double &y_collision, double &theta_collision,
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

// path_seq is the sequence of paths (turn sequences) we got via sliding
RowMatrixXd Yao2020Planner::plan(VectorXdRef_const start, VectorXdRef_const goal, double time_limit, int max_iter, std::vector<RowMatrixXd> &path_seq, bool just_return_rrt_trj, bool ignore_obstacles_in_normalization) {
  double twopirho = 2*M_PI*rho;

  RowMatrixXd pose_seq = spatial_planner.plan(start.head<3>(), goal.head<3>(), time_limit, max_iter);
  RowMatrixXd add_term = RowMatrixXd::Zero(rects.rows(), rects.cols());
  add_term(52, 3) += 1;
  spatial_planner = NoTimeDubinsPlanner(rho, occupancy, rects + add_term, map_lb, map_ub);

  if (std::isinf(pose_seq(0, 0))) {
    return std::numeric_limits<double>::infinity()*RowMatrixXd::Ones(1, 2);
  }
  std::vector<std::vector<RowVector2d>> subpaths;
  std::vector<RowVector2d> subpath;
  subpath.push_back(RowVector2d::Zero()); // Initial sequence of L circles
  subpath[0](0) = 1;
  subpath.push_back(RowVector2d::Zero()); // Initial sequence of R circles
  subpath[1](0) = -1;
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
  
  // The below code is fine, but it doesn't take advantage of the fact that the spatial planner
  // returns a sequence of Dubins paths, so we don't need to append dummy segments
  /*
  // We start by sliding the final C segment in the path.
  // When we do so, we will replace the final C segment and the subsequent
  // two segments by the Dubins path, then delete the segment
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
  */

  subpaths.push_back(subpath);

  if (debug) {
    int num_total_turns = 0;
    for (auto subpath : subpaths) {
      num_total_turns += subpath.size();
    }
    RowMatrixXd overall_turns(num_total_turns, 2);
    int overall_turn_idx = 0;
    for (auto subpath : subpaths) {
      for (int turn_idx = 0; turn_idx < subpath.size(); ++turn_idx) {
        overall_turns.row(overall_turn_idx) = subpath[turn_idx];
        ++overall_turn_idx;
      }
    }
    path_seq.push_back(overall_turns);
  }

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

  int num_collisions = 0;

  // Normalize
  for (int subpath_idx = subpaths.size() - 1; subpath_idx >= 0; --subpath_idx) {
    if (just_return_rrt_trj) {
      break;
    }
    if (subpath_idx != 0) {
      throw std::runtime_error("We should never have subpath_idx != 0");
    }
    // TODO: take out
    /*
    if (num_collisions == 1) {
      break;
    }
    */
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
      if (debug) {
        std::cout << "turn " << turn_idx << std::endl;
        std::cout << "We currently have " << subpaths.size() << " subpaths" << std::endl;
        std::cout << "printing current subpath" << std::endl;
        for (int turn_idx2 = 0; turn_idx2 < subpaths[subpath_idx].size(); ++turn_idx2) {
          std::cout << "turn " << turn_idx2 << ": " << subpaths[subpath_idx][turn_idx2](0) << " " << subpaths[subpath_idx][turn_idx2](1) << std::endl;
        }
        std::cout << std::endl;
      }

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

      if (debug) {
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

      ICSlidingClass sliding_obj(turns_before_slide, x_0_tmp, y_0_tmp, theta_0_tmp, x_f_subpath, y_f_subpath, theta_f_subpath, rho, false);

      while (slide_amount < prev_turn(1)) {
        slide_amount += step_size;
        if (slide_amount > prev_turn(1)) {
          slide_amount = prev_turn(1);
        }
        if (debug) {
          // std::cout << "sliding by " << slide_amount << std::endl;
        }

        RowMatrixXd turns_after_slide = sliding_obj.slide(slide_amount);

        /*
        std::ofstream turns_before_slide_file("/home/noopygbhat/catkin_ws/src/mapf/scripts/turns_before_slide.txt");
        turns_before_slide_file << turns_before_slide.format(IOFormat(FullPrecision)) << std::endl;
        turns_before_slide_file.close();

        std::ofstream turns_after_prev_slide_file("/home/noopygbhat/catkin_ws/src/mapf/scripts/turns_after_prev_slide.txt");
        turns_after_prev_slide_file << turns_after_prev_slide.format(IOFormat(FullPrecision)) << std::endl;
        turns_after_prev_slide_file.close();

        std::ofstream turns_after_slide_file("/home/noopygbhat/catkin_ws/src/mapf/scripts/turns_after_slide.txt");
        turns_after_slide_file << turns_after_slide.format(IOFormat(FullPrecision)) << std::endl;
        turns_after_slide_file.close();

        std::ofstream boundary_conditions_file("/home/noopygbhat/catkin_ws/src/mapf/scripts/boundary_conditions.txt");
        boundary_conditions_file << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << x_0_tmp << " " << y_0_tmp << " " << theta_0_tmp << std::endl;
        boundary_conditions_file << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << x_f_subpath << " " << y_f_subpath << " " << theta_f_subpath << std::endl;
        boundary_conditions_file.close();
        */

        /*
        if (debug) {
          std::cout << "turns_after_slide" << std::endl;
          std::cout << turns_after_slide << std::endl;
        }
        */

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

        /*
        if (debug) {
          std::cout << "did collision check" << std::endl;
        }
        */

        if (ignore_obstacles_in_normalization) {
          collision = false;
        }
        // local_collision_turn_idx indexes into turns_after_prev_slide
        if (collision) {
          ++num_collisions;
          /*
          std::cout << "turns_after_slide for collision" << std::endl;
          std::cout << turns_after_slide << std::endl;
          std::cout << local_collision_turn_idx << std::endl;
          */
          /*
          if (debug) {
            std::cout << "collision" << std::endl;

            std::cout << "local_collision_turn_idx" << std::endl;
            std::cout << local_collision_turn_idx << std::endl;
          }
          */

          if (local_collision_turn_idx == 0) {
            throw std::runtime_error("local_collision_turn_idx should not equal 0, since we are simply truncating local turn 0 (i.e. turn_idx - 1)");
          }

          // Convert collision_dist along specific turn in turns_after_slide to collision_dist along turns_after_slide overall
          collision_dist += turns_after_slide.col(1).head(local_collision_turn_idx).sum();
          //std::cout << local_collision_turn_idx << std::endl;

          /*
          if (debug) {
            std::cout << "updated collision dist" << std::endl;
          }
          */
          
          // It is possible that turns_after_prev_slide and turns_after_slide have a different sequence of segment types.
          // Compute the local_collision_turn_idx for turns_after_prev_slide using collision_dist
          double dist = 0.;
          bool found = false;
          for (local_collision_turn_idx = 0; local_collision_turn_idx < num_turns_after_slide; ++local_collision_turn_idx) {
            dist += turns_after_prev_slide(local_collision_turn_idx, 1);
            if (dist >= collision_dist) {
              found = true;
              break;
            }
          }
          /*
          std::cout << dist << " " << collision_dist << std::endl;
          std::cout << turns_after_prev_slide << std::endl << std::endl;
          std::cout << turns_after_slide << std::endl << std::endl;
          */
          if (!found) {
            std::cout << "collision dist = " << collision_dist << " dist of turns after prev slide " << dist << " second minus first: " << dist - collision_dist << " dist of turns after slide: " << turns_after_slide.col(1).sum() << std::endl;
            throw std::runtime_error("Could not compute local_collision_turn_idx for turns_after_prev_slide");
          }
          // Convert collision_dist along turns_after_prev_slide overall to collision_dist along specific turn in turns_after_prev_slide

          /*
          if (debug) {
            std::cout << "updated local_collision_turn_idx" << std::endl;
          }
          */

          collision_dist -= turns_after_prev_slide.col(1).head(local_collision_turn_idx).sum();

          /*
          if (debug) {
            std::cout << "updated collision_dist again" << std::endl;
          }
          */

          if (local_collision_turn_idx == 0) {
            throw std::runtime_error("local_collision_turn_idx should not equal 0, since we are simply truncating local turn 0 (i.e. turn_idx - 1)");
          }

          /*
          if (debug) {
            std::cout << "local_collision_turn_idx: " << local_collision_turn_idx << std::endl;
          }
          */

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
          subpaths[subpath_idx][turn_idx - 1 + local_collision_turn_idx](0) = turns_after_prev_slide(local_collision_turn_idx, 0);
          subpaths[subpath_idx][turn_idx - 1 + local_collision_turn_idx](1) = collision_dist;

          // If local_collision_turn_idx == 3, then the last three segments of subpaths[subpath_idx] comprise a Dubins path.
          // The third-to-last segment must therefore by type C, and we'll be sliding it next.

          /*
          std::cout << "num turns in next subpath " << turn_idx - 1 + local_collision_turn_idx << std::endl;
          for (auto turn : subpaths[subpath_idx]) {
            std::cout << turn << std::endl;
          }
          std::cout << "turns_after_prev_slide" << std::endl;
          std::cout << turns_after_prev_slide << std::endl;
          std::cout << "turns_after_slide" << std::endl;
          std::cout << turns_after_slide << std::endl;
          */

          if (local_collision_turn_idx == 2) {
            // If local_collision_turn_idx == 2, then the last two segments of subpaths[subpath_idx] comprise a Dubins path. 
            // The second-to-last segment must therefore be type C, and we'll be sliding it next. We need two segments in front, but we only have
            // one, so append one dummy segment
            if (subpaths[subpath_idx][subpaths[subpath_idx].size() - 2](0) == 0) {
              throw std::runtime_error("Second-to-last segment is not type C");
            }
            subpaths[subpath_idx].push_back(RowVector2d::Zero());
          } else if (local_collision_turn_idx == 1) {
            // If local_collision_turn_idx == 1, we don't want to slide segment 1 next, because that will result in a discontinuous change in the path between
            // segment 2 and segment 1 unless we were sliding segment 1 exactly such that the corresponding endpoints of the two segments were staying pointed at
            // one another. If segment 0 is type S, then we can slide the segment before segment 0 because it must be type C. If segment 0 is type C, then we slide segment 0.
            // We need two segments in front, but currently we only have one, so add a dummy
            if (subpaths[subpath_idx][subpaths[subpath_idx].size() - 1](0) == 0) {
              throw std::runtime_error("Last segment is not type C");
            }
            // Check if C segment
            if (turns_after_prev_slide(0, 0) != 0) {
              subpaths[subpath_idx].push_back(RowVector2d::Zero());
            }
          }

          // New subpath after collision
          subpaths[subpath_idx + 1].clear();
          
          // Insert the turns in turns_after_prev_slide that occur after the collision point
          for (int local_turn_idx = local_collision_turn_idx; local_turn_idx < num_turns_after_slide; ++local_turn_idx) {
            subpaths[subpath_idx + 1].push_back(turns_after_prev_slide.row(local_turn_idx));
          }
          subpaths[subpath_idx + 1][0](1) -= collision_dist; // Colliding turn

          // Don't set x_f_subpath = x_collision and so forth because we want to find the 
          // point at the same distance along the path in the sliding iteration before collision
          x_f_subpath = x_0_tmp;
          y_f_subpath = y_0_tmp;
          theta_f_subpath = theta_0_tmp;
          double ctheta = cos(theta_f_subpath);
          double stheta = sin(theta_f_subpath);
          // Technically the subpath may have size larger than turn_idx + local_collision_turn_idx now
          // but only because we pushed back dummy segments, so we don't need to integrate those
          for (int turn_idx2 = turn_idx - 1; turn_idx2 < turn_idx + local_collision_turn_idx; ++turn_idx2) {
            double turn_dir = subpaths[subpath_idx][turn_idx2](0);
            double turn_dist = subpaths[subpath_idx][turn_idx2](1);
            if (turn_dir == 0) {
              // S segment
              x_f_subpath += turn_dist*ctheta;
              y_f_subpath += turn_dist*stheta;
            } else {
              // C segment
              double rho_times_turn_dir = rho*turn_dir;
              theta_f_subpath += turn_dist/rho_times_turn_dir;
              double next_ctheta = cos(theta_f_subpath);
              double next_stheta = sin(theta_f_subpath);
              x_f_subpath += rho_times_turn_dir*(-stheta + next_stheta);
              y_f_subpath += rho_times_turn_dir*(ctheta - next_ctheta);
              ctheta = next_ctheta;
              stheta = next_stheta;
            }
          }

          /*
          x_f_subpath = x_0;
          y_f_subpath = y_0;
          theta_f_subpath = theta_0;
          ctheta = cos(theta_f_subpath);
          stheta = sin(theta_f_subpath);
          for (int subpath_idx2 = 0; subpath_idx2 < subpath_idx + 1; ++subpath_idx2) {
            for (int turn_idx2 = 0; turn_idx2 < subpaths[subpath_idx2].size(); ++turn_idx2) {
              double turn_dir = subpaths[subpath_idx2][turn_idx2](0);
              double turn_dist = subpaths[subpath_idx2][turn_idx2](1);
              double rho_times_turn_dir = rho*turn_dir;
              if (turn_dir == 0) {
                // S segment
                x_f_subpath += turn_dist*ctheta;
                y_f_subpath += turn_dist*stheta;
              } else {
                // C segment
                theta_f_subpath += turn_dist/rho_times_turn_dir;
                double next_ctheta = cos(theta_f_subpath);
                double next_stheta = sin(theta_f_subpath);
                x_f_subpath += rho_times_turn_dir*(-stheta + next_stheta);
                y_f_subpath += rho_times_turn_dir*(ctheta - next_ctheta);
                ctheta = next_ctheta;
                stheta = next_stheta;
              }
            }
          }
          */

          x_0_tmp = x_f_subpath;
          y_0_tmp = y_f_subpath;
          theta_0_tmp = theta_f_subpath;
          c_0_tmp_prev = ctheta;
          s_0_tmp_prev = stheta;

          // Since local_collision_turn_idx != 0, the second new subpath consists of the latter three or fewer segments of turns_after_prev_slide,
          // so it's a Dubins path and thereby a quintet path. Thus we move on to the first new subpath
          ++subpath_idx;
          turn_idx = 0; // So we exit the loop over turn_idx

          if (debug) {
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
              throw std::runtime_error("After splitting, integrating sequence of segments doesnt reach goal state");
            }
          }

          break;
        }

        turns_after_prev_slide = turns_after_slide;

        if (debug) {
          /*
          if (debug) {
            std::cout << "getting subpath" << std::endl;
          }
          */
          int num_total_turns = 0;
          for (auto subpath : subpaths) {
            num_total_turns += subpath.size();
          }
          RowMatrixXd overall_turns(num_total_turns, 2);
          int overall_turn_idx = 0;
          for (int subpath_idx2 = 0; subpath_idx2 < subpaths.size(); ++subpath_idx2) {
            for (int turn_idx2 = 0; turn_idx2 < subpaths[subpath_idx2].size(); ++turn_idx2) {
              if (subpath_idx2 == subpath_idx && turn_idx - 1 <= turn_idx2 && turn_idx2 <= turn_idx + 2) {
                overall_turns.row(overall_turn_idx) = turns_after_prev_slide.row(turn_idx2 - (turn_idx - 1));
              } else {
                overall_turns.row(overall_turn_idx) = subpaths[subpath_idx2][turn_idx2];
              }
              ++overall_turn_idx;
            }
          }
          path_seq.push_back(overall_turns);
        }
      }

      if (!collision) {
        if (slide_amount != prev_turn(1) || turns_after_prev_slide(0, 1) != 0) {
          throw std::runtime_error("slide_amount != prev_turn(1) || turns_after_prev_slide(0, 1) != 0");
        }
        /*
        if (debug) {
          std::cout << "updating subpath" << std::endl;
        }
        */
        // Replace cur_turn and next_turn with three turns of the Dubins path, and delete prev_turn
        subpaths[subpath_idx].erase(subpaths[subpath_idx].begin() + turn_idx - 1); // We reduced this segment to zero length
        subpaths[subpath_idx][turn_idx - 1] = turns_after_prev_slide.row(1);
        subpaths[subpath_idx][turn_idx] = turns_after_prev_slide.row(2);
        subpaths[subpath_idx][turn_idx + 1] = turns_after_prev_slide.row(3);

        if (debug) {
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

  int num_total_turns = 0;
  for (auto subpath : subpaths) {
    num_total_turns += subpath.size();
  }
  RowMatrixXd overall_turns(num_total_turns, 2);
  int overall_turn_idx = 0;
  for (auto subpath : subpaths) {
    for (int turn_idx = 0; turn_idx < subpath.size(); ++turn_idx) {
      overall_turns.row(overall_turn_idx) = subpath[turn_idx];
      ++overall_turn_idx;
    }
  }
  return overall_turns;
}

bool Yao2020Planner::is_state_valid(VectorXdRef_const state_vec) {
  return spatial_planner.is_state_valid(state_vec);
}
