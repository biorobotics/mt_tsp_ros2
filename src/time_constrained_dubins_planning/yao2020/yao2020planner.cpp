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
  /*
  RowMatrixXd add_term = RowMatrixXd::Zero(rects.rows(), rects.cols());
  add_term(52, 3) += 1;
  spatial_planner = NoTimeDubinsPlanner(rho, occupancy, rects + add_term, map_lb, map_ub);
  */
  spatial_planner = NoTimeDubinsPlanner(rho, occupancy, rects, map_lb, map_ub);

  if (std::isinf(pose_seq(0, 0))) {
    return std::numeric_limits<double>::infinity()*RowMatrixXd::Ones(1, 2);
  }
  std::vector<std::vector<RowVector2d>> subpaths;
  std::vector<RowVector2d> subpath;
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
      // Get rid of full circles
      if (turn_dir != 0 && turn_dist >= twopirho) {
        double turn_dist_mod_twopirho = fmod(turn_dist, twopirho);
        subpath.back()(1) = turn_dist_mod_twopirho;
      }
    }
  }
  
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
    // For the initial path from the RRT, we have a sequence of Dubins paths, so the third-to-last segment must be C and that's the final
    // C eligible for sliding. For a path resulting from some IC-sliding we already did, we append dummies such that the below is true
    int final_C_idx = subpaths[subpath_idx].size() - 3;
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

    // turn_idx is the turn we're sliding
    // loop guard is turn_idx > 0 because we don't want to slide turn_idx 0 itself
    for (int turn_idx = final_C_idx; turn_idx > 0; --turn_idx) {
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

        // Get rid of full circles
        if (subpaths[subpath_idx][turn_idx - 1](0) != 0 && subpaths[subpath_idx][turn_idx - 1](1) >= twopirho) {
          double combined_dist_mod_twopirho = fmod(combined_dist, twopirho);
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

      double slide_amount = 0.;
      RowMatrixXd turns_before_slide(num_turns_after_slide, 2);
      for (int local_turn_idx = 0; local_turn_idx < num_turns_after_slide; ++local_turn_idx) {
        turns_before_slide.row(local_turn_idx) = subpaths[subpath_idx][turn_idx - 1 + local_turn_idx];
      }

      RowMatrixXd turns_after_prev_slide = turns_before_slide;

      bool collision = false;

      ICSlidingClass sliding_obj(turns_before_slide, x_0_tmp, y_0_tmp, theta_0_tmp, x_f_subpath, y_f_subpath, theta_f_subpath, rho, false, false);

      while (slide_amount < prev_turn(1)) {
        slide_amount += slide_step_size;
        if (slide_amount > prev_turn(1)) {
          slide_amount = prev_turn(1);
        }

        RowMatrixXd turns_after_slide = sliding_obj.slide(slide_amount, turns_after_prev_slide.col(1).sum());

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

        if (ignore_obstacles_in_normalization) {
          collision = false;
        }
        // local_collision_turn_idx indexes into turns_after_prev_slide
        if (collision) {
          ++num_collisions;
          if (local_collision_turn_idx == 0) {
            throw std::runtime_error("local_collision_turn_idx should not equal 0, since we are simply truncating local turn 0 (i.e. turn_idx - 1)");
          }

          // Convert collision_dist along specific turn in turns_after_slide to collision_dist along turns_after_slide overall
          collision_dist += turns_after_slide.col(1).head(local_collision_turn_idx).sum();
          //std::cout << local_collision_turn_idx << std::endl;

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
          if (!found) {
            std::cout << "collision dist = " << collision_dist << " dist of turns after prev slide " << dist << " second minus first: " << dist - collision_dist << " dist of turns after slide: " << turns_after_slide.col(1).sum() << std::endl;
            throw std::runtime_error("Could not compute local_collision_turn_idx for turns_after_prev_slide");
          }
          // Convert collision_dist along turns_after_prev_slide overall to collision_dist along specific turn in turns_after_prev_slide

          collision_dist -= turns_after_prev_slide.col(1).head(local_collision_turn_idx).sum();

          if (local_collision_turn_idx == 0) {
            throw std::runtime_error("local_collision_turn_idx should not equal 0, since we are simply truncating local turn 0 (i.e. turn_idx - 1)");
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

          if (subpaths[subpath_idx + 1].size() == 1) {
            // One row, and it must be type C. Add dummy S and C afterward
            subpaths[subpath_idx + 1].push_back(RowVector2d(0, 0));
            subpaths[subpath_idx + 1].push_back(RowVector2d(1, 0));
          }

          if (subpaths[subpath_idx + 1].size() == 2) {
            // Two rows. They must be type SC or CC because we're at the tail end of a Dubins path.
            // Either way, insert dummy C at the beginning
            if (subpaths[subpath_idx + 1][0](0) == 1) {
              // If CC, make sure we don't have two consecutive Ls or Rs
              subpaths[subpath_idx + 1].insert(subpaths[subpath_idx + 1].begin(), RowVector2d(-1, 0));
            } else {
              subpaths[subpath_idx + 1].insert(subpaths[subpath_idx + 1].begin(), RowVector2d(1, 0));
            }
          }

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

  // Now we've normalized the path. Compute length intervals of each subpath
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
  std::vector<Vector3d> elongation_intervals_per_subpath;
  std::vector<Vector3d> start_qs_per_subpath;
  std::vector<Vector3d> end_qs_per_subpath;
  double cur_length = 0.;
  for (int subpath_idx = 0; subpath_idx < subpaths.size(); ++subpath_idx) {
    if (subpaths[subpath_idx].size() != 3) {
      throw std::runtime_error("All subpaths should have size 3 after normalization");
    }
    if (subpaths[subpath_idx][0](0) != 0 &&
        subpaths[subpath_idx][1](0) != 0 &&
        subpaths[subpath_idx][2](0) != 0) {
      throw std::runtime_error("CCC path after normalization");
    }

    double x_subpath_start = x;
    double y_subpath_start = y;
    double theta_subpath_start = theta;

    start_qs_per_subpath.push_back(Vector3d(x, y, theta));
    for (int turn_idx = 0; turn_idx < subpaths[subpath_idx].size(); ++turn_idx) {
      double turn_dir = subpaths[subpath_idx][turn_idx](0);
      double turn_dist = subpaths[subpath_idx][turn_idx](1);
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

      cur_length += turn_dist;
    }

    double x_subpath_end = x;
    double y_subpath_end = y;
    double theta_subpath_end = theta;
    end_qs_per_subpath.push_back(Vector3d(x, y, theta));

    Vector3d elongation_intervals = get_elongation_intervals(x_subpath_start, y_subpath_start, theta_subpath_start, x_subpath_end, y_subpath_end, theta_subpath_end, rho);
    elongation_intervals_per_subpath.push_back(elongation_intervals);
  }

  double des_length = vmax*(goal(3) - start(3));
  std::vector<bool> tried_adding_circles_to_subpath(subpaths.size(), false);

  std::vector<std::pair<double, double>> feasible_length_intervals;

  // int max_elongation_iter = 2000;
  int max_elongation_iter = 1; // TODO: replace with the above
  for (int elongation_iter = 0; elongation_iter < max_elongation_iter; ++elongation_iter) {
    // Incrementally compute Minkowski sum of the feasible length intervals of the individual subpaths
    feasible_length_intervals.clear();
    feasible_length_intervals.push_back(std::pair<double, double>(0., 0.));

    for (int subpath_idx = 0; subpath_idx < subpaths.size(); ++subpath_idx) {
      const Ref<const Vector3d> elongation_intervals = elongation_intervals_per_subpath[subpath_idx];
      std::vector<std::pair<double, double>> new_feasible_length_intervals;
      for (auto interval : feasible_length_intervals) {
        if (std::isinf(elongation_intervals(1))) {
          // Only one interval
          new_feasible_length_intervals.push_back(std::pair<double, double>(interval.first + elongation_intervals(0), std::numeric_limits<double>::infinity()));
        } else {
          // Two intervals
          new_feasible_length_intervals.push_back(std::pair<double, double>(interval.first + elongation_intervals(0), interval.first + elongation_intervals(1)));
          new_feasible_length_intervals.push_back(std::pair<double, double>(interval.first + elongation_intervals(2), std::numeric_limits<double>::infinity()));
        }
      }
      feasible_length_intervals = new_feasible_length_intervals;
    }

    bool feasible_ignoring_obstacles = false;
    for (auto interval : feasible_length_intervals) {
      if (interval.first <= des_length && des_length <= interval.second) {
        feasible_ignoring_obstacles = true;
        break;
      }
    }

    if (debug) {
      std::cout << "Is elongation of normalized subpaths to desired sum of lengths feasible, keeping subpath endpoints fixed and ignoring collisions during elongation? ";
      if (feasible_ignoring_obstacles) {
        std::cout << "yes" << std::endl;
      } else {
        std::cout << "no" << std::endl;
      }
    }

    if (!feasible_ignoring_obstacles) {
      return std::numeric_limits<double>::infinity()*RowMatrixXd::Ones(1, 2);
    }

    double remaining_length = des_length - cur_length;

    if (remaining_length >= 2*M_PI*rho) {
      for (int subpath_idx = 0; subpath_idx < subpaths.size(); ++subpath_idx) {
        if (!tried_adding_circles_to_subpath[subpath_idx]) {
          // Try adding full circles on the left
          x = start_qs_per_subpath[subpath_idx](0);
          y = start_qs_per_subpath[subpath_idx](1);
          theta = start_qs_per_subpath[subpath_idx](2);
          if (spatial_planner.collision_free(x, y, theta, 1, 2*M_PI*rho, x, y)) {
            subpaths[subpath_idx].insert(subpaths[subpath_idx].begin(), RowVector2d::Zero());
            subpaths[subpath_idx][0](0) = 1;
            double remaining_length_mod_twopirho = fmod(remaining_length, twopirho);
            subpaths[subpath_idx][0](1) = remaining_length - remaining_length_mod_twopirho;
            remaining_length = remaining_length_mod_twopirho;
            break;
          }

          // Try adding full circles on the right
          if (spatial_planner.collision_free(x, y, theta, -1, 2*M_PI*rho, x, y)) {
            subpaths[subpath_idx].insert(subpaths[subpath_idx].begin(), RowVector2d::Zero());
            subpaths[subpath_idx][0](0) = -1;
            double remaining_length_mod_twopirho = fmod(remaining_length, twopirho);
            subpaths[subpath_idx][0](1) = remaining_length - remaining_length_mod_twopirho;
            remaining_length = remaining_length_mod_twopirho;
            break;
          }
          tried_adding_circles_to_subpath[subpath_idx] = true;
        }
      }

      if (debug) {
        int num_total_turns = 0;
        for (auto subpath : subpaths) {
          num_total_turns += subpath.size();
        }
        RowMatrixXd overall_turns(num_total_turns, 2);
        int overall_turn_idx = 0;
        for (int subpath_idx2 = 0; subpath_idx2 < subpaths.size(); ++subpath_idx2) {
          for (int turn_idx2 = 0; turn_idx2 < subpaths[subpath_idx2].size(); ++turn_idx2) {
            overall_turns.row(overall_turn_idx) = subpaths[subpath_idx2][turn_idx2];
            ++overall_turn_idx;
          }
        }
        path_seq.push_back(overall_turns);
      }
    }

    for (int subpath_idx = 0; subpath_idx < subpaths.size(); ++subpath_idx) {
      // Elongate this subpath as close as we can get to the desired length via IC-sliding

      // Slide initial incomplete arc forward around left turning circle 
      RowMatrixXd turns_before_slide(num_turns_after_slide, 2);
      turns_before_slide(0, 0) = 1;
      turns_before_slide(0, 1) = 0;
      for (int local_turn_idx = 0; local_turn_idx < num_turns_after_slide; ++local_turn_idx) {
        if (subpaths[subpath_idx][0](1) >= 2*M_PI*rho) {
          if (subpaths[subpath_idx].size() != 4) {
            throw std::runtime_error("If first segment length >= 2*M_PI*rho, subpaths[subpath_idx].size() should equal 4");
          }
          // We have a full circle at the beginning
          turns_before_slide.row(1 + local_turn_idx) = subpaths[subpath_idx][1 + local_turn_idx];
        } else {
          // We don't have a full circle at the beginning
          if (subpaths[subpath_idx].size() != 3) {
            throw std::runtime_error("If first segment length < 2*M_PI*rho, subpaths[subpath_idx].size() should equal 3");
          }
          turns_before_slide.row(1 + local_turn_idx) = subpaths[subpath_idx][local_turn_idx];
        }
      }

      RowMatrixXd turns_after_prev_slide = turns_before_slide;

      double x_0_subpath = start_qs_per_subpath[subpath_idx](0);
      double y_0_subpath = start_qs_per_subpath[subpath_idx](1);
      double theta_0_subpath = start_qs_per_subpath[subpath_idx](2);

      double x_f_subpath = end_qs_per_subpath[subpath_idx](0);
      double y_f_subpath = end_qs_per_subpath[subpath_idx](1);
      double theta_f_subpath = end_qs_per_subpath[subpath_idx](2);

      ICSlidingClass sliding_obj(turns_before_slide, x_0_subpath, y_0_subpath, theta_0_subpath, x_f_subpath, y_f_subpath, theta_f_subpath, rho, false, true);
      double slide_amount = 0.;
      while (slide_amount < remaining_length) {
        slide_amount += slide_step_size;
        if (slide_amount > remaining_length) {
          slide_amount = remaining_length;
        }
        RowMatrixXd turns_after_slide = sliding_obj.slide(slide_amount, turns_after_prev_slide.col(1).sum());

        double x_collision;
        double y_collision;
        double theta_collision;
        int collision_turn_idx;
        double collision_dist;
        bool collision = check_path_collides(x_collision, y_collision, theta_collision,
                                             collision_turn_idx, collision_dist,
                                             x_0_subpath, y_0_subpath, theta_0_subpath,
                                             x_f_subpath, y_f_subpath, theta_f_subpath,
                                             turns_after_slide);

        if (collision) {
          // Split

          // Convert collision_dist along specific turn in turns_after_slide to collision_dist along turns_after_slide overall
          collision_dist += turns_after_slide.col(1).head(collision_turn_idx).sum();

          // It is possible that turns_after_prev_slide and turns_after_slide have a different sequence of segment types.
          // Compute the collision_turn_idx for turns_after_prev_slide using collision_dist
          double dist = 0.;
          bool found = false;
          for (collision_turn_idx = 0; collision_turn_idx < num_turns_after_slide; ++collision_turn_idx) {
            dist += turns_after_prev_slide(collision_turn_idx, 1);
            if (dist >= collision_dist) {
              found = true;
              break;
            }
          }
          if (!found) {
            std::cout << "collision dist = " << collision_dist << " dist of turns after prev slide " << dist << " second minus first: " << dist - collision_dist << " dist of turns after slide: " << turns_after_slide.col(1).sum() << std::endl;
            throw std::runtime_error("Could not compute collision_turn_idx for turns_after_prev_slide");
          }

          // Convert collision_dist along turns_after_prev_slide overall to collision_dist along specific turn in turns_after_prev_slide
          collision_dist -= turns_after_prev_slide.col(1).head(collision_turn_idx).sum();

          // Split subpath into two subpaths

          // New subpath before collision
          subpaths.insert(subpaths.begin() + subpath_idx, std::vector<RowVector2d>(collision_turn_idx + 1));
          for (int turn_idx = 0; turn_idx < collision_turn_idx; ++turn_idx) {
            subpaths[subpath_idx][turn_idx] = turns_after_prev_slide.row(turn_idx);
          }
          subpaths[subpath_idx].back()(1) = collision_dist;

          if (subpaths[subpath_idx + 1].size() == 4) {
            // We had a full circle at the beginning
            subpaths[subpath_idx].insert(subpaths[subpath_idx].begin(), subpaths[subpath_idx + 1][0]);
          }

          // New subpath after collision
          subpaths[subpath_idx + 1].clear();
          
          // Insert the turns in turns_after_prev_slide that occur after the collision point
          for (int turn_idx = collision_turn_idx; turn_idx < num_turns_after_slide; ++turn_idx) {
            subpaths[subpath_idx + 1].push_back(turns_after_prev_slide.row(turn_idx));
          }
          subpaths[subpath_idx + 1][0](1) -= collision_dist; // Colliding turn

          if (subpaths[subpath_idx + 1].size() == 1) {
            // One row, and it must be type C. Add dummy S and C afterward
            subpaths[subpath_idx + 1].push_back(RowVector2d(0, 0));
            subpaths[subpath_idx + 1].push_back(RowVector2d(1, 0));
          }

          if (subpaths[subpath_idx + 1].size() == 2) {
            // Two rows. They must be type SC or CC because we're at the tail end of a Dubins path.
            // Either way, insert dummy C at the beginning
            if (subpaths[subpath_idx + 1][0](0) == 1) {
              // If CC, make sure we don't have two consecutive Ls or Rs
              subpaths[subpath_idx + 1].insert(subpaths[subpath_idx + 1].begin(), RowVector2d(-1, 0));
            } else {
              subpaths[subpath_idx + 1].insert(subpaths[subpath_idx + 1].begin(), RowVector2d(1, 0));
            }
          }

          // Don't set x_f_subpath = x_collision and so forth because we want to find the 
          // point at the same distance along the path in the sliding iteration before collision
          x_f_subpath = x_0_subpath;
          y_f_subpath = y_0_subpath;
          theta_f_subpath = theta_0_subpath;
          double ctheta = cos(theta_f_subpath);
          double stheta = sin(theta_f_subpath);
          for (int turn_idx = 0; turn_idx < collision_turn_idx + 1; ++turn_idx) {
            double turn_dir = subpaths[subpath_idx][turn_idx](0);
            double turn_dist = subpaths[subpath_idx][turn_idx](1);
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

          start_qs_per_subpath.insert(start_qs_per_subpath.begin() + subpath_idx, start_qs_per_subpath[subpath_idx]);
          end_qs_per_subpath.insert(end_qs_per_subpath.begin() + subpath_idx, Vector3d(x_f_subpath, y_f_subpath, theta_f_subpath));
          start_qs_per_subpath[subpath_idx + 1] = end_qs_per_subpath[subpath_idx];

          
          Vector3d elongation_intervals = get_elongation_intervals(x_0_subpath, y_0_subpath, theta_0_subpath, x_f_subpath, y_f_subpath, theta_f_subpath, rho);
          elongation_intervals_per_subpath.insert(elongation_intervals_per_subpath.begin() + subpath_idx, elongation_intervals);

          tried_adding_circles_to_subpath.insert(tried_adding_circles_to_subpath.begin() + subpath_idx, false);

          // Don't work on the new segments until we've tried elongating all subsequent ones
          ++subpath_idx;
          break;
        }

        turns_after_prev_slide = turns_after_slide;
      }

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

      // Slide initial incomplete arc forward around right turning circle 

      // Slide final incomplete arc backward around left turning circle 

      // Slide final incomplete arc backward around right turning circle 

      // TODO: if we split, update the start and end points, the elongation intervals, 
      // and the tried_adding_circles per subpath
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
