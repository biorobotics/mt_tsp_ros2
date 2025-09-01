#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_motion_validator_rects.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/geometry_utils.h"
#include <ompl/base/spaces/SE2StateSpace.h>
#include <ompl/base/spaces/TimeStateSpace.h>

DubinsMotionValidatorRects::DubinsMotionValidatorRects(const ob::SpaceInformationPtr si, double vmax, double rho, RowMatrixXdRef_const &rects, Vector2dRef_const map_lb, Vector2dRef_const map_ub, const Ref<const VectorXd> &start, const Ref<const VectorXd> &goal, bool monte_carlo_prop) : DubinsMotionValidator(si, vmax, rho), rects(rects), map_lb(map_lb), map_ub(map_ub), start(start), goal(goal), monte_carlo_prop(monte_carlo_prop) {
  for (int row = 0; row < rects.rows(); ++row) {
    Point point1(rects(row, 0), rects(row, 1)); // xlow, ylow
    Point point2(rects(row, 0), rects(row, 3)); // xlow, yhigh
    Point point3(rects(row, 2), rects(row, 3)); // xhigh, yhigh
    Point point4(rects(row, 2), rects(row, 1)); // xhigh, ylow
    segments.push_back(Segment(point1, point2));
    segments.push_back(Segment(point2, point3));
    segments.push_back(Segment(point3, point4));
    segments.push_back(Segment(point4, point1));
  }
  // Add the map boundaries
  Point point1(map_lb(0), map_lb(1)); // xlow, ylow
  Point point2(map_lb(0), map_ub(1)); // xlow, yhigh
  Point point3(map_ub(0), map_ub(1)); // xhigh, yhigh
  Point point4(map_ub(0), map_lb(1)); // xhigh, ylow
  segments.push_back(Segment(point1, point2));
  segments.push_back(Segment(point2, point3));
  segments.push_back(Segment(point3, point4));
  segments.push_back(Segment(point4, point1));

  aabb_tree = Tree(segments.begin(), segments.end()); 

  rng_ = ompl::RNG(1);
}

const Tree &DubinsMotionValidatorRects::get_aabb_tree() {
  return aabb_tree;
}

bool DubinsMotionValidatorRects::collision_free(double x, double y, double theta, double turn_dir, double turn_dist, double next_x, double next_y) const {
  /*
  // Line segment
  for (int rect_idx = 0; rect_idx < rects.rows(); ++rect_idx) {
    if (turn_dir == 0 && line_segment_intersects_rect(x, y, next_x, next_y, 
                                                      rects(rect_idx, 0), rects(rect_idx, 1),
                                                      rects(rect_idx, 2), rects(rect_idx, 3))) {
      return false;
    } else if (turn_dir != 0) {
      if (arc_intersects_rect(x, y, theta, turn_dir, turn_dist/rho, rho, 
                              next_x, next_y,
                              rects(rect_idx, 0), rects(rect_idx, 1),
                              rects(rect_idx, 2), rects(rect_idx, 3))) {
        return false;
      }
    }
  }

  // Check if we exit the map
  if (!point_in_rect(x, y, 
                     map_lb(0), map_lb(1),
                     map_ub(0), map_ub(1)) ||
      !point_in_rect(next_x, next_y, 
                     map_lb(0), map_lb(1),
                     map_ub(0), map_ub(1))) {
    return false;
  }

  if (turn_dir != 0) {
    double xir = map_lb(0);
    double xfr = map_ub(0);
    double yir = map_lb(1);
    double yfr = map_ub(1);
    if (arc_intersects_line_segment(x, y, theta, turn_dir, turn_dist/rho, rho,
                                    xir, yir, xir, yfr) ||
        arc_intersects_line_segment(x, y, theta, turn_dir, turn_dist/rho, rho,
                                    xir, yfr, xfr, yfr) ||
        arc_intersects_line_segment(x, y, theta, turn_dir, turn_dist/rho, rho,
                                    xfr, yfr, xfr, yir) ||
        arc_intersects_line_segment(x, y, theta, turn_dir, turn_dist/rho, rho,
                                    xfr, yir, xir, yir)) {
      return false;
    }
  }

  return true;
  */

  // S segment
  if (turn_dir == 0) {
    Point point1(x, y);
    Point point2(next_x, next_y);
    Segment segment(point1, point2);
    return !aabb_tree.do_intersect(segment);
  }

  // C segment
  double c = cos(theta);
  double s = sin(theta);

  Vector2d dir(c, s);
  Vector2d perp(-s*turn_dir, c*turn_dir);
  Vector2d center = Vector2d(x, y) + perp*rho;

  CGAL::Bbox_2 bbox(center(0) - rho, center(1) - rho, center(0) + rho, center(1) + rho);
  std::list<Primitive_id> primitives;
  aabb_tree.all_intersected_primitives(bbox, std::back_inserter(primitives));

  for (Primitive_id primitive : primitives) {
    if (arc_intersects_line_segment(x, y, theta, turn_dir, turn_dist/rho, rho,
                                    primitive->source().x(), primitive->source().y(), primitive->target().x(), primitive->target().y())) {
      return false;
    }
  }
  return true;

  // Below code doesn't work. I think it's becauseCGAL always assumes the arc goes counterclockwise from point1 to point3
  /*
  if (turn_dist >= 2*M_PI*rho) {
    Vector2d dir(c, s);
    Vector2d perp(-s*turn_dir, c*turn_dir);
    Vector2d center = Vector2d(x, y) + perp*rho;

    Point_Circular_k center_CGAL(center(0), center(1));
    Circle_2 circle_CGAL(center_CGAL, rho*rho);

    std::list<Primitive_id> primitives;
    aabb_tree.all_intersected_primitives(circle_CGAL.bbox(), std::back_inserter(primitives));
    for (Primitive_id primitive : primitives) {
      Matrix2d intersections = line_segment_intersects_circle(primitive->source().x(), primitive->source().y(), primitive->target().x(), primitive->target().y(), center(0), center(1), rho);
      if (std::isfinite(intersections(0, 0))) {
        return false;
      }
    }
    return true;
  } else {
    double rho_times_turn_dir = rho*turn_dir;
    double theta_mid = theta + 0.5*turn_dist/rho_times_turn_dir;
    double x_mid = x + rho_times_turn_dir*(-s + sin(theta_mid));
    double y_mid = y + rho_times_turn_dir*(c - cos(theta_mid));

    Point_Circular_k point1_CGAL(x, y);
    Point_Circular_k point2_CGAL(x_mid, y_mid);
    Point_Circular_k point3_CGAL(next_x, next_y);

    Circular_arc_2 arc_CGAL(point1_CGAL, point2_CGAL, point3_CGAL);

    std::list<Primitive_id> primitives;
    aabb_tree.all_intersected_primitives(arc_CGAL.bbox(), std::back_inserter(primitives));
    for (Primitive_id primitive : primitives) {
      if (arc_intersects_line_segment(x, y, theta, turn_dir, turn_dist/rho, rho,
                                      primitive->source().x(), primitive->source().y(), primitive->target().x(), primitive->target().y())) {
        return false;
      }
    }
    return true;
  }
  */
}

// If we set validMotion = false, stopState isn't used so we don't have to populate it. Same deal with turns
RowMatrixXd DubinsMotionValidatorRects::checkMotionForward(const ob::State *s1, const ob::State *s2, double maxDuration, ob::State *stopState, bool &reach, bool &validMotion, int num_checks) const {
  return checkMotionForward_internal(s1, s2, maxDuration, stopState, reach, validMotion);
  /*
  RowMatrixXd turns1 = checkMotionForward_internal(s1, s2, maxDuration, stopState, reach, validMotion);
  bool reach2;
  bool valid2;
  RowMatrixXd turns2 = DubinsMotionValidator::checkMotionForward(s1, s2, maxDuration, stopState, reach2, valid2);
  if (validMotion != valid2) {
    std::cout << "mismatch, using finer collision check resolution" << std::endl;
    DubinsMotionValidator::checkMotionForward(s1, s2, maxDuration, stopState, reach2, valid2, 100000);
    if (validMotion != valid2) {
      throw std::runtime_error("Mismatch on valid forward");
    }
  }
  return turns1;
  */
}

RowMatrixXd DubinsMotionValidatorRects::checkMotionForward_internal(const ob::State *s1, const ob::State *s2, double maxDuration, ob::State *stopState, bool &reach, bool &validMotion) const {
  double x1 = s1->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getX();
  double y1 = s1->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getY();
  double theta1 = s1->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getYaw();
  double t1 = s1->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position;

  double x2 = s2->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getX();
  double y2 = s2->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getY();
  double theta2 = s2->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getYaw();
  double t2 = s2->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position;

  auto timer_start = std::chrono::high_resolution_clock::now();
  auto timer_stop = std::chrono::high_resolution_clock::now();
  auto nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();

  bool compute_elongated_path = !monte_carlo_prop; // false; // t2 - t1 <= maxDuration;
  RowMatrixXd turns(1, 2);
  if (compute_elongated_path) {
    timer_start = std::chrono::high_resolution_clock::now();
    turns = elongated_dubins_path(x1, y1, theta1, x2, y2, theta2, vmax*(t2 - t1), rho, false);
    timer_stop = std::chrono::high_resolution_clock::now();
    nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
    path_elongation_time += ((double)nanos)/1e9;
    if (std::isinf(turns(0, 0))) {
      reach = false;
      validMotion = false;
      // std::cout << "failed on initial forward elongation check" << std::endl;
      ++num_discarded_samples;
      return turns;
    }
  } else {
    double duration = rng_.uniformReal(0, maxDuration);
    turns(0, 1) = duration*vmax;
    double dir_tmp = rng_.uniformReal(0, 3);
    if (dir_tmp < 1) {
      turns(0, 0) = -1;
    } else if (dir_tmp < 2) {
      turns(0, 0) = 0;
    } else {
      turns(0, 0) = 1;
    }
  }

  double ctheta1 = cos(theta1);
  double stheta1 = sin(theta1);

  double x = x1;
  double y = y1;
  double theta = theta1;
  double t = t1;

  double ctheta = ctheta1;
  double stheta = stheta1;

  double next_x = x;
  double next_y = y;
  double next_theta = theta;
  double next_t = t;

  double next_ctheta = ctheta;
  double next_stheta = stheta;

  double valid_path_length = 0.;

  /*
  // Check if the state reached after maxDuration can get to the goal (assuming no obstacles)
  for (int turn_idx = 0; turn_idx < turns.rows(); ++turn_idx) {
    double turn_dir = turns(turn_idx, 0);
    double turn_dist = turns(turn_idx, 1);
    if (turn_dist == 0) {
      continue;
    }
    bool stop_on_this_turn = t2 - t1 > maxDuration + 1e-2 && valid_path_length + turn_dist >= maxDuration*vmax;
    if (stop_on_this_turn) {
      turn_dist = maxDuration*vmax - valid_path_length;
    }
    double rho_times_turn_dir = rho*turn_dir;

    if (turn_dir == 0) {
      // S segment
      next_theta = theta;
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

    next_t = t + turn_dist/vmax;

    valid_path_length += turn_dist;

    if (stop_on_this_turn) {
      if (!check_elongation_possible(next_x, next_y, next_theta, goal(0), goal(1), goal(2), vmax*(goal(3) - next_t), rho)) {
        reach = false;
        validMotion = false;
        // std::cout << "failed on forward elongation check" << std::endl;
        return turns;
      }
      break;
    }

    x = next_x;
    y = next_y;
    theta = next_theta;
    t = next_t;

    ctheta = next_ctheta;
    stheta = next_stheta;
  }

  x = x1;
  y = y1;
  theta = theta1;
  t = t1;

  ctheta = ctheta1;
  stheta = stheta1;

  next_x = x;
  next_y = y;
  next_theta = theta;
  next_t = t;

  next_ctheta = ctheta;
  next_stheta = stheta;

  valid_path_length = 0.;

  */

  // Now perform collision-checks
  timer_start = std::chrono::high_resolution_clock::now();
  for (int turn_idx = 0; turn_idx < turns.rows(); ++turn_idx) {
    double turn_dir = turns(turn_idx, 0);
    double turn_dist = turns(turn_idx, 1);
    if (turn_dist == 0) {
      continue;
    }
    bool stop_on_this_turn = valid_path_length + turn_dist >= maxDuration*vmax;
    if (stop_on_this_turn) {
      turn_dist = maxDuration*vmax - valid_path_length;
    }
    double rho_times_turn_dir = rho*turn_dir;

    if (turn_dir == 0) {
      // S segment
      next_theta = theta;
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

    next_t = t + turn_dist/vmax;

    if (!collision_free(x, y, theta, turn_dir, turn_dist, next_x, next_y)) {
      reach = false;
      validMotion = false;
      turns(turn_idx, 1) = turn_dist; // Not needed, I think

      timer_stop = std::chrono::high_resolution_clock::now();
      nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
      collision_check_time += ((double)nanos)/1e9;
      // std::cout << "forward collision check failed" << std::endl;
      return turns.topRows(turn_idx + 1);
    }

    valid_path_length += turn_dist;

    if (stop_on_this_turn) {
      reach = false;
      stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setX(next_x);
      stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setY(next_y);
      stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setYaw(next_theta);
      stopState->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position = next_t;
      validMotion = true;
      turns(turn_idx, 1) = turn_dist;

      timer_stop = std::chrono::high_resolution_clock::now();
      nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
      collision_check_time += ((double)nanos)/1e9;

      return turns.topRows(turn_idx + 1);
    }

    x = next_x;
    y = next_y;
    theta = next_theta;
    t = next_t;

    ctheta = next_ctheta;
    stheta = next_stheta;
  }

  reach = compute_elongated_path;
  stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setX(next_x);
  stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setY(next_y);
  stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setYaw(next_theta);
  stopState->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position = next_t;
  validMotion = true;

  timer_stop = std::chrono::high_resolution_clock::now();
  nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
  collision_check_time += ((double)nanos)/1e9;

  return turns;
}

RowMatrixXd DubinsMotionValidatorRects::checkMotionBackward(const ob::State *s1, const ob::State *s2, double maxDuration, ob::State *stopState, bool &reach, bool &validMotion, int num_checks) const {
  return checkMotionBackward_internal(s1, s2, maxDuration, stopState, reach, validMotion);
  /*
  RowMatrixXd turns1 = checkMotionBackward_internal(s1, s2, maxDuration, stopState, reach, validMotion);
  bool reach2;
  bool valid2;
  RowMatrixXd turns2 = DubinsMotionValidator::checkMotionBackward(s1, s2, maxDuration, stopState, reach2, valid2);
  if (validMotion != valid2) {
    DubinsMotionValidator::checkMotionBackward(s1, s2, maxDuration, stopState, reach2, valid2, 100000);
    if (validMotion != valid2) {
      std::cout << validMotion << " " << valid2 << std::endl;
      throw std::runtime_error("Mismatch on valid backward");
    }
  }
  return turns1;
  */
}

RowMatrixXd DubinsMotionValidatorRects::checkMotionBackward_internal(const ob::State *s1, const ob::State *s2, double maxDuration, ob::State *stopState, bool &reach, bool &validMotion) const {
  if (!si_->isValid(s2)) {
    // TODO: do I still need this?
    reach = false;
    validMotion = false;
    return std::numeric_limits<double>::infinity()*RowMatrixXd::Ones(1, 2);
  }
  double x1 = s1->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getX();
  double y1 = s1->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getY();
  double theta1 = s1->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getYaw();
  double t1 = s1->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position;

  double x2 = s2->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getX();
  double y2 = s2->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getY();
  double theta2 = s2->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getYaw();
  double t2 = s2->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position;

  auto timer_start = std::chrono::high_resolution_clock::now();
  auto timer_stop = std::chrono::high_resolution_clock::now();
  auto nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();

  RowMatrixXd turns(1, 2);
  bool compute_elongated_path = !monte_carlo_prop; // false; // t2 - t1 <= maxDuration;
  if (compute_elongated_path) {
    timer_start = std::chrono::high_resolution_clock::now();
    turns = elongated_dubins_path(x1, y1, theta1, x2, y2, theta2, vmax*(t2 - t1), rho, false);
    timer_stop = std::chrono::high_resolution_clock::now();
    nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
    path_elongation_time += ((double)nanos)/1e9;
    if (std::isinf(turns(0, 0))) {
      reach = false;
      validMotion = false;
      // std::cout << "failed on initial backward elongation check" << std::endl;
      ++num_discarded_samples;
      return turns;
    }
  } else {
    double duration = rng_.uniformReal(0, maxDuration);
    turns(0, 1) = duration*vmax;
    double dir_tmp = rng_.uniformReal(0, 3);
    if (dir_tmp < 1) {
      turns(0, 0) = -1;
    } else if (dir_tmp < 2) {
      turns(0, 0) = 0;
    } else {
      turns(0, 0) = 1;
    }
  }

  double ctheta2 = cos(theta2);
  double stheta2 = sin(theta2);

  double x = x2;
  double y = y2;
  double theta = theta2;
  double t = t2;

  double ctheta = ctheta2;
  double stheta = stheta2;

  double next_x = x;
  double next_y = y;
  double next_theta = theta;
  double next_t = t;

  double next_ctheta = ctheta;
  double next_stheta = stheta;

  double valid_path_length = 0.;

  /*
  // Check if the state reached after maxDuration can get to the start (assuming no obstacles)
  for (int turn_idx = turns.rows() - 1; turn_idx >= 0; --turn_idx) {
    double turn_dir = turns(turn_idx, 0);
    double turn_dist = turns(turn_idx, 1);
    if (turn_dist == 0) {
      continue;
    }
    bool stop_on_this_turn = t2 - t1 > maxDuration + 1e-2 && valid_path_length + turn_dist >= maxDuration*vmax;
    if (stop_on_this_turn) {
      turn_dist = maxDuration*vmax - valid_path_length;
    }
    double rho_times_turn_dir = rho*turn_dir;

    next_t = t - turn_dist/vmax;
    if (turn_dir == 0) {
      // S segment
      next_theta = theta;
      next_x = x - turn_dist*ctheta;
      next_y = y - turn_dist*stheta;
    } else {
      // C segment
      next_theta = theta - turn_dist/rho_times_turn_dir;
      next_ctheta = cos(next_theta);
      next_stheta = sin(next_theta);
      next_x = x - rho_times_turn_dir*(stheta - next_stheta);
      next_y = y - rho_times_turn_dir*(-ctheta + next_ctheta);
    }

    valid_path_length += turn_dist;

    if (stop_on_this_turn) {
      if (!check_elongation_possible(start(0), start(1), start(2), next_x, next_y, next_theta, vmax*(next_t - start(3)), rho)) {
        reach = false;
        validMotion = false;
        // std::cout << "failed on backward elongation check" << std::endl;
        return turns;
      }
    }
    x = next_x;
    y = next_y;
    theta = next_theta;
    t = next_t;

    ctheta = next_ctheta;
    stheta = next_stheta;
  }

  x = x2;
  y = y2;
  theta = theta2;
  t = t2;

  ctheta = ctheta2;
  stheta = stheta2;

  next_x = x;
  next_y = y;
  next_theta = theta;
  next_t = t;

  next_ctheta = ctheta;
  next_stheta = stheta;

  valid_path_length = 0.;
  */

  // Now perform collision-checks
  timer_start = std::chrono::high_resolution_clock::now();
  for (int turn_idx = turns.rows() - 1; turn_idx >= 0; --turn_idx) {
    double turn_dir = turns(turn_idx, 0);
    double turn_dist = turns(turn_idx, 1);
    if (turn_dist == 0) {
      continue;
    }

    bool stop_on_this_turn = valid_path_length + turn_dist >= maxDuration*vmax;
    if (stop_on_this_turn) {
      turn_dist = maxDuration*vmax - valid_path_length;
    }
    double rho_times_turn_dir = rho*turn_dir;

    next_t = t - turn_dist/vmax;
    if (turn_dir == 0) {
      // S segment
      next_theta = theta;
      next_x = x - turn_dist*ctheta;
      next_y = y - turn_dist*stheta;
    } else {
      // C segment
      next_theta = theta - turn_dist/rho_times_turn_dir;
      next_ctheta = cos(next_theta);
      next_stheta = sin(next_theta);
      next_x = x - rho_times_turn_dir*(stheta - next_stheta);
      next_y = y - rho_times_turn_dir*(-ctheta + next_ctheta);
    }

    // Use next_x etc because it's checkMotionBackward
    if (!collision_free(next_x, next_y, next_theta, turn_dir, turn_dist, x, y)) {
      reach = false;
      validMotion = false;
      turns(turn_idx, 1) = turn_dist; // Not needed, I think

      timer_stop = std::chrono::high_resolution_clock::now();
      nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
      collision_check_time += ((double)nanos)/1e9;

      // std::cout << "backward collision check failed" << std::endl;
      return turns.bottomRows(turns.rows() - turn_idx);
    }

    valid_path_length += turn_dist;

    if (stop_on_this_turn) {
      reach = false;
      stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setX(next_x);
      stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setY(next_y);
      stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setYaw(next_theta);
      stopState->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position = next_t;
      validMotion = true;
      turns(turn_idx, 1) = turn_dist;

      timer_stop = std::chrono::high_resolution_clock::now();
      nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
      collision_check_time += ((double)nanos)/1e9;

      return turns.bottomRows(turns.rows() - turn_idx);
    }
    x = next_x;
    y = next_y;
    theta = next_theta;
    t = next_t;

    ctheta = next_ctheta;
    stheta = next_stheta;
  }

  reach = compute_elongated_path;
  stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setX(next_x);
  stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setY(next_y);
  stopState->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setYaw(next_theta);
  stopState->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position = next_t;
  validMotion = true;

  timer_stop = std::chrono::high_resolution_clock::now();
  nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
  collision_check_time += ((double)nanos)/1e9;

  return turns;
}

void DubinsMotionValidatorRects::checkMotionForward_connect(std::vector<RowMatrixXd> &turns_vec, std::vector<Vector4d> &states_vec, const ob::State *s1, const ob::State *s2, double maxDuration, bool &reach) const {
  double x1 = s1->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getX();
  double y1 = s1->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getY();
  double theta1 = s1->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getYaw();
  double t1 = s1->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position;

  double x2 = s2->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getX();
  double y2 = s2->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getY();
  double theta2 = s2->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getYaw();
  double t2 = s2->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position;

  auto timer_start = std::chrono::high_resolution_clock::now();
  RowMatrixXd turns = elongated_dubins_path(x1, y1, theta1, x2, y2, theta2, vmax*(t2 - t1), rho, false);
  auto timer_stop = std::chrono::high_resolution_clock::now();
  auto nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
  path_elongation_time += ((double)nanos)/1e9;
  if (std::isinf(turns(0, 0))) {
    reach = false;
    return;
  }

  double ctheta1 = cos(theta1);
  double stheta1 = sin(theta1);

  double x = x1;
  double y = y1;
  double theta = theta1;
  double t = t1;

  double ctheta = ctheta1;
  double stheta = stheta1;

  double next_x = x;
  double next_y = y;
  double next_theta = theta;
  double next_t = t;

  double next_ctheta = ctheta;
  double next_stheta = stheta;

  double cur_segment_length = 0.;

  timer_start = std::chrono::high_resolution_clock::now();

  int start_turn_idx = 0;
  double maxDist = maxDuration*vmax;

  for (int turn_idx = 0; turn_idx < turns.rows(); ++turn_idx) {
    double turn_dir = turns(turn_idx, 0);
    double turn_dist = turns(turn_idx, 1);
    if (turn_dist == 0) {
      continue;
    }
    double rho_times_turn_dir = rho*turn_dir;

    while (cur_segment_length + turn_dist >= maxDist) {
      // Move only until we get to maxDist along current segment
      turn_dist = maxDist - cur_segment_length;

      if (turn_dir == 0) {
        // S segment
        next_theta = theta;
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

      next_t = t + turn_dist/vmax;

      if (!collision_free(x, y, theta, turn_dir, turn_dist, next_x, next_y)) {
        // If we collide, return the states and turns we've already accumulated
        reach = false;

        timer_stop = std::chrono::high_resolution_clock::now();
        nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
        collision_check_time += ((double)nanos)/1e9;
        return;
      }

      // If we don't collide, add states and turns
      turns_vec.push_back(turns.block(start_turn_idx, 0, turn_idx - start_turn_idx + 1, 2));
      turns_vec.back()(turn_idx - start_turn_idx, 1) = turn_dist; // The last turn is too long. Shorten it
      states_vec.push_back(Vector4d(next_x, next_y, next_theta, next_t));

      x = next_x;
      y = next_y;
      theta = next_theta;
      t = next_t;

      ctheta = next_ctheta;
      stheta = next_stheta;

      cur_segment_length = 0.; // Start a new segment
      turns(turn_idx, 1) -= turn_dist;
      turn_dist = turns(turn_idx, 1);

      start_turn_idx = turn_idx;
    }

    if (turn_dir == 0) {
      // S segment
      next_theta = theta;
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

    next_t = t + turn_dist/vmax;

    if (!collision_free(x, y, theta, turn_dir, turn_dist, next_x, next_y)) {
      reach = false;

      timer_stop = std::chrono::high_resolution_clock::now();
      nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
      collision_check_time += ((double)nanos)/1e9;
      return;
    }

    x = next_x;
    y = next_y;
    theta = next_theta;
    t = next_t;

    ctheta = next_ctheta;
    stheta = next_stheta;

    cur_segment_length += turn_dist;
  }

  reach = true;

  // If we reach, need to add the final state and the final sequence of turns, add states and turns
  turns_vec.push_back(turns.block(start_turn_idx, 0, turns.rows() - start_turn_idx, 2));
  states_vec.push_back(Vector4d(next_x, next_y, next_theta, next_t));

  timer_stop = std::chrono::high_resolution_clock::now();
  nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
  collision_check_time += ((double)nanos)/1e9;
}

void DubinsMotionValidatorRects::checkMotionBackward_connect(std::vector<RowMatrixXd> &turns_vec, std::vector<Vector4d> &states_vec, const ob::State *s1, const ob::State *s2, double maxDuration, bool &reach) const {
  if (!si_->isValid(s2)) {
    // TODO: do I still need this?
    reach = false;
    return;
  }
  double x1 = s1->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getX();
  double y1 = s1->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getY();
  double theta1 = s1->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getYaw();
  double t1 = s1->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position;

  double x2 = s2->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getX();
  double y2 = s2->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getY();
  double theta2 = s2->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getYaw();
  double t2 = s2->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position;

  auto timer_start = std::chrono::high_resolution_clock::now();
  RowMatrixXd turns = elongated_dubins_path(x1, y1, theta1, x2, y2, theta2, vmax*(t2 - t1), rho, false);
  auto timer_stop = std::chrono::high_resolution_clock::now();
  auto nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
  path_elongation_time += ((double)nanos)/1e9;
  if (std::isinf(turns(0, 0))) {
    reach = false;
    return;
  }

  double ctheta2 = cos(theta2);
  double stheta2 = sin(theta2);

  double x = x2;
  double y = y2;
  double theta = theta2;
  double t = t2;

  double ctheta = ctheta2;
  double stheta = stheta2;

  double next_x = x;
  double next_y = y;
  double next_theta = theta;
  double next_t = t;

  double next_ctheta = ctheta;
  double next_stheta = stheta;

  double cur_segment_length = 0.;

  timer_start = std::chrono::high_resolution_clock::now();
  
  int start_turn_idx = turns.rows() - 1;
  double maxDist = maxDuration*vmax;

  for (int turn_idx = turns.rows() - 1; turn_idx >= 0; --turn_idx) {
    double turn_dir = turns(turn_idx, 0);
    double turn_dist = turns(turn_idx, 1);
    if (turn_dist == 0) {
      continue;
    }
    double rho_times_turn_dir = rho*turn_dir;

    while (cur_segment_length + turn_dist >= maxDist) {
      // Move only until we get to maxDist along current segment
      turn_dist = maxDist - cur_segment_length;

      next_t = t - turn_dist/vmax;
      if (turn_dir == 0) {
        // S segment
        next_theta = theta;
        next_x = x - turn_dist*ctheta;
        next_y = y - turn_dist*stheta;
      } else {
        // C segment
        next_theta = theta - turn_dist/rho_times_turn_dir;
        next_ctheta = cos(next_theta);
        next_stheta = sin(next_theta);
        next_x = x - rho_times_turn_dir*(stheta - next_stheta);
        next_y = y - rho_times_turn_dir*(-ctheta + next_ctheta);
      }

      // Use next_x etc because it's checkMotionBackward
      if (!collision_free(next_x, next_y, next_theta, turn_dir, turn_dist, x, y)) {
        // If we collide, return the states and turns we've already accumulated
        reach = false;

        timer_stop = std::chrono::high_resolution_clock::now();
        nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
        collision_check_time += ((double)nanos)/1e9;

        return;
      }

      // If we don't collide, add states and turns
      turns_vec.push_back(turns.block(turn_idx, 0, start_turn_idx - turn_idx + 1, 2));
      turns_vec.back()(0, 1) = turn_dist; // The last turn is too long. Shorten it
      states_vec.push_back(Vector4d(next_x, next_y, next_theta, next_t));

      x = next_x;
      y = next_y;
      theta = next_theta;
      t = next_t;

      ctheta = next_ctheta;
      stheta = next_stheta;

      cur_segment_length = 0.; // Start a new segment
      turns(turn_idx, 1) -= turn_dist;
      turn_dist = turns(turn_idx, 1);

      start_turn_idx = turn_idx;
    }

    next_t = t - turn_dist/vmax;
    if (turn_dir == 0) {
      // S segment
      next_theta = theta;
      next_x = x - turn_dist*ctheta;
      next_y = y - turn_dist*stheta;
    } else {
      // C segment
      next_theta = theta - turn_dist/rho_times_turn_dir;
      next_ctheta = cos(next_theta);
      next_stheta = sin(next_theta);
      next_x = x - rho_times_turn_dir*(stheta - next_stheta);
      next_y = y - rho_times_turn_dir*(-ctheta + next_ctheta);
    }

    // Use next_x etc because it's checkMotionBackward
    if (!collision_free(next_x, next_y, next_theta, turn_dir, turn_dist, x, y)) {
      // If we collide, return the states and turns we've already accumulated
      reach = false;

      timer_stop = std::chrono::high_resolution_clock::now();
      nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
      collision_check_time += ((double)nanos)/1e9;

      return;
    }

    x = next_x;
    y = next_y;
    theta = next_theta;
    t = next_t;

    ctheta = next_ctheta;
    stheta = next_stheta;

    cur_segment_length += turn_dist;
  }

  reach = true;

  // If we reach, need to add the final state and the final sequence of turns, add states and turns
  turns_vec.push_back(turns.block(0, 0, start_turn_idx + 1, 2));
  states_vec.push_back(Vector4d(next_x, next_y, next_theta, next_t));

  timer_stop = std::chrono::high_resolution_clock::now();
  nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
  collision_check_time += ((double)nanos)/1e9;
}

void DubinsMotionValidatorRects::set_monte_carlo_prop(bool monte_carlo_prop) {
  this->monte_carlo_prop = monte_carlo_prop;
}
