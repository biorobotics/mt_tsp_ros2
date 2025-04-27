#pragma once

#include <ompl/control/DirectedControlSampler.h>
#include <ompl/control/StatePropagator.h>
#include <ompl/control/SpaceInformation.h>
#include <cmath>
#include <ompl/control/spaces/RealVectorControlSpace.h>
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_motion_validator_rects.h"

using namespace ompl::control;
 
class DubinsControlSampler : public DirectedControlSampler
{
public:
    DubinsControlSampler(const SpaceInformation *si, std::shared_ptr<DubinsMotionValidatorRects> motion_validator, double max_w, double max_delta_t) : DirectedControlSampler(si), motion_validator(motion_validator), max_w(max_w), max_delta_t(max_delta_t)
    {
      collision_check_time = 0.;
    }

    void reset_timing_info() {
      collision_check_time = 0.;
    }

    double get_collision_check_time() {
      return collision_check_time;
    }

    ~DubinsControlSampler() override = default;

    bool collision_free(double x, double y, double theta, double turn_dir, double turn_dist, double next_x, double next_y, double rho) const {
      // S segment
      if (turn_dir == 0) {
        Point point1(x, y);
        Point point2(next_x, next_y);
        Segment segment(point1, point2);
        return !motion_validator->get_aabb_tree().do_intersect(segment);
      }

      // C segment
      double c = cos(theta);
      double s = sin(theta);

      Vector2d dir(c, s);
      Vector2d perp(-s*turn_dir, c*turn_dir);
      Vector2d center = Vector2d(x, y) + perp*rho;

      CGAL::Bbox_2 bbox(center(0) - rho, center(1) - rho, center(0) + rho, center(1) + rho);
      std::list<Primitive_id> primitives;
      motion_validator->get_aabb_tree().all_intersected_primitives(bbox, std::back_inserter(primitives));

      for (Primitive_id primitive : primitives) {
        if (arc_intersects_line_segment(x, y, theta, turn_dir, turn_dist/rho, rho,
                                        primitive->source().x(), primitive->source().y(), primitive->target().x(), primitive->target().y())) {
          return false;
        }
      }
      return true;
    }

    unsigned int sampleTo(Control *control, const base::State *source, base::State *dest) override
    {
        RealVectorControlSpace::ControlType *real_vector_control = control->as<RealVectorControlSpace::ControlType>();
        double w = rng_.uniformReal(-max_w, max_w);
        double delta_t = rng_.uniformReal(0., max_delta_t);
        if (std::abs(w) < 1e-4) {
          // Avoiding potential numerical issues
          w = 0.;
        }
        (*real_vector_control)[0] = w;
        (*real_vector_control)[1] = delta_t;

        double x = source->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getX();
        double y = source->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getY();
        double theta = source->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getYaw();
        double t = source->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position;

        double ctheta = cos(theta);
        double stheta = sin(theta);

        auto timer_start = std::chrono::high_resolution_clock::now();

        double turn_dir = (0 < w) - (w < 0);
        double turn_dist = motion_validator->get_vmax()*delta_t;
        // rho here is the turning radius determined by the sampled
        // angular velocity, not the min turning radius
        double rho_times_turn_dir = motion_validator->get_vmax()/w;

        double next_x;
        double next_y;
        double next_theta;
        double next_ctheta;
        double next_stheta;
        double next_t;

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

        next_t = t + delta_t;

        if (!collision_free(x, y, theta, turn_dir, turn_dist, next_x, next_y, std::abs(rho_times_turn_dir))) {
          auto timer_stop = std::chrono::high_resolution_clock::now();
          auto nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
          collision_check_time += ((double)nanos)/1e9;
          return 0;
        }

        dest->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setX(next_x);
        dest->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setY(next_y);
        dest->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setYaw(next_theta);
        dest->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position = next_t;
        std::cout << "reached state with t = " << next_t << std::endl;
        auto timer_stop = std::chrono::high_resolution_clock::now();
        auto nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
        collision_check_time += ((double)nanos)/1e9;

        return 1;
    }

    unsigned int sampleTo(Control *control, const Control * /*previous*/, const base::State *source,
                          base::State *dest) override
    {
        return sampleTo(control, source, dest);
    }
private:
  ControlSamplerPtr cs_;
  double collision_check_time;
  std::shared_ptr<DubinsMotionValidatorRects> motion_validator;
  double max_w;
  double max_delta_t;
  ompl::RNG rng_;
};

struct DubinsControlSamplerAllocator {
  DubinsControlSamplerAllocator(std::shared_ptr<DubinsMotionValidatorRects> motion_validator, double max_w, double max_delta_t) : motion_validator(motion_validator), max_w(max_w), max_delta_t(max_delta_t) {
  }

  std::shared_ptr<DubinsControlSampler> operator() (const SpaceInformation *si) {
    return std::make_shared<DubinsControlSampler>(si, motion_validator, max_w, max_delta_t);
  }

  std::shared_ptr<DubinsMotionValidatorRects> motion_validator;

  double max_w;
  double max_delta_t;
};
