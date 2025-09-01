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
    DubinsControlSampler(const SpaceInformation *si, std::shared_ptr<DubinsMotionValidatorRects> motion_validator, double max_w, double max_delta_t);

    void reset_timing_info();

    double get_collision_check_time();

    ~DubinsControlSampler() override = default;

    bool collision_free(double x, double y, double theta, double turn_dir, double turn_dist, double next_x, double next_y, double rho) const;

    unsigned int sampleTo(Control *control, const ob::State *source, ob::State *dest) override;

    unsigned int sampleTo(Control *control, const Control * /*previous*/, const ob::State *source,
                          ob::State *dest) override;
private:
  ControlSamplerPtr cs_;
  double collision_check_time;
  std::shared_ptr<DubinsMotionValidatorRects> motion_validator;
  double max_w;
  double max_delta_t;
  ompl::RNG rng_;
};

struct DubinsControlSamplerAllocator {
  DubinsControlSamplerAllocator(std::shared_ptr<DubinsMotionValidatorRects> motion_validator, double max_w, double max_delta_t);

  std::shared_ptr<DubinsControlSampler> operator() (const SpaceInformation *si);

  std::shared_ptr<DubinsMotionValidatorRects> motion_validator;

  double max_w;
  double max_delta_t;
};
