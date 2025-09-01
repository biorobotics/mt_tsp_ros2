#pragma once
#include <ompl/base/SpaceInformation.h>
#include <ompl/base/MotionValidator.h>
#include "mt_tsp_ros2/elongate_dubins_path.h"

namespace ob = ompl::base;

class NoTimeDubinsMotionValidator : public ob::MotionValidator {
  public:
    explicit NoTimeDubinsMotionValidator(const ob::SpaceInformationPtr si, double rho);

    bool checkMotion(const ob::State *s1, const ob::State *s2, std::pair<ob::State*, double> &lastValid) const override;

    bool checkMotion(const ob::State *s1, const ob::State *s2) const override;

  protected:
    double rho;
};
