#pragma once
#include <ompl/base/spaces/SE2StateSpace.h>
#include <ompl/base/spaces/DubinsStateSpace.h>
#include <ompl/base/samplers/DeterministicStateSampler.h>

namespace ob = ompl::base;

class DubinsStateSpaceDeterministicSampling : public ob::DubinsStateSpace {
  public:
    DubinsStateSpaceDeterministicSampling(double rho) : ob::DubinsStateSpace(rho) {
    }

    ob::StateSamplerPtr allocStateSampler() const override {
      return std::make_shared<ob::SE2DeterministicStateSampler>(this);
    }
};
