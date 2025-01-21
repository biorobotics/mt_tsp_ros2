#pragma once
#include <ompl/base/spaces/SE2StateSpace.h>
#include <ompl/base/samplers/DeterministicStateSampler.h>

namespace ob = ompl::base;

class SE2StateSpaceDeterministicSampling : public ob::SE2StateSpace {
  public:
    ob::StateSamplerPtr allocStateSampler() const override {
      return std::make_shared<ob::SE2DeterministicStateSampler>(this);
    }
};
