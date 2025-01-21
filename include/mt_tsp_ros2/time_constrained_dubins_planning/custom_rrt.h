#include "ompl/geometric/planners/rrt/RRT.h"

using namespace ompl::geometric;
namespace base = ompl::base;

class CustomRRT : public RRT {
  public:
    CustomRRT(const base::SpaceInformationPtr &si) : RRT(si) {
      rng_ = ompl::RNG(1);
    }
};
