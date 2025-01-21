#include "ompl/geometric/planners/rrt/RRTstar.h"

using namespace ompl::geometric;
namespace base = ompl::base;

class CustomRRTstar : public RRTstar {
  public:
    CustomRRTstar(const base::SpaceInformationPtr &si) : RRTstar(si) {
      rng_ = ompl::RNG(1);
    }
};
