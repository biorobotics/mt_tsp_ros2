#pragma once
#include <ompl/base/SpaceInformation.h>
#include <ompl/base/spaces/SE2StateSpace.h>

using namespace Eigen;

namespace ob = ompl::base;

typedef Matrix<bool, Dynamic, Dynamic, RowMajor> RowMatrixXb;
typedef const Ref<const RowMatrixXb> &RowMatrixXbRef_const;
typedef const Ref<const Vector2d> &Vector2dRef_const;

class DubinsStateValidityChecker : public ob::StateValidityChecker {
  public:
    DubinsStateValidityChecker(const ob::SpaceInformationPtr &si, RowMatrixXbRef_const &occupancy, Vector2dRef_const map_lb, Vector2dRef_const map_ub);
 
    virtual bool isValid(const ob::State *state) const;

  private:
    RowMatrixXb occupancy;
    Vector2d map_lb;
    Vector2d map_ub;
    Vector2d cell_size;
};
