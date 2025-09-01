#pragma once
#include <ompl/base/SpaceInformation.h>
#include <ompl/base/MotionValidator.h>
#include "mt_tsp_ros2/elongate_dubins_path.h"
#include <chrono>

namespace ob = ompl::base;

class DubinsMotionValidator : public ob::MotionValidator {
  public:
    explicit DubinsMotionValidator(const ob::SpaceInformationPtr si, double vmax, double rho);

    virtual RowMatrixXd checkMotionForward(const ob::State *s1, const ob::State *s2, double maxDuration, ob::State *stopState, bool &reach, bool &validMotion, int num_checks = 1000) const;

    virtual RowMatrixXd checkMotionBackward(const ob::State *s1, const ob::State *s2, double maxDuration, ob::State *stopState, bool &reach, bool &validMotion, int num_checks = 1000) const;

    bool checkMotion(const ob::State *s1, const ob::State *s2, std::pair<ob::State*, double> &lastValid) const override;

    bool checkMotion(const ob::State *s1, const ob::State *s2) const override;

    double get_path_elongation_time() const;

    double get_collision_check_time() const;

    void reset_timing_info();
    
    double get_vmax() const;

    int get_num_discarded_samples();

    void reset_num_discarded_samples();

  protected:
    double vmax;
    double rho;
    mutable double path_elongation_time;
    mutable double collision_check_time;
    mutable int num_discarded_samples;
};
