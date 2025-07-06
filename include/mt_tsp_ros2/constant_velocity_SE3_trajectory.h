#pragma once
#include "mt_tsp_ros2/constant_velocity_rotation_trajectory.h"
#include "mt_tsp_ros2/linear_trajectory.h"
#include "mt_tsp_ros2/SE3_spline.h"

class ConstantVelocitySE3Trajectory {
  public:
    ConstantVelocitySE3Trajectory(const SE3Spline &spline, double t) : pos_trj(spline.get_pos_spline(), t), rot_trj(spline.get_rot_spline(), t) {
    }

    virtual Matrix4d operator()(double t) const {
      Matrix4d ret = Matrix4d::Identity();
      ret.topLeftCorner<3, 3>() = rot_trj(t);
      ret.topRightCorner<3, 1>() = pos_trj(t);
      return ret;
    }

  private:
    LinearTrajectory pos_trj;
    ConstantVelocityRotationTrajectory rot_trj;
};
