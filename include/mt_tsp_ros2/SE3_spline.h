#pragma once
#include "mt_tsp_ros2/rotation_spline.h"
#include "mt_tsp_ros2/cpp_spline_3d.h"
#include "mt_tsp_ros2/linear_trajectory.h"
#include "mt_tsp_ros2/constant_velocity_rotation_trajectory.h"

class SE3Spline {
  public:
    SE3Spline(const CppSpline3D &pos_spline, const RotationSpline &rot_spline) : pos_spline(pos_spline), rot_spline(rot_spline) {
      do_const_vel_pos_trj = false;
      do_const_vel_rot_trj = false;
    }

    void set_const_vel_pos_trj(bool do_const_vel_pos_trj, const Ref<const Vector3d> &p0, const Ref<const Vector3d> &v) {
      this->do_const_vel_pos_trj = do_const_vel_pos_trj;
      if (do_const_vel_pos_trj) {
        this->const_vel_pos_trj = LinearTrajectory(p0, v);
      }
    }

    void set_const_vel_pos_trj(bool do_const_vel_pos_trj, double t) {
      this->do_const_vel_pos_trj = do_const_vel_pos_trj;
      if (do_const_vel_pos_trj) {
        this->const_vel_pos_trj = LinearTrajectory(pos_spline, t);
      }
    }

    void set_const_vel_rot_trj(bool do_const_vel_rot_trj, const Ref<const Matrix3d> &R0, const Ref<const Vector3d> &axis, double w) {
      this->do_const_vel_rot_trj = do_const_vel_rot_trj;
      if (do_const_vel_rot_trj) {
        this->const_vel_rot_trj = ConstantVelocityRotationTrajectory(R0, axis, w);
      }
    }

    void set_const_vel_rot_trj(bool do_const_vel_rot_trj, double t) {
      this->do_const_vel_rot_trj = do_const_vel_rot_trj;
      if (do_const_vel_rot_trj) {
        this->const_vel_rot_trj = ConstantVelocityRotationTrajectory(rot_spline, t);
      }
    }

    virtual Matrix4d operator()(double t) const {
      Matrix4d ret = Matrix4d::Identity();
      ret.topLeftCorner<3, 3>() = do_const_vel_rot_trj ? const_vel_rot_trj(t) : rot_spline(t);
      ret.topRightCorner<3, 1>() = do_const_vel_pos_trj ? const_vel_pos_trj(t).head<3>() : pos_spline(t);
      return ret;
    }

    const CppSpline3D &get_pos_spline() const {
      return pos_spline;
    }

    const RotationSpline &get_rot_spline() const {
      return rot_spline;
    }

  private:
    CppSpline3D pos_spline;
    RotationSpline rot_spline;

    bool do_const_vel_pos_trj;
    LinearTrajectory const_vel_pos_trj;

    bool do_const_vel_rot_trj;
    ConstantVelocityRotationTrajectory const_vel_rot_trj;
;
};
