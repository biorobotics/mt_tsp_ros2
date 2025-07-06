#pragma once
#include "mt_tsp_ros2/rotation_spline.h"
#include "mt_tsp_ros2/cpp_spline_3d.h"

class SE3Spline {
  public:
    SE3Spline(const CppSpline3D &pos_spline, const RotationSpline &rot_spline) : pos_spline(pos_spline), rot_spline(rot_spline) {
      constant_orientation = false;
    }

    void set_constant_orientation(bool constant_orientation, const Ref<const Matrix3d> &rmat) {
      this->constant_orientation = constant_orientation;
      if (constant_orientation) {
        constant_rmat = rmat;
      }
    }

    virtual Matrix4d operator()(double t) const {
      Matrix4d ret = Matrix4d::Identity();
      ret.topLeftCorner<3, 3>() = constant_orientation ? constant_rmat : rot_spline(t);
      ret.topRightCorner<3, 1>() = pos_spline(t);
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

    bool constant_orientation;
    Matrix3d constant_rmat;
};
