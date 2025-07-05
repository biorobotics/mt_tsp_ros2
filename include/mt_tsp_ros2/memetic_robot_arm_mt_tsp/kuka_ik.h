#pragma once
#include <Eigen/Dense>

using namespace Eigen;

const int dim_q = 7;

const double tol = 1e-8;

Vector3d unit(const Ref<const Vector3d> &v) {
  return v/v.norm();
}

Matrix3d skew(const Ref<const Vector3d> &v) {
  Matrix3d ret;
  ret << 0, -v(2), v(1),
         v(2), 0, -v(0),
         -v(1), v(0), 0;
  return ret;
}

Matrix4d dh_calc(double a, double alpha, double d, double theta) {
  Matrix4d ret;
  double ctheta = cos(theta);
  double stheta = sin(theta);
  double calpha = cos(alpha);
  double salpha = sin(alpha);
  ret << ctheta, -stheta * calpha,  stheta * salpha, a*ctheta,
         stheta,  ctheta * calpha, -ctheta * salpha, a*stheta,
         0.0,     salpha,                    calpha,        d,
         0.0,        0.0,                       0.0,        1;
  return ret;
}

// Untested, converted from python
/*
// Returns false if pose is outside workspace or if robot is at singularity
bool ReferencePlane(Ref<Vector3d> ref_plane_vector, Ref<Matrix3d> rot_base_elbow, Ref<VectorXd> joints, const Ref<const Matrix4d> &pose, const Ref<const Vector3d> &elbow, const Ref<const Vector4d> &l, const Ref<const RowMatrixXd> &dh):
  // REFERENCEPLANE
  // 
  // Calculates the vector normal to the reference plane.
  // A virtual robotic manipulator is created with the same structure as the
  // real robot, but keeping its joint3 = 0.
  //
  // From the current end-effector position and robot configuration,
  // compute the set of joint positions for the virtual manipulator (i.e. theta_3=0)
  // With these we compute the plane of the virtual Shoulder-Elbow-Wrist, which
  // is the reference plane. The NSParam is calculated from the angle between
  // the actual robot Shoulder-Elbow-Wrist plane and the reference plane.
  //
  // Input:  pose            - Homogeneous Matrix size(4,4) of current end-effector pose
  //         rconf           - Robot Configuration 8-bit number
  // Output: ref_plan_vector - Vector normal to the reference plane.
  //         rot_base_elbow  - Rotation Matrix size(3,3) from base to elbow
  //         joints          - Joint values of the virtual robot (debug)
  //
  // The robot configuration parameter is considered because there are usually
  // 2 possible solutions: elbow 'up' and 'down'. Therefore and as the author
  // suggested we are selecting the configuration that matches the current
  // robot.

  if (joints.size() != dim_q) {
    throw std::runtime_error("Joints vector is not correct size");
  }

  // Joint values of virtual manipulator
  joints.setZero();

  const Ref<const Vector3d> &xend = pose.topRightCorner<3, 1>(); // end-effector position from base    
  Vector3d xs0(0, 0, dh(0,2)); // shoulder position from base 
  Vector3d xwt(0, 0, dh(-1,2)); // end-effector position from wrist
  Vector3d xw0 = xend - pose.topLeftCorner<3, 3>()*xwt; // wrist position from base
  Vector3d xsw = xw0 - xs0; // shoulder to wrist vector
  double norm_xsw = xsw.norm();

  double lbs = l(0);
  double lse = l(1); // upper arm length (shoulder to elbow)
  double lew = l(2); // lower arm length (elbow to wrist)

  // Check if pose is within arm+forearm reach
  if (!(norm_xsw < lse + lew && norm_xsw > lse - lew)) {
    return false;
  }

  // -- Joint 4 --
  // Elbow joint can be directly calculated since it only depends on the 
  // robot configuration and the xsw vector 
  if (!(std::abs((norm_xsw*norm_xsw - lse*lse - lew*lew)-(2*lse*lew)) > tol)) {
    // Singularity
    return false;
  }
  // Cosine law - According to our robot, joint 4 rotates backwards
  joints(3) = elbow * acos((norm_xsw*norm_xsw - lse*lse - lew*lew)/(2*lse*lew));

  // Shoulder Joints
  Matrix4d T34 = dh_calc(dh(3,0),dh(3,1),dh(3,2),joints(3));
  const Ref<const Matrix3d> &R34 = T34.topLeftCorner<3, 3>();

  // These are the vectors corresponding to our DH parameters
  Vector3d xse(0, lse, 0);
  Vector3d xew(0, 0, lew);
  // m = member between parentheses. Check equation (14)
  Vector3d m = xse + R34*xew;

  // -- Joint 1 --
  // Since joint3 is locked as 0, the only joint defining the orientation of
  // the xsw vector in the xy-plane is joint 1. Therefore and since we are
  // only interested in the transformation T03 (disregarding joint limits), we
  // chose to simply set joint 1 as the atan of xsw y and x coordinates 
  // (even if if goes beyond the joint limit).

  // Cannot be this because if x and y are 0, then it is not defined.
  if (xsw.cross(Vector3d(0, 0, 1)).norm() > tol) {
    joints(0) = atan2(xsw(1),xsw(0));
  } else {
    joints(0) = 0;
  }

  // -- Joint 2 --
  // Can be found through geometric relations
  // Let phi be the angle E-S-W, and theta2 the angle (z-axis)-S-E.
  // Then, theta2 = atan2(r,xsw(3)) -/+ phi.
  // phi can be calculated as a function of theta3:
  //   atan2(lew*sin(theta4),lse+lew*cos(theta4))
  // z-axis
  //   ^
  //   |  E O------------O W
  //   |   /        .  
  //   |  /      .
  //   | /    .    xsw
  //   |/  .
  // S O___________________ r-axis
  //
  double r = sqrt(xsw(0)*xsw(0) + xsw(1)*xsw(1));
  double dsw = xsw.norm();
  double phi = acos((lse*lse + dsw*dsw - lew*lew)/(2*lse*dsw));

  joints(1) = atan2(r, xsw(2)) + elbow * phi;

  // Lower arm transformation
  Matrix4d T01 = dh_calc(dh(0,0),dh(0,1),dh(0,2),joints(0));
  Matrix4d T12 = dh_calc(dh(1,0),dh(1,1),dh(1,2),joints(1));
  Matrix4d T23 = dh_calc(dh(2,0),dh(2,1),dh(2,2),0);
  Matrix4d T34 = dh_calc(dh(3,0),dh(3,1),dh(3,2),joints(3));
  Matrix4d T04 = T01*T12*T23*T34;

  rot_base_elbow = T01.topLeftCorner<3, 3>()*T12.topLeftCorner<3, 3>()*T23.topLeftCorner<3, 3>();

  // With T03 we can calculate the reference elbow position and with it the
  // vector normal to the reference plane.
  const Ref<const Vector3d> &x0e = T04.topRightCorner<3, 1>() // reference elbow position
  Vector3d v1 = unit(x0e - xs0); // unit vector from shoulder to elbow
  Vector3d v2 = unit(xw0 - xs0); // unit vector from shoulder to wrist

  ref_plane_vector = v1.cross(v2);
}

VectorXd kuka_ik(const Ref<const Matrix4d> &pose, double nsparam, int rconf, const Ref<const Vector4d> &l, const Ref<const RowMatrixXd> &dh) {
  int arm = rconf >> 3;
  int elbow = (rconf >> 2) & 1;
  int wrist = rconf & 1;

  VectorXd joints(dim_q);

  const Ref<const Vector3d> &xend = pose.topRightCorner<3, 1>(); // end-effector position from base    
  Vector3d xs(0, 0, dh(0,2)); // shoulder position from base 
  Vector3d xwt(0, 0, dh(-1,2)); // end-effector position from wrist
  Vector3d xw = xend - pose.topLeftCorner<3, 3>()*xwt; // wrist position from base
  Vector3d xsw = xw - xs; // shoulder to wrist vector
  double norm_xsw = xsw.norm();
  if (norm_xsw == 0) {
    return VectorXd::Zero(0);
  }
  Vector3d usw = xsw/norm_xsw;

  double lbs = l(0);
  double lse = l(1); // upper arm length (shoulder to elbow)
  double lew = l(2); // lower arm length (elbow to wrist)

  // Check if pose is within arm+forearm reach
  if (!(norm_xsw < lse + lew && norm_xsw > lse - lew)) {
    return VectorXd::Zero(0);
  }

  // -- Joint 4 --
  // Elbow joint can be directly calculated since it does only depend on the 
  // robot configuration and the xsw vector 
  if (!(std::abs((norm_xsw*norm_xsw - lse*lse - lew*lew)-(2*lse*lew)) > tol)) {
    return VectorXd::Zero(0);
  }
  // Cosine law - According to our robot, joint 4 rotates backwards
  joints(3) = elbow * acos((norm_xsw*norm_xsw - lse*lse - lew*lew)/(2*lse*lew));

  Matrix4d T34 = dh_calc(dh(3,0),dh(3,1),dh(3,2),joints(3));
  const Ref<const Matrix3d> &R34 = T34.topLeftCorner<3, 3>();

  // Shoulder Joints
  // First compute the reference joint angles when the arm angle is zero.
  Matrix3d R03_o;
  Vector3d tmp_ref_plane_vector;
  VectorXd tmp_joints;
  ReferencePlane(tmp_ref_plane_vector, R03_o, tmp_joints, pose, elbow, l, dh);

  Matrix3d skew_usw = skew(usw);
  // Following eq. (15), the auxiliary matrixes As Bs and Cs can be calculated 
  // by substituting eq. (6) into (9). 
  // R0psi = I3 + sin(psi)*skew_usw + (1-cos(psi))*skew_usw²    (6)
  // R03 = R0psi * R03_o                                         (9)
  // Substituting (distributive prop.) we get:
  // R03 = R03_o*skew_usw*sin(psi) + R03_o*(-skew_usw²)*cos(psi) + R03_o(I3 + skew_usw²)
  // R03 =      As       *sin(psi) +        Bs         *cos(psi) +          Cs
  Matrix3d As = skew_usw * R03_o;
  Matrix3d Bs = -skew_usw*skew_usw * R03_o;
  Matrix3d Cs = usw*usw.transpose() * R03_o;

  double psi = nsparam;
  Matrix3d  R03 = As*sin(psi) + Bs*cos(psi) + Cs;

  //  T03 transformation matrix (DH parameters)
  // [ cos(j1)*cos(j2)*cos(j3) - sin(j1)*sin(j3), cos(j1)*sin(j2), cos(j3)*sin(j1) + cos(j1)*cos(j2)*sin(j3), 0.4*cos(j1)*sin(j2)]
  // [ cos(j1)*sin(j3) + cos(j2)*cos(j3)*sin(j1), sin(j1)*sin(j2), cos(j2)*sin(j1)*sin(j3) - cos(j1)*cos(j3), 0.4*sin(j1)*sin(j2)]
  // [                          -cos(j3)*sin(j2),         cos(j2),                          -sin(j2)*sin(j3),  0.4*cos(j2) + 0.34]
  // [                                         0,               0,                                         0,                   1]
  joints(0) = atan2(arm * R03(1,1), arm * R03(0,1));
  joints(1) = arm * acos(R03(2,1));

  joints(2) = atan2(arm * -R03(2,2), arm * -R03(2,0));

  Matrix3d Aw = R34.transpose() * As.transpose() * pose.topLeftCorner<3, 3>();
  Matrix3d Bw = R34.transpose() * Bs.transpose() * pose.topLeftCorner<3, 3>();
  Matrix3d Cw = R34.transpose() * Cs.transpose() * pose.topLeftCorner<3, 3>();

  Matrix3d R47 = Aw*sin(psi) + Bw*cos(psi) + Cw;

  //  T47 transformation matrix (DH parameters)
  // [ cos(j5)*cos(j6)*cos(j7) - sin(j5)*sin(j7), - cos(j7)*sin(j5) - cos(j5)*cos(j6)*sin(j7), cos(j5)*sin(j6), (63*cos(j5)*sin(j6))/500]
  // [ cos(j5)*sin(j7) + cos(j6)*cos(j7)*sin(j5),   cos(j5)*cos(j7) - cos(j6)*sin(j5)*sin(j7), sin(j5)*sin(j6), (63*sin(j5)*sin(j6))/500]
  // [                          -cos(j7)*sin(j6),                             sin(j6)*sin(j7),         cos(j6),   (63*cos(j6))/500 + 2/5]
  // [                                         0,                                           0,               0,                        1]
  joints(4) = atan2(wrist * R47(1,2), wrist * R47(0,2));
  joints(5) = wrist * acos(R47(2,2));
  joints(6) = atan2(wrist * R47(2,1), wrist * -R47(2,0));

  return joints;
}
*/

const Matrix4d kuka_fk(const Ref<const VectorXd> &joints, const Ref<const Vector4d> &l, const Ref<const RowMatrixXd> &dh) {
  int nj = dh.rows();

  // Assign joint values to the theta column of the DH parameters
  MatrixXd dh_tmp = dh;
  dh_tmp.col(3) = joints;
  // Store transformations from the base reference frame to the index joint
  // e.g: tr(:,:,2) is the T02 -> transformation from base to joint 2 (DH table)
  std::vector<Matrix4d> tr(nj);
  // Rotation Matrix applied with Denavit-Hartenberg parameters [same as (3)]
  // R = [Xx,Yx,Zx,   --  Xx = cos(theta), Yx = -sin(theta) * cos(alpha), Zx =  sin(theta) * sin(alpha)
  //      Xy,YY,Zy,   --  Xy = sin(theta), Yy =  cos(theta) * cos(alpha), Zy = -cos(theta) * sin(alpha)
  //      Xz,Yz,Zz];  --  Xz = 0.0,        Yz =  sin(alpha),              Zz =  cos(alpha)   

  for (int i = 0; i < nj; ++i) {
    double a = dh_tmp(i,0);
    double alpha = dh_tmp(i,1);
    double d = dh_tmp(i,2);
    double theta = dh_tmp(i,3);
    Matrix4d tmp = dh_calc(a, alpha, d, theta);
    if (i == 0) {
      tr[0] = tmp;
    } else {
      tr[i] = tr[i - 1] * tmp;
    }
  }

  return tr.back(); // end-effector transformation from base
}
