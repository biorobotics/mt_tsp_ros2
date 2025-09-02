#include "mt_tsp_ros2/time_constrained_dubins_planning/yao2020/ic_sliding.h"

ICSlidingClass::ICSlidingClass(RowMatrixXdRef_const turns, double x_0, double y_0, double theta_0, double x_f, double y_f, double theta_f, double rho, bool direction) : turns(turns), x_0(x_0), y_0(y_0), theta_0(theta_0), x_f(x_f), y_f(y_f), theta_f(theta_f), rho(rho), direction(direction) {
  slide_idx = direction ? turns.rows() - 2 : 1;
  if (turns(slide_idx, 0) == 0) {
    throw std::runtime_error("Cannot slide an S segment");
  }

  /*
  if (turns.rows() < 2) {
    throw std::runtime_error("Cannot slide if we do not have at least two segments");
  }
  */

  if (turns.rows() != 4) {
    throw std::runtime_error("Cannot slide if we do not have at exactly four segments");
  }

  turn_to_shorten = direction ? slide_idx + 1 : slide_idx - 1;

  path_type_options.push_back(DubinsPathType::LSL);
  path_type_options.push_back(DubinsPathType::LSR);
  path_type_options.push_back(DubinsPathType::RSR);
  path_type_options.push_back(DubinsPathType::RSL);
  path_type_options.push_back(DubinsPathType::LRL);
  path_type_options.push_back(DubinsPathType::RLR);
}

RowMatrixXd ICSlidingClass::slide(double amount, double prev_length) {
  // Compute configuration after/before shortened segment (after if sliding backward, before if sliding forward)
  double x_mid;
  double y_mid;
  double theta_mid;
  double dist = turns(turn_to_shorten, 1) - amount;
  if (turns(turn_to_shorten, 0) == 0) {
    // S segment
    if (direction) {
      // Sliding forward
      double ctheta_f = cos(theta_f);
      double stheta_f = sin(theta_f);
      x_mid = x_f - ctheta_f*dist;
      y_mid = y_f - stheta_f*dist;
      theta_mid = theta_f;
    } else {
      // Sliding backward
      double ctheta_0 = cos(theta_0);
      double stheta_0 = sin(theta_0);
      x_mid = x_0 + ctheta_0*dist;
      y_mid = y_0 + stheta_0*dist;
      theta_mid = theta_0;
    }
  } else {
    // C segment
    double C_sign = turns(turn_to_shorten, 0);
    if (direction) {
      // Sliding forward
      double c_f = cos(theta_f);
      double s_f = sin(theta_f);
      theta_mid = theta_f - C_sign*dist/rho;
      double c_mid = cos(theta_mid);
      double s_mid = sin(theta_mid);
      x_mid = x_f - rho/C_sign*(-s_mid + s_f);
      y_mid = y_f - rho/C_sign*(c_mid - c_f);
    } else {
      // Sliding backward
      double c_0 = cos(theta_0);
      double s_0 = sin(theta_0);
      theta_mid = theta_0 + C_sign*dist/rho;
      double c_mid = cos(theta_mid);
      double s_mid = sin(theta_mid);
      x_mid = x_0 + rho/C_sign*(-s_0 + s_mid);
      y_mid = y_0 + rho/C_sign*(c_0 - c_mid);
    }
  }

  // Compute Dubins path
  RowMatrixXd ret(4, 2);
  ret(turn_to_shorten, 0) = turns(turn_to_shorten, 0);
  ret(turn_to_shorten, 1) = dist;
  double q_0_for_dubins[3];
  double q_f_for_dubins[3];

  DubinsPath closest_path_among_options;
  double closest_path_length_abs_diff_among_options = std::numeric_limits<double>::infinity();

  if (direction) {
    // Sliding forward
    q_0_for_dubins[0] = x_0;
    q_0_for_dubins[1] = y_0;
    q_0_for_dubins[2] = theta_0;

    q_f_for_dubins[0] = x_mid;
    q_f_for_dubins[1] = y_mid;
    q_f_for_dubins[2] = theta_mid;

    for (auto path_type : path_type_options) {
      DubinsPath path;
      int dubins_error = dubins_path(&path, q_0_for_dubins, q_f_for_dubins, rho, path_type);
      if (dubins_error) {
        continue;
      }
      double l = dubins_path_length(&path) + dist;
      double path_length_abs_diff = std::abs(l - prev_length);
      if (path_length_abs_diff < closest_path_length_abs_diff_among_options) {
        closest_path_length_abs_diff_among_options = path_length_abs_diff;
        closest_path_among_options = path;
      }

      if (path.type == DubinsPathType::LSL) {
        std::cout << "LSL ";
      } else if (path.type == DubinsPathType::LSR) {
        std::cout << "LSR ";
      } else if (path.type == DubinsPathType::RSR) {
        std::cout << "RSR ";
      } else if (path.type == DubinsPathType::RSL) {
        std::cout << "RSL ";
      } else if (path.type == DubinsPathType::RLR) {
        std::cout << "RLR ";
      } else if (path.type == DubinsPathType::LRL) {
        std::cout << "LRL ";
      }
      std::cout << l << " ";
    }
    std::cout << " prev length " << prev_length << std::endl;

    if (closest_path_among_options.type == DubinsPathType::LSL) {
      ret(0, 0) = 1.;
      ret(1, 0) = 0.;
      ret(2, 0) = 1.;
    } else if (closest_path_among_options.type == DubinsPathType::LSR) {
      ret(0, 0) = 1.;
      ret(1, 0) = 0.;
      ret(2, 0) = -1.;
    } else if (closest_path_among_options.type == DubinsPathType::RSR) {
      ret(0, 0) = -1.;
      ret(1, 0) = 0.;
      ret(2, 0) = -1.;
    } else if (closest_path_among_options.type == DubinsPathType::RSL) {
      ret(0, 0) = -1.;
      ret(1, 0) = 0.;
      ret(2, 0) = 1.;
    } else if (closest_path_among_options.type == DubinsPathType::RLR) {
      ret(0, 0) = -1.;
      ret(1, 0) = 1.;
      ret(2, 0) = -1.;
    } else if (closest_path_among_options.type == DubinsPathType::LRL) {
      ret(0, 0) = 1.;
      ret(1, 0) = -1.;
      ret(2, 0) = 1.;
    }
    ret(0, 1) = closest_path_among_options.param[0]*rho;
    ret(1, 1) = closest_path_among_options.param[1]*rho;
    ret(2, 1) = closest_path_among_options.param[2]*rho;
  } else {
    // Sliding backward
    q_0_for_dubins[0] = x_mid;
    q_0_for_dubins[1] = y_mid;
    q_0_for_dubins[2] = theta_mid;

    q_f_for_dubins[0] = x_f;
    q_f_for_dubins[1] = y_f;
    q_f_for_dubins[2] = theta_f;

    for (auto path_type : path_type_options) {
      DubinsPath path;
      int dubins_error = dubins_path(&path, q_0_for_dubins, q_f_for_dubins, rho, path_type);
      if (dubins_error) {
        continue;
      }
      double l = dubins_path_length(&path) + dist;
      double path_length_abs_diff = std::abs(l - prev_length);
      if (path_length_abs_diff < closest_path_length_abs_diff_among_options) {
        closest_path_length_abs_diff_among_options = path_length_abs_diff;
        closest_path_among_options = path;
      }
    }

    if (closest_path_among_options.type == DubinsPathType::LSL) {
      ret(1, 0) = 1.;
      ret(2, 0) = 0.;
      ret(3, 0) = 1.;
    } else if (closest_path_among_options.type == DubinsPathType::LSR) {
      ret(1, 0) = 1.;
      ret(2, 0) = 0.;
      ret(3, 0) = -1.;
    } else if (closest_path_among_options.type == DubinsPathType::RSR) {
      ret(1, 0) = -1.;
      ret(2, 0) = 0.;
      ret(3, 0) = -1.;
    } else if (closest_path_among_options.type == DubinsPathType::RSL) {
      ret(1, 0) = -1.;
      ret(2, 0) = 0.;
      ret(3, 0) = 1.;
    } else if (closest_path_among_options.type == DubinsPathType::RLR) {
      ret(1, 0) = -1.;
      ret(2, 0) = 1.;
      ret(3, 0) = -1.;
    } else if (closest_path_among_options.type == DubinsPathType::LRL) {
      ret(1, 0) = 1.;
      ret(2, 0) = -1.;
      ret(3, 0) = 1.;
    }
    ret(1, 1) = closest_path_among_options.param[0]*rho;
    ret(2, 1) = closest_path_among_options.param[1]*rho;
    ret(3, 1) = closest_path_among_options.param[2]*rho;
  }
  return ret;
}

void ICSlidingClass::batch_slide(std::vector<RowMatrixXd> &turns_seq, const Ref<const VectorXd> &amounts) {
  turns_seq.resize(amounts.size());
  double prev_length = turns.col(1).sum();
  for (int amount_idx = 0; amount_idx < amounts.size(); ++amount_idx) {
    turns_seq[amount_idx] = slide(amounts[amount_idx], prev_length);
    prev_length = turns_seq[amount_idx].col(1).sum();
  }
}
