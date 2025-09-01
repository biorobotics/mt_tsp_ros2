#include "mt_tsp_ros2/time_constrained_dubins_planning/yao2020/ic_sliding.h"

// Returns a sequence of (turn direction, dist) pairs.
// turns should have two or more rows.
// amount is how much to shorten the segment in the sliding direction.
// direction = true means slide forward, false means slide backward.
// Performs sliding on segment turns.rows() - 2 (zero-based indexing) if forward, and on segment 1 if backward
RowMatrixXd ic_sliding(RowMatrixXdRef_const turns, double x_0, double y_0, double theta_0, double x_f, double y_f, double theta_f, double rho, double amount, bool direction) {
  int slide_idx = direction ? turns.rows() - 2 : 1;
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

  int turn_to_shorten = direction ? slide_idx + 1 : slide_idx - 1;
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
  ret.col(0) = turns.col(0);
  ret(turn_to_shorten, 1) = dist;
  double q_0_for_dubins[3];
  double q_f_for_dubins[3];
  std::unordered_set<DubinsPathType> path_type_options;
  path_type_options.insert(DubinsPathType::LSL);
  path_type_options.insert(DubinsPathType::LSR);
  path_type_options.insert(DubinsPathType::RSR);
  path_type_options.insert(DubinsPathType::RSL);
  path_type_options.insert(DubinsPathType::LRL);
  path_type_options.insert(DubinsPathType::RLR);
  DubinsPath shortest_path_among_options;
  double shortest_path_length_among_options = std::numeric_limits<double>::infinity();

  if (direction) {
    // Sliding forward
    q_0_for_dubins[0] = x_0;
    q_0_for_dubins[1] = y_0;
    q_0_for_dubins[2] = theta_0;

    q_f_for_dubins[0] = x_mid;
    q_f_for_dubins[1] = y_mid;
    q_f_for_dubins[2] = theta_mid;

    if (std::abs(turns(0, 1)) > 1e-10) {
      if (turns(0, 0) != -1) {
        path_type_options.erase(DubinsPathType::RSR);
        path_type_options.erase(DubinsPathType::RSL);
        path_type_options.erase(DubinsPathType::RLR);
      }

      if (turns(0, 0) != 1) {
        path_type_options.erase(DubinsPathType::LSR);
        path_type_options.erase(DubinsPathType::LSL);
        path_type_options.erase(DubinsPathType::LRL);
      }
    }

    if (std::abs(turns(1, 1)) > 1e-10) {
      if (turns(1, 0) != -1) {
        path_type_options.erase(DubinsPathType::LRL);
      }

      if (turns(1, 0) != 0) {
        path_type_options.erase(DubinsPathType::LSL);
        path_type_options.erase(DubinsPathType::LSR);
        path_type_options.erase(DubinsPathType::RSR);
        path_type_options.erase(DubinsPathType::RSL);
      }

      if (turns(1, 0) != 1) {
        path_type_options.erase(DubinsPathType::RLR);
      }
    }

    if (std::abs(turns(2, 1)) > 1e-10) {
      if (turns(2, 0) != -1) {
        path_type_options.erase(DubinsPathType::RSR);
        path_type_options.erase(DubinsPathType::LSR);
        path_type_options.erase(DubinsPathType::RLR);
      }

      if (turns(2, 0) != 1) {
        path_type_options.erase(DubinsPathType::LSL);
        path_type_options.erase(DubinsPathType::RSL);
        path_type_options.erase(DubinsPathType::LRL);
      }
    }

    for (auto path_type : path_type_options) {
      DubinsPath path;
      int dubins_error = dubins_path(&path, q_0_for_dubins, q_f_for_dubins, rho, path_type);
      double l = dubins_path_length(&path);
      if (l < shortest_path_length_among_options) {
        shortest_path_length_among_options = l;
        shortest_path_among_options = path;
      }
    }

    if (shortest_path_among_options.type == DubinsPathType::LSL) {
      ret(0, 0) = 1.;
      ret(1, 0) = 0.;
      ret(2, 0) = 1.;
    } else if (shortest_path_among_options.type == DubinsPathType::LSR) {
      ret(0, 0) = 1.;
      ret(1, 0) = 0.;
      ret(2, 0) = -1.;
    } else if (shortest_path_among_options.type == DubinsPathType::RSR) {
      ret(0, 0) = -1.;
      ret(1, 0) = 0.;
      ret(2, 0) = -1.;
    } else if (shortest_path_among_options.type == DubinsPathType::RSL) {
      ret(0, 0) = -1.;
      ret(1, 0) = 0.;
      ret(2, 0) = 1.;
    } else if (shortest_path_among_options.type == DubinsPathType::RLR) {
      ret(0, 0) = -1.;
      ret(1, 0) = 1.;
      ret(2, 0) = -1.;
    } else if (shortest_path_among_options.type == DubinsPathType::LRL) {
      ret(0, 0) = 1.;
      ret(1, 0) = -1.;
      ret(2, 0) = 1.;
    }
    ret(0, 1) = shortest_path_among_options.param[0]*rho;
    ret(1, 1) = shortest_path_among_options.param[1]*rho;
    ret(2, 1) = shortest_path_among_options.param[2]*rho;
  } else {
    // Sliding backward
    q_0_for_dubins[0] = x_mid;
    q_0_for_dubins[1] = y_mid;
    q_0_for_dubins[2] = theta_mid;

    q_f_for_dubins[0] = x_f;
    q_f_for_dubins[1] = y_f;
    q_f_for_dubins[2] = theta_f;

    if (std::abs(turns(1, 1)) > 1e-10) {
      if (turns(1, 0) != -1) {
        path_type_options.erase(DubinsPathType::RSR);
        path_type_options.erase(DubinsPathType::RSL);
        path_type_options.erase(DubinsPathType::RLR);
      }

      if (turns(1, 0) != 1) {
        path_type_options.erase(DubinsPathType::LSR);
        path_type_options.erase(DubinsPathType::LSL);
        path_type_options.erase(DubinsPathType::LRL);
      }
    }

    if (std::abs(turns(2, 1)) > 1e-10) {
      if (turns(2, 0) != -1) {
        path_type_options.erase(DubinsPathType::LRL);
      }

      if (turns(2, 0) != 0) {
        path_type_options.erase(DubinsPathType::LSL);
        path_type_options.erase(DubinsPathType::LSR);
        path_type_options.erase(DubinsPathType::RSR);
        path_type_options.erase(DubinsPathType::RSL);
      }

      if (turns(2, 0) != 1) {
        path_type_options.erase(DubinsPathType::RLR);
      }
    }

    if (std::abs(turns(3, 1)) > 1e-10) {
      if (turns(3, 0) != -1) {
        path_type_options.erase(DubinsPathType::RSR);
        path_type_options.erase(DubinsPathType::LSR);
        path_type_options.erase(DubinsPathType::RLR);
      }

      if (turns(3, 0) != 1) {
        path_type_options.erase(DubinsPathType::LSL);
        path_type_options.erase(DubinsPathType::RSL);
        path_type_options.erase(DubinsPathType::LRL);
      }
    }

    for (auto path_type : path_type_options) {
      DubinsPath path;
      int dubins_error = dubins_path(&path, q_0_for_dubins, q_f_for_dubins, rho, path_type);
      double l = dubins_path_length(&path);
      if (l < shortest_path_length_among_options) {
        shortest_path_length_among_options = l;
        shortest_path_among_options = path;
      }
    }

    if (shortest_path_among_options.type == DubinsPathType::LSL) {
      ret(1, 0) = 1.;
      ret(2, 0) = 0.;
      ret(3, 0) = 1.;
    } else if (shortest_path_among_options.type == DubinsPathType::LSR) {
      ret(1, 0) = 1.;
      ret(2, 0) = 0.;
      ret(3, 0) = -1.;
    } else if (shortest_path_among_options.type == DubinsPathType::RSR) {
      ret(1, 0) = -1.;
      ret(2, 0) = 0.;
      ret(3, 0) = -1.;
    } else if (shortest_path_among_options.type == DubinsPathType::RSL) {
      ret(1, 0) = -1.;
      ret(2, 0) = 0.;
      ret(3, 0) = 1.;
    } else if (shortest_path_among_options.type == DubinsPathType::RLR) {
      ret(1, 0) = -1.;
      ret(2, 0) = 1.;
      ret(3, 0) = -1.;
    } else if (shortest_path_among_options.type == DubinsPathType::LRL) {
      ret(1, 0) = 1.;
      ret(2, 0) = -1.;
      ret(3, 0) = 1.;
    }
    ret(1, 1) = shortest_path_among_options.param[0]*rho;
    ret(2, 1) = shortest_path_among_options.param[1]*rho;
    ret(3, 1) = shortest_path_among_options.param[2]*rho;
  }
  return ret;
}

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

  slide_amount_where_near_turn_changes_direction = turns(turn_to_shorten, 1);
  double arclength_where_near_turn_changes_direction = 0.;
  if (turns(turn_to_shorten, 0) != 0) {
    std::unordered_set<DubinsPathType> path_type_options_zero_case;
    path_type_options_zero_case.insert(DubinsPathType::LSL);
    path_type_options_zero_case.insert(DubinsPathType::LSR);
    path_type_options_zero_case.insert(DubinsPathType::RSR);
    path_type_options_zero_case.insert(DubinsPathType::RSL);
    path_type_options_zero_case.insert(DubinsPathType::LRL);
    path_type_options_zero_case.insert(DubinsPathType::RLR);

    if (direction) {
      // Sliding forward
      if (std::abs(turns(0, 1)) > 1e-10) {
        if (turns(0, 0) != -1) {
          path_type_options_zero_case.erase(DubinsPathType::RSR);
          path_type_options_zero_case.erase(DubinsPathType::RSL);
          path_type_options_zero_case.erase(DubinsPathType::RLR);
        }

        if (turns(0, 0) != 1) {
          path_type_options_zero_case.erase(DubinsPathType::LSR);
          path_type_options_zero_case.erase(DubinsPathType::LSL);
          path_type_options_zero_case.erase(DubinsPathType::LRL);
        }
      }

      if (std::abs(turns(1, 1)) > 1e-10) {
        if (turns(1, 0) != -1) {
          path_type_options_zero_case.erase(DubinsPathType::LRL);
        }

        if (turns(1, 0) != 0) {
          path_type_options_zero_case.erase(DubinsPathType::LSL);
          path_type_options_zero_case.erase(DubinsPathType::LSR);
          path_type_options_zero_case.erase(DubinsPathType::RSR);
          path_type_options_zero_case.erase(DubinsPathType::RSL);
        }

        if (turns(1, 0) != 1) {
          path_type_options_zero_case.erase(DubinsPathType::RLR);
        }
      }

      if (std::abs(turns(3, 1)) > 1e-10) {
        if (turns(3, 0) != -1) {
          path_type_options_zero_case.erase(DubinsPathType::RSR);
          path_type_options_zero_case.erase(DubinsPathType::LSR);
          path_type_options_zero_case.erase(DubinsPathType::RLR);
        }

        if (turns(3, 0) != 1) {
          path_type_options_zero_case.erase(DubinsPathType::LSL);
          path_type_options_zero_case.erase(DubinsPathType::RSL);
          path_type_options_zero_case.erase(DubinsPathType::LRL);
        }
      }
    } else {
      // Sliding backward
      if (std::abs(turns(0, 1)) > 1e-10) {
        if (turns(0, 0) != -1) {
          path_type_options_zero_case.erase(DubinsPathType::RSR);
          path_type_options_zero_case.erase(DubinsPathType::RSL);
          path_type_options_zero_case.erase(DubinsPathType::RLR);
        }

        if (turns(0, 0) != 1) {
          path_type_options_zero_case.erase(DubinsPathType::LSR);
          path_type_options_zero_case.erase(DubinsPathType::LSL);
          path_type_options_zero_case.erase(DubinsPathType::LRL);
        }
      }

      if (std::abs(turns(2, 1)) > 1e-10) {
        if (turns(2, 0) != -1) {
          path_type_options_zero_case.erase(DubinsPathType::LRL);
        }

        if (turns(2, 0) != 0) {
          path_type_options_zero_case.erase(DubinsPathType::LSL);
          path_type_options_zero_case.erase(DubinsPathType::LSR);
          path_type_options_zero_case.erase(DubinsPathType::RSR);
          path_type_options_zero_case.erase(DubinsPathType::RSL);
        }

        if (turns(2, 0) != 1) {
          path_type_options_zero_case.erase(DubinsPathType::RLR);
        }
      }

      if (std::abs(turns(3, 1)) > 1e-10) {
        if (turns(3, 0) != -1) {
          path_type_options_zero_case.erase(DubinsPathType::RSR);
          path_type_options_zero_case.erase(DubinsPathType::LSR);
          path_type_options_zero_case.erase(DubinsPathType::RLR);
        }

        if (turns(3, 0) != 1) {
          path_type_options_zero_case.erase(DubinsPathType::LSL);
          path_type_options_zero_case.erase(DubinsPathType::RSL);
          path_type_options_zero_case.erase(DubinsPathType::LRL);
        }
      }
    }

    double q_0_for_dubins[3];
    double q_f_for_dubins[3];

    q_0_for_dubins[0] = x_0;
    q_0_for_dubins[1] = y_0;
    q_0_for_dubins[2] = theta_0;

    q_f_for_dubins[0] = x_f;
    q_f_for_dubins[1] = y_f;
    q_f_for_dubins[2] = theta_f;

    for (auto path_type : path_type_options_zero_case) {
      DubinsPath path;
      int dubins_error = dubins_path(&path, q_0_for_dubins, q_f_for_dubins, rho, path_type);
      double arclength;
      if (direction) {
        arclength = path.param[2]*rho;
      } else {
        arclength = path.param[0]*rho;
      }
      if (turns(turn_to_shorten, 1) > arclength && arclength > arclength_where_near_turn_changes_direction) {
        arclength_where_near_turn_changes_direction = arclength;
      }
    }

    // Compute the Dubins path where segment at slide_idx has length 0
    slide_amount_where_near_turn_changes_direction = turns(turn_to_shorten, 1) - arclength_where_near_turn_changes_direction;
    std::cout << "slide_amount_where_near_turn_changes_direction: " << slide_amount_where_near_turn_changes_direction << std::endl;
  }
}

RowMatrixXd ICSlidingClass::slide(double amount) {
  std::unordered_set<DubinsPathType> path_type_options;
  path_type_options.insert(DubinsPathType::LSL);
  path_type_options.insert(DubinsPathType::LSR);
  path_type_options.insert(DubinsPathType::RSR);
  path_type_options.insert(DubinsPathType::RSL);
  path_type_options.insert(DubinsPathType::LRL);
  path_type_options.insert(DubinsPathType::RLR);

  if (direction) {
    // Sliding forward
    if (std::abs(turns(0, 1)) > 1e-10) {
      if (turns(0, 0) != -1) {
        path_type_options.erase(DubinsPathType::RSR);
        path_type_options.erase(DubinsPathType::RSL);
        path_type_options.erase(DubinsPathType::RLR);
      }

      if (turns(0, 0) != 1) {
        path_type_options.erase(DubinsPathType::LSR);
        path_type_options.erase(DubinsPathType::LSL);
        path_type_options.erase(DubinsPathType::LRL);
      }
    }

    if (std::abs(turns(1, 1)) > 1e-10) {
      if (turns(1, 0) != -1) {
        path_type_options.erase(DubinsPathType::LRL);
      }

      if (turns(1, 0) != 0) {
        path_type_options.erase(DubinsPathType::LSL);
        path_type_options.erase(DubinsPathType::LSR);
        path_type_options.erase(DubinsPathType::RSR);
        path_type_options.erase(DubinsPathType::RSL);
      }

      if (turns(1, 0) != 1) {
        path_type_options.erase(DubinsPathType::RLR);
      }
    }

    if (std::abs(turns(2, 1)) > 1e-10) {
      if (amount > slide_amount_where_near_turn_changes_direction) {
        if (turns(2, 0) != 1) {
          path_type_options.erase(DubinsPathType::RSR);
          path_type_options.erase(DubinsPathType::LSR);
          path_type_options.erase(DubinsPathType::RLR);
        }

        if (turns(2, 0) != -1) {
          path_type_options.erase(DubinsPathType::LSL);
          path_type_options.erase(DubinsPathType::RSL);
          path_type_options.erase(DubinsPathType::LRL);
        }
      } else {
        if (turns(2, 0) != -1) {
          path_type_options.erase(DubinsPathType::RSR);
          path_type_options.erase(DubinsPathType::LSR);
          path_type_options.erase(DubinsPathType::RLR);
        }

        if (turns(2, 0) != 1) {
          path_type_options.erase(DubinsPathType::LSL);
          path_type_options.erase(DubinsPathType::RSL);
          path_type_options.erase(DubinsPathType::LRL);
        }
      }
    }
  } else {
    // Sliding backward
    if (std::abs(turns(1, 1)) > 1e-10) {
      if (amount > slide_amount_where_near_turn_changes_direction) {
        if (turns(1, 0) != 1) {
          path_type_options.erase(DubinsPathType::RSR);
          path_type_options.erase(DubinsPathType::RSL);
          path_type_options.erase(DubinsPathType::RLR);
        }

        if (turns(1, 0) != -1) {
          path_type_options.erase(DubinsPathType::LSR);
          path_type_options.erase(DubinsPathType::LSL);
          path_type_options.erase(DubinsPathType::LRL);
        }
      } else {
        if (turns(1, 0) != -1) {
          path_type_options.erase(DubinsPathType::RSR);
          path_type_options.erase(DubinsPathType::RSL);
          path_type_options.erase(DubinsPathType::RLR);
        }

        if (turns(1, 0) != 1) {
          path_type_options.erase(DubinsPathType::LSR);
          path_type_options.erase(DubinsPathType::LSL);
          path_type_options.erase(DubinsPathType::LRL);
        }
      }
    }

    if (std::abs(turns(2, 1)) > 1e-10) {
      if (turns(2, 0) != -1) {
        path_type_options.erase(DubinsPathType::LRL);
      }

      if (turns(2, 0) != 0) {
        path_type_options.erase(DubinsPathType::LSL);
        path_type_options.erase(DubinsPathType::LSR);
        path_type_options.erase(DubinsPathType::RSR);
        path_type_options.erase(DubinsPathType::RSL);
      }

      if (turns(2, 0) != 1) {
        path_type_options.erase(DubinsPathType::RLR);
      }
    }

    if (std::abs(turns(3, 1)) > 1e-10) {
      if (turns(3, 0) != -1) {
        path_type_options.erase(DubinsPathType::RSR);
        path_type_options.erase(DubinsPathType::LSR);
        path_type_options.erase(DubinsPathType::RLR);
      }

      if (turns(3, 0) != 1) {
        path_type_options.erase(DubinsPathType::LSL);
        path_type_options.erase(DubinsPathType::RSL);
        path_type_options.erase(DubinsPathType::LRL);
      }
    }
  }

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
  ret.col(0) = turns.col(0);
  if (amount > slide_amount_where_near_turn_changes_direction) {
    ret(slide_idx, 0) = -ret(slide_idx, 0);
  }
  ret(turn_to_shorten, 1) = dist;
  double q_0_for_dubins[3];
  double q_f_for_dubins[3];

  DubinsPath shortest_path_among_options;
  double shortest_path_length_among_options = std::numeric_limits<double>::infinity();

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
      double l = dubins_path_length(&path);
      if (l < shortest_path_length_among_options) {
        shortest_path_length_among_options = l;
        shortest_path_among_options = path;
      }
    }

    if (shortest_path_among_options.type == DubinsPathType::LSL) {
      ret(0, 0) = 1.;
      ret(1, 0) = 0.;
      ret(2, 0) = 1.;
    } else if (shortest_path_among_options.type == DubinsPathType::LSR) {
      ret(0, 0) = 1.;
      ret(1, 0) = 0.;
      ret(2, 0) = -1.;
    } else if (shortest_path_among_options.type == DubinsPathType::RSR) {
      ret(0, 0) = -1.;
      ret(1, 0) = 0.;
      ret(2, 0) = -1.;
    } else if (shortest_path_among_options.type == DubinsPathType::RSL) {
      ret(0, 0) = -1.;
      ret(1, 0) = 0.;
      ret(2, 0) = 1.;
    } else if (shortest_path_among_options.type == DubinsPathType::RLR) {
      ret(0, 0) = -1.;
      ret(1, 0) = 1.;
      ret(2, 0) = -1.;
    } else if (shortest_path_among_options.type == DubinsPathType::LRL) {
      ret(0, 0) = 1.;
      ret(1, 0) = -1.;
      ret(2, 0) = 1.;
    }
    ret(0, 1) = shortest_path_among_options.param[0]*rho;
    ret(1, 1) = shortest_path_among_options.param[1]*rho;
    ret(2, 1) = shortest_path_among_options.param[2]*rho;
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
      double l = dubins_path_length(&path);
      if (l < shortest_path_length_among_options) {
        shortest_path_length_among_options = l;
        shortest_path_among_options = path;
      }
    }

    if (shortest_path_among_options.type == DubinsPathType::LSL) {
      ret(1, 0) = 1.;
      ret(2, 0) = 0.;
      ret(3, 0) = 1.;
    } else if (shortest_path_among_options.type == DubinsPathType::LSR) {
      ret(1, 0) = 1.;
      ret(2, 0) = 0.;
      ret(3, 0) = -1.;
    } else if (shortest_path_among_options.type == DubinsPathType::RSR) {
      ret(1, 0) = -1.;
      ret(2, 0) = 0.;
      ret(3, 0) = -1.;
    } else if (shortest_path_among_options.type == DubinsPathType::RSL) {
      ret(1, 0) = -1.;
      ret(2, 0) = 0.;
      ret(3, 0) = 1.;
    } else if (shortest_path_among_options.type == DubinsPathType::RLR) {
      ret(1, 0) = -1.;
      ret(2, 0) = 1.;
      ret(3, 0) = -1.;
    } else if (shortest_path_among_options.type == DubinsPathType::LRL) {
      ret(1, 0) = 1.;
      ret(2, 0) = -1.;
      ret(3, 0) = 1.;
    }
    ret(1, 1) = shortest_path_among_options.param[0]*rho;
    ret(2, 1) = shortest_path_among_options.param[1]*rho;
    ret(3, 1) = shortest_path_among_options.param[2]*rho;
  }
  return ret;
}
