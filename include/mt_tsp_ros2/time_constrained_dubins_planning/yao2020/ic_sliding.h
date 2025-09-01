#pragma once
#include "mt_tsp_ros2/elongate_dubins_path.h"
#include <unordered_set>

// Returns a sequence of (turn direction, dist) pairs.
// turns should have two or more rows.
// amount is how much to shorten the segment in the sliding direction.
// direction = true means slide forward, false means slide backward.
// Performs sliding on segment turns.rows() - 2 (zero-based indexing) if forward, and on segment 1 if backward
RowMatrixXd ic_sliding(RowMatrixXdRef_const turns, double x_0, double y_0, double theta_0, double x_f, double y_f, double theta_f, double rho, double amount, bool direction);

class ICSlidingClass {
  public:
    ICSlidingClass(RowMatrixXdRef_const turns, double x_0, double y_0, double theta_0, double x_f, double y_f, double theta_f, double rho, bool direction);

    RowMatrixXd slide(double amount);

  private:
    RowMatrixXd turns;
    double x_0;
    double y_0;
    double theta_0;
    double x_f;
    double y_f;
    double theta_f;
    double rho;
    bool direction;
    int slide_idx;
    int turn_to_shorten;

    double slide_amount_where_near_turn_changes_direction;
};
