#pragma once
#include "mt_tsp_ros2/elongate_dubins_path.h"
#include <unordered_set>

class ICSlidingClass {
  public:
    ICSlidingClass(RowMatrixXdRef_const turns, double x_0, double y_0, double theta_0, double x_f, double y_f, double theta_f, double rho, bool final_turn, bool lengthen);

    RowMatrixXd slide(double amount, double prev_length);

    void batch_slide(std::vector<RowMatrixXd> &turns_seq, const Ref<const VectorXd> &amounts);

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
    int turn_to_shorten_or_lengthen;
    bool final_turn;
    bool lengthen;

    std::vector<DubinsPathType> path_type_options;
};
