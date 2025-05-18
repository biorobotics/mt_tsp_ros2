#include "mt_tsp_ros2/elongate_one_sided_dubins_path.h"
#include "opendubins/dubins.h"
#include <random>

using namespace opendubins;

int main() {
  std::mt19937 rng;
  std::uniform_real_distribution<double> x_dist(-50, 50);
  std::uniform_real_distribution<double> y_dist(-50, 50);
  std::uniform_real_distribution<double> theta_dist(0, 2*M_PI);

  double rho = 12.;

  for (int i = 0; i < 1000; ++i) {
    double x_0 = x_dist(rng);
    double y_0 = y_dist(rng);
    double theta_0 = y_dist(rng);
    double x_f = x_dist(rng);
    double y_f = y_dist(rng);
    RowMatrixXd turns = turns_for_one_sided_dubins_path(x_0, y_0, theta_0, x_f, y_f, rho);
    Dubins path(AngleInterval(Point(x_0, y_0), theta_0, 0), AngleInterval(Point(x_f, y_f), 0, 2*M_PI), rho);
    if ((turns(0, 0) < 0 && turns(1, 0) > 0 && path.getType() != DType::DIP_RLp) || 
        (turns(0, 0) > 0 && turns(1, 0) < 0 && path.getType() != DType::DIP_LRp) || 
        (turns(0, 0) < 0 && turns(1, 0) == 0 && path.getType() != DType::DIP_RS) || 
        (turns(0, 0) > 0 && turns(1, 0) == 0 && path.getType() != DType::DIP_LS)
       ) {
      throw std::runtime_error("Path type incorrect");
    }
    if (std::abs(turns(0, 1) - std::abs(rho*path.getLen1())) > 1e-4) {
      throw std::runtime_error("Large difference from GDIP on segment 1");
    }
    if ((turns(1, 0) != 0 && std::abs(turns(1, 1) - std::abs(rho*path.getLen3())) > 1e-4) || 
        (turns(1, 0) == 0 && std::abs(turns(1, 1) - path.getLen2()) > 1e-4)) {
      throw std::runtime_error("Large difference from GDIP on segment 2");
    }
  }

  return 0;
}
