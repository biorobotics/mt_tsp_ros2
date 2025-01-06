#include "mt_tsp_ros2/polyanya_wrapper.h"

#include "scenario.h"
#include "point.h"
#include "cfg.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>

#include <cstdio>
#include <memory>
#include <stdexcept>
#include <array>
#include <fstream>

std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), static_cast<int>(buffer.size()), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

PolyanyaWrapper::PolyanyaWrapper(const Ref<const Matrix<bool, Dynamic, Dynamic, RowMajor>> &occupancy, const std::string &polyanya_path, const std::string &map_save_path) {
  // Write grid map to file
  std::ofstream map_file(map_save_path + "/grid.map");

  map_file << "type octile\n";
  map_file << "height " << std::to_string(occupancy.rows()) << "\n";
  map_file << "width " << std::to_string(occupancy.cols()) << "\n";
  map_file << "map\n";
  for (int x_idx = 0; x_idx < occupancy.rows(); ++x_idx) {
    std::string row_str = "";
    for (int y_idx = 0; y_idx < occupancy.cols(); ++y_idx) {
      row_str = row_str + (occupancy(x_idx, y_idx) ? "T" : ".");
    }
    map_file << row_str << "\n";
  }
  map_file.close();

  // Run the shell script to convert map file to mesh
  std::ofstream mesh_file(map_save_path + "/grid.mesh");
  std::string cmd_str = polyanya_path + "/utils/bin/gridmap2rects < " + map_save_path + "/grid.map";
  mesh_file << exec(cmd_str.c_str()) << std::endl;
  mesh_file.close();

  std::ifstream mesh_file_in(map_save_path + "/grid.mesh");
  if (!mesh_file_in.is_open())
  {
    std::cerr << "Unable to open mesh" << std::endl;
    exit(1);
  }
  mesh = std::make_shared<polyanya::Mesh>(mesh_file_in);
  mesh_file_in.close();

  si = std::make_shared<polyanya::SearchInstance>(mesh.get());
}

MatrixXd PolyanyaWrapper::shortest_path(const Ref<const Vector2d>& start, const Ref<const Vector2d>& goal) {
  polyanya::Point start_point;
  start_point.x = start(0);
  start_point.y = start(1);
  polyanya::Point goal_point;
  goal_point.x = goal(0);
  goal_point.y = goal(1);
  si->set_start_goal(start_point, goal_point);
  si->search();
  std::vector<polyanya::Point> path;
  si->get_path_points(path);
  MatrixXd ret(path.size(), 2);
  for (int i = 0; i < path.size(); ++i) {
    ret(i, 0) = path[i].x;
    ret(i, 1) = path[i].y;
  }
  return ret;
}
