cmake_minimum_required(VERSION 3.8)
project(mt_tsp_ros2)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# find dependencies
find_package(ament_cmake REQUIRED)
# uncomment the following section in order to fill in
# further dependencies manually.
# find_package(<dependency> REQUIRED)

find_package(ament_cmake_ros REQUIRED)
find_package(pluginlib REQUIRED)
find_package(rviz_common REQUIRED)
find_package(pybind11 REQUIRED)

# The following is for TFTrajectoryDisplay

set(CMAKE_AUTOMOC ON)
qt5_wrap_cpp(MOC_FILES
  include/mt_tsp_ros2/tf_trajectory_display.hpp
)

add_library(tf_trajectory_display src/tf_trajectory_display.cpp ${MOC_FILES})
target_include_directories(tf_trajectory_display PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
ament_target_dependencies(tf_trajectory_display
  pluginlib
  rviz_common
)
install(TARGETS tf_trajectory_display
        EXPORT export_mt_tsp_ros2
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin
)
install(DIRECTORY include/
        DESTINATION include
)
install(FILES rviz_common_plugins.xml
        DESTINATION share/${PROJECT_NAME}
)
ament_export_include_directories(include)
ament_export_targets(export_mt_tsp_ros2)
pluginlib_export_plugin_description_file(rviz_common rviz_common_plugins.xml)

# Makes sure other packages can find meshes from this package
install(DIRECTORY meshes DESTINATION share/${PROJECT_NAME})
install(DIRECTORY irb1600_aggregated DESTINATION share/${PROJECT_NAME})

# Python bindings

include_directories(
  include
)

ament_python_install_package(${PROJECT_NAME})

pybind11_add_module(astar_planner_3d_wrapper
  src/astar_planner_3d.cpp
  src/astar_planner_3d_wrapper.cpp
  src/bind_astar_planner_3d_wrapper.cpp
)
install(TARGETS astar_planner_3d_wrapper
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)
message(WARNING "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}")

# There by default
if(BUILD_TESTING)
  find_package(ament_lint_auto REQUIRED)
  # the following line skips the linter which checks for copyrights
  # comment the line when a copyright and license is added to all source files
  set(ament_cmake_copyright_FOUND TRUE)
  # the following line skips cpplint (only works in a git repo)
  # comment the line when this package is in a git repo and when
  # a copyright and license is added to all source files
  set(ament_cmake_cpplint_FOUND TRUE)
  ament_lint_auto_find_test_dependencies()
endif()

ament_package()
