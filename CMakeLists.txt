cmake_minimum_required(VERSION 3.8)
project(mt_tsp_ros2)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

add_compile_options(-fopenmp)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake_modules)

# find dependencies
find_package(ament_cmake REQUIRED)
# uncomment the following section in order to fill in
# further dependencies manually.
# find_package(<dependency> REQUIRED)

find_package(ament_cmake_ros REQUIRED)
find_package(pluginlib REQUIRED)
find_package(rviz_common REQUIRED)
find_package(pybind11 REQUIRED)
find_package(OpenMP REQUIRED)
find_package(ompl REQUIRED)
find_package(CGAL REQUIRED)
# find_package(IPOPT 3.11.9 REQUIRED)
find_package(IPOPT REQUIRED)

set(IPOPT_INCLUDE_DIRS /usr/include/coin)

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
  $ENV{HOME}/Clarabel.cpp/include
  $ENV{HOME}/gurobi1103/linux64/include/
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

pybind11_add_module(gurobi_wrapper
  src/gurobi_wrapper.cpp
  src/bind_gurobi_wrapper.cpp
)
target_link_directories(gurobi_wrapper PRIVATE /$ENV{HOME}/gurobi1103/linux64/lib/)
target_link_libraries(gurobi_wrapper PRIVATE gurobi110)
install(TARGETS gurobi_wrapper
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)

pybind11_add_module(dtp_astar
  src/bind_dtp_astar.cpp
)
install(TARGETS dtp_astar
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)

pybind11_add_module(gtsp_wrapper
  src/gtsp.cpp
  src/bind_gtsp.cpp
)
target_link_directories(gtsp_wrapper PRIVATE /$ENV{HOME}/gurobi1103/linux64/lib/)
target_link_libraries(gtsp_wrapper PRIVATE gurobi110)
install(TARGETS gtsp_wrapper
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)

pybind11_add_module(gmdm_wrapper
  src/bind_gmdm.cpp
  src/dubins.cpp
)
target_include_directories(gmdm_wrapper PRIVATE /$ENV{HOME}/cppitertools/cppitertools/)
target_link_libraries(gmdm_wrapper PRIVATE OpenMP::OpenMP_CXX)
install(TARGETS gmdm_wrapper
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)

pybind11_add_module(polyanya_wrapper
  src/bind_polyanya_wrapper.cpp
  src/polyanya_wrapper.cpp
  $ENV{HOME}/polyanya/anyangle/polyanya/structs/mesh.cpp
  $ENV{HOME}/polyanya/anyangle/polyanya/helpers/geometry.cpp
  $ENV{HOME}/polyanya/anyangle/polyanya/helpers/timer.cpp
  $ENV{HOME}/polyanya/anyangle/polyanya/helpers/scenario.cpp
  $ENV{HOME}/polyanya/anyangle/polyanya/search/expansion.cpp
  $ENV{HOME}/polyanya/anyangle/polyanya/search/searchinstance.cpp
)
target_include_directories(polyanya_wrapper PRIVATE
  $ENV{HOME}/polyanya/anyangle/polyanya/structs
  $ENV{HOME}/polyanya/anyangle/polyanya/helpers
  $ENV{HOME}/polyanya/anyangle/polyanya/search
)
install(TARGETS polyanya_wrapper
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)

pybind11_add_module(implicit_grid_astar
  src/bind_implicit_grid_astar.cpp
)
install(TARGETS implicit_grid_astar
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)

pybind11_add_module(time_constrained_dubins_planning
  src/bind_time_constrained_dubins_planning.cpp
  src/dubins.cpp
)

target_include_directories(time_constrained_dubins_planning PRIVATE ${OMPL_INCLUDE_DIRS})
target_link_libraries(time_constrained_dubins_planning PRIVATE /usr/local/lib/libompl.so OpenMP::OpenMP_CXX CGAL::CGAL)
install(TARGETS time_constrained_dubins_planning
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)

pybind11_add_module(box_qp
  src/bind_box_qp.cpp
  src/box-qp.cpp
  src/box_qp_exception.cpp
)
install(TARGETS box_qp
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)

pybind11_add_module(visibility
  src/visibility_wrapper.cpp
  src/cgal_visibility_example.cpp
  src/bind_visibility.cpp
)
target_link_libraries(visibility PRIVATE CGAL::CGAL)
install(TARGETS visibility
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)

pybind11_add_module(mtvg_problem
  src/bind_min_time_segment_interception_astar_problem.cpp
)
install(TARGETS mtvg_problem
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)

pybind11_add_module(trajopt_gcs
  src/mt_ixg_star.cpp
  src/trajopt_through_convex_sets.cpp
  src/bind_trajopt_gcs.cpp
)
target_link_directories(trajopt_gcs PRIVATE $ENV{HOME}/Clarabel.cpp/rust_wrapper/target/release/)
target_link_libraries(trajopt_gcs PRIVATE clarabel_c)
install(TARGETS trajopt_gcs
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)

pybind11_add_module(parallel_squareform
  src/bind_parallel_squareform.cpp
)
target_link_libraries(parallel_squareform PRIVATE OpenMP::OpenMP_CXX)
install(TARGETS parallel_squareform
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)

pybind11_add_module(grid_2d_astar
  src/bind_grid_2d_astar.cpp
)
install(TARGETS grid_2d_astar
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)
message(WARNING "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}")

pybind11_add_module(mt_cbss
  src/bind_cbss.cpp
)
target_link_libraries(mt_cbss PRIVATE OpenMP::OpenMP_CXX)
install(TARGETS mt_cbss
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)
message(WARNING "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}")

pybind11_add_module(nearest_neighbor_edge_posterior
  src/bind_nearest_neighbor_edge_posterior.cpp
)
target_include_directories(nearest_neighbor_edge_posterior PRIVATE ${OMPL_INCLUDE_DIRS})
target_link_libraries(nearest_neighbor_edge_posterior PRIVATE /usr/local/lib/libompl.so OpenMP::OpenMP_CXX)
install(TARGETS nearest_neighbor_edge_posterior
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)
message(WARNING "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}")

pybind11_add_module(dag_dfs
  src/bind_dag_dfs.cpp
)
target_link_libraries(dag_dfs PRIVATE OpenMP::OpenMP_CXX)
install(TARGETS dag_dfs
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)
message(WARNING "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}")

pybind11_add_module(memetic_dubins_ce_mt_tsp
  src/bind_memetic_dubins_ce_mt_tsp.cpp
  src/dubins.cpp
)
target_link_libraries(memetic_dubins_ce_mt_tsp PRIVATE OpenMP::OpenMP_CXX)
install(TARGETS memetic_dubins_ce_mt_tsp
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)

pybind11_add_module(memetic_dubins_mt_tsp
  src/bind_memetic_dubins_mt_tsp.cpp
  src/dubins.cpp
)
target_link_libraries(memetic_dubins_mt_tsp PRIVATE OpenMP::OpenMP_CXX)
install(TARGETS memetic_dubins_mt_tsp
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)

pybind11_add_module(memetic_vs_dubins_mt_tsp
  src/bind_memetic_vs_dubins_mt_tsp.cpp
  src/dubins.cpp
)
target_link_libraries(memetic_vs_dubins_mt_tsp PRIVATE OpenMP::OpenMP_CXX)
install(TARGETS memetic_vs_dubins_mt_tsp
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)

pybind11_add_module(memetic_sampling
  src/bind_memetic_sampling.cpp
  src/dubins.cpp
)
target_link_libraries(memetic_sampling PRIVATE OpenMP::OpenMP_CXX)
install(TARGETS memetic_sampling
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)

pybind11_add_module(cpp_spline
  src/bind_cpp_spline.cpp
)
target_include_directories(cpp_spline PRIVATE)
install(TARGETS cpp_spline
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)

pybind11_add_module(one_sided_dubins_elongation
  src/bind_one_sided_dubins_elongation.cpp
  src/dubins.cpp
)
target_link_libraries(one_sided_dubins_elongation PRIVATE OpenMP::OpenMP_CXX)
install(TARGETS one_sided_dubins_elongation
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)

pybind11_add_module(memetic_pcg_utils
  src/bind_memetic_pcg_utils.cpp
  src/dubins.cpp
)
target_link_libraries(memetic_pcg_utils PRIVATE OpenMP::OpenMP_CXX)
install(TARGETS memetic_pcg_utils
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)

add_executable(compare_to_gdip 
  src/compare_to_gdip.cpp
  src/dubins.cpp
  $ENV{HOME}/gdip/gdip/core/src/opendubins/angleinterval.cpp
  $ENV{HOME}/gdip/gdip/core/src/opendubins/arc.cpp
  $ENV{HOME}/gdip/gdip/core/src/opendubins/circle.cpp
  $ENV{HOME}/gdip/gdip/core/src/opendubins/dip.cpp
  $ENV{HOME}/gdip/gdip/core/src/opendubins/dubins.cpp
  $ENV{HOME}/gdip/gdip/core/src/opendubins/dubinsmaneuver.cpp
  $ENV{HOME}/gdip/gdip/core/src/opendubins/dubinspath.cpp
  $ENV{HOME}/gdip/gdip/core/src/opendubins/gdip.cpp
  $ENV{HOME}/gdip/gdip/core/src/opendubins/line.cpp
  $ENV{HOME}/gdip/gdip/core/src/opendubins/path.cpp
  $ENV{HOME}/gdip/gdip/core/src/opendubins/state.cpp
)
target_include_directories(compare_to_gdip PRIVATE $ENV{HOME}/gdip/gdip/include)
target_link_libraries(compare_to_gdip PRIVATE OpenMP::OpenMP_CXX)

pybind11_add_module(adjacency_list
  src/bind_adjacency_list.cpp
)
target_link_libraries(adjacency_list PRIVATE OpenMP::OpenMP_CXX)
install(TARGETS adjacency_list
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)
message(WARNING "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}")

pybind11_add_module(efat
  src/bind_efat.cpp
)
install(TARGETS efat
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)
message(WARNING "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}")

pybind11_add_module(dubins_trj_through_seq_of_targets
  src/bind_dubins_trj_through_seq_of_targets.cpp
  src/dubins.cpp
)
target_link_libraries(dubins_trj_through_seq_of_targets PRIVATE OpenMP::OpenMP_CXX)
install(TARGETS dubins_trj_through_seq_of_targets
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)
message(WARNING "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}")

pybind11_add_module(vs_dubins_trj_through_seq_of_targets
  src/bind_vs_dubins_trj_through_seq_of_targets.cpp
  src/dubins.cpp
)
target_link_libraries(vs_dubins_trj_through_seq_of_targets PRIVATE OpenMP::OpenMP_CXX)
install(TARGETS vs_dubins_trj_through_seq_of_targets
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)
message(WARNING "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}")

pybind11_add_module(reopt_gtsp_tour
  src/bind_reopt_gtsp_tour.cpp
)
install(TARGETS reopt_gtsp_tour
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)
message(WARNING "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}")

pybind11_add_module(memetic_robot_arm
  src/bind_memetic_robot_arm.cpp
)
target_include_directories(memetic_robot_arm PRIVATE ${IPOPT_INCLUDE_DIRS})
target_link_libraries(memetic_robot_arm PRIVATE OpenMP::OpenMP_CXX ${IPOPT_LIBRARIES} ${IPOPT_LINK_FLAGS})
target_compile_definitions(memetic_robot_arm PRIVATE ${IPOPT_DEFINITIONS})
install(TARGETS memetic_robot_arm
  DESTINATION "${PYTHON_INSTALL_DIR}/${PROJECT_NAME}"
)

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
