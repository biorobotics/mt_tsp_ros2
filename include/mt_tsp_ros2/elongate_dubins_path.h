#pragma once
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include "mt_tsp_ros2/dubins.h"
#include <stdexcept>
#include <omp.h>
#include <chrono>

using namespace Eigen;
typedef Ref<Matrix<double, Dynamic, 1>> VectorXdRef;
typedef Matrix<double, Dynamic, Dynamic, RowMajor> RowMatrixXd;

typedef const Ref<const Vector2d>& Vector2dRef_const;
typedef Matrix<bool, Dynamic, 1> VectorXb;
typedef const Ref<const RowMatrixXd>& RowMatrixXdRef_const;
typedef const Ref<const VectorXd>& VectorXdRef_const;

double root_find(double rho, double s, double l_m);

double arclength(const Vector2d &v1, const Vector2d &v2, bool left, double rho);

RowMatrixXd turns_for_dubins_path(double x_0, double y_0, double theta_0, double x_f, double y_f, double theta_f, double rho);

bool check_elongation_possible(double x_0, double y_0, double theta_0, double x_f, double y_f, double theta_f, double s, double rho, bool verbose);

bool check_elongation_possible(double x_0, double y_0, double theta_0, double x_f, double y_f, double theta_f, double s, double rho);

bool check_elongation_possible_with_profiling(double x_0, double y_0, double theta_0, double x_f, double y_f, double theta_f, double s, double rho, double &CCC_time, double &CSC_time, double &check_time);

Vector3d get_elongation_intervals(double x_0, double y_0, double theta_0, double x_f, double y_f, double theta_f, double rho);

// Returns a sequence of (turn direction, dist) pairs
RowMatrixXd elongated_dubins_path(double x_0, double y_0, double theta_0, double x_f, double y_f, double theta_f, double s, double rho, bool verbose);

VectorXb batch_elongation_check_with_profiling(RowMatrixXdRef_const q0s, RowMatrixXdRef_const qfs, VectorXdRef_const ss, double rho, int num_openmp_threads, Ref<Vector3d> profiling_info);

VectorXb batch_elongation_check(RowMatrixXdRef_const q0s, RowMatrixXdRef_const qfs, VectorXdRef_const ss, double rho, int num_openmp_threads);

VectorXd batch_dubins_path_lengths(RowMatrixXdRef_const q0s, RowMatrixXdRef_const qfs, double rho, int num_openmp_threads);

RowMatrixXd get_ccc_path(double x_0, double y_0, double theta_0, double x_f, double y_f, double theta_f, double rho, bool short_path, bool lrl);
