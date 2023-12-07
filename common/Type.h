#pragma once

#include <Eigen/Sparse>
#include <Eigen/Geometry>

template <int d> using MatM = Eigen::Matrix<double, 3*d-3, 3*d-3>; ///< if d=2, 3*3; if d=3, 6*6;
template <int d> using VecM = Eigen::Matrix<double, 3*d-3, 1>; ///< if d=2, 3; if d=3, 6;
template <int d> using MatE = Eigen::Matrix<double, d+1, d+1>; ///< if d=2, 3*3; if d=3, 4*4;
template <int d> using VecE = Eigen::Matrix<double, d+1, 1>; ///< if d=2, 3; if d=3, 4;
template <int d> using Veci = Eigen::Matrix<int, d, 1>;
template <int d> using Vec = Eigen::Matrix<double, d, 1>;
template <int d> using Mat = Eigen::Matrix<double, d, d>;
typedef Eigen::Matrix<double, 6, 1> Vec6d;
typedef Eigen::Matrix<double, 6, 6> Mat6d;
typedef Eigen::Matrix<double, 3, 6> Mat3x6d;
typedef Eigen::Matrix<double, 6, 3> Mat6x3d;
