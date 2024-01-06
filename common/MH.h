#pragma once

#include "Type.h"

#include <Eigen/Sparse>
#include <Eigen/Geometry>

namespace MH{///< Math Helper Functions
    template<int d> MatE<d> InverseE(const MatE<d>& E); ///< inverse an affine transformation matrix
    Mat<3> Bracket3(const Vec<3>& a); ///< bracket3(a)*b = a.cross(b)
    Vec<3> Unbracket3(const Mat<3>& A);
    template<int d> MatM<d> Adjoint(const MatE<d>& E); ///< adjoint transformation of E
    template<int d> MatM<d> Addot(const MatE<d>& E, const VecM<d>& phi);
    Mat6d ad(const Vec6d& phi); ///< used to compute Addot and Coriolis force
    Mat<3> Exp3(const Vec<3>& w); ///< exp representation of rotation
};