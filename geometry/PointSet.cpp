//
// Created by zhangjie on 2024/1/6.
//

#include "PointSet.h"

template<int d>
void PointSet<d>::Update_Nbs()
{
    //update data structure nbs_searcher
    std::function<bool(const int)> valid = [&](const int idx)->bool {return points->I(idx) != -1; };
    nbs_searcher->Update_Points(points->XRef(), valid);
    //use updated nbs_searcher to update tang_nbs
    int pn = points->Size();
    tang_nbs.resize(pn);
#pragma omp parallel for
    for (int i = 0; i < pn; i++) {
        if (!Valid_Particle(i)) continue;
        const VectorD &pos = points->X(i);
        const MatrixD& local_frame = points->E(i);
        std::function<bool(const int)> filter_func = [&](const int p) {return Is_Tangential_Neighbor(pos, local_frame, p); };
        nbs_searcher->Find_Neighbors(pos, v_r, filter_func, tang_nbs[i]);
    }
}

template<int d>
bool PointSet<d>::Is_Tangential_Neighbor(const VectorD& pos, const MatrixD& local_frame, const int p) const
{
    ////check angle
    VectorD dir = local_frame.col(d - 1);
    VectorD dir_p = points->E(p).col(d - 1);
    real dot = dir.dot(dir_p);
    if (dot < t_dot) return false;	////skip the points with large angles

    ////check distance
    VectorD u = points->X(p) - pos;
    VectorT t; Project_To_TPlane(u, local_frame, t);
    return t.norm() < t_r;
}

template<int d>
ArraySlice<int> PointSet<d>::Record_Tangential_Nbs(const VectorD& pos, const MatrixD& local_frame) const
{
    std::function<bool(const int)> filter_func = [&](const int p) {return Is_Tangential_Neighbor(pos, local_frame, p); };
    return nbs_searcher->Record_Neighbors(pos, v_r, filter_func);
}

template<int d>
Array<int> PointSet<d>::Find_Tangential_Nbs(const VectorD& pos, const MatrixD& local_frame) const
{
    std::function<bool(const int)> filter_func = [&](const int p) {return Is_Tangential_Neighbor(pos, local_frame, p); };
    return nbs_searcher->Find_Neighbors(pos, v_r, filter_func);
}

template<int d>
void PointSet<d>::Update_Local_Frame(const real dt)
{
    int p_n = points->Size();
    for (int i = 0; i < p_n; i++) {
        if (!Valid_Particle(i))continue;

        if constexpr (d == 2) {
            VectorD e0 = points->E(i).col(0);
            VectorD e1 = points->E(i).col(1);

            ////return u \dot e1 in tangential space
            std::function<real(const int)> V_Normal_Val = [&](const int idx)->real
            {return points->V(idx).dot(e1); };

            VectorT dvdx = Grad_TPlane(i, V_Normal_Val);
            std::cout << dvdx[0] << ", ";
            e0 = Rotation2(dvdx[0] * dt) * e0;
            e1 = Rotation2(dvdx[0] * dt) * e1;
        }
        else if constexpr (d == 3) {
            ////TOIMPROVE: use rotation matrix or quaternion
            VectorD e0 = points->E(i).col(0);
            VectorD e1 = points->E(i).col(1);
            VectorD e2 = points->E(i).col(2);

            std::function<real(const int)> U0_Val = [&](const int idx)->real
            {return points->V(idx).dot(points->E(i).col(0)); };
            VectorT du0 = Grad_TPlane(i, U0_Val);

            std::function<real(const int)> U1_Val = [&](const int idx)->real
            {return points->V(idx).dot(points->E(i).col(1)); };
            VectorT du1 = Grad_TPlane(i, U1_Val);

            std::function<real(const int)> U2_Val = [&](const int idx)->real
            {return points->V(idx).dot(points->E(i).col(2)); };
            VectorT du2 = Grad_TPlane(i, U2_Val);

            VectorD de0 = (real).5 * (du1[0] - du0[1]) * e1 + (du2[0]) * e2;
            VectorD de1 = (real).5 * (du0[1] - du1[0]) * e0 + (du2[1]) * e2;
            VectorD de2 = (-du2[0]) * e0 + (-du2[1]) * e1;

            e0 += de0 * dt; e0.normalize();
            e1 += de1 * dt; e1.normalize();
            e2 += de2 * dt; e2.normalize();

            points->E(i).col(0) = e0;
            points->E(i).col(1) = e1;
            points->E(i).col(2) = e2;
        }
    }
}

template<int d>
void PointSet<d>::Initialize_Local_Frame(const VectorD& v, const Array<int>& nbs,MatrixD& local_frame){

    ////calculate the PCA normal
    VectorD xp = VectorD::Zero();
    real w = (real)0;
    for (size_t i = 0; i < nbs.size(); i++) {
        ////here we use volumetric distance instead of tangential distance
        real dis = (points->X(i) - v).norm();
        real w0 = W_PCA(dis);
        xp += w0 * points->X(i);
        w += w0;
    }
    if (w != (real)0)xp /= w;
    MatrixD C = MatrixD::Zero();
    real wc = (real)0;
    for (size_t i = 0; i < nbs.size(); i++) {
        const VectorD& xi = points->X(i);
        real dis = (xi - xp).norm();
        real w0 = W_PCA(dis);
        C += w0 * (xi - xp) * (xi - xp).transpose();
        wc += w0;
    }
    if (wc != (real)0)C /= wc;
    VectorD normal = AuxFunc::Min_Eigenvector(C);
    // if (normal.dot(points->Normal(p)) < (real)0)normal *= (real)-1;
    if (normal.dot(Normal(v)) < (real)0)normal *= (real)-1; ////TOTEST

    ////update local frame according to the PCA normal
    if constexpr (d == 2) {
        VectorD tang = -AuxFunc::Orthogonal_Vector(normal);
        local_frame.col(0) = tang.normalized();
        local_frame.col(1) = normal.normalized();
    }
    else if constexpr (d == 3) {
        VectorD t1 = -AuxFunc::Orthogonal_Vector(normal);
        VectorD t2 = t1.cross(normal);
        local_frame.col(0) = t1.normalized();
        local_frame.col(1) = t2.normalized();
        local_frame.col(2) = normal.normalized();
    }
}

template<int d>
void PointSet<d>::Reinitialize_Local_Frames()
{
    ////update local normals with PCA
    int pn = points->Size();
#pragma omp parallel for
    for (int p = 0; p < pn; p++) {
        if (!Valid_Particle(p))continue;

        ////calculate the PCA normal
        size_t nbs_num = tang_nbs[p].size();
        VectorD xp = VectorD::Zero();
        real w = (real)0;
        for (size_t i = 0; i < nbs_num; i++) {
            int q = tang_nbs[p][i];
            ////here we use volumetric distance instead of tangential distance
            real dis = (points->X(q) - points->X(p)).norm();
            real w0 = W_PCA(dis);
            xp += w0 * points->X(q);
            w += w0;
        }
        if (w != (real)0)xp /= w;
        MatrixD C = MatrixD::Zero();
        real wc = (real)0;
        for (int i = 0; i < nbs_num; i++) {
            int q = tang_nbs[p][i];
            const VectorD& xq = points->X(q);
            real dis = (xq - xp).norm();
            real w0 = W_PCA(dis);
            C += w0 * (xq - xp) * (xq - xp).transpose();
            wc += w0;
        }
        if (wc != (real)0)C /= wc;
        VectorD normal = AuxFunc::Min_Eigenvector(C);
        if (normal.dot(points->Normal(p)) < (real)0)normal *= (real)-1;

        ////update local frame according to the PCA normal
        if constexpr (d == 2) {
            VectorD tang = -AuxFunc::Orthogonal_Vector(normal);
            points->E(p).col(0) = tang.normalized();
            points->E(p).col(1) = normal.normalized();
        }
        else if constexpr (d == 3) {
            VectorD t1 = -AuxFunc::Orthogonal_Vector(normal);
            VectorD t2 = t1.cross(normal);
            points->E(p).col(0) = t1.normalized();
            points->E(p).col(1) = t2.normalized();
            points->E(p).col(2) = normal.normalized();
        }
    }

    //#pragma omp parallel for
    //for (int p = 0; p < pn; p++) {
    //
    //}

    ////Correct local frame to align with the tangent space
    const bool use_local_frame_correction = true;
    if(use_local_frame_correction) {////
#pragma omp parallel for
        for (int p = 0; p < pn; p++) {
            if (!Valid_Particle(p))continue;
            ////return h in tangential space
            std::function<real(const int)> H_Val = [=](const int idx)->real {
                VectorD u = points->X(idx) - points->X(p); real h = Project_To_TPlane_H(u, points->E(p)); return h; };

            VectorT dzdx = Grad_TPlane(p, H_Val);

            if constexpr (d == 2) {
                VectorD tang = (points->E(p) * VectorD(1, dzdx(0))).normalized();
                VectorD normal = -AuxFunc::Orthogonal_Vector(tang);

                points->E(p).col(0) = tang.normalized();
                points->E(p).col(1) = normal.normalized();
            }
            else if constexpr (d == 3) {
                VectorD t0 = points->E(p) * VectorD(1, 0, dzdx(0));
                VectorD t1 = points->E(p) * VectorD(0, 1, dzdx(1));
                VectorD normal = t1.cross(t0).normalized();

                t0 = -AuxFunc::Orthogonal_Vector(normal);
                t1 = t0.cross(normal);
                points->E(p).col(0) = t0.normalized();
                points->E(p).col(1) = t1.normalized();
                points->E(p).col(2) = normal.normalized();
            }
        }
    }
}

template<int d>
Vector<real,d> PointSet<d>::Normal(const VectorD& pos) const
{
    int closest_p = Closest_Point(pos);
    Array<int> nbs = Find_Tangential_Nbs(pos, points->E(closest_p));
    size_t nb_n = nbs.size();
    if (nb_n == 0)return VectorD::Zero();

    VectorD nml = VectorD::Zero();
    for (int i = 0; i < nb_n; i++) {
        int p = nbs[i];
        real dis = (pos - points->X(p)).norm();
        real w0 = W_PCA(dis);
        nml += w0 * Normal(p);
    }
    return nml.normalized();
}

template<int d>
Matrix<real, d> PointSet<d>::Local_Frame(const VectorD& pos) const
{
    VectorD normal = Normal(pos);
    MatrixD e;
    if constexpr (d == 2) {
        VectorD tang = -AuxFunc::Orthogonal_Vector(normal);
        e.col(0) = tang.normalized();
        e.col(1) = normal.normalized();
    }
    else if constexpr (d == 3) {
        VectorD t1 = -AuxFunc::Orthogonal_Vector(normal);
        VectorD t2 = t1.cross(normal);
        e.col(0) = t1.normalized();
        e.col(1) = t2.normalized();
        e.col(2) = normal.normalized();
    }
    return e;
}

template<int d>
void PointSet<d>::Update_Metric_Tensor()
{
    int p_n = points->Size();
    for (int i = 0; i < p_n; i++) {
        if (!Valid_Particle(i))continue;

        ////return h in tangential space
        std::function<real(const int)> H_Val = [=](const int idx)->real
        {VectorD u = points->X(idx) - points->X(i); real h = Project_To_TPlane_H(u, points->E(i)); return h; };

        VectorT dzdx = Grad_TPlane(i, H_Val);
        points->dH(i) = dzdx;
        points->G(i) = Metric_Tensor(dzdx);
    }
}

template<int d>
Matrix<real, 1> PointSet<d>::Metric_Tensor(const Vector1& dzdx) const
{
    Matrix<real, 1> mt; mt << 1 + pow(dzdx[0], 2);
    return mt;
}

template<int d>
Matrix<real, 2> PointSet<d>::Metric_Tensor(const Vector2& dzdx) const
{
    Matrix<real, 2> mt;
    real g11 = 1 + pow(dzdx[0], 2);
    real g22 = 1 + pow(dzdx[1], 2);
    real g12 = dzdx[0] * dzdx[1];
    mt << g11, g12, g12, g22;
    return mt;
}

template<int d>
Vector<real, d> PointSet<d>::Project_To_Surface(const VectorD& pos) const
{
    ////The current implementation only conducts one iteration. The function can be called for multiple times for an iterative projection.
    using namespace LeastSquares;
    MatrixD local_frame = Local_Frame(pos);
    Array<int> nbs = Find_Tangential_Nbs(pos, local_frame);
    LS<d - 1, 2> ls;
    size_t n = nbs.size(); Array<real> data(n * (size_t)d);
    for (size_t i = 0; i < n; i++) {
        int p = nbs[i];
        VectorD u = points->X(p) - pos;
        VectorD th; Project_To_TPlane(u, local_frame, th);
        for (size_t j = 0; j < d; j++)data[i * (size_t)d + j] = th[j];
    }
    ls.Fit(&data[0], (int)n);
    VectorD proj_pos = pos + local_frame.col(d - 1) * ls(VectorT::Zero());
    return proj_pos;
}

template<int d>
real PointSet<d>::Unsigned_Distance(const VectorD& pos) const
{
    VectorD proj_pos = Project_To_Surface(pos);
    VectorD v = proj_pos - pos;
    return v.norm();
}

template<int d>
real PointSet<d>::Signed_Distance(const VectorD& pos) const
{
    VectorD proj_pos = Project_To_Surface(pos);
    VectorD v = proj_pos - pos;
    VectorD normal = Normal(proj_pos);
    real sign = v.normalized().dot(normal) > (real)0 ? (real)1 : (real)-1;
    return sign * v.norm();
}

template<int d>
Vector<real, d - 1> PointSet<d>::Local_Dzdx(const int i, const int j)
{
    if constexpr (d == 2) {
        const VectorD ti = points->E(i).col(0);
        const VectorD ni = points->E(i).col(1);
        const VectorD nj = points->E(j).col(2);
        const real nj_ti = nj.dot(ti);
        real nj_ni = nj.dot(ni);

        if (abs(nj_ni) < 1e-8) {
            if (nj_ni < 0) nj_ni = 1e-8;
            else nj_ni = -1e-8;
        }
        return VectorT(nj_ti / nj_ni);
    }
    else if constexpr (d == 3) {
        const VectorD ti0 = points->E(i).col(0);
        const VectorD ti1 = points->E(i).col(1);
        const VectorD ni = points->E(i).col(2);
        const VectorD nj = points->E(j).col(2);
        const real nj_ti0 = nj.dot(ti0);
        const real nj_ti1 = nj.dot(ti1);
        real nj_ni = nj.dot(ni);

        if (abs(nj_ni) < 1e-8) {
            if (nj_ni < 0) nj_ni = 1e-8;
            else nj_ni = -1e-8;
        }
        return VectorT(nj_ti0 / nj_ni, nj_ti1 / nj_ni);
    }
}

template<int d>
Vector<real, d - 1> PointSet<d>::Local_Dzdx(const MatrixD& lf, const int j)
{
    if constexpr (d == 2) {
        const VectorD ti = lf.col(0);
        const VectorD ni = lf.col(1);
        const VectorD nj = points->E(j).col(2);
        const real nj_ti = nj.dot(ti);
        real nj_ni = nj.dot(ni);

        if (abs(nj_ni) < 1e-8) {
            if (nj_ni < 0) nj_ni = 1e-8;
            else nj_ni = -1e-8;}
        return VectorT(nj_ti / nj_ni);
    }
    else if constexpr (d == 3) {
        const VectorD ti0 = lf.col(0);
        const VectorD ti1 = lf.col(1);
        const VectorD ni = lf.col(2);
        const VectorD nj = points->E(j).col(2);
        const real nj_ti0 = nj.dot(ti0);
        const real nj_ti1 = nj.dot(ti1);
        real nj_ni = nj.dot(ni);

        if (abs(nj_ni) < 1e-8) {
            if (nj_ni < 0) nj_ni = 1e-8;
            else nj_ni = -1e-8;}
        return VectorT(nj_ti0 / nj_ni, nj_ti1 / nj_ni);
    }
}

template<int d>
real PointSet<d>::Calculate_Number_Density(const Array<VectorD>& X, Array<real>& nden)
{
    int pn = points->Size(); real avg = (real)0; int n = 0; nden.resize(pn); //n is the number of active particles
#pragma omp parallel for
    for (int i = 0; i < pn; i++) {
        if (!Valid_Particle(i))continue;
        real nd = (real)0;
        size_t nb_n = tang_nbs[i].size();
        for (size_t k = 0; k < nb_n; k++) {
            int j = tang_nbs[i][k];
            VectorT lr_ij; Project_To_TPlane(X[i] - X[j], points->E(i), lr_ij);
            nd += t_kernel.W_Poly6(lr_ij.norm());
        }
        nden[i] = nd; avg += nd; n++;
    }
    if (n > 0)avg /= (real)n; return avg;
}

template<int d>
real PointSet<d>::Calculate_Number_Density(const Array<VectorD>& X, const VectorD& pos) const
{
    MatrixD local_frame = Local_Frame(pos);
    Array<int> nbs_arr = Find_Tangential_Nbs(pos, local_frame);
    real nd = (real)0; size_t nb_n = nbs_arr.size();
    for (size_t k = 0; k < nb_n; k++) {
        size_t j = nbs_arr[k];
        VectorT lr_ij; Project_To_TPlane(pos - X[j], local_frame, lr_ij);
        nd += t_kernel.W_Poly6(lr_ij.norm());
    }
    return nd;
}

template<int d>
void PointSet<d>::Point_Reseeding()
{
    real reseeding_nden = (real)0.6 * init_avg_nden;
    real deleting_nden = (real)2 * init_avg_nden;
    int pn = points->Size();

    ////find all reseeding pairs
    Array<Vector2i> reseeding_pairs;
    //Array<int> deleting_idx;
    for (int i = 0; i < pn; i++) {
        if (!Valid_Particle(i))continue;
        size_t nb_n = tang_nbs[i].size();
        for (size_t k = 0; k < nb_n; k++) {
            int j = tang_nbs[i][k];
            if (j <= i)continue;//avoid counting itself
            VectorD mid_pos = (real).5 * (points->X(i) + points->X(j));
            real mid_nden = Calculate_Number_Density(points->XRef(), mid_pos);
            if (mid_nden < reseeding_nden) {
                reseeding_pairs.push_back(Vector2i(i, j));
            }
        }
        real i_nden = Calculate_Number_Density(points->XRef(), points->X(i));
        //if(i_nden>deleting_nden)
        //	deleting_idx.push_back(i);
    }

    ////reseed points according to the selected pairs
    reseeded_points.clear();
    int new_size = pn + (int)reseeding_pairs.size();
    points->Resize(new_size);
    for (int k = 0; k < reseeding_pairs.size(); k++) {
        int i = reseeding_pairs[k][0]; int j = reseeding_pairs[k][1];
        //int p=Add_Particle(); ////some bug with Add_Particles()
        int p = pn + k;
        reseeded_points.push_back(p);
        points->X(p) = (real).5 * (points->X(i) + points->X(j));
        points->V(p) = (real).5 * (points->V(i) + points->V(j));
        points->E(p) = Local_Frame(points->X(p));
        points->M(p) = (real)1;
        points->I(p) = 0;
    }
    if (reseeding_pairs.size() > 0) Update();
    ////deleting points with high density
    //for (int k=0; k<deleting_idx.size(); k++)
    //{Remove_Particle(k);
    //std::cout<<"particle removed"<<std::endl;
    //}
}

template<int d>
void PointSet<d>::Point_Relaxation()
{
    int pn = points->Size();
    Array<real> nden(pn);
    Array<VectorT> local_f(pn);
    Array<VectorD> relaxed_pos(pn);

    real kp = (real)1e3;
    real nden_0 = (real).8 * init_avg_nden;	////set the nden_0 to be smaller than the initial average to force points to push each other
    real one_over_m = (real)1;
    real vis = (real)1;
    real dt = (real).02;
    int relax_iter_num = 4;

    auto P = [&](const int idx)->real {return kp * (nden[idx] / nden_0 - (real)1.); };	////lambda function to access P
    auto Vol = [&](const int idx)->real {return (real)1 / nden[idx]; };
    relaxed_pos = points->XRef();	////initialize relaxed_pos by *copying* data from points

    for (int iter = 0; iter < relax_iter_num; iter++) {
        ////update number density
        for (int i = 0; i < pn; i++) {
            if (points->I(i) == -1)continue;
            real nd = (real)0;
            int nb_n = (int)tang_nbs[i].size();
            for (int k = 0; k < nb_n; k++) {
                int j = tang_nbs[i][k];
                VectorT lr_ij; Project_To_TPlane(relaxed_pos[i] - relaxed_pos[j], points->E(i), lr_ij);
                nd += t_kernel.W_Poly6(lr_ij.norm());
            }
            nden[i] = nd;
        }

        ////update forces
        for (int i = 0; i < pn; i++) {
            if (points->I(i) == -1)continue;
            VectorT lf = VectorT::Zero();
            size_t nb_n = tang_nbs[i].size();
            for (size_t k = 0; k < nb_n; k++) {
                int j = tang_nbs[i][k];
                VectorT lr_ij; Project_To_TPlane(relaxed_pos[i] - relaxed_pos[j], points->E(i), lr_ij);
                real lr2 = lr_ij.squaredNorm();
                VectorT lf_p = -(P(i) * pow(Vol(i), 2) + P(j) * pow(Vol(j), 2)) * t_kernel.Grad_Spiky(lr_ij);
                lf += lf_p;
            }
            local_f[i] = one_over_m * lf;
        }

        ////time integration
        for (int i = 0; i < pn; i++) {
            if (points->I(i) == -1)continue;
            VectorD delta_x; Unproject_To_World(VectorT(local_f[i] * dt * dt), points->E(i), delta_x);
            relaxed_pos[i] += delta_x;
            relaxed_pos[i] = Project_To_Surface(relaxed_pos[i]);
        }
    }

    ////update positions
    for (int i = 0; i < pn; i++) {
        if (points->I(i) == -1)continue;
        points->X(i) = relaxed_pos[i];
    }
}

template<int d>
void PointSet<d>::Print_Statistics()
{
    using namespace AuxFunc;
    Seperation_Line(2);
    int avg_nb_num = 0;
    int min_nb_num = std::numeric_limits<int>::max();
    int max_nb_num = 0;
    int n = 0;
    for (int i = 0; i < points->Size(); i++) {
        if (points->I(i) == -1)continue;
        int i_nb_num = (int)tang_nbs[i].size();
        avg_nb_num += i_nb_num; n++;
        if (i_nb_num < min_nb_num) min_nb_num = i_nb_num;
        if (i_nb_num > max_nb_num) max_nb_num = i_nb_num;
    }
    if (n != 0)avg_nb_num /= n;
    std::cout << "[#Point] active: " << n << ", total: " << points->Size() << std::endl;
    std::cout << "[#Nb] avg: " << avg_nb_num << ", min: " << min_nb_num << ", max: " << max_nb_num << std::endl;
    std::cout << "[Nden] avg: " << avg_nden << ", init: " << init_avg_nden << std::endl;
    Seperation_Line(2);
}

template class PointSet<2>;
template class PointSet<3>;