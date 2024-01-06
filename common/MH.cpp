#include "MH.h"

namespace MH{
    template<int d> MatE<d> InverseE(const MatE<d>& E){
        MatE<d> Einv = MatE<d>::Identity();
        Mat<d> R = E.template block<d,d>(0,0);
        Vec<d> p = E.template block<d,1>(0,d);
        Mat<d> Rt = R.transpose();
        Einv.template block<d,d>(0,0) = Rt;
        Einv.template block<d,1>(0,d) = -Rt*p;
        return Einv;
    };
    template MatE<2> InverseE<2>(const MatE<2>& E);
    template MatE<3> InverseE<3>(const MatE<3>& E);

    Mat<3> Bracket3(const Vec<3>& a){
        Mat<3> A = Mat<3>::Zero();
        A(0,1) = -a(2);
        A(0,2) =  a(1);
        A(1,0) =  a(2);
        A(1,2) = -a(0);
        A(2,0) = -a(1);
        A(2,1) =  a(0);
        return A;
    }

    Vec<3> Unbracket3(const Mat<3>& A){
        Vec<3> a;
        a(0) = A(2,1);
        a(1) = A(0,2);
        a(2) = A(1,0);
        return a;
    }

    template<int d> MatM<d> Adjoint(const MatE<d>& E){
        MatM<d> ad = MatM<d>::Zero();
        if constexpr (d==2){
            ad(0,0) = 1;
            ad(1,0) = E(1,2);//p_y
            ad(2,0) = -E(0,2);//-p_x
            ad.template block<2,2>(1,1) = E.template block<2,2>(0,0);//R
        }
        else if constexpr (d==3){
            ad.template block<3,3>(0,0) = E.template block<3,3>(0,0);//R
            ad.template block<3,3>(3,3) = E.template block<3,3>(0,0);//R
            ad.template block<3,3>(3,0) = Bracket3(E.template block<3,1>(0,3))*E.template block<3,3>(0,0);//[p]R
        }
        return ad;
    };
    template MatM<2> Adjoint<2>(const MatE<2>& E);
    template MatM<3> Adjoint<3>(const MatE<3>& E);

    template<int d> MatM<d> Addot(const MatE<d>& E, const VecM<d>& phi){
        MatM<d> addot = MatM<d>::Zero();
        if constexpr (d==2){
            // 0    0
            // -TRv TRw
            Mat<2> T;
            T << 0,-1,1,0;
            Mat<2> R = E.template block<2,2>(0,0);
            addot.template block<2,1>(1,0) = -T*R*phi.template segment<2>(1);//-TRv
            addot.template block<2,2>(1,1) = T*R*phi(0);//TRw
        }
        else if constexpr (d==3){
            Mat<3> R = E.template block<3,3>(0,0);
            Vec<3> p = E.template block<3,1>(0,3);
            Vec<3> omega = phi.template segment<3>(0);
            Vec<3> v = phi.template segment<3>(3);
            addot.template block<3,3>(0,0) = R*Bracket3(omega);
            addot.template block<3,3>(3,0) = R*Bracket3(v)+Bracket3(p)*R*Bracket3(omega);
            addot.template block<3,3>(3,3) = R*Bracket3(omega);
        }
        return addot;
    }
    template MatM<2> Addot<2>(const MatE<2>& E, const VecM<2>& phi);
    template MatM<3> Addot<3>(const MatE<3>& E, const VecM<3>& phi);

    Mat6d ad(const Vec6d& phi){
        Mat6d ad = Mat6d::Zero();
        Mat<3> wbrac = Bracket3(phi.segment<3>(0));
        Mat<3> vbrac = Bracket3(phi.segment<3>(3));
        ad.block<3,3>(0,0) = wbrac;
        ad.block<3,3>(3,0) = vbrac;
        ad.block<3,3>(3,3) = wbrac;
        return ad;
    }

    Mat<3> Exp3(const Vec<3>& w){
        Mat<3> R = Mat<3>::Identity();
        double wlen = w.norm();
        if(wlen > 1e-6){// rotation angle bigger than 1e-6
            Vec<3> w1 = w/wlen;
            double wX = w1(0);
            double wY = w1(1);
            double wZ = w1(2);
            double c = std::cos(wlen);
            double s = std::sin(wlen);
            double c1 = 1-c;
            R << c + wX*wX*c1, -wZ*s + wX*wY*c1, wY*s + wX*wZ*c1,
                    wZ*s + wX*wY*c1, c + wY*wY*c1, -wX*s + wY*wZ*c1,
                    -wY*s + wX*wZ*c1, wX*s + wY*wZ*c1, c + wZ*wZ*c1;
        }
        return R;
    }
}