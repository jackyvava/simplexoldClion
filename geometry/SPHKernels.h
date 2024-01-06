//
// Created by zhangjie on 2024/1/6.
//

#ifndef SIMPLEXOLDCLION_SPHKERNELS_H
#define SIMPLEXOLDCLION_SPHKERNELS_H

#include "Common.h"
#include "Constants.h"
#include <iostream>

template<int d> class KernelSPH
{
    Typedef_VectorDii(d);
public:

    const real gaussian_trunc = 3;//h=3*sigma, for gaussian kernel

    real h;
    real h2;
    real h3;
    real hd1;
    real h13;
    real h23;

    real coef_poly6;
    real coef_poly6_grad;
    real coef_poly6_lap;

    real coef_spiky;
    real coef_spiky_grad;
    real coef_spiky_lap;

    real coef_vis;
    real coef_vis_grad_1;
    real coef_vis_grad_2;
    real coef_vis_grad_3;
    real coef_vis_lap;

    real coef_cubic_spline;
    real coef_quintic;

    real coef_bell;
    real coef_bell_grad;

    real coef_gaussian;
    real coef_gaussian_grad;

    real coef_quadratic;

    real coef_new_quadratic;

    KernelSPH() { Precompute_Coefs(1.0); }
    KernelSPH(const real _h) { Precompute_Coefs(_h); }
    virtual void Initialize(const real _h) { Precompute_Coefs(_h); }
    virtual void Precompute_Coefs(const real _h)
    {
        h = _h;
        h2 = h * h;
        h3 = h * h * h;
        hd1 = (real)1 / h;
        h13 = h * one_third;
        h23 = h * two_thirds;

        if constexpr (d == 3) {
            ////poly6
            coef_poly6 = 315.0 / (64.0 * pi * pow(h, 9));
            coef_poly6_grad = -945.0 / (32.0 * pi * pow(h, 9));
            coef_poly6_lap = coef_poly6_grad;
            ////spiky
            coef_spiky = 15.0 / (pi * pow(h, 6));
            coef_spiky_grad = -45.0 / (pi * pow(h, 6));
            coef_spiky_lap = -90.0 / (pi * pow(h, 6));
            ////vis
            coef_vis = 15.0 / (2.0 * pi * h3);
            coef_vis_grad_1 = coef_vis;
            coef_vis_grad_2 = -1.5 / h3;
            coef_vis_grad_3 = 2.0 / h2;
            coef_vis_lap = 45.0 / (pi * h3);
            ////cubic
            // coef_cubic_spline=3.0/(2.0*pi*h3);
            coef_cubic_spline = 8.0 / (pi * h3);
            ////interpolation4
            coef_quintic = 81.0 / (359.0 * pi * h3);
            ////bell
            coef_bell = 105.0 / (16.0 * pi * h3);
            coef_bell_grad = coef_bell;
            ////Gaussian kernel
            coef_gaussian = pow(gaussian_trunc, 3) / (pow(pi, 1.5) * h3);
            coef_gaussian_grad = coef_gaussian * (-2) * pow(gaussian_trunc, 2) / pow(h, 2);
            ////Quadratic kernel
            coef_quadratic = 5.0 / (4.0 * pi * h3);
            ////New quadratic kernel
            coef_new_quadratic = 315.0 / (208 * pi * h3);
        }
        else if constexpr (d == 2) {
            ////poly6
            coef_poly6 = 4.0 / (pi * pow(h, 8));
            coef_poly6_grad = -24.0 / (pi * pow(h, 8));
            coef_poly6_lap = coef_poly6_grad;
            ////spiky
            coef_spiky = 10.0 / (pi * pow(h, 5));
            coef_spiky_grad = -30.0 / (pi * pow(h, 5));
            coef_spiky_lap = -60.0 / (pi * pow(h, 5));
            ////vis
            coef_vis = 10.0 / (3.0 * pi * h2);
            coef_vis_grad_1 = coef_vis;
            coef_vis_grad_2 = -1.5 / h3;
            coef_vis_grad_3 = 2.0 / h2;
            coef_vis_lap = 20.0 / (pi * h2);
            ////cubic
            // coef_cubic_spline=15.0/(7.0*pi*h2);
            coef_cubic_spline = 40.0 / (7.0 * pi * h2);
            ////interpolation4
            coef_quintic = 63.0 / (478.0 * pi * h2);
            ////bell
            coef_bell = 5.0 / (pi * h2);
            coef_bell_grad = coef_bell;
            ////Gaussian kernel
            coef_gaussian = pow(gaussian_trunc, 2) / (pi * h2);
            coef_gaussian_grad = coef_gaussian * (-2) * pow(gaussian_trunc, 2) / pow(h, 2);
            ////Quadratic kernel
            coef_quadratic = 2.0 / (pi * h2);
            ////New quadratic kernel
            coef_new_quadratic = 15.0 / (7.0 * pi * h2);
        }
        else if constexpr (d == 1) {
            ////poly6
            coef_poly6 = 35.0 / (32.0 * pow(h, 7));
            coef_poly6_grad = -105.0 / (16.0 * pow(h, 7));
            coef_poly6_lap = coef_poly6_grad;
            ////spiky
            coef_spiky = 2.0 / pow(h, 4);
            coef_spiky_grad = -6.0 / pow(h, 4);
            coef_spiky_lap = -12.0 / pow(h, 4);
            ////vis
            coef_vis = (real)0;
            coef_vis_grad_1 = coef_vis;
            coef_vis_grad_2 = -1.5 / h3;
            coef_vis_grad_3 = 2.0 / h2;
            coef_vis_lap = (real)0;
            ////cubic
            // coef_cubic_spline=1.0/h;
            coef_cubic_spline = 4.0 / (3.0 * h);
            ////interpolation4
            coef_quintic = 1.0 / (40 * h);
            ////bell
            coef_bell = 5.0 / (4.0 * h);
            coef_bell_grad = coef_bell;
            ////Gaussian kernel
            coef_gaussian = gaussian_trunc / (sqrt(pi) * h);
            coef_gaussian_grad = coef_gaussian * (-2) * pow(gaussian_trunc, 2) / pow(h, 2);
            ////Quadratic kernel
            coef_quadratic = hd1;
            ////New quadratic kernel
            coef_new_quadratic = hd1;
        }
    }

    //////////////////////////////////////////////////////////////////////////
    ////kernels

    inline real W_Poly6(const real r) const
    {
        if (r < h) return coef_poly6 * pow(h2 - r * r, 3);
        else return (real)0;
    }

    inline real W_Spiky(const real r) const
    {
        if (r < h) return coef_spiky * pow(h - r, 3);
        else return (real)0;
    }

    inline real W_Cubic(const real r) const
    {
        //Truncate to [0,h). So h here is 2h in mathematical formula
        real u = r / h;
        // if(u>=0&&u<1) return coef_cubic_spline*(two_thirds-pow(u,2)+(real).5*pow(u,3));
        // else if(u>=1&&u<2) return coef_cubic_spline*(one_sixth*pow((real)2.-u,3));
        if (u >= 0 && u <= 0.5) return coef_cubic_spline * ((real)6 * pow(u, 3) - (real)6 * pow(u, 2) + 1);
            //else if (u > 0.5 && u < 1) return coef_cubic_spline * (2 * (1 - pow((real)1. - u, 3))); --- wrong
        else if (u > 0.5 && u < 1) return coef_cubic_spline * (2 * pow(1 - u, 3));
        else return (real)0;
    }

    inline real W_Quintic(const real r) const
    {
        if (r >= 0 && r < h13) {
            return coef_quintic * (pow(3 - 3 * r * hd1, 5) - 6 * pow(2 - 3 * r * hd1, 5) + 15 * pow(1 - 3 * r * hd1, 5));
        }
        else if (r >= h13 && r < h23) {
            return coef_quintic * (pow(3 - 3 * r * hd1, 5) - 6 * pow(2 - 3 * r * hd1, 5));
        }
        else if (r >= h23 && r < h) {
            return coef_quintic * (pow(3 - 3 * r * hd1, 5));
        }
        else { return 0; }
    }

    inline real W_Bell(const real r) const
    {
        if (r <= h) return coef_bell * (1.0 + 3.0 * r * hd1) * pow(1.0 - r * hd1, 3);
        else return (real)0;
    }

    inline real W_Gaussian(const real r) const
    {
        return coef_gaussian * exp(-pow(gaussian_trunc * r / h, 2));
    }

    inline real W_Quadratic(const real r) const
    {
        if (r <= 2 * h) return coef_quadratic * (3.0 / 16.0 * r * r / h2 - 0.75 * r * hd1 + 0.75);
        else return (real)0;
    }

    inline real W_New_Quadratic(const real r) const
    {
        if (r <= 2 * h) return coef_new_quadratic * (two_thirds - 9.0 / 8.0 * r * r / h2 + 19.0 / 24.0 * pow(r / h, 3) - 5.0 / 32.0 * pow(r / h, 4));
        else return (real)0;
    }
    //////////////////////////////////////////////////////////////////////////
    ////gradients

    inline VectorD Grad_Poly6(const VectorD& vr) const
    {
        real r = vr.norm();
        //std::cout << "vr: " << vr << std::endl;
        //std::cout << "r: " << r << std::endl;
        //std::cout << "h: " << h << std::endl;
        if (r > 0 && r < h)
        {
            return coef_poly6_grad * pow(h2 - r * r, 2) * vr;
        }
        else { return VectorD::Zero(); }
    }

    inline VectorD Grad_Spiky(const VectorD& vr) const
    {
        real r = vr.norm();
        if (r > 0 && r < h) return coef_spiky_grad * pow(h - r, 2) * vr / r;
        else return VectorD::Zero();
    }

    inline VectorD Grad_Vis(const VectorD& vr) const
    {
        real r = vr.norm();
        if (r > 0 && r < h) return coef_vis_grad_1 * (coef_vis_grad_2 * r + coef_vis_grad_3 - h / (2 * r * r * r)) * vr;
        else return VectorD::Zero();
    }

    inline VectorD Grad_Cubic(const VectorD& vr) const
    {
        real r = vr.norm();
        real u = r / h;
        // if (u>0&&u<1) return coef_cubic_spline*(1.5*r/h3-2.0/h2)*vr;
        // else if (u>=1&&u<2) return coef_cubic_spline*(-0.5*hd1*pow(2.0-r*hd1,2))*vr/r;
        const VectorD gradu = vr / (r * h);
        if (u > 0 && u <= 0.5) return 6.0 * coef_cubic_spline * u * (3.0 * u - 2.0) * gradu;
        else if (u > 0.5 && u < 1) return -6.0 * coef_cubic_spline * pow(1.0 - u, 2) * gradu;
        else return VectorD::Zero();
    }

    inline VectorD Grad_Quintic(const VectorD& vr) const
    {
        real r = vr.norm();
        if (r > 0 && r < h13) {
            return vr / r * coef_quintic * hd1 * (-(real)15 * pow(3 - 3 * r * hd1, 4) + (real)90 * pow(2 - 3 * r * hd1, 4) - (real)225 * pow(1 - 3 * r * hd1, 4));
        }
        else if (r >= h13 && r < h23) {
            return vr / r * coef_quintic * hd1 * (-(real)15 * pow(3 - 3 * r * hd1, 4) + (real)90 * pow(2 - 3 * r * hd1, 4));
        }
        else if (r >= h23 && r < h) {
            return vr / r * coef_quintic * hd1 * (-(real)15 * pow(3 - 3 * r * hd1, 4));
        }
        else { return VectorD::Zero(); }
    }

    inline VectorD Grad_Bell(const VectorD& vr) const
    {
        real r = vr.norm();
        if (r > 0 && r <= h) return coef_bell_grad * (-12.0 / h2) * pow(1.0 - r * hd1, 2) * vr;
        else return VectorD::Zero();
    }

    inline VectorD Grad_Gaussian(const VectorD& vr) const
    {
        real r = vr.norm();
        //return coef_gaussian_grad * (-2 / h2) * exp(-pow(r / h, 2)) * vr;
        return coef_gaussian_grad * exp(-pow(gaussian_trunc * r / h, 2)) * vr;
    }

    // ---This kernel seems to be behaving erratically -- please check or delete -- (yitong)
    //inline VectorD Grad_Quadratic(const VectorD& vr) const
    //{
    //	real r = vr.norm();
    //	if (r > 0 && r <= 2 * h) return coef_quadratic * (3.0 * r / (8.0 * h2) - 0.75 * hd1) * vr / r;
    //	else return VectorD::Zero();
    //}

    inline VectorD Grad_New_Quadratic(const VectorD& vr) const
    {
        real r = vr.norm();
        if (r > 0 && r <= 2 * h) return coef_quadratic * (19.0 * r / (8.0 * h3) - 9.0 / (4.0 * h2) - 5.0 * r * r / (8.0 * h3 * h)) * vr;
        else return VectorD::Zero();
    }
    //////////////////////////////////////////////////////////////////////////
    ////Laplacian
    inline real Lap_Poly6(const VectorD& vr) const
    {
        real r = vr.norm();
        if (r <= h) return coef_poly6_lap * (h2 - r * r) * (3 * h2 - 7 * r * r);
        else return (real)0;
    }

    inline real Lap_Spiky(const VectorD& vr) const
    {
        real r = vr.norm();
        if (r > 0 && r <= h) return coef_spiky_lap * (h - r) * (h - 2.0 * r) / r;
        else return (real)0;
    }

    inline real Lap_Vis(const VectorD& vr) const
    {
        real r = vr.norm();
        if (r < h) return coef_vis_lap * (h - r) / h3;
        else return (real)0;
    }
};
#endif
