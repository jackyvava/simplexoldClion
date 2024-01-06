//
// Created by zhangjie on 2024/1/6.
//

#ifndef SIMPLEXOLDCLION_KERNEL_H
#define SIMPLEXOLDCLION_KERNEL_H
#include <cmath>
#include "Common.h"
#include "Constants.h"
namespace Kernel {
    class Kernel {
    public:
        virtual real calValue(real x) = 0;
    };

    class GaussianKernel: public Kernel {
        real sigma; real sigma_sq; real one_over_sigma_sqrt2pi;
    public:
        GaussianKernel(real sigma_ = 1) : sigma(sigma_) { sigma_sq = sigma * sigma; one_over_sigma_sqrt2pi = 1.0 / sigma / std::sqrt(2 * pi); }
        real calValue(real x) override { return one_over_sigma_sqrt2pi * std::exp(-0.5*x*x/sigma_sq); }
    };
}

#endif

