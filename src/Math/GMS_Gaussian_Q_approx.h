
/*
 * Copyright (C) Bernard Gingold, 2020-2026 
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#ifndef __GMS_GAUSSIAN_Q_APPROX_H__
#define __GMS_GAUSSIAN_Q_APPROX_H__ 160920261342

#include <cstdint>
#include <algorithm>
#include "GMS_config.h"
#include "GMS_fast_simd_funcs_approx.h"

namespace file_info 
{

     static const unsigned int GMS_GAUSSIAN_Q_APPROX_MAJOR = 1;
     static const unsigned int GMS_GAUSSIAN_Q_APPROX_MINOR = 1;
     static const unsigned int GMS_GAUSSIAN_Q_APPROX_MICRO = 0;
     static const unsigned int GMS_GAUSSIAN_Q_APPROX_FULLVER =
       1000U*GMS_GAUSSIAN_Q_APPROX_MAJOR+100U*GMS_GAUSSIAN_Q_APPROX_MINOR+
       10U*GMS_GAUSSIAN_Q_APPROX_MICRO;
     static const char GMS_GAUSSIAN_Q_APPROX_CREATION_DATE[] = "16-09-2026 13:45AM +00200 (WED 16 AUG 2026 GMT+2)";
     static const char GMS_GAUSSIAN_Q_APPROX_BUILD_DATE[]    = __DATE__; 
     static const char GMS_GAUSSIAN_Q_APPROX_BUILD_TIME[]    = __TIME__;
     static const char GMS_GAUSSIAN_Q_APPROX_SYNOPSIS[]      = "Gaussian Q function approximations, based on the paper: Gaussian Q Function Approximation in Wireless \
                                                                Communication System Design: A \
                                                                Gradient-Based Optimization Approach ";

}

namespace gms 
{

namespace math 
{

__ATTR_ALWAYS_INLINE__
static inline
float gaussian_Q_approx_chiani(const float x)
{
    constexpr float C0083333333333333333333333333 = 0.083333333333333333333333333;
    constexpr float C05                           = 0.5f;
    constexpr float C03                           = 0.333333333333333333333333333;
    constexpr float C025                          = 0.25f;
    const float xx            = x*x;
    const float left_exp_val  = C0083333333333333333333333333*expapprox(-C05*xx);
    const float right_exp_val = C025*expapprox(-C03*(2.0*xx));
    return (left_exp_val+right_exp_val);
}

__ATTR_ALWAYS_INLINE__
static inline
double gaussian_Q_approx_chiani(const double x)
{
    constexpr double C0083333333333333333333333333 = 0.083333333333333333333333333;
    constexpr double C05                           = 0.5;
    constexpr double C03                           = 0.333333333333333333333333333;
    constexpr double C025                          = 0.25;
    const double xx            = x*x;
    const double left_exp_val  = C0083333333333333333333333333*expapprox_d(-C05*xx);
    const double right_exp_val = C025*expapprox_d(-C03*(2.0*xx));
    return (left_exp_val+right_exp_val);
}

__ATTR_ALWAYS_INLINE__
static inline 
float gaussian_Q_approx_loskot_2T(const float x)
{
    constexpr float C0208  = 0.208f;
    constexpr float C0971  = 0.971f;
    constexpr float C0147  = 0.147f;
    constexpr float C0525  = 0.525f;
    const float xx = x*x;
    return (std::fma(C0208,expapprox(-C0971*xx),C0147*expapprox(-C0525*xx)));
}

__ATTR_ALWAYS_INLINE__
static inline 
double gaussian_Q_approx_loskot_2T(const double x)
{
    constexpr double C0208  = 0.208;
    constexpr double C0971  = 0.971;
    constexpr double C0147  = 0.147;
    constexpr double C0525  = 0.525;
    const double xx = x*x;
    return (std::fma(C0208,expapprox_d(-C0971*xx),C0147*expapprox_d(-C0525*xx)));
}

__ATTR_ALWAYS_INLINE__
static inline 
float gaussian_Q_approx_loskot_3T(const float x)
{
    constexpr float C0168 = 0.168f;
    constexpr float C0876 = 0.876f;
    constexpr float C0144 = 0.144f;
    constexpr float C0525 = 0.525f;
    constexpr float C0002 = 0.002f;
    constexpr float C0603 = 0.603f;
    const float xx = x*x;
    return (std::fma(C0168,expapprox(-C0876*xx),
                           std::fma(C0144,expapprox(-C0525*xx),C0002*expapprox(-C0603*xx))));
}

__ATTR_ALWAYS_INLINE__
static inline 
double gaussian_Q_approx_loskot_3T(const double x)
{
    constexpr double C0168 = 0.168;
    constexpr double C0876 = 0.876;
    constexpr double C0144 = 0.144;
    constexpr double C0525 = 0.525;
    constexpr double C0002 = 0.002;
    constexpr double C0603 = 0.603;
    const double xx = x*x;
    return (std::fma(C0168,expapprox_d(-C0876*xx),
                           std::fma(C0144,expapprox_d(-C0525*xx),C0002*expapprox_d(-C0603*xx))));
}

__ATTR_ALWAYS_INLINE__
static inline 
float gaussian_Q_approx_sadhwani_1T(const float x)
{
    constexpr float C025 = 0.25f;
    constexpr float C05  = 0.5f;
    const float xx = x*x;
    return (C025*expapprox(-C05*xx));
}

__ATTR_ALWAYS_INLINE__
static inline 
double gaussian_Q_approx_sadhwani_1T(const double x)
{
    constexpr double C025 = 0.25;
    constexpr double C05  = 0.5;
    const double xx = x*x;
    return (C025*expapprox_d(-C05*xx));
}

__ATTR_ALWAYS_INLINE__
static inline 
float gaussian_Q_approx_sadhwani_2T(const float x)
{
    constexpr float C0125 = 0.125f;
    constexpr float C025  = 0.25f;
    constexpr float C05   = 0.5f;
    const float xx  = x*x;
    return (std::fma(C0125,expapprox(-xx),C025*expapprox(-C05*xx)));
}

__ATTR_ALWAYS_INLINE__
static inline 
double gaussian_Q_approx_sadhwani_2T(const double x)
{
    constexpr double C0125 = 0.125;
    constexpr double C025  = 0.25;
    constexpr double C05   = 0.5;
    const double xx = x*x;
    return (std::fma(C0125,expapprox_d(-xx),C025*expapprox_d(-C05*xx)));
}



} //math

} //gms

#endif /*__GMS_GAUSSIAN_Q_APPROX_H__*/