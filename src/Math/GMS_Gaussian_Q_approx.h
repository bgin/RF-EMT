
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

} //math

} //gms

#endif /*__GMS_GAUSSIAN_Q_APPROX_H__*/