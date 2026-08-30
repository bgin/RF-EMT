/* Copyright (C) Bernard Gingold, 2020-2026 
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

#ifndef __GMS_KDE_FFT_H__
#define __GMS_KDE_FFT_H__


namespace file_info 
{

     static const unsigned int GMS_KDE_FFT_MAJOR = 1;
     static const unsigned int GMS_KDE_FFT_MINOR = 1;
     static const unsigned int GMS_KDE_FFT_MICRO = 0;
     static const unsigned int GMS_KDE_FFT_FULLVER =
       1000U*GMS_KDE_FFT_MAJOR+100U*GMS_KDE_FFT_MINOR+
       10U*GMS_KDE_FFT_MICRO;
     static const char GMS_KDE_FFT_CREATION_DATE[] = "10-12-2025 01:38PM +00200 (WED 10 DEC 2025 GMT+2)";
     static const char GMS_KDE_FFT_BUILD_DATE[]    = __DATE__; 
     static const char GMS_KDE_FFT_BUILD_TIME[]    = __TIME__;
     static const char GMS_KDE_FFT_SYNOPSIS[]      = "Algorithm AS 176: Kernel Density Estimation Using the Fast Fourier Transform (written by B.W. Silverman), ported to C++.";

}

#include <cstdint>
#include "GMS_config.h"

namespace gms
{

namespace math 
{
   
__ATTR_HOT__
__ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
__ATTR_OPTIMIZE_03__
#endif    
void denest(float * __restrict__,
            const int32_t,
            float,
            float,
            float,
            float * __restrict__,
            float * __restrict__,
            const int32_t,
            int32_t,
            int32_t * __restrict__);

} //math

} //gms

#endif /*__GMS_KDE_FFT_H__*/
