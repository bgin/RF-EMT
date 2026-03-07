/*MIT License
!Copyright (c) 2020 Bernard Gingold
!Permission is hereby granted, free of charge, to any person obtaining a copy
!of this software and associated documentation files (the "Software"), to deal
!in the Software without restriction, including without limitation the rights
!to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!copies of the Software, and to permit persons to whom the Software is
!furnished to do so, subject to the following conditions:
!The above copyright notice and this permission notice shall be included in all
!copies or substantial portions of the Software.
!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!SOFTWARE.
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