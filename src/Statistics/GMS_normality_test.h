/*
!        ALGORITHM AS R94 APPL. STATIST. (1995) VOL.44, NO.4
!        Calculates the Shapiro-Wilk W test and its significance level
! Fortran 90 version by Alan.Miller @ vic.cmis.csiro.au
! Latest revision - 4 December 1998
! C++ Port by Bernard Gingold beniekg@gmail.com On 13 October 2025.
! 
*/

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

#ifndef __GMS_NORMALITY_TEST_H__
#define __GMS_NORMALITY_TEST_H__ 131020250657

namespace file_info 
{

     static const unsigned int GMS_GMS_NORMALITY_TEST_MAJOR = 1;
     static const unsigned int GMS_GMS_NORMALITY_TEST_MINOR = 1;
     static const unsigned int GMS_GMS_NORMALITY_TEST_MICRO = 0;
     static const unsigned int GMS_GMS_NORMALITY_TEST_FULLVER =
       1000U*GMS_GMS_NORMALITY_TEST_MAJOR+100U*GMS_GMS_NORMALITY_TEST_MINOR+
       10U*GMS_GMS_NORMALITY_TEST_MICRO;
     static const char GMS_GMS_NORMALITY_TEST_CREATION_DATE[] =   "13-10-2025 06:59AM +00200 (MON 13 OCT 2025 GMT+2)";
     static const char GMS_AM_BB_CMPLX_SINE_SIGNAL_BUILD_DATE[]    = __DATE__; 
     static const char GMS_AM_BB_CMPLX_SINE_SIGNAL_BUILD_TIME[]    = __TIME__;
     static const char GMS_AM_BB_CMPLX_SINE_SIGNAL_SYNOPSIS[]      = "ALGORITHM AS R94 APPL. STATIST. (1995) VOL.44, NO.4 (C++ Port)";

}

#include <cstdint>
#include "GMS_config.h"

namespace gms 
{

namespace math 
{



__ATTR_HOT__
__ATTR_OPTIMIZE_03__
void         shapiro_wilk(bool,
                          float * __restrict,
                          const int32_t ,
                          const int32_t ,
                          const int32_t ,
                          float * __restrict ,
                          float &,
                          float &,
                          int32_t &); 



} // math

} // gms

#endif /*__GMS_NORMALITY_TEST_H__*/

