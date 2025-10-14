
/*
    Adapted from: 
    [freely-available stand-alone C++
    implementation](https://www.johndcook.com/blog/cpp_phi_inverse/).
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

#ifndef __GMS_NORMALITY_TEST_V2_H__
#define __GMS_NORMALITY_TEST_V2_H__ 141020251201

namespace file_info 
{

     static const unsigned int GMS_NORMALITY_TEST_V2_MAJOR = 1;
     static const unsigned int GMS_NORMALITY_TEST_V2_MINOR = 1;
     static const unsigned int GMS_NORMALITY_TEST_V2_MICRO = 0;
     static const unsigned int GMS_NORMALITY_TEST_V2_FULLVER =
       1000U*GMS_NORMALITY_TEST_V2_MAJOR+100U*GMS_NORMALITY_TEST_V2_MINOR+
       10U*GMS_NORMALITY_TEST_V2_MICRO;
     static const char GMS_NORMALITY_TEST_V2_CREATION_DATE[] =   "14-10-2025 12:02PM +00200 (TUE 14 OCT 2025 GMT+2)";
     static const char GMS_NORMALITY_TEST_V2_BUILD_DATE[]    = __DATE__; 
     static const char GMS_NORMALITY_TEST_V2_BUILD_TIME[]    = __TIME__;
     static const char GMS_NORMALITY_TEST_V2_SYNOPSIS[]      = "Adapted from: https://www.johndcook.com/blog/cpp_phi_inverse/";

}

#include <cstdint>
#include "GMS_config.h"
#include "GMS_dyn_array.h"

namespace gms 
{

namespace math 
{

__ATTR_HOT__
__ATTR_OPTIMIZE_03__
double 
calculate_W_statistics(const darray_r8_t &);

}

}

#endif /*__GMS_NORMALITY_TEST_V2_H__*/