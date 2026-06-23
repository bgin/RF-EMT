/*******************************************************************************
*
*
*  Copyright [2019 - 2023] [Intel Corporation]
* 
*  Licensed under the Apache License, Version 2.0 (the "License");
*  you may not use this file except in compliance with the License.
*  
*  You may obtain a copy of the License at
*  
*     http://www.apache.org/licenses/LICENSE-2.0 
*  
*  Unless required by applicable law or agreed to in writing, software 
*  distributed under the License is distributed on an "AS IS" BASIS, 
*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*  See the License for the specific language governing permissions and 
*  limitations under the License. 
*  
*  SPDX-License-Identifier: Apache-2.0 
*  
* 
*
*******************************************************************************/

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

#ifndef __GMS_MATRIX_INV_CHOLESKY_H__
#define __GMS_MATRIX_INV_CHOLESKY_H__  230620261010

#include <cstdlib>
#include <cstdio>
#include <cstdint>
#include <cmath>
#include <immintrin.h>
#include "GMS_config.h"

namespace file_info 
{

     static const unsigned int GMS_MATRIX_INV_CHOLESKY_MAJOR = 1;
     static const unsigned int GMS_MATRIX_INV_CHOLESKY_MINOR = 1;
     static const unsigned int GMS_MATRIX_INV_CHOLESKY_MICRO = 0;
     static const unsigned int GMS_MATRIX_INV_CHOLESKY_FULLVER =
       1000U*GMS_MATRIX_INV_CHOLESKY_MAJOR+100U*GMS_MATRIX_INV_CHOLESKY_MINOR+
       10U*GMS_MATRIX_INV_CHOLESKY_MICRO;
     static const char GMS_MATRIX_INV_CHOLESKY_CREATION_DATE[] = "23-06-2026 10:12AM +00200 (TUE 23 JUN 2026 GMT+2)";
     static const char GMS_MATRIX_INV_CHOLESKY_BUILD_DATE[]    = __DATE__; 
     static const char GMS_MATRIX_INV_CHOLESKY_BUILD_TIME[]    = __TIME__;
     static const char GMS_MATRIX_INV_CHOLESKY_SYNOPSIS[]      = "Further adaptation and optimization of avx512 cholesky inverse kernels (FlexRAN).";

}

namespace gms 
{

namespace math 
{

namespace mat_inv_chol_detail
{

constexpr std::int32_t NUM_MEM        = 100;
constexpr std::int32_t MAT_SQR_SIZE_2 = 2;
constexpr std::int32_t MAT_SQR_SIZE_3 = 3;
constexpr std::int32_t MAT_SQR_SIZE_4 = 4; 
constexpr std::int32_t MAT_SQR_SIZE_5 = 5; 
constexpr std::int32_t MAT_SQR_SIZE_6 = 6; 
constexpr std::int32_t MAT_SQR_SIZE_7 = 7;
constexpr std::int32_t MAT_SQR_SIZE_8 = 8;
constexpr std::int32_t MAT_SQR_SIZE_9 = 9;
constexpr std::int32_t MAT_SQR_SIZE_10= 10;
constexpr std::int32_t MAT_SQR_SIZE_11= 11;
constexpr std::int32_t MAT_SQR_SIZE_12= 12;
constexpr std::int32_t MAT_SQR_SIZE_13= 13;
constexpr std::int32_t MAT_SQR_SIZE_14= 14;
constexpr std::int32_t MAT_SQR_SIZE_15= 15;
constexpr std::int32_t MAT_SQR_SIZE_16= 16;

}

using namespace mat_inv_chol_detail;

__ATTR_HOT__
__ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
void 
mat_inv_cholesky_2x2_16xf32(const __m512 matBRe[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2],
                            const __m512 matBIm[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2],
                            __m512       matInvBRe[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2],
                            __m512       matInvBIm[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2]);

__ATTR_HOT__
__ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
void 
mat_inv_cholesky_3x3_16xf32(const __m512 matBRe[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3],
                            const __m512 matBIm[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3],
                            __m512       matInvBRe[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3],
                            __m512       matInvBIm[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3]);

__ATTR_HOT__
__ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
void 
mat_inv_cholesky_4x4_16xf32(const __m512 matBRe[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4],
                            const __m512 matBIm[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4],
                            __m512       matInvBRe[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4],
                            __m512       matInvBIm[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4]);

__ATTR_HOT__
__ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
void 
mat_inv_cholesky_5x5_16xf32(const __m512 matBRe[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5],
                            const __m512 matBIm[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5],
                            __m512       matInvBRe[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5],
                            __m512       matInvBIm[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5]);

__ATTR_HOT__
__ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
void 
mat_inv_cholesky_6x6_16xf32(const __m512 matBRe[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6],
                            const __m512 matBIm[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6],
                            __m512       matInvBRe[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6],
                            __m512       matInvBIm[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6]);

__ATTR_HOT__
__ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
void 
mat_inv_cholesky_7x7_16xf32(const __m512 matBRe[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7],
                            const __m512 matBIm[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7],
                            __m512       matInvBRe[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7],
                            __m512       matInvBIm[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7]);

__ATTR_HOT__
__ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
void 
mat_inv_cholesky_8x8_16xf32(const __m512 matBRe[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8],
                            const __m512 matBIm[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8],
                            __m512       matInvBRe[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8],
                            __m512       matInvBIm[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8]);

__ATTR_HOT__
__ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
void 
mat_inv_cholesky_9x9_16xf32(const __m512 matBRe[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9],
                            const __m512 matBIm[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9],
                            __m512       matInvBRe[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9],
                            __m512       matInvBIm[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9]);

__ATTR_HOT__
__ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
void 
mat_inv_cholesky_10x10_16xf32(const __m512 matBRe[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10],
                              const __m512 matBIm[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10],
                              __m512       matInvBRe[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10],
                              __m512       matInvBIm[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10]);

__ATTR_HOT__
__ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
void 
mat_inv_cholesky_11x11_16xf32(const __m512 matBRe[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11],
                              const __m512 matBIm[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11],
                              __m512       matInvBRe[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11],
                              __m512       matInvBIm[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11]);

__ATTR_HOT__
__ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
void 
mat_inv_cholesky_12x12_16xf32(const __m512 matBRe[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12],
                              const __m512 matBIm[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12],
                              __m512       matInvBRe[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12],
                              __m512       matInvBIm[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12]);

__ATTR_HOT__
__ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
void 
mat_inv_cholesky_13x13_16xf32(const __m512 matBRe[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13],
                              const __m512 matBIm[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13],
                              __m512       matInvBRe[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13],
                              __m512       matInvBIm[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13]);

__ATTR_HOT__
__ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
void 
mat_inv_cholesky_14x14_16xf32(const __m512 matBRe[MAT_SQR_SIZE_14][MAT_SQR_SIZE_14],
                              const __m512 matBIm[MAT_SQR_SIZE_14][MAT_SQR_SIZE_14],
                              __m512       matInvBRe[MAT_SQR_SIZE_14][MAT_SQR_SIZE_14],
                              __m512       matInvBIm[MAT_SQR_SIZE_14][MAT_SQR_SIZE_14]);

__ATTR_HOT__
__ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
void 
mat_inv_cholesky_15x15_16xf32(const __m512 matBRe[MAT_SQR_SIZE_15][MAT_SQR_SIZE_15],
                              const __m512 matBIm[MAT_SQR_SIZE_15][MAT_SQR_SIZE_15],
                              __m512       matInvBRe[MAT_SQR_SIZE_15][MAT_SQR_SIZE_15],
                              __m512       matInvBIm[MAT_SQR_SIZE_15][MAT_SQR_SIZE_15]);

__ATTR_HOT__
__ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
void 
mat_inv_cholesky_16x16_16xf32(const __m512 matBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                              const __m512 matBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                              __m512       matInvBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                              __m512       matInvBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16]);


} // math 

} // gms












#endif /*__GMS_MATRIX_INV_CHOLESKY_H__*/