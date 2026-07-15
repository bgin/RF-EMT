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

#include <cstdint>
#include <cmath>
#include <immintrin.h>
#include <string>
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

namespace mat_inv_chol_tsc_instr
{

static constexpr std::uint64_t RDTSCP_SKYLAKE_SERVER_LAT = 42ull; // accordingly to intrinsic guide

template<std::size_t MAT_SQR_SIZE_MAX>
struct alignas(64) tsc_instrumentation_block_t
{

  static_assert(MAT_SQR_SIZE_MAX<=16ull,"The allowable size (i.e. 16) was exceeded!!");
/* The members are laid-out accordingly to the instrumented regions, except the array members
   which are declared first.
*/
  __ATTR_ALIGN__(64) std::uint64_t 
  m_region_loop2D_s[MAT_SQR_SIZE_MAX] = {};
  __ATTR_ALIGN__(64)  std::uint64_t
  m_region_loop2D_e[MAT_SQR_SIZE_MAX] = {};
  __ATTR_ALIGN__(64)  std::uint64_t 
  m_region_loop2D_d[MAT_SQR_SIZE_MAX] = {};

  __ATTR_ALIGN__(64)  std::uint64_t 
  m_region_loop3D_s[MAT_SQR_SIZE_MAX] = {};
  __ATTR_ALIGN__(64)  std::uint64_t 
  m_region_loop3D_e[MAT_SQR_SIZE_MAX] = {};
  __ATTR_ALIGN__(64)  std::uint64_t 
  m_region_loop3D_d[MAT_SQR_SIZE_MAX] = {};
  /*Scalar members*/
  
   std::uint64_t m_region_prefetch_s{};
   std::uint64_t m_region_prefetch_e{};
   std::uint64_t m_region_prefetch_d{};
   std::uint64_t m_region_g00_s{};
   std::uint64_t m_region_g00_e{};
   std::uint64_t m_region_g00_d{}; 
   std::uint64_t m_region_g11_s{};
   std::uint64_t m_region_g11_e{};
   std::uint64_t m_region_g11_d{};
   std::uint64_t m_region_g22_s{};
   std::uint64_t m_region_g22_e{};
   std::uint64_t m_region_g22_d{};
   std::uint64_t m_region_g33_s{};
   std::uint64_t m_region_g33_e{};
   std::uint64_t m_region_g33_d{};   
   std::uint64_t m_region_g44_s{};
   std::uint64_t m_region_g44_e{};
   std::uint64_t m_region_g44_d{};
   std::uint64_t m_region_g55_s{};
   std::uint64_t m_region_g55_e{};
   std::uint64_t m_region_g55_d{};  
   std::uint64_t m_region_g66_s{};
   std::uint64_t m_region_g66_e{};
   std::uint64_t m_region_g66_d{}; 
   std::uint64_t m_region_g77_s{};
   std::uint64_t m_region_g77_e{};
   std::uint64_t m_region_g77_d{};
   std::uint64_t m_region_gii_even_col8_s{};
   std::uint64_t m_region_gii_even_col8_e{};
   std::uint64_t m_region_gii_even_col8_d{};  
   std::uint64_t m_region_gii_odd_col9_s{};
   std::uint64_t m_region_gii_odd_col9_e{};
   std::uint64_t m_region_gii_odd_col9_d{}; 
   std::uint64_t m_region_gii_even_col10_s{};
   std::uint64_t m_region_gii_even_col10_e{};
   std::uint64_t m_region_gii_even_col10_d{}; 
   std::uint64_t m_region_gii_odd_col11_s{};
   std::uint64_t m_region_gii_odd_col11_e{};
   std::uint64_t m_region_gii_odd_col11_d{}; 
   std::uint64_t m_region_gii_even_col12_s{};
   std::uint64_t m_region_gii_even_col12_e{};
   std::uint64_t m_region_gii_even_col12_d{};
   std::uint64_t m_region_gii_odd_col13_s{};
   std::uint64_t m_region_gii_odd_col13_e{};
   std::uint64_t m_region_gii_odd_col13_d{};  
   std::uint64_t m_region_gii_even_col14_s{};
   std::uint64_t m_region_gii_even_col14_e{};
   std::uint64_t m_region_gii_even_col14_d{};
   std::uint64_t m_region_gii_odd_col15_s{};
   std::uint64_t m_region_gii_odd_col15_e{};
   std::uint64_t m_region_gii_odd_col15_d{};  
  //////////////////////////////////////////
   std::uint64_t m_region_lii_col0_s{};
   std::uint64_t m_region_lii_col0_e{};
   std::uint64_t m_region_lii_col0_d{};
   std::uint64_t m_region_lii_col1_s{};
   std::uint64_t m_region_lii_col1_e{};
   std::uint64_t m_region_lii_col1_d{};
   std::uint64_t m_region_lii_col2_s{};
   std::uint64_t m_region_lii_col2_e{};
   std::uint64_t m_region_lii_col2_d{};
   std::uint64_t m_region_lii_col3_s{};
   std::uint64_t m_region_lii_col3_e{};
   std::uint64_t m_region_lii_col3_d{};
   std::uint64_t m_region_lii_col4_s{};
   std::uint64_t m_region_lii_col4_e{};
   std::uint64_t m_region_lii_col4_d{};
   std::uint64_t m_region_lii_col5_s{};
   std::uint64_t m_region_lii_col5_e{};
   std::uint64_t m_region_lii_col5_d{};
   std::uint64_t m_region_lii_col6_s{};
   std::uint64_t m_region_lii_col6_e{};
   std::uint64_t m_region_lii_col6_d{};
   std::uint64_t m_region_lii_col7_s{};
   std::uint64_t m_region_lii_col7_e{};
   std::uint64_t m_region_lii_col7_d{};
   std::uint64_t m_region_lii_col8_s{};
   std::uint64_t m_region_lii_col8_e{};
   std::uint64_t m_region_lii_col8_d{};
   std::uint64_t m_region_lii_col9_s{};
   std::uint64_t m_region_lii_col9_e{};
   std::uint64_t m_region_lii_col9_d{};
   std::uint64_t m_region_lii_col10_s{};
   std::uint64_t m_region_lii_col10_e{};
   std::uint64_t m_region_lii_col10_d{};
   std::uint64_t m_region_lii_col11_s{};
   std::uint64_t m_region_lii_col11_e{};
   std::uint64_t m_region_lii_col11_d{};
   std::uint64_t m_region_lii_col12_s{};
   std::uint64_t m_region_lii_col12_e{};
   std::uint64_t m_region_lii_col12_d{};
   std::uint64_t m_region_lii_col13_s{};
   std::uint64_t m_region_lii_col13_e{};
   std::uint64_t m_region_lii_col13_d{};
   std::uint64_t m_region_lii_col14_s{};
   std::uint64_t m_region_lii_col14_e{};
   std::uint64_t m_region_lii_col14_d{};
   std::uint64_t m_region_lii_col15_s{};
   std::uint64_t m_region_lii_col15_e{};
   std::uint64_t m_region_lii_col15_d{};

};

__ATTR_ALWAYS_INLINE__ 
static inline 
std::uint64_t compute_delta(const std::uint64_t tsc_start,
                            const std::uint64_t tsc_end,
                            const bool do_latency_correction) 
{
    std::uint64_t tsc_s_corrected{};
    std::uint64_t tsc_e_corrected{};
    std::uint64_t tmp_res{};
    if(do_latency_correction)
    {
       tsc_s_corrected = tsc_start-RDTSCP_SKYLAKE_SERVER_LAT;
       tsc_e_corrected = tsc_end-RDTSCP_SKYLAKE_SERVER_LAT;
       tmp_res = tsc_e_corrected-tsc_s_corrected;
    }
    else 
    {
       tmp_res = tsc_end-tsc_start;
    }
    std::uint64_t delta{};
    (tmp_res>0ull)?delta = tmp_res:delta = UINT64_MAX;
    return (delta);
}

__ATTR_ALWAYS_INLINE__ 
static inline 
std::int32_t print_results(const std::string &&reg_name,
                           const std::uint64_t tsc_start,
                           const std::uint64_t tsc_end,
                           const std::uint64_t tsc_delta)
{
    std::int32_t 
    printf_ret = std::printf("[PMC-TSC:] region=%s,START=%llu,END=%llu,DELTA=%llu\n",reg_name.c_str(),tsc_start,tsc_end,tsc_delta);
    return (printf_ret);
}

}

using namespace mat_inv_chol_detail;


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_prefetching,bool mitigate_nan>         
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
mat_inv_cholesky_2x2_16xf32(const __m512 matBRe[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2],
                            const __m512 matBIm[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2],
                            __m512       matInvBRe[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2],
                            __m512       matInvBIm[MAT_SQR_SIZE_2][MAT_SQR_SIZE_2]);


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_prefetching,bool mitigate_nan>
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
mat_inv_cholesky_3x3_16xf32(const __m512 matBRe[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3],
                            const __m512 matBIm[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3],
                            __m512       matInvBRe[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3],
                            __m512       matInvBIm[MAT_SQR_SIZE_3][MAT_SQR_SIZE_3]);


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_prefetching,bool mitigate_nan>
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
mat_inv_cholesky_4x4_16xf32(const __m512 matBRe[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4],
                            const __m512 matBIm[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4],
                            __m512       matInvBRe[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4],
                            __m512       matInvBIm[MAT_SQR_SIZE_4][MAT_SQR_SIZE_4]);


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_prefetching,bool mitigate_nan>
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
mat_inv_cholesky_5x5_16xf32(const __m512 matBRe[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5],
                            const __m512 matBIm[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5],
                            __m512       matInvBRe[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5],
                            __m512       matInvBIm[MAT_SQR_SIZE_5][MAT_SQR_SIZE_5]);


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_prefetching,bool mitigate_nan>
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
mat_inv_cholesky_6x6_16xf32(const __m512 matBRe[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6],
                            const __m512 matBIm[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6],
                            __m512       matInvBRe[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6],
                            __m512       matInvBIm[MAT_SQR_SIZE_6][MAT_SQR_SIZE_6]);


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_prefetching,bool mitigate_nan>
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
mat_inv_cholesky_7x7_16xf32(const __m512 matBRe[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7],
                            const __m512 matBIm[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7],
                            __m512       matInvBRe[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7],
                            __m512       matInvBIm[MAT_SQR_SIZE_7][MAT_SQR_SIZE_7]);


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_prefetching,bool mitigate_nan>
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
mat_inv_cholesky_8x8_16xf32(const __m512 matBRe[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8],
                            const __m512 matBIm[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8],
                            __m512       matInvBRe[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8],
                            __m512       matInvBIm[MAT_SQR_SIZE_8][MAT_SQR_SIZE_8]);


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_prefetching,bool mitigate_nan>
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
mat_inv_cholesky_9x9_16xf32(const __m512 matBRe[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9],
                            const __m512 matBIm[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9],
                            __m512       matInvBRe[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9],
                            __m512       matInvBIm[MAT_SQR_SIZE_9][MAT_SQR_SIZE_9]);


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_prefetching,bool mitigate_nan>
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
mat_inv_cholesky_10x10_16xf32(const __m512 matBRe[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10],
                              const __m512 matBIm[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10],
                              __m512       matInvBRe[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10],
                              __m512       matInvBIm[MAT_SQR_SIZE_10][MAT_SQR_SIZE_10]);


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_prefetching,bool mitigate_nan>
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
mat_inv_cholesky_11x11_16xf32(const __m512 matBRe[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11],
                              const __m512 matBIm[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11],
                              __m512       matInvBRe[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11],
                              __m512       matInvBIm[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11]);


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_prefetching,bool mitigate_nan>
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
mat_inv_cholesky_12x12_16xf32(const __m512 matBRe[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12],
                              const __m512 matBIm[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12],
                              __m512       matInvBRe[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12],
                              __m512       matInvBIm[MAT_SQR_SIZE_12][MAT_SQR_SIZE_12]);


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_prefetching,bool mitigate_nan>
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
mat_inv_cholesky_13x13_16xf32(const __m512 matBRe[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13],
                              const __m512 matBIm[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13],
                              __m512       matInvBRe[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13],
                              __m512       matInvBIm[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13]);

// An attempt to force a constexpr evaluation of the SIMD loop.
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_prefetching,bool mitigate_nan,
         std::int32_t N_VECS_13>
__ATTR_HOT__
__ATTR_ALIGN__(32)
constexpr
void        
mat_inv_cholesky_13x13_16xf32_v2(const __m512 matBRe[N_VECS_13][N_VECS_13],
                                 const __m512 matBIm[N_VECS_13][N_VECS_13],
                                 __m512       matInvBRe[N_VECS_13][N_VECS_13],
                                 __m512       matInvBIm[N_VECS_13][N_VECS_13]);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_prefetching,bool mitigate_nan>
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
mat_inv_cholesky_14x14_16xf32(const __m512 matBRe[MAT_SQR_SIZE_14][MAT_SQR_SIZE_14],
                              const __m512 matBIm[MAT_SQR_SIZE_14][MAT_SQR_SIZE_14],
                              __m512       matInvBRe[MAT_SQR_SIZE_14][MAT_SQR_SIZE_14],
                              __m512       matInvBIm[MAT_SQR_SIZE_14][MAT_SQR_SIZE_14]);


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_prefetching,bool mitigate_nan>
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
mat_inv_cholesky_15x15_16xf32(const __m512 matBRe[MAT_SQR_SIZE_15][MAT_SQR_SIZE_15],
                              const __m512 matBIm[MAT_SQR_SIZE_15][MAT_SQR_SIZE_15],
                              __m512       matInvBRe[MAT_SQR_SIZE_15][MAT_SQR_SIZE_15],
                              __m512       matInvBIm[MAT_SQR_SIZE_15][MAT_SQR_SIZE_15]);


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_prefetching,bool mitigate_nan>
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
mat_inv_cholesky_16x16_16xf32(const __m512 matBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                              const __m512 matBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                              __m512       matInvBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                              __m512       matInvBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16]);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_prefetching,bool mitigate_nan>
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
mat_inv_cholesky_16x16_16xf32_tsc_instr(const __m512 matBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                        const __m512 matBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                        __m512       matInvBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                        __m512       matInvBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                        mat_inv_chol_tsc_instr::tsc_instrumentation_block_t<MAT_SQR_SIZE_16> &);                             


} // math 

} // gms












#endif /*__GMS_MATRIX_INV_CHOLESKY_H__*/