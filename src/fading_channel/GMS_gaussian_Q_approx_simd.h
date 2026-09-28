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

#ifndef __GMS_GAUSSIAN_Q_APPROX_SIMD_H__
#define __GMS_GAUSSIAN_Q_APPROX_SIMD_H__ 280920260802

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_fast_simd_funcs_approx.h"
#include "GMS_simd_utils.h"

namespace file_info 
{

     static const unsigned int GMS_GAUSSIAN_Q_APPROX_SIMD_MAJOR = 1;
     static const unsigned int GMS_GAUSSIAN_Q_APPROX_SIMD_MINOR = 1;
     static const unsigned int GMS_GAUSSIAN_Q_APPROX_SIMD_MICRO = 0;
     static const unsigned int GMS_GAUSSIAN_Q_APPROX_SIMD_FULLVER =
       1000U*GMS_GAUSSIAN_Q_APPROX_SIMD_MAJOR+100U*GMS_GAUSSIAN_Q_APPROX_SIMD_MINOR+
       10U*GMS_GAUSSIAN_Q_APPROX_SIMD_MICRO;
     static const char GMS_GAUSSIAN_Q_APPROX_SIMD_CREATION_DATE[] = "28-09-2026 08:03AM +00200 (MON 28 SEP 2026 GMT+2)";
     static const char GMS_GAUSSIAN_Q_APPROX_SIMD_BUILD_DATE[]    = __DATE__; 
     static const char GMS_GAUSSIAN_Q_APPROX_SIMD_BUILD_TIME[]    = __TIME__;
     static const char GMS_GAUSSIAN_Q_APPROX_SIMD_SYNOPSIS[]      = "Gaussian Q function approximations, based on the paper: Gaussian Q Function Approximation in Wireless \
                                                                Communication System Design: A \
                                                                Gradient-Based Optimization Approach ";

}

#if !defined(GAUSSIAN_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE)
#define GAUSSIAN_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE 0
#endif

namespace gms 
{

namespace math 
{

#if (GAUSSIAN_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
__m128 
gaussian_Q_approx_chiani_4xf32(const __m128 x)
{
const __m128 C0083333333333333333333333333 = _mm_set1_ps(0.083333333333333333333333333f);
const __m128 CN05                          = _mm_set1_ps(-0.5f);
const __m128 CN03                          = _mm_set1_ps(-0.333333333333333333333333333f);
const __m128 C025                          = _mm_set1_ps(0.25f);
const __m128 xx                            = _mm_mul_ps(x,x);
const __m128 left_exp_val                  = simd_fast_exp_approx_4xf32(_mm_mul_ps(CN05,xx));
const __m128 right_exp_val                 = _mm_mul_ps(C025,simd_fast_exp_approx_4xf32(_mm_mul_ps(CN03,_mm_add_ps(xx,xx))));
return (_mm_add_ps(left_exp_val,right_exp_val));
}

#if (GAUSSIAN_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
__m128d 
gaussian_Q_approx_chiani_2xf64(const __m128d x)
{
const __m128d C0083333333333333333333333333 = _mm_set1_pd(0.083333333333333333333333333);
const __m128d CN05                          = _mm_set1_pd(-0.5);
const __m128d CN03                          = _mm_set1_pd(-0.333333333333333333333333333);
const __m128d C025                          = _mm_set1_pd(0.25);
const __m128d xx                            = _mm_mul_pd(x,x);
const __m128d left_exp_val                  = simd_fast_exp_approx_2xf64(_mm_mul_pd(CN05,xx));
const __m128d right_exp_val                 = _mm_mul_pd(C025,simd_fast_exp_approx_2xf64(_mm_mul_pd(CN03,_mm_add_pd(xx,xx))));
return (_mm_add_pd(left_exp_val,right_exp_val));
}

#if (GAUSSIAN_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX2
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx2")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
__m256 
gaussian_Q_approx_chiani_8xf32(const __m256 x)
{
const __m256 C0083333333333333333333333333 = _mm256_set1_ps(0.083333333333333333333333333f);
const __m256 CN05                          = _mm256_set1_ps(-0.5f);
const __m256 CN03                          = _mm256_set1_ps(-0.333333333333333333333333333f);
const __m256 C025                          = _mm256_set1_ps(0.25f);
const __m256 xx                            = _mm256_mul_ps(x,x);
const __m256 left_exp_val                  = _mm256_exp_ps(_mm256_mul_ps(CN05,xx));
const __m256 right_exp_val                 = _mm256_mul_ps(C025,_mm256_exp_ps(_mm256_mul_ps(CN03,_mm256_add_ps(xx,xx))));
return (_mm256_add_ps(left_exp_val,right_exp_val));
}

#if (GAUSSIAN_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX2
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx2")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
__m256d 
gaussian_Q_approx_chiani_4xf64(const __m256d x)
{
const __m256d C0083333333333333333333333333 = _mm256_set1_pd(0.083333333333333333333333333);
const __m256d CN05                          = _mm256_set1_pd(-0.5);
const __m256d CN03                          = _mm256_set1_pd(-0.333333333333333333333333333);
const __m256d C025                          = _mm256_set1_pd(0.25);
const __m256d xx                            = _mm256_mul_pd(x,x);
const __m256d left_exp_val                  = _mm256_exp_pd(_mm256_mul_pd(CN05,xx));
const __m256d right_exp_val                 = _mm256_mul_pd(C025,_mm256_exp_pd(_mm256_mul_pd(CN03,_mm256_add_pd(xx,xx))));
return (_mm256_add_pd(left_exp_val,right_exp_val));
}

#if (GAUSSIAN_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
__m512
gaussian_Q_approx_chiani_16xf32(const __m512 x)
{
const __m512 C0083333333333333333333333333 = _mm512_set1_ps(0.083333333333333333333333333f);
const __m512 CN05                          = _mm512_set1_ps(-0.5f);
const __m512 CN03                          = _mm512_set1_ps(-0.333333333333333333333333333f);
const __m512 C025                          = _mm512_set1_ps(0.25f);
const __m512 xx                            = _mm512_mul_ps(x,x);
const __m512 left_exp_val                  = _mm512_exp_ps(_mm512_mul_ps(CN05,xx));
const __m512 right_exp_val                 = _mm512_mul_ps(C025,_mm512_exp_ps(_mm512_mul_ps(CN03,_mm512_add_ps(xx,xx))));
return (_mm512_add_ps(left_exp_val,right_exp_val));
}



}

}
#endif /*__GMS_GAUSSIAN_Q_APPROX_SIMD_H__*/