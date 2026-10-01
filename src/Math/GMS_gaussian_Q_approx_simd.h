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

#if !defined(GAUSSIAN_Q_APPROX_SIMD_RECIPROCAL_BY_DIVISION)
#define GAUSSIAN_Q_APPROX_SIMD_RECIPROCAL_BY_DIVISION 1
#endif 

namespace gms 
{

namespace math 
{

/*
#include <immintrin.h>
#include <cmath>

__m512 gaussian_Q_approx_sadhwani_1T_vec(__m512 x);

void marcum_Q_approx_sadhwani_1T_vector(
    const float* __restrict mu, 
    const float* __restrict a, 
    const float* __restrict b, 
    float* __restrict result, 
    int n) 
{
    // Constant for the 0.5f subtraction
    __m512 v_half = _mm512_set1_ps(0.5f);
    // Constant for the 1.0f subtraction in the 'else' branch
    __m512 v_one  = _mm512_set1_ps(1.0f);

    int i = 0;
    // Process 16 float elements at a time
    for (; i <= n - 16; i += 16) {
        // 1. Load inputs into 512-bit registers
        __m512 v_mu = _mm512_loadu_ps(&mu[i]);
        __m512 v_a  = _mm512_loadu_ps(&a[i]);
        __m512 v_b  = _mm512_loadu_ps(&b[i]);

        // 2. Perform unconditional shared math: pow_term = (b/a)^(mu - 0.5)
        __m512 v_base = _mm512_div_ps(v_b, v_a);
        __m512 v_exp  = _mm512_sub_ps(v_mu, v_half);
        // Requires SVML or an AVX-512 compatible math library linkage
        __m512 v_pow_term = _mm512_pow_ps(v_base, v_exp); 

        // 3. Prepare both conditional code execution branches simultaneously
        
        // Path A: b > a -> pow_term * gaussian_Q_approx(b - a)
        __m512 v_diff_A = _mm512_sub_ps(v_b, v_a);
        __m512 v_gauss_A = gaussian_Q_approx_sadhwani_1T_vec(v_diff_A);
        __m512 v_res_A = _mm512_mul_ps(v_pow_term, v_gauss_A);

        // Path B: a > b -> 1.0 - pow_term * gaussian_Q_approx(a - b)
        __m512 v_diff_B = _mm512_sub_ps(v_a, v_b);
        __m512 v_gauss_B = gaussian_Q_approx_sadhwani_1T_vec(v_diff_B);
        __m512 v_res_B = _mm512_sub_ps(v_one, _mm512_mul_ps(v_pow_term, v_gauss_B));

        // 4. If-Conversion: Generate a 16-bit condition mask evaluating (b > a)
        __mmask16 mask_gt = _mm512_cmp_ps_mask(v_b, v_a, _CMP_GT_OS);

        // 5. Masked Blend: Merge the paths based on the boolean condition mask
        // If a mask bit is 1 (b > a is true), pick from v_res_A. 
        // If a mask bit is 0 (b > a is false / a >= b), pick from v_res_B.
        __m512 v_final_res = _mm512_mask_blend_ps(mask_gt, v_res_B, v_res_A);

        // 6. Store the merged output
        _mm512_storeu_ps(&result[i], v_final_res);
    }

    // 7. Cleanup loop for any remaining structural array elements (n % 16)
    if (i < n) {
        int remaining = n - i;
        __mmask16 tail_mask = (1U << remaining) - 1;

        __m512 v_mu = _mm512_maskz_loadu_ps(tail_mask, &mu[i]);
        __m512 v_a  = _mm512_maskz_loadu_ps(tail_mask, &a[i]);
        __m512 v_b  = _mm512_maskz_loadu_ps(tail_mask, &b[i]);

        __m512 v_base = _mm512_div_ps(v_b, v_a);
        __m512 v_exp  = _mm512_sub_ps(v_mu, v_half);
        __m512 v_pow_term = _mm512_pow_ps(v_base, v_exp);

        __m512 v_res_A = _mm512_mul_ps(v_pow_term, gaussian_Q_approx_sadhwani_1T_vec(_mm512_sub_ps(v_b, v_a)));
        __m512 v_res_B = _mm512_sub_ps(v_one, _mm512_mul_ps(v_pow_term, gaussian_Q_approx_sadhwani_1T_vec(_mm512_sub_ps(v_a, v_b))));

        __mmask16 mask_gt = _mm512_mask_cmp_ps_mask(tail_mask, v_b, v_a, _CMP_GT_OS);
        __m512 v_final_res = _mm512_mask_blend_ps(mask_gt, v_res_B, v_res_A);

        _mm512_mask_storeu_ps(&result[i], tail_mask, v_final_res);
    }
}

*/

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
__m512d
gaussian_Q_approx_chiani_8xf64(const __m512d x)
{
const __m512d C0083333333333333333333333333 = _mm512_set1_pd(0.083333333333333333333333333);
const __m512d CN05                          = _mm512_set1_pd(-0.5);
const __m512d CN03                          = _mm512_set1_pd(-0.333333333333333333333333333);
const __m512d C025                          = _mm512_set1_pd(0.25);
const __m512d xx                            = _mm512_mul_pd(x,x);
const __m512d left_exp_val                  = _mm512_exp_pd(_mm512_mul_pd(CN05,xx));
const __m512d right_exp_val                 = _mm512_mul_pd(C025,_mm512_exp_pd(_mm512_mul_pd(CN03,_mm512_add_pd(xx,xx))));
return (_mm512_add_pd(left_exp_val,right_exp_val));
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
__m128
gaussian_Q_approx_loskot_2T_4xf32(const __m128 x)
{
const __m128 xx     = _mm_mul_ps(x,x);
const __m128 C0208  = _mm_set1_ps(0.208f);
const __m128 C0971  = _mm_mul_ps(_mm_set1_ps(-0.971f),xx);
const __m128 C0147  = _mm_set1_ps(0.147f);
const __m128 C0525  = _mm_mul_ps(_mm_set1_ps(-0.525f),xx);
const __m128 term   = _mm_mul_ps(C0147,simd_fast_exp_approx_4xf32(C0525));
return _mm_fmadd_ps(C0208,simd_fast_exp_approx_4xf32(C0971),term);
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
gaussian_Q_approx_loskot_2T_2xf64(const __m128d x)
{
const __m128d xx     = _mm_mul_pd(x,x);
const __m128d C0208  = _mm_set1_pd(0.208);
const __m128d C0971  = _mm_mul_pd(_mm_set1_pd(-0.971),xx);
const __m128d C0147  = _mm_set1_pd(0.147);
const __m128d C0525  = _mm_mul_ps(_mm_set1_pd(-0.525),xx);
const __m128d term   = _mm_mul_pd(C0147,simd_fast_exp_approx_2xf64(C0525));
return _mm_fmadd_pd(C0208,simd_fast_exp_approx_2xf64(C0971),term);
}

#if (GAUSSIAN_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX2
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
__m256 
gaussian_Q_approx_loskot_2T_8xf32(const __m256 x) 
{
const __m256 xx     = _mm256_mul_ps(x,x);
const __m256 C0208  = _mm256_set1_ps(0.208f);
const __m256 C0971  = _mm256_mul_ps(_mm256_set1_ps(-0.971f),xx);
const __m256 C0147  = _mm256_set1_ps(0.147f);
const __m256 C0525  = _mm256_mul_ps(_mm256_set1_ps(-0.525f),xx);
const __m256 term   = _mm256_mul_ps(C0147,mm256_exp_ps(C0525));
return _mm256_fmadd_ps(C0208,mm256_exp_ps(C0971),term);
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
gaussian_Q_approx_loskot_2T_4xf64(const __m256d x) 
{
const __m256d xx     = _mm256_mul_pd(x,x);
const __m256d C0208  = _mm256_set1_pd(0.208);
const __m256d C0971  = _mm256_mul_pd(_mm256_set1_ps(-0.971),xx);
const __m256d C0147  = _mm256_set1_pd(0.147);
const __m256d C0525  = _mm256_mul_pd(_mm256_set1_pd(-0.525),xx);
const __m256d term   = _mm256_mul_pd(C0147,mm256_exp_pd(C0525));
return _mm256_fmadd_pd(C0208,mm256_exp_pd(C0971),term);
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
gaussian_Q_approx_loskot_2T_16xf32(const __m512 x) 
{
const __m512 xx     = _mm512_mul_ps(x,x);
const __m512 C0208  = _mm512_set1_ps(0.208f);
const __m512 C0971  = _mm512_mul_ps(_mm512_set1_ps(-0.971f),xx);
const __m512 C0147  = _mm512_set1_ps(0.147f);
const __m512 C0525  = _mm512_mul_ps(_mm512_set1_ps(-0.525f),xx);
const __m512 term   = _mm512_mul_ps(C0147,mm512_exp_ps(C0525));
return _mm512_fmadd_ps(C0208,mm512_exp_ps(C0971),term);
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
__m512d 
gaussian_Q_approx_loskot_2T_8xf64(const __m512d x) 
{
const __m512d xx     = _mm512_mul_pd(x,x);
const __m512d C0208  = _mm512_set1_pd(0.208);
const __m512d C0971  = _mm512_mul_pd(_mm512_set1_pd(-0.971),xx);
const __m512d C0147  = _mm512_set1_pd(0.147);
const __m512d C0525  = _mm512_mul_pd(_mm512_set1_pd(-0.525),xx);
const __m512d term   = _mm512_mul_pd(C0147,mm512_exp_pd(C0525));
return _mm512_fmadd_pd(C0208,mm512_exp_pd(C0971),term);
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
__m128
gaussian_Q_approx_loskot_3T_4xf32(const __m128 x)
{
const __m128 xx    = _mm_mul_ps(x,x);
const __m128 C0168 = _mm_set1_ps(0.168f);
const __m128 C0876 = _mm_mul_ps(_mm_set1_ps(-0.876f),xx);
const __m128 C0144 = _mm_set1_ps(0.144f);
const __m128 C0525 = _mm_mul_ps(_mm_set1_ps(-0.525f),xx);
const __m128 C0002 = _mm_set1_ps(0.002f);
const __m128 C0603 = _mm_mul_ps(_mm_set1_ps(-0.603f),xx);
return (_mm_fmadd_ps(C0168,simd_fast_exp_approx_4xf32(C0876),
                     _mm_fmadd_ps(C0144,simd_fast_exp_approx_4xf32(C0525),
                                 _mm_mul_ps(C0002,simd_fast_exp_approx_4xf32(C0603)))));
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
gaussian_Q_approx_loskot_3T_2xf64(const __m128d x)
{
const __m128d xx    = _mm_mul_pd(x,x);
const __m128d C0168 = _mm_set1_pd(0.168);
const __m128d C0876 = _mm_mul_pd(_mm_set1_pd(-0.876),xx);
const __m128d C0144 = _mm_set1_pd(0.144);
const __m128d C0525 = _mm_mul_pd(_mm_set1_pd(-0.525),xx);
const __m128d C0002 = _mm_set1_pd(0.002);
const __m128d C0603 = _mm_mul_pd(_mm_set1_pd(-0.603),xx);
return (_mm_fmadd_pd(C0168,simd_fast_exp_approx_2xf64(C0876),
                     _mm_fmadd_pd(C0144,simd_fast_exp_approx_2xf64(C0525),
                                 _mm_mul_pd(C0002,simd_fast_exp_approx_2xf64(C0603)))));
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
gaussian_Q_approx_loskot_3T_8xf32(const __m256 x)
{
const __m256 xx    = _mm256_mul_ps(x,x);
const __m256 C0168 = _mm256_set1_ps(0.168f);
const __m256 C0876 = _mm256_mul_ps(_mm256_set1_ps(-0.876f),xx);
const __m256 C0144 = _mm256_set1_ps(0.144f);
const __m256 C0525 = _mm256_mul_ps(_mm256_set1_ps(-0.525f),xx);
const __m256 C0002 = _mm256_set1_ps(0.002f);
const __m256 C0603 = _mm256_mul_ps(_mm256_set1_ps(-0.603f),xx);
return (_mm256_fmadd_ps(C0168,_mm256_exp_ps(C0876),
                     _mm256_fmadd_ps(C0144,_mm256_exp_ps(C0525),
                                 _mm256_mul_ps(C0002,_mm256_exp_ps(C0603)))));
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
gaussian_Q_approx_loskot_3T_4xf64(const __m128d x)
{
const __m256d xx    = _mm256_mul_pd(x,x);
const __m256d C0168 = _mm256_set1_pd(0.168);
const __m256d C0876 = _mm256_mul_pd(_mm256_set1_pd(-0.876),xx);
const __m256d C0144 = _mm256_set1_pd(0.144);
const __m256d C0525 = _mm256_mul_pd(_mm256_set1_pd(-0.525),xx);
const __m256d C0002 = _mm256_set1_pd(0.002);
const __m256d C0603 = _mm256_mul_pd(_mm256_set1_pd(-0.603),xx);
return (_mm256_fmadd_pd(C0168,_mm256_exp_pd(C0876),
                     _mm256_fmadd_pd(C0144,_mm256_exp_pd(C0525),
                                 _mm256_mul_pd(C0002,_mm256_exp_pd(C0603)))));
}

#if (GAUSSIAN_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
__m512
gaussian_Q_approx_loskot_3T_16xf32(const __m512 x)
{
const __m512 xx    = _mm512_mul_ps(x,x);
const __m512 C0168 = _mm512_set1_ps(0.168f);
const __m512 C0876 = _mm512_mul_ps(_mm512_set1_ps(-0.876f),xx);
const __m512 C0144 = _mm512_set1_ps(0.144f);
const __m512 C0525 = _mm512_mul_ps(_mm512_set1_ps(-0.525f),xx);
const __m512 C0002 = _mm512_set1_ps(0.002f);
const __m512 C0603 = _mm512_mul_ps(_mm512_set1_ps(-0.603f),xx);
return (_mm512_fmadd_ps(C0168,_mm512_exp_ps(C0876),
                     _mm512_fmadd_ps(C0144,_mm512_exp_ps(C0525),
                                 _mm512_mul_ps(C0002,_mm512_exp_ps(C0603)))));
}

#if (GAUSSIAN_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
__m512d
gaussian_Q_approx_loskot_3T_8xf32(const __m512d x)
{
const __m512d xx    = _mm512_mul_pd(x,x);
const __m512d C0168 = _mm512_set1_pd(0.168);
const __m512d C0876 = _mm512_mul_pd(_mm512_set1_pd(-0.876),xx);
const __m512d C0144 = _mm512_set1_pd(0.144);
const __m512d C0525 = _mm512_mul_pd(_mm512_set1_pd(-0.525),xx);
const __m512d C0002 = _mm512_set1_pd(0.002);
const __m512d C0603 = _mm512_mul_pd(_mm512_set1_pd(-0.603),xx);
return (_mm512_fmadd_pd(C0168,_mm512_exp_pd(C0876),
                     _mm512_fmadd_pd(C0144,_mm512_exp_pd(C0525),
                                 _mm512_mul_pd(C0002,_mm512_exp_pd(C0603)))));
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
__m128
gaussian_Q_approx_sadhwani_1T_4xf32(const __m128 x)
{
const __m128 xx    = _mm_mul_ps(x,x);
const __m128 C025  = _mm_set1_ps(0.25f);
const __m128 C05   = _mm_mul_ps(_mm_set1_ps(-0.5f),xx);
return (_mm_mul_ps(C025,simd_fast_exp_approx_4xf32(C05)));
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
gaussian_Q_approx_sadhwani_1T_2xf64(const __m128d x)
{
const __m128d xx    = _mm_mul_pd(x,x);
const __m128d C025  = _mm_set1_pd(0.25);
const __m128d C05   = _mm_mul_pd(_mm_set1_ps(-0.5),xx);
return (_mm_mul_pd(C025,simd_fast_exp_approx_2xf64(C05)));
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
gaussian_Q_approx_sadhwani_1T_8xf32(const __m256 x)
{
const __m256 xx    = _mm256_mul_ps(x,x);
const __m256 C025  = _mm256_set1_ps(0.25f);
const __m256 C05   = _mm256_mul_ps(_mm256_set1_ps(-0.5f),xx);
return (_mm256_mul_ps(C025,_mm256_exp_ps(C05)));
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
gaussian_Q_approx_sadhwani_1T_4xf64(const __m256d x)
{
const __m256d xx    = _mm256_mul_pd(x,x);
const __m256d C025  = _mm256_set1_pd(0.25);
const __m256d C05   = _mm256_mul_pd(_mm256_set1_ps(-0.5),xx);
return (_mm256_mul_pd(C025,_mm256_exp_pd(C05)));
}

#if (GAUSSIAN_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
__m512
gaussian_Q_approx_sadhwani_1T_16xf32(const __m512 x)
{
const __m512 xx    = _mm512_mul_ps(x,x);
const __m512 C025  = _mm512_set1_ps(0.25f);
const __m512 C05   = _mm512_mul_ps(_mm512_set1_ps(-0.5f),xx);
return (_mm512_mul_ps(C025,_mm512_exp_ps(C05)));
}

#if (GAUSSIAN_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
__m512d
gaussian_Q_approx_sadhwani_1T_8xf64(const __m512d x)
{
const __m512d xx    = _mm512_mul_pd(x,x);
const __m512d C025  = _mm512_set1_pd(0.25);
const __m512d C05   = _mm512_mul_pd(_mm512_set1_ps(-0.5),xx);
return (_mm512_mul_pd(C025,_mm512_exp_pd(C05)));
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
__m128 
gaussian_Q_approx_sadhwani_2T_4xf32(const __m128 x)
{
const __m128 xx    = _mm_mul_ps(x,x);
const __m128 C0125 = _mm_set1_ps(0.125f);
const __m128 C025  = _mm_set1_ps(0.25f);
const __m128 negxx = _mm_sub_ps(_mm_setzero_ps(),xx);
const __m128 C05   = _mm_mul_ps(_mm_set1_ps(-0.5f),xx);
return (_mm_fmadd_ps(C0125,simd_fast_exp_approx_4xf32(negxx),
                           _mm_mul_ps(C025,simd_fast_exp_approx_4xf32(C05))));
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
gaussian_Q_approx_sadhwani_2T_2xf64(const __m128d x)
{
const __m128d xx    = _mm_mul_pd(x,x);
const __m128d C0125 = _mm_set1_pd(0.125);
const __m128d C025  = _mm_set1_pd(0.25);
const __m128d negxx = _mm_sub_pd(_mm_setzero_pd(),xx);
const __m128d C05   = _mm_mul_pd(_mm_set1_pd(-0.5),xx);
return (_mm_fmadd_pd(C0125,simd_fast_exp_approx_2xf64(negxx),
                           _mm_mul_pd(C025,simd_fast_exp_approx_2xf64(C05))));
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
gaussian_Q_approx_sadhwani_2T_8xf32(const __m256 x)
{
const __m256 xx    = _mm256_mul_ps(x,x);
const __m256 C0125 = _mm256_set1_ps(0.125f);
const __m256 C025  = _mm256_set1_ps(0.25f);
const __m256 negxx = _mm256_sub_ps(_mm256_setzero_ps(),xx);
const __m256 C05   = _mm256_mul_ps(_mm256_set1_ps(-0.5f),xx);
return (_mm256_fmadd_ps(C0125,_mm256_exp_ps(negxx),
                           _mm256_mul_ps(C025,_mm256_exp_ps(C05))));
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
gaussian_Q_approx_sadhwani_2T_4xf64(const __m256d x)
{
const __m256d xx    = _mm256_mul_pd(x,x);
const __m256d C0125 = _mm256_set1_pd(0.125);
const __m256d C025  = _mm256_set1_pd(0.25);
const __m256d negxx = _mm256_sub_pd(_mm256_setzero_pd(),xx);
const __m256d C05   = _mm256_mul_pd(_mm256_set1_pd(-0.5),xx);
return (_mm256_fmadd_pd(C0125,_mm256_exp_pd(negxx),
                           _mm256_mul_pd(C025,_mm256_exp_pd(C05))));
}

#if (GAUSSIAN_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
__m512
gaussian_Q_approx_sadhwani_2T_16xf32(const __m512 x)
{
const __m512 xx    = _mm512_mul_ps(x,x);
const __m512 C0125 = _mm512_set1_ps(0.125f);
const __m512 C025  = _mm512_set1_ps(0.25f);
const __m512 negxx = _mm512_sub_ps(_mm512_setzero_ps(),xx);
const __m512 C05   = _mm512_mul_ps(_mm512_set1_ps(-0.5f),xx);
return (_mm512_fmadd_ps(C0125,_mm512_exp_ps(negxx),
                           _mm512_mul_ps(C025,_mm512_exp_ps(C05))));
}

#if (GAUSSIAN_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
__m512d
gaussian_Q_approx_sadhwani_2T_8xf64(const __m512d x)
{
const __m512d xx    = _mm512_mul_pd(x,x);
const __m512d C0125 = _mm512_set1_pd(0.125);
const __m512d C025  = _mm512_set1_pd(0.25);
const __m512d negxx = _mm512_sub_pd(_mm512_setzero_pd(),xx);
const __m512d C05   = _mm512_mul_pd(_mm512_set1_pd(-0.5),xx);
return (_mm512_fmadd_pd(C0125,_mm512_exp_pd(negxx),
                           _mm512_mul_pd(C025,_mm512_exp_pd(C05))));
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
__m128 
gaussian_Q_approx_sadhwani_4T_4xf32(const __m128 x)
{
const __m128 xx  = _mm_mul_ps(x,x);
const __m128 C1  = _mm_set1_ps(0.0625f);
const __m128 C2  = _mm_set1_ps(0.125f);
const __m128 negxx = _mm_sub_ps(_mm_setzero_ps(),xx);
const __m128 C3  = _mm_mul_ps(_mm_set1_ps(-0.5f),xx);
const __m128 C4  = _mm_mul_ps(_mm_set1_ps(-0.33333333333333333333333f),xx);
const __m128 C5  = _mm_mul_ps(_mm_set1_ps(-0.058823529411764705882352941f),xx);
const __m128 tenxx = _mm_mul_ps(_mm_set1_ps(10.0f),xx);
return (_mm_fmadd_ps(C1,simd_fast_exp_approx_4xf32(C3),
                    _mm_fmadd_ps(C2,simd_fast_exp_approx_4xf32(negxx),
                                 _mm_fmadd_ps(C2,simd_fast_exp_approx_4xf32(C4),
                                                _mm_mul_ps(C2,simd_fast_exp_approx_4xf32(C5))))));
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
gaussian_Q_approx_sadhwani_4T_2xf64(const __m128d x)
{
const __m128d xx  = _mm_mul_pd(x,x);
const __m128d C1  = _mm_set1_pd(0.0625);
const __m128d C2  = _mm_set1_ps(0.125);
const __m128d negxx = _mm_sub_pd(_mm_setzero_pd(),xx);
const __m128d C3  = _mm_mul_pd(_mm_set1_pd(-0.5),xx);
const __m128d C4  = _mm_mul_pd(_mm_set1_pd(-0.33333333333333333333333),xx);
const __m128d C5  = _mm_mul_pd(_mm_set1_pd(-0.058823529411764705882352941),xx);
const __m128d tenxx = _mm_mul_pd(_mm_set1_pd(10.0),xx);
return (_mm_fmadd_pd(C1,simd_fast_exp_approx_2xf64(C3),
                    _mm_fmadd_pd(C2,simd_fast_exp_approx_2xf64(negxx),
                                 _mm_fmadd_pd(C2,simd_fast_exp_approx_2xf64(C4),
                                                _mm_mul_pd(C2,simd_fast_exp_approx_2xf64(C5))))));
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
gaussian_Q_approx_sadhwani_4T_8xf32(const __m256 x)
{
const __m256 xx  = _mm256_mul_ps(x,x);
const __m256 C1  = _mm256_set1_ps(0.0625f);
const __m256 C2  = _mm256_set1_ps(0.125f);
const __m256 negxx = _mm256_sub_ps(_mm256_setzero_ps(),xx);
const __m256 C3  = _mm256_mul_ps(_mm256_set1_ps(-0.5f),xx);
const __m256 C4  = _mm256_mul_ps(_mm256_set1_ps(-0.33333333333333333333333f),xx);
const __m256 C5  = _mm256_mul_ps(_mm256_set1_ps(-0.058823529411764705882352941f),xx);
const __m256 tenxx = _mm256_mul_ps(_mm256_set1_ps(10.0f),xx);
return (_mm256_fmadd_ps(C1,_mm256_exp_ps(C3),
                    _mm256_fmadd_ps(C2,_mm256_exp_ps(negxx),
                                 _mm256_fmadd_ps(C2,_mm256_exp_ps(C4),
                                                _mm256_mul_ps(C2,_mm256_exp_ps(C5))))));
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
gaussian_Q_approx_sadhwani_4T_4xf64(const __m256d x)
{
const __m256d xx  = _mm256_mul_pd(x,x);
const __m256d C1  = _mm256_set1_pd(0.0625);
const __m256d C2  = _mm256_set1_pd(0.125);
const __m256d negxx = _mm256_sub_pd(_mm256_setzero_pd(),xx);
const __m256d C3  = _mm256_mul_pd(_mm256_set1_pd(-0.5),xx);
const __m256d C4  = _mm256_mul_pd(_mm256_set1_pd(-0.33333333333333333333333),xx);
const __m256d C5  = _mm256_mul_pd(_mm256_set1_pd(-0.058823529411764705882352941),xx);
const __m256d tenxx = _mm256_mul_pd(_mm256_set1_pd(10.0),xx);
return (_mm256_fmadd_pd(C1,_mm256_exp_pd(C3),
                    _mm256_fmadd_pd(C2,_mm256_exp_pd(negxx),
                                 _mm256_fmadd_pd(C2,_mm256_exp_pd(C4),
                                                _mm256_mul_pd(C2,_mm256_exp_pd(C5))))));
}

#if (GAUSSIAN_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
__m512
gaussian_Q_approx_sadhwani_4T_16xf32(const __m512 x)
{
const __m512 xx  = _mm512_mul_ps(x,x);
const __m512 C1  = _mm512_set1_ps(0.0625f);
const __m512 C2  = _mm512_set1_ps(0.125f);
const __m512 negxx = _mm512_sub_ps(_mm512_setzero_ps(),xx);
const __m512 C3  = _mm512_mul_ps(_mm512_set1_ps(-0.5f),xx);
const __m512 C4  = _mm512_mul_ps(_mm512_set1_ps(-0.33333333333333333333333f),xx);
const __m512 C5  = _mm512_mul_ps(_mm512_set1_ps(-0.058823529411764705882352941f),xx);
const __m512 tenxx = _mm512_mul_ps(_mm512_set1_ps(10.0f),xx);
return (_mm512_fmadd_ps(C1,_mm512_exp_ps(C3),
                    _mm512_fmadd_ps(C2,_mm512_exp_ps(negxx),
                                 _mm512_fmadd_ps(C2,_mm512_exp_ps(C4),
                                                _mm512_mul_ps(C2,_mm512_exp_ps(C5))))));
}

#if (GAUSSIAN_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
__m512d
gaussian_Q_approx_sadhwani_4T_8xf64(const __m512d x)
{
const __m512d xx  = _mm512_mul_pd(x,x);
const __m512d C1  = _mm512_set1_pd(0.0625);
const __m512d C2  = _mm512_set1_pd(0.125);
const __m512d negxx = _mm512_sub_pd(_mm512_setzero_pd(),xx);
const __m512d C3  = _mm512_mul_pd(_mm512_set1_pd(-0.5),xx);
const __m512d C4  = _mm512_mul_pd(_mm512_set1_pd(-0.33333333333333333333333),xx);
const __m512d C5  = _mm512_mul_pd(_mm512_set1_pd(-0.058823529411764705882352941),xx);
const __m512d tenxx = _mm512_mul_pd(_mm512_set1_pd(10.0f),xx);
return (_mm512_fmadd_pd(C1,_mm512_exp_pd(C3),
                    _mm512_fmadd_pd(C2,_mm512_exp_pd(negxx),
                                 _mm512_fmadd_pd(C2,_mm512_exp_pd(C4),
                                                _mm512_mul_pd(C2,_mm512_exp_pd(C5))))));
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
__m128 
gaussian_Q_approx_borjesson_4xf32(const __m128 x)
{
const __m128 C0339     = _mm_set1_ps(0.339f);
const __m128 xx        = _mm_mul_ps(x,x);
const __m128 C5510     = _mm_set1_ps(5.510f);
const __m128 C1        = _mm_mul_ps(_mm_set1_ps(-0.5f),xx);
const __m128 C0661     = _mm_set1_ps(0.661f);
const __m128 inv2PI    = _mm_set1_ps(0.3989422804014326779399460599344f);
const __m128 exp_val   = _mm_mul_ps(inv2PI,simd_fast_exp_approx_4xf32(C1));
const __m128 one       = _mm_set1_ps(1.0f);
const __m128 sqrt_term = _mm_mul_ps(C0339,_mm_sqrt_ps(_mm_mul_ps(xx,C5510)));
const __m128 right_term= _mm_rcp_ps(_mm_fmadd_ps(C0661,x,sqrt_term));
return (_mm_mul_ps(exp_val,right_term));
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
gaussian_Q_approx_borjesson_2xf64(const __m128d x)
{
const __m128d C0339     = _mm_set1_pd(0.339);
const __m128d xx        = _mm_mul_pd(x,x);
const __m128d C5510     = _mm_set1_pd(5.510);
const __m128d C1        = _mm_mul_pd(_mm_set1_pd(-0.5),xx);
const __m128d C0661     = _mm_set1_pd(0.661);
const __m128d inv2PI    = _mm_set1_pd(0.3989422804014326779399460599344);
const __m128d exp_val   = _mm_mul_pd(inv2PI,simd_fast_exp_approx_2xf64(C1));
const __m128d one       = _mm_set1_pd(1.0);
const __m128d sqrt_term = _mm_mul_pd(C0339,_mm_sqrt_pd(_mm_mul_pd(xx,C5510)));
#if (GAUSSIAN_Q_APPROX_SIMD_RECIPROCAL_BY_DIVISION) == 1
const __m128d right_term= _mm_div_pd(one,_mm_fmadd_pd(C0661,x,sqrt_term));
#else 
const __m128d right_term= _mm_rcp14_pd(_mm_fmadd_pd(C0661,x,sqrt_term));
#endif 
return (_mm_mul_pd(exp_val,right_term));
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
gaussian_Q_approx_borjesson_8xf32(const __m256 x)
{
const __m256 C0339     = _mm256_set1_ps(0.339f);
const __m256 xx        = _mm256_mul_ps(x,x);
const __m256 C5510     = _mm256_set1_ps(5.510f);
const __m256 C1        = _mm256_mul_ps(_mm256_set1_ps(-0.5f),xx);
const __m256 C0661     = _mm256_set1_ps(0.661f);
const __m256 inv2PI    = _mm256_set1_ps(0.3989422804014326779399460599344f);
const __m256 exp_val   = _mm256_mul_ps(inv2PI,_mm256_exp_ps(C1));
const __m256 one       = _mm256_set1_ps(1.0f);
const __m256 sqrt_term = _mm256_mul_ps(C0339,_mm256_sqrt_ps(_mm256_mul_ps(xx,C5510)));
const __m256 right_term= _mm256_rcp_ps(_mm256_fmadd_ps(C0661,x,sqrt_term));
return (_mm256_mul_ps(exp_val,right_term));
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
gaussian_Q_approx_borjesson_4xf64(const __m256d x)
{
const __m256d C0339     = _mm256_set1_pd(0.339);
const __m256d xx        = _mm256_mul_pd(x,x);
const __m256d C5510     = _mm256_set1_pd(5.510);
const __m256d C1        = _mm256_mul_pd(_mm256_set1_pd(-0.5),xx);
const __m256d C0661     = _mm256_set1_pd(0.661);
const __m256d inv2PI    = _mm256_set1_pd(0.3989422804014326779399460599344);
const __m256d exp_val   = _mm256_mul_pd(inv2PI,_mm256_exp_pd(C1));
const __m256d one       = _mm256_set1_pd(1.0f);
const __m256d sqrt_term = _mm256_mul_pd(C0339,_mm256_sqrt_pd(_mm256_mul_pd(xx,C5510)));
#if (GAUSSIAN_Q_APPROX_SIMD_RECIPROCAL_BY_DIVISION) == 1
const __m256d right_term= _mm256_div_pd(one,_mm256_fmadd_pd(C0661,x,sqrt_term));
#else 
const __m256d right_term= _mm256_rcp14_pd(_mm256_fmadd_pd(C0661,x,sqrt_term));
#endif 
return (_mm256_mul_pd(exp_val,right_term));
}

#if (GAUSSIAN_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
__m512
gaussian_Q_approx_borjesson_16xf32(const __m512 x)
{
const __m512 C0339     = _mm512_set1_ps(0.339f);
const __m512 xx        = _mm512_mul_ps(x,x);
const __m512 C5510     = _mm512_set1_ps(5.510f);
const __m512 C1        = _mm512_mul_ps(_mm512_set1_ps(-0.5f),xx);
const __m512 C0661     = _mm512_set1_ps(0.661f);
const __m512 inv2PI    = _mm512_set1_ps(0.3989422804014326779399460599344f);
const __m512 exp_val   = _mm512_mul_ps(inv2PI,_mm512_exp_ps(C1));
const __m512 one       = _mm512_set1_ps(1.0f);
const __m512 sqrt_term = _mm512_mul_ps(C0339,_mm512_sqrt_ps(_mm512_mul_ps(xx,C5510)));
const __m512 right_term= _mm512_rcp_ps(_mm512_fmadd_ps(C0661,x,sqrt_term));
return (_mm512_mul_ps(exp_val,right_term));
}

#if (GAUSSIAN_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
__m512d
gaussian_Q_approx_borjesson_8xf64(const __m512d x)
{
const __m512d C0339     = _mm512_set1_pd(0.339);
const __m512d xx        = _mm512_mul_pd(x,x);
const __m512d C5510     = _mm512_set1_pd(5.510);
const __m512d C1        = _mm512_mul_pd(_mm512_set1_pd(-0.5),xx);
const __m512d C0661     = _mm512_set1_pd(0.661);
const __m512d inv2PI    = _mm512_set1_pd(0.3989422804014326779399460599344f);
const __m512d exp_val   = _mm512_mul_pd(inv2PI,_mm512_exp_pd(C1));
const __m512d one       = _mm512_set1_pd(1.0f);
const __m512d sqrt_term = _mm512_mul_pd(C0339,_mm512_sqrt_pd(_mm512_mul_pd(xx,C5510)));
#if (GAUSSIAN_Q_APPROX_SIMD_RECIPROCAL_BY_DIVISION) == 1
const __m512d right_term= _mm512_div_pd(one,_mm512_fmadd_pd(C0661,x,sqrt_term));
#else 
const __m512d right_term= _mm512_rcp14_pd(_mm512_fmadd_pd(C0661,x,sqrt_term));
#endif 
return (_mm512_mul_pd(exp_val,right_term));
}

/* Currently unused*/
#if !defined(GAUSSIAN_Q_APPROX_SIMD_PREVENT_BROADCAST_INSERTION)
#define GAUSSIAN_Q_APPROX_SIMD_PREVENT_BROADCAST_INSERTION 0
#endif 

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
gaussian_Q_approx_sadhwani_summed_4xf32(const __m128 x,const std::int32_t n) 
{
if (__builtin_expect(n>128,0)) {return (_mm_set1_ps(-1.0));}
__ATTR_ALIGN__(16) 
const __m128 theta_lut[128] = 
{
_mm_set1_ps(26561.07370058031301596),
_mm_set1_ps(2951.52672978394593883),
_mm_set1_ps(1062.76301068144402961),
_mm_set1_ps(542.38938216885594557),
_mm_set1_ps(328.24391161267163852),
_mm_set1_ps(219.84388820771997075),
_mm_set1_ps(157.49790211713116150),
_mm_set1_ps(118.38163379128440056),
_mm_set1_ps(92.23973597198599350),
_mm_set1_ps(73.90969946516314337),
_mm_set1_ps(60.56287935539454281),
_mm_set1_ps(50.54400438065793111),
_mm_set1_ps(42.83209252697486136),
_mm_set1_ps(36.76965720093080847),
_mm_set1_ps(31.91778586717441257),
_mm_set1_ps(27.97440695989148551),
_mm_set1_ps(24.72611246534243534),
_mm_set1_ps(22.01866773467314076),
_mm_set1_ps(19.73836125017757581),
_mm_set1_ps(17.79986866518916599),
_mm_set1_ps(16.13815476591664932),
_mm_set1_ps(14.70294783067476274),
_mm_set1_ps(13.45489318988675542),
_mm_set1_ps(12.36282713537543643),
_mm_set1_ps(11.40181318331400817),
_mm_set1_ps(10.55170643021994792),
_mm_set1_ps(9.79608972017488178),
_mm_set1_ps(9.12147551075540264),
_mm_set1_ps(8.51670021771579755),
_mm_set1_ps(7.97245976035782089),
_mm_set1_ps(7.48094990150482975),
_mm_set1_ps(7.03558520545476096),
_mm_set1_ps(6.63077757000577961),
_mm_set1_ps(6.26176032547270278),
_mm_set1_ps(5.92444749260821357),
_mm_set1_ps(5.61532039144371531),
_mm_set1_ps(5.33133569098266236),
_mm_set1_ps(5.06985038851206493),
_mm_set1_ps(4.82856024770480818),
_mm_set1_ps(4.60544900516174227),
_mm_set1_ps(4.39874624527744107),
_mm_set1_ps(4.20689229309963242),
_mm_set1_ps(4.02850882009587341),
_mm_set1_ps(3.86237412456220319),
_mm_set1_ps(3.70740225596646766),
_mm_set1_ps(3.56262531497762280),
_mm_set1_ps(3.42717838884132009),
_mm_set1_ps(3.30028668303544404),
_mm_set1_ps(3.18125449075313638),
_mm_set1_ps(3.06945570625878084),
_mm_set1_ps(2.96432564001850762),
_mm_set1_ps(2.86535393539436001),
_mm_set1_ps(2.77207842067833932),
_mm_set1_ps(2.68407975793799203),
_mm_set1_ps(2.60097677280667572),
_mm_set1_ps(2.52242236796785813),
_mm_set1_ps(2.44809993843393503),
_mm_set1_ps(2.37772021942528200),
_mm_set1_ps(2.31101850820766996),
_mm_set1_ps(2.24775221004035153),
_mm_set1_ps(2.18769866574001037),
_mm_set1_ps(2.13065322453304118),
_mm_set1_ps(2.07642753105719757),
_mm_set1_ps(2.02484799975171637),
_mm_set1_ps(1.97575445357971491),
_mm_set1_ps(1.92899890717015055),
_mm_set1_ps(1.88444447714107866),
_mm_set1_ps(1.84196440464714883),
_mm_set1_ps(1.80144117714498853),
_mm_set1_ps(1.76276573804219461),
_mm_set1_ps(1.72583677433227445),
_mm_set1_ps(1.69056007355504945),
_mm_set1_ps(1.65684794248971023),
_mm_set1_ps(1.62461868091121531),
_mm_set1_ps(1.59379610454104137),
_mm_set1_ps(1.56430911201842116),
_mm_set1_ps(1.53609129132299893),
_mm_set1_ps(1.50908056160717741),
_mm_set1_ps(1.48321884685699712),
_mm_set1_ps(1.45845177820343341),
_mm_set1_ps(1.43472842205929152),
_mm_set1_ps(1.41200103156710211),
_mm_set1_ps(1.39022481911630291),
_mm_set1_ps(1.36935774792838338),
_mm_set1_ps(1.34936034092082524),
_mm_set1_ps(1.33019550524813823),
_mm_set1_ps(1.31182837108428330),
_mm_set1_ps(1.29422614335786368),
_mm_set1_ps(1.27735796528209100),
_mm_set1_ps(1.26119479263764811),
_mm_set1_ps(1.24570927786994146),
_mm_set1_ps(1.23087566315441310),
_mm_set1_ps(1.21666968166584666),
_mm_set1_ps(1.20306846636115816),
_mm_set1_ps(1.19005046565099137),
_mm_set1_ps(1.17759536539444132),
_mm_set1_ps(1.16568401670416888),
_mm_set1_ps(1.15429836909675965),
_mm_set1_ps(1.14342140856596108),
_mm_set1_ps(1.13303710019499548),
_mm_set1_ps(1.12313033495892256),
_mm_set1_ps(1.11368688039941000),
_mm_set1_ps(1.10469333488267329),
_mm_set1_ps(1.09613708517704822),
_mm_set1_ps(1.08800626710994641),
_mm_set1_ps(1.08028972908511101),
_mm_set1_ps(1.07297699826029014),
_mm_set1_ps(1.06605824920296688),
_mm_set1_ps(1.05952427485770095),
_mm_set1_ps(1.05336645967322906),
_mm_set1_ps(1.04757675475074774),
_mm_set1_ps(1.04214765488702521),
_mm_set1_ps(1.03707217739714941),
_mm_set1_ps(1.03234384261201928),
_mm_set1_ps(1.02795665595514785),
_mm_set1_ps(1.02390509151210085),
_mm_set1_ps(1.02018407701397673),
_mm_set1_ps(1.01678898016385832),
_mm_set1_ps(1.01371559624214091),
_mm_set1_ps(1.01096013693317843),
_mm_set1_ps(1.00851922032179320),
_mm_set1_ps(1.00638986201394287),
_mm_set1_ps(1.00456946734126218),
_mm_set1_ps(1.00305582461434328),
_mm_set1_ps(1.00184709939451322),
_mm_set1_ps(1.00094182975856194),
_mm_set1_ps(1.00033892253539203),
_mm_set1_ps(1.00003765049793447)
};
const __m128 inv2n = _mm_rcp_ps(_mm_set1_ps(static_cast<float>(n+n)));
const __m128 xx    = _mm_mul_ps(x,x);
const __m128 C1    = _mm_set1_ps(-0.5f);
__m128 sum = _mm_setzero_ps();
for(std::int32_t j = 0; j<n; ++j)  
{
    __m128 xmm0    = theta_lut[j];
    __m128 exp_arg = _mm_mul_ps(C1,_mm_mul_ps(xx,xmm0));
    __m128 exp_val = simd_fast_exp_approx_4xf32(exp_arg);
    sum            = _mm_add_ps(sum,exp_val); 
}
return (_mm_mul_ps(inv2n,sum));
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
gaussian_Q_approx_sadhwani_summed_2xf64(const __m128d x,const std::int32_t n) 
{
if (__builtin_expect(n>128,0)) {return (_mm_set1_pd(-1.0));}
__ATTR_ALIGN__(16) 
const __m128d theta_lut[128] = 
{
_mm_set1_pd(26561.07370058031301596),
_mm_set1_pd(2951.52672978394593883),
_mm_set1_pd(1062.76301068144402961),
_mm_set1_pd(542.38938216885594557),
_mm_set1_pd(328.24391161267163852),
_mm_set1_pd(219.84388820771997075),
_mm_set1_pd(157.49790211713116150),
_mm_set1_pd(118.38163379128440056),
_mm_set1_pd(92.23973597198599350),
_mm_set1_pd(73.90969946516314337),
_mm_set1_pd(60.56287935539454281),
_mm_set1_pd(50.54400438065793111),
_mm_set1_pd(42.83209252697486136),
_mm_set1_pd(36.76965720093080847),
_mm_set1_pd(31.91778586717441257),
_mm_set1_pd(27.97440695989148551),
_mm_set1_pd(24.72611246534243534),
_mm_set1_pd(22.01866773467314076),
_mm_set1_pd(19.73836125017757581),
_mm_set1_pd(17.79986866518916599),
_mm_set1_pd(16.13815476591664932),
_mm_set1_pd(14.70294783067476274),
_mm_set1_pd(13.45489318988675542),
_mm_set1_pd(12.36282713537543643),
_mm_set1_pd(11.40181318331400817),
_mm_set1_pd(10.55170643021994792),
_mm_set1_pd(9.79608972017488178),
_mm_set1_pd(9.12147551075540264),
_mm_set1_pd(8.51670021771579755),
_mm_set1_pd(7.97245976035782089),
_mm_set1_pd(7.48094990150482975),
_mm_set1_pd(7.03558520545476096),
_mm_set1_pd(6.63077757000577961),
_mm_set1_pd(6.26176032547270278),
_mm_set1_pd(5.92444749260821357),
_mm_set1_pd(5.61532039144371531),
_mm_set1_pd(5.33133569098266236),
_mm_set1_pd(5.06985038851206493),
_mm_set1_pd(4.82856024770480818),
_mm_set1_pd(4.60544900516174227),
_mm_set1_pd(4.39874624527744107),
_mm_set1_pd(4.20689229309963242),
_mm_set1_pd(4.02850882009587341),
_mm_set1_pd(3.86237412456220319),
_mm_set1_pd(3.70740225596646766),
_mm_set1_pd(3.56262531497762280),
_mm_set1_pd(3.42717838884132009),
_mm_set1_pd(3.30028668303544404),
_mm_set1_pd(3.18125449075313638),
_mm_set1_pd(3.06945570625878084),
_mm_set1_pd(2.96432564001850762),
_mm_set1_pd(2.86535393539436001),
_mm_set1_pd(2.77207842067833932),
_mm_set1_pd(2.68407975793799203),
_mm_set1_pd(2.60097677280667572),
_mm_set1_pd(2.52242236796785813),
_mm_set1_pd(2.44809993843393503),
_mm_set1_pd(2.37772021942528200),
_mm_set1_pd(2.31101850820766996),
_mm_set1_pd(2.24775221004035153),
_mm_set1_pd(2.18769866574001037),
_mm_set1_pd(2.13065322453304118),
_mm_set1_pd(2.07642753105719757),
_mm_set1_pd(2.02484799975171637),
_mm_set1_pd(1.97575445357971491),
_mm_set1_pd(1.92899890717015055),
_mm_set1_pd(1.88444447714107866),
_mm_set1_pd(1.84196440464714883),
_mm_set1_pd(1.80144117714498853),
_mm_set1_pd(1.76276573804219461),
_mm_set1_pd(1.72583677433227445),
_mm_set1_pd(1.69056007355504945),
_mm_set1_pd(1.65684794248971023),
_mm_set1_pd(1.62461868091121531),
_mm_set1_pd(1.59379610454104137),
_mm_set1_pd(1.56430911201842116),
_mm_set1_pd(1.53609129132299893),
_mm_set1_pd(1.50908056160717741),
_mm_set1_pd(1.48321884685699712),
_mm_set1_pd(1.45845177820343341),
_mm_set1_pd(1.43472842205929152),
_mm_set1_pd(1.41200103156710211),
_mm_set1_pd(1.39022481911630291),
_mm_set1_pd(1.36935774792838338),
_mm_set1_pd(1.34936034092082524),
_mm_set1_pd(1.33019550524813823),
_mm_set1_pd(1.31182837108428330),
_mm_set1_pd(1.29422614335786368),
_mm_set1_pd(1.27735796528209100),
_mm_set1_pd(1.26119479263764811),
_mm_set1_pd(1.24570927786994146),
_mm_set1_pd(1.23087566315441310),
_mm_set1_pd(1.21666968166584666),
_mm_set1_pd(1.20306846636115816),
_mm_set1_pd(1.19005046565099137),
_mm_set1_pd(1.17759536539444132),
_mm_set1_pd(1.16568401670416888),
_mm_set1_pd(1.15429836909675965),
_mm_set1_pd(1.14342140856596108),
_mm_set1_pd(1.13303710019499548),
_mm_set1_pd(1.12313033495892256),
_mm_set1_pd(1.11368688039941000),
_mm_set1_pd(1.10469333488267329),
_mm_set1_pd(1.09613708517704822),
_mm_set1_pd(1.08800626710994641),
_mm_set1_pd(1.08028972908511101),
_mm_set1_pd(1.07297699826029014),
_mm_set1_pd(1.06605824920296688),
_mm_set1_pd(1.05952427485770095),
_mm_set1_pd(1.05336645967322906),
_mm_set1_pd(1.04757675475074774),
_mm_set1_pd(1.04214765488702521),
_mm_set1_pd(1.03707217739714941),
_mm_set1_pd(1.03234384261201928),
_mm_set1_pd(1.02795665595514785),
_mm_set1_pd(1.02390509151210085),
_mm_set1_pd(1.02018407701397673),
_mm_set1_pd(1.01678898016385832),
_mm_set1_pd(1.01371559624214091),
_mm_set1_pd(1.01096013693317843),
_mm_set1_pd(1.00851922032179320),
_mm_set1_pd(1.00638986201394287),
_mm_set1_pd(1.00456946734126218),
_mm_set1_pd(1.00305582461434328),
_mm_set1_pd(1.00184709939451322),
_mm_set1_pd(1.00094182975856194),
_mm_set1_pd(1.00033892253539203),
_mm_set1_pd(1.00003765049793447)
};
#if (GAUSSIAN_Q_APPROX_SIMD_RECIPROCAL_BY_DIVISION) == 1
const __m128d one   = _mm_set1_pd(1.0);
const __m128d inv2n = _mm_div_pd(one,_mm_set1_pd(static_cast<double>(n+n)));
#else 
const __m128d inv2n = _mm_rcp14_pd(_mm_set1_pd(static_cast<double>(n+n)));
#endif 
const __m128d xx    = _mm_mul_pd(x,x);
const __m128d C1    = _mm_set1_pd(-0.5f);
__m128d sum = _mm_setzero_pd();
for(std::int32_t j = 0; j<n; ++j)  
{
    __m128d xmm0    = theta_lut[j];
    __m128d exp_arg = _mm_mul_pd(C1,_mm_mul_pd(xx,xmm0));
    __m128d exp_val = simd_fast_exp_approx_2xf64(exp_arg);
    sum            = _mm_add_pd(sum,exp_val); 
}
return (_mm_mul_pd(inv2n,sum));
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
gaussian_Q_approx_sadhwani_summed_8xf32(const __m256 x,const std::int32_t n) 
{
if (__builtin_expect(n>128,0)) {return (_mm256_set1_ps(-1.0));}
__ATTR_ALIGN__(32) 
const __m256 theta_lut[128] = 
{
_mm256_set1_ps(26561.07370058031301596),
_mm256_set1_ps(2951.52672978394593883),
_mm256_set1_ps(1062.76301068144402961),
_mm256_set1_ps(542.38938216885594557),
_mm256_set1_ps(328.24391161267163852),
_mm256_set1_ps(219.84388820771997075),
_mm256_set1_ps(157.49790211713116150),
_mm256_set1_ps(118.38163379128440056),
_mm256_set1_ps(92.23973597198599350),
_mm256_set1_ps(73.90969946516314337),
_mm256_set1_ps(60.56287935539454281),
_mm256_set1_ps(50.54400438065793111),
_mm256_set1_ps(42.83209252697486136),
_mm256_set1_ps(36.76965720093080847),
_mm256_set1_ps(31.91778586717441257),
_mm256_set1_ps(27.97440695989148551),
_mm256_set1_ps(24.72611246534243534),
_mm256_set1_ps(22.01866773467314076),
_mm256_set1_ps(19.73836125017757581),
_mm256_set1_ps(17.79986866518916599),
_mm256_set1_ps(16.13815476591664932),
_mm256_set1_ps(14.70294783067476274),
_mm256_set1_ps(13.45489318988675542),
_mm256_set1_ps(12.36282713537543643),
_mm256_set1_ps(11.40181318331400817),
_mm256_set1_ps(10.55170643021994792),
_mm256_set1_ps(9.79608972017488178),
_mm256_set1_ps(9.12147551075540264),
_mm256_set1_ps(8.51670021771579755),
_mm256_set1_ps(7.97245976035782089),
_mm256_set1_ps(7.48094990150482975),
_mm256_set1_ps(7.03558520545476096),
_mm256_set1_ps(6.63077757000577961),
_mm256_set1_ps(6.26176032547270278),
_mm256_set1_ps(5.92444749260821357),
_mm256_set1_ps(5.61532039144371531),
_mm256_set1_ps(5.33133569098266236),
_mm256_set1_ps(5.06985038851206493),
_mm256_set1_ps(4.82856024770480818),
_mm256_set1_ps(4.60544900516174227),
_mm256_set1_ps(4.39874624527744107),
_mm256_set1_ps(4.20689229309963242),
_mm256_set1_ps(4.02850882009587341),
_mm256_set1_ps(3.86237412456220319),
_mm256_set1_ps(3.70740225596646766),
_mm256_set1_ps(3.56262531497762280),
_mm256_set1_ps(3.42717838884132009),
_mm256_set1_ps(3.30028668303544404),
_mm256_set1_ps(3.18125449075313638),
_mm256_set1_ps(3.06945570625878084),
_mm256_set1_ps(2.96432564001850762),
_mm256_set1_ps(2.86535393539436001),
_mm256_set1_ps(2.77207842067833932),
_mm256_set1_ps(2.68407975793799203),
_mm256_set1_ps(2.60097677280667572),
_mm256_set1_ps(2.52242236796785813),
_mm256_set1_ps(2.44809993843393503),
_mm256_set1_ps(2.37772021942528200),
_mm256_set1_ps(2.31101850820766996),
_mm256_set1_ps(2.24775221004035153),
_mm256_set1_ps(2.18769866574001037),
_mm256_set1_ps(2.13065322453304118),
_mm256_set1_ps(2.07642753105719757),
_mm256_set1_ps(2.02484799975171637),
_mm256_set1_ps(1.97575445357971491),
_mm256_set1_ps(1.92899890717015055),
_mm256_set1_ps(1.88444447714107866),
_mm256_set1_ps(1.84196440464714883),
_mm256_set1_ps(1.80144117714498853),
_mm256_set1_ps(1.76276573804219461),
_mm256_set1_ps(1.72583677433227445),
_mm256_set1_ps(1.69056007355504945),
_mm256_set1_ps(1.65684794248971023),
_mm256_set1_ps(1.62461868091121531),
_mm256_set1_ps(1.59379610454104137),
_mm256_set1_ps(1.56430911201842116),
_mm256_set1_ps(1.53609129132299893),
_mm256_set1_ps(1.50908056160717741),
_mm256_set1_ps(1.48321884685699712),
_mm256_set1_ps(1.45845177820343341),
_mm256_set1_ps(1.43472842205929152),
_mm256_set1_ps(1.41200103156710211),
_mm256_set1_ps(1.39022481911630291),
_mm256_set1_ps(1.36935774792838338),
_mm256_set1_ps(1.34936034092082524),
_mm256_set1_ps(1.33019550524813823),
_mm256_set1_ps(1.31182837108428330),
_mm256_set1_ps(1.29422614335786368),
_mm256_set1_ps(1.27735796528209100),
_mm256_set1_ps(1.26119479263764811),
_mm256_set1_ps(1.24570927786994146),
_mm256_set1_ps(1.23087566315441310),
_mm256_set1_ps(1.21666968166584666),
_mm256_set1_ps(1.20306846636115816),
_mm256_set1_ps(1.19005046565099137),
_mm256_set1_ps(1.17759536539444132),
_mm256_set1_ps(1.16568401670416888),
_mm256_set1_ps(1.15429836909675965),
_mm256_set1_ps(1.14342140856596108),
_mm256_set1_ps(1.13303710019499548),
_mm256_set1_ps(1.12313033495892256),
_mm256_set1_ps(1.11368688039941000),
_mm256_set1_ps(1.10469333488267329),
_mm256_set1_ps(1.09613708517704822),
_mm256_set1_ps(1.08800626710994641),
_mm256_set1_ps(1.08028972908511101),
_mm256_set1_ps(1.07297699826029014),
_mm256_set1_ps(1.06605824920296688),
_mm256_set1_ps(1.05952427485770095),
_mm256_set1_ps(1.05336645967322906),
_mm256_set1_ps(1.04757675475074774),
_mm256_set1_ps(1.04214765488702521),
_mm256_set1_ps(1.03707217739714941),
_mm256_set1_ps(1.03234384261201928),
_mm256_set1_ps(1.02795665595514785),
_mm256_set1_ps(1.02390509151210085),
_mm256_set1_ps(1.02018407701397673),
_mm256_set1_ps(1.01678898016385832),
_mm256_set1_ps(1.01371559624214091),
_mm256_set1_ps(1.01096013693317843),
_mm256_set1_ps(1.00851922032179320),
_mm256_set1_ps(1.00638986201394287),
_mm256_set1_ps(1.00456946734126218),
_mm256_set1_ps(1.00305582461434328),
_mm256_set1_ps(1.00184709939451322),
_mm256_set1_ps(1.00094182975856194),
_mm256_set1_ps(1.00033892253539203),
_mm256_set1_ps(1.00003765049793447)
};
const __m256 inv2n = _mm256_rcp_ps(_mm256_set1_ps(static_cast<float>(n+n)));
const __m256 xx    = _mm256_mul_ps(x,x);
const __m256 C1    = _mm256_set1_ps(-0.5f);
__m256 sum = _mm256_setzero_ps();
for(std::int32_t j = 0; j<n; ++j)  
{
    __m256 xmm0    = theta_lut[j];
    __m256 exp_arg = _mm256_mul_ps(C1,_mm256_mul_ps(xx,xmm0));
    __m256 exp_val = _mm256_exp_ps(exp_arg);
    sum            = _mm256_add_ps(sum,exp_val); 
}
return (_mm256_mul_ps(inv2n,sum));
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
gaussian_Q_approx_sadhwani_summed_4xf64(const __m256d x,const std::int32_t n) 
{
if (__builtin_expect(n>128,0)) {return (_mm256_set1_pd(-1.0));}
__ATTR_ALIGN__(32) 
const __m256d theta_lut[128] = 
{
_mm256_set1_pd(26561.07370058031301596),
_mm256_set1_pd(2951.52672978394593883),
_mm256_set1_pd(1062.76301068144402961),
_mm256_set1_pd(542.38938216885594557),
_mm256_set1_pd(328.24391161267163852),
_mm256_set1_pd(219.84388820771997075),
_mm256_set1_pd(157.49790211713116150),
_mm256_set1_pd(118.38163379128440056),
_mm256_set1_pd(92.23973597198599350),
_mm256_set1_pd(73.90969946516314337),
_mm256_set1_pd(60.56287935539454281),
_mm256_set1_pd(50.54400438065793111),
_mm256_set1_pd(42.83209252697486136),
_mm256_set1_pd(36.76965720093080847),
_mm256_set1_pd(31.91778586717441257),
_mm256_set1_pd(27.97440695989148551),
_mm256_set1_pd(24.72611246534243534),
_mm256_set1_pd(22.01866773467314076),
_mm256_set1_pd(19.73836125017757581),
_mm256_set1_pd(17.79986866518916599),
_mm256_set1_pd(16.13815476591664932),
_mm256_set1_pd(14.70294783067476274),
_mm256_set1_pd(13.45489318988675542),
_mm256_set1_pd(12.36282713537543643),
_mm256_set1_pd(11.40181318331400817),
_mm256_set1_pd(10.55170643021994792),
_mm256_set1_pd(9.79608972017488178),
_mm256_set1_pd(9.12147551075540264),
_mm256_set1_pd(8.51670021771579755),
_mm256_set1_pd(7.97245976035782089),
_mm256_set1_pd(7.48094990150482975),
_mm256_set1_pd(7.03558520545476096),
_mm256_set1_pd(6.63077757000577961),
_mm256_set1_pd(6.26176032547270278),
_mm256_set1_pd(5.92444749260821357),
_mm256_set1_pd(5.61532039144371531),
_mm256_set1_pd(5.33133569098266236),
_mm256_set1_pd(5.06985038851206493),
_mm256_set1_pd(4.82856024770480818),
_mm256_set1_pd(4.60544900516174227),
_mm256_set1_pd(4.39874624527744107),
_mm256_set1_pd(4.20689229309963242),
_mm256_set1_pd(4.02850882009587341),
_mm256_set1_pd(3.86237412456220319),
_mm256_set1_pd(3.70740225596646766),
_mm256_set1_pd(3.56262531497762280),
_mm256_set1_pd(3.42717838884132009),
_mm256_set1_pd(3.30028668303544404),
_mm256_set1_pd(3.18125449075313638),
_mm256_set1_pd(3.06945570625878084),
_mm256_set1_pd(2.96432564001850762),
_mm256_set1_pd(2.86535393539436001),
_mm256_set1_pd(2.77207842067833932),
_mm256_set1_pd(2.68407975793799203),
_mm256_set1_pd(2.60097677280667572),
_mm256_set1_pd(2.52242236796785813),
_mm256_set1_pd(2.44809993843393503),
_mm256_set1_pd(2.37772021942528200),
_mm256_set1_pd(2.31101850820766996),
_mm256_set1_pd(2.24775221004035153),
_mm256_set1_pd(2.18769866574001037),
_mm256_set1_pd(2.13065322453304118),
_mm256_set1_pd(2.07642753105719757),
_mm256_set1_pd(2.02484799975171637),
_mm256_set1_pd(1.97575445357971491),
_mm256_set1_pd(1.92899890717015055),
_mm256_set1_pd(1.88444447714107866),
_mm256_set1_pd(1.84196440464714883),
_mm256_set1_pd(1.80144117714498853),
_mm256_set1_pd(1.76276573804219461),
_mm256_set1_pd(1.72583677433227445),
_mm256_set1_pd(1.69056007355504945),
_mm256_set1_pd(1.65684794248971023),
_mm256_set1_pd(1.62461868091121531),
_mm256_set1_pd(1.59379610454104137),
_mm256_set1_pd(1.56430911201842116),
_mm256_set1_pd(1.53609129132299893),
_mm256_set1_pd(1.50908056160717741),
_mm256_set1_pd(1.48321884685699712),
_mm256_set1_pd(1.45845177820343341),
_mm256_set1_pd(1.43472842205929152),
_mm256_set1_pd(1.41200103156710211),
_mm256_set1_pd(1.39022481911630291),
_mm256_set1_pd(1.36935774792838338),
_mm256_set1_pd(1.34936034092082524),
_mm256_set1_pd(1.33019550524813823),
_mm256_set1_pd(1.31182837108428330),
_mm256_set1_pd(1.29422614335786368),
_mm256_set1_pd(1.27735796528209100),
_mm256_set1_pd(1.26119479263764811),
_mm256_set1_pd(1.24570927786994146),
_mm256_set1_pd(1.23087566315441310),
_mm256_set1_pd(1.21666968166584666),
_mm256_set1_pd(1.20306846636115816),
_mm256_set1_pd(1.19005046565099137),
_mm256_set1_pd(1.17759536539444132),
_mm256_set1_pd(1.16568401670416888),
_mm256_set1_pd(1.15429836909675965),
_mm256_set1_pd(1.14342140856596108),
_mm256_set1_pd(1.13303710019499548),
_mm256_set1_pd(1.12313033495892256),
_mm256_set1_pd(1.11368688039941000),
_mm256_set1_pd(1.10469333488267329),
_mm256_set1_pd(1.09613708517704822),
_mm256_set1_pd(1.08800626710994641),
_mm256_set1_pd(1.08028972908511101),
_mm256_set1_pd(1.07297699826029014),
_mm256_set1_pd(1.06605824920296688),
_mm256_set1_pd(1.05952427485770095),
_mm256_set1_pd(1.05336645967322906),
_mm256_set1_pd(1.04757675475074774),
_mm256_set1_pd(1.04214765488702521),
_mm256_set1_pd(1.03707217739714941),
_mm256_set1_pd(1.03234384261201928),
_mm256_set1_pd(1.02795665595514785),
_mm256_set1_pd(1.02390509151210085),
_mm256_set1_pd(1.02018407701397673),
_mm256_set1_pd(1.01678898016385832),
_mm256_set1_pd(1.01371559624214091),
_mm256_set1_pd(1.01096013693317843),
_mm256_set1_pd(1.00851922032179320),
_mm256_set1_pd(1.00638986201394287),
_mm256_set1_pd(1.00456946734126218),
_mm256_set1_pd(1.00305582461434328),
_mm256_set1_pd(1.00184709939451322),
_mm256_set1_pd(1.00094182975856194),
_mm256_set1_pd(1.00033892253539203),
_mm256_set1_pd(1.00003765049793447)
};
#if (GAUSSIAN_Q_APPROX_SIMD_RECIPROCAL_BY_DIVISION) == 1
const __m256d one   = _mm256_set1_pd(1.0);
const __m256d inv2n = _mm256_div_pd(one,_mm256_set1_pd(static_cast<float>(n+n)));
#else 
const __m256d inv2n = _mm256_rcp14_pd(_mm256_set1_pd(static_cast<float>(n+n)));
#endif 
const __m256d xx    = _mm256_mul_pd(x,x);
const __m256d C1    = _mm256_set1_pd(-0.5f);
__m256d sum = _mm256_setzero_pd();
for(std::int32_t j = 0; j<n; ++j)  
{
    __m256d xmm0    = theta_lut[j];
    __m256d exp_arg = _mm256_mul_pd(C1,_mm256_mul_pd(xx,xmm0));
    __m256d exp_val = _mm256_exp_pd(exp_arg);
    sum             = _mm256_add_pd(sum,exp_val); 
}
return (_mm256_mul_pd(inv2n,sum));
}

#if (GAUSSIAN_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
__m512
gaussian_Q_approx_sadhwani_summed_16xf32(const __m512 x,const std::int32_t n) 
{
if (__builtin_expect(n>128,0)) {return (_mm512_set1_ps(-1.0));}
__ATTR_ALIGN__(64) 
const __m512 theta_lut[128] = 
{
_mm512_set1_ps(26561.07370058031301596),
_mm512_set1_ps(2951.52672978394593883),
_mm512_set1_ps(1062.76301068144402961),
_mm512_set1_ps(542.38938216885594557),
_mm512_set1_ps(328.24391161267163852),
_mm512_set1_ps(219.84388820771997075),
_mm512_set1_ps(157.49790211713116150),
_mm512_set1_ps(118.38163379128440056),
_mm512_set1_ps(92.23973597198599350),
_mm512_set1_ps(73.90969946516314337),
_mm512_set1_ps(60.56287935539454281),
_mm512_set1_ps(50.54400438065793111),
_mm512_set1_ps(42.83209252697486136),
_mm512_set1_ps(36.76965720093080847),
_mm512_set1_ps(31.91778586717441257),
_mm512_set1_ps(27.97440695989148551),
_mm512_set1_ps(24.72611246534243534),
_mm512_set1_ps(22.01866773467314076),
_mm512_set1_ps(19.73836125017757581),
_mm512_set1_ps(17.79986866518916599),
_mm512_set1_ps(16.13815476591664932),
_mm512_set1_ps(14.70294783067476274),
_mm512_set1_ps(13.45489318988675542),
_mm512_set1_ps(12.36282713537543643),
_mm512_set1_ps(11.40181318331400817),
_mm512_set1_ps(10.55170643021994792),
_mm512_set1_ps(9.79608972017488178),
_mm512_set1_ps(9.12147551075540264),
_mm512_set1_ps(8.51670021771579755),
_mm512_set1_ps(7.97245976035782089),
_mm512_set1_ps(7.48094990150482975),
_mm512_set1_ps(7.03558520545476096),
_mm512_set1_ps(6.63077757000577961),
_mm512_set1_ps(6.26176032547270278),
_mm512_set1_ps(5.92444749260821357),
_mm512_set1_ps(5.61532039144371531),
_mm512_set1_ps(5.33133569098266236),
_mm512_set1_ps(5.06985038851206493),
_mm512_set1_ps(4.82856024770480818),
_mm512_set1_ps(4.60544900516174227),
_mm512_set1_ps(4.39874624527744107),
_mm512_set1_ps(4.20689229309963242),
_mm512_set1_ps(4.02850882009587341),
_mm512_set1_ps(3.86237412456220319),
_mm512_set1_ps(3.70740225596646766),
_mm512_set1_ps(3.56262531497762280),
_mm512_set1_ps(3.42717838884132009),
_mm512_set1_ps(3.30028668303544404),
_mm512_set1_ps(3.18125449075313638),
_mm512_set1_ps(3.06945570625878084),
_mm512_set1_ps(2.96432564001850762),
_mm512_set1_ps(2.86535393539436001),
_mm512_set1_ps(2.77207842067833932),
_mm512_set1_ps(2.68407975793799203),
_mm512_set1_ps(2.60097677280667572),
_mm512_set1_ps(2.52242236796785813),
_mm512_set1_ps(2.44809993843393503),
_mm512_set1_ps(2.37772021942528200),
_mm512_set1_ps(2.31101850820766996),
_mm512_set1_ps(2.24775221004035153),
_mm512_set1_ps(2.18769866574001037),
_mm512_set1_ps(2.13065322453304118),
_mm512_set1_ps(2.07642753105719757),
_mm512_set1_ps(2.02484799975171637),
_mm512_set1_ps(1.97575445357971491),
_mm512_set1_ps(1.92899890717015055),
_mm512_set1_ps(1.88444447714107866),
_mm512_set1_ps(1.84196440464714883),
_mm512_set1_ps(1.80144117714498853),
_mm512_set1_ps(1.76276573804219461),
_mm512_set1_ps(1.72583677433227445),
_mm512_set1_ps(1.69056007355504945),
_mm512_set1_ps(1.65684794248971023),
_mm512_set1_ps(1.62461868091121531),
_mm512_set1_ps(1.59379610454104137),
_mm512_set1_ps(1.56430911201842116),
_mm512_set1_ps(1.53609129132299893),
_mm512_set1_ps(1.50908056160717741),
_mm512_set1_ps(1.48321884685699712),
_mm512_set1_ps(1.45845177820343341),
_mm512_set1_ps(1.43472842205929152),
_mm512_set1_ps(1.41200103156710211),
_mm512_set1_ps(1.39022481911630291),
_mm512_set1_ps(1.36935774792838338),
_mm512_set1_ps(1.34936034092082524),
_mm512_set1_ps(1.33019550524813823),
_mm512_set1_ps(1.31182837108428330),
_mm512_set1_ps(1.29422614335786368),
_mm512_set1_ps(1.27735796528209100),
_mm512_set1_ps(1.26119479263764811),
_mm512_set1_ps(1.24570927786994146),
_mm512_set1_ps(1.23087566315441310),
_mm512_set1_ps(1.21666968166584666),
_mm512_set1_ps(1.20306846636115816),
_mm512_set1_ps(1.19005046565099137),
_mm512_set1_ps(1.17759536539444132),
_mm512_set1_ps(1.16568401670416888),
_mm512_set1_ps(1.15429836909675965),
_mm512_set1_ps(1.14342140856596108),
_mm512_set1_ps(1.13303710019499548),
_mm512_set1_ps(1.12313033495892256),
_mm512_set1_ps(1.11368688039941000),
_mm512_set1_ps(1.10469333488267329),
_mm512_set1_ps(1.09613708517704822),
_mm512_set1_ps(1.08800626710994641),
_mm512_set1_ps(1.08028972908511101),
_mm512_set1_ps(1.07297699826029014),
_mm512_set1_ps(1.06605824920296688),
_mm512_set1_ps(1.05952427485770095),
_mm512_set1_ps(1.05336645967322906),
_mm512_set1_ps(1.04757675475074774),
_mm512_set1_ps(1.04214765488702521),
_mm512_set1_ps(1.03707217739714941),
_mm512_set1_ps(1.03234384261201928),
_mm512_set1_ps(1.02795665595514785),
_mm512_set1_ps(1.02390509151210085),
_mm512_set1_ps(1.02018407701397673),
_mm512_set1_ps(1.01678898016385832),
_mm512_set1_ps(1.01371559624214091),
_mm512_set1_ps(1.01096013693317843),
_mm512_set1_ps(1.00851922032179320),
_mm512_set1_ps(1.00638986201394287),
_mm512_set1_ps(1.00456946734126218),
_mm512_set1_ps(1.00305582461434328),
_mm512_set1_ps(1.00184709939451322),
_mm512_set1_ps(1.00094182975856194),
_mm512_set1_ps(1.00033892253539203),
_mm512_set1_ps(1.00003765049793447)
};
#if (GAUSSIAN_Q_APPROX_SIMD_RECIPROCAL_BY_DIVISION) == 1
const __m512 one   = _mm512_set1_ps(1.0);
const __m512 inv2n = _mm512_div_ps(one,_mm512_set1_ps(static_cast<float>(n+n)));
#else
const __m512 inv2n = _mm512_rcp14_ps(_mm512_set1_ps(static_cast<float>(n+n)));
#endif 
const __m512 xx    = _mm512_mul_ps(x,x);
const __m512 C1    = _mm512_set1_ps(-0.5f);
__m512 sum = _mm512_setzero_ps();
for(std::int32_t j = 0; j<n; ++j)  
{
    __m512 xmm0    = theta_lut[j];
    __m512 exp_arg = _mm512_mul_ps(C1,_mm512_mul_ps(xx,xmm0));
    __m512 exp_val = _mm512_exp_ps(exp_arg);
    sum            = _mm512_add_ps(sum,exp_val); 
}
return (_mm512_mul_ps(inv2n,sum));
}

#if (GAUSSIAN_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
__m512d
gaussian_Q_approx_sadhwani_summed_8xf64(const __m512d x,const std::int32_t n) 
{
if (__builtin_expect(n>128,0)) {return (_mm512_set1_pd(-1.0));}
__ATTR_ALIGN__(64) 
const __m512d theta_lut[128] = 
{
_mm512_set1_pd(26561.07370058031301596),
_mm512_set1_pd(2951.52672978394593883),
_mm512_set1_pd(1062.76301068144402961),
_mm512_set1_pd(542.38938216885594557),
_mm512_set1_pd(328.24391161267163852),
_mm512_set1_pd(219.84388820771997075),
_mm512_set1_pd(157.49790211713116150),
_mm512_set1_pd(118.38163379128440056),
_mm512_set1_pd(92.23973597198599350),
_mm512_set1_pd(73.90969946516314337),
_mm512_set1_pd(60.56287935539454281),
_mm512_set1_pd(50.54400438065793111),
_mm512_set1_pd(42.83209252697486136),
_mm512_set1_pd(36.76965720093080847),
_mm512_set1_pd(31.91778586717441257),
_mm512_set1_pd(27.97440695989148551),
_mm512_set1_pd(24.72611246534243534),
_mm512_set1_pd(22.01866773467314076),
_mm512_set1_pd(19.73836125017757581),
_mm512_set1_pd(17.79986866518916599),
_mm512_set1_pd(16.13815476591664932),
_mm512_set1_pd(14.70294783067476274),
_mm512_set1_pd(13.45489318988675542),
_mm512_set1_pd(12.36282713537543643),
_mm512_set1_pd(11.40181318331400817),
_mm512_set1_pd(10.55170643021994792),
_mm512_set1_pd(9.79608972017488178),
_mm512_set1_pd(9.12147551075540264),
_mm512_set1_pd(8.51670021771579755),
_mm512_set1_pd(7.97245976035782089),
_mm512_set1_pd(7.48094990150482975),
_mm512_set1_pd(7.03558520545476096),
_mm512_set1_pd(6.63077757000577961),
_mm512_set1_pd(6.26176032547270278),
_mm512_set1_pd(5.92444749260821357),
_mm512_set1_pd(5.61532039144371531),
_mm512_set1_pd(5.33133569098266236),
_mm512_set1_pd(5.06985038851206493),
_mm512_set1_pd(4.82856024770480818),
_mm512_set1_pd(4.60544900516174227),
_mm512_set1_pd(4.39874624527744107),
_mm512_set1_pd(4.20689229309963242),
_mm512_set1_pd(4.02850882009587341),
_mm512_set1_pd(3.86237412456220319),
_mm512_set1_pd(3.70740225596646766),
_mm512_set1_pd(3.56262531497762280),
_mm512_set1_pd(3.42717838884132009),
_mm512_set1_pd(3.30028668303544404),
_mm512_set1_pd(3.18125449075313638),
_mm512_set1_pd(3.06945570625878084),
_mm512_set1_pd(2.96432564001850762),
_mm512_set1_pd(2.86535393539436001),
_mm512_set1_pd(2.77207842067833932),
_mm512_set1_pd(2.68407975793799203),
_mm512_set1_pd(2.60097677280667572),
_mm512_set1_pd(2.52242236796785813),
_mm512_set1_pd(2.44809993843393503),
_mm512_set1_pd(2.37772021942528200),
_mm512_set1_pd(2.31101850820766996),
_mm512_set1_pd(2.24775221004035153),
_mm512_set1_pd(2.18769866574001037),
_mm512_set1_pd(2.13065322453304118),
_mm512_set1_pd(2.07642753105719757),
_mm512_set1_pd(2.02484799975171637),
_mm512_set1_pd(1.97575445357971491),
_mm512_set1_pd(1.92899890717015055),
_mm512_set1_pd(1.88444447714107866),
_mm512_set1_pd(1.84196440464714883),
_mm512_set1_pd(1.80144117714498853),
_mm512_set1_pd(1.76276573804219461),
_mm512_set1_pd(1.72583677433227445),
_mm512_set1_pd(1.69056007355504945),
_mm512_set1_pd(1.65684794248971023),
_mm512_set1_pd(1.62461868091121531),
_mm512_set1_pd(1.59379610454104137),
_mm512_set1_pd(1.56430911201842116),
_mm512_set1_pd(1.53609129132299893),
_mm512_set1_pd(1.50908056160717741),
_mm512_set1_pd(1.48321884685699712),
_mm512_set1_pd(1.45845177820343341),
_mm512_set1_pd(1.43472842205929152),
_mm512_set1_pd(1.41200103156710211),
_mm512_set1_pd(1.39022481911630291),
_mm512_set1_pd(1.36935774792838338),
_mm512_set1_pd(1.34936034092082524),
_mm512_set1_pd(1.33019550524813823),
_mm512_set1_pd(1.31182837108428330),
_mm512_set1_pd(1.29422614335786368),
_mm512_set1_pd(1.27735796528209100),
_mm512_set1_pd(1.26119479263764811),
_mm512_set1_pd(1.24570927786994146),
_mm512_set1_pd(1.23087566315441310),
_mm512_set1_pd(1.21666968166584666),
_mm512_set1_pd(1.20306846636115816),
_mm512_set1_pd(1.19005046565099137),
_mm512_set1_pd(1.17759536539444132),
_mm512_set1_pd(1.16568401670416888),
_mm512_set1_pd(1.15429836909675965),
_mm512_set1_pd(1.14342140856596108),
_mm512_set1_pd(1.13303710019499548),
_mm512_set1_pd(1.12313033495892256),
_mm512_set1_pd(1.11368688039941000),
_mm512_set1_pd(1.10469333488267329),
_mm512_set1_pd(1.09613708517704822),
_mm512_set1_pd(1.08800626710994641),
_mm512_set1_pd(1.08028972908511101),
_mm512_set1_pd(1.07297699826029014),
_mm512_set1_pd(1.06605824920296688),
_mm512_set1_pd(1.05952427485770095),
_mm512_set1_pd(1.05336645967322906),
_mm512_set1_pd(1.04757675475074774),
_mm512_set1_pd(1.04214765488702521),
_mm512_set1_pd(1.03707217739714941),
_mm512_set1_pd(1.03234384261201928),
_mm512_set1_pd(1.02795665595514785),
_mm512_set1_pd(1.02390509151210085),
_mm512_set1_pd(1.02018407701397673),
_mm512_set1_pd(1.01678898016385832),
_mm512_set1_pd(1.01371559624214091),
_mm512_set1_pd(1.01096013693317843),
_mm512_set1_pd(1.00851922032179320),
_mm512_set1_pd(1.00638986201394287),
_mm512_set1_pd(1.00456946734126218),
_mm512_set1_pd(1.00305582461434328),
_mm512_set1_pd(1.00184709939451322),
_mm512_set1_pd(1.00094182975856194),
_mm512_set1_pd(1.00033892253539203),
_mm512_set1_pd(1.00003765049793447)
};
#if (GAUSSIAN_Q_APPROX_SIMD_RECIPROCAL_BY_DIVISION) == 1
const __m512d one   = _mm512_set1_pd(1.0);
const __m512d inv2n = _mm512_div_pd(one,_mm512_set1_pd(static_cast<double>(n+n)));
#else 
const __m512d inv2n = _mm512_rcp14_pd(_mm512_set1_pd(static_cast<double>(n+n)));
#endif 
const __m512d xx    = _mm512_mul_pd(x,x);
const __m512d C1    = _mm512_set1_pd(-0.5f);
__m512d sum = _mm512_setzero_pd();
for(std::int32_t j = 0; j<n; ++j)  
{
    __m512d xmm0    = theta_lut[j];
    __m512d exp_arg = _mm512_mul_pd(C1,_mm512_mul_pd(xx,xmm0));
    __m512d exp_val = _mm512_exp_pd(exp_arg);
    sum             = _mm512_add_pd(sum,exp_val); 
}
return (_mm512_mul_pd(inv2n,sum));
}

#if (GAUSSIAN_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=sse
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
__m128 
gaussian_Q_approx_karagiannidis_lioumpas_4xf32(const __m128 x)
{
if(__builtin_expect(_mm_cmp_ps_mask(x,_mm_setzero_ps(),_CMP_EQ_OQ),0))
{
   return (_mm_set1_ps(std::numeric_limits<float>::quiet_NaN());)
}
const __m128 a       = _mm_mul_ps(x,_mm_set1_ps(-1.98f));
const __m128 xx      = _mm_mul_ps(x,x);
const __m128 sqrt2PI = _mm_set1_ps(2.506628274631000502415765284811f);
const __m128 halfxx  = _mm_mul_ps(_mm_set1_ps(0.5f),xx);
const __m128 den     = _mm_mul_ps(_mm_set1_ps(1.135f,_mm_mul_ps(x,sqrt2PI)));
const __m128 exp_val1= _mm_mul_ps(_mm_set1_ps(1.0),simd_fast_exp_approx_4xf32(_mm_mul_ps(a,x)));
const __m128 exp_val2= simd_fast_exp_approx_4xf32(_mm_sub_ps(_mm_setzero_ps(),xx));
return (_mm_div_ps(_mm_mul_ps(exp_val1,exp_val2),den));
}

#if (GAUSSIAN_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=sse
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
__m128d 
gaussian_Q_approx_karagiannidis_lioumpas_2xf64(const __m128d x)
{
if(__builtin_expect(_mm_cmp_pd_mask(x,_mm_setzero_pd(),_CMP_EQ_OQ),0))
{
   return (_mm_set1_pd(std::numeric_limits<double>::quiet_NaN());)
}
const __m128d a       = _mm_mul_pd(x,_mm_set1_pd(-1.98));
const __m128d xx      = _mm_mul_pd(x,x);
const __m128d sqrt2PI = _mm_set1_pd(2.506628274631000502415765284811);
const __m128d halfxx  = _mm_mul_pd(_mm_set1_pd(0.5),xx);
const __m128d den     = _mm_mul_pd(_mm_set1_pd(1.135,_mm_mul_pd(x,sqrt2PI)));
const __m128d exp_val1= _mm_mul_pd(_mm_set1_pd(1.0),simd_fast_exp_approx_2xf64(_mm_mul_pd(a,x)));
const __m128d exp_val2= simd_fast_exp_approx_2xf64(_mm_sub_pd(_mm_setzero_pd(),xx));
return (_mm_div_pd(_mm_mul_pd(exp_val1,exp_val2),den));
}

}
}
#endif /*__GMS_GAUSSIAN_Q_APPROX_SIMD_H__*/