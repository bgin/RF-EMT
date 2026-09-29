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

}

}
#endif /*__GMS_GAUSSIAN_Q_APPROX_SIMD_H__*/