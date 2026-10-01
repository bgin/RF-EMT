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

#ifndef __GMS_MARCUM_Q_APPROX_SIMD_H__
#define __GMS_MARCUM_Q_APPROX_SIMD_H__ 3020260826 


#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_gaussian_Q_approx_simd.h"

namespace file_info 
{

     static const unsigned int GMS_MARCUM_Q_APPROX_SIMD_MAJOR = 1;
     static const unsigned int GMS_MARCUM_Q_APPROX_SIMD_MINOR = 1;
     static const unsigned int GMS_MARCUM_Q_APPROX_SIMD_MICRO = 0;
     static const unsigned int GMS_MARCUM_Q_APPROX_SIMD_FULLVER =
       1000U*GMS_MARCUM_Q_APPROX_SIMD_MAJOR+100U*GMS_MARCUM_Q_APPROX_SIMD_MINOR+
       10U*GMS_MARCUM_Q_APPROX_SIMD_MICRO;
     static const char GMS_MARCUM_Q_APPROX_SIMD_CREATION_DATE[] = "30-09-2026 08:27AM +00200 (WED 30 AUG 2026 GMT+2)";
     static const char GMS_MARCUM_Q_APPROX_SIMD_BUILD_DATE[]    = __DATE__; 
     static const char GMS_MARCUM_Q_APPROX_SIMD_BUILD_TIME[]    = __TIME__;
     static const char GMS_MARCUM_Q_APPROX_SIMD_AUTHOR[]        = "Bernard Gingold, beniekg@gmail.com";
     static const char GMS_MARCUM_Q_APPROX_SIMD_SYNOPSIS[]      = "Marcum-Q function approximations SIMD vectorized, based on the Wikipedia article: https://en.wikipedia.org/wiki/Marcum_Q-function#Asymptotic_forms ";

}

#if !defined(MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE)
#define MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE 0
#endif

#if !defined(MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY)
#define MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY 0
#endif 

namespace gms 
{

namespace math
{


//__m128 
//marcum_Q_approx_chiani_4xf32(const __m128 mu,const __m128 a,
//                             const __m128 b)
//{
//#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
//const __mmask8 is_a_eq_b = _mm_cmp_ps_mask(a,b,_CMP_EQ_OQ);
//if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm_set1_ps(-1.0));}
//#endif
/*
const __m128   one      = _mm_set1_ps(1.0f);
const __m128   mu_half  = _mm_sub_ps(mu,_mm_set1_ps(0.5f));
const __mmask8 is_b_gt_a= _mm_cmp_ps_mask(b,a,_CMP_EQ_OQ);
const __m128   pow_term = _mm_pow_ps(_mm_div_ps(b,a),mu_half);
const __m128   b_sub_a  = _mm_sub_ps(b,a);
const __m128   br_b_gt_a= _mm_mul_ps(pow_term,gaussian_Q_approx_chiani_4xf32(b_sub_a));
const __m128   a_sub_b  = _mm_sub_ps(a,b);
const __m128   br_a_gt_b= _mm_sub_ps(one,_mm_mul_ps(pow_term,gaussian_Q_approx_chiani_4xf32(a_sub_b)));
*/

//COMMON_BODY_FUNC_BLOCK_4XF32(gaussian_Q_approx_chiani_4xf32(b_sub_a),gaussian_Q_approx_chiani_4xf32(a_sub_b))
//return (_mm_mask_blend_ps(is_b_gt_a,br_a_gt_b,br_b_gt_a));
//}

#define COMMON_BODY_FUNC_BLOCK_4XF32(name1,name2)\
    const __m128   one      = _mm_set1_ps(1.0f);\
    const __m128   mu_half  = _mm_sub_ps(mu,_mm_set1_ps(0.5f));\
    const __mmask8 is_b_gt_a= _mm_cmp_ps_mask(b,a,_CMP_EQ_OQ);\
    const __m128   pow_term = _mm_pow_ps(_mm_div_ps(b,a),mu_half);\
    const __m128   b_sub_a  = _mm_sub_ps(b,a);\
    const __m128   br_b_gt_a= _mm_mul_ps(pow_term,name1);\
    const __m128   a_sub_b  = _mm_sub_ps(a,b);\
    const __m128   br_a_gt_b= _mm_sub_ps(one,_mm_mul_ps(pow_term,name2));

#define COMMON_BODY_FUNC_BLOCK_2XF64(name1,name2)\
    const __m128d   one      = _mm_set1_pd(1.0);\
    const __m128d   mu_half  = _mm_sub_pd(mu,_mm_set1_pd(0.5));\
    const __mmask8  is_b_gt_a= _mm_cmp_pd_mask(b,a,_CMP_EQ_OQ);
    const __m128d   pow_term = _mm_pow_pd(_mm_div_pd(b,a),mu_half);\
    const __m128d   b_sub_a  = _mm_sub_pd(b,a);\
    const __m128d   br_b_gt_a= _mm_mul_pd(pow_term,name1);\
    const __m128d   a_sub_b  = _mm_sub_pd(a,b);\
    const __m128d   br_a_gt_b= _mm_sub_pd(one,_mm_mul_ps(pow_term,name2));

#define COMMON_BODY_FUNC_BLOCK_8XF32(name1,name2)
    const __m256   one      = _mm256_set1_ps(1.0f);\
    const __m256   mu_half  = _mm256_sub_ps(mu,_mm256_set1_ps(0.5f));\
    const __mmask8 is_b_gt_a= _mm256_cmp_ps_mask(b,a,_CMP_EQ_OQ);\
    const __m256   pow_term = _mm256_pow_ps(_mm256_div_ps(b,a),mu_half);\
    const __m256   b_sub_a  = _mm256_sub_ps(b,a);\
    const __m256   br_b_gt_a= _mm256_mul_ps(pow_term,name1);\
    const __m256   a_sub_b  = _mm256_sub_ps(a,b);\
    const __m256   br_a_gt_b= _mm256_sub_ps(one,_mm256_mul_ps(pow_term,name2));

#define COMMON_BODY_FUNC_BLOCK_4XF64(name1,name2)\
    const __m256d   one      = _mm256_set1_pd(1.0);\
    const __m256d   mu_half  = _mm256_sub_pd(mu,_mm256_set1_pd(0.5));\
    const __mmask8  is_b_gt_a= _mm256_cmp_pd_mask(b,a,_CMP_EQ_OQ);\
    const __m256d   pow_term = _mm256_pow_pd(_mm256_div_pd(b,a),mu_half);\
    const __m256d   b_sub_a  = _mm256_sub_pd(b,a);\
    const __m256d   br_b_gt_a= _mm256_mul_pd(pow_term,name1);\
    const __m256d   a_sub_b  = _mm256_sub_pd(a,b);\
    const __m256d   br_a_gt_b= _mm256_sub_pd(one,_mm256_mul_pd(pow_term,name2));

#define COMMON_BODY_FUNC_BLOCK_16XF32(name1,name2)\
    const __m512   one      = _mm512_set1_ps(1.0f);\
    const __m512   mu_half  = _mm512_sub_ps(mu,_mm512_set1_ps(0.5f));\
    const __mmask16 is_b_gt_a= _mm512_cmp_ps_mask(b,a,_CMP_EQ_OQ);\
    const __m512   pow_term = _mm512_pow_ps(_mm512_div_ps(b,a),mu_half);\
    const __m512   b_sub_a  = _mm512_sub_ps(b,a);\
    const __m512   br_b_gt_a= _mm512_mul_ps(pow_term,name1);\
    const __m512   a_sub_b  = _mm512_sub_ps(a,b);\
    const __m512   br_a_gt_b= _mm512_sub_ps(one,_mm512_mul_ps(pow_term,name2)); 

#define COMMON_BODY_FUNC_BLOCK_8XF64(name1,name2)\
    const __m512d   one      = _mm512_set1_pd(1.0);\
    const __m512d   mu_half  = _mm512_sub_pd(mu,_mm512_set1_pd(0.5));\
    const __mmask8 is_b_gt_a= _mm512_cmp_pd_mask(b,a,_CMP_EQ_OQ);\
    const __m512d   pow_term = _mm512_pow_pd(_mm512_div_pd(b,a),mu_half);\
    const __m512d   b_sub_a  = _mm512_sub_pd(b,a);\
    const __m512d   br_b_gt_a= _mm512_mul_pd(pow_term,name1);\
    const __m512d   a_sub_b  = _mm512_sub_pd(a,b);\
    const __m512d   br_a_gt_b= _mm512_sub_pd(one,_mm512_mul_pd(pow_term,name2));

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_chiani_4xf32(const __m128 mu,const __m128 a,
                             const __m128 b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm_cmp_ps_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm_set1_ps(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_4XF32(gaussian_Q_approx_chiani_4xf32(b_sub_a),gaussian_Q_approx_chiani_4xf32(a_sub_b))
return (_mm_mask_blend_ps(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_chiani_2xf64(const __m128d mu,const __m128d a,
                             const __m128d b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm_cmp_pd_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm_set1_pd(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_2XF64(gaussian_Q_approx_chiani_2xf64(b_sub_a),gaussian_Q_approx_chiani_2xf64(a_sub_b))
return (_mm_mask_blend_pd(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
__m256 marcum_Q_approx_chiani_8xf32(const __m256 mu,const __m256 a,
                                    const __m256 b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm256_cmp_ps_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm256_set1_ps(-1.0));}
#endif
COMMON_BODY_FUNC_BLOCK_8XF32(gaussian_Q_approx_chiani_8xf32(b_sub_a),gaussian_Q_approx_chiani_8xf32(a_sub_b))
return (_mm256_mask_blend_ps(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
__m256d marcum_Q_approx_chiani_4xf64(const __m256d mu,const __m256d a,
                                    const __m256d b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm256_cmp_pd_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm256_set1_pd(-1.0));}
#endif
COMMON_BODY_FUNC_BLOCK_4XF64(gaussian_Q_approx_chiani_4xf64(b_sub_a),gaussian_Q_approx_chiani_4xf64(a_sub_b))
return (_mm256_mask_blend_pd(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
__m512 marcum_Q_approx_chiani_16xf32(const __m512 mu,const __m512 a,
                                    const __m512 b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask16 is_a_eq_b = _mm512_cmp_ps_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFFFF==is_a_eq_b,0)) { return (_mm512_set1_ps(-1.0));}
#endif
COMMON_BODY_FUNC_BLOCK_16XF32(gaussian_Q_approx_chiani_16xf32(b_sub_a),gaussian_Q_approx_chiani_16xf32(a_sub_b))
return (_mm512_mask_blend_ps(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
__m512d marcum_Q_approx_chiani_8xf64(const __m512d mu,const __m512d a,
                                    const __m512d b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm512_cmp_pd_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm512_set1_pd(-1.0));}
#endif
COMMON_BODY_FUNC_BLOCK_8XF64(gaussian_Q_approx_chiani_8xf64(b_sub_a),gaussian_Q_approx_chiani_8xf64(a_sub_b))
return (_mm512_mask_blend_pd(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_loskot_2T_4xf32(const __m128 mu,const __m128 a,
                                const __m128 b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm_cmp_ps_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm_set1_ps(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_4XF32(gaussian_Q_approx_loskot_2T_4xf32(b_sub_a),gaussian_Q_approx_Loskot_2T_4xf32(a_sub_b))
return (_mm_mask_blend_ps(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_loskot_2T_2xf64(const __m128d mu,const __m128d a,
                                const __m128d b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm_cmp_pd_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm_set1_pd(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_2XF64(gaussian_Q_approx_loskot_2T_2xf64(b_sub_a),gaussian_Q_approx_Loskot_2T_2xf64(a_sub_b))
return (_mm_mask_blend_pd(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_loskot_2T_8xf32(const __m256 mu,const __m256 a,
                                const __m256 b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm256_cmp_ps_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm256_set1_ps(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_8XF32(gaussian_Q_approx_loskot_2T_8xf32(b_sub_a),gaussian_Q_approx_loskot_2T_8xf32(a_sub_b))
return (_mm256_mask_blend_ps(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_loskot_2T_4xf64(const __m256d mu,const __m256d a,
                                const __m256d b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm256_cmp_pd_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm256_set1_pd(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_4XF64(gaussian_Q_approx_loskot_2T_4xf64(b_sub_a),gaussian_Q_approx_loskot_2T_4xf64(a_sub_b))
return (_mm256_mask_blend_pd(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_loskot_2T_16xf32(const __m512 mu,const __m512 a,
                                const __m512 b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask16 is_a_eq_b = _mm512_cmp_ps_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFFFF==is_a_eq_b,0)) { return (_mm512_set1_ps(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_16XF32(gaussian_Q_approx_loskot_2T_16xf32(b_sub_a),gaussian_Q_approx_loskot_2T_16xf32(a_sub_b))
return (_mm256_mask_blend_ps(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_loskot_2T_8xf64(const __m512d mu,const __m512d a,
                                const __m512d b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm512_cmp_pd_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm512_set1_pd(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_8XF64(gaussian_Q_approx_loskot_2T_8xf64(b_sub_a),gaussian_Q_approx_loskot_2T_8xf64(a_sub_b))
return (_mm512_mask_blend_pd(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_loskot_3T_4xf32(const __m128 mu,const __m128 a,
                                const __m128 b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm_cmp_ps_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm_set1_ps(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_4XF32(gaussian_Q_approx_loskot_3T_4xf32(b_sub_a),gaussian_Q_approx_loskot_3T_4xf32(a_sub_b))
return (_mm_mask_blend_ps(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_loskot_3T_2xf64(const __m128d mu,const __m128d a,
                                const __m128d b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm_cmp_pd_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm_set1_pd(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_2XF64(gaussian_Q_approx_loskot_3T_2xf64(b_sub_a),gaussian_Q_approx_loskot_3T_2xf64(a_sub_b))
return (_mm_mask_blend_pd(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_loskot_3T_8xf32(const __m256 mu,const __m256 a,
                                const __m256 b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm256_cmp_ps_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm256_set1_ps(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_8XF32(gaussian_Q_approx_loskot_3T_8xf32(b_sub_a),gaussian_Q_approx_loskot_3T_8xf32(a_sub_b))
return (_mm256_mask_blend_ps(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_loskot_3T_4xf64(const __m256d mu,const __m256d a,
                                const __m256d b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm256_cmp_pd_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm256_set1_pd(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_4XF64(gaussian_Q_approx_loskot_3T_4xf64(b_sub_a),gaussian_Q_approx_loskot_3T_4xf64(a_sub_b))
return (_mm256_mask_blend_pd(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_loskot_3T_16xf32(const __m512 mu,const __m512 a,
                                const __m512 b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask16 is_a_eq_b = _mm512_cmp_ps_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFFFF==is_a_eq_b,0)) { return (_mm512_set1_ps(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_16XF32(gaussian_Q_approx_loskot_3T_16xf32(b_sub_a),gaussian_Q_approx_loskot_3T_16xf32(a_sub_b))
return (_mm512_mask_blend_ps(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_loskot_3T_8xf64(const __m512d mu,const __m512d a,
                                const __m512d b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm512_cmp_pd_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm512_set1_pd(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_8XF64(gaussian_Q_approx_loskot_3T_8xf64(b_sub_a),gaussian_Q_approx_loskot_3T_8xf64(a_sub_b))
return (_mm512_mask_blend_pd(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_sadhwani_1T_4xf32(const __m128 mu,const __m128 a,
                                const __m128 b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm_cmp_ps_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm_set1_ps(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_4XF32(gaussian_Q_approx_sadhwani_1T_4xf32(b_sub_a),gaussian_Q_approx_sadhwani_1T_4xf32(a_sub_b))
return (_mm_mask_blend_ps(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_sadhwani_1T_2xf64(const __m128d mu,const __m128d a,
                                const __m128d b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm_cmp_pd_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm_set1_pd(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_2XF64(gaussian_Q_approx_sadhwani_1T_2xf64(b_sub_a),gaussian_Q_approx_sadhwani_1T_2xf64(a_sub_b))
return (_mm_mask_blend_pd(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_sadhwani_1T_8xf32(const __m256 mu,const __m256 a,
                                const __m256 b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm256_cmp_ps_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm256_set1_ps(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_8XF32(gaussian_Q_approx_sadhwani_1T_8xf32(b_sub_a),gaussian_Q_approx_sadhwani_1T_8xf32(a_sub_b))
return (_mm256_mask_blend_ps(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_sadhwani_1T_4xf64(const __m256d mu,const __m256d a,
                                const __m256d b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm256_cmp_pd_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm256_set1_pd(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_4XF64(gaussian_Q_approx_sadhwani_1T_4xf64(b_sub_a),gaussian_Q_approx_sadhwani_1T_4xf64(a_sub_b))
return (_mm256_mask_blend_pd(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_sadhwani_1T_16xf32(const __m512 mu,const __m512 a,
                                const __m512 b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask16 is_a_eq_b = _mm512_cmp_ps_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFFFF==is_a_eq_b,0)) { return (_mm512_set1_ps(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_16XF32(gaussian_Q_approx_sadhwani_1T_16xf32(b_sub_a),gaussian_Q_approx_sadhwani_1T_16xf32(a_sub_b))
return (_mm512_mask_blend_ps(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_sadhwani_1T_8xf64(const __m512d mu,const __m512d a,
                                const __m512d b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm512_cmp_pd_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm512_set1_pd(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_8XF64(gaussian_Q_approx_sadhwani_1T_8xf64(b_sub_a),gaussian_Q_approx_sadhwani_1T_8xf64(a_sub_b))
return (_mm512_mask_blend_pd(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_sadhwani_2T_4xf32(const __m128 mu,const __m128 a,
                                const __m128 b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm_cmp_ps_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm_set1_ps(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_4XF32(gaussian_Q_approx_sadhwani_2T_4xf32(b_sub_a),gaussian_Q_approx_sadhwani_2T_4xf32(a_sub_b))
return (_mm_mask_blend_ps(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_sadhwani_2T_2xf64(const __m128d mu,const __m128d a,
                                const __m128d b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm_cmp_pd_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm_set1_pd(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_2XF64(gaussian_Q_approx_sadhwani_2T_2xf64(b_sub_a),gaussian_Q_approx_sadhwani_2T_2xf64(a_sub_b))
return (_mm_mask_blend_pd(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_sadhwani_2T_8xf32(const __m256 mu,const __m256 a,
                                const __m256 b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm256_cmp_ps_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm256_set1_ps(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_8XF32(gaussian_Q_approx_sadhwani_2T_8xf32(b_sub_a),gaussian_Q_approx_sadhwani_2T_8xf32(a_sub_b))
return (_mm256_mask_blend_ps(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_sadhwani_2T_4xf64(const __m256d mu,const __m256d a,
                                const __m256d b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm256_cmp_pd_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm256_set1_pd(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_4XF64(gaussian_Q_approx_sadhwani_2T_4xf64(b_sub_a),gaussian_Q_approx_sadhwani_2T_4xf64(a_sub_b))
return (_mm256_mask_blend_pd(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_sadhwani_2T_16xf32(const __m512 mu,const __m512 a,
                                const __m512 b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask16 is_a_eq_b = _mm512_cmp_ps_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFFFF==is_a_eq_b,0)) { return (_mm512_set1_ps(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_16XF32(gaussian_Q_approx_sadhwani_2T_16xf32(b_sub_a),gaussian_Q_approx_sadhwani_2T_16xf32(a_sub_b))
return (_mm512_mask_blend_ps(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_sadhwani_2T_8xf64(const __m512d mu,const __m512d a,
                                const __m512d b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm512_cmp_pd_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm512_set1_pd(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_8XF64(gaussian_Q_approx_sadhwani_2T_8xf64(b_sub_a),gaussian_Q_approx_sadhwani_2T_8xf64(a_sub_b))
return (_mm512_mask_blend_pd(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_sadhwani_4T_4xf32(const __m128 mu,const __m128 a,
                                const __m128 b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm_cmp_ps_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm_set1_ps(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_4XF32(gaussian_Q_approx_sadhwani_4T_4xf32(b_sub_a),gaussian_Q_approx_sadhwani_4T_4xf32(a_sub_b))
return (_mm_mask_blend_ps(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_sadhwani_4T_2xf64(const __m128d mu,const __m128d a,
                                const __m128d b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm_cmp_pd_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm_set1_pd(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_2XF64(gaussian_Q_approx_sadhwani_4T_2xf64(b_sub_a),gaussian_Q_approx_sadhwani_4T_2xf64(a_sub_b))
return (_mm_mask_blend_pd(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_sadhwani_4T_8xf32(const __m256 mu,const __m256 a,
                                const __m256 b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm256_cmp_ps_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm256_set1_ps(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_8XF32(gaussian_Q_approx_sadhwani_4T_8xf32(b_sub_a),gaussian_Q_approx_sadhwani_4T_8xf32(a_sub_b))
return (_mm256_mask_blend_ps(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_sadhwani_4T_4xf64(const __m256d mu,const __m256d a,
                                const __m256d b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm256_cmp_pd_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm256_set1_pd(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_4XF64(gaussian_Q_approx_sadhwani_4T_4xf64(b_sub_a),gaussian_Q_approx_sadhwani_4T_4xf64(a_sub_b))
return (_mm256_mask_blend_pd(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_sadhwani_4T_16xf32(const __m512 mu,const __m512 a,
                                const __m512 b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask16 is_a_eq_b = _mm512_cmp_ps_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFFFF==is_a_eq_b,0)) { return (_mm512_set1_ps(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_16XF32(gaussian_Q_approx_sadhwani_4T_16xf32(b_sub_a),gaussian_Q_approx_sadhwani_4T_16xf32(a_sub_b))
return (_mm512_mask_blend_ps(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_sadhwani_4T_8xf64(const __m512d mu,const __m512d a,
                                const __m512d b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm512_cmp_pd_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm512_set1_pd(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_8XF64(gaussian_Q_approx_sadhwani_4T_8xf64(b_sub_a),gaussian_Q_approx_sadhwani_4T_8xf64(a_sub_b))
return (_mm512_mask_blend_pd(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_borjesson_4xf32(const __m128 mu,const __m128 a,
                                const __m128 b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm_cmp_ps_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm_set1_ps(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_4XF32(gaussian_Q_approx_borjesson_4xf32(b_sub_a),gaussian_Q_approx_borjesson_4xf32(a_sub_b))
return (_mm_mask_blend_ps(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_borjesson_2xf64(const __m128d mu,const __m128d a,
                                const __m128d b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm_cmp_pd_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm_set1_pd(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_2XF64(gaussian_Q_approx_borjesson_2xf64(b_sub_a),gaussian_Q_approx_borjesson_2xf64(a_sub_b))
return (_mm_mask_blend_pd(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_borjesson_8xf32(const __m256 mu,const __m256 a,
                                const __m256 b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm256_cmp_ps_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm256_set1_ps(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_8XF32(gaussian_Q_approx_borjesson_8xf32(b_sub_a),gaussian_Q_approx_borjesson_8xf32(a_sub_b))
return (_mm256_mask_blend_ps(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_borjesson_4xf64(const __m256d mu,const __m256d a,
                                const __m256d b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm256_cmp_pd_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm256_set1_pd(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_4XF64(gaussian_Q_approx_borjesson_4xf64(b_sub_a),gaussian_Q_approx_borjesson_4xf64(a_sub_b))
return (_mm256_mask_blend_pd(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_borjesson_16xf32(const __m512 mu,const __m512 a,
                                const __m512 b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask16 is_a_eq_b = _mm512_cmp_ps_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFFFF==is_a_eq_b,0)) { return (_mm512_set1_ps(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_16XF32(gaussian_Q_approx_borjesson_16xf32(b_sub_a),gaussian_Q_approx_borjesson_16xf32(a_sub_b))
return (_mm512_mask_blend_ps(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

#if (MARCUM_Q_APPROX_SIMD_OVERRIDE_COMPILER_CMD_LINE) == 1
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
marcum_Q_approx_borjesson_8xf64(const __m512d mu,const __m512d a,
                                const __m512d b)
{
#if (MARCUM_Q_APPROX_SIMD_HANDLE_ARGS_A_B_EQUALITY) == 1
const __mmask8 is_a_eq_b = _mm512_cmp_pd_mask(a,b,_CMP_EQ_OQ);
if(__builtin_expect(0xFF==is_a_eq_b,0)) { return (_mm512_set1_pd(-1.0));}
#endif 
COMMON_BODY_FUNC_BLOCK_8XF64(gaussian_Q_approx_borjesson_8xf64(b_sub_a),gaussian_Q_approx_borjesson_8xf64(a_sub_b))
return (_mm512_mask_blend_pd(is_b_gt_a,br_a_gt_b,br_b_gt_a));
}

}

}
#endif /*__GMS_MARCUM_Q_APPROX_SIMD_H__*/