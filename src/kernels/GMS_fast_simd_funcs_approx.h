/*

The MIT License (MIT)

Copyright (c) 2015 Jacques-Henri Jourdan <jourgun@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

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

#ifndef __GMS_FAST_SIMD_FUNCS_APPROX_H__
#define __GMS_FAST_SIMD_FUNCS_APPROX_H__ 1901620260612

namespace file_info 
{

     static const unsigned int GMS_FAST_SIMD_FUNCS_APPROX_MAJOR = 1;
     static const unsigned int GMS_FAST_SIMD_FUNCS_APPROX_MINOR = 1;
     static const unsigned int GMS_FAST_SIMD_FUNCS_APPROX_MICRO = 0;
     static const unsigned int GMS_FAST_SIMD_FUNCS_APPROX_FULLVER =
       1000U*GMS_FAST_SIMD_FUNCS_APPROX_MAJOR+100U*GMS_FAST_SIMD_FUNCS_APPROX_MINOR+
       10U*GMS_FAST_SIMD_FUNCS_APPROX_MICRO;
     static const char GMS_FAST_SIMD_FUNCS_APPROX_CREATION_DATE[] = "19-06-2026 06:14AM +00200 (FRI 19 JUN 2026 GMT+2)";
     static const char GMS_FAST_SIMD_FUNCS_APPROX_BUILD_DATE[]    = __DATE__; 
     static const char GMS_FAST_SIMD_FUNCS_APPROX_BUILD_TIME[]    = __TIME__;
     static const char GMS_FAST_SIMD_FUNCS_APPROX_SYNOPSIS[]      = "Manual vectorization of 'simd_math_prims.h' scalar implementation.";

}

#include <immintrin.h>
#include <cmath>
#include <cstdint>
#include "GMS_config.h"

namespace gms 
{


namespace math 
{
   

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_ALWAYS_INLINE__
static inline 
__m128 
simd_fast_exp_approx_4xf32(const __m128 val)
{
const __m128  exp_cst1{_mm_set1_ps(2139095040.f)};
const __m128  exp_cst2{_mm_setzero_ps()};
const __m128  exp_cst3{_mm_set1_ps(12102203.1615614f)};
const __m128  exp_cst4{_mm_set1_ps(1065353216.f)};
const __m128  C0509871020{_mm_set1_ps(0.509871020f)};
const __m128  C0312146713{_mm_set1_ps(0.312146713f)};
const __m128  C0166617139f{_mm_set1_ps(0.166617139f)};
const __m128  C0166617139{_mm_set1_ps(0.166617139f)};
const __m128  CN2190619930e3{_mm_set1_ps(-2.190619930e-3f)};
const __m128  C13555747234e2{_mm_set1_ps(1.3555747234e-2f)};
const __m128i C0x7F800000{_mm_set1_epi32(0x7F800000)};
const __m128i C0x7FFFFF{_mm_set1_epi32(0x7FFFFF)};
const __m128i C0x3F800000{_mm_set1_epi32(0x3F800000)};
////////////////////////////////////////////////////////////////
union {__m128i v_i; __m128 v_f;} xu,xu2;
__m128  val2,val3,val4,b;
__m128i val4i;
val2    = _mm_fmadd_ps(exp_cst3,val,exp_cst4);
const __mmask8 val2_lt_exp_cst1 = _mm_cmp_ps_mask(val2,exp_cst1,_CMP_LT_OQ);
val3    = _mm_mask_mov_ps(exp_cst1,val2_lt_exp_cst1,val2);
const __mmask8 val3_gt_exp_cst2 = _mm_cmp_ps_mask(val3,exp_cst2,_CMP_GT_OQ);
val4    = _mm_mask_mov_ps(exp_cst2,val3_gt_exp_cst2,val3);
val4i   = _mm_cvtps_epi32(val4);
xu.v_i  = _mm_and_si128(val4i,C0x7F800000);
xu2.v_i = _mm_or_si128(_mm_and_si128(val4i,C0x7FFFFF),C0x3F800000);
b       = xu2.v_f;
__m128  factor = _mm_fmadd_ps(b,C13555747234e2,CN2190619930e3);
__m128  horner_scheme = _mm_fmadd_ps(
                                _mm_fmadd_ps(
                                         _mm_fmadd_ps(factor,b,C0166617139),b,C0312146713),b,C0509871020);
return (_mm_mul_ps(xu.v_f,horner_scheme));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_ALWAYS_INLINE__
static inline 
__m128d
simd_fast_exp_approx_2xf64(const __m128d val)
{
const __m128d exp_cst1{_mm_set1_pd(9218868437227405312.0)};
const __m128d exp_cst2{_mm_setzero_pd()};
const __m128d C6497320848556798092{_mm_set1_pd(6497320848556798.092)};
const __m128d C46071824188000174080{_mm_set1_pd(4607182418800017408.0)};
const __m128d C217150255054231565039e4{_mm_set1_pd(2.17150255054231565039e-4)};
const __m128d CN58813825553864185693e5{_mm_set1_pd(-5.8813825553864185693e-5)};
const __m128d C67148791796145116483e3{_mm_set1_pd(6.7148791796145116483e-3)};
const __m128d C24869768911673215212e2{_mm_set1_pd(2.4869768911673215212e-2)};
const __m128d C01226618159032020501{_mm_set1_pd(0.1226618159032020501)};
const __m128d C03453457447382168695{_mm_set1_pd(0.3453457447382168695)};
const __m128d C05002494548377929861{_mm_set1_pd(0.5002494548377929861)};
const __m128i C0x7FF0000000000000{_mm_set1_epi64x(0x7FF0000000000000)};
const __m128i C0xFFFFFFFFFFFFF{_mm_set1_epi64x(0xFFFFFFFFFFFFF)};
const __m128i C0x3FF0000000000000{_mm_set1_epi64x(0x3FF0000000000000)};
////////////////////////////////////////////////////////////////////////////
union {__m128i v_i; __m128d v_f;} xu,xu2;
__m128d  val2,val3,val4,b;
__m128i val4i;
val2    = _mm_fmadd_pd(C6497320848556798092,val,C46071824188000174080);
const __mmask8 val2_lt_exp_cst1 = _mm_cmp_pd_mask(val2,exp_cst1,_CMP_LT_OQ);
val3    = _mm_mask_mov_pd(exp_cst1,val2_lt_exp_cst1,val2);
const __mmask8 val3_gt_exp_cst2 = _mm_cmp_pd_mask(val3,exp_cst2,_CMP_GT_OQ);
val4    = _mm_mask_mov_pd(exp_cst2,val3_gt_exp_cst2,val3);
val4i   = _mm_cvtpd_epi64(val4);
xu.v_i  = _mm_and_si128(val4i,C0x7FF0000000000000);
xu2.v_i = _mm_or_si128(_mm_and_si128(val4i,C0xFFFFFFFFFFFFF),C0x3FF0000000000000);
b       = xu2.v_f;
__m128d factor = _mm_fmadd_pd(b,C217150255054231565039e4,CN58813825553864185693e5);
__m128d horner_scheme = _mm_fmadd_pd(
                           _mm_fmadd_pd(
                              _mm_fmadd_pd(
                                 _mm_fmadd_pd(
                                    _mm_fmadd_pd(factor,b,C67148791796145116483e3),b,C24869768911673215212e2),
                                                        b,C01226618159032020501),b,C03453457447382168695),b,C05002494548377929861);
return (_mm_mul_pd(xu.v_f,horner_scheme));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_ALWAYS_INLINE__
static inline 
__m128 
simd_fast_cos_approx_4xf32(const __m128 val)
{
const __m128 C18929864824e5{_mm_set1_ps(1.8929864824e-5f)};
const __m128 CN13422947025e3f{_mm_set1_ps(-1.3422947025e-3f)};
const __m128 CN41518035216e2f{_mm_set1_ps(4.1518035216e-2f)};
const __m128 CN04998515820f{_mm_set1_ps(-0.4998515820f)};
const __m128 C1{_mm_set1_ps(1.0f)};
const __m128 val2 = _mm_mul_ps(val,val);
__m128 cos_horner_scheme = _mm_fmadd_ps(
                          _mm_fmadd_ps(
                             _mm_fmadd_ps(
                                _mm_fmadd_ps(C18929864824e5,val2,CN13422947025e3f),
                                                            val2,CN41518035216e2f),val2,CN04998515820f),val2,C1);
return (cos_horner_scheme);
}

} // math


} // gms



#endif /*__GMS_FAST_SIMD_FUNCS_APPROX__*/