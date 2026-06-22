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
#include <cstdio>
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
simd_fast_log_approx_4xf32(const __m128 val) 
{
const __m128  C06931471805{_mm_set1_ps(0.6931471805f)};
const __m128  C3110401639em2{_mm_set1_ps(3.110401639e-2f)};
const __m128  CN0288739945{_mm_set1_ps(-0.288739945f)};
const __m128  C1130626167{_mm_set1_ps(1.130626167f)};
const __m128  CN2461222105{_mm_set1_ps(-2.461222105f)};
const __m128  C3529304993{_mm_set1_ps(3.529304993f)};
const __m128  CN89970756366{_mm_set1_ps(-89.970756366f)};
const __m128  CNINFINITY{_mm_set1_ps(-(float)INFINITY)};
const __m128i C0x7FFFFF{_mm_set1_epi32(0x7FFFFF)};
const __m128i C0x3F800000{_mm_set1_epi32(0x3F800000)};
///////////////////////////////////////////////////////
union {__m128 v_f; __m128i v_i;} valu;
__m128 exp,addcst,x;
valu.v_f = val;
__mmask8 val_gt_0 = _mm_cmp_ps_mask(val,_mm_setzero_ps(),_CMP_GT_OQ);
exp      = _mm_cvtepi32_ps(_mm_srai_epi32(valu.v_i,23));
addcst   = _mm_mask_mov_ps(CNINFINITY,val_gt_0,CN89970756366);
valu.v_i = _mm_or_si128(_mm_and_si128(valu.v_i,C0x7FFFFF),C0x3F800000);
x        = valu.v_f;
__m128 factor = _mm_fmadd_ps(C06931471805,exp,addcst);
__m128 horner_scheme = _mm_fmadd_ps(
                          _mm_fmadd_ps(
                             _mm_fmadd_ps(
                                _mm_fmadd_ps(C3110401639em2,x,CN0288739945),
                                                            x,C1130626167),x,CN2461222105),x,C3529304993);
return (_mm_fmadd_ps(x,horner_scheme,factor));
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
simd_fast_log_approx_2xf64(const __m128d val) 
{
const __m128d C0693147180559945309{_mm_set1_pd(0.693147180559945309)};
const __m128d C38452996567427771004em3{_mm_set1_pd(3.8452996567427771004e-3)};
const __m128d CN58063046913513248806em2{_mm_set1_pd(-5.8063046913513248806e-2)};
const __m128d C03933681393710522449{_mm_set1_pd(0.3933681393710522449)};
const __m128d CN15768742873140365359{_mm_set1_pd(-1.5768742873140365359)};
const __m128d C41526922366643242018{_mm_set1_pd(4.1526922366643242018)};
const __m128d CN75448441898968122224{_mm_set1_pd(-7.5448441898968122224)};
const __m128d C96834305233551590662{_mm_set1_pd(9.6834305233551590662)};
const __m128d CN89283961901995856725{_mm_set1_pd(-8.9283961901995856725)};
const __m128d C63599808385417079442{_mm_set1_pd(6.3599808385417079442)};
const __m128d CN71157470503155781879{_mm_set1_pd(-711.57470503155781879)};
const __m128d CNINFINITY{_mm_set1_pd(-(double)INFINITY)};
const __m128i C0xFFFFFFFFFFFFF{_mm_set1_epi64x(0xFFFFFFFFFFFFF)};
const __m128i C0x3FF0000000000000{_mm_set1_epi64x(0x3FF0000000000000)};
/////////////////////////////////////////////////////////////////////////
union {__m128d v_f; __m128i v_i;} valu;
__m128d exp,addcst,x;
valu.v_f = val;
exp      = _mm_cvtepi64_pd(_mm_srai_epi64(valu.v_i,52));
__mmask8 val_gt_0 = _mm_cmp_pd_mask(val,_mm_setzero_pd(),_CMP_GT_OQ);
valu.v_i = _mm_or_epi64(_mm_and_si128(valu.v_i,C0xFFFFFFFFFFFFF),C0x3FF0000000000000);
addcst   = _mm_mask_mov_pd(CNINFINITY,val_gt_0,CN71157470503155781879);
x        = valu.v_f;
__m128d factor{_mm_fmadd_pd(C0693147180559945309,exp,addcst)};
__m128d horner_scheme = _mm_fmadd_pd(
                           _mm_fmadd_pd(
                              _mm_fmadd_pd(
                                 _mm_fmadd_pd(
                                    _mm_fmadd_pd(
                                       _mm_fmadd_pd(
                                          _mm_fmadd_pd(
                                              _mm_fmadd_pd(C38452996567427771004em3,x,CN58063046913513248806em2),x,C03933681393710522449),
                                                                  x,CN15768742873140365359),x,C41526922366643242018),x,CN75448441898968122224),
                                                                           x,C96834305233551590662),x,CN89283961901995856725),x,C63599808385417079442);
return (_mm_fmadd_pd(x,horner_scheme,factor));
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
   /* Correct only in [-pi, pi]
   Absolute error bounded by 5e-5
   Continuous error */
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

/*
    return 1. + val2 * (-0.4999999989854595 + val2 * (4.1666664031946494e-2 +
               val2 * (-1.3888865732154298e-3 + val2 * (2.4800621164404954e-5 +
                val2 * (-2.7535709220609759e-7 + val2 * (2.0609775337948301e-9 +
                 val2 * (-9.7393350917248365e-12)))))));
*/

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
simd_fast_cos_approx_2xf64(const __m128d val)
{
   /* Correct only in [-pi, pi]
   Absolute error bounded by 2e-10
   Continuous error */
const __m128d CN97393350917248365em12{_mm_set1_pd(-9.7393350917248365e-12)};
const __m128d C20609775337948301em9{2.0609775337948301e-9};
const __m128d CN27535709220609759em7{_mm_set1_pd(-2.7535709220609759e-7)};
const __m128d C24800621164404954em5{_mm_set1_pd(2.4800621164404954e-5)};
const __m128d CN13888865732154298em3{_mm_set1_pd(-1.3888865732154298e-3)};
const __m128d C41666664031946494em2{_mm_set1_pd(4.1666664031946494e-2)};
const __m128d CN04999999989854595{_mm_set1_pd(-0.4999999989854595)};
const __m128d C1{_mm_set1_pd(1.0)};

__m128d val2 = _mm_mul_pd(val,val);
__m128d cos_horner_scheme = _mm_fmadd_pd(
                               _mm_fmadd_pd(
                                  _mm_fmadd_pd(
                                     _mm_fmadd_pd(
                                        _mm_fmadd_pd(
                                           _mm_fmadd_pd(
                                              _mm_fmadd_pd(CN97393350917248365em12,val2,C20609775337948301em9),
                                                               val2,CN27535709220609759em7),val2,C24800621164404954em5),
                                                                     val2,CN13888865732154298em3),val2,C41666664031946494em2),
                                                                                                val2,CN04999999989854595),val2,C1);
return (cos_horner_scheme);                                                                        
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
simd_fast_sin_approx_4xf32(const __m128 val)
{
   /* Correct only in [-pi, pi]
   Absolute error bounded by 6e-6
   Continuous error */
const __m128 C21478401777em6{_mm_set1_ps(2.1478401777e-6f)};
const __m128 CN19264918228em4{_mm_set1_ps(-1.9264918228e-4f)};
const __m128 C83089787513em3{_mm_set1_ps(8.3089787513e-3f)};
const __m128 CN01666243672{_mm_set1_ps(-0.1666243672f)};
const __m128 C09999793767{_mm_set1_ps(0.9999793767f)};
__m128 val2 = _mm_mul_ps(val,val);
__m128 sin_horner_scheme = _mm_fmadd_ps(
                              _mm_fmadd_ps(
                                 _mm_fmadd_ps(
                                    _mm_fmadd_ps(C21478401777em6,val2,CN19264918228em4),
                                                     val2,C83089787513em3),val2,CN01666243672),
                                                                               val2,C09999793767);
return (_mm_mul_ps(val,sin_horner_scheme));
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
simd_fast_sin_approx_2xf64(const __m128d val)
{
   /* Correct only in [-pi, pi]
   Absolute error bounded by 2e-9
   Continuous error */
__m128d val2 = _mm_mul_pd(val,val);
const __m128d C0{_mm_set1_pd(1.3451481340051383601e-10)};
const __m128d C1{_mm_set1_pd(-2.4676970823046321831e-8)};
const __m128d C2{_mm_set1_pd(2.7529454331962521774e-6)};
const __m128d C3{_mm_set1_pd(-1.9840155355055654144e-4)};
const __m128d C4{_mm_set1_pd(8.3333103922589284663e-3)};
const __m128d C5{_mm_set1_pd(-0.1666666458182987439)};
const __m128d C6{_mm_set1_pd(0.9999999945159759653)};
__m128d sin_horner_scheme = _mm_fmadd_pd(_mm_fmadd_pd(_mm_fmadd_pd(_mm_fmadd_pd(
                                     _mm_fmadd_pd(_mm_fmadd_pd(C0,val2,C1),val2,C2),val2,C3),val2,C4),val2,C5),val2,C6);
return (_mm_mul_pd(val,sin_horner_scheme));
}


/*
    Original scalar implementation, used mainly or the scalar loop remainder.
*/

/* Workaround a lack of optimization in gcc */


/* Relative error bounded by 1e-5 for normalized outputs
   Returns invalid outputs for nan inputs
   Continuous error */
__ATTR_ALWAYS_INLINE__
static inline
float expapprox(float val) {
  union { int32_t i; float f; } xu, xu2;
  float exp_cst1_f = 2139095040.f;
  float exp_cst2_f = 0.f;
  float val2, val3, val4, b;
  int32_t val4i;
  val2 = 12102203.1615614f*val+1065353216.f;
  val3 = val2 < exp_cst1_f ? val2 : exp_cst1_f;
  val4 = val3 > exp_cst2_f ? val3 : exp_cst2_f;
  val4i = (int32_t) val4;
  xu.i = val4i & 0x7F800000;
  xu2.i = (val4i & 0x7FFFFF) | 0x3F800000;
  b = xu2.f;

  /* Generated in Sollya with:
     > f=remez(1-x*exp(-(x-1)*log(2)),
               [|(x-1)*(x-2), (x-1)*(x-2)*x, (x-1)*(x-2)*x*x|],
               [1.000001,1.999999], exp(-(x-1)*log(2)));
     > plot(exp((x-1)*log(2))/(f+x)-1, [1,2]);
     > f+x;
  */
  return
    xu.f * (0.509871020f + b * (0.312146713f + b * (0.166617139f + b *
              (-2.190619930e-3f + b * 1.3555747234e-2f))));
}

/* Same code with better precision, less performance and using double */


/* Relative error bounded by 3e-9 for normalized outputs
   Returns invalid outputs for nan inputs
   Continuous error

   Vectorizable only with AVX512dq extensions because of the
   double->int64 cast. On GCC, use option -mavx512dq. */
__ATTR_ALWAYS_INLINE__
static inline
double expapprox_d(double val) {
  union { int64_t i; double f; } xu, xu2;
  double exp_cst1_d = 9218868437227405312.;
  double exp_cst2_d = 0.;
  double val2, val3, val4, b;
  int64_t val4i;
  val2 = 6497320848556798.092*val+4607182418800017408.;
  val3 = val2 < exp_cst1_d ? val2 : exp_cst1_d;
  val4 = val3 > exp_cst2_d ? val3 : exp_cst2_d;
  val4i = (int64_t) val4;
  xu.i = val4i & 0x7FF0000000000000;
  xu2.i = (val4i & 0xFFFFFFFFFFFFF) | 0x3FF0000000000000;
  b = xu2.f;

  /* Generated in Sollya with:
     > f=remez(1-x*exp(-(x-1)*log(2)),
               [|(x-1)*(x-2), (x-1)*(x-2)*x, (x-1)*(x-2)*x*x,
                 (x-1)*(x-2)*x*x*x, (x-1)*(x-2)*x*x*x*x|],
                 [1.000001,1.999999], exp(-(x-1)*log(2)));
     > plot(exp((x-1)*log(2))/(f+x)-1, [1,2]);
     > f+x;
  */
  return
    xu.f * (0.5002494548377929861 + b * (0.3453457447382168695 + b *
            (0.1226618159032020501 + b * (2.4869768911673215212e-2 + b *
             (6.7148791796145116483e-3 + b * (-5.8813825553864185693e-5 + b *
              2.17150255054231565039e-4))))));
}

/* Absolute error bounded by 1e-5 for normalized inputs
   Returns a finite number for +inf input
   Returns -inf for nan and <= 0 inputs.
   Continuous error. */
__ATTR_ALWAYS_INLINE__
static inline
float logapprox(float val) {
  union { float f; int32_t i; } valu;
  float exp, addcst, x;
  valu.f = val;
  exp = valu.i >> 23;
  /* -89.970756366f = -127 * log(2) + constant term of polynomial bellow. */
  addcst = val > 0 ? -89.970756366f : -(float)INFINITY;
  valu.i = (valu.i & 0x7FFFFF) | 0x3F800000;
  x = valu.f;


  /* Generated in Sollya using:
    > f = remez(log(x)-(x-1)*log(2),
            [|1,(x-1)*(x-2), (x-1)*(x-2)*x, (x-1)*(x-2)*x*x,
              (x-1)*(x-2)*x*x*x|], [1,2], 1, 1e-8);
    > plot(f+(x-1)*log(2)-log(x), [1,2]);
    > f+(x-1)*log(2)
 */
  return
    x * (3.529304993f + x * (-2.461222105f + x * (1.130626167f +
      x * (-0.288739945f + x * 3.110401639e-2f))))
    + (addcst + 0.6931471805f*exp);
}

/* Absolute error bounded by 5e-9 for normalized inputs
   Returns a finite number for +inf input
   Returns -inf for nan and <= 0 inputs.
   Continuous error.

   Vectorizable only with AVX512dq extensions because of the
   int64 left shift. On GCC, use option -mavx512dq. */
__ATTR_ALWAYS_INLINE__
static inline
double logapprox_d(double val) {
  union { double f; int64_t i; } valu;
  double exp, addcst, x;
  valu.f = val;
  exp = valu.i >> 52;
  /* -711.57470503155781879 = -1023 * log(2) + constant term of polynomial bellow. */
  addcst = val > 0 ? -711.57470503155781879 : -(double)INFINITY;
  valu.i = (valu.i & 0xFFFFFFFFFFFFF) | 0x3FF0000000000000;
  x = valu.f;


  /* Generated in Sollya using:
    > f = remez(log(x)-(x-1)*log(2),
            [|1,(x-1)*(x-2), (x-1)*(x-2)*x, (x-1)*(x-2)*x*x,
              (x-1)*(x-2)*x*x*x|], [1,2], 1, 1e-8);
    > plot(f+(x-1)*log(2)-log(x), [1,2]);
    > f+(x-1)*log(2)
 */
  return x * (6.3599808385417079442 + x * (-8.9283961901995856725 +
          x * (9.6834305233551590662 + x * (-7.5448441898968122224 +
           x * (4.1526922366643242018 + x * (-1.5768742873140365359 +
            x * (0.3933681393710522449 + x * (-5.8063046913513248806e-2 +
             x * 3.8452996567427771004e-3))))))))
    + (addcst + 0.693147180559945309*exp);
}

/* Correct only in [-pi, pi]
   Absolute error bounded by 5e-5
   Continuous error */
__ATTR_ALWAYS_INLINE__
static inline
float cosapprox(float val) {
  float val2 = val*val;
  /* Generated in Sollya using:
     > f = remez(cos(x)-1, [|x*x, x*x*x*x, x*x*x*x*x*x, x*x*x*x*x*x*x*x|],
                            [0.000001, pi], 1, 1e-8);
     > plot(f-cos(x)+1, [0, pi]);
     > f+1
  */
  return
    1.f + val2 * (-0.4998515820f + val2 * (4.1518035216e-2f + val2 *
      (-1.3422947025e-3f + val2 * 1.8929864824e-5f)));
}

/* Correct only in [-pi, pi]
   Absolute error bounded by 2e-10
   Continuous error */
__ATTR_ALWAYS_INLINE__
static inline
double cosapprox_d(double val) {
  double val2 = val*val;
  /* Generated in Sollya using:
     > f = remez(cos(x)-1, [|x*x, x*x*x*x, x*x*x*x*x*x, x*x*x*x*x*x*x*x,
                             x*x*x*x*x*x*x*x*x*x, x*x*x*x*x*x*x*x*x*x*x*x,
                             x*x*x*x*x*x*x*x*x*x*x*x*x*x|], [1e-10, pi], 1, 1e-10);
     > plot(f-cos(x)+1, [0, pi]);
     > f+1
 */
  return 1. + val2 * (-0.4999999989854595 + val2 * (4.1666664031946494e-2 +
               val2 * (-1.3888865732154298e-3 + val2 * (2.4800621164404954e-5 +
                val2 * (-2.7535709220609759e-7 + val2 * (2.0609775337948301e-9 +
                 val2 * (-9.7393350917248365e-12)))))));
}

/* Correct only in [-pi, pi]
   Absolute error bounded by 6e-6
   Continuous error */
__ATTR_ALWAYS_INLINE__
static inline
float sinapprox(float val) {
  float val2 = val*val;
  return
    val * (0.9999793767f + val2 * (-0.1666243672f + val2 *
          (8.3089787513e-3f + val2 * (-1.9264918228e-4f + val2 *
           2.1478401777e-6f))));
}

/* Correct only in [-pi, pi]
   Absolute error bounded by 2e-9
   Continuous error */
__ATTR_ALWAYS_INLINE__
static inline
double sinapprox_d(double val) {
  double val2 = val*val;
  return
    val * (0.9999999945159759653 + val2 * (-0.1666666458182987439 +
     val2 * (8.3333103922589284663e-3 + val2 * (-1.9840155355055654144e-4 +
      val2 * (2.7529454331962521774e-6 + val2 * (-2.4676970823046321831e-8 +
       val2 * 1.3451481340051383601e-10))))));
}


} // math


} // gms



#endif /*__GMS_FAST_SIMD_FUNCS_APPROX__*/