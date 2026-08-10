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

#ifndef __GMS_INTEGRANDS_FUNC_CH4_SSE_H__
#define __GMS_INTEGRANDS_FUNC_CH4_SSE_H__ 040820261033

#include <cstdint>
#include <immintrin.h>
#include <cmath>
#include "GMS_config.h"


namespace file_info 
{

     static const unsigned int GMS_INTEGRANDS_FUNC_CH4_SSE_MAJOR = 1;
     static const unsigned int GMS_INTEGRANDS_FUNC_CH4_SSE_MINOR = 1;
     static const unsigned int GMS_INTEGRANDS_FUNC_CH4_SSE_MICRO = 0;
     static const unsigned int GMS_INTEGRANDS_FUNC_CH4_SSE_FULLVER =
       1000U*GMS_INTEGRANDS_FUNC_CH4_SSE_MAJOR+100U*GMS_INTEGRANDS_FUNC_CH4_SSE_MINOR+
       10U*GMS_INTEGRANDS_FUNC_CH4_SSE_MICRO;
     static const char GMS_INTEGRANDS_FUNC_CH4_SSE_CREATION_DATE[] = "04-08-2026 10:33AM +00200 (TUE 04 AUG 2026 GMT+2)";
     static const char GMS_INTEGRANDS_FUNC_CH4_SSE_BUILD_DATE[]    = __DATE__; 
     static const char GMS_INTEGRANDS_FUNC_CH4_SSE_BUILD_TIME[]    = __TIME__;
     static const char GMS_INTEGRANDS_FUNC_CH4_SSE_SYNOPSIS[]      = "Integrands of the functional form SSE-accelerated used for fading channel and BER,BEP,SER,SEP quadrature computation(QUADPACK).\
	                                                          Based on the M.K Simon, M.S. Alouini: Digital Communication over Fading Channels 1st ed\
															  ISBN-13 978-0471317791";

}


#if !defined(INTEGRANDS_FUNC_CH4_SSE_DO_WARMUP_CALL)
#define INTEGRANDS_FUNC_CH4_SSE_DO_WARMUP_CALL 1
#endif

#if !defined(INTEGRANDS_FUNC_CH4_SSE_DEFINE_FUNC_BODY)
#define INTEGRANDS_FUNC_CH4_SSE_DEFINE_FUNC_BODY 1
#endif 

namespace gms 
{

namespace fading_channel
{

namespace 
{
  /*
      Workaround for the weird ICPC error.
  */
const __m128d  NZ128DP = _mm_set1_pd(-0.0);
__ATTR_ALWAYS_INLINE__		     
static inline 
__m128d
negate_xmm2r8(const __m128d v) 
{
	return (_mm_xor_pd(v,NZ128DP));
}
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
integrand_4_1_gauss_Q_func_sse_pd(const __m128d y)
{
  const __m128d vC03989422804014326779399460599344{_mm_set1_pd(0.3989422804014326779399460599344)};
  __m128d exp_arg{_mm_mul_pd(_mm_set1_pd(0.5),_mm_mul_pd(y,y))};
  __m128d exp_val{_mm_exp_pd(negate_xmm2r8(exp_arg))};
  return (_mm_mul_pd(vC03989422804014326779399460599344,exp_val));
}

/* Parametrized macro block for insertion into the loop 
   body*/

#if (INTEGRANDS_FUNC_CH4_SSE_DEFINE_FUNC_BODY) == 1

#define INTEGRAND_4_1_GAUSS_FUNC_SSE_BODY(result,y)\
{\
const __m128d vC03989422804014326779399460599344{_mm_set1_pd(0.3989422804014326779399460599344)};\
__m128d exp_arg{_mm_mul_pd(_mm_set1_pd(0.5),_mm_mul_pd((y),(y)))};\
__m128d exp_val{_mm_exp_pd(negate_xmm2r8(exp_arg))};\
(result) = _mm_mul_pd(vC03989422804014326779399460599344,exp_val);\
}

#endif

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
integrand_4_2_gauss_Q_func_sse_pd(const __m128d x,
                                  const __m128d theta)
{
  const __m128d xsqr{_mm_mul_pd(x,x)};
  const __m128d sin_val{_mm_sin_pd(theta)};
  const __m128d tmp_sin_mul{_mm_mul_pd(sin_val,sin_val)};
  const __m128d sin_denom{_mm_add_pd(tmp_sin_mul,tmp_sin_mul)};
  const __m128d exp_arg{_mm_div_pd(xsqr,sin_denom)};
  const __m128d exp_val{_mm_exp_pd(negate_xmm2r8(exp_arg))};
  return (exp_val);
}

#if (INTEGRANDS_FUNC_CH4_SSE_DEFINE_FUNC_BODY) == 1

#define INTEGRAND_4_2_GAUSS_Q_FUNC_SSE_BODY(result,x,theta)\
{\
const __m128d xsqr{_mm_mul_pd((x),(x))};\
const __m128d sin_val{_mm_sin_pd((theta))};\
const __m128d tmp_sin_mul{_mm_mul_pd(sin_val,sin_val)};\
const __m128d sin_denom{_mm_add_pd(tmp_sin_mul,tmp_sin_mul)};\
const __m128d exp_arg{_mm_div_pd(xsqr,sin_denom)};\
const __m128d exp_val{_mm_exp_pd(negate_xmm2r8(exp_arg))};\
(result) = exp_val;\
}

#endif 

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
integrand_4_6_sin_gauss_Q_func_sse_pd(const __m128d x1,const __m128d y1,
                                      const __m128d rho,const __m128d theta,
                                      const __m128d phi_s /*an upper limit of integration, computed one per integration full step*/)
{
  const __m128d one{_mm_set1_pd(1.0)};
  const __m128d squared_sum{_mm_fmadd_pd(x1,x1,_mm_mul_pd(y1,y1))};
  const __m128d tsin{_mm_sin_pd(theta)};
  const __m128d sin_sqr_tht{_mm_mul_pd(tsin,tsin)};
  const __m128d S_hat{_mm_sqrt_pd(squared_sum)};
  const __m128d one_m_rho{_mm_sub_pd(one,_mm_mul_pd(rho,rho))};
  const __m128d sin2theta{_mm_sin_pd(_mm_add_pd(theta,theta))};
  const __m128d one_m_sin2theta{_mm_sub_pd(one,_mm_mul_pd(rho,sin2theta))};
  const __m128d sqr_1m_rho{_mm_sqrt_pd(one_m_rho)};
  const __m128d S_hat_half{_mm_mul_pd(_mm_set1_pd(0.5),_mm_mul_pd(S_hat,S_hat))};
  const __m128d ratio1{_mm_div_pd(sqr_1m_rho,one_m_sin2theta)};
  const __m128d sin_phi_s{_mm_sin_pd(phi_s)};
  const __m128d sin_phi_s_sqr{_mm_mul_pd(sin_phi_s,sin_phi_s)};
  const __m128d ratio_sin{_mm_div_pd(sin_phi_s_sqr,sin_sqr_tht)};
  const __m128d ratio2{_mm_div_pd(one_m_sin2theta,one_m_rho)};
  const __m128d exp_arg{_mm_mul_pd(S_hat_half,_mm_mul_pd(ratio2,ratio_sin))};
  const __m128d exp_val{_mm_exp_pd(negate_xmm2r8(exp_arg))};
  return (_mm_mul_pd(ratio1,exp_val));
}

#if (INTEGRANDS_FUNC_CH4_SSE_DEFINE_FUNC_BODY) == 1

#define INTEGRAND_4_6_SIN_GAUSS_Q_FUNC_SSE_BODY(result,x1,y1,rho,theta,phi_s)\
{\
  const __m128d one{_mm_set1_pd(1.0)};\
  const __m128d squared_sum{_mm_fmadd_pd((x1),(x1),_mm_mul_pd((y1),(y1)))};\
  const __m128d tsin{_mm_sin_pd((theta))};\
  const __m128d sin_sqr_tht{_mm_mul_pd(tsin,tsin)};\
  const __m128d S_hat{_mm_sqrt_pd(squared_sum)};\
  const __m128d one_m_rho{_mm_sub_pd(one,_mm_mul_pd((rho),(rho)))};\
  const __m128d sin2theta{_mm_sin_pd(_mm_add_pd((theta),(theta)))};\
  const __m128d one_m_sin2theta{_mm_sub_pd(one,_mm_mul_pd(rho,sin2theta))};\
  const __m128d sqr_1m_rho{_mm_sqrt_pd(one_m_rho)};\
  const __m128d S_hat_half{_mm_mul_pd(_mm_set1_pd(0.5),_mm_mul_pd(S_hat,S_hat))};\
  const __m128d ratio1{_mm_div_pd(sqr_1m_rho,one_m_sin2theta)};\
  const __m128d sin_phi_s{_mm_sin_pd((phi_s))};\
  const __m128d sin_phi_s_sqr{_mm_mul_pd(sin_phi_s,sin_phi_s)};\
  const __m128d ratio_sin{_mm_div_pd(sin_phi_s_sqr,sin_sqr_tht)};\
  const __m128d ratio2{_mm_div_pd(one_m_sin2theta,one_m_rho)};\
  const __m128d exp_arg{_mm_mul_pd(S_hat_half,_mm_mul_pd(ratio2,ratio_sin))};\
  const __m128d exp_val{_mm_exp_pd(negate_xmm2r8(exp_arg))};\
  (result) = _mm_mul_pd(ratio1,exp_val);\
}

#endif 

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
integrand_4_6_cos_gauss_Q_func_sse_pd(const __m128d x1,const __m128d y1,
                                      const __m128d rho,const __m128d theta,
                                      const __m128d phi_s /*an upper limit of integration, computed one per integration full step*/)
{
  const __m128d one{_mm_set1_pd(1.0)};
  const __m128d squared_sum{_mm_fmadd_pd(x1,x1,_mm_mul_pd(y1,y1))};
  const __m128d tsin{_mm_sin_pd(theta)};
  const __m128d sin_sqr_tht{_mm_mul_pd(tsin,tsin)};
  const __m128d S_hat{_mm_sqrt_pd(squared_sum)};
  const __m128d one_m_rho{_mm_sub_pd(one,_mm_mul_pd(rho,rho))};
  const __m128d sin2theta{_mm_sin_pd(_mm_add_pd(theta,theta))};
  const __m128d one_m_sin2theta{_mm_sub_pd(one,_mm_mul_pd(rho,sin2theta))};
  const __m128d sqr_1m_rho{_mm_sqrt_pd(one_m_rho)};
  const __m128d S_hat_half{_mm_mul_pd(_mm_set1_pd(0.5),_mm_mul_pd(S_hat,S_hat))};
  const __m128d ratio1{_mm_div_pd(sqr_1m_rho,one_m_sin2theta)};
  const __m128d cos_phi_s{_mm_cos_pd(phi_s)};
  const __m128d cos_phi_s_sqr{_mm_mul_pd(cos_phi_s,cos_phi_s)};
  const __m128d ratio_cos{_mm_div_pd(cos_phi_s_sqr,sin_sqr_tht)};
  const __m128d ratio2{_mm_div_pd(one_m_sin2theta,one_m_rho)};
  const __m128d exp_arg{_mm_mul_pd(S_hat_half,_mm_mul_pd(ratio2,ratio_cos))};
  const __m128d exp_val{_mm_exp_pd(negate_xmm2r8(exp_arg))};
  return (_mm_mul_pd(ratio1,exp_val));
}

#if (INTEGRANDS_FUNC_CH4_SSE_DEFINE_FUNC_BODY) == 1

#define INTEGRAND_4_6_COS_GAUSS_Q_FUNC_SSE_BODY(result,x1,y1,rho,theta,phi_s)\
{\
  const __m128d one{_mm_set1_pd(1.0)};\
  const __m128d squared_sum{_mm_fmadd_pd((x1),(x1),_mm_mul_pd((y1),(y1)))};\
  const __m128d tsin{_mm_sin_pd((theta))};\
  const __m128d sin_sqr_tht{_mm_mul_pd(tsin,tsin)};\
  const __m128d S_hat{_mm_sqrt_pd(squared_sum)};\
  const __m128d one_m_rho{_mm_sub_pd(one,_mm_mul_pd((rho),(rho)))};\
  const __m128d sin2theta{_mm_sin_pd(_mm_add_pd((theta),(theta)))};\
  const __m128d one_m_sin2theta{_mm_sub_pd(one,_mm_mul_pd((rho),sin2theta))};\
  const __m128d sqr_1m_rho{_mm_sqrt_pd(one_m_rho)};\
  const __m128d S_hat_half{_mm_mul_pd(_mm_set1_pd(0.5),_mm_mul_pd(S_hat,S_hat))};\
  const __m128d ratio1{_mm_div_pd(sqr_1m_rho,one_m_sin2theta)};\
  const __m128d cos_phi_s{_mm_cos_pd((phi_s))};\
  const __m128d cos_phi_s_sqr{_mm_mul_pd(cos_phi_s,cos_phi_s)};\
  const __m128d ratio_cos{_mm_div_pd(cos_phi_s_sqr,sin_sqr_tht)};\
  const __m128d ratio2{_mm_div_pd(one_m_sin2theta,one_m_rho)};\
  const __m128d exp_arg{_mm_mul_pd(S_hat_half,_mm_mul_pd(ratio2,ratio_cos))};\
  const __m128d exp_val{_mm_exp_pd(negate_xmm2r8(exp_arg))};\
  (result) = _mm_mul_pd(ratio1,exp_val);\
}

#endif 

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
integrand_4_7_x1_gauss_Q_func_sse_pd(const __m128d x1,
                                      const __m128d rho,const __m128d theta)
{
  const __m128d one{_mm_set1_pd(1.0)};
  const __m128d half{_mm_set1_pd(0.5)};
  const __m128d rho_m_sin2tht{_mm_mul_pd(rho,_mm_sin_pd(_mm_add_pd(theta,theta)))};
  const __m128d one_m_rho_sqr{_mm_sub_pd(one,_mm_mul_pd(rho,rho))};
  const __m128d tmp_sin{_mm_sin_pd(theta)};
  const __m128d sqr1m_rhosqr{_mm_sqrt_pd(one_m_rho_sqr)};
  const __m128d one_m_rhosintht{_mm_sub_pd(one,rho_m_sin2tht)};
  const __m128d sinthtpow2{_mm_mul_pd(tmp_sin,tmp_sin)};
  const __m128d lead_factor{_mm_div_pd(sqr1m_rhosqr,one_m_rhosintht)};
  const __m128d half_x1{_mm_mul_pd(half,_mm_mul_pd(x1,x1))};
  const __m128d ratio{_mm_div_pd(one_m_rhosintht,_mm_mul_pd(one_m_rho_sqr,sinthtpow2))};
  const __m128d exp_arg{_mm_mul_pd(half_x1,ratio)};
  const __m128d exp_val{_mm_exp_pd(negate_xmm2r8(exp_arg))};
  return (_mm_mul_pd(lead_factor,exp_val));
}

#if (INTEGRANDS_FUNC_CH4_SSE_DEFINE_FUNC_BODY) == 1

#define INTEGRAND_4_7_X1_GAUSS_Q_FUNC_SSE_BODY(result,x1,rho,theta)\
{\
  const __m128d one{_mm_set1_pd(1.0)};\
  const __m128d half{_mm_set1_pd(0.5)};\
  const __m128d rho_m_sin2tht{_mm_mul_pd(rho,_mm_sin_pd(_mm_add_pd((theta),(theta))))};\
  const __m128d one_m_rho_sqr{_mm_sub_pd(one,_mm_mul_pd((rho),(rho)))};\
  const __m128d tmp_sin{_mm_sin_pd((theta))};\
  const __m128d sqr1m_rhosqr{_mm_sqrt_pd(one_m_rho_sqr)};\
  const __m128d one_m_rhosintht{_mm_sub_pd(one,rho_m_sin2tht)};\
  const __m128d sinthtpow2{_mm_mul_pd(tmp_sin,tmp_sin)};\
  const __m128d lead_factor{_mm_div_pd(sqr1m_rhosqr,one_m_rhosintht)};\
  const __m128d half_x1{_mm_mul_pd(half,_mm_mul_pd((x1),(x1)))};\
  const __m128d ratio{_mm_div_pd(one_m_rhosintht,_mm_mul_pd(one_m_rho_sqr,sinthtpow2))};\
  const __m128d exp_arg{_mm_mul_pd(half_x1,ratio)};\
  const __m128d exp_val{_mm_exp_pd(negate_xmm2r8(exp_arg))};\
  (result) = _mm_mul_pd(lead_factor,exp_val);\
}

#endif 

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
integrand_4_7_y1_gauss_Q_func_sse_pd(const __m128d y1,
                                      const __m128d rho,const __m128d theta)
{
  const __m128d one{_mm_set1_pd(1.0)};
  const __m128d half{_mm_set1_pd(0.5)};
  const __m128d rho_m_sin2tht{_mm_mul_pd(rho,_mm_sin_pd(_mm_add_pd(theta,theta)))};
  const __m128d one_m_rho_sqr{_mm_sub_pd(one,_mm_mul_pd(rho,rho))};
  const __m128d tmp_sin{_mm_sin_pd(theta)};
  const __m128d sqr1m_rhosqr{_mm_sqrt_pd(one_m_rho_sqr)};
  const __m128d one_m_rhosintht{_mm_sub_pd(one,rho_m_sin2tht)};
  const __m128d sinthtpow2{_mm_mul_pd(tmp_sin,tmp_sin)};
  const __m128d lead_factor{_mm_div_pd(sqr1m_rhosqr,one_m_rhosintht)};
  const __m128d half_y1{_mm_mul_pd(half,_mm_mul_pd(y1,y1))};
  const __m128d ratio{_mm_div_pd(one_m_rhosintht,_mm_mul_pd(one_m_rho_sqr,sinthtpow2))};
  const __m128d exp_arg{_mm_mul_pd(half_y1,ratio)};
  const __m128d exp_val{_mm_exp_pd(negate_xmm2r8(exp_arg))};
  return (_mm_mul_pd(lead_factor,exp_val));
}

#if (INTEGRANDS_FUNC_CH4_SSE_DEFINE_FUNC_BODY) == 1

#define INTEGRAND_4_7_Y1_GAUSS_Q_FUNC_SSE_BODY(result,y1,rho,theta)\
{\
  const __m128d one{_mm_set1_pd(1.0)};\
  const __m128d half{_mm_set1_pd(0.5)};\
  const __m128d rho_m_sin2tht{_mm_mul_pd(rho,_mm_sin_pd(_mm_add_pd((theta),(theta))))};\
  const __m128d one_m_rho_sqr{_mm_sub_pd(one,_mm_mul_pd((rho),(rho)))};\
  const __m128d tmp_sin{_mm_sin_pd((theta))};\
  const __m128d sqr1m_rhosqr{_mm_sqrt_pd(one_m_rho_sqr)};\
  const __m128d one_m_rhosintht{_mm_sub_pd(one,rho_m_sin2tht)};\
  const __m128d sinthtpow2{_mm_mul_pd(tmp_sin,tmp_sin)};\
  const __m128d lead_factor{_mm_div_pd(sqr1m_rhosqr,one_m_rhosintht)};\
  const __m128d half_y1{_mm_mul_pd(half,_mm_mul_pd((y1),(y1)))};\
  const __m128d ratio{_mm_div_pd(one_m_rhosintht,_mm_mul_pd(one_m_rho_sqr,sinthtpow2))};\
  const __m128d exp_arg{_mm_mul_pd(half_y1,ratio)};\
  const __m128d exp_val{_mm_exp_pd(negate_xmm2r8(exp_arg))};\
  (result) = _mm_mul_pd(lead_factor,exp_val);\
}

#endif 

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
integrand_4_8_x1_gauss_Q_func_sse_pd(const __m128d x1,const __m128d theta)
{
  const __m128d x1_sqr{_mm_mul_pd(x1,x1)};
  const __m128d tmp_sin{_mm_sin_pd(theta)};
  const __m128d two_sinthtsqr{_mm_add_pd(tmp_sin,tmp_sin)};
  const __m128d exp_arg{_mm_div_pd(x1_sqr,two_sinthtsqr)};
  const __m128d exp_val{_mm_exp_pd(negate_xmm2r8(exp_arg))};
  return (exp_val);
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
integrand_4_8_y1_gauss_Q_func_sse_pd(const __m128d y1,const __m128d theta)
{
  const __m128d y1_sqr{_mm_mul_pd(y1,y1)};
  const __m128d tmp_sin{_mm_sin_pd(theta)};
  const __m128d two_sinthtsqr{_mm_add_pd(tmp_sin,tmp_sin)};
  const __m128d exp_arg{_mm_div_pd(y1_sqr,two_sinthtsqr)};
  const __m128d exp_val{_mm_exp_pd(negate_xmm2r8(exp_arg))};
  return (exp_val);
}

#include "GMS_bessel_i0_sse.h"

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
integrand_4_10_marcum_Q_func_sse_pd(const __m128d x,const __m128d s)
{
  using namespace gms::math;
  const __m128d half{_mm_set1_pd(0.5)};
  const __m128d C15{_mm_set1_pd(15)};
  const __m128d sx{_mm_mul_pd(s,x)};
  const __mmask8 is_sx_ge15{_mm_cmp_pd_mask(sx,C15,_CMP_GE_OQ)};
  const __m128d xspow2{_mm_fmadd_pd(x,x,_mm_mul_pd(s,s))};
  const __m128d besi0_val{_mm_mask_blend_pd(is_sx_ge15,bessel_i0_le15_sse_pd(sx),bessel_i0_ge15_sse_pd(sx))};
  const __m128d exp_arg{_mm_mul_pd(half,xspow2)};
  const __m128d exp_val{_mm_exp_pd(negate_xmm2r8(exp_arg))};
  return (_mm_mul_pd(x,_mm_mul_pd(exp_val,besi0_val)));
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
integrand_4_16_marcum_Q_func_sse_pd(const __m128d beta,const __m128d psi,
                                    const __m128d theta)
{
  const __m128d one{_mm_set1_pd(1.0)};
  const __m128d half{_mm_set1_pd(0.5)};
  const __m128d psip2{_mm_mul_pd(psi,psi)};
  const __m128d sintht{_mm_sin_pd(theta)};
  const __m128d half_beta{_mm_mul_pd(half,_mm_mul_pd(beta,beta))};
  const __m128d lead_fac_num{_mm_add_pd(one,sintht)};
  const __m128d lead_fac_den{_mm_add_pd(_mm_add_pd(one,_mm_add_pd(sintht,sintht)),psip2)};
  const __m128d lead_fac_rat{_mm_div_pd(lead_fac_num,lead_fac_den)};
  const __m128d exp_arg{_mm_mul_pd(half_beta,lead_fac_den)};
  const __m128d exp_val{_mm_exp_pd(negate_xmm2r8(exp_arg))};
  return (_mm_mul_pd(lead_fac_rat,exp_val));
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
integrand_4_20_marcum_Q_func_lo_sse_pd(const __m128d beta,const __m128d psi,
                                       const __m128d theta)
{
  const __m128d one{_mm_set1_pd(1.0)};
  const __m128d half{_mm_set1_pd(0.5)};
  const __m128d psip2{_mm_mul_pd(psi,psi)};
  const __m128d costht{_mm_cos_pd(theta)};
  const __m128d lead_fac_num{_mm_sub_pd(one,costht)};
  const __m128d half_beta{_mm_mul_pd(half,_mm_mul_pd(beta,beta))};
  const __m128d lead_fac_den{_mm_add_pd(_mm_sub_pd(one,_mm_add_pd(costht,costht)),psip2)};
  const __m128d lead_fac_rat{_mm_div_pd(lead_fac_num,lead_fac_den)};
  const __m128d exp_arg{_mm_mul_pd(half_beta,lead_fac_den)};
  const __m128d exp_val{_mm_exp_pd(negate_xmm2r8(exp_arg))};
  return (_mm_mul_pd(lead_fac_rat,exp_val));
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
integrand_4_20_marcum_Q_func_hi_sse_pd(const __m128d beta,const __m128d psi,
                                       const __m128d theta)
{
  const __m128d one{_mm_set1_pd(1.0)};
  const __m128d half{_mm_set1_pd(0.5)};
  const __m128d psip2{_mm_mul_pd(psi,psi)};
  const __m128d costht{_mm_cos_pd(theta)};
  const __m128d lead_fac_num{_mm_add_pd(one,costht)};
  const __m128d half_beta{_mm_mul_pd(half,_mm_mul_pd(beta,beta))};
  const __m128d lead_fac_den{_mm_add_pd(_mm_add_pd(one,_mm_add_pd(costht,costht)),psip2)};
  const __m128d lead_fac_rat{_mm_div_pd(lead_fac_num,lead_fac_den)};
  const __m128d exp_arg{_mm_mul_pd(half_beta,lead_fac_den)};
  const __m128d exp_val{_mm_exp_pd(negate_xmm2r8(exp_arg))};
  return (_mm_mul_pd(lead_fac_rat,exp_val));
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
integrand_4_26_marcum_Q_func_sse_pd(const __m128d beta,const __m128d psi,
                                    const __m128d theta)
{
  const __m128d one{_mm_set1_pd(1.0)};
  const __m128d half{_mm_set1_pd(0.5)};
  const __m128d psip2{_mm_mul_pd(psi,psi)};
  const __m128d half_beta{_mm_mul_pd(half,_mm_mul_pd(beta,beta))};
  const __m128d psim_sintht{_mm_mul_pd(psi,_mm_sin_pd(theta))};
  const __m128d tmp{_mm_sub_pd(one,psip2)};
  const __m128d num{_mm_mul_pd(tmp,tmp)};
  const __m128d denom{_mm_add_pd(_mm_add_pd(one,_mm_add_pd(psim_sintht,psim_sintht)),psip2)};
  const __m128d exp_arg1{_mm_mul_pd(half_beta,denom)};
  const __m128d ratio{_mm_div_pd(num,denom)};
  const __m128d exp_arg2{_mm_mul_pd(half_beta,ratio)};
  const __m128d left_exp_val{_mm_exp_pd(negate_xmm2r8(exp_arg1))};
  const __m128d right_exp_val{_mm_exp_pd(negate_xmm2r8(exp_arg2))};
  return (_mm_add_pd(left_exp_val,right_exp_val));
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
integrand_4_27_marcum_Q_func_sse_pd(const __m128d alpha,const __m128d psi,
                                    const __m128d theta)
{
  const __m128d one{_mm_set1_pd(1.0)};
  const __m128d half{_mm_set1_pd(0.5)};
  const __m128d psip2{_mm_mul_pd(psi,psi)};
  const __m128d half_alpha{_mm_mul_pd(half,_mm_mul_pd(alpha,alpha))};
  const __m128d psim_sintht{_mm_mul_pd(psi,_mm_sin_pd(theta))};
  const __m128d tmp{_mm_sub_pd(one,psip2)};
  const __m128d num{_mm_mul_pd(tmp,tmp)};
  const __m128d denom{_mm_add_pd(_mm_add_pd(one,_mm_add_pd(psim_sintht,psim_sintht)),psip2)};
  const __m128d exp_arg1{_mm_mul_pd(half_alpha,denom)};
  const __m128d ratio{_mm_div_pd(num,denom)};
  const __m128d exp_arg2{_mm_mul_pd(half_alpha,ratio)};
  const __m128d left_exp_val{_mm_exp_pd(negate_xmm2r8(exp_arg1))};
  const __m128d right_exp_val{_mm_exp_pd(negate_xmm2r8(exp_arg2))};
  return (_mm_add_pd(left_exp_val,right_exp_val));
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
integrand_4_32_marcum_Q_m_func_sse_pd(const __m128d x,const __m128d s,
                                      const __m128d d_m,const std::int32_t i_m)
{
  /*
      This kernel relies on the scalar version of the modified cylindrical
      Bessel function of the first kind of the order n, hence substantial
      performance impact due to scalarization shall be expected. 
  */
  const __m128d half{_mm_set1_pd(0.5)};
  const __m128d sx{_mm_mul_pd(s,x)};
  const __m128d xx{_mm_mul_pd(x,x)};
  const double * __restrict__ p_sx{reinterpret_cast<const double * __restrict__>(&sx)};
  const __m128d ss{_mm_mul_pd(s,s)};
  const __m128d x_to_m{_mm_pow_pd(x,d_m)};
  const __m128d cyl_bes_val{_mm_setr_pd(std::cyl_bessel_i(i_m,p_sx[0]),
                                        std::cyl_bessel_i(i_m,p_sx[1]))};
  const __m128d exp_arg{_mm_mul_pd(half,_mm_add_pd(xx,ss))};
  const __m128d exp_val{_mm_exp_pd(negate_xmm2r8(exp_arg))};
  return (_mm_mul_pd(x_to_m,_mm_mul_pd(exp_val,cyl_bes_val)));
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
integrand_4_42_marcum_Q_m_func(const __m128d beta,const __m128d psi,
                               const __m128d m,const __m128d theta)
{
  const __m128d half{_mm_set1_pd(0.5)};
  const __m128d one{_mm_set1_pd(1.0)};
  const __m128d C15707963267948966192313216916398{_mm_set1_pd(1.5707963267948966192313216916398)};
  const __m128d psip2{_mm_mul_pd(psi,psi)};
  const __m128d halfbeta{_mm_mul_pd(half,_mm_mul_pd(beta,beta))};
  const __m128d psinpm{_mm_pow_pd(psi,_mm_sub_pd(m,one))};
  const __m128d tht_p_pi2{_mm_add_pd(theta,C15707963267948966192313216916398)};
  const __m128d cos_arg_left{_mm_sub_pd(m,_mm_mul_pd(one,tht_p_pi2))};
  const __m128d cos_val_left{_mm_cos_pd(cos_arg_left)};
  const __m128d cos_arg_right{_mm_mul_pd(m,tht_p_pi2)};
  const __m128d cos_val_right{_mm_cos_pd(cos_arg_right)};
  const __m128d psi_m_sintht{_mm_mul_pd(psi,_mm_sin_pd(theta))};
  const __m128d sin_factor{_mm_add_pd(one,_mm_add_pd(_mm_add_pd(psi_m_sintht,psi_m_sintht),psip2))};
  const __m128d exp_arg{_mm_mul_pd(halfbeta,sin_factor)};
  const __m128d cos_num_fac{_mm_mul_pd(psinpm,_mm_sub_pd(cos_val_left,_mm_mul_pd(psi,cos_val_right)))};
  const __m128d exp_val{_mm_exp_pd(negate_xmm2r8(exp_arg))};
  const __m128d ratio_factor{_mm_div_pd(cos_num_fac,sin_factor)};
  return (_mm_mul_pd(ratio_factor,exp_val));
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
integrand_4_45_marcum_Q_m_func_sse_pd(const __m128d theta,const __m128d beta,
                                      const __m128d m)
{
  const __m128d one{_mm_set1_pd(1.0)};
  const __m128d pow_arg{_mm_add_pd(one,_mm_add_pd(m,m))};
  const __m128d betap2{_mm_mul_pd(beta,beta)};
  const __m128d tmp_sin{_mm_sin_pd(theta)};
  const __m128d sinp2{_mm_mul_pd(tmp_sin,tmp_sin)};
  const __m128d costht{_mm_cos_pd(theta)};
  const __m128d denom{_mm_add_pd(sinp2,sinp2)};
  const __m128d sinpm{_mm_pow_pd(tmp_sin,pow_arg)};
  const __m128d exp_arg{_mm_div_pd(betap2,denom)};
  const __m128d cos_ratio{_mm_div_pd(costht,sinpm)};
  const __m128d exp_val{_mm_exp_pd(negate_xmm2r8(exp_arg))};
  return (_mm_mul_pd(cos_ratio,exp_val));
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
integrand_4_66_pawula_func_sse_pd(const __m128d psi,const __m128d delphi,
                                  const __m128d t,const __m128d A,
                                  const __m128d sigmasqr)
{
  const __m128d one{_mm_set1_pd(1.0)};
  const __m128d phase_diff{_mm_sub_pd(delphi,psi)};
  const __m128d A_factor{_mm_div_pd(_mm_mul_pd(A,A),_mm_add_pd(sigmasqr,sigmasqr))};
  const __m128d cost{_mm_cos_pd(t)};
  const __m128d cosphdif{_mm_cos_pd(phase_diff)};
  const __m128d cosfdiff_fac{_mm_sub_pd(one,_mm_mul_pd(cosphdif,cost))};
  const __m128d exp_arg{_mm_mul_pd(A_factor,cosfdiff_fac)};
  const __m128d inv_cosdiff{_mm_div_pd(one,cosfdiff_fac)};
  const __m128d exp_val{_mm_exp_pd(negate_xmm2r8(exp_arg))};
  return (_mm_mul_pd(inv_cosdiff,exp_val));
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
integrand_4_67_pawula_func_sse_pd(const __m128d psi,
                                  const __m128d t,const __m128d A,
                                  const __m128d sigmasqr)
{
  const __m128d one{_mm_set1_pd(1.0)};
  const __m128d A_factor{_mm_div_pd(_mm_mul_pd(A,A),_mm_add_pd(sigmasqr,sigmasqr))};
  const __m128d cost{_mm_cos_pd(t)};
  const __m128d cospsi{_mm_cos_pd(psi)};
  const __m128d cos_term{_mm_sub_pd(one,_mm_mul_pd(cospsi,cost))};
  const __m128d exp_arg{_mm_mul_pd(A_factor,cos_term)};
  const __m128d lead_factor{_mm_div_pd(one,cos_term)};
  const __m128d exp_val{_mm_exp_pd(negate_xmm2r8(exp_arg))};
  return (_mm_mul_pd(lead_factor,exp_val));
}


} // fading_channel

} // gms

#endif /*__GMS_INTEGRANDS_FUNC_CH4_SSE_H__*/