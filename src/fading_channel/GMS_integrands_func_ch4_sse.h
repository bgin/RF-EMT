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
#include "GMS_config.h"
#include "GMS_simd_utils.h"

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
#define INTEGRANDS_func_CH4_SSE_DO_WARMUP_CALL 1
#endif

namespace gms 
{

namespace fading_channel
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
__m128d
integrand_4_1_gauss_Q_func_sse_pd(const __m128d y)
{
  const __m128 vC03989422804014326779399460599344{_mm_set1_pd(0.3989422804014326779399460599344)};
  __m128d exp_arg{_mm_mul_pd(_mm_set1_pd(0.5),_mm_mul_pd(y,y))};
  __m128d exp_val{_mm_exp_pd(gms::math::negate_xmm2r8(exp_arg))};
  return (_mm_mul_ps(vC03989422804014326779399460599344,exp_val));
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
integrand_4_2_gauss_Q_func_sse_pd(const __m128d x,
                                  const __m128d theta)
{
  const __m128d xsqr{_mm_mul_pd(x,x)};
  const __m128d sin_val{_mm_sin_pd(theta)};
  const __m128d tmp_sin_mul{_mm_mul_pd(sin_val,sin_val)};
  const __m128d sin_denom{_mm_add_ps(tmp_sin_mul,tmp_sin_mul)};
  const __m128d exp_arg{_mm_div_pd(xsqr,sin_denom)};
  const __m128d exp_val{_mm_exp_pd(gms::math::negate_xmm2r8(exp_arg))};
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
integrand_4_6_sin_gauss_Q_func_sse_pd(const __m128d x1,const __m128d y1,
                                      const __m128d rho,const __m128d theta,
                                      const phi_s /*an upper limit of integration, computed one per integration full step*/)
{
  const __m128d one{_mm_set1_pd(1.0)};
  const __m128d squared_sum{_mm_fmadd_pd(x1,x1,_mm_mul_pd(y1,y1))};
  const __m128d tsin{_mm_sin_pd(theta)};
  const __m128d sin_sqr_tht{_mm_mul_pd(tsin,tsin)};
  const __m128d S_hat{_mm_sqrt_pd(squared_sum)};
  const __m128d one_m_rho{_mm_sub_pd(one,_mm_mul_pd(rho,rho))};
  const __m128d sin2theta{_mm_sin_pd(_mm_add_pd(theta,theta))};
  const __m128d one_m_sin2theta{_mm_sub_pd(one,_mm_mul_pd(rho,sin2theta))};
  const __m128d sqr_1m_rho{_mm_sqrt_pd(on_m_rho)};
  const __m128d S_hat_half{_mm_mul_pd(_mm_set1_pd(0.5),_mm_mul_pd(S_hat,S_hat))};
  const __m128d ratio1{_mm_div_pd(sqr_1m_rho,one_m_sin2theta)};
  const __m128d sin_phi_s{_mm_sin_pd(phi_s)};
  const __m128d sin_phi_s_sqr{_mm_mul_pd(sin_phi_s,sin_phi_s)};
  const __m128d ratio_sin{_mm_div_pd(sin_phi_s_sqr,sin_sqr_tht)};
  const __m128d ratio2{_mm_div_pd(one_m_sin2theta,one_m_rho)};
  const __m128d exp_arg{_mm_mul_pd(S_hat_half,_mm_mul_pd(ratio2,ratio_sin))};
  const __m128d exp_val{_mm_exp_pd(gms::math::negate_xmm2r8(exp_arg))};
  return (_mm_mul_pd(ratio1,exp_val));
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
integrand_4_6_cos_gauss_Q_func_sse_pd(const __m128d x1,const __m128d y1,
                                      const __m128d rho,const __m128d theta,
                                      const phi_s /*an upper limit of integration, computed one per integration full step*/)
{
  const __m128d one{_mm_set1_pd(1.0)};
  const __m128d squared_sum{_mm_fmadd_pd(x1,x1,_mm_mul_pd(y1,y1))};
  const __m128d tsin{_mm_sin_pd(theta)};
  const __m128d sin_sqr_tht{_mm_mul_pd(tsin,tsin)};
  const __m128d S_hat{_mm_sqrt_pd(squared_sum)};
  const __m128d one_m_rho{_mm_sub_pd(one,_mm_mul_pd(rho,rho))};
  const __m128d sin2theta{_mm_sin_pd(_mm_add_pd(theta,theta))};
  const __m128d one_m_sin2theta{_mm_sub_pd(one,_mm_mul_pd(rho,sin2theta))};
  const __m128d sqr_1m_rho{_mm_sqrt_pd(on_m_rho)};
  const __m128d S_hat_half{_mm_mul_pd(_mm_set1_pd(0.5),_mm_mul_pd(S_hat,S_hat))};
  const __m128d ratio1{_mm_div_pd(sqr_1m_rho,one_m_sin2theta)};
  const __m128d cos_phi_s{_mm_cos_pd(phi_s)};
  const __m128d cos_phi_s_sqr{_mm_mul_pd(cos_phi_s,cos_phi_s)};
  const __m128d ratio_cos{_mm_div_pd(cos_phi_s_sqr,sin_sqr_tht)};
  const __m128d ratio2{_mm_div_pd(one_m_sin2theta,one_m_rho)};
  const __m128d exp_arg{_mm_mul_pd(S_hat_half,_mm_mul_pd(ratio2,ratio_cos))};
  const __m128d exp_val{_mm_exp_pd(gms::math::negate_xmm2r8(exp_arg))};
  return (_mm_mul_pd(ratio1,exp_val));
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
  const __m128d exp_val{_mm_exp_pd(gms::math::negate_xmm2r8(exp_arg))};
  return (_mm_mul_pd(lead_factor,exp_val));
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
  const __m128d exp_val{_mm_exp_pd(gms::math::negate_xmm2r8(exp_arg))};
  return (_mm_mul_pd(lead_factor,exp_val));
}


}

}


#endif /*__GMS_INTEGRANDS_FUNC_CH4_SSE_H__*/