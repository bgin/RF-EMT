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

#ifndef __GMS_INTEGRANDS_FUNC_CH5_H__
#define __GMS_INTEGRANDS_FUNC_CH5_H__ 250820260640

#include <cstdint>
#include "GMS_config.h"

namespace file_info 
{

     static const unsigned int GMS_INTEGRANDS_FUNC_CH5_MAJOR = 1;
     static const unsigned int GMS_INTEGRANDS_FUNC_CH5_MINOR = 1;
     static const unsigned int GMS_INTEGRANDS_FUNC_CH5_MICRO = 0;
     static const unsigned int GMS_INTEGRANDS_FUNC_CH5_FULLVER =
       1000U*GMS_INTEGRANDS_FUNC_CH5_MAJOR+100U*GMS_INTEGRANDS_FUNC_CH5_MINOR+
       10U*GMS_INTEGRANDS_FUNC_CH5_MICRO;
     static const char GMS_INTEGRANDS_FUNC_CH5_CREATION_DATE[] = "25-08-2026 07:09AM +00200 (TUE 25 AUG 2026 GMT+2)";
     static const char GMS_INTEGRANDS_FUNC_CH5_BUILD_DATE[]    = __DATE__; 
     static const char GMS_INTEGRANDS_FUNC_CH5_BUILD_TIME[]    = __TIME__;
     static const char GMS_INTEGRANDS_FUNC_CH5_SYNOPSIS[]      = "Integrands of the functional form used for fading channel and BER,BEP,SER,SEP quadrature computation(QUADPACK).\
	                                                          Based on the M.K Simon, M.S. Alouini: Digital Communication over Fading Channels 1st ed\
															  ISBN-13 978-0471317791";

}

#if !defined(INTEGRANDS_FUNC_CH5_USE_CEPHES_DOUBLE)
#define INTEGRANDS_FUNC_CH5_USE_CEPHES_DOUBLE 1
#endif 

#if !defined(INTEGRANDS_FUNC_CH5_DO_WARMUP_CALL)
#define INTEGRANDS_FUNC_CH5_DO_WARMUP_CALL 1
#endif 

#if !defined(INTEGRANDS_FUNC_CH5_SPLIT_MATH_IMPL_PERF_TEST)
#define INTEGRANDS_FUNC_CH5_SPLIT_MATH_IMPL_PERF_TEST 1
#endif 

#if !defined(INTEGRANDS_FUNC_CH5_USE_PARAM_MACRO)
#define INTEGRANDS_FUNC_CH5_USE_PARAM_MACRO 0
#endif


namespace gms
{

namespace fading_channel
{


struct alignas(64) func_args_ch5_payload_t
{
      double arg1d;
      double arg2d;
      double arg3d;
      double arg4d;
      double arg5d;
      double arg6d;
      double arg7d;
      std::int32_t arg8i;
};

namespace detail
{

/*
   stdlibc++ warmup calls
*/
__ATTR_HOT__
__ATTR_ALIGN__(32)
double warmup_sin(const double);

__ATTR_HOT__
__ATTR_ALIGN__(32)
double warmup_cos(const double);

__ATTR_HOT__
__ATTR_ALIGN__(32)
double warmup_tan(const double);

__ATTR_HOT__
__ATTR_ALIGN__(32)
double warmup_exp(const double);

__ATTR_HOT__
__ATTR_ALIGN__(32)
double warmup_cyl_bess_j(const std::int32_t,const double);

}// detail

/*
   Chapter 5, formula: 5.6
   Rayleigh Fading Channel SNR PDF passed through the MGF (Moment Generating Function)
   i.e. the Laplace Transform and further integrated.
   I(a,gamma), where a is constant modulation-detection dependent, gamma is instantenous SNR.
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_Rayleigh_chan_5_6(const double,const double,
                                   const double);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_Rayleigh_chan_5_6_iface(const double,void * __restrict__);

/*
    Chapter 5, formula: 5.9
    Nakagami-q (Hoyt) Fading Channel
    I(a,q,gamma), where a is constant modulation-detection dependent, gamma is instantenous SNR, q - PDF parameter.
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_Hoyt_chan_5_9(const double,const double,
                               const double,const double);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_Hoyt_chan_5_9_iface(const double,void * __restrict__);

/*
   Chapter 5, formula: 5.12
   Nakagami-n (Rice) Fading Channel
   I(a,n,gamma), where a is constant modulation-detection dependent, gamma is instantenous SNR, n -  PDF parameter.
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_Rice_chan_5_12(const double,const double,
                                const double,const double);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_Rice_chan_5_12_iface(const double,void * __restrict__);

/*
   Chapter 5, formula: 5.16
   Nakagami-m Fading Channel
   I(a,m,gamma), where a is constant modulation-detection dependent, gamma is instantenous SNR, n -  PDF parameter.
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_Nakagami_m_chan_5_16(const double,const double,
                                      const double,const double);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_Nakagami_m_chan_5_16_iface(const double,void * __restrict__);

/*
    Chapter 5, formula: 5.20,5.21
    Log-Normal Shadowing channel
    I(a,mu,sigma), where a is constant modulation-detection dependent,mu is a 10*log10(gamma),sigma - logarithmic standard deviation of shadowing
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_Log_Norm_Shadow_chan_5_20(const double,const double,
                                           const double, const double,
                                           const double);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_Log_Norm_Shadow_chan_5_20_iface(const double,void * __restrict__);

/*
    Chapter 5, formula: 5.25
    Composite Log-Normal Shadowing and Nakagami-m channel
    I(a,mu,sigma,m), where a is constant modulation-detection dependent,mu is a 10*log10(gamma),sigma - logarithmic standard deviation of shadowing
                     m -- parameter of the Nakagami-m PDF
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_LNSh_Nakagami_m_chan_5_25(const double,const double,
                                          const double,const double,
                                          const double);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_LNSh_Nakagami_m_chan_5_25_iface(const double,void * __restrict__);

/*
    Chapter 5, formula: 5.39
    Rayleigh Fading Channel (Laplace Transform of SNR (per-bit) PDF)
    I(b,psi,gamma,l), described by the 5.38
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_Rayleigh_lt_chan_5_39(const double,const double,
                                       const double,const double,
                                       const double);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_Rayleigh_lt_chan_5_39_iface(const double,void * __restrict__);

/*
    Chapter 5, formula: 5.40
    Nakagami-q (Hoyt) Fading Channel (Laplace Transform of SNR (per-bit) PDF)
    I(b,psi,n,gamma,l), described by the 5.38
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_Hoyt_lt_chan_5_40(const double,const double,
                                   const double,const double,
                                   const double,const double);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_Hoyt_lt_chan_5_40_iface(const double,void * __restrict__);

/*
    Chapter 5, formula: 5.41
    Nakagami-n (Rice) Fading Channel (Laplace Transform of SNR (per-bit) PDF)
    I(b,psi,n,gamma,l), described by the 5.38
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_Rice_lt_chan_5_41(const double,const double,
                                   const double,const double,
                                   const double,const double);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_Rice_lt_chan_5_41_iface(const double,void * __restrict__);

/*
    Chapter 5, formula: 5.43
    Nakagami-m Fading Channel (Laplace Transform of SNR (per-bit) PDF)
    I(b,psi,n,gamma,l), described by the 5.38
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_Nakagami_m_lt_chan_5_43(const double,const double,
                                         const double,const double,
                                         const double,const double);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_Nakagami_m_lt_chan_5_43_iface(const double,void * __restrict__);

/*
    Chapter 5, formula: 5.44
    Log-Normal Shadowing Fading Channel (Laplace Transform of SNR (per-bit) PDF)
    I(b,psi,mu,sigma,l), described by the 5.38
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_Log_Norm_Shadow_lt_chan_5_44(const double,const double,
                                              const double,const double,
                                              const double,const double);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_Log_Norm_Shadow_lt_chan_5_44_iface(const double,void * __restrict__);

/*
    Chapter 5, formula: 5.47
    Log-Normal Shadowing composed with the Nakagami-m Fading Channel (Laplace Transform of SNR (per-bit) PDF)
    I(b,psi,mu,sigma,m,l), described by the 5.38
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_lnsh_Nakagami_m_lt_cha_5_47(const double,const double,
                                             const double,const double,
                                             const double,const double,
                                             const double);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_lnsh_Nakagami_m_lt_chan_5_47_iface(const double,void * __restrict__);

/*
   ai -- SNR for the ith signal
   thetai -- angle associated with the correct signal detection decision region (the ith signal)
   phii   -- same as above (the ith signal)
   The integral will be computed for the each ith signal in the 2D signal set.
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_avg_err_prob_QAM_5_71(const double,const double,const double);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_avg_err_prob_QAM_5_71_iface(const double,void * __restrict__);

/*
   ai - SNR for the ith signal
   thetai -- angle associated with the correct signal detection decision region (the ith signal)
   psii   -- same as above (the ith signal)
   gamma -- SNR
   m - Nakagami-m PDF parameter
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_avg_err_prob_Nakagami_m_QAM_5_76(const double,const double,
                                                  const double,const double,
                                                  const double,const double);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_avg_err_prob_Nakagami_m_QAM_5_76_iface(const double,void * __restrict__);


} // fading_channel

} // gms

#endif /*__GMS_INTEGRANDS_FUNC_CH5_H__*/