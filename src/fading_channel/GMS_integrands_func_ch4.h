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

#ifndef __GMS_INTEGRANDS_FUNC_CH4_H__
#define __GMS_INTEGRANDS_FUNC_CH4_H__ 280720260622

#include <cstdint>
#include "GMS_config.h"


namespace file_info 
{

     static const unsigned int GMS_INTEGRANDS_FUNC_CH4_MAJOR = 1;
     static const unsigned int GMS_INTEGRANDS_FUNC_CH4_MINOR = 1;
     static const unsigned int GMS_INTEGRANDS_FUNC_CH4_MICRO = 0;
     static const unsigned int GMS_INTEGRANDS_FUNC_CH4_FULLVER =
       1000U*GMS_INTEGRANDS_FUNC_CH4_MAJOR+100U*GMS_INTEGRANDS_FUNC_CH4_MINOR+
       10U*GMS_INTEGRANDS_FUNC_CH4_MICRO;
     static const char GMS_INTEGRANDS_FUNC_CH4_CREATION_DATE[] = "28-07-2026 11:20AM +00200 (MON 27 JUL 2026 GMT+2)";
     static const char GMS_INTEGRANDS_FUNC_CH4_BUILD_DATE[]    = __DATE__; 
     static const char GMS_INTEGRANDS_FUNC_CH4_BUILD_TIME[]    = __TIME__;
     static const char GMS_INTEGRANDS_FUNC_CH4_SYNOPSIS[]      = "Integrands of the functional form used for fading channel and BER,BEP,SER,SEP quadrature computation(QUADPACK).\
	                                                          Based on the M.K Simon, M.S. Alouini: Digital Communication over Fading Channels 1st ed\
															  ISBN-13 978-0471317791";

}

namespace gms
{

namespace fading_channel
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




/*
  Chapter 4, formula: 4.1
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_std_lib,bool do_funcs_warmup>
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_4_1_gauss_Q_func(const double);

/*
  Chapter 4, formula: 4.2
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_std_lib,bool do_funcs_warmup>
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_4_2_gauss_Q_func(const double,const double);

/*
  Chapter 4, formula: 4.6
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_std_lib,bool do_funcs_warmup,
         std::int32_t choose_sin_or_cos>
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_4_6_gauss_Q_func(const double,const double,
                                  const double,const double);
                                  
/*
  Chapter 4, formula: 4.7
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_std_lib,bool do_funcs_warmup,
         std::int32_t choose_x1_or_y1>
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_4_7_gauss_Q_func(const double,const double,
                             const double,const double);

/*
  Chapter 4, formula: 4.8 (case: rho = 0)
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_std_lib,bool do_funcs_warmup,
         std::int32_t choose_x1_or_y1>
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_4_8_gauss_Q_func(const double,const double,const double);

/*
   Chapter 4, formula: 4.10
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_std_lib,bool do_funcs_warmup>
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_4_10_marcum_Q_func(const double,const double);

/*
   Chapter 4, formula: 4.16
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_std_lib,bool do_funcs_warmup>
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_4_16_marcum_Q_func(const double,
                                    const double,const double);

/*
   Chapter 4, formula: 4.20
   Range: -PI,0
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_std_lib,bool do_funcs_warmup>
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_4_20_marcum_Q_func_lo(const double,const double,
                                       const double);

/*
   Chapter 4, formula: 4.20
   Range: 0,PI
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_std_lib,bool do_funcs_warmup>
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_4_20_marcum_Q_func_hi(const double,const double,
                                       const double);

/*
   Chapter 4, formula: 4.26
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_std_lib,bool do_funcs_warmup>
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_4_26_marcum_Q_func(const double,const double,
                               const double,const double);

/*
   Chapter 4, formula: 4.27
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_std_lib,bool do_funcs_warmup>
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_4_27_marcum_Q_func(const double,const double,
                               const double,const double);

/*
   Chapter 4, formula: 4.32
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_std_lib,bool do_funcs_warmup>
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_4_32_marcum_Q_m_func(const double,const double,const std::int32_t);

/*
   Chapter 4, formula: 4.42
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_std_lib,bool do_funcs_warmup>
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_4_42_marcum_Q_m_func(const double,const double,
                                      const double,std::int32_t);

/*
   Chapter 4, formula: 4.45
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_std_lib,bool do_funcs_warmup>
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_4_45_marcum_Q_m_func(const double,const double,std::int32_t);

/*
   Chapter 4, formula: 4.50
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_std_lib,bool do_funcs_warmup>
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_4_50_marcum_Q_m_func(const double,const double,
                                      const double,const double);

/*
   Chapter 4, formula: 4.66
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_std_lib,bool do_funcs_warmup>
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_4_66_pawula_func(const double,const double,const double);

/*
   Chapter 4, formula: 4.67
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_std_lib,bool do_funcs_warmup>
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_4_67_pawula_func(const double,const double);

/*
   Chapter 4, formula: 4.68
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_std_lib,bool do_funcs_warmup>
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_4_68_pawula_func(const double,const double);


                            
} // fading_channel

} // gms

#endif /*__GMS_INTEGRANDS_FUNC_CH4_H__*/