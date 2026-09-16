
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

#ifndef __GMS_INTEGRANDS_FUNC_CH8_H__
#define __GMS_INTEGRANDS_FUNC_CH8_H__ 160920261158

#include <cstdint>
#include "GMS_config.h"
#include "GMS_cephes_double.h"

namespace file_info 
{

     static const unsigned int GMS_INTEGRANDS_FUNC_CH8_MAJOR = 1;
     static const unsigned int GMS_INTEGRANDS_FUNC_CH8_MINOR = 1;
     static const unsigned int GMS_INTEGRANDS_FUNC_CH8_MICRO = 0;
     static const unsigned int GMS_INTEGRANDS_FUNC_CH8_FULLVER =
       1000U*GMS_INTEGRANDS_FUNC_CH8_MAJOR+100U*GMS_INTEGRANDS_FUNC_CH8_MINOR+
       10U*GMS_INTEGRANDS_FUNC_CH8_MICRO;
     static const char GMS_INTEGRANDS_FUNC_CH8_CREATION_DATE[] = "16-09-2026 11:58AM +00200 (WED 16 AUG 2026 GMT+2)";
     static const char GMS_INTEGRANDS_FUNC_CH8_BUILD_DATE[]    = __DATE__; 
     static const char GMS_INTEGRANDS_FUNC_CH8_BUILD_TIME[]    = __TIME__;
     static const char GMS_INTEGRANDS_FUNC_CH8_SYNOPSIS[]      = "Integrands of the functional form used for fading channel and BER,BEP,SER,SEP quadrature computation(QUADPACK).\
	                                                          Based on the chapter 8 of M.K Simon, M.S. Alouini: Digital Communication over Fading Channels 1st ed\
															  ISBN-13 978-0471317791";

}

#if !defined(INTEGRANDS_FUNC_CH8_USE_CEPHES_DOUBLE)
#define INTEGRANDS_FUNC_CH8_USE_CEPHES_DOUBLE 1
#endif 

#if !defined(INTEGRANDS_FUNC_CH8_DO_WARMUP_CALL)
#define INTEGRANDS_FUNC_CH8_DO_WARMUP_CALL 1
#endif 

#if !defined(INTEGRANDS_FUNC_CH8_SPLIT_MATH_IMPL_PERF_TEST)
#define INTEGRANDS_FUNC_CH8_SPLIT_MATH_IMPL_PERF_TEST 1
#endif 

#if !defined(INTEGRANDS_FUNC_CH8_USE_PARAM_MACRO)
#define INTEGRANDS_FUNC_CH8_USE_PARAM_MACRO 0
#endif

#if !defined(INTEGRANDS_FUNC_CH8_OVERRIDE_COMPILER_CMD_LINE)
#define INTEGRANDS_FUNC_CH8_OVERRIDE_COMPILER_CMD_LINE 0
#endif 


namespace gms
{

namespace fading_channel
{

struct alignas(64) func_args_ch8_payload_t
{
      double arg1d;
      double arg2d;
      double arg3d;
      double arg4d;
      double arg5d;
};

#if (INTEGRANDS_FUNC_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
__ATTR_HOT__
__ATTR_ALIGN__(32)
double analytic_SEP_MAM_8_1(const double,const double,
                            const double,const double);

/*
   Symbol Error Probability (SEP) of multiple amplitude shift keying
*/
#if (INTEGRANDS_FUNC_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_SEP_MAM_8_5(const double,const double,
                             const double,const double,
                             const double);

#if (INTEGRANDS_FUNC_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_SEP_MAM_8_5_iface(const double,void * __restrict__);

#if (INTEGRANDS_FUNC_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
__ATTR_HOT__
__ATTR_ALIGN__(32)
double analytic_SEP_QAM_8_10(const double,const double,
                             const double,const double);

#if (INTEGRANDS_FUNC_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
__ATTR_HOT__
__ATTR_ALIGN__(32)
double analytic_SEP_QAM4_8_11(const double,const double,
                              const double,const double);

#if (INTEGRANDS_FUNC_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_SEP_QAM_8_10(const double,const double,
                              const double,const double,
                              const double);

#if (INTEGRANDS_FUNC_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_SEP_QAM_8_10_iface(const double,void * __restrict__);

#if (INTEGRANDS_FUNC_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_SEP_QAM4_8_11(const double,const double,
                              const double,const double,
                              const double);

#if (INTEGRANDS_FUNC_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_SEP_QAM4_8_11_iface(const double,void * __restrict__);


#if (INTEGRANDS_FUNC_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_SEP_MPSK_8_17(const double,const double,
                               const double,const double,
                               const double);

#if (INTEGRANDS_FUNC_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
__ATTR_HOT__
__ATTR_ALIGN__(32)
double integrand_SEP_MPSK_8_17_iface(const double,void * __restrict__);



}

}

#endif /*__GMS_INTEGRANDS_FUNC_CH8_H__*/