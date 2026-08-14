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

#ifndef __GMS_FUNCTIONALS_CH4_SSE_H__
#define __GMS_FUNCTIONALS_CH4_SSE_H__ 1108020260108

#include <cstdint>
#include "GMS_config.h"


namespace file_info 
{

     static const unsigned int GMS_FUNCTIONALS_CH4_SSE_MAJOR = 1;
     static const unsigned int GMS_FUNCTIONALS_CH4_SSE_MINOR = 1;
     static const unsigned int GMS_FUNCTIONALS_CH4_SSE_MICRO = 0;
     static const unsigned int GMS_FUNCTIONALS_CH4_SSE_FULLVER =
       1000U*GMS_FUNCTIONALS_CH4_SSE_MAJOR+100U*GMS_FUNCTIONALS_CH4_SSE_MINOR+
       10U*GMS_FUNCTIONALS_CH4_SSE_MICRO;
     static const char GMS_FUNCTIONALS_CH4_SSE_CREATION_DATE[] = "11-08-2026 01:07AM +00200 (TUE 11 AUG 2026 GMT+2)";
     static const char GMS_FUNCTIONALS_CH4_SSE_BUILD_DATE[]    = __DATE__; 
     static const char GMS_FUNCTIONALS_CH4_SSE_BUILD_TIME[]    = __TIME__;
     static const char GMS_FUNCTIONALS_CH4_SSE_SYNOPSIS[]      = "Numerical integration (computing of Functionals) of the integrands (ch.4) used for the fading channel PDF/CDF and BER,BEP,SER,SEP quadrature computation(QUADPACK).\
	                                                          Based on the M.K Simon, M.S. Alouini: Digital Communication over Fading Channels 1st ed\
															  ISBN-13 978-0471317791";

}

#if !defined(FUNCTIONALS_CH4_SSE_INTERMIX_LOAD_COMPUTE)
#define FUNCTIONALS_CH4_SSE_INTERMIX_LOAD_COMPUTE 1 
#endif 

// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (FUNCTIONALS_CH4_SSE_USE_PMC_INSTRUMENTATION)
#define FUNCTIONALS_CH4_SSE_USE_PMC_INSTRUMENTATION 0
#endif 

#if (FUNCTIONALS_CH4_SSE_USE_PMC_INSTRUMENTATION) == 1
#include "GMS_hw_perf_macros.h"

#define PMC_VARS                      \
uint64_t prog_counters_start[4] = {}; \
uint64_t prog_counters_stop[4]  = {}; \
uint64_t tsc_start,tsc_stop;          \
uint64_t act_cyc_start,act_cyc_stop;  \
uint64_t ref_cyc_start,ref_cyc_stop;  \
[[maybe_unused]] uint64_t dummy1;     \
[[maybe_unused]] uint64_t dummy2;     \
[[maybe_unused]] uint64_t dummy3;     \
int32_t core_counter_width;           \
double utilization,nom_ghz,avg_ghz;
#endif 

namespace gms 
{

namespace fading_channel
{


struct alignas(64) integrator_payload_sse_t
{
    double * __restrict__ __ATTR_ALIGN__(16) in_buf{nullptr};
    double * __restrict__ __ATTR_ALIGN__(16) out_buf{nullptr};
    double * __restrict__ __ATTR_ALIGN__(16) functional{nullptr};
    double * __restrict__ __ATTR_ALIGN__(16) functional_arg1{nullptr};
    double * __restrict__ __ATTR_ALIGN__(16) functional_arg2{nullptr};
    double                lo1;
    double                hi1;
    double                lo2;
    double                hi2;
    double                lo3;
    double                hi3;
    double                lo4;
    double                hi4;
    std::int32_t          n_func_args;
    std::int32_t          n_integrand_vals;
    std::int32_t          which_integrator; /*Currently unimplemented, will be done later*/
};  

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
std::int32_t 
integrate_4_1_gauss_Q_func_sse(integrator_payload_sse_t * __restrict__);


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
std::int32_t 
integrate_4_2_gauss_Q_func_sse(integrator_payload_sse_t * __restrict__);



} // fading_channel

} // gms

#endif /*__GMS_FUNCTIONALS_CH4_SSE_H__*/