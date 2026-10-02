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

#ifndef __GMS_ANALYTIC_BEP_SEP_CH8_H__
#define __GMS_ANALYTIC_BEP_SEP_CH8_H__ 011020260105

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_gaussian_Q_approx_simd.h"
#include "GMS_marcum_Q_approx_simd.h"



namespace file_info 
{

     static const unsigned int GMS_ANALYTIC_BEP_SEP_CH8_SSE_MAJOR = 1;
     static const unsigned int GMS_ANALYTIC_BEP_SEP_CH8_SSE_MINOR = 1;
     static const unsigned int GMS_ANALYTIC_BEP_SEP_CH8_SSE_MICRO = 0;
     static const unsigned int GMS_ANALYTIC_BEP_SEP_CH8_SSE_FULLVER =
       1000U*GMS_ANALYTIC_BEP_SEP_CH8_SSE_MAJOR+100U*GMS_ANALYTIC_BEP_SEP_CH8_SSE_MINOR+
       10U*GMS_ANALYTIC_BEP_SEP_CH8_SSE_MICRO;
     static const char GMS_ANALYTIC_BEP_SEP_CH8_SSE_CREATION_DATE[] = "01-10-2026 01:05PM +00200 (THR 01 OCT 2026 GMT+2)";
     static const char GMS_ANALYTIC_BEP_SEP_CH8_SSE_BUILD_DATE[]    = __DATE__; 
     static const char GMS_ANALYTIC_BEP_SEP_CH8_SSE_BUILD_TIME[]    = __TIME__;
     static const char GMS_ANALYTIC_BEP_SEP_CH8_SSE_AUTHOR[]        = "Programmer: Bernard Gingold, beniekg@gmail.com";
     static const char GMS_ANALYTIC_BEP_SEP_CH8_SSE_SYNOPSIS[]      = "Analytic Formulae of the Bit and Symbol Error Probabilities based on the approximated Gaussian-Q function.\
	                                                          Based on the chapter 8 of M.K Simon, M.S. Alouini: Digital Communication over Fading Channels 1st ed\
															  ISBN-13 978-0471317791";

}

#if !defined(ANALYTIC_BEP_SEP_CH8_SSE_OVERRIDE_COMPILER_CMD_LINE)
#define ANALYTIC_BEP_SEP_CH8_SSE_OVERRIDE_COMPILER_CMD_LINE 0
#endif

namespace gms 
{

namespace fading_channel
{

enum class Gaussian_Q_approxmations_sse_t : std::int32_t 
{
    Gaussian_Q_approx_chiani,
    Gaussian_Q_approx_loskot_2T,
    Gaussian_Q_approx_loskot_3T,
    Gaussian_Q_approx_sadhwani_1T,
    Gaussian_Q_approx_sadhwani_2T,
    Gaussian_Q_approx_sadhwani_4T,
    Gaussian_Q_approx_borjesson,
    Gaussian_Q_approx_sadhwani_summed,
    Gaussian_Q_approx_karagiannidis_lioumpas
};

#if (ANALYTIC_BEP_SEP_CH8_SSE_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=sse
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
__m128 
avg_symbol_E_to_carrier_A_4xf32(const __m128 Ac,const __m128 Tc,
                                const __m128 M)
{
const __m128 third    = _mm_set1_ps(0.3333333333333333333333333333333333333333f);
const __m128 AA       = _mm_mul_ps(Ac,Ac);
const __m128 MM       = _mm_mul_ps(M,M); 
const __m128 one      = _mm_mul_ps(1.0f);
const __m128 MM_ratio = _mm_sub_ps(_mm_mul_ps(MM,one),third);
const __m128 lead_fac = _mm_mul_ps(AA,Ts);
return (_mm_mul_ps(lead_fac,MM_ratio)); 
}

#if (ANALYTIC_BEP_SEP_CH8_SSE_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=sse
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
__m128d 
avg_symbol_E_to_carrier_A_2xf64(const __m128d Ac,const __m128d Tc,
                                const __m128d M)
{
const __m128d third    = _mm_set1_pd(0.3333333333333333333333333333333333333333);
const __m128d AA       = _mm_mul_pd(Ac,Ac);
const __m128d MM       = _mm_mul_pd(M,M); 
const __m128d one      = _mm_mul_pd(1.0);
const __m128d MM_ratio = _mm_sub_pd(_mm_mul_pd(MM,one),third);
const __m128d lead_fac = _mm_mul_pd(AA,Ts);
return (_mm_mul_pd(lead_fac,MM_ratio)); 
}

#if (ANALYTIC_BEP_SEP_CH8_SSE_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=sse
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif
__m128 
avg_bit_E_to_carrier_A_4xf32(const __m128 Ac,const __m128 Ts, 
                             const __m128 M)
{
const __m128 num   = _mm_fmsub_ps(M,M,_mm_set1_ps(1.0));
const __m128 Eg    = _mm_mul_ps(_mm_mul_ps(Ac,Ac,),Ts);
const __m128 log2M = _mm_mul_ps(_mm_set1_ps(6.0),_mm_log2_ps(M));
return (_mm_mul_ps(Eg,_mm_div_ps(num,log2M))); 
}

#if (ANALYTIC_BEP_SEP_CH8_SSE_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=sse
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif
__m128d 
avg_bit_E_to_carrier_A_2xf64(const __m128d Ac,const __m128d Ts, 
                             const __m128d M)
{
const __m128d num   = _mm_fmsub_pd(M,M,_mm_set1_pd(1.0));
const __m128d Eg    = _mm_mul_pd(_mm_mul_pd(Ac,Ac,),Ts);
const __m128d log2M = _mm_mul_pd(_mm_set1_pd(6.0),_mm_log2_pd(M));
return (_mm_mul_pd(Eg,_mm_div_pd(num,log2M))); 
}


}

}

#endif /*__GMS_ANALYTIC_BEP_SEP_CH8_SSE_H__*/