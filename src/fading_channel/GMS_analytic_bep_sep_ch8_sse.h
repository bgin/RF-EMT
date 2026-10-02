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
#include "GMS_bessel_i0_sse.h"



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
tikhonov_phase_err_pdf_4xf32(const __m128 rho_eq,const __m128 phi_c)
{
const __m128 negPI         = _mm_set1_ps(-3.1415926535897932384626433832795f);
const float * __restrict__ p_phi_c = reinterpret_cast<const float * __restrict__>(&phi_c);
const __m128 posPI         = _mm_set1_ps(+3.1415926535897932384626433832795);
const __m128  bessi0       = _mm_setr_ps(std::cyl_bessel_i(0,p_phi_c[0],std::cyl_bessel_i(0,p_phi_c[1]),
                                         std::cyl_bessel_i(0,p_phi_c[2]),std::cyl_bessel_i(0,p_phi_c[3])));
const __mmask8 phic_le_pi   = _mm_cmp_ps_mask(phi_c,posPI,_CMP_LE_OQ);
const __m128  cosarg        = _mm_add_ps(phi_c,phi_c);
const __mmask8 phic_gt_pi   = _mm_cmp_ps_mask(phi_c,negPI,_CMP_GE_OQ);
const __mmask8 mask_pi      = _kand_mask8(phic_le_pi,phic_gt_pi);
const __m128  cos2phic      = _mm_mask_blend_ps(mask_pi,_mm_cos_ps(cosarg),gms::math::simd_fast_cos_approx_4xf32(cosarg));
const __m128  exp_val       = gms::math::simd_fast_exp_approx_4xf32(_mm_mul_ps(rho_eq,cos2phic));
return (_mm_div_ps(exp_val,_mm_mul_ps(posPI,bessi0)));
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
tikhonov_phase_err_pdf_2xf64(const __m128d rho_eq,__m128d phi_c)
{
const __m128d negPI        = _mm_set1_pd(-3.1415926535897932384626433832795);
const __m128d posPI        = _mm_set1_pd(+3.1415926535897932384626433832795);
const __m128d bessi0       = _mm_mask_blend_pd(_mm_cmp_pd_mask(rho_eq,_mm_set1_pd(15.0),_CMP_GE_OQ),
                                               gms::math::bessel_i0_le15_sse_pd(rho_eq),
                                               gms::math::bessel_i0_ge15_sse_pd(rho_eq));
const __mmask8 phic_le_pi   = _mm_cmp_pd_mask(phi_c,posPI,_CMP_LE_OQ);
const __m128d  cosarg       = _mm_add_pd(phi_c,phi_c);
const __mmask8 phic_gt_pi   = _mm_cmp_pd_mask(phi_c,negPI,_CMP_GE_OQ);
const __mmask8 mask_pi      = _kand_mask8(phic_le_pi,phic_gt_pi);
const __m128d  cos2phic     = _mm_mask_blend_pd(mask_pi,_mm_cos_pd(cosarg),gms::math::simd_fast_cos_approx_2xf64(cosarg));
const __m128d  exp_val      = gms::math::simd_fast_exp_approx_2xf64(_mm_mul_pd(rho_eq,cos2phic));
return (_mm_div_pd(exp_val,_mm_mul_pd(posPI,bessi0)));
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
G_PLL_bandwidth_4xf32(const __m128 Bl,const __m128 Tb)
{
return (_mm_rcp_ps(_mm_mul_ps(Bl,Tb)));
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
template<bool use_div_as_reciprocal>
__m128d 
G_PLL_bandwidth_2xf64(const __m128d Bl,const __m128d Tb)
{
if constexpr (use_div_as_reciprocal)
{
return (_mm_div_pd(_mm_set1_pd(1.0),_mm_mul_pd(Bl,Tb)));
}
else
{ 
return (_mm_rcp14_pd(_mm_mul_pd(Bl,Tb)));
}
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
G_PLL_bandwidth_2xf64_v2(const __m128d Bl,const __m128d Tb)
{
return (_mm_div_pd(_mm_set1_pd(1.0),_mm_mul_pd(Bl,Tb)));
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
BPSK_param_a_8_61_4xf32(const __m128 Ac,const __m128 Ts, 
                        const __m128 M, const __m128 N0,
                        const __m128 Bl,const __m128 Tb)
{
const __m128 Eb = avg_bit_E_to_carrier_A_4xf32(Ac,Ts,M);
const __m128 G  = G_PLL_bandwidth_4xf32(Bl,Tb);
const __m128 ratio = _mm_div_ps(Eb,_mm_add_ps(N0,N0));
const __m128 sqrtG = _mm_sub_ps(_mm_sqrt_ps(G),_mm_set1_ps(1.0));
return (_mm_sub_ps(ratio,_mm_mul_ps(sqrtG,sqrtG)));
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
BPSK_param_a_8_61_2xf64(const __m128d Ac,const __m128d Ts, 
                        const __m128d M, const __m128d N0,
                        const __m128d Bl,const __m128d Tb)
{
const __m128d Eb = avg_bit_E_to_carrier_A_2xf64(Ac,Ts,M);
const __m128d G  = G_PLL_bandwidth_2xf64_v2(Bl,Tb);
const __m128d ratio = _mm_div_pd(Eb,_mm_add_pd(N0,N0));
const __m128d sqrtG = _mm_sub_pd(_mm_sqrt_pd(G),_mm_set1_pd(1.0));
return (_mm_sub_pd(ratio,_mm_mul_pd(sqrtG,sqrtG)));
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
BPSK_param_b_8_61_4xf32(const __m128 Ac,const __m128 Ts, 
                        const __m128 M, const __m128 N0,
                        const __m128 Bl,const __m128 Tb)
{
const __m128 Eb = avg_bit_E_to_carrier_A_4xf32(Ac,Ts,M);
const __m128 G  = G_PLL_bandwidth_4xf32(Bl,Tb);
const __m128 ratio = _mm_div_ps(Eb,_mm_add_ps(N0,N0));
const __m128 sqrtG = _mm_add_ps(_mm_sqrt_ps(G),_mm_set1_ps(1.0));
return (_mm_sub_ps(ratio,_mm_mul_ps(sqrtG,sqrtG)));
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
BPSK_param_b_8_61_2xf64(const __m128d Ac,const __m128d Ts, 
                        const __m128d M, const __m128d N0,
                        const __m128d Bl,const __m128d Tb)
{
const __m128d Eb = avg_bit_E_to_carrier_A_2xf64(Ac,Ts,M);
const __m128d G  = G_PLL_bandwidth_2xf64_v2(Bl,Tb);
const __m128d ratio = _mm_div_pd(Eb,_mm_add_pd(N0,N0));
const __m128d sqrtG = _mm_add_pd(_mm_sqrt_pd(G),_mm_set1_pd(1.0));
return (_mm_sub_pd(ratio,_mm_mul_pd(sqrtG,sqrtG)));
}

}
}
#endif /*__GMS_ANALYTIC_BEP_SEP_CH8_SSE_H__*/