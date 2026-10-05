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

enum class Gaussian_Q_approximations_sse_t : std::int32_t 
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
OQPSK_param_a1_8_61_4xf32(const __m128 Ac,const __m128 Ts, 
                        const __m128 M, const __m128 N0,
                        const __m128 Bl,const __m128 Tb)
{
return (BPSK_param_a_8_61_4xf32(Ac,Ts,M,N0,Bl,Tb));
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
OQPSK_param_b1_8_61_2xf64(const __m128d Ac,const __m128d Ts, 
                          const __m128d M, const __m128d N0,
                          const __m128d Bl,const __m128d Tb)
{
return (BPSK_param_a_8_61_2xf64(Ac,Ts,M,N0,Bl,Tb));
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
OQPSK_param_a2_8_64_4xf32(const __m128 Ac,const __m128 Ts, 
                        const __m128 M, const __m128 N0,
                        const __m128 Bl,const __m128 Tb)
{
const __m128 Eb = avg_bit_E_to_carrier_A_4xf32(Ac,Ts,M);
const __m128 G  = G_PLL_bandwidth_4xf32(Bl,Tb);
const __m128 sqrt2G = _mm_sub_ps(_mm_set1_ps(1.0),_mm_sqrt_ps(_mm_add_ps(G,H)));
const __m128 factor = _mm_add_ps(G,sqrt2G);
return (_mm_mul_ps(_mm_div_ps(Eb,N0),factor));
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
OQPSK_param_a2_8_64_2xf64(const __m128d Ac,const __m128d Ts, 
                          const __m128d M, const __m128d N0,
                          const __m128d Bl,const __m128d Tb)
{
const __m128d Eb = avg_bit_E_to_carrier_A_2xf64(Ac,Ts,M);
const __m128d G  = G_PLL_bandwidth_2xf64_v2(Bl,Tb);
const __m128d sqrt2G = _mm_sub_pd(_mm_set1_pd(1.0),_mm_sqrt_pd(_mm_add_pd(G,H)));
const __m128d factor = _mm_add_pd(G,sqrt2G);
return (_mm_mul_pd(_mm_div_pd(Eb,N0),factor));
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
OQPSK_param_b2_8_64_4xf32(const __m128 Ac,const __m128 Ts, 
                          const __m128 M, const __m128 N0,
                          const __m128 Bl,const __m128 Tb)
{
return (OQPSK_param_a2_8_64_4xf32(Ac,Ts,M,N0,Bl,Tb));
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
OQPSK_param_b2_8_64_2xf64(const __m128d Ac,const __m128d Ts, 
                          const __m128d M, const __m128d N0,
                          const __m128d Bl,const __m128d Tb)
{
return (OQPSK_param_a2_8_64_2xf64(Ac,Ts,M,N0,Bl,Tb));
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
MSK_param_a1_8_65_4xf32(const __m128 Ac,const __m128 Ts, 
                        const __m128 M, const __m128 N0,
                        const __m128 Bl,const __m128 Tb)
{
return (BPSK_param_a_8_61_4xf32(Ac,Ts,M,N0,Bl,Tb));
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
MSK_param_a1_8_65_2xf64(const __m128d Ac,const __m128d Ts, 
                        const __m128d M, const __m128d N0,
                        const __m128d Bl,const __m128d Tb)
{
return (BPSK_param_a_8_61_2xf64(Ac,Ts,M,N0,Bl,Tb));
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
MSK_param_b1_8_65_4xf32(const __m128 Ac,const __m128 Ts, 
                        const __m128 M, const __m128 N0,
                        const __m128 Bl,const __m128 Tb)
{
return (BPSK_param_b_8_61_4xf32(Ac,Ts,M,N0,Bl,Tb));
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
MSK_param_b1_8_65_2xf64(const __m128d Ac,const __m128d Ts, 
                        const __m128d M, const __m128d N0,
                        const __m128d Bl,const __m128d Tb)
{
return (BPSK_param_b_8_61_2xf64(Ac,Ts,M,N0,Bl,Tb));
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
MSK_param_a2_8_65_4xf32(const __m128 Ac,const __m128 Ts, 
                        const __m128 M, const __m128 N0,
                        const __m128 Bl,const __m128 Tb)
{
const __m128 Eb = avg_bit_E_to_carrier_A_4xf32(Ac,Ts,M);
const __m128 G  = G_PLL_bandwidth_4xf32(Bl,Tb);
const __m128 sqrt2G = _mm_sub_ps(_mm_set1_ps(0.7026423672846755428877589264195f),_mm_sqrt_ps(_mm_add_ps(G,G)));
const __m128 factor = _mm_add_ps(G,sqrt2G);
return (_mm_mul_ps(_mm_div_ps(Eb,N0),factor));
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
MSK_param_a2_8_65_2xf64(const __m128d Ac,const __m128d Ts, 
                        const __m128d M, const __m128d N0,
                        const __m128d Bl,const __m128d Tb)
{
const __m128d Eb = avg_bit_E_to_carrier_A_2xf64(Ac,Ts,M);
const __m128d G  = G_PLL_bandwidth_2xf64(Bl,Tb);
const __m128d sqrt2G = _mm_sub_pd(_mm_set1_pd(0.7026423672846755428877589264195),_mm_sqrt_pd(_mm_add_pd(G,G)));
const __m128d factor = _mm_add_pd(G,sqrt2G);
return (_mm_mul_pd(_mm_div_pd(Eb,N0),factor));
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
MSK_param_b2_8_65_4xf32(const __m128 Ac,const __m128 Ts, 
                        const __m128 M, const __m128 N0,
                        const __m128 Bl,const __m128 Tb)
{
const __m128 Eb = avg_bit_E_to_carrier_A_4xf32(Ac,Ts,M);
const __m128 G  = G_PLL_bandwidth_4xf32(Bl,Tb);
const __m128 sqrt2G = _mm_add_ps(_mm_set1_ps(0.7026423672846755428877589264195f),_mm_sqrt_ps(_mm_add_ps(G,G)));
const __m128 factor = _mm_add_ps(G,sqrt2G);
return (_mm_mul_ps(_mm_div_ps(Eb,N0),factor));
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
MSK_param_b2_8_65_2xf64(const __m128d Ac,const __m128d Ts, 
                        const __m128d M, const __m128d N0,
                        const __m128d Bl,const __m128d Tb)
{
const __m128d Eb = avg_bit_E_to_carrier_A_2xf64(Ac,Ts,M);
const __m128d G  = G_PLL_bandwidth_2xf64(Bl,Tb);
const __m128d sqrt2G = _mm_add_pd(_mm_set1_pd(0.7026423672846755428877589264195),_mm_sqrt_pd(_mm_add_pd(G,G)));
const __m128d factor = _mm_add_pd(G,sqrt2G);
return (_mm_mul_pd(_mm_div_pd(Eb,N0),factor));
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
template<Gaussian_Q_approximations_sse_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline 
__m128 analytic_SEP_MAM_8_1_4xf32(const __m128 Ac,const __m128 Ts, 
                                  const __m128 M, const __m128 N0
                                  const std::int32_t n) 
{
__m128 result;
__m128 Q_func_val;
__m128 M_ratio   = _mm_div_ps(_mm_sub_ps(M,_mm_set1_ps(1.0)),M);
__m128 left_term = _mm_add_ps(M_ratio,M_ratio);
__m128 tmp       = _mm_mul_ps(Ac,Ac);
__m128 num       = _mm_mul_ps(_mm_add_ps(tmp,tmp),Ts);
__m128 sqr_rat   = _mm_div_ps(num,N0);
__m128 Q_func_arg= _mm_sqrt_ps(sqr_rat);
    if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_chiani)
    {
        Q_func_val = gms::math::gaussian_Q_approx_chiani_4xf32(Q_func_arg);
        result = _mm_mul_ps(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T_4xf32(Q_func_arg);
        result = _mm_mul_ps(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T_4xf32(Q_func_arg);
        result = _mm_mul_ps(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T_4xf32(Q_func_arg);
        result = _mm_mul_ps(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T_4xf32(Q_func_arg);
        result = _mm_mul_ps(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T_4xf32(Q_func_arg);
        result = _mm_mul_ps(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_borjesson)
    {
        Q_func_val = gms::math::gaussian_Q_approx_borjesson_4xf32(Q_func_arg);
        result = _mm_mul_ps(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed_4xf32(Q_func_arg,n);
        result = _mm_mul_ps(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_karagiannidis_lioumpas)
    {
        Q_func_val = gms::math::gaussian_Q_approx_karagiannidis_lioumpas_4xf32(Q_func_arg);
        result = _mm_mul_ps(left_term,Q_func_val);
    }
    return (result);
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
template<Gaussian_Q_approximations_sse_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline 
__m128d analytic_SEP_MAM_8_1_2xf64(const __m128d Ac,const __m128d Ts, 
                                  const __m128d M, const __m128d N0
                                  const std::int32_t n) 
{
__m128d result;
__m128d Q_func_val;
__m128d M_ratio   = _mm_div_pd(_mm_sub_pd(M,_mm_set1_pd(1.0)),M);
__m128d left_term = _mm_add_pd(M_ratio,M_ratio);
__m128d tmp       = _mm_mul_pd(Ac,Ac);
__m128d num       = _mm_mul_pd(_mm_add_pd(tmp,tmp),Ts);
__m128d sqr_rat   = _mm_div_pd(num,N0);
__m128d Q_func_arg= _mm_sqrt_pd(sqr_rat);
    if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_chiani)
    {
        Q_func_val = gms::math::gaussian_Q_approx_chiani_2xf64(Q_func_arg);
        result = _mm_mul_pd(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T_2xf64(Q_func_arg);
        result = _mm_mul_pd(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T_2xf64(Q_func_arg);
        result = _mm_mul_pd(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T_2xf64(Q_func_arg);
        result = _mm_mul_pd(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T_2xf64(Q_func_arg);
        result = _mm_mul_pd(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T_2xf64(Q_func_arg);
        result = _mm_mul_pd(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_borjesson)
    {
        Q_func_val = gms::math::gaussian_Q_approx_borjesson_2xf64(Q_func_arg);
        result = _mm_mul_pd(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed_2xf64(Q_func_arg,n);
        result = _mm_mul_pd(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_karagiannidis_lioumpas)
    {
        Q_func_val = gms::math::gaussian_Q_approx_karagiannidis_lioumpas_2xf64(Q_func_arg);
        result = _mm_mul_pd(left_term,Q_func_val);
    }
    return (result);
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
template<Gaussian_Q_approximations_sse_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline 
__m128 analytic_SEP_MAM_8_3_4xf32(const __m128 Ac,const __m128 Ts, 
                                  const __m128 M, const __m128 N0
                                  const std::int32_t n) 
{
__m128 result;
__m128 Q_func_val;
__m128 one   = _mm_set1_ps(1.0);
__m128 Es    = avg_symbol_E_to_carrier_A_4xf32(Ac,Ts,M);
__m128 sixEs = _mm_mul_ps(_mm_set1_ps(6.0),Es);
__m128 MM_sub_1 = _mm_fmsub_ps(M,M,one);
__m128 tmp      = _mm_div_ps(_mm_sub_ps(M,one),M);
__m128 left_term= _mm_add_ps(tmp,tmp);
__m128 ratio    = _mm_div_ps(Es,_mm_mul_ps(N0,MM_sub_1));
__m128 Q_func_arg = _mm_sqrt_ps(ratio);
    if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_chiani)
    {
        Q_func_val = gms::math::gaussian_Q_approx_chiani_4xf32(Q_func_arg);
        result = _mm_mul_ps(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T_4xf32(Q_func_arg);
        result = _mm_mul_ps(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T_4xf32(Q_func_arg);
        result = _mm_mul_ps(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T_4xf32(Q_func_arg);
        result = _mm_mul_ps(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T_4xf32(Q_func_arg);
        result = _mm_mul_ps(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T_4xf32(Q_func_arg);
        result = _mm_mul_ps(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_borjesson)
    {
        Q_func_val = gms::math::gaussian_Q_approx_borjesson_4xf32(Q_func_arg);
        result = _mm_mul_ps(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed_4xf32(Q_func_arg,n);
        result = _mm_mul_ps(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_karagiannidis_lioumpas)
    {
        Q_func_val = gms::math::gaussian_Q_approx_karagiannidis_lioumpas_4xf32(Q_func_arg);
        result = _mm_mul_ps(left_term,Q_func_val);
    }
    return (result);
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
template<Gaussian_Q_approximations_sse_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline 
__m128d analytic_SEP_MAM_8_3_2xf64(const __m128d Ac,const __m128d Ts, 
                                   const __m128d M, const __m128d N0
                                   const std::int32_t n) 
{
__m128d result;
__m128d Q_func_val;
__m128d one   = _mm_set1_pd(1.0);
__m128d Es    = avg_symbol_E_to_carrier_A_2xf64(Ac,Ts,M);
__m128d sixEs = _mm_mul_pd(_mm_set1_pd(6.0),Es);
__m128d MM_sub_1 = _mm_fmsub_pd(M,M,one);
__m128d tmp      = _mm_div_pd(_mm_sub_pd(M,one),M);
__m128d left_term= _mm_add_pd(tmp,tmp);
__m128d ratio    = _mm_div_pd(Es,_mm_mul_pd(N0,MM_sub_1));
__m128d Q_func_arg = _mm_sqrt_pd(ratio);
    if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_chiani)
    {
        Q_func_val = gms::math::gaussian_Q_approx_chiani_2xf64(Q_func_arg);
        result = _mm_mul_pd(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T_2xf64(Q_func_arg);
        result = _mm_mul_pd(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T_2xf64(Q_func_arg);
        result = _mm_mul_pd(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T_2xf64(Q_func_arg);
        result = _mm_mul_pd(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T_2xf64(Q_func_arg);
        result = _mm_mul_pd(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T_2xf64(Q_func_arg);
        result = _mm_mul_pd(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_borjesson)
    {
        Q_func_val = gms::math::gaussian_Q_approx_borjesson_2xf64(Q_func_arg);
        result = _mm_mul_pd(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed_2xf64(Q_func_arg,n);
        result = _mm_mul_pd(left_term,Q_func_val);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_karagiannidis_lioumpas)
    {
        Q_func_val = gms::math::gaussian_Q_approx_karagiannidis_lioumpas_2xf64(Q_func_arg);
        result = _mm_mul_pd(left_term,Q_func_val);
    }
    return (result);
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
template<Gaussian_Q_approximations_sse_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline 
__m128 analytic_SEP_MAM2_8_4_4xf32(const __m128 Ac,const __m128 Ts, 
                                   const __m128 M, const __m128 N0
                                   const std::int32_t n) 
{
__m128 result;
__m128 Q_func_val;
__m128 Eb    = avg_bit_E_to_carrier_A_4xf32(Ac,Ts,M);
__m128 ratio = _mm_div_ps(_mm_add_ps(Eb,Eb),N0);
__m128 Q_func_arg = _mm_sqrt_ps(ratio);
    if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_chiani)
    {
        Q_func_val = gms::math::gaussian_Q_approx_chiani_4xf32(Q_func_arg);
        result = Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T_4xf32(Q_func_arg);
        result = Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T_4xf32(Q_func_arg);
        result = Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T_4xf32(Q_func_arg);
        result = Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T_4xf32(Q_func_arg);
        result = Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T_4xf32(Q_func_arg);
        result = Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_borjesson)
    {
        Q_func_val = gms::math::gaussian_Q_approx_borjesson_4xf32(Q_func_arg);
        result = Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed_4xf32(Q_func_arg,n);
        result = Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_karagiannidis_lioumpas)
    {
        Q_func_val = gms::math::gaussian_Q_approx_karagiannidis_lioumpas_4xf32(Q_func_arg);
        result = Q_func_val;
    }
    return (result);
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
template<Gaussian_Q_approximations_sse_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline 
__m128d analytic_SEP_MAM2_8_4_2xf64(const __m128d Ac,const __m128d Ts, 
                                   const __m128d M, const __m128d N0
                                   const std::int32_t n) 
{
__m128d result;
__m128d Q_func_val;
__m128d Eb    = avg_bit_E_to_carrier_A_2xf64(Ac,Ts,M);
__m128d ratio = _mm_div_pd(_mm_add_pd(Eb,Eb),N0);
__m128d Q_func_arg = _mm_sqrt_pd(ratio);
    if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_chiani)
    {
        Q_func_val = gms::math::gaussian_Q_approx_chiani_2xf64(Q_func_arg);
        result = Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T_2xf64(Q_func_arg);
        result = Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T_2xf64(Q_func_arg);
        result = Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T_2xf64(Q_func_arg);
        result = Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T_2xf64(Q_func_arg);
        result = Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T_2xf64(Q_func_arg);
        result = Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_borjesson)
    {
        Q_func_val = gms::math::gaussian_Q_approx_borjesson_2xf64(Q_func_arg);
        result = Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed_2xf64(Q_func_arg,n);
        result = Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_karagiannidis_lioumpas)
    {
        Q_func_val = gms::math::gaussian_Q_approx_karagiannidis_lioumpas_2xf64(Q_func_arg);
        result = Q_func_val;
    }
    return (result);
}

#define ANALYTIC_SEP_QAM_8_10_4XF32_COMMON_BLOCK\
    Q_func_val2 = _mm_mul_ps(Q_func_val,Q_func_val);\
    left_term_p2= _mm_mul_ps(left_term,left_term);\
    diff        = _mm_fmsub_ps(left_term,Q_func_val,_mm_mul_ps(left_term_p2,Q_func_val2));\
    result      = diff;

#define ANALYTIC_SEP_QAM_8_10_2XF64_COMMON_BLOCK\
    Q_func_val2 = _mm_mul_pd(Q_func_val,Q_func_val);\
    left_term_p2= _mm_mul_pd(left_term,left_term);\
    diff        = _mm_fmsub_pd(left_term,Q_func_val,_mm_mul_pd(left_term_p2,Q_func_val2));\
    result      = diff;

#if (ANALYTIC_BEP_SEP_CH8_SSE_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=sse
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif
template<Gaussian_Q_approximations_sse_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline 
__m128 analytic_SEP_QAM_8_10_4xf32(const __m128 Ac,const __m128 Ts,
                                   const __m128 M,const __m128 N0,
                                   const std::int32_t n)
{
__m128 result;
__m128 Q_func_val;
__m128 Q_func_val2;
__m128 left_term_p2;
__m128 diff;
__m128 one   = _mm_set1_ps(1.0f);
__m128 Es    = avg_symbol_E_to_carrier_A_4xf32(Ac,Ts,M);
__m128 M_m_1 = _mm_sub_ps(M,one);
__m128 Es3   = _mm_add_ps(Es,_mm_add_ps(Es,Es));
__m128 sqrtM = _mm_sqrt_ps(M);
__m128 tmp   = _mm_div_ps(_mm_sub_ps(sqrtM,one),sqrtM);
__m128 left_term = _mm_mul_ps(_mm_set1_ps(4.0f),tmp);
__m128 ratio = _mm_div_ps(Es3,_mm_mul_ps(N0,M_m_1));
__m128 Q_func_arg = _mm_sqrt_ps(ratio);
    if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_chiani)
    {
        Q_func_val = gms::math::gaussian_Q_approx_chiani_4xf32(Q_func_arg);
        ANALYTIC_SEP_QAM_8_10_4XF32_COMMON_BLOCK
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T_4xf32(Q_func_arg);
        ANALYTIC_SEP_QAM_8_10_4XF32_COMMON_BLOCK
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T_4xf32(Q_func_arg);
        ANALYTIC_SEP_QAM_8_10_4XF32_COMMON_BLOCK
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T_4xf32(Q_func_arg);
        ANALYTIC_SEP_QAM_8_10_4XF32_COMMON_BLOCK
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T_4xf32(Q_func_arg);
        ANALYTIC_SEP_QAM_8_10_4XF32_COMMON_BLOCK
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T_4xf32(Q_func_arg);
        ANALYTIC_SEP_QAM_8_10_4XF32_COMMON_BLOCK
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_borjesson)
    {
        Q_func_val = gms::math::gaussian_Q_approx_borjesson_4xf32(Q_func_arg);
        ANALYTIC_SEP_QAM_8_10_4XF32_COMMON_BLOCK
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed_4xf32(Q_func_arg,n);
        ANALYTIC_SEP_QAM_8_10_4XF32_COMMON_BLOCK
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_karagiannidis_lioumpas)
    {
        Q_func_val = gms::math::gaussian_Q_approx_karagiannidis_lioumpas_4xf32(Q_func_arg);
        ANALYTIC_SEP_QAM_8_10_4XF32_COMMON_BLOCK
    }
    return (result);
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
template<Gaussian_Q_approximations_sse_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline 
__m128d analytic_SEP_QAM_8_10_2xf64(const __m128d Ac,const __m128d Ts,
                                   const __m128d M,const __m128d N0,
                                   const std::int32_t n)
{
__m128d result;
__m128d Q_func_val;
__m128d Q_func_val2;
__m128d left_term_p2;
__m128d diff;
__m128d one   = _mm_set1_pd(1.0);
__m128d Es    = avg_symbol_E_to_carrier_A_2xf64(Ac,Ts,M);
__m128d M_m_1 = _mm_sub_pd(M,one);
__m128d Es3   = _mm_add_pd(Es,_mm_add_pd(Es,Es));
__m128d sqrtM = _mm_sqrt_pd(M);
__m128d tmp   = _mm_div_pd(_mm_sub_pd(sqrtM,one),sqrtM);
__m128d left_term = _mm_mul_pd(_mm_set1_pd(4.0),tmp);
__m128d ratio = _mm_div_pd(Es3,_mm_mul_pd(N0,M_m_1));
__m128d Q_func_arg = _mm_sqrt_pd(ratio);
    if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_chiani)
    {
        Q_func_val = gms::math::gaussian_Q_approx_chiani_2xf64(Q_func_arg);
        ANALYTIC_SEP_QAM_8_10_2XF64_COMMON_BLOCK
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T_2xf64(Q_func_arg);
        ANALYTIC_SEP_QAM_8_10_2XF64_COMMON_BLOCK
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T_2xf64(Q_func_arg);
        ANALYTIC_SEP_QAM_8_10_2XF64_COMMON_BLOCK
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T_2xf64(Q_func_arg);
        ANALYTIC_SEP_QAM_8_10_2XF64_COMMON_BLOCK
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T_2xf64(Q_func_arg);
        ANALYTIC_SEP_QAM_8_10_2XF64_COMMON_BLOCK
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T_2xf64(Q_func_arg);
        ANALYTIC_SEP_QAM_8_10_2XF64_COMMON_BLOCK
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_borjesson)
    {
        Q_func_val = gms::math::gaussian_Q_approx_borjesson_2xf64(Q_func_arg);
        ANALYTIC_SEP_QAM_8_10_2XF64_COMMON_BLOCK
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed_2xf64(Q_func_arg,n);
        ANALYTIC_SEP_QAM_8_10_2XF64_COMMON_BLOCK
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_karagiannidis_lioumpas)
    {
        Q_func_val = gms::math::gaussian_Q_approx_karagiannidis_lioumpas_2xf64(Q_func_arg);
        ANALYTIC_SEP_QAM_8_10_2XF64_COMMON_BLOCK
    }
    return (result);
}

#define COMMON_IF_CONSTEXPR_FUNC_4XF32_BLOCK_SINGLE_LINE\
    if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_chiani)\
    {\
        Q_func_val = gms::math::gaussian_Q_approx_chiani_4xf32(Q_func_arg);\
    }\
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_loskot_2T)\
    {\
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T_4xf32(Q_func_arg);\
    }\
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_loskot_3T)\
    {\
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T_4xf32(Q_func_arg);\
    }\
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_1T)\
    {\
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T_4xf32(Q_func_arg);\
    }\
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_2T)\
    {\
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T_4xf32(Q_func_arg);\
    }\
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_4T)\
    {\
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T_4xf32(Q_func_arg);\
    }\
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_borjesson)\
    {\
        Q_func_val = gms::math::gaussian_Q_approx_borjesson_4xf32(Q_func_arg);\
    }\
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_summed)\
    {\
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed_4xf32(Q_func_arg,n);\
    }\
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_karagiannidis_lioumpas)\
    {\
        Q_func_val = gms::math::gaussian_Q_approx_karagiannidis_lioumpas_4xf32(Q_func_arg);\
    }  

#define COMMON_IF_CONSTEXPR_FUNC_2XF64_BLOCK_SINGLE_LINE\
    if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_chiani)\
    {\
        Q_func_val = gms::math::gaussian_Q_approx_chiani_2xf64(Q_func_arg);\
    }\
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_loskot_2T)\
    {\
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T_2xf64(Q_func_arg);\
    }\
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_loskot_3T)\
    {\
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T_2xf64(Q_func_arg);\
    }\
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_1T)\
    {\
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T_2xf64(Q_func_arg);\
    }\
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_2T)\
    {\
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T_2xf64(Q_func_arg);\
    }\
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_4T)\
    {\
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T_2xf64(Q_func_arg);\
    }\
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_borjesson)\
    {\
        Q_func_val = gms::math::gaussian_Q_approx_borjesson_2xf64(Q_func_arg);\
    }\
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_summed)\
    {\
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed_2xf64(Q_func_arg,n);\
    }\
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_karagiannidis_lioumpas)\
    {\
        Q_func_val = gms::math::gaussian_Q_approx_karagiannidis_lioumpas_2xf64(Q_func_arg);\
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
template<Gaussian_Q_approximations_sse_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline 
__m128 analytic_SEP_QAM4_8_11_4xf32(const __m128 Ac,const __m128 Ts,
                                    const __m128 M, const __m128 N0,
                                    const std::int32_t n)  
{
__m128 result;
__m128 Q_func_val;
__m128 Q_func_val_a2;
__m128 Q_func_val_p2;
__m128 diff_result;
__m128 Es         = avg_symbol_E_to_carrier_A_4xf32(Ac,Ts,M);
__m128 sqrt_arg   = _mm_div_ps(Es,N0);
__m128 Q_func_arg = _mm_sqrt_ps(sqrt_arg);
    COMMON_IF_CONSTEXPR_FUNC_4XF32_BLOCK_SINGLE_LINE
    Q_func_val_a2 = _mm_add_ps(Q_func_val,Q_func_val);
    Q_func_val_p2 = _mm_mul_ps(Q_func_val,Q_func_val);
    diff_result   = _mm_sub_ps(Q_func_val_a2,Q_func_val_p2);
    result = diff_result;
    return (result);
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
template<Gaussian_Q_approximations_sse_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline 
__m128d analytic_SEP_QAM4_8_11_2xf64(const __m128d Ac,const __m128d Ts,
                                    const __m128d M, const __m128d N0,
                                    const std::int32_t n)  
{
__m128d result;
__m128d Q_func_val;
__m128d Q_func_val_a2;
__m128d Q_func_val_p2;
__m128d diff_result;
__m128d Es         = avg_symbol_E_to_carrier_A_2xf64(Ac,Ts,M);
__m128d sqrt_arg   = _mm_div_pd(Es,N0);
__m128d Q_func_arg = _mm_sqrt_pd(sqrt_arg);
    COMMON_IF_CONSTEXPR_FUNC_2XF64_BLOCK_SINGLE_LINE
    Q_func_val_a2 = _mm_add_pd(Q_func_val,Q_func_val);
    Q_func_val_p2 = _mm_mul_pd(Q_func_val,Q_func_val);
    diff_result   = _mm_sub_pd(Q_func_val_a2,Q_func_val_p2);
    result = diff_result;
    return (result);
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
template<Gaussian_Q_approximations_sse_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline 
__m128 analytic_BEP_QAM_AWGN_8_14_4xf32(const __m128 Ac,const __m128 Ts, 
                                        const __m128 M, const __m128 N0,
                                        const std::int32_t n)  
{
__m128 result;
__m128 Q_func_val  = _mm_setzero_ps(); 
__m128 one         = _mm_set1_ps(1.0f);
__m128 Eb          = avg_bit_E_to_carrier_A_4xf32(Ac,Ts,M);
__m128 log2M       = _mm_log2_ps(M);
__m128 invlog2M    = _mm_rcp_ps(log2M);
__m128 sqrtM       = _mm_sqrt_ps(M); 
__m128 sqrtM_sub1  = _mm_sub_ps(sqrtM,one);
__m128i itmp       = _mm_cvtps_epi32(_mm_mul_ps(sqrtM,_mm_set1_ps(0.5f)));
std::int32_t up_lim= _mm_extract_epi32(itmp,0);
if(__builtin_expect(0==up_lim,0)) { return _mm_set1_ps(std::numeric_limits<float>::quiet_NaN());}
__m128 left_term   = _mm_mul_ps(_mm_set1_ps(4.0f),_mm_div_ps(sqrtM_sub1,sqrtM));
__m128 num         = _mm_mul_ps(_mm_set1_ps(3.0),_mm_mul_ps(Eb,log2M));
__m128 den         = _mm_mul_ps(N0,sqrtM_sub1);
__m128 Q_func_arg  = _mm_sqrt_ps(_mm_div_ps(num,den));
std::int32_t rem   = up_lim%4;
if(0==rem)
{
    __ATTR_ALIGN__(16) float lut_idx[16] = 
    {
         1.0f,2.0f,3.0f,4.0f,
         5.0f,6.0f,7.0f,8.0f,
         9.0f,10.0f,11.0f,12.0f,
         13.0f,14.0f,15.0f,16.0f
    };
    if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_chiani)
    {
        for(std::int32_t i = 1; i<=ROUND_TO_FOUR(up_lim,3); i+= 4)
        {
            __m128 indices = _mm_load_ps(&lut_idx[i]);
            __m128 mul_fac = _mm_sub_ps(_mm_add_ps(indices,indices),one);
            __m128 arg     = _mm_mul_ps(mul_fac,Q_func_arg);
            Q_func_val     = _mm_add_ps(Q_func_val,gms::math::gaussian_Q_approx_chiani_4xf32(arg));
        }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_loskot_2T)
    {
        for(std::int32_t i = 1; i<=ROUND_TO_FOUR(up_lim,3); i+= 4)
        {
            __m128 indices = _mm_load_ps(&lut_idx[i]);
            __m128 mul_fac = _mm_sub_ps(_mm_add_ps(indices,indices),one);
            __m128 arg     = _mm_mul_ps(mul_fac,Q_func_arg);
            Q_func_val     = _mm_add_ps(Q_func_val,gms::math::gaussian_Q_approx_loskot_2T_4xf32(arg));
        }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_loskot_3T)
    {
        for(std::int32_t i = 1; i<=ROUND_TO_FOUR(up_lim,3); i+= 4)
        {
            __m128 indices = _mm_load_ps(&lut_idx[i]);
            __m128 mul_fac = _mm_sub_ps(_mm_add_ps(indices,indices),one);
            __m128 arg     = _mm_mul_ps(mul_fac,Q_func_arg);
            Q_func_val     = _mm_add_ps(Q_func_val,gms::math::gaussian_Q_approx_loskot_3T_4xf32(arg));
        }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_1T)
    {
        for(std::int32_t i = 1; i<=ROUND_TO_FOUR(up_lim,3); i+= 4)
        {
            __m128 indices = _mm_load_ps(&lut_idx[i]);
            __m128 mul_fac = _mm_sub_ps(_mm_add_ps(indices,indices),one);
            __m128 arg     = _mm_mul_ps(mul_fac,Q_func_arg);
            Q_func_val     = _mm_add_ps(Q_func_val,gms::math::gaussian_Q_approx_sadhwani_1T_4xf32(arg));
        }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_2T)
    {
        for(std::int32_t i = 1; i<=ROUND_TO_FOUR(up_lim,3); i+= 4)
        {
            __m128 indices = _mm_load_ps(&lut_idx[i]);
            __m128 mul_fac = _mm_sub_ps(_mm_add_ps(indices,indices),one);
            __m128 arg     = _mm_mul_ps(mul_fac,Q_func_arg);
            Q_func_val     = _mm_add_ps(Q_func_val,gms::math::gaussian_Q_approx_sadhwani_2T_4xf32(arg));
        }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_4T)
    {
        for(std::int32_t i = 1; i<=ROUND_TO_FOUR(up_lim,3); i+= 4)
        {
            __m128 indices = _mm_load_ps(&lut_idx[i]);
            __m128 mul_fac = _mm_sub_ps(_mm_add_ps(indices,indices),one);
            __m128 arg     = _mm_mul_ps(mul_fac,Q_func_arg);
            Q_func_val     = _mm_add_ps(Q_func_val,gms::math::gaussian_Q_approx_sadhwani_4T_4xf32(arg));
        }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_borjesson)
    {
        for(std::int32_t i = 1; i<=ROUND_TO_FOUR(up_lim,3); i+= 4)
        {
            __m128 indices = _mm_load_ps(&lut_idx[i]);
            __m128 mul_fac = _mm_sub_ps(_mm_add_ps(indices,indices),one);
            __m128 arg     = _mm_mul_ps(mul_fac,Q_func_arg);
            Q_func_val     = _mm_add_ps(Q_func_val,gms::math::gaussian_Q_approx_borjesson_4xf32(arg));
        }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_summed)
    {
        for(std::int32_t i = 1; i<=ROUND_TO_FOUR(up_lim,3); i+= 4)
        {
            __m128 indices = _mm_load_ps(&lut_idx[i]);
            __m128 mul_fac = _mm_sub_ps(_mm_add_ps(indices,indices),one);
            __m128 arg     = _mm_mul_ps(mul_fac,Q_func_arg);
            Q_func_val     = _mm_add_ps(Q_func_val,gms::math::gaussian_Q_approx_sadhwani_summed_4xf32(arg,n));
        }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_karagiannidis_lioumpas)
    {
        for(std::int32_t i = 1; i<=ROUND_TO_FOUR(up_lim,3); i+= 4)
        {
            __m128 indices = _mm_load_ps(&lut_idx[i]);
            __m128 mul_fac = _mm_sub_ps(_mm_add_ps(indices,indices),one);
            __m128 arg     = _mm_mul_ps(mul_fac,Q_func_arg);
            Q_func_val     = _mm_add_ps(Q_func_val,gms::math::gaussian_Q_approx_karagiannidis_lioumpas_4xf32(arg));
        }
    }
    result = _mm_mul_ps(left_term,_mm_mul_ps(invlog2M,Q_func_arg));
    return (result);
}
else
{
    result = _mm_set1_ps(std::numeric_limits<float>::quiet_NaN());
    return (result);
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
template<Gaussian_Q_approximations_sse_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline 
__m128d analytic_BEP_QAM_AWGN_8_14_2xf64(const __m128d Ac,const __m128d Ts, 
                                         const __m128d M, const __m128d N0,
                                         const std::int32_t n)  
{
__m128d result;
__m128d Q_func_val  = _mm_setzero_pd(); 
__m128d one         = _mm_set1_pd(1.0f);
__m128d Eb          = avg_bit_E_to_carrier_A_2xf64(Ac,Ts,M);
__m128d log2M       = _mm_log2_pd(M);
__m128d invlog2M    = _mm_rcp14_pd(log2M);
__m128d sqrtM       = _mm_sqrt_pd(M); 
__m128d sqrtM_sub1  = _mm_sub_pd(sqrtM,one);
__m128i itmp       = _mm_cvtpd_epi32(_mm_mul_pd(sqrtM,_mm_set1_pd(0.5)));
std::int32_t up_lim= _mm_extract_epi32(itmp,0);
if(__builtin_expect(0==up_lim,0)) { return _mm_set1_pd(std::numeric_limits<double>::quiet_NaN());}
__m128d left_term   = _mm_mul_pd(_mm_set1_pd(4.0),_mm_div_pd(sqrtM_sub1,sqrtM));
__m128d num         = _mm_mul_pd(_mm_set1_pd(3.0),_mm_mul_pd(Eb,log2M));
__m128d den         = _mm_mul_pd(N0,sqrtM_sub1);
__m128d Q_func_arg  = _mm_sqrt_pd(_mm_div_pd(num,den));
std::int32_t rem   = up_lim%2;
if(0==rem)
{
    __ATTR_ALIGN__(16) double lut_idx[16] = 
    {
         1.0,2.0,3.0,4.0,
         5.0,6.0,7.0,8.0,
         9.0,10.0,11.0,12.0,
         13.0,14.0,15.0,16.0
    };
    _mm_prefetch((const char*)&lut_idx[0],_MM_HINT_T0);
    _mm_prefetch((const char*)&lut_idx[64],_MM_HINT_T0);
    if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_chiani)
    {
        for(std::int32_t i = 1; i<=ROUND_DOWN(up_lim,1); i+= 2)
        {
            __m128d indices = _mm_load_pd(&lut_idx[i]);
            __m128d mul_fac = _mm_sub_pd(_mm_add_ps(indices,indices),one);
            __m128d arg     = _mm_mul_pd(mul_fac,Q_func_arg);
            Q_func_val     = _mm_add_pd(Q_func_val,gms::math::gaussian_Q_approx_chiani_2xf64(arg));
        }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_loskot_2T)
    {
        for(std::int32_t i = 1; i<=ROUND_DOWN(up_lim,1); i+= 2)
        {
            __m128d indices = _mm_load_pd(&lut_idx[i]);
            __m128d mul_fac = _mm_sub_pd(_mm_add_pd(indices,indices),one);
            __m128d arg     = _mm_mul_pd(mul_fac,Q_func_arg);
            Q_func_val      = _mm_add_pd(Q_func_val,gms::math::gaussian_Q_approx_loskot_2T_2xf64(arg));
        }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_loskot_3T)
    {
        for(std::int32_t i = 1; i<=ROUND_DOWN(up_lim,1); i+= 2)
        {
            __m128d indices = _mm_load_pd(&lut_idx[i]);
            __m128d mul_fac = _mm_sub_pd(_mm_add_pd(indices,indices),one);
            __m128d arg     = _mm_mul_pd(mul_fac,Q_func_arg);
            Q_func_val      = _mm_add_pd(Q_func_val,gms::math::gaussian_Q_approx_loskot_3T_2xf64(arg));
        }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_1T)
    {
        for(std::int32_t i = 1; i<=ROUND_DOWN(up_lim,1); i+= 2)
        {
            __m128d indices = _mm_load_pd(&lut_idx[i]);
            __m128d mul_fac = _mm_sub_pd(_mm_add_pd(indices,indices),one);
            __m128d arg     = _mm_mul_pd(mul_fac,Q_func_arg);
            Q_func_val     = _mm_add_pd(Q_func_val,gms::math::gaussian_Q_approx_sadhwani_1T_2xf64(arg));
        }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_2T)
    {
        for(std::int32_t i = 1; i<=ROUND_DOWN(up_lim,1); i+= 2)
        {
            __m128d indices = _mm_load_pd(&lut_idx[i]);
            __m128d mul_fac = _mm_sub_pd(_mm_add_ps(indices,indices),one);
            __m128d arg     = _mm_mul_pd(mul_fac,Q_func_arg);
            Q_func_val     = _mm_add_pd(Q_func_val,gms::math::gaussian_Q_approx_sadhwani_2T_2xf64(arg));
        }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_4T)
    {
        for(std::int32_t i = 1; i<=ROUND_DOWN(up_lim,1); i+= 2)
        {
            __m128d indices = _mm_load_pd(&lut_idx[i]);
            __m128d mul_fac = _mm_sub_pd(_mm_add_ps(indices,indices),one);
            __m128d arg     = _mm_mul_pd(mul_fac,Q_func_arg);
            Q_func_val     = _mm_add_pd(Q_func_val,gms::math::gaussian_Q_approx_sadhwani_4T_2xf64(arg));
        }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_borjesson)
    {
        for(std::int32_t i = 1; i<=ROUND_DOWN(up_lim,1); i+= 2)
        {
            __m128d indices = _mm_load_pd(&lut_idx[i]);
            __m128d mul_fac = _mm_sub_pd(_mm_add_pd(indices,indices),one);
            __m128d arg     = _mm_mul_pd(mul_fac,Q_func_arg);
            Q_func_val     = _mm_add_pd(Q_func_val,gms::math::gaussian_Q_approx_borjesson_2xf64(arg));
        }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_summed)
    {
        for(std::int32_t i = 1; i<=ROUND_DOWN(up_lim,1); i+= 2)
        {
            __m128d indices = _mm_load_pd(&lut_idx[i]);
            __m128d mul_fac = _mm_sub_pd(_mm_add_ps(indices,indices),one);
            __m128d arg     = _mm_mul_pd(mul_fac,Q_func_arg);
            Q_func_val     = _mm_add_pd(Q_func_val,gms::math::gaussian_Q_approx_sadhwani_summed_2xf64(arg,n));
        }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_karagiannidis_lioumpas)
    {
        for(std::int32_t i = 1; i<=ROUND_DOWN(up_lim,1); i+= 2)
        {
            __m128d indices = _mm_load_pd(&lut_idx[i]);
            __m128d mul_fac = _mm_sub_pd(_mm_add_pd(indices,indices),one);
            __m128d arg     = _mm_mul_pd(mul_fac,Q_func_arg);
            Q_func_val     = _mm_add_pd(Q_func_val,gms::math::gaussian_Q_approx_karagiannidis_lioumpas_2xf64(arg));
        }
    }
    result = _mm_mul_pd(left_term,_mm_mul_pd(invlog2M,Q_func_arg));
    return (result);
}
else
{
    result = _mm_set1_ps(std::numeric_limits<double>::quiet_NaN());
    return (result);
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
template<Gaussian_Q_approximations_sse_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline 
__m128 analytic_BEP_MPSK2_8_18_4xf32(const __m128 Ac,const __m128 Ts,
                                     const __m128 M, const __m128 N0,
                                     const std::int32_t n)  
{
__m128 result;
__m128 Q_func_val;
__m128 Eb   = avg_bit_E_to_carrier_A_4xf32(Ac,Ts,M);
__m128 ratio= _mm_div_ps(_mm_add_ps(Eb,Eb),N0);
__m128 Q_func_arg = _mm_sqrt_ps(ratio);
    COMMON_IF_CONSTEXPR_FUNC_4XF32_BLOCK_SINGLE_LINE
result = Q_func_val;
return (result);
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
template<Gaussian_Q_approximations_sse_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline 
__m128d analytic_BEP_MPSK2_8_18_2xf64(const __m128d Ac,const __m128d Ts,
                                     const __m128d M, const __m128d N0,
                                     const std::int32_t n)  
{
__m128d result;
__m128d Q_func_val;
__m128d Eb   = avg_bit_E_to_carrier_A_2xf64(Ac,Ts,M);
__m128d ratio= _mm_div_pd(_mm_add_pd(Eb,Eb),N0);
__m128d Q_func_arg = _mm_sqrt_pd(ratio);
    COMMON_IF_CONSTEXPR_FUNC_2XF64_BLOCK_SINGLE_LINE
result = Q_func_val;
return (result);
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
template<Gaussian_Q_approximations_sse_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline 
__m128 analytic_SEP_QPSKM4_8_19_4xf32(const __m128 Ac,const __m128 Ts,
                                      const __m128 M,const __m128 N0,
                                      const std::int32_t n) 
{
__m128 result;
__m128 Q_func_val;
__m128 sqrQ;
__m128 Es  = avg_symbol_E_to_carrier_A_4xf32(Ac,Ts,M);
__m128 snr = Es/N0;
__m128 Q_func_arg = _mm_sqrt_ps(snr);
    if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_chiani)
    {
        Q_func_val = gms::math::gaussian_Q_approx_chiani_4xf32(Q_func_arg);
        sqrQ       = _mm_mul_ps(Q_func_val,Q_func_val);
        result     = _mm_sub_ps(_mm_add_ps(Q_func_val,Q_func_val),sqrQ);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T_4xf32(Q_func_arg);
        sqrQ       = _mm_mul_ps(Q_func_val,Q_func_val);
        result     = _mm_sub_ps(_mm_add_ps(Q_func_val,Q_func_val),sqrQ);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T_4xf32(Q_func_arg);
        sqrQ       = _mm_mul_ps(Q_func_val,Q_func_val);
        result     = _mm_sub_ps(_mm_add_ps(Q_func_val,Q_func_val),sqrQ);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T_4xf32(Q_func_arg);
        sqrQ       = _mm_mul_ps(Q_func_val,Q_func_val);
        result     = _mm_sub_ps(_mm_add_ps(Q_func_val,Q_func_val),sqrQ);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T_4xf32(Q_func_arg);
        sqrQ       = _mm_mul_ps(Q_func_val,Q_func_val);
        result     = _mm_sub_ps(_mm_add_ps(Q_func_val,Q_func_val),sqrQ);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T_4xf32(Q_func_arg);
        sqrQ       = _mm_mul_ps(Q_func_val,Q_func_val);
        result     = _mm_sub_ps(_mm_add_ps(Q_func_val,Q_func_val),sqrQ);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_borjesson)
    {
        Q_func_val = gms::math::gaussian_Q_approx_borjesson_4xf32(Q_func_arg);
        sqrQ       = _mm_mul_ps(Q_func_val,Q_func_val);
        result     = _mm_sub_ps(_mm_add_ps(Q_func_val,Q_func_val),sqrQ);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed_4xf32(Q_func_arg,n);
        sqrQ       = _mm_mul_ps(Q_func_val,Q_func_val);
        result     = _mm_sub_ps(_mm_add_ps(Q_func_val,Q_func_val),sqrQ);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approximations_sse_t::Gaussian_Q_approx_karagiannidis_lioumpas)
    {
        Q_func_val = gms::math::gaussian_Q_approx_karagiannidis_lioumpas_4xf32(Q_func_arg);
        sqrQ       = _mm_mul_ps(Q_func_val,Q_func_val);
        result     = _mm_sub_ps(_mm_add_ps(Q_func_val,Q_func_val),sqrQ);
    }
    return (result);
}

}
}
#endif /*__GMS_ANALYTIC_BEP_SEP_CH8_SSE_H__*/