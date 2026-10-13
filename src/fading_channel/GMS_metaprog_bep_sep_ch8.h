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

#ifndef __GMS_METAPROG_BEP_SEP_CH8_H__
#define __GMS_METAPROG_BEP_SEP_CH8_H__ 131020260825

#include <cstdint>
#include <immintrin.h>
#include <cmath>
#include "GMS_config.h"
#include "GMS_fast_simd_funcs_approx.h"

namespace file_info 
{

     static const unsigned int GMS_METAPROG_BEP_SEP_CH8_MAJOR = 1;
     static const unsigned int GMS_METAPROG_BEP_SEP_CH8_MINOR = 1;
     static const unsigned int GMS_METAPROG_BEP_SEP_CH8_MICRO = 0;
     static const unsigned int GMS_METAPROG_BEP_SEP_CH8_FULLVER =
       1000U*GMS_METAPROG_BEP_SEP_CH8_MAJOR+100U*GMS_METAPROG_BEP_SEP_CH8_MINOR+
       10U*GMS_METAPROG_BEP_SEP_CH8_MICRO;
     static const char GMS_METAPROG_BEP_SEP_CH8_CREATION_DATE[] = "13-10-2026 08:25PM +00200 (TUE 13 OCT 2026 GMT+2)";
     static const char GMS_METAPROG_BEP_SEP_CH8_BUILD_DATE[]    = __DATE__; 
     static const char GMS_METAPROG_BEP_SEP_CH8_BUILD_TIME[]    = __TIME__;
     static const char GMS_METAPROG_BEP_SEP_CH8_AUTHOR[]        = "Programmer: Bernard Gingold, beniekg@gmail.com";
     static const char GMS_METAPROG_BEP_SEP_CH8_SYNOPSIS[]      = "Analytic Formulae of the Bit and Symbol Error Probabilities based on the approximated Gaussian-Q function.\
	                                                          Based on the chapter 8 of M.K Simon, M.S. Alouini: Digital Communication over Fading Channels 1st ed\
															  ISBN-13 978-0471317791";

}

namespace gms
{

namespace fading_channel
{

template<std::size_t Index, std::size_t BatchSize>
struct Avg_symbol_E_to_carrier_A_4xf32_t
{
static inline void 
compile_time_execute(const float * __restrict__ Ac,
                     const float * __restrict__ Ts,
                     const float * __restrict__ M,
                     float       * __restrict__ result)
{
    if constexpr(Index < BatchSize)
    {
        constexpr std::size_t idx = Index*4ull;
        __m128 vAc   = _mm_loadu_ps(Ac+idx);
        __m128 vTs   = _mm_loadu_ps(Ts+idx);
        __m128 vM    = _mm_loadu_ps(M+idx);
        __m128 third = _mm_set1_ps(0.3333333333f);
        __m128 AA    = _mm_mul_ps(vAc,vAc);
        __m128 one   = _mm_set1_ps(1.0f);
        __m128 MM    = _mm_mul_ps(vM,vM);
        __m128 MM_rat= _mm_sub_ps(_mm_mul_ps(MM,one),third);
        __m128 factor= _mm_mul_ps(AA,vTs); 
        __m128 vres  = _mm_mul_ps(factor,MM_rat);
        _mm_storeu_ps(result+idx,vres);

        Avg_symbol_E_to_carrier_A_4xf32_t<Index+1ull,BatchSize>::compile_time_execute(Ac,Ts,M,result);
    }
}

};

template<std::size_t Index, std::size_t BatchSize>
struct Avg_symbol_E_to_carrier_A_2xf64_t
{
static inline void 
compile_time_execute(const double * __restrict__ Ac,
                     const double * __restrict__ Ts,
                     const double * __restrict__ M,
                     float       * __restrict__ result)
{
    if constexpr(Index < BatchSize)
    {
        constexpr std::size_t idx = Index*2ull;
        __m128d vAc   = _mm_loadu_pd(Ac+idx);
        __m128d vTs   = _mm_loadu_pd(Ts+idx);
        __m128d vM    = _mm_loadu_pd(M+idx);
        __m128d third = _mm_set1_pd(0.333333333333333333333);
        __m128d AA    = _mm_mul_pd(vAc,vAc);
        __m128d one   = _mm_set1_pd(1.0);
        __m128d MM    = _mm_mul_pd(vM,vM);
        __m128d MM_rat= _mm_sub_pd(_mm_mul_pd(MM,one),third);
        __m128d factor= _mm_mul_pd(AA,vTs); 
        __m128d vres  = _mm_mul_pd(factor,MM_rat);
        _mm_storeu_pd(result+idx,vres);

        Avg_symbol_E_to_carrier_A_2xf64_t<Index+1ull,BatchSize>::compile_time_execute(Ac,Ts,M,result);
    }
}

};

template<std::size_t Index, std::size_t BatchSize>
struct Avg_bit_E_to_carrier_A_4xf32_t
{
static inline void 
compile_time_execute(const float * __restrict__ Ac,
                     const float * __restrict__ Ts,
                     const float * __restrict__ M,
                     float       * __restrict__ result)
{
    if constexpr(Index < BatchSize)
    {
        constexpr std::size_t idx = Index*4ull;
        __m128 vAc   = _mm_loadu_ps(Ac+idx);
        __m128 vTs   = _mm_loadu_ps(Ts+idx);
        __m128 vM    = _mm_loadu_ps(M+idx);
        __m128 one   = _mm_set1_ps(1.0f);
        __m128 num   = _mm_fmsub_ps(vM,vM,one);
        __m128 Eg    = _mm_mul_ps(_mm_mul_ps(vAc,vAc),vTs);
        __m128 log2M = _mm_mul_ps(_mm_set1_ps(6.0f),_mm_log2_ps(vM));
        __m128 vres  = _mm_mul_ps(Eg,_mm_div_ps(num,log2M));
        _mm_storeu_ps(result+idx,vres);

        Avg_bit_E_to_carrier_A_4xf32_t<Index+1ull,BatchSize>::compile_time_execute(Ac,Ts,M,result);
    }
}
};

template<std::size_t Index, std::size_t BatchSize>
struct Avg_bit_E_to_carrier_A_2xf64_t
{
static inline void 
compile_time_execute(const double * __restrict__ Ac,
                     const double * __restrict__ Ts,
                     const double * __restrict__ M,
                     double       * __restrict__ result)
{
    if constexpr(Index < BatchSize)
    {
        constexpr std::size_t idx = Index*2ull;
        __m128d vAc   = _mm_loadu_pd(Ac+idx);
        __m128d vTs   = _mm_loadu_pd(Ts+idx);
        __m128d vM    = _mm_loadu_pd(M+idx);
        __m128d one   = _mm_set1_pd(1.0);
        __m128d num   = _mm_fmsub_pd(vM,vM,one);
        __m128d Eg    = _mm_mul_pd(_mm_mul_pd(vAc,vAc),vTs);
        __m128d log2M = _mm_mul_pd(_mm_set1_pd(6.0f),_mm_log2_pd(vM));
        __m128d vres  = _mm_mul_pd(Eg,_mm_div_pd(num,log2M));
        _mm_storeu_pd(result+idx,vres);

        Avg_bit_E_to_carrier_A_2xf64_t<Index+1ull,BatchSize>::compile_time_execute(Ac,Ts,M,result);
    }
}
};

template<std::size_t Index, std::size_t BatchSize>
struct Tikhonov_phase_err_pdf_4xf33_t
{
static inline void 
compile_time_execute(const float * __restrict__ rho_eq,
                     const float * __restrict__ phi_c,
                     float       * __restric__ result) 
{
    if constexpr(Index < BatchSize)
    {
        constexpr std::size_t idx = Index*4ull;
        __m128 vrho_eq = _mm_loadu_ps(rho_eq+idx);
        __m128 vphi_c  = _mm_loadu_ps(phi_c+idx);
        __m128 negPI   = _mm_set1_ps(-3.1415926535897932384626433832795f);
        const float * __restrict__ p_phi_c = reinterpret_cast<const float * __restrict__>(&vphi_c);
        __m128 posPI   = _mm_set1_ps(+3.1415926535897932384626433832795);
        __m128 bessi0  = _mm_setr_ps(std::cyl_bessel_i(0,p_phi_c[0],std::cyl_bessel_i(0,p_phi_c[1]),
                                     std::cyl_bessel_i(0,p_phi_c[2]),std::cyl_bessel_i(0,p_phi_c[3])));
        __mmask8 phic_le_pi    = _mm_cmp_ps_mask(vphi_c,posPI,_CMP_LE_OQ);
        __m128  cosarg         = _mm_add_ps(vphi_c,vphi_c);
        __mmask8 phic_gt_pi    = _mm_cmp_ps_mask(vphi_c,negPI,_CMP_GE_OQ);
        __mmask8 mask_pi       = _kand_mask8(phic_le_pi,phic_gt_pi);
        __m128  cos2phic       = _mm_mask_blend_ps(mask_pi,_mm_cos_ps(cosarg),gms::math::simd_fast_cos_approx_4xf32(cosarg));
        __m128  exp_val        = gms::math::simd_fast_exp_approx_4xf32(_mm_mul_ps(vrho_eq,cos2phic));
        __m128 vres            = _mm_div_ps(exp_val,_mm_mul_ps(posPI,bessi0));
        _mm_storeu_ps(result+idx,vres);

        Tikhonov_phase_err_pdf_4xf33_t<Index+1ull,BatchSize>::compile_time_execute(rho_eq,phi_c,result);
    }
}
};



}

}

#endif /*__GMS_METAPROG_BEP_SEP_CH8_H__*/