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
#include "GMS_bessel_i0_sse.h"
#include "GMS_simd_cephes_cos.h"
#include "GMS_simd_cephes_sin.h"

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

namespace 
{
static inline 
__m128 
avg_bit_E_to_carrier_A_4xf32(const __m128 Ac,const __m128 Ts, 
                             const __m128 M)
{
const __m128 num   = _mm_fmsub_ps(M,M,_mm_set1_ps(1.0));
const __m128 Eg    = _mm_mul_ps(_mm_mul_ps(Ac,Ac,),Ts);
const __m128 log2M = _mm_mul_ps(_mm_set1_ps(6.0),_mm_log2_ps(M));
return (_mm_mul_ps(Eg,_mm_div_ps(num,log2M))); 
}

static inline 
__m128d 
avg_bit_E_to_carrier_A_2xf64(const __m128d Ac,const __m128d Ts, 
                             const __m128d M)
{
const __m128d num   = _mm_fmsub_pd(M,M,_mm_set1_pd(1.0));
const __m128d Eg    = _mm_mul_pd(_mm_mul_pd(Ac,Ac,),Ts);
const __m128d log2M = _mm_mul_pd(_mm_set1_pd(6.0),_mm_log2_pd(M));
return (_mm_mul_pd(Eg,_mm_div_pd(num,log2M))); 
}

static inline 
__m128 
avg_symbol_E_to_carrier_A_4xf32(const __m128 Ac,const __m128 Tc,
                                const __m128 M)
{
const __m128 third    = _mm_set1_ps(0.3333333333333333333333333333333333333333f);
const __m128 AA       = _mm_mul_ps(Ac,Ac);
const __m128 MM       = _mm_mul_ps(M,M); 
const __m128 one      = _mm_set1_ps(1.0f);
const __m128 MM_ratio = _mm_sub_ps(_mm_mul_ps(MM,one),third);
const __m128 lead_fac = _mm_mul_ps(AA,Ts);
return (_mm_mul_ps(lead_fac,MM_ratio)); 
}

static inline 
__m128d 
avg_symbol_E_to_carrier_A_2xf64(const __m128d Ac,const __m128d Tc,
                                const __m128d M)
{
const __m128d third    = _mm_set1_pd(0.3333333333333333333333333333333333333333);
const __m128d AA       = _mm_mul_pd(Ac,Ac);
const __m128d MM       = _mm_mul_pd(M,M); 
const __m128d one      = _mm_set1_pd(1.0);
const __m128d MM_ratio = _mm_sub_pd(_mm_mul_pd(MM,one),third);
const __m128d lead_fac = _mm_mul_pd(AA,Ts);
return (_mm_mul_pd(lead_fac,MM_ratio)); 
}

static inline 
__m128
G_PLL_bandwidth_4xf32(const __m128 Bl,const __m128 Tb)
{
return (_mm_rcp_ps(_mm_mul_ps(Bl,Tb)));
}


template<std::int32_t use_div_as_reciprocal>
static inline 
__m128d 
G_PLL_bandwidth_2xf64(const __m128d Bl,const __m128d Tb)
{
if constexpr (use_div_as_reciprocal==1)
{
return (_mm_div_pd(_mm_set1_pd(1.0),_mm_mul_pd(Bl,Tb)));
}
else if constexpr(use_div_as_reciprocal==2)
{ 
return (_mm_rcp14_pd(_mm_mul_pd(Bl,Tb)));
}
}

static inline 
__m128d 
G_PLL_bandwidth_2xf64_v2(const __m128d Bl,const __m128d Tb)
{
return (_mm_div_pd(_mm_set1_pd(1.0),_mm_mul_pd(Bl,Tb)));
}

static inline 
__m128 spec_pow_to_rand_pow_8_145_4xf32(const __m128 m_I,const __m128 m_Q,
                                        const __m128 sigma)
{
__m128 sigma_sqr = _mm_mul_ps(sigma,sigma);
__m128 one       = _mm_set1_ps(1.0f);
__m128 ricc_fac  = _mm_add_ps(one,riccian_factor_8_144_4xf32(m_I,m_Q,sigma));
return (_mm_mul_ps(_mm_add_ps(sigma_sqr,sigma_sqr),ricc_fac));
}

static inline 
__m128d spec_pow_to_rand_pow_8_145_2xf64(const __m128d m_I,const __m128d m_Q,
                                        const __m128d sigma)
{
__m128d sigma_sqr = _mm_mul_pd(sigma,sigma);
__m128d one       = _mm_set1_pd(1.0);
__m128d ricc_fac  = _mm_add_pd(one,riccian_factor_8_144_2xf64(m_I,m_Q,sigma));
return (_mm_mul_pd(_mm_add_pd(sigma_sqr,sigma_sqr),ricc_fac));
}

}

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
struct Tikhonov_phase_err_pdf_4xf32_t
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

        Tikhonov_phase_err_pdf_4xf32_t<Index+1ull,BatchSize>::compile_time_execute(rho_eq,phi_c,result);
    }
}

};

template<std::size_t Index, std::size_t BatchSize>
struct Tikhonov_phase_err_pdf_2xf64_t
{
static inline void 
compile_time_execute(const double * __restrict__ rho_eq,
                     const double * __restrict__ phi_c,
                     double       * __restric__ result)
{
    if constexpr(Index < BatchSize)
    {
        constexpr std::size_t idx = Index*2ull;
        __m128 vrho_eq = _mm_loadu_pd(rho_eq+idx);
        __m128 vphi_c  = _mm_loadu_pd(phi_c+idx);
        __m128d negPI        = _mm_set1_pd(-3.1415926535897932384626433832795);
        __m128d posPI        = _mm_set1_pd(+3.1415926535897932384626433832795);
        __m128d bessi0       = _mm_mask_blend_pd(_mm_cmp_pd_mask(vrho_eq,_mm_set1_pd(15.0),_CMP_GE_OQ),
                                               gms::math::bessel_i0_le15_sse_pd(vrho_eq),
                                               gms::math::bessel_i0_ge15_sse_pd(vrho_eq));
        __mmask8 phic_le_pi   = _mm_cmp_pd_mask(vphi_c,posPI,_CMP_LE_OQ);
        __m128d  cosarg       = _mm_add_pd(vphi_c,vphi_c);
        __mmask8 phic_gt_pi   = _mm_cmp_pd_mask(vphi_c,negPI,_CMP_GE_OQ);
        __mmask8 mask_pi      = _kand_mask8(phic_le_pi,phic_gt_pi);
        __m128d  cos2phic     = _mm_mask_blend_pd(mask_pi,_mm_cos_pd(cosarg),gms::math::simd_fast_cos_approx_2xf64(cosarg));
        __m128d  exp_val      = gms::math::simd_fast_exp_approx_2xf64(_mm_mul_pd(vrho_eq,cos2phic));
        __m128d  vres         = _mm_div_pd(exp_val,_mm_mul_pd(posPI,bessi0));
        _mm_storeu_pd(result+idx,vres);

        Tikhonov_phase_err_pdf_2xf64_t<Index+1ull,BatchSize>::compile_time_execute(rho_eq,phi_c,result);
    }
}

};

template<std::size_t Index,std::size_t BatchSize>
struct G_PLL_Bandwith_4xf32_t
{
static inline void
compile_time_execute(const float * __restrict__ Bl,
                     const float * __restrict__ Tb,
                     float       * __restrict__ result)
{
    if constexpr(Index < BatchSize)
    {
        constexpr std::size_t idx = Index*4ull;
        __m128 vBl  = _mm_loadu_ps(Bl+idx);
        __m128 vTb  = _mm_loadu_ps(Tb+idx);
        __m128 vres = _mm_rcp_ps(_mm_mul_ps(vBl,vTb));
        _mm_storeu_ps(result+idx,vres);

        G_PLL_Bandwith_4xf32_t<Index+1ull,BatchSize>::compile_time_execute(Bl,Tb,result); 
    }
}

};

template<std::size_t Index,std::size_t BatchSize,
         std::int32_t reciprocal_by_div>
struct G_PLL_Bandwith_2xf64_t
{
static inline void
compile_time_execute(const double * __restrict__ Bl,
                     const double * __restrict__ Tb,
                     double       * __restrict__ result)
{
    if constexpr(Index < BatchSize)
    {
        constexpr std::size_t idx = Index*2ull;
        __m128d vBl  = _mm_loadu_pd(Bl+idx);
        __m128d vTb  = _mm_loadu_pd(Tb+idx);
        if constexpr(reciprocal_by_div==1)
        {
            __m128d vres = _mm_div_pd(_mm_set1_pd(1.0),
                                      _mm_mul_pd(vBl,vTb));
        }
        else if constexpr(reciprocal_by_div==2)
        {
        __m128d vres = _mm_rcp14_pd(_mm_mul_pd(vBl,vTb));
        }
        _mm_storeu_pd(result+idx,vres);

        G_PLL_Bandwith_2xf64_t<Index+1ull,BatchSize,reciprocal_by_div>::compile_time_execute(Bl,Tb,result); 
    }
}

};

template<std::size_t Index,std::size_t BatchSize,
         std::int32_t copy_pointers_to_stack>
struct BPSK_Param_a_8_61_4xf32_t
{
static inline void
compile_time_execute(const float * __restrict__ Ac,
                     const float * __restrict__ Ts,
                     const float * __restrict__ M,
                     const float * __restrict__ N0,
                     const float * __restrict__ Bl,
                     const float * __restrict__ Tb,
                     float       * __restrict__ result)
{
    if constexpr(Index < BatchSize)
    {
        constexpr std::size_t idx = Index*4ull;
        if constexpr(copy_pointers_to_stack==1)
        {
            const float * __restrict__ pAc = Ac;
            const float * __restrict__ pTs = Ts;
            const float * __restrict__ pM  = M;
            const float * __restrict__ pN0 = N0;
            const float * __restrict__ pBl = Bl;
            const float * __restrict__ pTb = Tb;
            float       * __restrict__ presult = result;
            __m128 vAc   = _mm_loadu_ps(pAc+idx);
            __m128 vTs   = _mm_loadu_ps(pTs+idx);
            __m128 vM    = _mm_loadu_ps(pM+idx);
            __m128 vN0   = _mm_loadu_ps(pN0+idx);
            __m128 vBl   = _mm_loadu_ps(pBl+idx);
            __m128 vTb   = _mm_loadu_ps(pTb+idx);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            __m128 vAc   = _mm_loadu_ps(Ac+idx);
            __m128 vTs   = _mm_loadu_ps(Ts+idx);
            __m128 vM    = _mm_loadu_ps(M+idx);
            __m128 vN0   = _mm_loadu_ps(N0+idx);
            __m128 vBl   = _mm_loadu_ps(Bl+idx);
            __m128 vTb   = _mm_loadu_ps(Tb+idx);
        }
        __m128 Eb    = avg_bit_E_to_carrier_A_4xf32(vAc,vTs,vM);
        __m128 ratio = _mm_div_ps(Eb,_mm_add_ps(vN0,vN0));
        __m128 G     = G_PLL_bandwidth_4xf32(vBl,vTb);
        __m128 sqrtG = _mm_sub_ps(_mm_sqrt_ps(G),_mm_set1_ps(1.0f));
        __m128 vres  = _mm_sub_ps(ratio,_mm_mul_ps(sqrtG,sqrtG));
        _mm_storeu_ps(result+idx,vres);
        if constexpr(copy_pointers_to_stack==1)
        {
           BPSK_Param_a_8_61_4xf32_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(pAc,pTs,pM,pN0,pBl,pTb,presult);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            BPSK_Param_a_8_61_4xf32_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(Ac,Ts,M,N0,Bl,Tb,result);
        }
    }
}

};

template<std::size_t Index,std::size_t BatchSize,
         std::int32_t copy_pointers_to_stack>
struct BPSK_Param_a_8_61_2xf64_t
{
static inline void
compile_time_execute(const double * __restrict__ Ac,
                     const double * __restrict__ Ts,
                     const double * __restrict__ M,
                     const double * __restrict__ N0,
                     const double * __restrict__ Bl,
                     const double * __restrict__ Tb,
                     double       * __restrict__ result)
{
    if constexpr(Index < BatchSize)
    {
        constexpr std::size_t idx = Index*2ull;
        if constexpr(copy_pointers_to_stack==1)
        {
            const double * __restrict__ pAc = Ac;
            const double * __restrict__ pTs = Ts;
            const double * __restrict__ pM  = M;
            const double * __restrict__ pN0 = N0;
            const double * __restrict__ pBl = Bl;
            const double * __restrict__ pTb = Tb;
            double       * __restrict__ presult = result;
            __m128d vAc   = _mm_loadu_pd(pAc+idx);
            __m128d vTs   = _mm_loadu_pd(pTs+idx);
            __m128d vM    = _mm_loadu_pd(pM+idx);
            __m128d vN0   = _mm_loadu_pd(pN0+idx);
            __m128d vBl   = _mm_loadu_pd(pBl+idx);
            __m128d vTb   = _mm_loadu_pd(pTb+idx);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            __m128d vAc   = _mm_loadu_pd(Ac+idx);
            __m128d vTs   = _mm_loadu_pd(Ts+idx);
            __m128d vM    = _mm_loadu_pd(M+idx);
            __m128d vN0   = _mm_loadu_pd(N0+idx);
            __m128d vBl   = _mm_loadu_pd(Bl+idx);
            __m128d vTb   = _mm_loadu_pd(Tb+idx);
        }
        __m128d Eb    = avg_bit_E_to_carrier_A_2xf64(vAc,vTs,vM);
        __m128d ratio = _mm_div_pd(Eb,_mm_add_pd(vN0,vN0));
        __m128d G     = G_PLL_bandwidth_2xf64_v2(vBl,vTb);
        __m128d sqrtG = _mm_sub_pd(_mm_sqrt_pd(G),_mm_set1_pd(1.0f));
        __m128d vres  = _mm_sub_pd(ratio,_mm_mul_pd(sqrtG,sqrtG));
        
        if constexpr(copy_pointers_to_stack==1)
        {
           _mm_storeu_pd(presult+idx,vres);
           BPSK_Param_a_8_61_2xf64_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(pAc,pTs,pM,pN0,pBl,pTb,presult);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {   
            _mm_storeu_pd(result+idx,vres);
            BPSK_Param_a_8_61_2xf64_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(Ac,Ts,M,N0,Bl,Tb,result);
        }
    }
}

};

template<std::size_t Index,std::size_t BatchSize,
         std::int32_t copy_pointers_to_stack>
struct BPSK_Param_b_8_61_4xf32_t
{
static inline void
compile_time_execute(const float * __restrict__ Ac,
                     const float * __restrict__ Ts,
                     const float * __restrict__ M,
                     const float * __restrict__ N0,
                     const float * __restrict__ Bl,
                     const float * __restrict__ Tb,
                     float       * __restrict__ result)
{
    if constexpr(Index < BatchSize)
    {
        constexpr std::size_t idx = Index*4ull;
        if constexpr(copy_pointers_to_stack==1)
        {
            const float * __restrict__ pAc = Ac;
            const float * __restrict__ pTs = Ts;
            const float * __restrict__ pM  = M;
            const float * __restrict__ pN0 = N0;
            const float * __restrict__ pBl = Bl;
            const float * __restrict__ pTb = Tb;
            float       * __restrict__ presult = result;
            __m128 vAc   = _mm_loadu_ps(pAc+idx);
            __m128 vTs   = _mm_loadu_ps(pTs+idx);
            __m128 vM    = _mm_loadu_ps(pM+idx);
            __m128 vN0   = _mm_loadu_ps(pN0+idx);
            __m128 vBl   = _mm_loadu_ps(pBl+idx);
            __m128 vTb   = _mm_loadu_ps(pTb+idx);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            __m128 vAc   = _mm_loadu_ps(Ac+idx);
            __m128 vTs   = _mm_loadu_ps(Ts+idx);
            __m128 vM    = _mm_loadu_ps(M+idx);
            __m128 vN0   = _mm_loadu_ps(N0+idx);
            __m128 vBl   = _mm_loadu_ps(Bl+idx);
            __m128 vTb   = _mm_loadu_ps(Tb+idx);
        }
        __m128 Eb    = avg_bit_E_to_carrier_A_4xf32(vAc,vTs,vM);
        __m128 ratio = _mm_div_ps(Eb,_mm_add_ps(vN0,vN0));
        __m128 G     = G_PLL_bandwidth_4xf32(vBl,vTb);
        __m128 sqrtG = _mm_add_ps(_mm_sqrt_ps(G),_mm_set1_ps(1.0f));
        __m128 vres  = _mm_sub_ps(ratio,_mm_mul_ps(sqrtG,sqrtG));
        _mm_storeu_ps(result+idx,vres);
        if constexpr(copy_pointers_to_stack==1)
        {
           BPSK_Param_b_8_61_4xf32_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(pAc,pTs,pM,pN0,pBl,pTb,presult);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            BPSK_Param_b_8_61_4xf32_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(Ac,Ts,M,N0,Bl,Tb,result);
        }
    }
}

};

template<std::size_t Index,std::size_t BatchSize,
         std::int32_t copy_pointers_to_stack>
struct BPSK_Param_b_8_61_2xf64_t
{
static inline void
compile_time_execute(const double * __restrict__ Ac,
                     const double * __restrict__ Ts,
                     const double * __restrict__ M,
                     const double * __restrict__ N0,
                     const double * __restrict__ Bl,
                     const double * __restrict__ Tb,
                     double       * __restrict__ result)
{
    if constexpr(Index < BatchSize)
    {
        constexpr std::size_t idx = Index*2ull;
        if constexpr(copy_pointers_to_stack==1)
        {
            const double * __restrict__ pAc = Ac;
            const double * __restrict__ pTs = Ts;
            const double * __restrict__ pM  = M;
            const double * __restrict__ pN0 = N0;
            const double * __restrict__ pBl = Bl;
            const double * __restrict__ pTb = Tb;
            double       * __restrict__ presult = result;
            __m128d vAc   = _mm_loadu_pd(pAc+idx);
            __m128d vTs   = _mm_loadu_pd(pTs+idx);
            __m128d vM    = _mm_loadu_pd(pM+idx);
            __m128d vN0   = _mm_loadu_pd(pN0+idx);
            __m128d vBl   = _mm_loadu_pd(pBl+idx);
            __m128d vTb   = _mm_loadu_pd(pTb+idx);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            __m128d vAc   = _mm_loadu_pd(Ac+idx);
            __m128d vTs   = _mm_loadu_pd(Ts+idx);
            __m128d vM    = _mm_loadu_pd(M+idx);
            __m128d vN0   = _mm_loadu_pd(N0+idx);
            __m128d vBl   = _mm_loadu_pd(Bl+idx);
            __m128d vTb   = _mm_loadu_pd(Tb+idx);
        }
        __m128d Eb    = avg_bit_E_to_carrier_A_2xf64(vAc,vTs,vM);
        __m128d ratio = _mm_div_pd(Eb,_mm_add_pd(vN0,vN0));
        __m128d G     = G_PLL_bandwidth_2xf64_v2(vBl,vTb);
        __m128d sqrtG = _mm_add_pd(_mm_sqrt_pd(G),_mm_set1_pd(1.0));
        __m128d vres  = _mm_sub_pd(ratio,_mm_mul_pd(sqrtG,sqrtG));
        _mm_storeu_pd(result+idx,vres);
        if constexpr(copy_pointers_to_stack==1)
        {
           _mm_storeu_pd(presult+idx,vres);
           BPSK_Param_b_8_61_2xf64_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(pAc,pTs,pM,pN0,pBl,pTb,presult);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            _mm_storeu_pd(result+idx,vres);
            BPSK_Param_b_8_61_2xf64_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(Ac,Ts,M,N0,Bl,Tb,result);
        }
    }
}

};

template<std::size_t Index,std::size_t BatchSize,
         std::int32_t copy_pointers_to_stack>
struct OQPSK_Param_a2_8_64_4xf32_t
{
static inline void
compile_time_execute(const float * __restrict__ Ac,
                     const float * __restrict__ Ts,
                     const float * __restrict__ M,
                     const float * __restrict__ N0,
                     const float * __restrict__ Bl,
                     const float * __restrict__ Tb,
                     float       * __restrict__ result)
{
    if constexpr(Index < BatchSize)
    {
        constexpr std::size_t idx = Index*4ull;
        if constexpr(copy_pointers_to_stack==1)
        {
            const float * __restrict__ pAc = Ac;
            const float * __restrict__ pTs = Ts;
            const float * __restrict__ pM  = M;
            const float * __restrict__ pN0 = N0;
            const float * __restrict__ pBl = Bl;
            const float * __restrict__ pTb = Tb;
            float       * __restrict__ presult = result;
            __m128 vAc   = _mm_loadu_ps(pAc+idx);
            __m128 vTs   = _mm_loadu_ps(pTs+idx);
            __m128 vM    = _mm_loadu_ps(pM+idx);
            __m128 vN0   = _mm_loadu_ps(pN0+idx);
            __m128 vBl   = _mm_loadu_ps(pBl+idx);
            __m128 vTb   = _mm_loadu_ps(pTb+idx);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            __m128 vAc   = _mm_loadu_ps(Ac+idx);
            __m128 vTs   = _mm_loadu_ps(Ts+idx);
            __m128 vM    = _mm_loadu_ps(M+idx);
            __m128 vN0   = _mm_loadu_ps(N0+idx);
            __m128 vBl   = _mm_loadu_ps(Bl+idx);
            __m128 vTb   = _mm_loadu_ps(Tb+idx);
        }
        __m128 Eb     = avg_bit_E_to_carrier_A_4xf32(vAc,vTs,vM);
        __m128 G      = G_PLL_bandwidth_4xf32(vBl,vTb);
        __m128 sqrt2G = _mm_sub_ps(_mm_set1_ps(1.0f),_mm_sqrt_ps(_mm_add_ps(G,G)));
        __m128 factor = _mm_add_ps(G,sqrt2G);
        __m128 vres   = _mm_mul_ps(_mm_div_ps(Eb,N0),factor);
        _mm_storeu_ps(result+idx,vres);
        if constexpr(copy_pointers_to_stack==1)
        {
           OQPSK_Param_a2_8_64_4xf32_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(pAc,pTs,pM,pN0,pBl,pTb,presult);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            OQPSK_Param_a2_8_64_4xf32_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(Ac,Ts,M,N0,Bl,Tb,result);
        }
    }
}

};

template<std::size_t Index,std::size_t BatchSize,
         std::int32_t copy_pointers_to_stack>
struct OQPSK_Param_a2_8_64_2xf64_t
{
static inline void
compile_time_execute(const double * __restrict__ Ac,
                     const double * __restrict__ Ts,
                     const double * __restrict__ M,
                     const double * __restrict__ N0,
                     const double * __restrict__ Bl,
                     const double * __restrict__ Tb,
                     double       * __restrict__ result)
{
    if constexpr(Index < BatchSize)
    {
        constexpr std::size_t idx = Index*2ull;
        if constexpr(copy_pointers_to_stack==1)
        {
            const double * __restrict__ pAc = Ac;
            const double * __restrict__ pTs = Ts;
            const double * __restrict__ pM  = M;
            const double * __restrict__ pN0 = N0;
            const double * __restrict__ pBl = Bl;
            const double * __restrict__ pTb = Tb;
            double       * __restrict__ presult = result;
            __m128d vAc   = _mm_loadu_pd(pAc+idx);
            __m128d vTs   = _mm_loadu_pd(pTs+idx);
            __m128d vM    = _mm_loadu_pd(pM+idx);
            __m128d vN0   = _mm_loadu_pd(pN0+idx);
            __m128d vBl   = _mm_loadu_pd(pBl+idx);
            __m128d vTb   = _mm_loadu_pd(pTb+idx);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            __m128d vAc   = _mm_loadu_pd(Ac+idx);
            __m128d vTs   = _mm_loadu_pd(Ts+idx);
            __m128d vM    = _mm_loadu_pd(M+idx);
            __m128d vN0   = _mm_loadu_pd(N0+idx);
            __m128d vBl   = _mm_loadu_pd(Bl+idx);
            __m128d vTb   = _mm_loadu_pd(Tb+idx);
        }
        __m128d Eb     = avg_bit_E_to_carrier_A_2xf64(vAc,vTs,vM);
        __m128d G      = G_PLL_bandwidth_2xf64_v2(vBl,vTb);
        __m128d sqrt2G = _mm_sub_pd(_mm_set1_pd(1.0),_mm_sqrt_pd(_mm_add_pd(G,G)));
        __m128d factor = _mm_add_pd(G,sqrt2G);
        __m128d vres   = _mm_mul_pd(_mm_div_pd(Eb,N0),factor);
        _mm_storeu_pd(result+idx,vres);
        if constexpr(copy_pointers_to_stack==1)
        {
           OQPSK_Param_a2_8_64_2xf64_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(pAc,pTs,pM,pN0,pBl,pTb,presult);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            OQPSK_Param_a2_8_64_2xf64_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(Ac,Ts,M,N0,Bl,Tb,result);
        }
    }
}

};

template<std::size_t Index,std::size_t BatchSize,
         std::int32_t copy_pointers_to_stack>
struct MSK_Param_a2_8_65_4xf32_t
{
static inline void
compile_time_execute(const float * __restrict__ Ac,
                     const float * __restrict__ Ts,
                     const float * __restrict__ M,
                     const float * __restrict__ N0,
                     const float * __restrict__ Bl,
                     const float * __restrict__ Tb,
                     float       * __restrict__ result)
{
    if constexpr(Index < BatchSize)
    {
        constexpr std::size_t idx = Index*4ull;
        if constexpr(copy_pointers_to_stack==1)
        {
            const float * __restrict__ pAc = Ac;
            const float * __restrict__ pTs = Ts;
            const float * __restrict__ pM  = M;
            const float * __restrict__ pN0 = N0;
            const float * __restrict__ pBl = Bl;
            const float * __restrict__ pTb = Tb;
            float       * __restrict__ presult = result;
            __m128 vAc   = _mm_loadu_ps(pAc+idx);
            __m128 vTs   = _mm_loadu_ps(pTs+idx);
            __m128 vM    = _mm_loadu_ps(pM+idx);
            __m128 vN0   = _mm_loadu_ps(pN0+idx);
            __m128 vBl   = _mm_loadu_ps(pBl+idx);
            __m128 vTb   = _mm_loadu_ps(pTb+idx);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            __m128 vAc   = _mm_loadu_ps(Ac+idx);
            __m128 vTs   = _mm_loadu_ps(Ts+idx);
            __m128 vM    = _mm_loadu_ps(M+idx);
            __m128 vN0   = _mm_loadu_ps(N0+idx);
            __m128 vBl   = _mm_loadu_ps(Bl+idx);
            __m128 vTb   = _mm_loadu_ps(Tb+idx);
        }
        __m128 Eb     = avg_bit_E_to_carrier_A_4xf32(vAc,vTs,vM);
        __m128 G      = G_PLL_bandwidth_4xf32(vBl,vTb);
        __m128 sqrt2G = _mm_sub_ps(_mm_set1_ps(0.7026423672846755428877589264195f),_mm_sqrt_ps(_mm_add_ps(G,G)));
        __m128 factor = _mm_add_ps(G,sqrt2G);
        __m128 vres   = _mm_mul_ps(_mm_div_ps(Eb,N0),factor);
        _mm_storeu_ps(result+idx,vres);
        if constexpr(copy_pointers_to_stack==1)
        {
           MSK_Param_a2_8_65_4xf32_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(pAc,pTs,pM,pN0,pBl,pTb,presult);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            MSK_Param_a2_8_65_4xf32_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(Ac,Ts,M,N0,Bl,Tb,result);
        }
    }
}

};

template<std::size_t Index,std::size_t BatchSize,
         std::int32_t copy_pointers_to_stack>
struct MSK_Param_a2_8_65_2xf64_t
{
static inline void
compile_time_execute(const double * __restrict__ Ac,
                     const double * __restrict__ Ts,
                     const double * __restrict__ M,
                     const double * __restrict__ N0,
                     const double * __restrict__ Bl,
                     const double * __restrict__ Tb,
                     double       * __restrict__ result)
{
    if constexpr(Index < BatchSize)
    {
        constexpr std::size_t idx = Index*2ull;
        if constexpr(copy_pointers_to_stack==1)
        {
            const double * __restrict__ pAc = Ac;
            const double * __restrict__ pTs = Ts;
            const double * __restrict__ pM  = M;
            const double * __restrict__ pN0 = N0;
            const double * __restrict__ pBl = Bl;
            const double * __restrict__ pTb = Tb;
            double       * __restrict__ presult = result;
            __m128d vAc   = _mm_loadu_pd(pAc+idx);
            __m128d vTs   = _mm_loadu_pd(pTs+idx);
            __m128d vM    = _mm_loadu_pd(pM+idx);
            __m128d vN0   = _mm_loadu_pd(pN0+idx);
            __m128d vBl   = _mm_loadu_pd(pBl+idx);
            __m128d vTb   = _mm_loadu_pd(pTb+idx);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            __m128d vAc   = _mm_loadu_pd(Ac+idx);
            __m128d vTs   = _mm_loadu_pd(Ts+idx);
            __m128d vM    = _mm_loadu_pd(M+idx);
            __m128d vN0   = _mm_loadu_pd(N0+idx);
            __m128d vBl   = _mm_loadu_pd(Bl+idx);
            __m128d vTb   = _mm_loadu_pd(Tb+idx);
        }
        __m128d Eb     = avg_bit_E_to_carrier_A_2xf64(vAc,vTs,vM);
        __m128d G      = G_PLL_bandwidth_2xf64_v2(vBl,vTb);
        __m128d sqrt2G = _mm_sub_pd(_mm_set1_pd(0.7026423672846755428877589264195),_mm_sqrt_pd(_mm_add_pd(G,G)));
        __m128d factor = _mm_add_pd(G,sqrt2G);
        __m128d vres   = _mm_mul_pd(_mm_div_pd(Eb,N0),factor);
       
        if constexpr(copy_pointers_to_stack==1)
        {
            _mm_storeu_pd(presult+idx,vres);
           MSK_Param_a2_8_65_2xf64_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(pAc,pTs,pM,pN0,pBl,pTb,presult);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {   
             _mm_storeu_pd(result+idx,vres);
            MSK_Param_a2_8_65_2xf64_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(Ac,Ts,M,N0,Bl,Tb,result);
        }
    }
}

};

template<std::size_t Index,std::size_t BatchSize,
         std::int32_t copy_pointers_to_stack>
struct MSK_Param_b2_8_65_4xf32_t
{
static inline void
compile_time_execute(const float * __restrict__ Ac,
                     const float * __restrict__ Ts,
                     const float * __restrict__ M,
                     const float * __restrict__ N0,
                     const float * __restrict__ Bl,
                     const float * __restrict__ Tb,
                     float       * __restrict__ result)
{
    if constexpr(Index < BatchSize)
    {
        constexpr std::size_t idx = Index*4ull;
        if constexpr(copy_pointers_to_stack==1)
        {
            const float * __restrict__ pAc = Ac;
            const float * __restrict__ pTs = Ts;
            const float * __restrict__ pM  = M;
            const float * __restrict__ pN0 = N0;
            const float * __restrict__ pBl = Bl;
            const float * __restrict__ pTb = Tb;
            float       * __restrict__ presult = result;
            __m128 vAc   = _mm_loadu_ps(pAc+idx);
            __m128 vTs   = _mm_loadu_ps(pTs+idx);
            __m128 vM    = _mm_loadu_ps(pM+idx);
            __m128 vN0   = _mm_loadu_ps(pN0+idx);
            __m128 vBl   = _mm_loadu_ps(pBl+idx);
            __m128 vTb   = _mm_loadu_ps(pTb+idx);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            __m128 vAc   = _mm_loadu_ps(Ac+idx);
            __m128 vTs   = _mm_loadu_ps(Ts+idx);
            __m128 vM    = _mm_loadu_ps(M+idx);
            __m128 vN0   = _mm_loadu_ps(N0+idx);
            __m128 vBl   = _mm_loadu_ps(Bl+idx);
            __m128 vTb   = _mm_loadu_ps(Tb+idx);
        }
        __m128 Eb     = avg_bit_E_to_carrier_A_4xf32(vAc,vTs,vM);
        __m128 G      = G_PLL_bandwidth_4xf32(vBl,vTb);
        __m128 sqrt2G = _mm_add_ps(_mm_set1_ps(0.7026423672846755428877589264195f),_mm_sqrt_ps(_mm_add_ps(G,G)));
        __m128 factor = _mm_add_ps(G,sqrt2G);
        __m128 vres   = _mm_mul_ps(_mm_div_ps(Eb,N0),factor);
        _mm_storeu_ps(result+idx,vres);
        if constexpr(copy_pointers_to_stack==1)
        {
           MSK_Param_b2_8_65_4xf32_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(pAc,pTs,pM,pN0,pBl,pTb,presult);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            MSK_Param_b2_8_65_4xf32_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(Ac,Ts,M,N0,Bl,Tb,result);
        }
    }
}

};

template<std::size_t Index,std::size_t BatchSize,
         std::int32_t copy_pointers_to_stack>
struct MSK_Param_b2_8_65_2xf64_t
{
static inline void
compile_time_execute(const double * __restrict__ Ac,
                     const double * __restrict__ Ts,
                     const double * __restrict__ M,
                     const double * __restrict__ N0,
                     const double * __restrict__ Bl,
                     const double * __restrict__ Tb,
                     double       * __restrict__ result)
{
    if constexpr(Index < BatchSize)
    {
        constexpr std::size_t idx = Index*2ull;
        if constexpr(copy_pointers_to_stack==1)
        {
            const double * __restrict__ pAc = Ac;
            const double * __restrict__ pTs = Ts;
            const double * __restrict__ pM  = M;
            const double * __restrict__ pN0 = N0;
            const double * __restrict__ pBl = Bl;
            const double * __restrict__ pTb = Tb;
            double       * __restrict__ presult = result;
            __m128d vAc   = _mm_loadu_pd(pAc+idx);
            __m128d vTs   = _mm_loadu_pd(pTs+idx);
            __m128d vM    = _mm_loadu_pd(pM+idx);
            __m128d vN0   = _mm_loadu_pd(pN0+idx);
            __m128d vBl   = _mm_loadu_pd(pBl+idx);
            __m128d vTb   = _mm_loadu_pd(pTb+idx);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            __m128d vAc   = _mm_loadu_pd(Ac+idx);
            __m128d vTs   = _mm_loadu_pd(Ts+idx);
            __m128d vM    = _mm_loadu_pd(M+idx);
            __m128d vN0   = _mm_loadu_pd(N0+idx);
            __m128d vBl   = _mm_loadu_pd(Bl+idx);
            __m128d vTb   = _mm_loadu_pd(Tb+idx);
        }
        __m128d Eb     = avg_bit_E_to_carrier_A_2xf64(vAc,vTs,vM);
        __m128d G      = G_PLL_bandwidth_2xf64_v2(vBl,vTb);
        __m128d sqrt2G = _mm_add_pd(_mm_set1_pd(0.7026423672846755428877589264195f),_mm_sqrt_pd(_mm_add_pd(G,G)));
        __m128d factor = _mm_add_pd(G,sqrt2G);
        __m128d vres   = _mm_mul_pd(_mm_div_pd(Eb,N0),factor);
        if constexpr(copy_pointers_to_stack==1)
        {  
            _mm_storeu_pd(presult+idx,vres);
           MSK_Param_b2_8_65_2xf64_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(pAc,pTs,pM,pN0,pBl,pTb,presult);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {   
            _mm_storeu_pd(result+idx,vres);
            MSK_Param_b2_8_65_2xf64_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(Ac,Ts,M,N0,Bl,Tb,result);
        }
    }
}

};

template<std::size_t Index,std::size_t BatchSize,
         std::int32_t copy_pointers_to_stack>
struct Analytic_BEP_DPSK_chernoff_bound_8_9_4xf32_t
{
static inline void 
compile_time_execute(const float * __restrict__ Ac,
                     const float * __restrict__ Ts,
                     const float * __restrict__ M,
                     const float * __restrict__ N0,
                     float       * __restrict__ result) 
{
    if constexpr(Index < BatchSize)
    {
        constexpr std::size_t idx = Index*4ull;
        if constexpr(copy_pointers_to_stack==1)
        {
            const float * __restrict__ pAc = Ac;
            const float * __restrict__ pTs = Ts;
            const float * __restrict__ pM  = M;
            const float * __restrict__ pN0 = N0;
            float       * __restrict__ presult = result;
            __m128 vAc   = _mm_loadu_ps(pAc+idx);
            __m128 vTs   = _mm_loadu_ps(pTs+idx);
            __m128 vM    = _mm_loadu_ps(pM+idx);
            __m128 vN0   = _mm_loadu_ps(pN0+idx);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            __m128 vAc   = _mm_loadu_ps(Ac+idx);
            __m128 vTs   = _mm_loadu_ps(Ts+idx);
            __m128 vM    = _mm_loadu_ps(M+idx);
            __m128 vN0   = _mm_loadu_ps(N0+idx);
        }
        __m128 one      = _mm_set1_ps(1.0f);
        __m128 Pi_div_M = _mm_div_ps(_mm_set1_ps(3.1415926535897932384626433832795f),vM);
        __m128 Es       = avg_symbol_E_to_carrier_A_4xf32(vAc,vTs,vM);
        __m128 snr      = _mm_div_ps(_mm_add_ps(Es,Es),vN0);
        __m128 cosPiM   = gms::math::_mm_ceph_cosf_ps(Pi_div_M);
        __m128 exp_arg  = _mm_sqrt_ps(_mm_mul_ps(snr,_mm_sub_ps(one,cosPiM)));
        __m128 cos_ratio= _mm_div_ps(_mm_add_ps(one,cosPiM),_mm_add_ps(cosPiM,cosPiM));
        __m128 left_term= _mm_mul_ps(_mm_set1_ps(1.03f),_mm_sqrt_ps(cos_ratio));  
        __m128 vres     = _mm_mul_ps(left_term,
                              gms::math::simd_fast_exp_approx_4xf32(_mm_sub_ps(_mm_setzero_ps(),exp_arg)));
        
        if constexpr(copy_pointers_to_stack==1)
        {
            _mm_storeu_ps(presult+idx,vres);
            Analytic_BEP_DPSK_chernoff_bound_8_9_4xf32_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(pAc,pTs,pM,pN0,presult);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            _mm_storeu_ps(result+idx,vres);
            Analytic_BEP_DPSK_chernoff_bound_8_9_4xf32_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(Ac,Ts,M,N0,result);
        }
    }
}

};

template<std::size_t Index,std::size_t BatchSize,
         std::int32_t copy_pointers_to_stack>
struct Analytic_BEP_DPSK_chernoff_bound_8_9_2xf64_t
{
static inline void 
compile_time_execute(const double * __restrict__ Ac,
                     const double * __restrict__ Ts,
                     const double * __restrict__ M,
                     const double * __restrict__ N0,
                     double       * __restrict__ result) 
{
    if constexpr(Index < BatchSize)
    {
        constexpr std::size_t idx = Index*2ull;
        if constexpr(copy_pointers_to_stack==1)
        {
            const double * __restrict__ pAc = Ac;
            const double * __restrict__ pTs = Ts;
            const double * __restrict__ pM  = M;
            const double * __restrict__ pN0 = N0;
            double       * __restrict__ presult = result;
            __m128d vAc   = _mm_loadu_pd(pAc+idx);
            __m128d vTs   = _mm_loadu_pd(pTs+idx);
            __m128d vM    = _mm_loadu_pd(pM+idx);
            __m128d vN0   = _mm_loadu_pd(pN0+idx);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            __m128d vAc   = _mm_loadu_pd(Ac+idx);
            __m128d vTs   = _mm_loadu_pd(Ts+idx);
            __m128d vM    = _mm_loadu_pd(M+idx);
            __m128d vN0   = _mm_loadu_pd(N0+idx);
        }
        __m128d one      = _mm_set1_pd(1.0);
        __m128d Pi_div_M = _mm_div_pd(_mm_set1_pd(3.1415926535897932384626433832795),vM);
        __m128d Es       = avg_symbol_E_to_carrier_A_2xf64(vAc,vTs,vM);
        __m128d snr      = _mm_div_pd(_mm_add_pd(Es,Es),vN0);
        __m128d cosPiM   = _mm_cos_pd(Pi_div_M);
        __m128d exp_arg  = _mm_sqrt_pd(_mm_mul_pd(snr,_mm_sub_pd(one,cosPiM)));
        __m128d cos_ratio= _mm_div_pd(_mm_add_pd(one,cosPiM),_mm_add_pd(cosPiM,cosPiM));
        __m128d left_term= _mm_mul_pd(_mm_set1_pd(1.03),_mm_sqrt_pd(cos_ratio));  
        __m128d vres     = _mm_mul_pd(left_term,
                              gms::math::simd_fast_exp_approx_2xf64(_mm_sub_pd(_mm_setzero_pd(),exp_arg)));
        
        if constexpr(copy_pointers_to_stack==1)
        {
            _mm_storeu_pd(presult+idx,vres);
            Analytic_BEP_DPSK_chernoff_bound_8_9_2xf64_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(pAc,pTs,pM,pN0,presult);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            _mm_storeu_pd(result+idx,vres);
            Analytic_BEP_DPSK_chernoff_bound_8_9_2xf64_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(Ac,Ts,M,N0,result);
        }
    }
}

};

template<std::size_t Index,std::size_t BatchSize,
         std::int32_t copy_pointers_to_stack>
struct Analytic_SEP_QAM4_8_107_4xf32_t
{
static inline void 
compile_time_execute(const float * __restrict__ gamma,
                     float *       __restrict__ result)
{
    if constexpr(Index < BatchSize)
    {
        constexpr std::size_t idx = Index+4ull;
        if constexpr(copy_pointers_to_stack==1)
        {
            const float * __restrict__ pgamma  = gamma;
            float       * __restrict__ presult = result;
            __m128 vgamma = _mm_loadu_ps(pgamma+idx);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            __m128 vgamma = _mm_loadu_ps(gamma+idx);
        }
        __m128 one           = _mm_set1_ps(1.0f);
        __m128 quarter       = _mm_set1_ps(0.25f);
        __m128 one_add_gamma = _mm_add_ps(one,vgamma); 
        __m128 sqrt_arg1     = _mm_div_ps(vgamma,one_add_gamma);
        __m128 sqrt_val1     = _mm_sqrt_ps(sqrt_arg1);
        __m128 sqrt_arg2     = _mm_div_ps(one_add_gamma,vgamma);
        __m128 sqrt_val2     = _mm_sqrt_ps(sqrt_arg2);
        __m128 term1         = _mm_sub_ps(one,sqrt_val1);
        __m128 atan_val      = _mm_mul_ps(_mm_set1_ps(1.2732395447351626861510701069801f),sqrt_val2);
        __m128 term2         = _mm_mul_ps(quarter,_mm_mul_ps(term1,atan_val));
        __m128 vres          = _mm_sub_ps(term1,term2);

        if constexpr(copy_pointers_to_stack==1)
        {
            _mm_storeu_ps(presult+idx,vres);
            Analytic_SEP_QAM4_8_107_4xf32_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(pgamma,presult);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            _mm_storeu_ps(result+idx,vres);
            Analytic_SEP_QAM4_8_107_4xf32_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(gamma,result);
        }
    }
}

};


template<std::size_t Index,std::size_t BatchSize,
         std::int32_t copy_pointers_to_stack>
struct Analytic_SEP_QAM4_8_107_2xf64_t
{
static inline void 
compile_time_execute(const double * __restrict__ gamma,
                     double *       __restrict__ result)
{
    if constexpr(Index < BatchSize)
    {
        constexpr std::size_t idx = Index+2ull;
        if constexpr(copy_pointers_to_stack==1)
        {
            const double * __restrict__ pgamma  = gamma;
            double       * __restrict__ presult = result;
            __m128d vgamma = _mm_loadu_pd(pgamma+idx);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            __m128d vgamma = _mm_loadu_pd(gamma+idx);
        }
        __m128d one           = _mm_set1_pd(1.0);
        __m128d quarter       = _mm_set1_pd(0.25);
        __m128d one_add_gamma = _mm_add_pd(one,vgamma); 
        __m128d sqrt_arg1     = _mm_div_pd(vgamma,one_add_gamma);
        __m128d sqrt_val1     = _mm_sqrt_pd(sqrt_arg1);
        __m128d sqrt_arg2     = _mm_div_pd(one_add_gamma,vgamma);
        __m128d sqrt_val2     = _mm_sqrt_pd(sqrt_arg2);
        __m128d term1         = _mm_sub_pd(one,sqrt_val1);
        __m128d atan_val      = _mm_mul_pd(_mm_set1_pd(1.2732395447351626861510701069801),sqrt_val2);
        __m128d term2         = _mm_mul_pd(quarter,_mm_mul_pd(term1,atan_val));
        __m128d vres          = _mm_sub_pd(term1,term2);

        if constexpr(copy_pointers_to_stack==1)
        {
            _mm_storeu_pd(presult+idx,vres);
            Analytic_SEP_QAM4_8_107_2xf64_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(pgamma,presult);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            _mm_storeu_ps(result+idx,vres);
            Analytic_SEP_QAM4_8_107_2xf64_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(gamma,result);
        }
    }
}

};

template<std::size_t Index,std::size_t BatchSize,
         std::int32_t copy_pointers_to_stack>
struct Analytic_SEP_MPSK_rayleigh_chan_8_112_4xf32_t
{
static inline void 
compile_time_execute(const float * __restrict__ gamma,
                     const float * __restrict__ M,
                     float       * __restrict__ result)
{
    if constexpr(Index < BatchSize)
    {
        constexpr std::size_t idx = Index*4ull;
        if constexpr(copy_pointers_to_stack==1)
        {
            const float  * __restrict__ pgamma  = gamma;
            const float  * __restrict__ pM      = M;
            float        * __restrict__ presult = result;
            __m128 vgamma = _mm_loadu_ps(pgamma+idx);
            __m128 vM     = _mm_loadu_ps(pM+idx);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            __m128 vgamma = _mm_loadu_ps(gamma+idx);
            __m128 vM     = _mm_loadu_ps(M+idx);
        }
        __m128 Pi       = _mm_set1_ps(3.1415926535897932384626433832795f);
        __m128 Pi_div_M = _mm_div_ps(Pi,vM);
        __m128 sin_tmp  = gms::math::_mm_ceph_sinf_ps(Pi_div_M);
        __m128 halfPi   = _mm_set1_ps(1.5707963267948966192313216916398f);
        __m128 one      = _mm_set1_ps(1.0f);
        __m128 g_psk    = _mm_mul_ps(sin_tmp,sin_tmp);
        __m128 M_sub_1  = _mm_sub_ps(vM,one);
        __m128 term1    = _mm_div_ps(M_sub_1,vM);
        __m128 gpskgamm = _mm_mul_ps(g_psk,vgamma);
        __m128 term2    = _mm_div_ps(gpskgamm,_mm_div_ps(one,gpskgamm));
        __m128 sqrt_val = _mm_sqrt_ps(term2);
        __m128 tmp      = _mm_mul_ps(M_sub_1,Pi);
        __m128 term3    = _mm_div_ps(vM,term3);
        __m128 cot_val  = _mm_tan_ps(_mm_mul_ps(halfPi,Pi_div_M));
        __m128 term4    = _mm_sub_ps(one,_mm_mul_ps(sqrt_val,term3));
        __m128 atan_val = _mm_mul_ps(halfPi,_mm_atan_ps(_mm_mul_ps(sqrt_val,cot_val)));   
        __m128 vres     = _mm_mul_ps(term1,_mm_mul_ps(term4,atan_val));

        if constexpr(copy_pointers_to_stack==1)
        {
            _mm_storeu_ps(presult+idx,vres);
            Analytic_SEP_MPSK_rayleigh_chan_8_112_4xf32_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(pgamma,M,presult);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            _mm_storeu_ps(result+idx,vres);
            Analytic_SEP_MPSK_rayleigh_chan_8_112_4xf32_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(gamma,M,result);
        }
    }
}

};

template<std::size_t Index,std::size_t BatchSize,
         std::int32_t copy_pointers_to_stack>
struct Analytic_SEP_MPSK_rayleigh_chan_8_112_2xf64_t
{
static inline void 
compile_time_execute(const double * __restrict__ gamma,
                     const double * __restrict__ M,
                     double       * __restrict__ result)
{
    if constexpr(Index < BatchSize)
    {
        constexpr std::size_t idx = Index*2ull;
        if constexpr(copy_pointers_to_stack==1)
        {
            const double  * __restrict__ pgamma  = gamma;
            const double  * __restrict__ pM      = M;
            double        * __restrict__ presult = result;
            __m128d vgamma = _mm_loadu_pd(pgamma+idx);
            __m128d vM     = _mm_loadu_pd(pM+idx);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            __m128d vgamma = _mm_loadu_pd(gamma+idx);
            __m128d vM     = _mm_loadu_pd(M+idx);
        }
        __m128d Pi       = _mm_set1_pd(3.1415926535897932384626433832795);
        __m128d Pi_div_M = _mm_div_pd(Pi,vM);
        __m128d sin_tmp  = _mm_sin_pd(Pi_div_M);
        __m128d halfPi   = _mm_set1_pd(1.5707963267948966192313216916398);
        __m128d one      = _mm_set1_pd(1.0f);
        __m128d g_psk    = _mm_mul_pd(sin_tmp,sin_tmp);
        __m128d M_sub_1  = _mm_sub_pd(vM,one);
        __m128d term1    = _mm_div_pd(M_sub_1,vM);
        __m128d gpskgamm = _mm_mul_pd(g_psk,vgamma);
        __m128d term2    = _mm_div_pd(gpskgamm,_mm_div_pd(one,gpskgamm));
        __m128d sqrt_val = _mm_sqrt_pd(term2);
        __m128d tmp      = _mm_mul_pd(M_sub_1,Pi);
        __m128d term3    = _mm_div_pd(vM,term3);
        __m128d cot_val  = _mm_tan_pd(_mm_mul_pd(halfPi,Pi_div_M));
        __m128d term4    = _mm_sub_pd(one,_mm_mul_pd(sqrt_val,term3));
        __m128d atan_val = _mm_mul_pd(halfPi,_mm_atan_pd(_mm_mul_pd(sqrt_val,cot_val)));   
        __m128d vres     = _mm_mul_pd(term1,_mm_mul_pd(term4,atan_val));

        if constexpr(copy_pointers_to_stack==1)
        {
            _mm_storeu_pd(presult+idx,vres);
            Analytic_SEP_MPSK_rayleigh_chan_8_112_2xf64_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(pgamma,M,presult);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            _mm_storeu_pd(result+idx,vres);
            Analytic_SEP_MPSK_rayleigh_chan_8_112_2xf64_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(gamma,M,result);
        }
    }
}

};

template<std::size_t Index,std::size_t BatchSize,
         std::int32_t copy_pointers_to_stack>
struct Riccian_factor_8_144_4xf32_t
{
static inline void 
compile_time_execute(const float * __restrict__ m_I,
                     const float * __restrict__ m_Q,
                     const float * __restrict__ sigma,
                     float       * __restrict__ result)
{
    if constexpr (Index < BatchSize)
    {
        constexpr std::size_t idx = Index*4ull;
        if constexpr(copy_pointers_to_stack==1)
        {
            const float  * __restrict__ pm_I    = m_I;
            const float  * __restrict__ pm_Q    = M_Q;
            const float  * __restrict__ psigma  = sigma;
            float        * __restrict__ presult = result;
            __m128 vm_I   = _mm_loadu_ps(pm_I+idx);
            __m128 vm_Q   = _mm_loadu_ps(pm_Q+idx);
            __m128 vsigma = _mm_loadu_ps(psigma+idx);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            __m128 vm_I   = _mm_loadu_ps(m_I+idx);
            __m128 vm_Q   = _mm_loadu_ps(m_Q+idx);
            __m128 vsigma = _mm_loadu_ps(sigma+idx);
        }

        __m128 m_I_sqr   = _mm_mul_ps(vm_I,vm_I);
        __m128 sigma_sqr = _mm_mul_ps(vsigma,vsigma);
        __m128 m_Q_sqr   = _mm_mul_ps(vm_Q,vm_Q);
        __m128 vres      = _mm_div_ps(_mm_add_ps(m_I_sqr,m_Q_sqr),_mm_add_ps(sigma_sqr,sigma_sqr));

        if constexpr(copy_pointers_to_stack==1)
        {
            _mm_storeu_ps(presult+idx,vres);
            Riccian_factor_8_144_4xf32_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(pm_I,pm_Q,psigma,presult);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            _mm_storeu_ps(result+idx,vres);
            Riccian_factor_8_144_4xf32_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(m_I,m_Q,sigma,result);
        }
    }
}

};

template<std::size_t Index,std::size_t BatchSize,
         std::int32_t copy_pointers_to_stack>
struct Riccian_factor_8_144_2xf64_t
{
static inline void 
compile_time_execute(const double * __restrict__ m_I,
                     const double * __restrict__ m_Q,
                     const double * __restrict__ sigma,
                     double       * __restrict__ result)
{
    if constexpr (Index < BatchSize)
    {
        constexpr std::size_t idx = Index*2ull;
        if constexpr(copy_pointers_to_stack==1)
        {
            const double  * __restrict__ pm_I    = m_I;
            const double  * __restrict__ pm_Q    = M_Q;
            const double  * __restrict__ psigma  = sigma;
            double        * __restrict__ presult = result;
            __m128d vm_I   = _mm_loadu_pd(pm_I+idx);
            __m128d vm_Q   = _mm_loadu_pd(pm_Q+idx);
            __m128d vsigma = _mm_loadu_pd(psigma+idx);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            __m128d vm_I   = _mm_loadu_pd(m_I+idx);
            __m128d vm_Q   = _mm_loadu_pd(m_Q+idx);
            __m128d vsigma = _mm_loadu_pd(sigma+idx);
        }

        __m128d m_I_sqr   = _mm_mul_pd(vm_I,vm_I);
        __m128d sigma_sqr = _mm_mul_pd(vsigma,vsigma);
        __m128d m_Q_sqr   = _mm_mul_pd(vm_Q,vm_Q);
        __m128d vres      = _mm_div_pd(_mm_add_pd(m_I_sqr,m_Q_sqr),_mm_add_pd(sigma_sqr,sigma_sqr));

        if constexpr(copy_pointers_to_stack==1)
        {
            _mm_storeu_pd(presult+idx,vres);
            Riccian_factor_8_144_2xf64_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(pm_I,pm_Q,psigma,presult);
        }
        else if constexpr(copy_pointers_to_stack==2)
        {
            _mm_storeu_ps(result+idx,vres);
            Riccian_factor_8_144_2xf64_t<Index+1ull,BatchSize,copy_pointers_to_stack>::compile_time_execute(m_I,m_Q,sigma,result);
        }
    }
}

};


}

}

#endif /*__GMS_METAPROG_BEP_SEP_CH8_H__*/