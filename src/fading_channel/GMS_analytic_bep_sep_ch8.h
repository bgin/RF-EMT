
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
#define __GMS_ANALYTIC_BEP_SEP_CH8_H__ 160920261158

#include <cstdint>
#include <complex>
#include "GMS_config.h"
#include "GMS_gaussian_Q_approx.h"
#include "GMS_marcum_Q_approx.h"
#include "GMS_cephes_double.h"


namespace file_info 
{

     static const unsigned int GMS_ANALYTIC_BEP_SEP_CH8_MAJOR = 1;
     static const unsigned int GMS_ANALYTIC_BEP_SEP_CH8_MINOR = 1;
     static const unsigned int GMS_ANALYTIC_BEP_SEP_CH8_MICRO = 0;
     static const unsigned int GMS_ANALYTIC_BEP_SEP_CH8_FULLVER =
       1000U*GMS_ANALYTIC_BEP_SEP_CH8_MAJOR+100U*GMS_ANALYTIC_BEP_SEP_CH8_MINOR+
       10U*GMS_ANALYTIC_BEP_SEP_CH8_MICRO;
     static const char GMS_ANALYTIC_BEP_SEP_CH8_CREATION_DATE[] = "18-09-2026 06:19AM +00200 (FRI 18 AUG 2026 GMT+2)";
     static const char GMS_ANALYTIC_BEP_SEP_CH8_BUILD_DATE[]    = __DATE__; 
     static const char GMS_ANALYTIC_BEP_SEP_CH8_BUILD_TIME[]    = __TIME__;
     static const char GMS_ANALYTIC_BEP_SEP_CH8_SYNOPSIS[]      = "Analytic Formulae of the Bit and Symbol Error Probabilities based on the approximated Gaussian-Q function.\
	                                                          Based on the chapter 8 of M.K Simon, M.S. Alouini: Digital Communication over Fading Channels 1st ed\
															  ISBN-13 978-0471317791";

}

#if !defined(ANALYTIC_BEP_SEP_CH8_USE_CEPHES_DOUBLE)
#define ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE 1
#endif 

#if !defined(ANALYTIC_BEP_SEP_CH8_DO_WARMUP_CALL)
#define ANALYTIC_BEP_SEP_CH8_DO_WARMUP_CALL 1
#endif 

#if !defined(ANALYTIC_BEP_SEP_CH8_SPLIT_MATH_IMPL_PERF_TEST)
#define ANALYTIC_BEP_SEP_CH8_SPLIT_MATH_IMPL_PERF_TEST 1
#endif 

#if !defined(ANALYTIC_BEP_SEP_CH8_USE_PARAM_MACRO)
#define ANALYTIC_BEP_SEP_CH8_USE_PARAM_MACRO 0
#endif

#if !defined(ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE)
#define ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE 0
#endif 

namespace gms 
{

namespace fading_channel
{

enum class Gaussian_Q_approxmations_t : std::int32_t 
{
    Gaussian_Q_approx_chiani,
    Gaussian_Q_approx_loskot_2T,
    Gaussian_Q_approx_loskot_3T,
    Gaussian_Q_approx_sadhwani_1T,
    Gaussian_Q_approx_sadhwani_2T,
    Gaussian_Q_approx_sadhwani_4T,
    Gaussian_Q_approx_cooper,
    Gaussian_Q_approx_borjesson,
    Gaussian_Q_approx_sadhwani_summed
};

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif
template<typename _Tp>
__ATTR_ALWAYS_INLINE__
static inline 
_Tp avg_symbol_E_to_carrier_A(const _Tp Ac,const _Tp Ts,
                                     const _Tp M) 
{
    constexpr _Tp third = static_cast<_Tp>(0.3333333333333333333333333333333333333333);
    constexpr _Tp one   = 1.0;
    const _Tp     AA    = Ac*Ac;
    const _Tp     MM    = M*M;
    const _Tp     MM_rat= (MM-one)*third;
    const _Tp     term  = AA*Ts;
    return (term*MM_rat);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif
template<typename _Tp>
__ATTR_ALWAYS_INLINE__
static inline 
_Tp avg_bit_E_to_carrier_A(const _Tp Ac,const _Tp Ts,
                                     const _Tp M)
{
    const _Tp Eg    = Ac*Ac*Ts;
    const _Tp num   = M*M-1.0;
    const _Tp log2M = 6.0*std::log2<_Tp>(M);
    return (Eg*(num/log2M));
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif
template<typename _Tp>
__ATTR_ALWAYS_INLINE__
static inline 
_Tp tikhonov_phase_err_pdf(const _Tp rho_eq,const _Tp phic)
{
    const _Tp cos2phic = std::cos<_Tp>(phic+phic);
    const _Tp exp_val  = std::exp<_Tp>(rho_eq*cos2phic);
    const _Tp bess_I0  = 3.1415926535897932384626433832795*std::cyl_bessel_i<_Tp>(0,rho_eq);
    return (exp_val/bess_I0);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif
template<typename _Tp>
__ATTR_ALWAYS_INLINE__
static inline 
_Tp G_PLL_bandwidth(const _Tp Bl,const _Tp Tb)
{
    return (static_cast<_Tp>(1.0)/(Bl*Tb));
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif
template<typename _Tp>
__ATTR_ALWAYS_INLINE__
static inline 
_Tp BPSK_param_a_8_61(const _Tp Ac,const _Tp Ts,
                      const _Tp M, const _Tp N0,
                      const _Tp Bl,const _Tp Tb)
{
    const _Tp Eb = avg_bit_E_to_carrier_A<_Tp>(Ac,Ts,M);
    const _Tp G  = G_PLL_bandwidth<_Tp>(Bl,Tb);
    const _Tp ratio = Eb/(N0+N0);
    const _Tp sqrtG = std::sqrt(G)-static_cast<_Tp>(1.0);
    return (ratio-(sqrtG*sqrtG));
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif
template<typename _Tp>
__ATTR_ALWAYS_INLINE__
static inline 
_Tp BPSK_param_b_8_61(const _Tp Ac,const _Tp Ts,
                      const _Tp M, const _Tp N0,
                      const _Tp Bl,const _Tp Tb)
{
    const _Tp Eb = avg_bit_E_to_carrier_A<_Tp>(Ac,Ts,M);
    const _Tp G  = G_PLL_bandwidth<_Tp>(Bl,Tb);
    const _Tp ratio = Eb/(N0+N0);
    const _Tp sqrtG = std::sqrt(G)+static_cast<_Tp>(1.0);
    return (ratio-(sqrtG*sqrtG));
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif
template<typename _Tp>
__ATTR_ALWAYS_INLINE__
static inline 
_Tp OQPSK_param_a1_8_64(const _Tp Ac,const _Tp Ts,
                        const _Tp M, const _Tp N0,
                        const _Tp Bl,const _Tp Tb)
{
    return (BPSK_param_a_8_61<_Tp>(Ac,Ts,M,N0,Bl,Tb));
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif
template<typename _Tp>
__ATTR_ALWAYS_INLINE__
static inline 
_Tp OQPSK_param_b1_8_64(const _Tp Ac,const _Tp Ts,
                        const _Tp M, const _Tp N0,
                        const _Tp Bl,const _Tp Tb)
{
    return (BPSK_param_b_8_61<_Tp>(Ac,Ts,M,N0,Bl,Tb));
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif
template<typename _Tp>
__ATTR_ALWAYS_INLINE__
static inline 
_Tp OQPSK_param_a2_8_64(const _Tp Ac,const _Tp Ts,
                        const _Tp M, const _Tp N0,
                        const _Tp Bl,const _Tp Tb)
{
    const _Tp Eb = avg_bit_E_to_carrier_A<_Tp>(Ac,Ts,M);
    const _Tp G  = G_PLL_bandwidth<_Tp>(Bl,Tb);
    const _Tp fac= (G+static_cast<_Tp>(1.0))-std::sqrt<_Tp>(G+G);
    return ((Eb/N0)*fac);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif
template<typename _Tp>
__ATTR_ALWAYS_INLINE__
static inline 
_Tp OQPSK_param_b2_8_64(const _Tp Ac,const _Tp Ts,
                        const _Tp M, const _Tp N0,
                        const _Tp Bl,const _Tp Tb)
{
    const _Tp Eb = avg_bit_E_to_carrier_A<_Tp>(Ac,Ts,M);
    const _Tp G  = G_PLL_bandwidth<_Tp>(Bl,Tb);
    const _Tp fac= (G+static_cast<_Tp>(1.0))+std::sqrt<_Tp>(G+G);
    return ((Eb/N0)*fac);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif
template<typename _Tp>
__ATTR_ALWAYS_INLINE__
static inline 
_Tp MSK_param_a1_8_65(  const _Tp Ac,const _Tp Ts,
                        const _Tp M, const _Tp N0,
                        const _Tp Bl,const _Tp Tb)
{
    return (BPSK_param_a_8_61<_Tp>(Ac,Ts,M,N0,Bl,Tb));
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif
template<typename _Tp>
__ATTR_ALWAYS_INLINE__
static inline 
_Tp MSK_param_b1_8_65(  const _Tp Ac,const _Tp Ts,
                        const _Tp M, const _Tp N0,
                        const _Tp Bl,const _Tp Tb)
{
    return (BPSK_param_b_8_61<_Tp>(Ac,Ts,M,N0,Bl,Tb));
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif
template<typename _Tp>
__ATTR_ALWAYS_INLINE__
static inline 
_Tp MSK_param_a2_8_65(  const _Tp Ac,const _Tp Ts,
                        const _Tp M, const _Tp N0,
                        const _Tp Bl,const _Tp Tb)
{
    const _Tp Eb          = avg_bit_E_to_carrier_A<_Tp>(Ac,Ts,M);
    const _Tp G           = G_PLL_bandwidth<_Tp>(Bl,Tb);
    constexpr _Tp pi_term = static_cast<_Tp>(0.7026423672846755428877589264195);
    const _Tp fac         = G+pi_term-std::sqrt<_Tp>(G+G);
    return ((Eb/N0)*fac);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif
template<typename _Tp>
__ATTR_ALWAYS_INLINE__
static inline 
_Tp MSK_param_b2_8_65(  const _Tp Ac,const _Tp Ts,
                        const _Tp M, const _Tp N0,
                        const _Tp Bl,const _Tp Tb)
{
    const _Tp Eb          = avg_bit_E_to_carrier_A<_Tp>(Ac,Ts,M);
    const _Tp G           = G_PLL_bandwidth<_Tp>(Bl,Tb);
    constexpr _Tp pi_term = static_cast<_Tp>(0.7026423672846755428877589264195);
    const _Tp fac         = G+pi_term+std::sqrt<_Tp>(G+G);
    return ((Eb/N0)*fac);
}


#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline 
double analytic_SEP_MAM_8_1(const double Ac,const double Ts,
                            const double M,const double N0,
                            const std::int32_t n) 
{   
    double result{};
    double Q_func_val = -1.0;
    const double M_ratio = (M-1.0)/M;
    const double left_term = M_ratio+M_ratio;
    const double num = 2.0*(Ac*Ac)*Ts;
    const double sqr_ratio = num/N0;
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double Q_func_arg = gms::math::cephes_d::sqrt(sqr_ratio);
#else 
    const double Q_func_arg = std::sqrt(sqr_ratio);
#endif
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
       Q_func_val = gms::math::gaussian_Q_approx_chiani(Q_func_arg);
       result = left_term*Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T(Q_func_arg);
        result = left_term*Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T(Q_func_arg);
        result = left_term*Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T(Q_func_arg);
        result = left_term*Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T(Q_func_arg);
        result = left_term*Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T(Q_func_arg);
        result = left_term*Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_cooper)
    {
        Q_func_val = gms::math::gaussian_Q_approx_cooper(Q_func_arg);
        result = left_term*Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson)
    {
        Q_func_val = gms::math::gaussian_Q_approx_borjesson(Q_func_arg);
        result = left_term*Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed(Q_func_arg,n);
        result = left_term*Q_func_val;
    }
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline 
float analytic_SEP_MAM_8_1(const float Ac,const float Ts,
                           const float M,const float N0,
                           const std::int32_t n) 
{   
    float result{};
    float Q_func_val = -1.0f;
    const float M_ratio = (M-1.0f)/M;
    const float left_term = M_ratio+M_ratio;
    const float num = 2.0f*(Ac*Ac)*Ts;
    const float sqr_ratio = num/N0;
    const float Q_func_arg = std::sqrt(sqr_ratio);
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
       Q_func_val = gms::math::gaussian_Q_approx_chiani(Q_func_arg);
       result = left_term*Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T(Q_func_arg);
        result = left_term*Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T(Q_func_arg);
        result = left_term*Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T(Q_func_arg);
        result = left_term*Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T(Q_func_arg);
        result = left_term*Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T(Q_func_arg);
        result = left_term*Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed(Q_func_arg,n);
        result = left_term*Q_func_val;
    }
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
double analytic_SEP_MAM_8_3(const double Ac,const double Ts,
                            const double M,const double N0,
                            const std::int32_t n)
{
    double result;
    double Q_func_val     = -1.0;
    double Es             = avg_symbol_E_to_carrier_A<double>(Ac,Ts,M);
    const double sixEs    = 6.0*Es;
    const double MM_m_one = M*M-1.0;
    const double left_term= 2.0*((M-1.0)/M);
    const double ratio    = sixEs/(N0*MM_m_one);
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double Q_func_arg = gms::math::cephes_d::sqrt(ratio);
#else 
    const double Q_func_arg = std::sqrt(ratio);
#endif 
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
       Q_func_val = gms::math::gaussian_Q_approx_chiani(Q_func_arg);
       result = left_term*Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T(Q_func_arg);
        result = left_term*Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T(Q_func_arg);
        result = left_term*Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T(Q_func_arg);
        result = left_term*Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T(Q_func_arg);
        result = left_term*Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T(Q_func_arg);
        result = left_term*Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed(Q_func_arg,n);
        result = left_term*Q_func_val;
    }
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
float analytic_SEP_MAM_8_3(const float Ac,const float Ts,
                            const float M,const float N0,
                            const std::int32_t n)
{
    float result;
    float Q_func_val     = -1.0f;
    float Es             = avg_symbol_E_to_carrier_A<float>(Ac,Ts,M);
    const float sixEs    = 6.0f*Es;
    const float MM_m_one = M*M-1.0f;
    const float left_term= 2.0f*((M-1.0f)/M);
    const float ratio    = sixEs/(N0*MM_m_one); 
    const double Q_func_arg = std::sqrt(ratio);
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
       Q_func_val = gms::math::gaussian_Q_approx_chiani(Q_func_arg);
       result = left_term*Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T(Q_func_arg);
        result = left_term*Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T(Q_func_arg);
        result = left_term*Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T(Q_func_arg);
        result = left_term*Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T(Q_func_arg);
        result = left_term*Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T(Q_func_arg);
        result = left_term*Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed(Q_func_arg,n);
        result = left_term*Q_func_val;
    }
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
double analytic_BEP_MAM2_8_4(const double Ac,const double Ts,
                            const double M,const double N0,
                            const std::int32_t n)
{
    double result;
    double Q_func_val{-1.0};
    double Eb    = avg_bit_E_to_carrier_A<double>(Ac,Ts,M);
    double ratio = (Eb+Eb)/N0; 
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double Q_func_arg = gms::math::cephes_d::sqrt(ratio);
#else 
    const double Q_func_arg = std::sqrt(ratio);
#endif
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
       Q_func_val = gms::math::gaussian_Q_approx_chiani(Q_func_arg);
       result = Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T(Q_func_arg);
        result = Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T(Q_func_arg);
        result = Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T(Q_func_arg);
        result = Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T(Q_func_arg);
        result = Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T(Q_func_arg);
        result = Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed(Q_func_arg,n);
        result = Q_func_val;
    }
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
float analytic_BEP_MAM2_8_4(const float Ac,const float Ts,
                            const float M,const float N0,
                            const std::int32_t n)
{
    float result;
    float Q_func_val{-1.0};
    float Eb    = avg_bit_E_to_carrier_A<float>(Ac,Ts,M);
    float ratio = (Eb+Eb)/N0; 
    const double Q_func_arg = std::sqrt(ratio);
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
       Q_func_val = gms::math::gaussian_Q_approx_chiani(Q_func_arg);
       result = Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T(Q_func_arg);
        result = Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T(Q_func_arg);
        result = Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T(Q_func_arg);
        result = Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T(Q_func_arg);
        result = Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T(Q_func_arg);
        result = Q_func_val;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed(Q_func_arg,n);
        result = Q_func_val;
    }
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
double analytic_SEP_QAM_8_10(const double Ac,const double Ts,
                            const double M,const double N0,
                            const std::int32_t n)
{
    double result;
    double Q_func_val    = -1.0;
    const double Es      = avg_symbol_E_to_carrier_A<double>(Ac,Ts,M);
    const double M_m_1   = M-1.0;
    const double threeEs = Es+Es+Es;
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double sqrtM   = gms::math::cephes_d::sqrt(M);
#else 
    const double sqrtM   = std::sqrt(M);
#endif 
    const double left_term= 4.0*((sqrtM-1.0)/sqrtM);
    const double ratio    = threeEs/(N0*M_m_1);
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double Q_func_arg= gms::math::cephes_d::sqrt(ratio);
#else 
    const double Q_func_arg= std::sqrt(ratio);
#endif 
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
       Q_func_val        = gms::math::gaussian_Q_approx_chiani(Q_func_arg);
       const double Q_func_val_p2 = Q_func_val*Q_func_val;
       const double left_term_p2  = left_term*left_term;
       const double diff_result   = (left_term*Q_func_val)-(left_term_p2*Q_func_val_p2); 
       result = diff_result;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T(Q_func_arg);
        const double Q_func_val_p2 = Q_func_val*Q_func_val;
        const double left_term_p2  = left_term*left_term;
        const double diff_result   = (left_term*Q_func_val)-(left_term_p2*Q_func_val_p2); 
        result = diff_result;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T(Q_func_arg);
        const double Q_func_val_p2 = Q_func_val*Q_func_val;
        const double left_term_p2  = left_term*left_term;
        const double diff_result   = (left_term*Q_func_val)-(left_term_p2*Q_func_val_p2); 
        result = diff_result;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T(Q_func_arg);
        const double Q_func_val_p2 = Q_func_val*Q_func_val;
        const double left_term_p2  = left_term*left_term;
        const double diff_result   = (left_term*Q_func_val)-(left_term_p2*Q_func_val_p2); 
        result = diff_result;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T(Q_func_arg);
        const double Q_func_val_p2 = Q_func_val*Q_func_val;
        const double left_term_p2  = left_term*left_term;
        const double diff_result   = (left_term*Q_func_val)-(left_term_p2*Q_func_val_p2); 
        result = diff_result;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T(Q_func_arg);
        const double Q_func_val_p2 = Q_func_val*Q_func_val;
        const double left_term_p2  = left_term*left_term;
        const double diff_result   = (left_term*Q_func_val)-(left_term_p2*Q_func_val_p2); 
        result = diff_result;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed(Q_func_arg,n);
        const double Q_func_val_p2 = Q_func_val*Q_func_val;
        const double left_term_p2  = left_term*left_term;
        const double diff_result   = (left_term*Q_func_val)-(left_term_p2*Q_func_val_p2); 
        result = diff_result;
    }
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
float analytic_SEP_QAM_8_10(const float Ac,const float Ts,
                            const float M,const float N0,
                            const std::int32_t n)
{
    float result;
    float Q_func_val    = -1.0f;
    const float Es      = avg_symbol_E_to_carrier_A<float>(Ac,Ts,M);
    const float M_m_1   = M-1.0f;
    const float threeEs = Es+Es+Es;
    const float sqrtM   = std::sqrt(M);
    const float left_term= 4.0f*((sqrtM-1.0f)/sqrtM);
    const float ratio    = threeEs/(N0*M_m_1);
    const float Q_func_arg= std::sqrt(ratio);
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
       Q_func_val        = gms::math::gaussian_Q_approx_chiani(Q_func_arg);
       const float Q_func_val_p2 = Q_func_val*Q_func_val;
       const float left_term_p2  = left_term*left_term;
       const float diff_result   = (left_term*Q_func_val)-(left_term_p2*Q_func_val_p2); 
       result = diff_result;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T(Q_func_arg);
        const float Q_func_val_p2 = Q_func_val*Q_func_val;
        const float left_term_p2  = left_term*left_term;
        const float diff_result   = (left_term*Q_func_val)-(left_term_p2*Q_func_val_p2); 
        result = diff_result;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T(Q_func_arg);
        const float Q_func_val_p2 = Q_func_val*Q_func_val;
        const float left_term_p2  = left_term*left_term;
        const float diff_result   = (left_term*Q_func_val)-(left_term_p2*Q_func_val_p2); 
        result = diff_result;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T(Q_func_arg);
        const float Q_func_val_p2 = Q_func_val*Q_func_val;
        const float left_term_p2  = left_term*left_term;
        const float diff_result   = (left_term*Q_func_val)-(left_term_p2*Q_func_val_p2); 
        result = diff_result;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T(Q_func_arg);
        const float Q_func_val_p2 = Q_func_val*Q_func_val;
        const float left_term_p2  = left_term*left_term;
        const float diff_result   = (left_term*Q_func_val)-(left_term_p2*Q_func_val_p2); 
        result = diff_result;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T(Q_func_arg);
        const float Q_func_val_p2 = Q_func_val*Q_func_val;
        const float left_term_p2  = left_term*left_term;
        const float diff_result   = (left_term*Q_func_val)-(left_term_p2*Q_func_val_p2); 
        result = diff_result;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed(Q_func_arg,n);
        const float Q_func_val_p2 = Q_func_val*Q_func_val;
        const float left_term_p2  = left_term*left_term;
        const float diff_result   = (left_term*Q_func_val)-(left_term_p2*Q_func_val_p2); 
        result = diff_result;
    }
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
double analytic_SEP_QAM4_8_11(const double Ac,const double Ts,
                              const double M,const double N0,
                              const std::int32_t n)
{
    double result;
    double Q_func_val = -1.0;
    const double Es   = avg_symbol_E_to_carrier_A<double>(Ac,Ts,M);
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double Q_func_arg = gms::math::cephes_d::sqrt(Es/N0);
#else 
    const double Q_func_arg = std::sqrt(Es/N0);
#endif 
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
       Q_func_val        = gms::math::gaussian_Q_approx_chiani(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed(Q_func_arg,n);
    }
    const double Q_func_val_a2  = Q_func_val+Q_func_val;
    const double Q_func_val_p2  = Q_func_val*Q_func_val;
    const double diff_result   =  Q_func_val_a2-Q_func_val_p2; 
    result = diff_result
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
float analytic_SEP_QAM4_8_11(const float Ac,const float Ts,
                              const float M,const float N0,
                              const std::int32_t n)
{
    float result;
    float Q_func_val = -1.0;
    const float Es   = avg_symbol_E_to_carrier_A<float>(Ac,Ts,M); 
    const float Q_func_arg = std::sqrt(Es/N0);
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
       Q_func_val        = gms::math::gaussian_Q_approx_chiani(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_cooper)
    {
        Q_func_val = gms::math::gaussian_Q_approx_cooper(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson)
    {
        Q_func_val = gms::math::gaussian_Q_approx_borjesson(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed(Q_func_arg,n);
    }
    const float Q_func_val_a2  = Q_func_val+Q_func_val;
    const float Q_func_val_p2  = Q_func_val*Q_func_val;
    const float diff_result   =  Q_func_val_a2-Q_func_val_p2; 
    result = diff_result;
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
double analytic_BEP_QAM_AWGN_8_14(const double Ac,const double Ts,
                                  const double M,const double N0,
                                  const std::int32_t n)
{
    double result;
    double Q_func_val;
    const double Eb       = avg_bit_E_to_carrier_A<double>(Ac,Ts,M);
    const double log2M    = std::log2<double>(M);
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double sqrtM    =  gms::math::cephes_d::sqrt(M);
#else 
    const double sqrtM    =  std::sqrt(M);
#endif 
    const std::int32_t upper_lim = static_cast<std::int32_t>(sqrtM*0.5);
    const double invlog2M  = 1.0/log2M;
    if(__builtin_expect(upper_lim==0,0)) { return (std::numeric_limits<double>::quiet_NaN());}
    const double sqrtM_m_1 = sqrtM-1.0;
    const double left_term = 4.0*(sqrtM_m_1/sqrtM);
    const double num       = 3.0*Eb*log2M;
    const double den       = N0*sqrtM_m_1;
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double Q_func_arg= gms::math::cephes_d::sqrt(num/den);
#else 
    const double Q_func_arg= std::sqrt(num/den);
#endif 
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
       for(std::int32_t i=1; i<=upper_lim; ++i) 
       {
           const double d_i = static_cast<double>(i);
           const double mul_fac = (d_i+d_i)-1.0;
           Q_func_val  += gms::math::gaussian_Q_approx_chiani(mul_fac*Q_func_arg);
       }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
       for(std::int32_t i=1; i<=upper_lim; ++i) 
       {
           const double d_i = static_cast<double>(i);
           const double mul_fac = (d_i+d_i)-1.0;
           Q_func_val += gms::math::gaussian_Q_approx_loskot_2T(mul_fac*Q_func_arg);
       }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        for(std::int32_t i=1; i<=upper_lim; ++i) 
       {
           const double d_i = static_cast<double>(i);
           const double mul_fac = (d_i+d_i)-1.0;
           Q_func_val += gms::math::gaussian_Q_approx_loskot_3T(mul_fac*Q_func_arg);
       }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        for(std::int32_t i=1; i<=upper_lim; ++i) 
       {
           const double d_i = static_cast<double>(i);
           const double mul_fac = (d_i+d_i)-1.0;
           Q_func_val += gms::math::gaussian_Q_approx_sadhwani_1T(mul_fac*Q_func_arg);
       }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
       for(std::int32_t i=1; i<=upper_lim; ++i) 
       {
           const double d_i = static_cast<double>(i);
           const double mul_fac = (d_i+d_i)-1.0;
           Q_func_val += gms::math::gaussian_Q_approx_sadhwani_2T(mul_fac*Q_func_arg);
       }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        for(std::int32_t i=1; i<=upper_lim; ++i) 
        { 
           const double d_i = static_cast<double>(i);
           const double mul_fac = (d_i+d_i)-1.0;
           Q_func_val += gms::math::gaussian_Q_approx_sadhwani_4T(mul_fac*Q_func_arg);
        }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_cooper)
    {
        for(std::int32_t i=1; i<=upper_lim; ++i) 
        { 
           const double d_i = static_cast<double>(i);
           const double mul_fac = (d_i+d_i)-1.0;
           Q_func_val += gms::math::gaussian_Q_approx_cooper(mul_fac*Q_func_arg);
        }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson)
    {
        for(std::int32_t i=1; i<=upper_lim; ++i) 
        { 
           const double d_i = static_cast<double>(i);
           const double mul_fac = (d_i+d_i)-1.0;
           Q_func_val += gms::math::gaussian_Q_approx_borjesson(mul_fac*Q_func_arg);
        }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        for(std::int32_t i=1; i<=upper_lim; ++i) 
        { 
           const double d_i = static_cast<double>(i);
           const double mul_fac = (d_i+d_i)-1.0;
           Q_func_val += gms::math::gaussian_Q_approx_sadhwani_summed(mul_fac*Q_func_arg,n);
        }
    }
    result = left_term*invlog2M*Q_func_val;
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
float analytic_BEP_QAM_AWGN_8_14(const float Ac,const float Ts,
                                  const float M,const float N0,
                                  const std::int32_t n)
{
    float result;
    float Q_func_val     = 0.0f;
    const float Eb       = avg_bit_E_to_carrier_A<float>(Ac,Ts,M);
    const float log2M    = std::log2<float>(M);
    const float sqrtM    =  std::sqrt(M);
    const std::int32_t upper_lim = static_cast<std::int32_t>(sqrtM*0.5f);
    const float invlog2M  = 1.0f/log2M;
    if(__builtin_expect(upper_lim==0,0)) { return (std::numeric_limits<float>::quiet_NaN());}
    const float sqrtM_m_1 = sqrtM-1.0f;
    const float left_term = 4.0f*(sqrtM_m_1/sqrtM);
    const float num       = 3.0f*Eb*log2M;
    const float den       = N0*sqrtM_m_1;
    const float Q_func_arg= std::sqrt(num/den);
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
       for(std::int32_t i=1; i<=upper_lim; ++i) 
       {
           const float d_i = static_cast<double>(i);
           const float mul_fac = (d_i+d_i)-1.0f;
           Q_func_val  += gms::math::gaussian_Q_approx_chiani(mul_fac*Q_func_arg);
       }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
       for(std::int32_t i=1; i<=upper_lim; ++i) 
       {
           const float d_i = static_cast<double>(i);
           const float mul_fac = (d_i+d_i)-1.0f;
           Q_func_val += gms::math::gaussian_Q_approx_loskot_2T(mul_fac*Q_func_arg);
       }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        for(std::int32_t i=1; i<=upper_lim; ++i) 
       {
           const float d_i = static_cast<double>(i);
           const float mul_fac = (d_i+d_i)-1.0f;
           Q_func_val += gms::math::gaussian_Q_approx_loskot_3T(mul_fac*Q_func_arg);
       }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        for(std::int32_t i=1; i<=upper_lim; ++i) 
       {
           const float d_i = static_cast<double>(i);
           const float mul_fac = (d_i+d_i)-1.0f;
           Q_func_val += gms::math::gaussian_Q_approx_sadhwani_1T(mul_fac*Q_func_arg);
       }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
       for(std::int32_t i=1; i<=upper_lim; ++i) 
       {
           const float d_i = static_cast<double>(i);
           const float mul_fac = (d_i+d_i)-1.0f;
           Q_func_val += gms::math::gaussian_Q_approx_sadhwani_2T(mul_fac*Q_func_arg);
       }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        for(std::int32_t i=1; i<=upper_lim; ++i) 
        { 
           const float d_i = static_cast<double>(i);
           const float mul_fac = (d_i+d_i)-1.0f;
           Q_func_val += gms::math::gaussian_Q_approx_sadhwani_4T(mul_fac*Q_func_arg);
        }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_cooper)
    {
        for(std::int32_t i=1; i<=upper_lim; ++i) 
        { 
           const float d_i = static_cast<double>(i);
           const float mul_fac = (d_i+d_i)-1.0f;
           Q_func_val += gms::math::gaussian_Q_approx_cooper(mul_fac*Q_func_arg);
        }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson)
    {
        for(std::int32_t i=1; i<=upper_lim; ++i) 
        { 
           const float d_i = static_cast<double>(i);
           const float mul_fac = (d_i+d_i)-1.0f;
           Q_func_val += gms::math::gaussian_Q_approx_borjesson(mul_fac*Q_func_arg);
        }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        for(std::int32_t i=1; i<=upper_lim; ++i) 
        { 
           const float d_i = static_cast<double>(i);
           const float mul_fac = (d_i+d_i)-1.0f;
           Q_func_val += gms::math::gaussian_Q_approx_sadhwani_summed(mul_fac*Q_func_arg,n);
        }
    }
    result = left_term*invlog2M*Q_func_val;
    return (result);
}


#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
double analytic_BEP_MPSK2_8_18(const double Ac,const double Ts,
                               const double M,const double N0,
                               const std::int32_t n)
{
    double result;
    double Q_func_val  = -1.0;
    const double Eb    = avg_bit_E_to_carrier_A<double>(Ac,Ts,M);
    const double ratio = (Eb+Eb)/N0;
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double Q_func_arg = gms::math::cephes_d::sqrt(ratio);
#else 
    const double Q_func_arg = std::sqrt(ratio);
#endif 
     if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
       Q_func_val        = gms::math::gaussian_Q_approx_chiani(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_cooper)
    {
        Q_func_val = gms::math::gaussian_Q_approx_cooper(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson)
    {
        Q_func_val = gms::math::gaussian_Q_approx_borjesson(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed(Q_func_arg,n);
    }
    result = Q_func_val;
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
float analytic_BEP_MPSK2_8_18(const float Ac,const float Ts,
                               const float M,const float N0,
                               const std::int32_t n)
{
    float result;
    float Q_func_val  = -1.0f;
    const float Eb    = avg_bit_E_to_carrier_A<float>(Ac,Ts,M);
    const float ratio = (Eb+Eb)/N0;
    const float Q_func_arg = std::sqrt(ratio);
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
       Q_func_val        = gms::math::gaussian_Q_approx_chiani(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_cooper)
    {
        Q_func_val = gms::math::gaussian_Q_approx_cooper(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson)
    {
        Q_func_val = gms::math::gaussian_Q_approx_borjesson(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed(Q_func_arg,n);
    }
    result = Q_func_val;
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
double analytic_SEP_QPSKM4_8_19(const double Ac,const double Ts,
                                const double M,const double N0,
                                const std::int32_t n)
{
    double result;
    double Q_func_val;
    double sqrQ;
    const double Es = avg_symbol_E_to_carrier_A<double>(Ac,Ts,M);
    const double ratio = Es/N0;
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double Q_func_arg = gms::math::cephes_d::sqrt(ratio);
#else 
    const double Q_func_arg = std::sqrt(ratio);
#endif 
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
       Q_func_val        = gms::math::gaussian_Q_approx_chiani(Q_func_arg);
       sqrQ              = Q_func_val*Q_func_val;
       result            = (Q_func_val+Q_func_val)-sqrQ;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T(Q_func_arg);
        sqrQ              = Q_func_val*Q_func_val;
        result            = (Q_func_val+Q_func_val)-sqrQ;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T(Q_func_arg);
        sqrQ              = Q_func_val*Q_func_val;
        result            = (Q_func_val+Q_func_val)-sqrQ;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T(Q_func_arg);
        sqrQ              = Q_func_val*Q_func_val;
        result            = (Q_func_val+Q_func_val)-sqrQ;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T(Q_func_arg);
        sqrQ              = Q_func_val*Q_func_val;
        result            = (Q_func_val+Q_func_val)-sqrQ;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T(Q_func_arg);
        sqrQ              = Q_func_val*Q_func_val;
        result            = (Q_func_val+Q_func_val)-sqrQ;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_cooper)
    {
        Q_func_val = gms::math::gaussian_Q_approx_cooper(Q_func_arg);
        sqrQ              = Q_func_val*Q_func_val;
        result            = (Q_func_val+Q_func_val)-sqrQ;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson)
    {
        Q_func_val = gms::math::gaussian_Q_approx_borjesson(Q_func_arg);
        sqrQ              = Q_func_val*Q_func_val;
        result            = (Q_func_val+Q_func_val)-sqrQ;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed(Q_func_arg,n);
        sqrQ              = Q_func_val*Q_func_val;
        result            = (Q_func_val+Q_func_val)-sqrQ;
    }
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
float analytic_SEP_QPSKM4_8_19(const float Ac,const float Ts,
                                const float M,const float N0,
                                const std::int32_t n)
{
    float result;
    float Q_func_val;
    float sqrQ;
    const float Es = avg_symbol_E_to_carrier_A<float>(Ac,Ts,M);
    const float ratio = Es/N0; 
    const float Q_func_arg = std::sqrt(ratio);
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
       Q_func_val        = gms::math::gaussian_Q_approx_chiani(Q_func_arg);
       sqrQ              = Q_func_val*Q_func_val;
       result            = (Q_func_val+Q_func_val)-sqrQ;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T(Q_func_arg);
        sqrQ              = Q_func_val*Q_func_val;
        result            = (Q_func_val+Q_func_val)-sqrQ;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T(Q_func_arg);
        sqrQ              = Q_func_val*Q_func_val;
        result            = (Q_func_val+Q_func_val)-sqrQ;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T(Q_func_arg);
        sqrQ              = Q_func_val*Q_func_val;
        result            = (Q_func_val+Q_func_val)-sqrQ;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T(Q_func_arg);
        sqrQ              = Q_func_val*Q_func_val;
        result            = (Q_func_val+Q_func_val)-sqrQ;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T(Q_func_arg);
        sqrQ              = Q_func_val*Q_func_val;
        result            = (Q_func_val+Q_func_val)-sqrQ;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_cooper)
    {
        Q_func_val = gms::math::gaussian_Q_approx_cooper(Q_func_arg);
        sqrQ              = Q_func_val*Q_func_val;
        result            = (Q_func_val+Q_func_val)-sqrQ;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson)
    {
        Q_func_val = gms::math::gaussian_Q_approx_borjesson(Q_func_arg);
        sqrQ              = Q_func_val*Q_func_val;
        result            = (Q_func_val+Q_func_val)-sqrQ;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed(Q_func_arg,n);
        sqrQ              = Q_func_val*Q_func_val;
        result            = (Q_func_val+Q_func_val)-sqrQ;
    }
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
double analytic_SEP_MPSK_coherent_8_25(const double Ac,const double Ts,
                                       const double M,const double N0,
                                       const std::int32_t n)
{
    double result;
    double Q_func_val;
    const double Es = avg_symbol_E_to_carrier_A<double>(Ac,Ts,M);
    const double PI_over_M = 3.1415926535897932384626433832795/M;
    const double ratio     = (Es+Es)/N0;
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double Q_func_arg = gms::math::cephes_d::sqrt(ratio);
    const double sin_PIM    = gms::math::cephes_d::sin(PI_over_M);
#else 
    const double Q_func_arg = std::sqrt(ratio);
    const double sin_PIM    = std::sin(PI_over_M);
#endif 
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
       Q_func_val        = gms::math::gaussian_Q_approx_chiani(Q_func_arg*sin_PIM);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T(Q_func_arg*sin_PIM);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T(Q_func_arg*sin_PIM);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T(Q_func_arg*sin_PIM);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T(Q_func_arg*sin_PIM);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T(Q_func_arg*sin_PIM);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_cooper)
    {
        Q_func_val = gms::math::gaussian_Q_approx_cooper(Q_func_arg*sin_PIM);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson)
    {
        Q_func_val = gms::math::gaussian_Q_approx_borjesson(Q_func_arg*sin_PIM);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed(Q_func_arg*sin_PIM,n);
    }
    result = Q_func_val+Q_func_val;
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
float analytic_SEP_MPSK_coherent_8_25(const float Ac,const float Ts,
                                       const float M,const float N0,
                                       const std::int32_t n)
{
    float result;
    float Q_func_val;
    const float Es = avg_symbol_E_to_carrier_A<float>(Ac,Ts,M);
    const float PI_over_M = 3.1415926535897932384626433832795f/M;
    const float ratio     = (Es+Es)/N0;
    const double Q_func_arg = std::sqrt(ratio);
    const double sin_PIM    = std::sin(PI_over_M);
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
       Q_func_val        = gms::math::gaussian_Q_approx_chiani(Q_func_arg*sin_PIM);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T(Q_func_arg*sin_PIM);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T(Q_func_arg*sin_PIM);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T(Q_func_arg*sin_PIM);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T(Q_func_arg*sin_PIM);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T(Q_func_arg*sin_PIM);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_cooper)
    {
        Q_func_val = gms::math::gaussian_Q_approx_cooper(Q_func_arg*sin_PIM);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson)
    {
        Q_func_val = gms::math::gaussian_Q_approx_borjesson(Q_func_arg*sin_PIM);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed(Q_func_arg*sin_PIM,n);
    }
    result = Q_func_val+Q_func_val;
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
double analytic_BEP_MPSK_AWGN_8_31(const double Ac,const double Ts,
                                   const double M,const double N0,
                                   const std::int32_t n)
{   
    if(__builtin_expect(M>1024.0,0)) { return (-1.0);}
    constexpr std::int32_t max_M = 256;
constexpr const double sin_arg_lut[m_Max] = 
{
3.1415926535897931159979634685442,
9.4247779607693793479938904056326,
15.7079632679489655799898173427209,
21.9911485751285518119857442798093,
28.2743338823081380439816712168977,
34.5575191894877207232639193534851,
40.8407044966673140606872038915753,
47.1238898038468931872557732276618,
53.4070751110264865246790577657521,
59.6902604182060656512476271018386,
65.9734457253856589886709116399288,
72.2566310325652381152394809760153,
78.5398163397448314526627655141056,
84.8230016469244105792313348501921,
91.1061869541040039166546193882823,
97.3893722612835830432231887243688,
103.6725575684631763806464732624590,
109.9557428756427555072150425985456,
116.2389281828223488446383271366358,
122.5221134900019279712068964727223,
128.8052987971815070977754658088088,
135.0884841043611004351987503468990,
141.3716694115406937726220348849893,
147.6548547187202871100453194230795,
153.9380400258998520257591735571623,
160.2212253330794453631824580952525,
166.5044106402590387006057426333427,
172.7875959474386320380290271714330,
179.0707812546181969537428813055158,
185.3539665617977902911661658436060,
191.6371518689773836285894503816962,
197.9203371761569769660127349197865,
204.2035224833365418817265890538692,
210.4867077905161352191498735919595,
216.7698930976957285565731581300497,
223.0530784048753218939964426681399,
229.3362637120548868097102968022227,
235.6194490192344801471335813403130,
241.9026343264140734845568658784032,
248.1858196335936668219801504164934,
254.4690049407732317376940045505762,
260.7521902479528534968267194926739,
267.0353755551324184125405736267567,
273.3185608623119833282544277608395,
279.6017461694916050873871427029371,
285.8849314766711700031009968370199,
292.1681167838507349188148509711027,
298.4513020910303566779475659132004,
304.7344873982099215936614200472832,
311.0176727053895433527941349893808,
317.3008580125691082685079891234636,
323.5840433197486731842218432575464,
329.8672286269282949433545581996441,
336.1504139341078598590684123337269,
342.4335992412874247747822664678097,
348.7167845484670465339149814099073,
354.9999698556466114496288355439901,
361.2831551628262332087615504860878,
367.5663404700057981244754046201706,
373.8495257771853630401892587542534,
380.1327110843649847993219736963511,
386.4158963915445497150358278304338,
392.6990816987241146307496819645166,
398.9822670059037363898823969066143,
405.2654523130833013055962510406971,
411.5486376202629230647289659827948,
417.8318229274424879804428201168776,
424.1150082346220528961566742509604,
430.3981935418016746552893891930580,
436.6813788489812395710032433271408,
442.9645641561608044867170974612236,
449.2477494633404262458498124033213,
455.5309347705199911615636665374041,
461.8141200776996129206963814795017,
468.0973053848791778364102356135845,
474.3804906920587427521240897476673,
480.6636759992383645112568046897650,
486.9468613064179294269706588238478,
493.2300466135974943426845129579306,
499.5132319207771161018172279000282,
505.7964172279566810175310820341110,
512.0796025351362459332449361681938,
518.3627878423158108489587903022766,
524.6459731494954894515103660523891,
530.9291584566750543672242201864719,
537.2123437638546192829380743205547,
543.4955290710341841986519284546375,
549.7787143782137491143657825887203,
556.0618996853934277169173583388329,
562.3450849925729926326312124729156,
568.6282702997525575483450666069984,
574.9114556069321224640589207410812,
581.1946409141116873797727748751640,
587.4778262212913659823243506252766,
593.7610115284709308980382047593594,
600.0441968356504958137520588934422,
606.3273821428300607294659130275249,
612.6105674500096256451797671616077,
618.8937527571891905608936212956905,
625.1769380643688691634451970458031,
631.4601233715484340791590511798859,
637.7433086787279989948729053139687,
644.0264939859075639105867594480515,
650.3096792930871288263006135821342,
656.5928646002668074288521893322468,
662.8760499074463723445660434663296,
669.1592352146259372602798976004124,
675.4424205218055021759937517344952,
681.7256058289850670917076058685780,
688.0087911361647456942591816186905,
694.2919764433443106099730357527733,
700.5751617505238755256868898868561,
706.8583470577034404414007440209389,
713.1415323648830053571145981550217,
719.4247176720625702728284522891045,
725.7079029792422488753800280392170,
731.9910882864218137910938821732998,
738.2742735936013787068077363073826,
744.5574589007809436225215904414654,
750.8406442079605085382354445755482,
757.1238295151401871407870203256607,
763.4070148223197520565008744597435,
769.6902001294993169722147285938263,
775.9733854366788818879285827279091,
782.2565707438584468036424368619919,
788.5397560510381254061940126121044,
794.8229413582176903219078667461872,
801.1061266653972552376217208802700,
807.3893119725768201533355750143528,
813.6724972797563850690494291484356,
819.9556825869359499847632832825184,
826.2388678941156285873148590326309,
832.5220532012951935030287131667137,
838.8052385084747584187425673007965,
845.0884238156543233344564214348793,
851.3716091228338882501702755689621,
857.6547944300135668527218513190746,
863.9379797371931317684357054531574,
870.2211650443726966841495595872402,
876.5043503515522615998634137213230,
882.7875356587318265155772678554058,
889.0707209659115051181288436055183,
895.3539062730910700338426977396011,
901.6370915802706349495565518736839,
907.9202768874501998652704060077667,
914.2034621946297647809842601418495,
920.4866475018093296966981142759323,
926.7698328089890082992496900260448,
933.0530181161685732149635441601276,
939.3362034233481381306773982942104,
945.6193887305277030463912524282932,
951.9025740377072679621051065623760,
958.1857593448869465646566823124886,
964.4689446520665114803705364465714,
970.7521299592460763960843905806541,
977.0353152664256413117982447147369,
983.3185005736052062275120988488197,
989.6016858807848848300636745989323,
995.8848711879644497457775287330151,
1002.1680564951440146614913828670979,
1008.4512418023235795772052370011806,
1014.7344271095031444929190911352634,
1021.0176124166827094086329452693462,
1027.3007977238623880111845210194588,
1033.5839830310419529268983751535416,
1039.8671683382215178426122292876244,
1046.1503536454010827583260834217072,
1052.4335389525806476740399375557899,
1058.7167242597602125897537916898727,
1064.9999095669397775054676458239555,
1071.2830948741193424211814999580383,
1077.5662801812991347105707973241806,
1083.8494654884786996262846514582634,
1090.1326507956582645419985055923462,
1096.4158361028378294577123597264290,
1102.6990214100173943734262138605118,
1108.9822067171969592891400679945946,
1115.2653920243765242048539221286774,
1121.5485773315560891205677762627602,
1127.8317626387356540362816303968430,
1134.1149479459152189519954845309258,
1140.3981332530950112413847818970680,
1146.6813185602745761570986360311508,
1152.9645038674541410728124901652336,
1159.2476891746337059885263442993164,
1165.5308744818132709042401984333992,
1171.8140597889928358199540525674820,
1178.0972450961724007356679067015648,
1184.3804304033519656513817608356476,
1190.6636157105315305670956149697304,
1196.9468010177110954828094691038132,
1203.2299863248906603985233232378960,
1209.5131716320704526879126206040382,
1215.7963569392500176036264747381210,
1222.0795422464295825193403288722038,
1228.3627275536091474350541830062866,
1234.6459128607887123507680371403694,
1240.9290981679682772664818912744522,
1247.2122834751478421821957454085350,
1253.4954687823274070979095995426178,
1259.7786540895069720136234536767006,
1266.0618393966865369293373078107834,
1272.3450247038661018450511619448662,
1278.6282100110458941344404593110085,
1284.9113953182254590501543134450912,
1291.1945806254050239658681675791740,
1297.4777659325845888815820217132568,
1303.7609512397641537972958758473396,
1310.0441365469437187130097299814224,
1316.3273218541232836287235841155052,
1322.6105071613028485444374382495880,
1328.8936924684824134601512923836708,
1335.1768777756619783758651465177536,
1341.4600630828417706652544438838959,
1347.7432483900213355809682980179787,
1354.0264336972009004966821521520615,
1360.3096190043804654123960062861443,
1366.5928043115600303281098604202271,
1372.8759896187395952438237145543098,
1379.1591749259191601595375686883926,
1385.4423602330987250752514228224754,
1391.7255455402782899909652769565582,
1398.0087308474578549066791310906410,
1404.2919161546374198223929852247238,
1410.5751014618172121117822825908661,
1416.8582867689967770274961367249489,
1423.1414720761763419432099908590317,
1429.4246573833559068589238449931145,
1435.7078426905354717746376991271973,
1441.9910279977150366903515532612801,
1448.2742133048946016060654073953629,
1454.5573986120741665217792615294456,
1460.8405839192537314374931156635284,
1467.1237692264332963532069697976112,
1473.4069545336128612689208239316940,
1479.6901398407926535583101212978363,
1485.9733251479722184740239754319191,
1492.2565104551517833897378295660019,
1498.5396957623313483054516837000847,
1504.8228810695109132211655378341675,
1511.1060663766904781368793919682503,
1517.3892516838700430525932461023331,
1523.6724369910496079683071002364159,
1529.9556222982291728840209543704987,
1536.2388076054087377997348085045815,
1542.5219929125885300891241058707237,
1548.8051782197680950048379600048065,
1555.0883635269476599205518141388893,
1561.3715488341272248362656682729721,
1567.6547341413067897519795224070549,
1573.9379194484863546676933765411377,
1580.2211047556659195834072306752205,
1586.5042900628454844991210848093033,
1592.7874753700250494148349389433861,
1599.0706606772046143305487930774689,
1605.3538459843841792462626472115517
};
    double result;
    double Q_func_val;
    const double Eb              = avg_bit_E_to_carrier_A<double>(Ac,Ts,M);
    const double log2M           = std::log2<double>(M); 
    const std::int32_t upper_lim = static_cast<std::int32_t>(std::max<double>(M*0.25,1.0));
    if(__builtin_expect(upper_lim==0,0)) { return (std::numeric_limits<double>::quiet_NaN());}
    const double left_term = 2.0/(std::max<double>(log2M,2.0));
    const double num       = ((Eb+Eb)*log2M)/N0;
    const double invM      = 1.0/M;
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double Q_func_arg= gms::math::cephes_d::sqrt(num);
#else 
    const double Q_func_arg= std::sqrt(num);
#endif 
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
       for(std::int32_t i=1; i<=upper_lim; ++i) 
       {
           const double lut_val = sin_arg_lut[i];
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
           const double sin_val = gms::math::cephes_d::sin(sin_val*invM);
#else 
           const double sin_val = std::sin(sin_val*invM);
#endif 
           Q_func_val  += gms::math::gaussian_Q_approx_chiani(Q_func_arg*sin_val);
       }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
       for(std::int32_t i=1; i<=upper_lim; ++i) 
       {
const double lut_val = sin_arg_lut[i];
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
           const double sin_val = gms::math::cephes_d::sin(sin_val*invM);
#else 
           const double sin_val = std::sin(sin_val*invM);
#endif 
           Q_func_val += gms::math::gaussian_Q_approx_loskot_2T(Q_func_arg*sin_val);
       }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        for(std::int32_t i=1; i<=upper_lim; ++i) 
       {
const double lut_val = sin_arg_lut[i];
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
           const double sin_val = gms::math::cephes_d::sin(sin_val*invM);
#else 
           const double sin_val = std::sin(sin_val*invM);
#endif 
           Q_func_val += gms::math::gaussian_Q_approx_loskot_3T(Q_func_arg*sin_val);
       }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        for(std::int32_t i=1; i<=upper_lim; ++i) 
       {
const double lut_val = sin_arg_lut[i];
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
           const double sin_val = gms::math::cephes_d::sin(sin_val*invM);
#else 
           const double sin_val = std::sin(sin_val*invM);
#endif 
           Q_func_val += gms::math::gaussian_Q_approx_sadhwani_1T(Q_func_arg*sin_val);
       }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
       for(std::int32_t i=1; i<=upper_lim; ++i) 
       {
const double lut_val = sin_arg_lut[i];
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
           const double sin_val = gms::math::cephes_d::sin(sin_val*invM);
#else 
           const double sin_val = std::sin(sin_val*invM);
#endif 
           Q_func_val += gms::math::gaussian_Q_approx_sadhwani_2T(Q_func_arg*sin_val);
       }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        for(std::int32_t i=1; i<=upper_lim; ++i) 
        { 
const double lut_val = sin_arg_lut[i];
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
           const double sin_val = gms::math::cephes_d::sin(sin_val*invM);
#else 
           const double sin_val = std::sin(sin_val*invM);
#endif 
           Q_func_val += gms::math::gaussian_Q_approx_sadhwani_4T(Q_func_arg*sin_val);
        }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_cooper)
    {
        for(std::int32_t i=1; i<=upper_lim; ++i) 
        { 
const double lut_val = sin_arg_lut[i];
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
           const double sin_val = gms::math::cephes_d::sin(sin_val*invM);
#else 
           const double sin_val = std::sin(sin_val*invM);
#endif 
           Q_func_val += gms::math::gaussian_Q_approx_cooper(Q_func_arg*sin_val);
        }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson)
    {
        for(std::int32_t i=1; i<=upper_lim; ++i) 
        { 
const double lut_val = sin_arg_lut[i];
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
           const double sin_val = gms::math::cephes_d::sin(sin_val*invM);
#else 
           const double sin_val = std::sin(sin_val*invM);
#endif 
           Q_func_val += gms::math::gaussian_Q_approx_borjesson(Q_func_arg*sin_val);
        }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        for(std::int32_t i=1; i<=upper_lim; ++i) 
        { 
const double lut_val = sin_arg_lut[i];
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
           const double sin_val = gms::math::cephes_d::sin(sin_val*invM);
#else 
           const double sin_val = std::sin(sin_val*invM);
#endif 
           Q_func_val += gms::math::gaussian_Q_approx_sadhwani_summed(Q_func_arg*sin_val,n);
        }
    }
    result = left_term*Q_func_val;
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
float analytic_BEP_MPSK_AWGN_8_31(const float Ac,const float Ts,
                                   const float M,const float N0,
                                   const std::int32_t n)
{   
    if(__builtin_expect(M>1024.0,0)) { return (-1.0f);}
    constexpr std::int32_t max_M = 256;
constexpr const float sin_arg_lut[max_M] = 
{
3.1415926535897931159979634685442,
9.4247779607693793479938904056326,
15.7079632679489655799898173427209,
21.9911485751285518119857442798093,
28.2743338823081380439816712168977,
34.5575191894877207232639193534851,
40.8407044966673140606872038915753,
47.1238898038468931872557732276618,
53.4070751110264865246790577657521,
59.6902604182060656512476271018386,
65.9734457253856589886709116399288,
72.2566310325652381152394809760153,
78.5398163397448314526627655141056,
84.8230016469244105792313348501921,
91.1061869541040039166546193882823,
97.3893722612835830432231887243688,
103.6725575684631763806464732624590,
109.9557428756427555072150425985456,
116.2389281828223488446383271366358,
122.5221134900019279712068964727223,
128.8052987971815070977754658088088,
135.0884841043611004351987503468990,
141.3716694115406937726220348849893,
147.6548547187202871100453194230795,
153.9380400258998520257591735571623,
160.2212253330794453631824580952525,
166.5044106402590387006057426333427,
172.7875959474386320380290271714330,
179.0707812546181969537428813055158,
185.3539665617977902911661658436060,
191.6371518689773836285894503816962,
197.9203371761569769660127349197865,
204.2035224833365418817265890538692,
210.4867077905161352191498735919595,
216.7698930976957285565731581300497,
223.0530784048753218939964426681399,
229.3362637120548868097102968022227,
235.6194490192344801471335813403130,
241.9026343264140734845568658784032,
248.1858196335936668219801504164934,
254.4690049407732317376940045505762,
260.7521902479528534968267194926739,
267.0353755551324184125405736267567,
273.3185608623119833282544277608395,
279.6017461694916050873871427029371,
285.8849314766711700031009968370199,
292.1681167838507349188148509711027,
298.4513020910303566779475659132004,
304.7344873982099215936614200472832,
311.0176727053895433527941349893808,
317.3008580125691082685079891234636,
323.5840433197486731842218432575464,
329.8672286269282949433545581996441,
336.1504139341078598590684123337269,
342.4335992412874247747822664678097,
348.7167845484670465339149814099073,
354.9999698556466114496288355439901,
361.2831551628262332087615504860878,
367.5663404700057981244754046201706,
373.8495257771853630401892587542534,
380.1327110843649847993219736963511,
386.4158963915445497150358278304338,
392.6990816987241146307496819645166,
398.9822670059037363898823969066143,
405.2654523130833013055962510406971,
411.5486376202629230647289659827948,
417.8318229274424879804428201168776,
424.1150082346220528961566742509604,
430.3981935418016746552893891930580,
436.6813788489812395710032433271408,
442.9645641561608044867170974612236,
449.2477494633404262458498124033213,
455.5309347705199911615636665374041,
461.8141200776996129206963814795017,
468.0973053848791778364102356135845,
474.3804906920587427521240897476673,
480.6636759992383645112568046897650,
486.9468613064179294269706588238478,
493.2300466135974943426845129579306,
499.5132319207771161018172279000282,
505.7964172279566810175310820341110,
512.0796025351362459332449361681938,
518.3627878423158108489587903022766,
524.6459731494954894515103660523891,
530.9291584566750543672242201864719,
537.2123437638546192829380743205547,
543.4955290710341841986519284546375,
549.7787143782137491143657825887203,
556.0618996853934277169173583388329,
562.3450849925729926326312124729156,
568.6282702997525575483450666069984,
574.9114556069321224640589207410812,
581.1946409141116873797727748751640,
587.4778262212913659823243506252766,
593.7610115284709308980382047593594,
600.0441968356504958137520588934422,
606.3273821428300607294659130275249,
612.6105674500096256451797671616077,
618.8937527571891905608936212956905,
625.1769380643688691634451970458031,
631.4601233715484340791590511798859,
637.7433086787279989948729053139687,
644.0264939859075639105867594480515,
650.3096792930871288263006135821342,
656.5928646002668074288521893322468,
662.8760499074463723445660434663296,
669.1592352146259372602798976004124,
675.4424205218055021759937517344952,
681.7256058289850670917076058685780,
688.0087911361647456942591816186905,
694.2919764433443106099730357527733,
700.5751617505238755256868898868561,
706.8583470577034404414007440209389,
713.1415323648830053571145981550217,
719.4247176720625702728284522891045,
725.7079029792422488753800280392170,
731.9910882864218137910938821732998,
738.2742735936013787068077363073826,
744.5574589007809436225215904414654,
750.8406442079605085382354445755482,
757.1238295151401871407870203256607,
763.4070148223197520565008744597435,
769.6902001294993169722147285938263,
775.9733854366788818879285827279091,
782.2565707438584468036424368619919,
788.5397560510381254061940126121044,
794.8229413582176903219078667461872,
801.1061266653972552376217208802700,
807.3893119725768201533355750143528,
813.6724972797563850690494291484356,
819.9556825869359499847632832825184,
826.2388678941156285873148590326309,
832.5220532012951935030287131667137,
838.8052385084747584187425673007965,
845.0884238156543233344564214348793,
851.3716091228338882501702755689621,
857.6547944300135668527218513190746,
863.9379797371931317684357054531574,
870.2211650443726966841495595872402,
876.5043503515522615998634137213230,
882.7875356587318265155772678554058,
889.0707209659115051181288436055183,
895.3539062730910700338426977396011,
901.6370915802706349495565518736839,
907.9202768874501998652704060077667,
914.2034621946297647809842601418495,
920.4866475018093296966981142759323,
926.7698328089890082992496900260448,
933.0530181161685732149635441601276,
939.3362034233481381306773982942104,
945.6193887305277030463912524282932,
951.9025740377072679621051065623760,
958.1857593448869465646566823124886,
964.4689446520665114803705364465714,
970.7521299592460763960843905806541,
977.0353152664256413117982447147369,
983.3185005736052062275120988488197,
989.6016858807848848300636745989323,
995.8848711879644497457775287330151,
1002.1680564951440146614913828670979,
1008.4512418023235795772052370011806,
1014.7344271095031444929190911352634,
1021.0176124166827094086329452693462,
1027.3007977238623880111845210194588,
1033.5839830310419529268983751535416,
1039.8671683382215178426122292876244,
1046.1503536454010827583260834217072,
1052.4335389525806476740399375557899,
1058.7167242597602125897537916898727,
1064.9999095669397775054676458239555,
1071.2830948741193424211814999580383,
1077.5662801812991347105707973241806,
1083.8494654884786996262846514582634,
1090.1326507956582645419985055923462,
1096.4158361028378294577123597264290,
1102.6990214100173943734262138605118,
1108.9822067171969592891400679945946,
1115.2653920243765242048539221286774,
1121.5485773315560891205677762627602,
1127.8317626387356540362816303968430,
1134.1149479459152189519954845309258,
1140.3981332530950112413847818970680,
1146.6813185602745761570986360311508,
1152.9645038674541410728124901652336,
1159.2476891746337059885263442993164,
1165.5308744818132709042401984333992,
1171.8140597889928358199540525674820,
1178.0972450961724007356679067015648,
1184.3804304033519656513817608356476,
1190.6636157105315305670956149697304,
1196.9468010177110954828094691038132,
1203.2299863248906603985233232378960,
1209.5131716320704526879126206040382,
1215.7963569392500176036264747381210,
1222.0795422464295825193403288722038,
1228.3627275536091474350541830062866,
1234.6459128607887123507680371403694,
1240.9290981679682772664818912744522,
1247.2122834751478421821957454085350,
1253.4954687823274070979095995426178,
1259.7786540895069720136234536767006,
1266.0618393966865369293373078107834,
1272.3450247038661018450511619448662,
1278.6282100110458941344404593110085,
1284.9113953182254590501543134450912,
1291.1945806254050239658681675791740,
1297.4777659325845888815820217132568,
1303.7609512397641537972958758473396,
1310.0441365469437187130097299814224,
1316.3273218541232836287235841155052,
1322.6105071613028485444374382495880,
1328.8936924684824134601512923836708,
1335.1768777756619783758651465177536,
1341.4600630828417706652544438838959,
1347.7432483900213355809682980179787,
1354.0264336972009004966821521520615,
1360.3096190043804654123960062861443,
1366.5928043115600303281098604202271,
1372.8759896187395952438237145543098,
1379.1591749259191601595375686883926,
1385.4423602330987250752514228224754,
1391.7255455402782899909652769565582,
1398.0087308474578549066791310906410,
1404.2919161546374198223929852247238,
1410.5751014618172121117822825908661,
1416.8582867689967770274961367249489,
1423.1414720761763419432099908590317,
1429.4246573833559068589238449931145,
1435.7078426905354717746376991271973,
1441.9910279977150366903515532612801,
1448.2742133048946016060654073953629,
1454.5573986120741665217792615294456,
1460.8405839192537314374931156635284,
1467.1237692264332963532069697976112,
1473.4069545336128612689208239316940,
1479.6901398407926535583101212978363,
1485.9733251479722184740239754319191,
1492.2565104551517833897378295660019,
1498.5396957623313483054516837000847,
1504.8228810695109132211655378341675,
1511.1060663766904781368793919682503,
1517.3892516838700430525932461023331,
1523.6724369910496079683071002364159,
1529.9556222982291728840209543704987,
1536.2388076054087377997348085045815,
1542.5219929125885300891241058707237,
1548.8051782197680950048379600048065,
1555.0883635269476599205518141388893,
1561.3715488341272248362656682729721,
1567.6547341413067897519795224070549,
1573.9379194484863546676933765411377,
1580.2211047556659195834072306752205,
1586.5042900628454844991210848093033,
1592.7874753700250494148349389433861,
1599.0706606772046143305487930774689,
1605.3538459843841792462626472115517
};
    float result;
    float Q_func_val;
    const float Eb              = avg_bit_E_to_carrier_A<float>(Ac,Ts,M);
    const float log2M           = std::log2<float>(M); 
    const std::int32_t upper_lim = static_cast<std::int32_t>(std::max<float>(M*0.25,1.0));
    if(__builtin_expect(upper_lim==0,0)) { return (std::numeric_limits<float>::quiet_NaN());}
    const float left_term = 2.0/(std::max<float>(log2M,2.0));
    const float num       = ((Eb+Eb)*log2M)/N0;
    const float invM      = 1.0/M;
    const float Q_func_arg= std::sqrt(num);
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
       for(std::int32_t i=1; i<=upper_lim; ++i) 
       {
           const float lut_val = sin_arg_lut[i];
           const float sin_val = std::sin(sin_val*invM);
           Q_func_val  += gms::math::gaussian_Q_approx_chiani(Q_func_arg*sin_val);
       }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
       for(std::int32_t i=1; i<=upper_lim; ++i) 
       {
           const float lut_val = sin_arg_lut[i];
           const float sin_val = std::sin(sin_val*invM); 
           Q_func_val += gms::math::gaussian_Q_approx_loskot_2T(Q_func_arg*sin_val);
       }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        for(std::int32_t i=1; i<=upper_lim; ++i) 
       {
            const float lut_val = sin_arg_lut[i];
            const float sin_val = std::sin(sin_val*invM);
            Q_func_val += gms::math::gaussian_Q_approx_loskot_3T(Q_func_arg*sin_val);
       }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        for(std::int32_t i=1; i<=upper_lim; ++i) 
       {
            const float lut_val = sin_arg_lut[i];
            const float sin_val = std::sin(sin_val*invM);
            Q_func_val += gms::math::gaussian_Q_approx_sadhwani_1T(Q_func_arg*sin_val);
       }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
       for(std::int32_t i=1; i<=upper_lim; ++i) 
       {
            const float lut_val = sin_arg_lut[i];
            const float sin_val = std::sin(sin_val*invM);
            Q_func_val += gms::math::gaussian_Q_approx_sadhwani_2T(Q_func_arg*sin_val);
       }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        for(std::int32_t i=1; i<=upper_lim; ++i) 
        { 
            const float lut_val = sin_arg_lut[i];
            const float sin_val = std::sin(sin_val*invM); 
            Q_func_val += gms::math::gaussian_Q_approx_sadhwani_4T(Q_func_arg*sin_val);
        }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_cooper)
    {
        for(std::int32_t i=1; i<=upper_lim; ++i) 
        { 
            const float lut_val = sin_arg_lut[i];
            const float sin_val = std::sin(sin_val*invM); 
            Q_func_val += gms::math::gaussian_Q_approx_cooper(Q_func_arg*sin_val);
        }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson)
    {
        for(std::int32_t i=1; i<=upper_lim; ++i) 
        { 
            const float lut_val = sin_arg_lut[i];
            const float sin_val = std::sin(sin_val*invM); 
            Q_func_val += gms::math::gaussian_Q_approx_borjesson(Q_func_arg*sin_val);
        }
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        for(std::int32_t i=1; i<=upper_lim; ++i) 
        { 
            const float lut_val = sin_arg_lut[i];
            const float sin_val = std::sin(sin_val*invM); 
            Q_func_val += gms::math::gaussian_Q_approx_sadhwani_summed(Q_func_arg*sin_val,n);
        }
    }
    result = left_term*Q_func_val;
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
double analytic_BEP_BPSK_8_47(const double Ac,const double Ts,
                              const double M,const double N0,
                              const double Phic,const std::int32_t n)
{
    double result;
    double Q_func_val;
    const double Eb         = avg_bit_E_to_carrier_A<double>(Ac,Ts,M);
    const double sqrt_arg   = (Eb+Eb)/N0;
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double cos_phic   = gms::math::cephes_d::cos(Phic);
    const double Q_func_arg = gms::math::cephes_d::sqrt(sqrt_arg);
#else 
    const double cos_phic   = std::cos(Phic);
    const double Q_func_arg = std::sqrt(sqrt_arg);
#endif 
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
       Q_func_val  = gms::math::gaussian_Q_approx_chiani(Q_func_arg*cos_phic);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T(Q_func_arg*cos_phic);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T(Q_func_arg*cos_phic);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T(Q_func_arg*cos_phic);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T(Q_func_arg*cos_phic);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T(Q_func_arg*cos_phic);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_cooper)
    {
        Q_func_val = gms::math::gaussian_Q_approx_cooper(Q_func_arg*cos_phic);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson)
    {
        Q_func_val = gms::math::gaussian_Q_approx_borjesson(Q_func_arg*cos_phic);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed(Q_func_arg*cos_phic,n);
    }
    result = Q_func_val;
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
float analytic_BEP_BPSK_8_47(const float Ac,const float Ts,
                              const float M,const float N0,
                              const float Phic,const std::int32_t n)
{
    float result;
    float Q_func_val;
    const float Eb         = avg_bit_E_to_carrier_A<float>(Ac,Ts,M);
    const float sqrt_arg   = (Eb+Eb)/N0;
    const float cos_phic   = std::cos(Phic);
    const float Q_func_arg = std::sqrt(sqrt_arg);
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
       Q_func_val  = gms::math::gaussian_Q_approx_chiani(Q_func_arg*cos_phic);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_2T(Q_func_arg*cos_phic);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_loskot_3T(Q_func_arg*cos_phic);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_1T(Q_func_arg*cos_phic);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_2T(Q_func_arg*cos_phic);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_4T(Q_func_arg*cos_phic);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_cooper)
    {
        Q_func_val = gms::math::gaussian_Q_approx_cooper(Q_func_arg*cos_phic);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson)
    {
        Q_func_val = gms::math::gaussian_Q_approx_borjesson(Q_func_arg*cos_phic);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = gms::math::gaussian_Q_approx_sadhwani_summed(Q_func_arg*cos_phic,n);
    }
    result = Q_func_val;
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
double analytic_BEP_QPSK_8_56(const double Ac,const double Ts,
                              const double M,const double N0,
                              const double Phic,const std::int32_t n)
{
    double result;
    double Q_func_val;
    const double Eb = avg_bit_E_to_carrier_A<double>(Ac,Ts,M);
    double sinphic;
    double cosphic;
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    __attribute__((unused)) std::int32_t ret_val = gms::math::cephes_d::sincos(Phic,&sinphic,&cosphic,0);
    const double sqrt_val = gms::math::cephes_d::sqrt((Eb+Eb)/N0);
#else 
    sinphic = std::sin(phic);
    cosphic = std::cos(phic);
    const double sqrt_val = std::sqrt((Eb+Eb)/N0); 
#endif 
    const double Q_func_arg1 = sqrt_val*(cosphic-sinphic);
    const double Q_func_arg2 = sqrt_val*(cosphic+sinphic);
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
        const double tmp1  = 0.5*gms::math::gaussian_Q_approx_chiani(Q_func_arg1);
        const double tmp2  = 0.5*gms::math::gaussian_Q_approx_chiani(Q_func_arg2);
        Q_func_val         = tmp1+tmp2;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        const double tmp1 = 0.5*gms::math::gaussian_Q_approx_loskot_2T(Q_func_arg1);
        const double tmp2 = 0.5*gms::math::gaussian_Q_approx_loskot_2T(Q_func_arg2);
        Q_func_val        = tmp1+tmp2;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        const double tmp1 = 0.5*gms::math::gaussian_Q_approx_loskot_3T(Q_func_arg1);
        const double tmp2 = 0.5*gms::math::gaussian_Q_approx_loskot_3T(Q_func_arg2);
        Q_func_val        = tmp1+tmp2;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        const double tmp1 = gms::math::gaussian_Q_approx_sadhwani_1T(Q_func_arg1);
        const double tmp2 = gms::math::gaussian_Q_approx_sadhwani_1T(Q_func_arg2);
        Q_func_val        = tmp1+tmp2;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        const double tmp1 = gms::math::gaussian_Q_approx_sadhwani_2T(Q_func_arg1);
        const double tmp2 = gms::math::gaussian_Q_approx_sadhwani_2T(Q_func_arg2);
        Q_func_val        = tmp1+tmp2;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        const double tmp1 = gms::math::gaussian_Q_approx_sadhwani_4T(Q_func_arg1);
        const double tmp2 = gms::math::gaussian_Q_approx_sadhwani_4T(Q_func_arg2);
        Q_func_val        = tmp1+tmp2;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_cooper)
    {
        const double tmp1 = gms::math::gaussian_Q_approx_cooper(Q_func_arg1);
        const double tmp2 = gms::math::gaussian_Q_approx_cooper(Q_func_arg2);
        Q_func_val        = tmp1+tmp2;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson)
    {
        const double tmp1 = gms::math::gaussian_Q_approx_borjesson(Q_func_arg1);
        const double tmp2 = gms::math::gaussian_Q_approx_borjesson(Q_func_arg2);
        Q_func_val        = tmp1+tmp2;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        const double tmp1 = gms::math::gaussian_Q_approx_sadhwani_summed(Q_func_arg1,n);
        const double tmp2 = gms::math::gaussian_Q_approx_sadhwani_summed(Q_func_arg2,n);
        Q_func_val        = tmp1+tmp2;
    }
    result = Q_func_val;
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
float analytic_BEP_QPSK_8_56(const float Ac,const float Ts,
                              const float M,const float N0,
                              const float Phic,const std::int32_t n)
{
    float result;
    float Q_func_val;
    const float Eb = avg_bit_E_to_carrier_A<float>(Ac,Ts,M);
    float sinphic;
    float cosphic;
    sinphic = std::sin(phic);
    cosphic = std::cos(phic);
    const float sqrt_val = std::sqrt((Eb+Eb)/N0);  
    const float Q_func_arg1 = sqrt_val*(cosphic-sinphic);
    const float Q_func_arg2 = sqrt_val*(cosphic+sinphic);
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
        const float tmp1  = 0.5*gms::math::gaussian_Q_approx_chiani(Q_func_arg1);
        const float tmp2  = 0.5*gms::math::gaussian_Q_approx_chiani(Q_func_arg2);
        Q_func_val         = tmp1+tmp2;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        const float tmp1 = 0.5*gms::math::gaussian_Q_approx_loskot_2T(Q_func_arg1);
        const float tmp2 = 0.5*gms::math::gaussian_Q_approx_loskot_2T(Q_func_arg2);
        Q_func_val        = tmp1+tmp2;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        const float tmp1 = 0.5*gms::math::gaussian_Q_approx_loskot_3T(Q_func_arg1);
        const float tmp2 = 0.5*gms::math::gaussian_Q_approx_loskot_3T(Q_func_arg2);
        Q_func_val        = tmp1+tmp2;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        const float tmp1 = gms::math::gaussian_Q_approx_sadhwani_1T(Q_func_arg1);
        const float tmp2 = gms::math::gaussian_Q_approx_sadhwani_1T(Q_func_arg2);
        Q_func_val        = tmp1+tmp2;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        const float tmp1 = gms::math::gaussian_Q_approx_sadhwani_2T(Q_func_arg1);
        const float tmp2 = gms::math::gaussian_Q_approx_sadhwani_2T(Q_func_arg2);
        Q_func_val        = tmp1+tmp2;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        const float tmp1 = gms::math::gaussian_Q_approx_sadhwani_4T(Q_func_arg1);
        const float tmp2 = gms::math::gaussian_Q_approx_sadhwani_4T(Q_func_arg2);
        Q_func_val        = tmp1+tmp2;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_cooper)
    {
        const float tmp1 = gms::math::gaussian_Q_approx_cooper(Q_func_arg1);
        const float tmp2 = gms::math::gaussian_Q_approx_cooper(Q_func_arg2);
        Q_func_val        = tmp1+tmp2;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson)
    {
        const float tmp1 = gms::math::gaussian_Q_approx_borjesson(Q_func_arg1);
        const float tmp2 = gms::math::gaussian_Q_approx_borjesson(Q_func_arg2);
        Q_func_val        = tmp1+tmp2;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        const float tmp1 = gms::math::gaussian_Q_approx_sadhwani_summed(Q_func_arg1,n);
        const float tmp2 = gms::math::gaussian_Q_approx_sadhwani_summed(Q_func_arg2,n);
        Q_func_val        = tmp1+tmp2;
    }
    result = Q_func_val;
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
double analytic_BEP_FSK2_noncoh_8_68(const double Ac,const double Ts,
                                     const double M,const double N0)
{
    const double Eb = avg_bit_E_to_carrier_A<double>(Ac,Ts,M);
    const double exp_arg = Eb/(N0+N+);
    return (0.5*gms::math::expapprox_d(-exp_arg));
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
float analytic_BEP_FSK2_noncoh_8_68(const float Ac,const float Ts,
                                     const float M,const float N0)
{
    const float Eb = avg_bit_E_to_carrier_A<float>(Ac,Ts,M);
    const float exp_arg = Eb/(N0+N+);
    return (0.5f*gms::math::expapprox(-exp_arg));
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
double analytic_BEP_BPSK_8_60(const double Ac,const double Ts,
                              const double M, const double N0,
                              const double Bl,const double Tb,
                              const double mu,const std::int32_t n)  
{
    double result;
    double Q1_func_val;
    const double a = BPSK_param_a_8_61<double>(Ac,Ts,M,N0,Bl,Tb);
    const double b = BPSK_param_b_8_61<double>(Ac,Ts,M,N0,Bl,Tb);
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double a_sqrt = gms::math::cephes_d::sqrt(a);
    const double b_sqrt = gms::math::cephes_d::sqrt(b);
#else 
    const double a_sqrt = std::sqrt(a);
    const double b_sqrt = std::sqrt(b);
#endif 
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
        const double tmp1  = 1.0-gms::math::marcum_Q_approx_chiani(mu,b_sqrt,a_sqrt);
        const double tmp2  = gms::math::marcum_Q_approx_chiani(mu,a_sqrt,b_sqrt);
        Q1_func_val        = 0.5*(tmp1+tmp2);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        const double tmp1 = 1.0-gms::math::marcum_Q_approx_loskot_2T(mu,b_sqrt,a_sqrt);
        const double tmp2 = gms::math::marcum_Q_approx_loskot_2T(mu,a_sqrt,b_sqrt);
        Q1_func_val       = 0.5*(tmp1+tmp2);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        const double tmp1 = 1.0-gms::math::marcum_Q_approx_loskot_3T(mu,b_sqrt,a_sqrt);
        const double tmp2 = gms::math::marcum_Q_approx_loskot_3T(mu,a_sqrt,b_sqrt);
        Q1_func_val       = 0.5*(tmp1+tmp2);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        const double tmp1 = 1.0-gms::math::marcum_Q_approx_sadhwani_1T(mu,b_sqrt,a_sqrt);
        const double tmp2 = gms::math::marcum_Q_approx_sadhwani_1T(mu,a_sqrt,b_sqrt);
        Q1_func_val       = 0.5*(tmp1+tmp2);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        const double tmp1 = 1.0-gms::math::marcum_Q_approx_sadhwani_2T(mu,b_sqrt,a_sqrt);
        const double tmp2 = gms::math::marcum_Q_approx_sadhwani_2T(mu,a_sqrt,b_sqrt);
        Q1_func_val       = 0.5*(tmp1+tmp2);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        const double tmp1 = 1.0-gms::math::marcum_Q_approx_sadhwani_4T(mu,b_sqrt,a_sqrt);
        const double tmp2 = gms::math::marcum_Q_approx_sadhwani_4T(mu,a_sqrt,b_sqrt);
        Q1_func_val        = 0.5*(tmp1+tmp2);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_cooper)
    {
        const double tmp1 = 1.0-gms::math::marcum_Q_approx_cooper(mu,b_sqrt,a_sqrt);
        const double tmp2 = gms::math::marcum_Q_approx_cooper(mu,a_sqrt,b_sqrt);
        Q1_func_val        = 0.5*(tmp1+tmp2);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson)
    {
        const double tmp1 = 1.0-gms::math::marcum_Q_approx_borjesson(mu,b_sqrt,a_sqrt);
        const double tmp2 = gms::math::marcum_Q_approx_borjesson(mu,a_sqrt,b_sqrt);
        Q1_func_val        = 0.5*(tmp1+tmp2);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        const double tmp1 = 1.0-gms::math::marcum_Q_approx_sadhwani_summed(mu,b_sqrt,a_sqrt,n);
        const double tmp2 = gms::math::marcum_Q_approx_sadhwani_summed(mu,a_sqrt,b_sqrt,n);
        Q1_func_val        = 0.5*(tmp1+tmp2);
    }
    result = Q1_func_val;
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
float analytic_BEP_BPSK_8_60(const float Ac,const float Ts,
                              const float M, const float N0,
                              const float Bl,const float Tb,
                              const float mu,const std::int32_t n)  
{
    float result;
    float Q1_func_val;
    const float a = BPSK_param_a_8_61<float>(Ac,Ts,M,N0,Bl,Tb);
    const float b = BPSK_param_b_8_61<float>(Ac,Ts,M,N0,Bl,Tb);
    const float a_sqrt = std::sqrt(a);
    const float b_sqrt = std::sqrt(b); 
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
        const float tmp1  = 1.0f-gms::math::marcum_Q_approx_chiani(mu,b_sqrt,a_sqrt);
        const float tmp2  = gms::math::marcum_Q_approx_chiani(mu,a_sqrt,b_sqrt);
        Q1_func_val        = 0.5f*(tmp1+tmp2);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        const float tmp1 = 1.0f-gms::math::marcum_Q_approx_loskot_2T(mu,b_sqrt,a_sqrt);
        const float tmp2 = gms::math::marcum_Q_approx_loskot_2T(mu,a_sqrt,b_sqrt);
        Q1_func_val       = 0.5f*(tmp1+tmp2);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        const float tmp1 = 1.0f-gms::math::marcum_Q_approx_loskot_3T(mu,b_sqrt,a_sqrt);
        const float tmp2 = gms::math::marcum_Q_approx_loskot_3T(mu,a_sqrt,b_sqrt);
        Q1_func_val       = 0.5f*(tmp1+tmp2);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        const float tmp1 = 1.0f-gms::math::marcum_Q_approx_sadhwani_1T(mu,b_sqrt,a_sqrt);
        const float tmp2 = gms::math::marcum_Q_approx_sadhwani_1T(mu,a_sqrt,b_sqrt);
        Q1_func_val       = 0.5f*(tmp1+tmp2);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        const float tmp1 = 1.0f-gms::math::marcum_Q_approx_sadhwani_2T(mu,b_sqrt,a_sqrt);
        const float tmp2 = gms::math::marcum_Q_approx_sadhwani_2T(mu,a_sqrt,b_sqrt);
        Q1_func_val       = 0.5f*(tmp1+tmp2);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        const float tmp1 = 1.0f-gms::math::marcum_Q_approx_sadhwani_4T(mu,b_sqrt,a_sqrt);
        const float tmp2 = gms::math::marcum_Q_approx_sadhwani_4T(mu,a_sqrt,b_sqrt);
        Q1_func_val        = 0.5f*(tmp1+tmp2);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_cooper)
    {
        const float tmp1 = 1.0f-gms::math::marcum_Q_approx_cooper(mu,b_sqrt,a_sqrt);
        const float tmp2 = gms::math::marcum_Q_approx_cooper(mu,a_sqrt,b_sqrt);
        Q1_func_val        = 0.5f*(tmp1+tmp2);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson)
    {
        const float tmp1 = 1.0f-gms::math::marcum_Q_approx_borjesson(mu,b_sqrt,a_sqrt);
        const float tmp2 = gms::math::marcum_Q_approx_borjesson(mu,a_sqrt,b_sqrt);
        Q1_func_val        = 0.5f*(tmp1+tmp2);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        const float tmp1 = 1.0f-gms::math::marcum_Q_approx_sadhwani_summed(mu,b_sqrt,a_sqrt,n);
        const float tmp2 = gms::math::marcum_Q_approx_sadhwani_summed(mu,a_sqrt,b_sqrt,n);
        Q1_func_val        = 0.5f*(tmp1+tmp2);
    }
    result = Q1_func_val;
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
double analytic_avg_BEP_QPSK_MSK_8_63(const double Ac,const double Ts,
                                      const double M, const double N0,
                                      const double Bl,const double Tb,
                                      const double mu,const std::int32_t n)
{
    // Marcum-Q order=1, is only allowed
    if(__builtin_expect(mu!=1.0,0)) { return (std::numeric_limits<double>::quiet_NaN());}
    double result;
    double Q1_func_val;
    double tmp1;
    double tmp2;
    double tmp3;
    double tmp4;
    double left;
    double right;
    const double a1q = OQPSK_param_a1_8_64<double>(Ac,Ts,M,N0,Bl,Tb);
    const double b1q = OQPSK_param_b1_8_64<double>(Ac,Ts,M,N0,Bl,Tb);
    const double a2q = OQPSK_param_a2_8_64<double>(Ac,Ts,M,N0,Bl,Tb);
    const double b2q = OQPSK_param_b2_8_64<double>(Ac,Ts,M,N0,Bl,Tb);
    const double a1m = MSK_param_a1_8_65<double>(Ac,Ts,M,N0,Bl,Tb);
    const double b1m = MSK_param_b1_8_65<double>(Ac,Ts,M,N0,Bl,Tb);
    const double a2m = MSK_param_a2_8_65<double>(Ac,Ts,M,N0,Bl,Tb);
    const double b2m = MSK_param_b2_8_65<double>(Ac,Ts,M,N0,Bl,Tb);
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double a1q_sqrt = gms::math::cephes_d::sqrt(a1q);
    const double b1q_sqrt = gms::math::cephes_d::sqrt(b1q);
    const double a2q_sqrt = gms::math::cephes_d::sqrt(a2q);
    const double b2q_sqrt = gms::math::cephes_d::sqrt(b2q);
    const double a1m_sqrt = gms::math::cephes_d::sqrt(a1m);
    const double b1m_sqrt = gms::math::cephes_d::sqrt(b1m);
    const double a2m_sqrt = gms::math::cephes_d::sqrt(a2m);
    const double b2m_sqrt = gms::math::cephes_d::sqrt(b2m);
#else 
    const double a1q_sqrt = std::sqrt(a1q);
    const double b1q_sqrt = std::sqrt(b1q);
    const double a2q_sqrt = std::sqrt(a2q);
    const double b2q_sqrt = std::sqrt(b2q);
    const double a1m_sqrt = std::sqrt(a1m);
    const double b1m_sqrt = std::sqrt(b1m);
    const double a2m_sqrt = std::sqrt(a2m);
    const double b2m_sqrt = std::sqrt(b2m);
#endif 
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
        tmp1 = 1.0-gms::math::marcum_Q_approx_chiani(mu,b1q_sqrt,a1q_sqrt);
        tmp2 = gms::math::marcum_Q_approx_chiani(mu,a1m_sqrt,b1m_sqrt);
        tmp3 = 1.0-gms::math::marcum_Q_approx_chiani(mu,b2q_sqrt,a2q_sqrt);
        tmp4 = gms::math::marcum_Q_approx_chiani(mu,a2m_sqrt,b2m_sqrt);
        left = 0.25*(tmp1+tmp2);
        right= 0.25*(tmp3+tmp4);
        Q1_func_val       = left+right;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        tmp1 = 1.0-gms::math::marcum_Q_approx_loskot_2T(mu,b1q_sqrt,a1q_sqrt);
        tmp2 = gms::math::marcum_Q_approx_loskot_2T(mu,a1m_sqrt,b1m_sqrt);
        tmp3 = 1.0-gms::math::marcum_Q_approx_loskot_2T(mu,b2q_sqrt,a2q_sqrt);
        tmp4 = gms::math::marcum_Q_approx_loskot_2T(mu,a2m_sqrt,b2m_sqrt);
        left = 0.25*(tmp1+tmp2);
        right= 0.25*(tmp3+tmp4);
        Q1_func_val       = left+right;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        tmp1 = 1.0-gms::math::marcum_Q_approx_loskot_3T(mu,b1q_sqrt,a1q_sqrt);
        tmp2 = gms::math::marcum_Q_approx_loskot_3T(mu,a1m_sqrt,b1m_sqrt);
        tmp3 = 1.0-gms::math::marcum_Q_approx_loskot_3T(mu,b2q_sqrt,a2q_sqrt);
        tmp4 = gms::math::marcum_Q_approx_loskot_3T(mu,a2m_sqrt,b2m_sqrt);
        left = 0.25*(tmp1+tmp2);
        right= 0.25*(tmp3+tmp4);
        Q1_func_val       = left+right;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        tmp1 = 1.0-gms::math::marcum_Q_approx_sadhwani_1T(mu,b1q_sqrt,a1q_sqrt);
        tmp2 = gms::math::marcum_Q_approx_sadhwani_1T(mu,a1m_sqrt,b1m_sqrt);
        tmp3 = 1.0-gms::math::marcum_Q_approx_sadhwani_1T(mu,b2q_sqrt,a2q_sqrt);
        tmp4 = gms::math::marcum_Q_approx_sadhwani_1T(mu,a2m_sqrt,b2m_sqrt);
        left = 0.25*(tmp1+tmp2);
        right= 0.25*(tmp3+tmp4);
        Q1_func_val       = left+right;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        tmp1 = 1.0-gms::math::marcum_Q_approx_sadhwani_2T(mu,b1q_sqrt,a1q_sqrt);
        tmp2 = gms::math::marcum_Q_approx_sadhwani_2T(mu,a1m_sqrt,b1m_sqrt);
        tmp3 = 1.0-gms::math::marcum_Q_approx_sadhwani_2T(mu,b2q_sqrt,a2q_sqrt);
        tmp4 = gms::math::marcum_Q_approx_sadhwani_2T(mu,a2m_sqrt,b2m_sqrt);
        left = 0.25*(tmp1+tmp2);
        right= 0.25*(tmp3+tmp4);
        Q1_func_val       = left+right;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        tmp1 = 1.0-gms::math::marcum_Q_approx_sadhwani_4T(mu,b1q_sqrt,a1q_sqrt);
        tmp2 = gms::math::marcum_Q_approx_sadhwani_4T(mu,a1m_sqrt,b1m_sqrt);
        tmp3 = 1.0-gms::math::marcum_Q_approx_sadhwani_4T(mu,b2q_sqrt,a2q_sqrt);
        tmp4 = gms::math::marcum_Q_approx_sadhwani_4T(mu,a2m_sqrt,b2m_sqrt);
        left = 0.25*(tmp1+tmp2);
        right= 0.25*(tmp3+tmp4);
        Q1_func_val       = left+right;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_cooper)
    {
        tmp1 = 1.0-gms::math::marcum_Q_approx_cooper(mu,b1q_sqrt,a1q_sqrt);
        tmp2 = gms::math::marcum_Q_approx_cooper(mu,a1m_sqrt,b1m_sqrt);
        tmp3 = 1.0-gms::math::marcum_Q_approx_cooper(mu,b2q_sqrt,a2q_sqrt);
        tmp4 = gms::math::marcum_Q_approx_cooper(mu,a2m_sqrt,b2m_sqrt);
        left = 0.25*(tmp1+tmp2);
        right= 0.25*(tmp3+tmp4);
        Q1_func_val       = left+right;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson)
    {
        tmp1 = 1.0-gms::math::marcum_Q_approx_borjesson(mu,b1q_sqrt,a1q_sqrt);
        tmp2 = gms::math::marcum_Q_approx_borjesson(mu,a1m_sqrt,b1m_sqrt);
        tmp3 = 1.0-gms::math::marcum_Q_approx_borjesson(mu,b2q_sqrt,a2q_sqrt);
        tmp4 = gms::math::marcum_Q_approx_borjesson(mu,a2m_sqrt,b2m_sqrt);
        left = 0.25*(tmp1+tmp2);
        right= 0.25*(tmp3+tmp4);
        Q1_func_val       = left+right;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        tmp1 = 1.0-gms::math::marcum_Q_approx_sadhwani_summed(mu,b1q_sqrt,a1q_sqrt,n);
        tmp2 = gms::math::marcum_Q_approx_sadhwani_summed(mu,a1m_sqrt,b1m_sqrt,n);
        tmp3 = 1.0-gms::math::marcum_Q_approx_sadhwani_summed(mu,b2q_sqrt,a2q_sqrt,n);
        tmp4 = gms::math::marcum_Q_approx_sadhwani_summed(mu,a2m_sqrt,b2m_sqrt,n);
        left = 0.25*(tmp1+tmp2);
        right= 0.25*(tmp3+tmp4);
        Q1_func_val       = left+right;
    }
    result = Q1_func_val;
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
float analytic_avg_BEP_QPSK_MSK_8_63(const float Ac,const float Ts,
                                      const float M, const float N0,
                                      const float Bl,const float Tb,
                                      const float mu,const std::int32_t n)
{
    // Marcum-Q order=1, is only allowed
    if(__builtin_expect(mu!=1.0,0)) { return (std::numeric_limits<float>::quiet_NaN());}
    float result;
    float Q1_func_val;
    float tmp1;
    float tmp2;
    float tmp3;
    float tmp4;
    float left;
    float right;
    const float a1q = OQPSK_param_a1_8_64<float>(Ac,Ts,M,N0,Bl,Tb);
    const float b1q = OQPSK_param_b1_8_64<float>(Ac,Ts,M,N0,Bl,Tb);
    const float a2q = OQPSK_param_a2_8_64<float>(Ac,Ts,M,N0,Bl,Tb);
    const float b2q = OQPSK_param_b2_8_64<float>(Ac,Ts,M,N0,Bl,Tb);
    const float a1m = MSK_param_a1_8_65<float>(Ac,Ts,M,N0,Bl,Tb);
    const float b1m = MSK_param_b1_8_65<float>(Ac,Ts,M,N0,Bl,Tb);
    const float a2m = MSK_param_a2_8_65<float>(Ac,Ts,M,N0,Bl,Tb);
    const float b2m = MSK_param_b2_8_65<float>(Ac,Ts,M,N0,Bl,Tb);
    const float a1q_sqrt = std::sqrt(a1q);
    const float b1q_sqrt = std::sqrt(b1q);
    const float a2q_sqrt = std::sqrt(a2q);
    const float b2q_sqrt = std::sqrt(b2q);
    const float a1m_sqrt = std::sqrt(a1m);
    const float b1m_sqrt = std::sqrt(b1m);
    const float a2m_sqrt = std::sqrt(a2m);
    const float b2m_sqrt = std::sqrt(b2m);
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
        tmp1 = 1.0-gms::math::marcum_Q_approx_chiani(mu,b1q_sqrt,a1q_sqrt);
        tmp2 = gms::math::marcum_Q_approx_chiani(mu,a1m_sqrt,b1m_sqrt);
        tmp3 = 1.0-gms::math::marcum_Q_approx_chiani(mu,b2q_sqrt,a2q_sqrt);
        tmp4 = gms::math::marcum_Q_approx_chiani(mu,a2m_sqrt,b2m_sqrt);
        left = 0.25*(tmp1+tmp2);
        right= 0.25*(tmp3+tmp4);
        Q1_func_val       = left+right;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        tmp1 = 1.0-gms::math::marcum_Q_approx_loskot_2T(mu,b1q_sqrt,a1q_sqrt);
        tmp2 = gms::math::marcum_Q_approx_loskot_2T(mu,a1m_sqrt,b1m_sqrt);
        tmp3 = 1.0-gms::math::marcum_Q_approx_loskot_2T(mu,b2q_sqrt,a2q_sqrt);
        tmp4 = gms::math::marcum_Q_approx_loskot_2T(mu,a2m_sqrt,b2m_sqrt);
        left = 0.25*(tmp1+tmp2);
        right= 0.25*(tmp3+tmp4);
        Q1_func_val       = left+right;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        tmp1 = 1.0-gms::math::marcum_Q_approx_loskot_3T(mu,b1q_sqrt,a1q_sqrt);
        tmp2 = gms::math::marcum_Q_approx_loskot_3T(mu,a1m_sqrt,b1m_sqrt);
        tmp3 = 1.0-gms::math::marcum_Q_approx_loskot_3T(mu,b2q_sqrt,a2q_sqrt);
        tmp4 = gms::math::marcum_Q_approx_loskot_3T(mu,a2m_sqrt,b2m_sqrt);
        left = 0.25*(tmp1+tmp2);
        right= 0.25*(tmp3+tmp4);
        Q1_func_val       = left+right;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        tmp1 = 1.0-gms::math::marcum_Q_approx_sadhwani_1T(mu,b1q_sqrt,a1q_sqrt);
        tmp2 = gms::math::marcum_Q_approx_sadhwani_1T(mu,a1m_sqrt,b1m_sqrt);
        tmp3 = 1.0-gms::math::marcum_Q_approx_sadhwani_1T(mu,b2q_sqrt,a2q_sqrt);
        tmp4 = gms::math::marcum_Q_approx_sadhwani_1T(mu,a2m_sqrt,b2m_sqrt);
        left = 0.25*(tmp1+tmp2);
        right= 0.25*(tmp3+tmp4);
        Q1_func_val       = left+right;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        tmp1 = 1.0-gms::math::marcum_Q_approx_sadhwani_2T(mu,b1q_sqrt,a1q_sqrt);
        tmp2 = gms::math::marcum_Q_approx_sadhwani_2T(mu,a1m_sqrt,b1m_sqrt);
        tmp3 = 1.0-gms::math::marcum_Q_approx_sadhwani_2T(mu,b2q_sqrt,a2q_sqrt);
        tmp4 = gms::math::marcum_Q_approx_sadhwani_2T(mu,a2m_sqrt,b2m_sqrt);
        left = 0.25*(tmp1+tmp2);
        right= 0.25*(tmp3+tmp4);
        Q1_func_val       = left+right;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        tmp1 = 1.0-gms::math::marcum_Q_approx_sadhwani_4T(mu,b1q_sqrt,a1q_sqrt);
        tmp2 = gms::math::marcum_Q_approx_sadhwani_4T(mu,a1m_sqrt,b1m_sqrt);
        tmp3 = 1.0-gms::math::marcum_Q_approx_sadhwani_4T(mu,b2q_sqrt,a2q_sqrt);
        tmp4 = gms::math::marcum_Q_approx_sadhwani_4T(mu,a2m_sqrt,b2m_sqrt);
        left = 0.25*(tmp1+tmp2);
        right= 0.25*(tmp3+tmp4);
        Q1_func_val       = left+right;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_cooper)
    {
        tmp1 = 1.0-gms::math::marcum_Q_approx_cooper(mu,b1q_sqrt,a1q_sqrt);
        tmp2 = gms::math::marcum_Q_approx_cooper(mu,a1m_sqrt,b1m_sqrt);
        tmp3 = 1.0-gms::math::marcum_Q_approx_cooper(mu,b2q_sqrt,a2q_sqrt);
        tmp4 = gms::math::marcum_Q_approx_cooper(mu,a2m_sqrt,b2m_sqrt);
        left = 0.25*(tmp1+tmp2);
        right= 0.25*(tmp3+tmp4);
        Q1_func_val       = left+right;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson)
    {
        tmp1 = 1.0-gms::math::marcum_Q_approx_borjesson(mu,b1q_sqrt,a1q_sqrt);
        tmp2 = gms::math::marcum_Q_approx_borjesson(mu,a1m_sqrt,b1m_sqrt);
        tmp3 = 1.0-gms::math::marcum_Q_approx_borjesson(mu,b2q_sqrt,a2q_sqrt);
        tmp4 = gms::math::marcum_Q_approx_borjesson(mu,a2m_sqrt,b2m_sqrt);
        left = 0.25*(tmp1+tmp2);
        right= 0.25*(tmp3+tmp4);
        Q1_func_val       = left+right;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        tmp1 = 1.0-gms::math::marcum_Q_approx_sadhwani_summed(mu,b1q_sqrt,a1q_sqrt,n);
        tmp2 = gms::math::marcum_Q_approx_sadhwani_summed(mu,a1m_sqrt,b1m_sqrt,n);
        tmp3 = 1.0-gms::math::marcum_Q_approx_sadhwani_summed(mu,b2q_sqrt,a2q_sqrt,n);
        tmp4 = gms::math::marcum_Q_approx_sadhwani_summed(mu,a2m_sqrt,b2m_sqrt,n);
        left = 0.25*(tmp1+tmp2);
        right= 0.25*(tmp3+tmp4);
        Q1_func_val       = left+right;
    }
    result = Q1_func_val;
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
double analytic_BEP_MFSK2_nonorthog_8_69(const double Ac,const double Ts,
                                         const double M, const double N0,
                                         const double rho,const double mu,
                                         const std::int32_t n)
{
    if(__builtin_expect(mu!=1.0,0)) { return (std::numeric_limits<double>::quiet_NaN());}
    const double Eb = avg_bit_E_to_carrier_A<double>(Ac,Ts,M);
    const double Eb_div_2N0 = Eb/(N0+N0);
    if(__builtin_expect(rho==0.0,0))
    {
       return (0.5*gms::math::expapprox_d(-Eb_div_2N0));
    }
    double result;
    double Q1_func_val;
    double a;
    double b;
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double rho_term = gms::math::cephes_d::sqrt(1.0-rho*rho);
#else
    const double rho_term = std::sqrt(1.0-rho*rho);
#endif 
    a = Eb_div_2N0*(1.0-rho_term);
    b = Eb_div_2N0*(1.0+rho_term);
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double bessI0_arg = gms::math::cephes_d::sqrt(a*b);
#else 
    const double bessI0_arg = std::sqrt(a*b);
#endif 
    const double exp_val    = 0.5*gms::math::expapprox_d(0.5*(a+b));
    const double bessI0_val = std::cyl_bessel_i<double>(0,bessI0_arg);
    const double right_fac  = exp_val*bessI0_val;
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double a_sqrt      = gms::math::cephes_d::sqrt(a);
    const double b_sqrt      = gms::math::cephes_d::sqrt(b);
#else 
    const double a_sqrt      = std::sqrt(a);
    const double b_sqrt      = std::sqrt(b);
#endif 
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
        Q1_func_val = gms::math::marcum_Q_approx_chiani(mu,a_sqrt,b_sqrt)-right_fac;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        Q1_func_val = gms::math::marcum_Q_approx_loskot_2T(mu,a_sqrt,b_sqrt)-right_fac;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        Q1_func_val = gms::math::marcum_Q_approx_loskot_3T(mu,a_sqrt,b_sqrt)-right_fac;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q1_func_val = gms::math::marcum_Q_approx_sadhwani_1T(mu,a_sqrt,b_sqrt)-right_fac;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q1_func_val = gms::math::marcum_Q_approx_sadhwani_2T(mu,a_sqrt,b_sqrt)-right_fac;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q1_func_val = gms::math::marcum_Q_approx_sadhwani_4T(mu,a_sqrt,b_sqrt)-right_fac;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_cooper)
    {
        Q1_func_val = gms::math::marcum_Q_approx_cooper(mu,a_sqrt,b_sqrt)-right_fac;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson)
    {
        Q1_func_val = gms::math::marcum_Q_approx_borjesson(mu,a_sqrt,b_sqrt)-right_fac;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q1_func_val = gms::math::marcum_Q_approx_sadhwani_summed(mu,a_sqrt,b_sqrt,n)-right_fac;
    }
    result = Q1_func_val;
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
float analytic_BEP_MFSK2_nonorthog_8_69(const float Ac,const float Ts,
                                         const float M, const float N0,
                                         const float rho,const float mu,
                                         const std::int32_t n)
{
    if(__builtin_expect(mu!=1.0,0)) { return (std::numeric_limits<float>::quiet_NaN());}
    const float Eb = avg_bit_E_to_carrier_A<float>(Ac,Ts,M);
    const float Eb_div_2N0 = Eb/(N0+N0);
    if(__builtin_expect(rho==0.0f,0))
    {
       return (0.5f*gms::math::expapprox(-Eb_div_2N0));
    }
    float result;
    float Q1_func_val;
    float a;
    float b;
    const float rho_term = std::sqrt(1.0-rho*rho);
    a = Eb_div_2N0*(1.0f-rho_term);
    b = Eb_div_2N0*(1.0f+rho_term);
    const float bessI0_arg = std::sqrt(a*b);
    const float exp_val    = 0.5*gms::math::expapprox(0.5f*(a+b));
    const float bessI0_val = std::cyl_bessel_i<float>(0,bessI0_arg);
    const float right_fac  = exp_val*bessI0_val;
    const float a_sqrt      = std::sqrt(a);
    const float b_sqrt      = std::sqrt(b);
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
        Q1_func_val = gms::math::marcum_Q_approx_chiani(mu,a_sqrt,b_sqrt)-right_fac;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        Q1_func_val = gms::math::marcum_Q_approx_loskot_2T(mu,a_sqrt,b_sqrt)-right_fac;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        Q1_func_val = gms::math::marcum_Q_approx_loskot_3T(mu,a_sqrt,b_sqrt)-right_fac;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q1_func_val = gms::math::marcum_Q_approx_sadhwani_1T(mu,a_sqrt,b_sqrt)-right_fac;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q1_func_val = gms::math::marcum_Q_approx_sadhwani_2T(mu,a_sqrt,b_sqrt)-right_fac;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q1_func_val = gms::math::marcum_Q_approx_sadhwani_4T(mu,a_sqrt,b_sqrt)-right_fac;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_cooper)
    {
        Q1_func_val = gms::math::marcum_Q_approx_cooper(mu,a_sqrt,b_sqrt)-right_fac;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson)
    {
        Q1_func_val = gms::math::marcum_Q_approx_borjesson(mu,a_sqrt,b_sqrt)-right_fac;
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q1_func_val = gms::math::marcum_Q_approx_sadhwani_summed(mu,a_sqrt,b_sqrt,n)-right_fac;
    }
    result = Q1_func_val;
    return (result);
}


#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
double analytic_BEP_MPSK_pair_err_prob_8_80(const double * __restrict__ betak,
                                            const double * __restrict__ betak_hat,
                                            const double Ac,const double Ts,
                                            const double M, const double N0,
                                            const double rhoc,const double phic,
                                            const double mu,const std::int32_t Ns,
                                            const std::int32_t n)                                
{
    const double * __restrict__ p_betak     = betak;
    const double * __restrict__ p_betak_hat = betak_hat;
    double Q1_func_val;
    double result;
    double nu;
    double abs_delta;
    double a;
    double b;
    constexpr const std::complex<double> j(0.0,1.0);
    const double d_Ns = static_cast<double>(Ns);
    const double Es   = avg_symbol_E_to_carrier_A<double>(Ac,Ts,M);
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double cosphic = gms::math::cephes_d::cos(phic);
#else 
    const double cosphic = std::cos(phic);
#endif 
    const double lead_fac = Es/(N0+N0);
    const double snr      = Es/N0;
    const double rho_snr  = rhoc/snr;
    const double invNs    = 1.0/d_Ns;
    const double term1    = = __builtin_fma(cosphic*rho_snr,invNs,1.0);
    std::complex<double> delta_sum(0.0,0.0);
    for(std::int32_t i = 0; i<Ns; ++i) 
    {
        const double bk     = p_betak[i];
        const double bk_hat = p_betak_hat[i];
        const double bkdiff = bk-bk_hat;
        const std::complex<double> cexp_arg = j*bkdiff; 
        delta_sum += std::exp(cexp_arg);
    }
    abs_delta = std::abs(delta_sum);
    nu        = std::arg(abs_delta);
    const double rho_snrp2 = rho_snr*rho_snr;
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double cosnu   = gms::math::cephes_d::cos(nu);
    const double cosnuph = gms::math::cephes_d::cos(nu+phic); 
#else 
    const double cosnu = std::cos(nu);
    const double cosnuph = std::cos(nu+phic);
#endif 
    const double diff_sqr   = (d_Ns*d_Ns)-(abs_delta*abs_delta);
    const double ratio1     = d_Ns-abs_delta*cosnu/(d_Ns+d_Ns*diff_sqr);
    const double term2      = __builtin_fma(rho_snrp2,ratio1,term1);
    const double tmp1       = lead_fac*term2;
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double s_diff_sqr = gms::math::cephes_d::sqrt(diff_sqr);
#else 
    const double s_diff_sqr = std::sqrt(diff_sqr);
#endif 
    const double num        = d_Ns*cosphic-(abs_delta*cosnuph);
    const double ratio2     = num/s_diff_sqr;
    const double term3      = __builtin_fma(s_diff_sqr,ratio2,rho_snr);
    const double tmp2       = snr*term3;
    b                       = tmp1+tmp2;
    a                       = tmp1-tmp2;
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double a_sqrt     = gms::math::cephes_d::sqrt(a);
    const double b_sqrt     = gms::math::cephes_d::sqrt(b);
#else 
    const double a_sqrt     = std::sqrt(a);
    const double b_sqrt     = std::sqrt(b);
#endif 
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
        Q1_func_val = 0.5*(1.0-gms::math::marcum_Q_approx_chiani(mu,b_sqrt,a_sqrt)+
                               gms::math::marcum_Q_approx_chiani(mu,a_sqrt,b_sqrt));
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        Q1_func_val = 0.5*(1.0-gms::math::marcum_Q_approx_loskot_2T(mu,b_sqrt,a_sqrt)+
                               gms::math::marcum_Q_approx_loskot_2T(mu,a_sqrt,b_sqrt));
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        Q1_func_val = 0.5*(1.0-gms::math::marcum_Q_approx_loskot_3T(mu,b_sqrt,a_sqrt)+
                               gms::math::marcum_Q_approx_loskot_3T(mu,a_sqrt,b_sqrt));
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q1_func_val = 0.5*(1.0-gms::math::marcum_Q_approx_sadhwani_1T(mu,b_sqrt,a_sqrt)+
                               gms::math::marcum_Q_approx_sadhwani_1T(mu,a_sqrt,b_sqrt));
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q1_func_val = 0.5*(1.0-gms::math::marcum_Q_approx_sadhwani_2T(mu,b_sqrt,a_sqrt)+
                               gms::math::marcum_Q_approx_sadhwani_2T(mu,a_sqrt,b_sqrt));
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q1_func_val = 0.5*(1.0-gms::math::marcum_Q_approx_sadhwani_4T(mu,b_sqrt,a_sqrt)+
                               gms::math::marcum_Q_approx_sadhwani_4T(mu,a_sqrt,b_sqrt));
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_cooper)
    {
        Q1_func_val = 0.5*(1.0-gms::math::marcum_Q_approx_cooper(mu,b_sqrt,a_sqrt)+
                               gms::math::marcum_Q_approx_cooper(mu,a_sqrt,b_sqrt));
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson)
    {
        Q1_func_val = 0.5*(1.0-gms::math::marcum_Q_approx_borjesson(mu,b_sqrt,a_sqrt)+
                               gms::math::marcum_Q_approx_borjesson(mu,a_sqrt,b_sqrt));
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q1_func_val = 0.5*(1.0-gms::math::marcum_Q_approx_sadhwani_summed(mu,b_sqrt,a_sqrt,n)+
                               gms::math::marcum_Q_approx_sadhwani_summed(mu,a_sqrt,b_sqrt,n));
    }
    result = Q1_func_val;
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
float analytic_BEP_MPSK_pair_err_prob_8_80( const float * __restrict__ betak,
                                            const float * __restrict__ betak_hat,
                                            const float Ac,const float Ts,
                                            const float M, const float N0,
                                            const float rhoc,const float phic,
                                            const float mu,const std::int32_t Ns,
                                            const std::int32_t n)                                
{
    const float * __restrict__ p_betak     = betak;
    const float * __restrict__ p_betak_hat = betak_hat;
    float Q1_func_val;
    float result;
    float nu;
    float abs_delta;
    float a;
    float b;
    constexpr const std::complex<float> j(0.0,1.0);
    const float d_Ns = static_cast<double>(Ns);
    const float Es   = avg_symbol_E_to_carrier_A<float>(Ac,Ts,M);
    const float cosphic = std::cos(phic);
    const float lead_fac = Es/(N0+N0);
    const float snr      = Es/N0;
    const float rho_snr  = rhoc/snr;
    const float invNs    = 1.0f/d_Ns;
    const float term1    = = __builtin_fmaf(cosphic*rho_snr,invNs,1.0);
    std::complex<float> delta_sum(0.0,0.0);
    for(std::int32_t i = 0; i<Ns; ++i) 
    {
        const float bk     = p_betak[i];
        const float bk_hat = p_betak_hat[i];
        const float bkdiff = bk-bk_hat;
        const std::complex<float> cexp_arg = j*bkdiff; 
        delta_sum += std::exp(cexp_arg);
    }
    abs_delta = std::abs(delta_sum);
    nu        = std::arg(abs_delta);
    const float rho_snrp2 = rho_snr*rho_snr;
    const float cosnu = std::cos(nu);
    const float cosnuph = std::cos(nu+phic);
    const float diff_sqr   = (d_Ns*d_Ns)-(abs_delta*abs_delta);
    const float ratio1     = d_Ns-abs_delta*cosnu/(d_Ns+d_Ns*diff_sqr);
    const float term2      = __builtin_fmaf(rho_snrp2,ratio1,term1);
    const float tmp1       = lead_fac*term2;
    const float s_diff_sqr = std::sqrt(diff_sqr); 
    const float num        = d_Ns*cosphic-(abs_delta*cosnuph);
    const float ratio2     = num/s_diff_sqr;
    const float term3      = __builtin_fmaf(s_diff_sqr,ratio2,rho_snr);
    const float tmp2       = snr*term3;
    b                       = tmp1+tmp2;
    a                       = tmp1-tmp2;
    const float a_sqrt     = std::sqrt(a);
    const float b_sqrt     = std::sqrt(b);
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
        Q1_func_val = 0.5*(1.0-gms::math::marcum_Q_approx_chiani(mu,b_sqrt,a_sqrt)+
                               gms::math::marcum_Q_approx_chiani(mu,a_sqrt,b_sqrt));
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        Q1_func_val = 0.5*(1.0-gms::math::marcum_Q_approx_loskot_2T(mu,b_sqrt,a_sqrt)+
                               gms::math::marcum_Q_approx_loskot_2T(mu,a_sqrt,b_sqrt));
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        Q1_func_val = 0.5*(1.0-gms::math::marcum_Q_approx_loskot_3T(mu,b_sqrt,a_sqrt)+
                               gms::math::marcum_Q_approx_loskot_3T(mu,a_sqrt,b_sqrt));
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q1_func_val = 0.5*(1.0-gms::math::marcum_Q_approx_sadhwani_1T(mu,b_sqrt,a_sqrt)+
                               gms::math::marcum_Q_approx_sadhwani_1T(mu,a_sqrt,b_sqrt));
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q1_func_val = 0.5*(1.0-gms::math::marcum_Q_approx_sadhwani_2T(mu,b_sqrt,a_sqrt)+
                               gms::math::marcum_Q_approx_sadhwani_2T(mu,a_sqrt,b_sqrt));
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q1_func_val = 0.5*(1.0-gms::math::marcum_Q_approx_sadhwani_4T(mu,b_sqrt,a_sqrt)+
                               gms::math::marcum_Q_approx_sadhwani_4T(mu,a_sqrt,b_sqrt));
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_cooper)
    {
        Q1_func_val = 0.5*(1.0-gms::math::marcum_Q_approx_cooper(mu,b_sqrt,a_sqrt)+
                               gms::math::marcum_Q_approx_cooper(mu,a_sqrt,b_sqrt));
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson)
    {
        Q1_func_val = 0.5*(1.0-gms::math::marcum_Q_approx_borjesson(mu,b_sqrt,a_sqrt)+
                               gms::math::marcum_Q_approx_borjesson(mu,a_sqrt,b_sqrt));
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q1_func_val = 0.5*(1.0-gms::math::marcum_Q_approx_sadhwani_summed(mu,b_sqrt,a_sqrt,n)+
                               gms::math::marcum_Q_approx_sadhwani_summed(mu,a_sqrt,b_sqrt,n));
    }
    result = Q1_func_val;
    return (result);
}    

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
double analytic_BEP_DPSK_up_bound_8_95(const double Ac,const double Ts,
                                       const double M, const double N0,
                                       const std::int32_t n)
{
    double result;
    double Q_func_val;
    const double PI_div_M = 3.1415926535897932384626433832795/M;
    const double Es       = avg_symbol_E_to_carrier_A<double>(Ac,Ts,M);
    const double snr      = (Es+Es)/N0;
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double cosPIM   = gms::math::cephes_d::cos(PI_div_M);
    const double Q_func_arg = gms::math::cephes_d::sqrt(snr*(1-cosPIM));
#else 
    const double cosPIM   = std::cos(PI_div_M);
    const double Q_func_arg = std::sqrt(snr*(1-cosPIM));
#endif
    const double cos_ratio = (1.0+cosPIM)/(cosPIM+cosPIM);
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double left_term = 2.06*gms::math::cephes_d::sqrt(cos_ratio);
#else 
    const double left_term = 2.06*std::sqrt(cos_ratio);
#endif 
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
       Q_func_val  = left_term*gms::math::gaussian_Q_approx_chiani(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = left_term*gms::math::gaussian_Q_approx_loskot_2T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = left_term*gms::math::gaussian_Q_approx_loskot_3T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = left_term*gms::math::gaussian_Q_approx_sadhwani_1T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = left_term*gms::math::gaussian_Q_approx_sadhwani_2T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = left_term*gms::math::gaussian_Q_approx_sadhwani_4T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_cooper)
    {
        Q_func_val = left_term*gms::math::gaussian_Q_approx_cooper(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson)
    {
        Q_func_val = left_term*gms::math::gaussian_Q_approx_borjesson(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = left_term*gms::math::gaussian_Q_approx_sadhwani_summed(Q_func_arg,n);
    }
    result = Q_func_val;
    return (result);
}
                                         
#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<Gaussian_Q_approxmations_t Q_func_approx>
__ATTR_ALWAYS_INLINE__
static inline
float analytic_BEP_DPSK_up_bound_8_95(const float Ac,const float Ts,
                                       const float M, const float N0,
                                       const std::int32_t n)
{
    float result;
    float Q_func_val;
    const float PI_div_M = 3.1415926535897932384626433832795f/M;
    const float Es       = avg_symbol_E_to_carrier_A<float>(Ac,Ts,M);
    const float snr      = (Es+Es)/N0;
    const float cosPIM   = std::cos(PI_div_M);
    const float Q_func_arg = std::sqrt(snr*(1-cosPIM));
    const float cos_ratio = (1.0+cosPIM)/(cosPIM+cosPIM);
    const float left_term = 2.06f*std::sqrt(cos_ratio);
    if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani)
    {
       Q_func_val  = left_term*gms::math::gaussian_Q_approx_chiani(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T)
    {
        Q_func_val = left_term*gms::math::gaussian_Q_approx_loskot_2T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T)
    {
        Q_func_val = left_term*gms::math::gaussian_Q_approx_loskot_3T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T)
    {
        Q_func_val = left_term*gms::math::gaussian_Q_approx_sadhwani_1T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T)
    {
        Q_func_val = left_term*gms::math::gaussian_Q_approx_sadhwani_2T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T)
    {
        Q_func_val = left_term*gms::math::gaussian_Q_approx_sadhwani_4T(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_cooper)
    {
        Q_func_val = left_term*gms::math::gaussian_Q_approx_cooper(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson)
    {
        Q_func_val = left_term*gms::math::gaussian_Q_approx_borjesson(Q_func_arg);
    }
    else if constexpr(Q_func_approx==Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed)
    {
        Q_func_val = left_term*gms::math::gaussian_Q_approx_sadhwani_summed(Q_func_arg,n);
    }
    result = Q_func_val;
    return (result);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
double analytic_BEP_DPSK_chernoff_bound_8_96(const double Ac,const double Ts,
                                             const double M, const double N0)
{
    const double PI_div_M = 3.1415926535897932384626433832795/M;
    const double Es       = avg_symbol_E_to_carrier_A<double>(Ac,Ts,M);
    const double snr      = (Es+Es)/N0;
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double cosPIM   = gms::math::cephes_d::cos(PI_div_M);
    const double exp_func_arg = gms::math::cephes_d::sqrt(snr*(1-cosPIM));
#else 
    const double cosPIM   = std::cos(PI_div_M);
    const double exp_func_arg = std::sqrt(snr*(1-cosPIM));
#endif
    const double cos_ratio = (1.0+cosPIM)/(cosPIM+cosPIM);
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double left_term = 1.03*gms::math::cephes_d::sqrt(cos_ratio);
#else 
    const double left_term = 1.03*std::sqrt(cos_ratio);
#endif 
    return (left_term*gms::math::expapprox_d(-exp_func_arg));
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
float analytic_BEP_DPSK_chernoff_bound_8_96(const float Ac,const float Ts,
                                             const float M, const float N0)
{
    const float PI_div_M = 3.1415926535897932384626433832795f/M;
    const float Es       = avg_symbol_E_to_carrier_A<float>(Ac,Ts,M);
    const float snr      = (Es+Es)/N0;
    const float cosPIM   = std::cos(PI_div_M);
    const float exp_func_arg = std::sqrt(snr*(1-cosPIM));
    const float cos_ratio = (1.0f+cosPIM)/(cosPIM+cosPIM); 
    const float left_term = 1.03f*std::sqrt(cos_ratio);
    return (left_term*gms::math::expapprox(-exp_func_arg));
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
double analytic_SEP_QAM4_8_107(const double gamma) //snr average
{
    constexpr const double four_div_pi = 1.2732395447351626861510701069801;
    const double sqrt_arg1              = gamma/(1.0+gamma);
    const double sqrt_arg2              = (1.0+gamma)/gamma;
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double sqrt_val1              = gms::math::cephes_d::sqrt(sqrt_arg1);
    const double sqrt_val2              = gms::math::cephes_d::sqrt(sqrt_arg2);
#else 
    const double sqrt_val1              = std::sqrt(sqrt_arg1);
    const double sqrt_val2              = std::sqrt(sqrt_arg2);
#endif 
    const double term1                  = 1.0-sqrt_val1;
    const double atan_val               = four_div_pi*std::atan(sqrt_val2);
    const double term2                  = 0.25*term1*atan_val;
    return (term1-term2);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
float analytic_SEP_QAM4_8_107(const float gamma) //snr average
{
    constexpr const float four_div_pi = 1.2732395447351626861510701069801f;
    const float sqrt_arg1              = gamma/(1.0f+gamma);
    const float sqrt_arg2              = (1.0f+gamma)/gamma;
    const float sqrt_val1              = std::sqrt(sqrt_arg1);
    const float sqrt_val2              = std::sqrt(sqrt_arg2);
    const float term1                  = 1.0f-sqrt_val1;
    const float atan_val               = four_div_pi*std::atan(sqrt_val2);
    const float term2                  = 0.25f*term1*atan_val;
    return (term1-term2);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
double analytic_SEP_MPSK_rayleigh_chan_8_112(const double gamma,
                                             const double M) 
{
    const double pi_div_M = 3.1415926535897932384626433832795/M;
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double sin_tmp  = gms::math::cephes_d::sin(pi_div_M);
#else 
    const double sin_tmp  = std::sin(pi_div_M);
#endif
    const double term1    = (M-1.0)/M;
    const double g_psk    = sin_tmp*sin_tmp;
    const double g_pskgamm= g_psk*gamma;
    const double term2    = g_pskgamm/(1.0+g_pskgamm);
#if (ANALYTIC_BEP_SEP_CH8_CEPHES_DOUBLE) == 1
    const double sqrt_val = gms::math::cephes_d::sqrt(term2);
#else 
    const double sqrt_val = std::sqrt(term2);
#endif 
    const double term3    = M/((M-1.0)*3.1415926535897932384626433832795);
    const double cot_val  = std::tan(1.5707963267948966192313216916398-pi_div_M);
    const double term4    = 1.0-sqrt_val*term3;
    const double atan_val = 1.5707963267948966192313216916398+std::atan(sqrt_val*cot_val);
    return (term1*term4*atan_val);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
float analytic_SEP_MPSK_rayleigh_chan_8_112(const float gamma,
                                            const float M) 
{
    const float pi_div_M = 3.1415926535897932384626433832795f/M;
    const float sin_tmp  = std::sin(pi_div_M);
    const float term1    = (M-1.0f)/M;
    const float g_psk    = sin_tmp*sin_tmp;
    const float g_pskgamm= g_psk*gamma;
    const float term2    = g_pskgamm/(1.0f+g_pskgamm);
    const float sqrt_val = std::sqrt(term2);
    const float term3    = M/((M-1.0f)*3.1415926535897932384626433832795f);
    const float cot_val  = std::tan(1.5707963267948966192313216916398f-pi_div_M);
    const float term4    = 1.0f-sqrt_val*term3;
    const float atan_val = 1.5707963267948966192313216916398f+std::atan(sqrt_val*cot_val);
    return (term1*term4*atan_val);
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
double riccian_factor_8_144(const double m_I,const double m_Q,
                            const double sigma)
{
    const double m_I_p2   = m_I*m_I;
    const double m_Q_p2   = m_Q*m_Q;
    const double sigma_p2 = sigma*sigma;
    return (m_I_p2+m_Q_p2/(sigma_p2+sigma_p2)); 
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
float riccian_factor_8_144(const float m_I,const float m_Q,
                            const float sigma)
{
    const float m_I_p2   = m_I*m_I;
    const float m_Q_p2   = m_Q*m_Q;
    const float sigma_p2 = sigma*sigma;
    return (m_I_p2+m_Q_p2/(sigma_p2+sigma_p2)); 
}

#if (ANALYTIC_BEP_SEP_CH8_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
__ATTR_ALWAYS_INLINE__
static inline
double spec_pow_to_rand_pow_8_145(const float m_I,const float m_Q,
                                  const float sigma)
{
    const double sigma_p2 = sigma*sigma;
    return ((sigma_p2+sigma_p2)*(1.0+riccian_factor_8_144(m_I,m_Q,sigma)));
}

}

}
#endif /*__GMS_ANALYTIC_BEP_SEP_CH8_H__*/