
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
#include "GMS_config.h"
#include "GMS_gaussian_Q_approx.h"
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
    double Q_func_val     = 0.0;
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

}

}
#endif /*__GMS_INTEGRANDS_FUNC_CH8_V2_H__*/