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

#ifndef __GMS_ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_H__
#define __GMS_ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_H__ 061020260556

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_analytic_bep_sep_ch8.h"
#include "GMS_analytic_bep_sep_ch8_sse.h"



namespace file_info 
{

     static const unsigned int GMS_ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F32_MAJOR = 1;
     static const unsigned int GMS_ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F32_MINOR = 1;
     static const unsigned int GMS_ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F32_MICRO = 0;
     static const unsigned int GMS_ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F32_FULLVER =
       1000U*GMS_ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F32_MAJOR+100U*GMS_ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F32_MINOR+
       10U*GMS_ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F32_MICRO;
     static const char GMS_ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F32_CREATION_DATE[] = "06-10-2026 05:55PM +00200 (TUE 06 OCT 2026 GMT+2)";
     static const char GMS_ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F32_BUILD_DATE[]    = __DATE__; 
     static const char GMS_ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F32_BUILD_TIME[]    = __TIME__;
     static const char GMS_ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F32_AUTHOR[]        = "Programmer: Bernard Gingold, beniekg@gmail.com";
     static const char GMS_ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F32_SYNOPSIS[]      = "Analytic Formulae of the Bit and Symbol Error Probabilities based on the approximated Gaussian-Q function.\
	                                                          Based on the chapter 8 of M.K Simon, M.S. Alouini: Digital Communication over Fading Channels 1st ed\
															  ISBN-13 978-0471317791";

}

#if !defined(ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F32_OVERRIDE_COMPILER_CMD_LINE)
#define ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F32_OVERRIDE_COMPILER_CMD_LINE 0
#endif

#if !defined(ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_USE_PEEL_LOOP)
#define ANALYTIC_BEP_SEP_CH8_ARRAY1D_USE_F32_PEEL_LOOP 1
#endif 

#if !defined(ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH)
#define ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH 0
#endif 
 
// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F32_USE_PMC_INSTRUMENTATION)
#define ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F32_USE_PMC_INSTRUMENTATION 0
#endif 

#if (ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F32_USE_PMC_INSTRUMENTATION) == 1
#include "GMS_hw_perf_macros.h"

#define PMC_VARS\
    uint64_t prog_counters_start[4] = {};\
    uint64_t prog_counters_stop[4]  = {};\
    uint64_t tsc_start,tsc_stop;\
    uint64_t act_cyc_start,act_cyc_stop;\
    uint64_t ref_cyc_start,ref_cyc_stop;\
    [[maybe_unused]] uint64_t dummy1;\
    [[maybe_unused]] uint64_t dummy2;\
    [[maybe_unused]] uint64_t dummy3;\
    int32_t core_counter_width;\
    double utilization,nom_ghz,avg_ghz;

#endif


namespace gms 
{

namespace fading_channel
{

#if (ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F32_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=sse
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<typename callable_sse_kernel,
         typename callable_scalar_kernel>
__ATTR_ALWAYS_INLINE__ 
static inline 
std::int32_t 
analytic_SEP_BEP_generic_4xf32_array1d_u16x(const float * __restrict__ Ac, 
                                        const float * __restrict__ Ts,
                                        const float * __restrict__ M,
                                        const float * __restrict__ N0,
                                        float       * __restrict__ sep_mam,
                                        const std::int32_t n,
                                        std::size_t sz,
                                        const callable_sse_kernel     & sse_kernel,
                                        const callable_scalar_kernel  & scalar_kernel) 
{
    if(__builtin_expect(sz<1ull,0)) {return (-1);}
    __m128 vAc;
    __m128 vTs;
    __m128 vM;
    __m128 vN0;
    const float * __restrict__ p_Ac      = Ac;
    const float * __restrict__ p_Ts      = Ts;
    const float * __restrict__ p_M       = M;
    const float * __restrict__ p_N0      = N0;
    float       * __restrict__ p_sep_mam = sep_mam;
    std::size_t             a_sz         = sz;
    if(a_sz>=1ull && a_sz<4ull)
    {
       for(std::uint32_t i = 0; i<4; ++i) 
       {
           const float valAc = p_Ac[i];
           const float valTs = p_Ts[i];
           const float valM  = p_M[i];
           const float valN0 = p_N0[i];
           p_sep_mam[i]      = scalar_kernel(valAc,valTs,valM,valN0,n);
       }
       return (0);
    }
    else if(a_sz==4ull)
    {
        vAc = _mm_loadu_ps(&p_Ac[0]);
        vTs = _mm_loadu_ps(&p_Ts[0]);
        vM  = _mm_loadu_ps(&p_M[0]);
        vN0 = _mm_loadu_ps(&p_N0[0]);
        _mm_storeu_ps(&p_sep_mam[0], sse_kernel(vAc,vTs,vM,vN0,n));
        return (0);
    }
    else if(a_sz>4ull && a_sz<=64ull)
    {
        std::size_t i,j;
        for(i = 0ull; i!=ROUND_DOWN(a_sz,3ull); i += 4ull)
        {
            vAc = _mm_loadu_ps(&p_Ac[i]);
            vTs = _mm_loadu_ps(&p_Ts[i]);
            vM  = _mm_loadu_ps(&p_M[i]);
            vN0 = _mm_loadu_ps(&p_N0[i]);
            _mm_storeu_ps(&p_sep_mam[i], sse_kernel(vAc,vTs,vM,vN0,n));
        }
        for(j = i; j<a_sz; ++j)  
        {
            const float valAc = p_Ac[j];
            const float valTs = p_Ts[j];
            const float valM  = p_M[j];
            const float valN0 = p_N0[j];
            p_sep_mam[j]      = scalar_kernel(valAc,valTs,valM,valN0,n);
        }
        return (0);
    }
    else if(a_sz>64ull)
    {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_USE_PEEL_LOOP) == 1
       while(((std::uintptr_t)&p_sep_mam & 15ull) && a_sz)
       {
           const float valAc = *p_Ac;
           const float valTs = *p_Ts;
           const float valM  = *p_M;
           const float valN0 = *p_N0;
           *p_sep_mam        =  scalar_kernel(valAc,valTs,valM,valN0,n);
           p_Ac++;
           p_Ts++;
           p_M++;
           p_N0++;
           p_sep_mam++;
           a_sz--;
       }
#else 
#error "***[FATAL]*** -- The loop peeling procedure must be always equal to 1!!"
#endif
       __m128 xmm0;
       __m128 xmm1;
       __m128 xmm2;
       __m128 xmm3;
       __m128 xmm4;
       __m128 xmm5;
       __m128 xmm6;
       __m128 xmm7;
       __m128 xmm8;
       __m128 xmm9;
       __m128 xmm10;
       __m128 xmm11;
       __m128 xmm12;
       __m128 xmm13;
       __m128 xmm14;
       std::size_t i,j;
       constexpr std::size_t L2_offset = 8000ull;
       constexpr std::size_t L1D_boundary = L2_offset;
    if(a_sz <= L1D_boundary)
    {
       for(i = 0ull; (i+63ull) < a_sz; i += 64ull) 
       {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
#endif 
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+8ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+8ull]);
            xmm12 = _mm_load_ps(&p_M[i+8ull]);
            xmm13 = _mm_load_ps(&p_N0[i+8ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+8ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+12ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+12ull]);
            xmm2 = _mm_load_ps(&p_M[i+12ull]);
            xmm3 = _mm_load_ps(&p_N0[i+12ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+12ull],xmm4);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+16ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+16ull],_MM_HINT_T0);
#endif
            xmm5 = _mm_load_ps(&p_Ac[i+16ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+16ull]);
            xmm7 = _mm_load_ps(&p_M[i+16ull]);
            xmm8 = _mm_load_ps(&p_N0[i+16ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+16ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+20ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+20ull]);
            xmm12 = _mm_load_ps(&p_M[i+20ull]);
            xmm13 = _mm_load_ps(&p_N0[i+20ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+20ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+24ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+24ull]);
            xmm2 = _mm_load_ps(&p_M[i+24ull]);
            xmm3 = _mm_load_ps(&p_N0[i+24ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+24ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+28ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+28ull]);
            xmm7 = _mm_load_ps(&p_M[i+28ull]);
            xmm8 = _mm_load_ps(&p_N0[i+28ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+28ull],xmm9);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+32ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+32ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+32ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+32ull],_MM_HINT_T0);
#endif

            xmm10 = _mm_load_ps(&p_Ac[i+32ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+32ull]);
            xmm12 = _mm_load_ps(&p_M[i+32ull]);
            xmm13 = _mm_load_ps(&p_N0[i+32ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+32ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+36ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+36ull]);
            xmm2 = _mm_load_ps(&p_M[i+36ull]);
            xmm3 = _mm_load_ps(&p_N0[i+36ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+36ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+40ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+40ull]);
            xmm7 = _mm_load_ps(&p_M[i+40ull]);
            xmm8 = _mm_load_ps(&p_N0[i+40ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+40ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+44ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+44ull]);
            xmm12 = _mm_load_ps(&p_M[i+44ull]);
            xmm13 = _mm_load_ps(&p_N0[i+44ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+44ull],xmm14);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+48ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+48ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+48ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+48ull],_MM_HINT_T0);
#endif
            
            xmm0 = _mm_load_ps(&p_Ac[i+48ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+48ull]);
            xmm2 = _mm_load_ps(&p_M[i+48ull]);
            xmm3 = _mm_load_ps(&p_N0[i+48ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+48ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+52ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+52ull]);
            xmm7 = _mm_load_ps(&p_M[i+52ull]);
            xmm8 = _mm_load_ps(&p_N0[i+52ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+52ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+56ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+56ull]);
            xmm12 = _mm_load_ps(&p_M[i+56ull]);
            xmm13 = _mm_load_ps(&p_N0[i+56ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+56ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+60ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+60ull]);
            xmm2 = _mm_load_ps(&p_M[i+60ull]);
            xmm3 = _mm_load_ps(&p_N0[i+60ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+60ull],xmm4);
       }
 
       for(; (i+47ull) < a_sz; i += 48ull)
       {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
#endif 
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+8ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+8ull]);
            xmm12 = _mm_load_ps(&p_M[i+8ull]);
            xmm13 = _mm_load_ps(&p_N0[i+8ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+8ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+12ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+12ull]);
            xmm2 = _mm_load_ps(&p_M[i+12ull]);
            xmm3 = _mm_load_ps(&p_N0[i+12ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+12ull],xmm4);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+16ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+16ull],_MM_HINT_T0);
#endif
            xmm5 = _mm_load_ps(&p_Ac[i+16ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+16ull]);
            xmm7 = _mm_load_ps(&p_M[i+16ull]);
            xmm8 = _mm_load_ps(&p_N0[i+16ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+16ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+20ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+20ull]);
            xmm12 = _mm_load_ps(&p_M[i+20ull]);
            xmm13 = _mm_load_ps(&p_N0[i+20ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+20ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+24ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+24ull]);
            xmm2 = _mm_load_ps(&p_M[i+24ull]);
            xmm3 = _mm_load_ps(&p_N0[i+24ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+24ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+28ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+28ull]);
            xmm7 = _mm_load_ps(&p_M[i+28ull]);
            xmm8 = _mm_load_ps(&p_N0[i+28ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+28ull],xmm9);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+32ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+32ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+32ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+32ull],_MM_HINT_T0);
#endif

            xmm10 = _mm_load_ps(&p_Ac[i+32ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+32ull]);
            xmm12 = _mm_load_ps(&p_M[i+32ull]);
            xmm13 = _mm_load_ps(&p_N0[i+32ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+32ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+36ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+36ull]);
            xmm2 = _mm_load_ps(&p_M[i+36ull]);
            xmm3 = _mm_load_ps(&p_N0[i+36ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+36ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+40ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+40ull]);
            xmm7 = _mm_load_ps(&p_M[i+40ull]);
            xmm8 = _mm_load_ps(&p_N0[i+40ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+40ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+44ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+44ull]);
            xmm12 = _mm_load_ps(&p_M[i+44ull]);
            xmm13 = _mm_load_ps(&p_N0[i+44ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+44ull],xmm14);
       }

       for(; (i+31ull) < a_sz; i += 32ull) 
       {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
#endif 
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+8ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+8ull]);
            xmm12 = _mm_load_ps(&p_M[i+8ull]);
            xmm13 = _mm_load_ps(&p_N0[i+8ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+8ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+12ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+12ull]);
            xmm2 = _mm_load_ps(&p_M[i+12ull]);
            xmm3 = _mm_load_ps(&p_N0[i+12ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+12ull],xmm4);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+16ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+16ull],_MM_HINT_T0);
#endif
            xmm5 = _mm_load_ps(&p_Ac[i+16ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+16ull]);
            xmm7 = _mm_load_ps(&p_M[i+16ull]);
            xmm8 = _mm_load_ps(&p_N0[i+16ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+16ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+20ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+20ull]);
            xmm12 = _mm_load_ps(&p_M[i+20ull]);
            xmm13 = _mm_load_ps(&p_N0[i+20ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+20ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+24ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+24ull]);
            xmm2 = _mm_load_ps(&p_M[i+24ull]);
            xmm3 = _mm_load_ps(&p_N0[i+24ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+24ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+28ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+28ull]);
            xmm7 = _mm_load_ps(&p_M[i+28ull]);
            xmm8 = _mm_load_ps(&p_N0[i+28ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+28ull],xmm9);
       }

       for(; (i+19ull) < a_sz; i += 20ull)
       {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
#endif 
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+8ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+8ull]);
            xmm12 = _mm_load_ps(&p_M[i+8ull]);
            xmm13 = _mm_load_ps(&p_N0[i+8ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+8ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+12ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+12ull]);
            xmm2 = _mm_load_ps(&p_M[i+12ull]);
            xmm3 = _mm_load_ps(&p_N0[i+12ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+12ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+16ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+16ull]);
            xmm7 = _mm_load_ps(&p_M[i+16ull]);
            xmm8 = _mm_load_ps(&p_N0[i+16ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+16ull],xmm9);
       }

       for(; (i+7ull) < a_sz; i += 8ull)
       {
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);
       }

       for(; (i+3ull) < a_sz; i += 4ull)
       {
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);
       }

       for(; (i+0ull) < a_sz; i += 1ull)
       {
            const float valAc = p_Ac[i];
            const float valTs = p_Ts[i];
            const float valM  = p_M[i];
            const float valN0 = p_N0[i];
            p_sep_mam[i]      = scalar_kernel(valAc,valTs,valM,valN0,n);
       }
       return (0);
    }
    else 
    {
       for(i = 0ull; (i+63ull) < a_sz; i += 64ull) 
       {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+0ull], _MM_HIN1_T0);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+0ull],_MM_HINT_T1);
#endif 
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+8ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+8ull]);
            xmm12 = _mm_load_ps(&p_M[i+8ull]);
            xmm13 = _mm_load_ps(&p_N0[i+8ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+8ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+12ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+12ull]);
            xmm2 = _mm_load_ps(&p_M[i+12ull]);
            xmm3 = _mm_load_ps(&p_N0[i+12ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+12ull],xmm4);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+16ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+16ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+16ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+16ull], _MM_HINT_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+16ull],_MM_HINT_T1);
#endif
            xmm5 = _mm_load_ps(&p_Ac[i+16ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+16ull]);
            xmm7 = _mm_load_ps(&p_M[i+16ull]);
            xmm8 = _mm_load_ps(&p_N0[i+16ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+16ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+20ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+20ull]);
            xmm12 = _mm_load_ps(&p_M[i+20ull]);
            xmm13 = _mm_load_ps(&p_N0[i+20ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+20ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+24ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+24ull]);
            xmm2 = _mm_load_ps(&p_M[i+24ull]);
            xmm3 = _mm_load_ps(&p_N0[i+24ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+24ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+28ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+28ull]);
            xmm7 = _mm_load_ps(&p_M[i+28ull]);
            xmm8 = _mm_load_ps(&p_N0[i+28ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+28ull],xmm9);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+32ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+32ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+32ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+32ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+32ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+32ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+32ull], _MM_HINT_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+32ull],_MM_HINT_T1);
#endif

            xmm10 = _mm_load_ps(&p_Ac[i+32ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+32ull]);
            xmm12 = _mm_load_ps(&p_M[i+32ull]);
            xmm13 = _mm_load_ps(&p_N0[i+32ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+32ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+36ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+36ull]);
            xmm2 = _mm_load_ps(&p_M[i+36ull]);
            xmm3 = _mm_load_ps(&p_N0[i+36ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+36ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+40ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+40ull]);
            xmm7 = _mm_load_ps(&p_M[i+40ull]);
            xmm8 = _mm_load_ps(&p_N0[i+40ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+40ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+44ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+44ull]);
            xmm12 = _mm_load_ps(&p_M[i+44ull]);
            xmm13 = _mm_load_ps(&p_N0[i+44ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+44ull],xmm14);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+48ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+48ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+48ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+48ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+48ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+48ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+48ull], _MM_HINT_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+48ull],_MM_HINT_T1);
#endif
            
            xmm0 = _mm_load_ps(&p_Ac[i+48ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+48ull]);
            xmm2 = _mm_load_ps(&p_M[i+48ull]);
            xmm3 = _mm_load_ps(&p_N0[i+48ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+48ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+52ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+52ull]);
            xmm7 = _mm_load_ps(&p_M[i+52ull]);
            xmm8 = _mm_load_ps(&p_N0[i+52ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+52ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+56ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+56ull]);
            xmm12 = _mm_load_ps(&p_M[i+56ull]);
            xmm13 = _mm_load_ps(&p_N0[i+56ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+56ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+60ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+60ull]);
            xmm2 = _mm_load_ps(&p_M[i+60ull]);
            xmm3 = _mm_load_ps(&p_N0[i+60ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+60ull],xmm4);
       }
 
       for(; (i+47ull) < a_sz; i += 48ull)
       {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+0ull], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+0ull],_MM_HINT_T1);
#endif 
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+8ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+8ull]);
            xmm12 = _mm_load_ps(&p_M[i+8ull]);
            xmm13 = _mm_load_ps(&p_N0[i+8ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+8ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+12ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+12ull]);
            xmm2 = _mm_load_ps(&p_M[i+12ull]);
            xmm3 = _mm_load_ps(&p_N0[i+12ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+12ull],xmm4);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+16ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+16ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+16ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+16ull], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+16ull],_MM_HINT_T1);
#endif
            xmm5 = _mm_load_ps(&p_Ac[i+16ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+16ull]);
            xmm7 = _mm_load_ps(&p_M[i+16ull]);
            xmm8 = _mm_load_ps(&p_N0[i+16ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+16ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+20ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+20ull]);
            xmm12 = _mm_load_ps(&p_M[i+20ull]);
            xmm13 = _mm_load_ps(&p_N0[i+20ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+20ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+24ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+24ull]);
            xmm2 = _mm_load_ps(&p_M[i+24ull]);
            xmm3 = _mm_load_ps(&p_N0[i+24ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+24ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+28ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+28ull]);
            xmm7 = _mm_load_ps(&p_M[i+28ull]);
            xmm8 = _mm_load_ps(&p_N0[i+28ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+28ull],xmm9);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+32ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+32ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+32ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+32ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+32ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+32ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+32ull], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+32ull],_MM_HINT_T1);
#endif

            xmm10 = _mm_load_ps(&p_Ac[i+32ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+32ull]);
            xmm12 = _mm_load_ps(&p_M[i+32ull]);
            xmm13 = _mm_load_ps(&p_N0[i+32ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+32ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+36ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+36ull]);
            xmm2 = _mm_load_ps(&p_M[i+36ull]);
            xmm3 = _mm_load_ps(&p_N0[i+36ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+36ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+40ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+40ull]);
            xmm7 = _mm_load_ps(&p_M[i+40ull]);
            xmm8 = _mm_load_ps(&p_N0[i+40ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+40ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+44ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+44ull]);
            xmm12 = _mm_load_ps(&p_M[i+44ull]);
            xmm13 = _mm_load_ps(&p_N0[i+44ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+44ull],xmm14);
       }

       for(; (i+31ull) < a_sz; i += 32ull) 
       {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+0ull], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+0ull],_MM_HINT_T1);
#endif 
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+8ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+8ull]);
            xmm12 = _mm_load_ps(&p_M[i+8ull]);
            xmm13 = _mm_load_ps(&p_N0[i+8ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+8ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+12ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+12ull]);
            xmm2 = _mm_load_ps(&p_M[i+12ull]);
            xmm3 = _mm_load_ps(&p_N0[i+12ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+12ull],xmm4);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+16ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+16ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+16ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+16ull], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+16ull],_MM_HINT_T1);
#endif
            xmm5 = _mm_load_ps(&p_Ac[i+16ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+16ull]);
            xmm7 = _mm_load_ps(&p_M[i+16ull]);
            xmm8 = _mm_load_ps(&p_N0[i+16ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+16ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+20ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+20ull]);
            xmm12 = _mm_load_ps(&p_M[i+20ull]);
            xmm13 = _mm_load_ps(&p_N0[i+20ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+20ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+24ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+24ull]);
            xmm2 = _mm_load_ps(&p_M[i+24ull]);
            xmm3 = _mm_load_ps(&p_N0[i+24ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+24ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+28ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+28ull]);
            xmm7 = _mm_load_ps(&p_M[i+28ull]);
            xmm8 = _mm_load_ps(&p_N0[i+28ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+28ull],xmm9);
       }

       for(; (i+19ull) < a_sz; i += 20ull)
       {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+0ull], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+0ull],_MM_HINT_T1);
#endif 
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+8ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+8ull]);
            xmm12 = _mm_load_ps(&p_M[i+8ull]);
            xmm13 = _mm_load_ps(&p_N0[i+8ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+8ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+12ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+12ull]);
            xmm2 = _mm_load_ps(&p_M[i+12ull]);
            xmm3 = _mm_load_ps(&p_N0[i+12ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+12ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+16ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+16ull]);
            xmm7 = _mm_load_ps(&p_M[i+16ull]);
            xmm8 = _mm_load_ps(&p_N0[i+16ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+16ull],xmm9);
       }

       for(; (i+7ull) < a_sz; i += 8ull)
       {
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);
       }

       for(; (i+3ull) < a_sz; i += 4ull)
       {
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);
       }

       for(; (i+0ull) < a_sz; i += 1ull)
       {
            const float valAc = p_Ac[i];
            const float valTs = p_Ts[i];
            const float valM  = p_M[i];
            const float valN0 = p_N0[i];
            p_sep_mam[i]      = scalar_kernel(valAc,valTs,valM,valN0,n);
       }
       return (0);
    }

    }
}

//////////////////////////////////////////////////////////////////////////////////

#if (ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F32_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=sse
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<typename callable_sse_kernel,
         typename callable_scalar_kernel>
__ATTR_ALWAYS_INLINE__ 
static inline 
std::int32_t 
analytic_SEP_BEP_generic_4xf32_array1d_u12x(const float * __restrict__ Ac, 
                                        const float * __restrict__ Ts,
                                        const float * __restrict__ M,
                                        const float * __restrict__ N0,
                                        float       * __restrict__ sep_mam,
                                        const std::int32_t n,
                                        std::size_t sz,
                                        const callable_sse_kernel     & sse_kernel,
                                        const callable_scalar_kernel  & scalar_kernel) 
{
    if(__builtin_expect(sz<1ull,0)) {return (-1);}
    __m128 vAc;
    __m128 vTs;
    __m128 vM;
    __m128 vN0;
    const float * __restrict__ p_Ac      = Ac;
    const float * __restrict__ p_Ts      = Ts;
    const float * __restrict__ p_M       = M;
    const float * __restrict__ p_N0      = N0;
    float       * __restrict__ p_sep_mam = sep_mam;
    std::size_t             a_sz         = sz;
    if(a_sz>=1ull && a_sz<4ull)
    {
       for(std::uint32_t i = 0; i<4; ++i) 
       {
           const float valAc = p_Ac[i];
           const float valTs = p_Ts[i];
           const float valM  = p_M[i];
           const float valN0 = p_N0[i];
           p_sep_mam[i]      = scalar_kernel(valAc,valTs,valM,valN0,n);
       }
       return (0);
    }
    else if(a_sz==4ull)
    {
        vAc = _mm_loadu_ps(&p_Ac[0]);
        vTs = _mm_loadu_ps(&p_Ts[0]);
        vM  = _mm_loadu_ps(&p_M[0]);
        vN0 = _mm_loadu_ps(&p_N0[0]);
        _mm_storeu_ps(&p_sep_mam[0], sse_kernel(vAc,vTs,vM,vN0,n));
        return (0);
    }
    else if(a_sz>4ull && a_sz<=48ull)
    {
        std::size_t i,j;
        for(i = 0ull; i!=ROUND_DOWN(a_sz,3ull); i += 4ull)
        {
            vAc = _mm_loadu_ps(&p_Ac[i]);
            vTs = _mm_loadu_ps(&p_Ts[i]);
            vM  = _mm_loadu_ps(&p_M[i]);
            vN0 = _mm_loadu_ps(&p_N0[i]);
            _mm_storeu_ps(&p_sep_mam[i], sse_kernel(vAc,vTs,vM,vN0,n));
        }
        for(j = i; j<a_sz; ++j)  
        {
            const float valAc = p_Ac[j];
            const float valTs = p_Ts[j];
            const float valM  = p_M[j];
            const float valN0 = p_N0[j];
            p_sep_mam[j]      = scalar_kernel(valAc,valTs,valM,valN0,n);
        }
        return (0);
    }
    else if(a_sz>48ull)
    {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_USE_PEEL_LOOP) == 1
       while(((std::uintptr_t)&p_sep_mam & 15ull) && a_sz)
       {
           const float valAc = *p_Ac;
           const float valTs = *p_Ts;
           const float valM  = *p_M;
           const float valN0 = *p_N0;
           *p_sep_mam        =  scalar_kernel(valAc,valTs,valM,valN0,n);
           p_Ac++;
           p_Ts++;
           p_M++;
           p_N0++;
           p_sep_mam++;
           a_sz--;
       }
#else 
#error "***[FATAL]*** -- The loop peeling procedure must be always equal to 1!!"
#endif
       __m128 xmm0;
       __m128 xmm1;
       __m128 xmm2;
       __m128 xmm3;
       __m128 xmm4;
       __m128 xmm5;
       __m128 xmm6;
       __m128 xmm7;
       __m128 xmm8;
       __m128 xmm9;
       __m128 xmm10;
       __m128 xmm11;
       __m128 xmm12;
       __m128 xmm13;
       __m128 xmm14;
       std::size_t i,j;
       constexpr std::size_t L2_offset = 8000ull;
       constexpr std::size_t L1D_boundary = L2_offset;
    if(a_sz <= L1D_boundary)
    {
       for(i = 0ull; (i+47ull) < a_sz; i += 48ull) 
       {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
#endif 
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+8ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+8ull]);
            xmm12 = _mm_load_ps(&p_M[i+8ull]);
            xmm13 = _mm_load_ps(&p_N0[i+8ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+8ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+12ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+12ull]);
            xmm2 = _mm_load_ps(&p_M[i+12ull]);
            xmm3 = _mm_load_ps(&p_N0[i+12ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+12ull],xmm4);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+16ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+16ull],_MM_HINT_T0);
#endif
            xmm5 = _mm_load_ps(&p_Ac[i+16ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+16ull]);
            xmm7 = _mm_load_ps(&p_M[i+16ull]);
            xmm8 = _mm_load_ps(&p_N0[i+16ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+16ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+20ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+20ull]);
            xmm12 = _mm_load_ps(&p_M[i+20ull]);
            xmm13 = _mm_load_ps(&p_N0[i+20ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+20ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+24ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+24ull]);
            xmm2 = _mm_load_ps(&p_M[i+24ull]);
            xmm3 = _mm_load_ps(&p_N0[i+24ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+24ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+28ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+28ull]);
            xmm7 = _mm_load_ps(&p_M[i+28ull]);
            xmm8 = _mm_load_ps(&p_N0[i+28ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+28ull],xmm9);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+32ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+32ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+32ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+32ull],_MM_HINT_T0);
#endif

            xmm10 = _mm_load_ps(&p_Ac[i+32ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+32ull]);
            xmm12 = _mm_load_ps(&p_M[i+32ull]);
            xmm13 = _mm_load_ps(&p_N0[i+32ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+32ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+36ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+36ull]);
            xmm2 = _mm_load_ps(&p_M[i+36ull]);
            xmm3 = _mm_load_ps(&p_N0[i+36ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+36ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+40ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+40ull]);
            xmm7 = _mm_load_ps(&p_M[i+40ull]);
            xmm8 = _mm_load_ps(&p_N0[i+40ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+40ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+44ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+44ull]);
            xmm12 = _mm_load_ps(&p_M[i+44ull]);
            xmm13 = _mm_load_ps(&p_N0[i+44ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+44ull],xmm14);
       }
 
       for(; (i+31ull) < a_sz; i += 32ull) 
       {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
#endif 
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+8ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+8ull]);
            xmm12 = _mm_load_ps(&p_M[i+8ull]);
            xmm13 = _mm_load_ps(&p_N0[i+8ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+8ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+12ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+12ull]);
            xmm2 = _mm_load_ps(&p_M[i+12ull]);
            xmm3 = _mm_load_ps(&p_N0[i+12ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+12ull],xmm4);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+16ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+16ull],_MM_HINT_T0);
#endif
            xmm5 = _mm_load_ps(&p_Ac[i+16ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+16ull]);
            xmm7 = _mm_load_ps(&p_M[i+16ull]);
            xmm8 = _mm_load_ps(&p_N0[i+16ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+16ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+20ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+20ull]);
            xmm12 = _mm_load_ps(&p_M[i+20ull]);
            xmm13 = _mm_load_ps(&p_N0[i+20ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+20ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+24ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+24ull]);
            xmm2 = _mm_load_ps(&p_M[i+24ull]);
            xmm3 = _mm_load_ps(&p_N0[i+24ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+24ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+28ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+28ull]);
            xmm7 = _mm_load_ps(&p_M[i+28ull]);
            xmm8 = _mm_load_ps(&p_N0[i+28ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+28ull],xmm9);
       }

       for(; (i+19ull) < a_sz; i += 20ull)
       {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
#endif 
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+8ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+8ull]);
            xmm12 = _mm_load_ps(&p_M[i+8ull]);
            xmm13 = _mm_load_ps(&p_N0[i+8ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+8ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+12ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+12ull]);
            xmm2 = _mm_load_ps(&p_M[i+12ull]);
            xmm3 = _mm_load_ps(&p_N0[i+12ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+12ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+16ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+16ull]);
            xmm7 = _mm_load_ps(&p_M[i+16ull]);
            xmm8 = _mm_load_ps(&p_N0[i+16ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+16ull],xmm9);
       }

       for(; (i+7ull) < a_sz; i += 8ull)
       {
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);
       }

       for(; (i+3ull) < a_sz; i += 4ull)
       {
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);
       }

       for(; (i+0ull) < a_sz; i += 1ull)
       {
            const float valAc = p_Ac[i];
            const float valTs = p_Ts[i];
            const float valM  = p_M[i];
            const float valN0 = p_N0[i];
            p_sep_mam[i]      = scalar_kernel(valAc,valTs,valM,valN0,n);
       }
       return (0);
    }
    else 
    {
       for(i = 0ull; (i+47ull) < a_sz; i += 48ull) 
       {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+0ull], _MM_HIN1_T0);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+0ull],_MM_HINT_T1);
#endif 
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+8ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+8ull]);
            xmm12 = _mm_load_ps(&p_M[i+8ull]);
            xmm13 = _mm_load_ps(&p_N0[i+8ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+8ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+12ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+12ull]);
            xmm2 = _mm_load_ps(&p_M[i+12ull]);
            xmm3 = _mm_load_ps(&p_N0[i+12ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+12ull],xmm4);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+16ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+16ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+16ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+16ull], _MM_HINT_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+16ull],_MM_HINT_T1);
#endif
            xmm5 = _mm_load_ps(&p_Ac[i+16ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+16ull]);
            xmm7 = _mm_load_ps(&p_M[i+16ull]);
            xmm8 = _mm_load_ps(&p_N0[i+16ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+16ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+20ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+20ull]);
            xmm12 = _mm_load_ps(&p_M[i+20ull]);
            xmm13 = _mm_load_ps(&p_N0[i+20ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+20ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+24ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+24ull]);
            xmm2 = _mm_load_ps(&p_M[i+24ull]);
            xmm3 = _mm_load_ps(&p_N0[i+24ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+24ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+28ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+28ull]);
            xmm7 = _mm_load_ps(&p_M[i+28ull]);
            xmm8 = _mm_load_ps(&p_N0[i+28ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+28ull],xmm9);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+32ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+32ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+32ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+32ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+32ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+32ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+32ull], _MM_HINT_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+32ull],_MM_HINT_T1);
#endif

            xmm10 = _mm_load_ps(&p_Ac[i+32ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+32ull]);
            xmm12 = _mm_load_ps(&p_M[i+32ull]);
            xmm13 = _mm_load_ps(&p_N0[i+32ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+32ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+36ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+36ull]);
            xmm2 = _mm_load_ps(&p_M[i+36ull]);
            xmm3 = _mm_load_ps(&p_N0[i+36ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+36ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+40ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+40ull]);
            xmm7 = _mm_load_ps(&p_M[i+40ull]);
            xmm8 = _mm_load_ps(&p_N0[i+40ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+40ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+44ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+44ull]);
            xmm12 = _mm_load_ps(&p_M[i+44ull]);
            xmm13 = _mm_load_ps(&p_N0[i+44ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+44ull],xmm14);
       }
 
       for(; (i+31ull) < a_sz; i += 32ull) 
       {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+0ull], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+0ull],_MM_HINT_T1);
#endif 
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+8ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+8ull]);
            xmm12 = _mm_load_ps(&p_M[i+8ull]);
            xmm13 = _mm_load_ps(&p_N0[i+8ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+8ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+12ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+12ull]);
            xmm2 = _mm_load_ps(&p_M[i+12ull]);
            xmm3 = _mm_load_ps(&p_N0[i+12ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+12ull],xmm4);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+16ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+16ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+16ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+16ull], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+16ull],_MM_HINT_T1);
#endif
            xmm5 = _mm_load_ps(&p_Ac[i+16ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+16ull]);
            xmm7 = _mm_load_ps(&p_M[i+16ull]);
            xmm8 = _mm_load_ps(&p_N0[i+16ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+16ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+20ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+20ull]);
            xmm12 = _mm_load_ps(&p_M[i+20ull]);
            xmm13 = _mm_load_ps(&p_N0[i+20ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+20ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+24ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+24ull]);
            xmm2 = _mm_load_ps(&p_M[i+24ull]);
            xmm3 = _mm_load_ps(&p_N0[i+24ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+24ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+28ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+28ull]);
            xmm7 = _mm_load_ps(&p_M[i+28ull]);
            xmm8 = _mm_load_ps(&p_N0[i+28ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+28ull],xmm9);
       }

       for(; (i+19ull) < a_sz; i += 20ull)
       {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+0ull], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+0ull],_MM_HINT_T1);
#endif 
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+8ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+8ull]);
            xmm12 = _mm_load_ps(&p_M[i+8ull]);
            xmm13 = _mm_load_ps(&p_N0[i+8ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+8ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+12ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+12ull]);
            xmm2 = _mm_load_ps(&p_M[i+12ull]);
            xmm3 = _mm_load_ps(&p_N0[i+12ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+12ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+16ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+16ull]);
            xmm7 = _mm_load_ps(&p_M[i+16ull]);
            xmm8 = _mm_load_ps(&p_N0[i+16ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+16ull],xmm9);
       }

       for(; (i+7ull) < a_sz; i += 8ull)
       {
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);
       }

       for(; (i+3ull) < a_sz; i += 4ull)
       {
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);
       }

       for(; (i+0ull) < a_sz; i += 1ull)
       {
            const float valAc = p_Ac[i];
            const float valTs = p_Ts[i];
            const float valM  = p_M[i];
            const float valN0 = p_N0[i];
            p_sep_mam[i]      = scalar_kernel(valAc,valTs,valM,valN0,n);
       }
       return (0);
    }

    }
}

///////////////////////////////////////////////////////////////////////////////////

#if (ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F32_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=sse
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<typename callable_sse_kernel,
         typename callable_scalar_kernel>
__ATTR_ALWAYS_INLINE__ 
static inline 
std::int32_t 
analytic_SEP_BEP_generic_4xf32_array1d_u8x(const float * __restrict__ Ac, 
                                        const float * __restrict__ Ts,
                                        const float * __restrict__ M,
                                        const float * __restrict__ N0,
                                        float       * __restrict__ sep_mam,
                                        const std::int32_t n,
                                        std::size_t sz,
                                        const callable_sse_kernel     & sse_kernel,
                                        const callable_scalar_kernel  & scalar_kernel) 
{
    if(__builtin_expect(sz<1ull,0)) {return (-1);}
    __m128 vAc;
    __m128 vTs;
    __m128 vM;
    __m128 vN0;
    const float * __restrict__ p_Ac      = Ac;
    const float * __restrict__ p_Ts      = Ts;
    const float * __restrict__ p_M       = M;
    const float * __restrict__ p_N0      = N0;
    float       * __restrict__ p_sep_mam = sep_mam;
    std::size_t             a_sz         = sz;
    if(a_sz>=1ull && a_sz<4ull)
    {
       for(std::uint32_t i = 0; i<4; ++i) 
       {
           const float valAc = p_Ac[i];
           const float valTs = p_Ts[i];
           const float valM  = p_M[i];
           const float valN0 = p_N0[i];
           p_sep_mam[i]      = scalar_kernel(valAc,valTs,valM,valN0,n);
       }
       return (0);
    }
    else if(a_sz==4ull)
    {
        vAc = _mm_loadu_ps(&p_Ac[0]);
        vTs = _mm_loadu_ps(&p_Ts[0]);
        vM  = _mm_loadu_ps(&p_M[0]);
        vN0 = _mm_loadu_ps(&p_N0[0]);
        _mm_storeu_ps(&p_sep_mam[0], sse_kernel(vAc,vTs,vM,vN0,n));
        return (0);
    }
    else if(a_sz>4ull && a_sz<=32ull)
    {
        std::size_t i,j;
        for(i = 0ull; i!=ROUND_DOWN(a_sz,3ull); i += 4ull)
        {
            vAc = _mm_loadu_ps(&p_Ac[i]);
            vTs = _mm_loadu_ps(&p_Ts[i]);
            vM  = _mm_loadu_ps(&p_M[i]);
            vN0 = _mm_loadu_ps(&p_N0[i]);
            _mm_storeu_ps(&p_sep_mam[i], sse_kernel(vAc,vTs,vM,vN0,n));
        }
        for(j = i; j<a_sz; ++j)  
        {
            const float valAc = p_Ac[j];
            const float valTs = p_Ts[j];
            const float valM  = p_M[j];
            const float valN0 = p_N0[j];
            p_sep_mam[j]      = scalar_kernel(valAc,valTs,valM,valN0,n);
        }
        return (0);
    }
    else if(a_sz>32ull)
    {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_USE_PEEL_LOOP) == 1
       while(((std::uintptr_t)&p_sep_mam & 15ull) && a_sz)
       {
           const float valAc = *p_Ac;
           const float valTs = *p_Ts;
           const float valM  = *p_M;
           const float valN0 = *p_N0;
           *p_sep_mam        =  scalar_kernel(valAc,valTs,valM,valN0,n);
           p_Ac++;
           p_Ts++;
           p_M++;
           p_N0++;
           p_sep_mam++;
           a_sz--;
       }
#else 
#error "***[FATAL]*** -- The loop peeling procedure must be always equal to 1!!"
#endif
       __m128 xmm0;
       __m128 xmm1;
       __m128 xmm2;
       __m128 xmm3;
       __m128 xmm4;
       __m128 xmm5;
       __m128 xmm6;
       __m128 xmm7;
       __m128 xmm8;
       __m128 xmm9;
       __m128 xmm10;
       __m128 xmm11;
       __m128 xmm12;
       __m128 xmm13;
       __m128 xmm14;
       std::size_t i,j;
       constexpr std::size_t L2_offset = 8000ull;
       constexpr std::size_t L1D_boundary = L2_offset;
    if(a_sz <= L1D_boundary)
    {
       for(i = 0ull; (i+31ull) < a_sz; i += 32ull) 
       {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
#endif 
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+8ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+8ull]);
            xmm12 = _mm_load_ps(&p_M[i+8ull]);
            xmm13 = _mm_load_ps(&p_N0[i+8ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+8ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+12ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+12ull]);
            xmm2 = _mm_load_ps(&p_M[i+12ull]);
            xmm3 = _mm_load_ps(&p_N0[i+12ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+12ull],xmm4);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+16ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+16ull],_MM_HINT_T0);
#endif
            xmm5 = _mm_load_ps(&p_Ac[i+16ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+16ull]);
            xmm7 = _mm_load_ps(&p_M[i+16ull]);
            xmm8 = _mm_load_ps(&p_N0[i+16ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+16ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+20ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+20ull]);
            xmm12 = _mm_load_ps(&p_M[i+20ull]);
            xmm13 = _mm_load_ps(&p_N0[i+20ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+20ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+24ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+24ull]);
            xmm2 = _mm_load_ps(&p_M[i+24ull]);
            xmm3 = _mm_load_ps(&p_N0[i+24ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+24ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+28ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+28ull]);
            xmm7 = _mm_load_ps(&p_M[i+28ull]);
            xmm8 = _mm_load_ps(&p_N0[i+28ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+28ull],xmm9);
       }
 
       for(; (i+19ull) < a_sz; i += 20ull)
       {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
#endif 
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+8ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+8ull]);
            xmm12 = _mm_load_ps(&p_M[i+8ull]);
            xmm13 = _mm_load_ps(&p_N0[i+8ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+8ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+12ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+12ull]);
            xmm2 = _mm_load_ps(&p_M[i+12ull]);
            xmm3 = _mm_load_ps(&p_N0[i+12ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+12ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+16ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+16ull]);
            xmm7 = _mm_load_ps(&p_M[i+16ull]);
            xmm8 = _mm_load_ps(&p_N0[i+16ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+16ull],xmm9);
       }

       for(; (i+7ull) < a_sz; i += 8ull)
       {
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);
       }

       for(; (i+3ull) < a_sz; i += 4ull)
       {
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);
       }

       for(; (i+0ull) < a_sz; i += 1ull)
       {
            const float valAc = p_Ac[i];
            const float valTs = p_Ts[i];
            const float valM  = p_M[i];
            const float valN0 = p_N0[i];
            p_sep_mam[i]      = scalar_kernel(valAc,valTs,valM,valN0,n);
       }
       return (0);
    }
    else 
    {
       for(i = 0ull; (i+31ull) < a_sz; i += 32ull) 
       {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+0ull], _MM_HIN1_T0);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+0ull],_MM_HINT_T1);
#endif 
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+8ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+8ull]);
            xmm12 = _mm_load_ps(&p_M[i+8ull]);
            xmm13 = _mm_load_ps(&p_N0[i+8ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+8ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+12ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+12ull]);
            xmm2 = _mm_load_ps(&p_M[i+12ull]);
            xmm3 = _mm_load_ps(&p_N0[i+12ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+12ull],xmm4);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+16ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+16ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+16ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+16ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+16ull], _MM_HINT_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+16ull],_MM_HINT_T1);
#endif
            xmm5 = _mm_load_ps(&p_Ac[i+16ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+16ull]);
            xmm7 = _mm_load_ps(&p_M[i+16ull]);
            xmm8 = _mm_load_ps(&p_N0[i+16ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+16ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+20ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+20ull]);
            xmm12 = _mm_load_ps(&p_M[i+20ull]);
            xmm13 = _mm_load_ps(&p_N0[i+20ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+20ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+24ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+24ull]);
            xmm2 = _mm_load_ps(&p_M[i+24ull]);
            xmm3 = _mm_load_ps(&p_N0[i+24ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+24ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+28ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+28ull]);
            xmm7 = _mm_load_ps(&p_M[i+28ull]);
            xmm8 = _mm_load_ps(&p_N0[i+28ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+28ull],xmm9);
       }
 
       for(; (i+19ull) < a_sz; i += 20ull)
       {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+0ull], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+0ull],_MM_HINT_T1);
#endif 
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+8ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+8ull]);
            xmm12 = _mm_load_ps(&p_M[i+8ull]);
            xmm13 = _mm_load_ps(&p_N0[i+8ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+8ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+12ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+12ull]);
            xmm2 = _mm_load_ps(&p_M[i+12ull]);
            xmm3 = _mm_load_ps(&p_N0[i+12ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+12ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+16ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+16ull]);
            xmm7 = _mm_load_ps(&p_M[i+16ull]);
            xmm8 = _mm_load_ps(&p_N0[i+16ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+16ull],xmm9);
       }

       for(; (i+7ull) < a_sz; i += 8ull)
       {
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);
       }

       for(; (i+3ull) < a_sz; i += 4ull)
       {
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);
       }

       for(; (i+0ull) < a_sz; i += 1ull)
       {
            const float valAc = p_Ac[i];
            const float valTs = p_Ts[i];
            const float valM  = p_M[i];
            const float valN0 = p_N0[i];
            p_sep_mam[i]      = scalar_kernel(valAc,valTs,valM,valN0,n);
       }
       return (0);
    }

    }
}

/////////////////////////////////////////////////////////////////////////////////

#if (ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F32_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=sse
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<typename callable_sse_kernel,
         typename callable_scalar_kernel>
__ATTR_ALWAYS_INLINE__ 
static inline 
std::int32_t 
analytic_SEP_BEP_generic_4xf32_array1d_u4x(const float * __restrict__ Ac, 
                                        const float * __restrict__ Ts,
                                        const float * __restrict__ M,
                                        const float * __restrict__ N0,
                                        float       * __restrict__ sep_mam,
                                        const std::int32_t n,
                                        std::size_t sz,
                                        const callable_sse_kernel     & sse_kernel,
                                        const callable_scalar_kernel  & scalar_kernel) 
{
    if(__builtin_expect(sz<1ull,0)) {return (-1);}
    __m128 vAc;
    __m128 vTs;
    __m128 vM;
    __m128 vN0;
    const float * __restrict__ p_Ac      = Ac;
    const float * __restrict__ p_Ts      = Ts;
    const float * __restrict__ p_M       = M;
    const float * __restrict__ p_N0      = N0;
    float       * __restrict__ p_sep_mam = sep_mam;
    std::size_t             a_sz         = sz;
    if(a_sz>=1ull && a_sz<4ull)
    {
       for(std::uint32_t i = 0; i<4; ++i) 
       {
           const float valAc = p_Ac[i];
           const float valTs = p_Ts[i];
           const float valM  = p_M[i];
           const float valN0 = p_N0[i];
           p_sep_mam[i]      = scalar_kernel(valAc,valTs,valM,valN0,n);
       }
       return (0);
    }
    else if(a_sz==4ull)
    {
        vAc = _mm_loadu_ps(&p_Ac[0]);
        vTs = _mm_loadu_ps(&p_Ts[0]);
        vM  = _mm_loadu_ps(&p_M[0]);
        vN0 = _mm_loadu_ps(&p_N0[0]);
        _mm_storeu_ps(&p_sep_mam[0], sse_kernel(vAc,vTs,vM,vN0,n));
        return (0);
    }
    else if(a_sz>4ull && a_sz<=16ull)
    {
        std::size_t i,j;
        for(i = 0ull; i!=ROUND_DOWN(a_sz,3ull); i += 4ull)
        {
            vAc = _mm_loadu_ps(&p_Ac[i]);
            vTs = _mm_loadu_ps(&p_Ts[i]);
            vM  = _mm_loadu_ps(&p_M[i]);
            vN0 = _mm_loadu_ps(&p_N0[i]);
            _mm_storeu_ps(&p_sep_mam[i], sse_kernel(vAc,vTs,vM,vN0,n));
        }
        for(j = i; j<a_sz; ++j)  
        {
            const float valAc = p_Ac[j];
            const float valTs = p_Ts[j];
            const float valM  = p_M[j];
            const float valN0 = p_N0[j];
            p_sep_mam[j]      = scalar_kernel(valAc,valTs,valM,valN0,n);
        }
        return (0);
    }
    else if(a_sz>16ull)
    {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_USE_PEEL_LOOP) == 1
       while(((std::uintptr_t)&p_sep_mam & 15ull) && a_sz)
       {
           const float valAc = *p_Ac;
           const float valTs = *p_Ts;
           const float valM  = *p_M;
           const float valN0 = *p_N0;
           *p_sep_mam        =  scalar_kernel(valAc,valTs,valM,valN0,n);
           p_Ac++;
           p_Ts++;
           p_M++;
           p_N0++;
           p_sep_mam++;
           a_sz--;
       }
#else 
#error "***[FATAL]*** -- The loop peeling procedure must be always equal to 1!!"
#endif
       __m128 xmm0;
       __m128 xmm1;
       __m128 xmm2;
       __m128 xmm3;
       __m128 xmm4;
       __m128 xmm5;
       __m128 xmm6;
       __m128 xmm7;
       __m128 xmm8;
       __m128 xmm9;
       __m128 xmm10;
       __m128 xmm11;
       __m128 xmm12;
       __m128 xmm13;
       __m128 xmm14;
       std::size_t i,j;
       constexpr std::size_t L2_offset = 8000ull;
       constexpr std::size_t L1D_boundary = L2_offset;
    if(a_sz <= L1D_boundary)
    {
       for(i = 0ull; (i+15ull) < a_sz; i += 16ull) 
       {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
#endif 
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+8ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+8ull]);
            xmm12 = _mm_load_ps(&p_M[i+8ull]);
            xmm13 = _mm_load_ps(&p_N0[i+8ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+8ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+12ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+12ull]);
            xmm2 = _mm_load_ps(&p_M[i+12ull]);
            xmm3 = _mm_load_ps(&p_N0[i+12ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+12ull],xmm4);

       }
 
       for(; (i+7ull) < a_sz; i += 8ull)
       {
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);
       }

       for(; (i+3ull) < a_sz; i += 4ull)
       {
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);
       }

       for(; (i+0ull) < a_sz; i += 1ull)
       {
            const float valAc = p_Ac[i];
            const float valTs = p_Ts[i];
            const float valM  = p_M[i];
            const float valN0 = p_N0[i];
            p_sep_mam[i]      = scalar_kernel(valAc,valTs,valM,valN0,n);
       }
       return (0);
    }
    else 
    {
       for(i = 0ull; (i+15ull) < a_sz; i += 16ull) 
       {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+0ull], _MM_HIN1_T0);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+0ull],_MM_HINT_T1);
#endif 
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);

            xmm10 = _mm_load_ps(&p_Ac[i+8ull]);
            xmm11 = _mm_load_ps(&p_Ts[i+8ull]);
            xmm12 = _mm_load_ps(&p_M[i+8ull]);
            xmm13 = _mm_load_ps(&p_N0[i+8ull]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_ps(&p_sep_mam[i+8ull],xmm14);

            xmm0 = _mm_load_ps(&p_Ac[i+12ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+12ull]);
            xmm2 = _mm_load_ps(&p_M[i+12ull]);
            xmm3 = _mm_load_ps(&p_N0[i+12ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+12ull],xmm4);
       }

       for(; (i+7ull) < a_sz; i += 8ull)
       {
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);
       }

       for(; (i+3ull) < a_sz; i += 4ull)
       {
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);
       }

       for(; (i+0ull) < a_sz; i += 1ull)
       {
            const float valAc = p_Ac[i];
            const float valTs = p_Ts[i];
            const float valM  = p_M[i];
            const float valN0 = p_N0[i];
            p_sep_mam[i]      = scalar_kernel(valAc,valTs,valM,valN0,n);
       }
       return (0);
    }

    }
}

//////////////////////////////////////////////////////////////////////////////////

#if (ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F32_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=sse
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
#endif 
template<typename callable_sse_kernel,
         typename callable_scalar_kernel>
__ATTR_ALWAYS_INLINE__ 
static inline 
std::int32_t 
analytic_SEP_BEP_generic_4xf32_array1d_u2x(const float * __restrict__ Ac, 
                                        const float * __restrict__ Ts,
                                        const float * __restrict__ M,
                                        const float * __restrict__ N0,
                                        float       * __restrict__ sep_mam,
                                        const std::int32_t n,
                                        std::size_t sz,
                                        const callable_sse_kernel     & sse_kernel,
                                        const callable_scalar_kernel  & scalar_kernel) 
{
    if(__builtin_expect(sz<1ull,0)) {return (-1);}
    __m128 vAc;
    __m128 vTs;
    __m128 vM;
    __m128 vN0;
    const float * __restrict__ p_Ac      = Ac;
    const float * __restrict__ p_Ts      = Ts;
    const float * __restrict__ p_M       = M;
    const float * __restrict__ p_N0      = N0;
    float       * __restrict__ p_sep_mam = sep_mam;
    std::size_t             a_sz         = sz;
    if(a_sz>=1ull && a_sz<4ull)
    {
       for(std::uint32_t i = 0; i<4; ++i) 
       {
           const float valAc = p_Ac[i];
           const float valTs = p_Ts[i];
           const float valM  = p_M[i];
           const float valN0 = p_N0[i];
           p_sep_mam[i]      = scalar_kernel(valAc,valTs,valM,valN0,n);
       }
       return (0);
    }
    else if(a_sz==4ull)
    {
        vAc = _mm_loadu_ps(&p_Ac[0]);
        vTs = _mm_loadu_ps(&p_Ts[0]);
        vM  = _mm_loadu_ps(&p_M[0]);
        vN0 = _mm_loadu_ps(&p_N0[0]);
        _mm_storeu_ps(&p_sep_mam[0], sse_kernel(vAc,vTs,vM,vN0,n));
        return (0);
    }
    else if(a_sz>4ull && a_sz<=8ull)
    {
        std::size_t i,j;
        for(i = 0ull; i!=ROUND_DOWN(a_sz,3ull); i += 4ull)
        {
            vAc = _mm_loadu_ps(&p_Ac[i]);
            vTs = _mm_loadu_ps(&p_Ts[i]);
            vM  = _mm_loadu_ps(&p_M[i]);
            vN0 = _mm_loadu_ps(&p_N0[i]);
            _mm_storeu_ps(&p_sep_mam[i], sse_kernel(vAc,vTs,vM,vN0,n));
        }
        for(j = i; j<a_sz; ++j)  
        {
            const float valAc = p_Ac[j];
            const float valTs = p_Ts[j];
            const float valM  = p_M[j];
            const float valN0 = p_N0[j];
            p_sep_mam[j]      = scalar_kernel(valAc,valTs,valM,valN0,n);
        }
        return (0);
    }
    else if(a_sz>8ull)
    {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_USE_PEEL_LOOP) == 1
       while(((std::uintptr_t)&p_sep_mam & 15ull) && a_sz)
       {
           const float valAc = *p_Ac;
           const float valTs = *p_Ts;
           const float valM  = *p_M;
           const float valN0 = *p_N0;
           *p_sep_mam        =  scalar_kernel(valAc,valTs,valM,valN0,n);
           p_Ac++;
           p_Ts++;
           p_M++;
           p_N0++;
           p_sep_mam++;
           a_sz--;
       }
#else 
#error "***[FATAL]*** -- The loop peeling procedure must be always equal to 1!!"
#endif
       __m128 xmm0;
       __m128 xmm1;
       __m128 xmm2;
       __m128 xmm3;
       __m128 xmm4;
       __m128 xmm5;
       __m128 xmm6;
       __m128 xmm7;
       __m128 xmm8;
       __m128 xmm9;
       std::size_t i,j;
       constexpr std::size_t L2_offset = 8000ull;
       constexpr std::size_t L1D_boundary = L2_offset;
    if(a_sz <= L1D_boundary)
    {
       for(i = 0ull; (i+7ull) < a_sz; i += 8ull) 
       {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
#endif 
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);

       }
 
       for(; (i+3ull) < a_sz; i += 4ull)
       {
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);
       }

       for(; (i+0ull) < a_sz; i += 1ull)
       {
            const float valAc = p_Ac[i];
            const float valTs = p_Ts[i];
            const float valM  = p_M[i];
            const float valN0 = p_N0[i];
            p_sep_mam[i]      = scalar_kernel(valAc,valTs,valM,valN0,n);
       }
       return (0);
    }
    else 
    {
       for(i = 0ull; (i+7ull) < a_sz; i += 8ull) 
       {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+0ull], _MM_HIN1_T0);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+0ull],_MM_HINT_T1);
#endif 
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_ps(&p_Ac[i+4ull]);
            xmm6 = _mm_load_ps(&p_Ts[i+4ull]);
            xmm7 = _mm_load_ps(&p_M[i+4ull]);
            xmm8 = _mm_load_ps(&p_N0[i+4ull]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_ps(&p_sep_mam[i+4ull],xmm9);
       }

       for(; (i+3ull) < a_sz; i += 4ull)
       {
            xmm0 = _mm_load_ps(&p_Ac[i+0ull]);
            xmm1 = _mm_load_ps(&p_Ts[i+0ull]);
            xmm2 = _mm_load_ps(&p_M[i+0ull]);
            xmm3 = _mm_load_ps(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_ps(&p_sep_mam[i+0ull],xmm4);
       }

       for(; (i+0ull) < a_sz; i += 1ull)
       {
            const float valAc = p_Ac[i];
            const float valTs = p_Ts[i];
            const float valM  = p_M[i];
            const float valN0 = p_N0[i];
            p_sep_mam[i]      = scalar_kernel(valAc,valTs,valM,valN0,n);
       }
       return (0);
    }

    }
}

}// fading_channel

}// gms
#endif /*__GMS_ANALYTIC_BEP_SEP_CH8_ARRAY1D_F32_H__*/