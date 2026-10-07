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

#ifndef __GMS_ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_H__
#define __GMS_ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_H__ 061020260556

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_analytic_bep_sep_ch8.h"
#include "GMS_analytic_bep_sep_ch8_sse.h"



namespace file_info 
{

     static const unsigned int GMS_ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F64_MAJOR = 1;
     static const unsigned int GMS_ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F64_MINOR = 1;
     static const unsigned int GMS_ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F64_MICRO = 0;
     static const unsigned int GMS_ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F64_FULLVER =
       1000U*GMS_ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F64_MAJOR+100U*GMS_ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F64_MINOR+
       10U*GMS_ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F64_MICRO;
     static const char GMS_ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F64_CREATION_DATE[] = "06-10-2026 05:55PM +00200 (TUE 06 OCT 2026 GMT+2)";
     static const char GMS_ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F64_BUILD_DATE[]    = __DATE__; 
     static const char GMS_ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F64_BUILD_TIME[]    = __TIME__;
     static const char GMS_ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F64_AUTHOR[]        = "Programmer: Bernard Gingold, beniekg@gmail.com";
     static const char GMS_ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F64_SYNOPSIS[]      = "Analytic Formulae of the Bit and Symbol Error Probabilities based on the approximated Gaussian-Q function.\
	                                                          Based on the chapter 8 of M.K Simon, M.S. Alouini: Digital Communication over Fading Channels 1st ed\
															  ISBN-13 978-0471317791";

}

#if !defined(ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F64_OVERRIDE_COMPILER_CMD_LINE)
#define ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F64_OVERRIDE_COMPILER_CMD_LINE 0
#endif

#if !defined(ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_USE_PEEL_LOOP)
#define ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_USE_PEEL_LOOP 1
#endif 

#if !defined(ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH)
#define ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH 0
#endif 
 
// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F64_USE_PMC_INSTRUMENTATION)
#define ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F64_USE_PMC_INSTRUMENTATION 0
#endif 

#if (ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F64_USE_PMC_INSTRUMENTATION) == 1
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

#if (ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F64_OVERRIDE_COMPILER_CMD_LINE) == 1
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
analytic_SEP_BEP_generic_2xf64_array1d_u16x(const double * __restrict__ Ac, 
                                        const double * __restrict__ Ts,
                                        const double * __restrict__ M,
                                        const double * __restrict__ N0,
                                        double       * __restrict__ sep_mam,
                                        const std::int32_t n,
                                        std::size_t sz,
                                        const callable_sse_kernel     & sse_kernel,
                                        const callable_scalar_kernel  & scalar_kernel) 
{
    if(__builtin_expect(sz<1ull,0)) {return (-1);}
    __m128d vAc;
    __m128d vTs;
    __m128d vM;
    __m128d vN0;
    const double * __restrict__ p_Ac      = Ac;
    const double * __restrict__ p_Ts      = Ts;
    const double * __restrict__ p_M       = M;
    const double * __restrict__ p_N0      = N0;
    double       * __restrict__ p_sep_mam = sep_mam;
    std::size_t             a_sz         = sz;
    if(a_sz>=1ull && a_sz<2ull)
    {
       for(std::uint32_t i = 0; i<2; ++i) 
       {
           const double valAc = p_Ac[i];
           const double valTs = p_Ts[i];
           const double valM  = p_M[i];
           const double valN0 = p_N0[i];
           p_sep_mam[i]      = scalar_kernel(valAc,valTs,valM,valN0,n);
       }
       return (0);
    }
    else if(a_sz==2ull)
    {
        vAc = _mm_loadu_pd(&p_Ac[0]);
        vTs = _mm_loadu_pd(&p_Ts[0]);
        vM  = _mm_loadu_pd(&p_M[0]);
        vN0 = _mm_loadu_pd(&p_N0[0]);
        _mm_storeu_pd(&p_sep_mam[0], sse_kernel(vAc,vTs,vM,vN0,n));
        return (0);
    }
    else if(a_sz>2ull && a_sz<=32ull)
    {
        std::size_t i,j;
        for(i = 0ull; i!=ROUND_DOWN(a_sz,1ull); i += 2ull)
        {
            vAc = _mm_loadu_pd(&p_Ac[i]);
            vTs = _mm_loadu_pd(&p_Ts[i]);
            vM  = _mm_loadu_pd(&p_M[i]);
            vN0 = _mm_loadu_pd(&p_N0[i]);
            _mm_storeu_pd(&p_sep_mam[i], sse_kernel(vAc,vTs,vM,vN0,n));
        }
        for(j = i; j<a_sz; ++j)  
        {
            const double valAc = p_Ac[j];
            const double valTs = p_Ts[j];
            const double valM  = p_M[j];
            const double valN0 = p_N0[j];
            p_sep_mam[j]      = scalar_kernel(valAc,valTs,valM,valN0,n);
        }
        return (0);
    }
    else if(a_sz>32ull)
    {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_USE_PEEL_LOOP) == 1
       while(((std::uintptr_t)&p_sep_mam & 15ull) && a_sz)
       {
           const double valAc = *p_Ac;
           const double valTs = *p_Ts;
           const double valM  = *p_M;
           const double valN0 = *p_N0;
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
       __m128d xmm0;
       __m128d xmm1;
       __m128d xmm2;
       __m128d xmm3;
       __m128d xmm4;
       __m128d xmm5;
       __m128d xmm6;
       __m128d xmm7;
       __m128d xmm8;
       __m128d xmm9;
       __m128d xmm10;
       __m128d xmm11;
       __m128d xmm12;
       __m128d xmm13;
       __m128d xmm14;
       std::size_t i,j;
       constexpr std::size_t L2_offset = 4000ull;
       constexpr std::size_t L1D_boundary = L2_offset;
    if(a_sz <= L1D_boundary)
    {
       for(i = 0ull; (i+31ull) < a_sz; i += 32ull) 
       {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
#endif 
            const std::size_t idx2  = i+2ull;
            const std::size_t idx4  = i+4ull;
            const std::size_t idx6  = i+6ull;
            const std::size_t idx8  = i+8ull;
            const std::size_t idx10 = i+10ull;
            const std::size_t idx12 = i+12ull;
            const std::size_t idx14 = i+14ull;
            const std::size_t idx16 = i+16ull;
            const std::size_t idx18 = i+18ull;
            const std::size_t idx20 = i+20ull;
            const std::size_t idx22 = i+22ull;
            const std::size_t idx24 = i+24ull;
            const std::size_t idx26 = i+26ull;
            const std::size_t idx28 = i+28ull;
            const std::size_t idx30 = i+30ull;

            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx2]);
            xmm6 = _mm_load_pd(&p_Ts[idx2]);
            xmm7 = _mm_load_pd(&p_M[idx2]);
            xmm8 = _mm_load_pd(&p_N0[idx2]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx2],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx4]);
            xmm11 = _mm_load_pd(&p_Ts[idx4]);
            xmm12 = _mm_load_pd(&p_M[idx4]);
            xmm13 = _mm_load_pd(&p_N0[idx4]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx4],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx6]);
            xmm1 = _mm_load_pd(&p_Ts[idx6]);
            xmm2 = _mm_load_pd(&p_M[idx6]);
            xmm3 = _mm_load_pd(&p_N0[idx6]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx6],xmm4);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[idx8],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[idx8],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[idx8], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[idx8],_MM_HINT_T0);
#endif
            xmm5 = _mm_load_pd(&p_Ac[idx8]);
            xmm6 = _mm_load_pd(&p_Ts[idx8]);
            xmm7 = _mm_load_pd(&p_M[idx8]);
            xmm8 = _mm_load_pd(&p_N0[idx8]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx8],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx10]);
            xmm11 = _mm_load_pd(&p_Ts[idx10]);
            xmm12 = _mm_load_pd(&p_M[idx10]);
            xmm13 = _mm_load_pd(&p_N0[idx10]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx10],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx12]);
            xmm1 = _mm_load_pd(&p_Ts[idx12]);
            xmm2 = _mm_load_pd(&p_M[idx12]);
            xmm3 = _mm_load_pd(&p_N0[idx12]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx12],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx14]);
            xmm6 = _mm_load_pd(&p_Ts[idx14]);
            xmm7 = _mm_load_pd(&p_M[idx14]);
            xmm8 = _mm_load_pd(&p_N0[idx14]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx14],xmm9);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[idx16],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[idx16],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[idx16], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[idx16],_MM_HINT_T0);
#endif

            xmm10 = _mm_load_pd(&p_Ac[idx16]);
            xmm11 = _mm_load_pd(&p_Ts[idx16]);
            xmm12 = _mm_load_pd(&p_M[idx16]);
            xmm13 = _mm_load_pd(&p_N0[idx16]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx16],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx18]);
            xmm1 = _mm_load_pd(&p_Ts[idx18]);
            xmm2 = _mm_load_pd(&p_M[idx18]);
            xmm3 = _mm_load_pd(&p_N0[idx18]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx18],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx20]);
            xmm6 = _mm_load_pd(&p_Ts[idx20]);
            xmm7 = _mm_load_pd(&p_M[idx20]);
            xmm8 = _mm_load_pd(&p_N0[idx20]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx20],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx22]);
            xmm11 = _mm_load_pd(&p_Ts[idx22]);
            xmm12 = _mm_load_pd(&p_M[idx22]);
            xmm13 = _mm_load_pd(&p_N0[idx22]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx22],xmm14);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[idx24],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[idx24],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[idx24], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[idx24],_MM_HINT_T0);
#endif
            
            xmm0 = _mm_load_pd(&p_Ac[idx24]);
            xmm1 = _mm_load_pd(&p_Ts[idx24]);
            xmm2 = _mm_load_pd(&p_M[idx24]);
            xmm3 = _mm_load_pd(&p_N0[idx24]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx24],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx26]);
            xmm6 = _mm_load_pd(&p_Ts[idx26]);
            xmm7 = _mm_load_pd(&p_M[idx26]);
            xmm8 = _mm_load_pd(&p_N0[idx26]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx26],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx28]);
            xmm11 = _mm_load_pd(&p_Ts[idx28]);
            xmm12 = _mm_load_pd(&p_M[idx28]);
            xmm13 = _mm_load_pd(&p_N0[idx28]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx28],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx30]);
            xmm1 = _mm_load_pd(&p_Ts[idx30]);
            xmm2 = _mm_load_pd(&p_M[idx30]);
            xmm3 = _mm_load_pd(&p_N0[idx30]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx30],xmm4);
       }
 
       for(; (i+23ull) < a_sz; i += 24ull)
       {
            const std::size_t idx2  = i+2ull;
            const std::size_t idx4  = i+4ull;
            const std::size_t idx6  = i+6ull;
            const std::size_t idx8  = i+8ull;
            const std::size_t idx10 = i+10ull;
            const std::size_t idx12 = i+12ull;
            const std::size_t idx14 = i+14ull;
            const std::size_t idx16 = i+16ull;
            const std::size_t idx18 = i+18ull;
            const std::size_t idx20 = i+20ull;
            const std::size_t idx22 = i+22ull;
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
#endif 
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx2]);
            xmm6 = _mm_load_pd(&p_Ts[idx2]);
            xmm7 = _mm_load_pd(&p_M[idx2]);
            xmm8 = _mm_load_pd(&p_N0[idx2]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx2],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx4]);
            xmm11 = _mm_load_pd(&p_Ts[idx4]);
            xmm12 = _mm_load_pd(&p_M[idx4]);
            xmm13 = _mm_load_pd(&p_N0[idx4]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx4],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx6]);
            xmm1 = _mm_load_pd(&p_Ts[idx6]);
            xmm2 = _mm_load_pd(&p_M[idx6]);
            xmm3 = _mm_load_pd(&p_N0[idx6]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx6],xmm4);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[idx8],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[idx8],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[idx8], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[idx8],_MM_HINT_T0);
#endif
            xmm5 = _mm_load_pd(&p_Ac[idx8]);
            xmm6 = _mm_load_pd(&p_Ts[idx8]);
            xmm7 = _mm_load_pd(&p_M[idx8]);
            xmm8 = _mm_load_pd(&p_N0[idx8]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx8],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx10]);
            xmm11 = _mm_load_pd(&p_Ts[idx10]);
            xmm12 = _mm_load_pd(&p_M[idx10]);
            xmm13 = _mm_load_pd(&p_N0[idx10]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx10],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx12]);
            xmm1 = _mm_load_pd(&p_Ts[idx12]);
            xmm2 = _mm_load_pd(&p_M[idx12]);
            xmm3 = _mm_load_pd(&p_N0[idx12]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx12],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx14]);
            xmm6 = _mm_load_pd(&p_Ts[idx14]);
            xmm7 = _mm_load_pd(&p_M[idx14]);
            xmm8 = _mm_load_pd(&p_N0[idx14]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx14],xmm9);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[idx16],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[idx16],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[idx16], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[idx16],_MM_HINT_T0);
#endif

            xmm10 = _mm_load_pd(&p_Ac[idx16]);
            xmm11 = _mm_load_pd(&p_Ts[idx16]);
            xmm12 = _mm_load_pd(&p_M[idx16]);
            xmm13 = _mm_load_pd(&p_N0[idx16]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx16],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx18]);
            xmm1 = _mm_load_pd(&p_Ts[idx18]);
            xmm2 = _mm_load_pd(&p_M[idx18]);
            xmm3 = _mm_load_pd(&p_N0[idx18]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx18],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx20]);
            xmm6 = _mm_load_pd(&p_Ts[idx20]);
            xmm7 = _mm_load_pd(&p_M[idx20]);
            xmm8 = _mm_load_pd(&p_N0[idx20]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx20],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx22]);
            xmm11 = _mm_load_pd(&p_Ts[idx22]);
            xmm12 = _mm_load_pd(&p_M[idx22]);
            xmm13 = _mm_load_pd(&p_N0[idx22]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx22],xmm14);
       }

       for(; (i+15ull) < a_sz; i += 16ull) 
       {
            const std::size_t idx2  = i+2ull;
            const std::size_t idx4  = i+4ull;
            const std::size_t idx6  = i+6ull;
            const std::size_t idx8  = i+8ull;
            const std::size_t idx10 = i+10ull;
            const std::size_t idx12 = i+12ull;
            const std::size_t idx14 = i+14ull;
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
#endif 
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx2]);
            xmm6 = _mm_load_pd(&p_Ts[idx2]);
            xmm7 = _mm_load_pd(&p_M[idx2]);
            xmm8 = _mm_load_pd(&p_N0[idx2]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx2],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx4]);
            xmm11 = _mm_load_pd(&p_Ts[idx4]);
            xmm12 = _mm_load_pd(&p_M[idx4]);
            xmm13 = _mm_load_pd(&p_N0[idx4]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx4],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx6]);
            xmm1 = _mm_load_pd(&p_Ts[idx6]);
            xmm2 = _mm_load_pd(&p_M[idx6]);
            xmm3 = _mm_load_pd(&p_N0[idx6]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx6],xmm4);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[idx8],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[idx8],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[idx8], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[idx8],_MM_HINT_T0);
#endif
            xmm5 = _mm_load_pd(&p_Ac[idx8]);
            xmm6 = _mm_load_pd(&p_Ts[idx8]);
            xmm7 = _mm_load_pd(&p_M[idx8]);
            xmm8 = _mm_load_pd(&p_N0[idx8]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx8],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx10]);
            xmm11 = _mm_load_pd(&p_Ts[idx10]);
            xmm12 = _mm_load_pd(&p_M[idx10]);
            xmm13 = _mm_load_pd(&p_N0[idx10]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx10],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx12]);
            xmm1 = _mm_load_pd(&p_Ts[idx12]);
            xmm2 = _mm_load_pd(&p_M[idx12]);
            xmm3 = _mm_load_pd(&p_N0[idx12]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx12],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx14]);
            xmm6 = _mm_load_pd(&p_Ts[idx14]);
            xmm7 = _mm_load_pd(&p_M[idx14]);
            xmm8 = _mm_load_pd(&p_N0[idx14]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx14],xmm9);
       }

       for(; (i+9ull) < a_sz; i += 10ull)
       {
            const std::size_t idx2  = i+2ull;
            const std::size_t idx4  = i+4ull;
            const std::size_t idx6  = i+6ull;
            const std::size_t idx8  = i+8ull;
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
#endif 
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx2]);
            xmm6 = _mm_load_pd(&p_Ts[idx2]);
            xmm7 = _mm_load_pd(&p_M[idx2]);
            xmm8 = _mm_load_pd(&p_N0[idx2]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx2],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx4]);
            xmm11 = _mm_load_pd(&p_Ts[idx4]);
            xmm12 = _mm_load_pd(&p_M[idx4]);
            xmm13 = _mm_load_pd(&p_N0[idx4]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx4],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx6]);
            xmm1 = _mm_load_pd(&p_Ts[idx6]);
            xmm2 = _mm_load_pd(&p_M[idx6]);
            xmm3 = _mm_load_pd(&p_N0[idx6]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx6],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx8]);
            xmm6 = _mm_load_pd(&p_Ts[idx8]);
            xmm7 = _mm_load_pd(&p_M[idx8]);
            xmm8 = _mm_load_pd(&p_N0[idx8]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx8],xmm9);
       }

       for(; (i+3ull) < a_sz; i += 4ull)
       {
            const std::size_t idx2  = i+2ull;
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
#endif 
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx2]);
            xmm6 = _mm_load_pd(&p_Ts[idx2]);
            xmm7 = _mm_load_pd(&p_M[idx2]);
            xmm8 = _mm_load_pd(&p_N0[idx2]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx2],xmm9);
       }

       for(; (i+1ull) < a_sz; i += 2ull)
       {
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);
       }

       for(; (i+0ull) < a_sz; i += 1ull)
       {
            const double valAc = p_Ac[i];
            const double valTs = p_Ts[i];
            const double valM  = p_M[i];
            const double valN0 = p_N0[i];
            p_sep_mam[i]      = scalar_kernel(valAc,valTs,valM,valN0,n);
       }
       return (0);
    }
    else 
    {
       for(i = 0ull; (i+31ull) < a_sz; i += 32ull) 
       {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
             _mm_prefetch((const char*)&p_Ac[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+0ull], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+0ull],_MM_HINT_T1);
#endif 
            const std::size_t idx2  = i+2ull;
            const std::size_t idx4  = i+4ull;
            const std::size_t idx6  = i+6ull;
            const std::size_t idx8  = i+8ull;
            const std::size_t idx10 = i+10ull;
            const std::size_t idx12 = i+12ull;
            const std::size_t idx14 = i+14ull;
            const std::size_t idx16 = i+16ull;
            const std::size_t idx18 = i+18ull;
            const std::size_t idx20 = i+20ull;
            const std::size_t idx22 = i+22ull;
            const std::size_t idx24 = i+24ull;
            const std::size_t idx26 = i+26ull;
            const std::size_t idx28 = i+28ull;
            const std::size_t idx30 = i+30ull;

            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx2]);
            xmm6 = _mm_load_pd(&p_Ts[idx2]);
            xmm7 = _mm_load_pd(&p_M[idx2]);
            xmm8 = _mm_load_pd(&p_N0[idx2]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx2],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx4]);
            xmm11 = _mm_load_pd(&p_Ts[idx4]);
            xmm12 = _mm_load_pd(&p_M[idx4]);
            xmm13 = _mm_load_pd(&p_N0[idx4]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx4],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx6]);
            xmm1 = _mm_load_pd(&p_Ts[idx6]);
            xmm2 = _mm_load_pd(&p_M[idx6]);
            xmm3 = _mm_load_pd(&p_N0[idx6]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx6],xmm4);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[idx8],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[idx8],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[idx8], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[idx8],_MM_HINT_T0);
             _mm_prefetch((const char*)&p_Ac[i+L2_offset+idx8],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+idx8],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+idx8], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+idx8],_MM_HINT_T1);
#endif
            xmm5 = _mm_load_pd(&p_Ac[idx8]);
            xmm6 = _mm_load_pd(&p_Ts[idx8]);
            xmm7 = _mm_load_pd(&p_M[idx8]);
            xmm8 = _mm_load_pd(&p_N0[idx8]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx8],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx10]);
            xmm11 = _mm_load_pd(&p_Ts[idx10]);
            xmm12 = _mm_load_pd(&p_M[idx10]);
            xmm13 = _mm_load_pd(&p_N0[idx10]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx10],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx12]);
            xmm1 = _mm_load_pd(&p_Ts[idx12]);
            xmm2 = _mm_load_pd(&p_M[idx12]);
            xmm3 = _mm_load_pd(&p_N0[idx12]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx12],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx14]);
            xmm6 = _mm_load_pd(&p_Ts[idx14]);
            xmm7 = _mm_load_pd(&p_M[idx14]);
            xmm8 = _mm_load_pd(&p_N0[idx14]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx14],xmm9);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[idx16],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[idx16],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[idx16], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[idx16],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+idx16],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+idx16],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+idx16], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+idx16],_MM_HINT_T1);
#endif

            xmm10 = _mm_load_pd(&p_Ac[idx16]);
            xmm11 = _mm_load_pd(&p_Ts[idx16]);
            xmm12 = _mm_load_pd(&p_M[idx16]);
            xmm13 = _mm_load_pd(&p_N0[idx16]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx16],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx18]);
            xmm1 = _mm_load_pd(&p_Ts[idx18]);
            xmm2 = _mm_load_pd(&p_M[idx18]);
            xmm3 = _mm_load_pd(&p_N0[idx18]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx18],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx20]);
            xmm6 = _mm_load_pd(&p_Ts[idx20]);
            xmm7 = _mm_load_pd(&p_M[idx20]);
            xmm8 = _mm_load_pd(&p_N0[idx20]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx20],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx22]);
            xmm11 = _mm_load_pd(&p_Ts[idx22]);
            xmm12 = _mm_load_pd(&p_M[idx22]);
            xmm13 = _mm_load_pd(&p_N0[idx22]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx22],xmm14);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[idx24],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[idx24],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[idx24], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[idx24],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+idx24],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+idx24],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+idx24], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+idx24],_MM_HINT_T1);
#endif
            
            xmm0 = _mm_load_pd(&p_Ac[idx24]);
            xmm1 = _mm_load_pd(&p_Ts[idx24]);
            xmm2 = _mm_load_pd(&p_M[idx24]);
            xmm3 = _mm_load_pd(&p_N0[idx24]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx24],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx26]);
            xmm6 = _mm_load_pd(&p_Ts[idx26]);
            xmm7 = _mm_load_pd(&p_M[idx26]);
            xmm8 = _mm_load_pd(&p_N0[idx26]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx26],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx28]);
            xmm11 = _mm_load_pd(&p_Ts[idx28]);
            xmm12 = _mm_load_pd(&p_M[idx28]);
            xmm13 = _mm_load_pd(&p_N0[idx28]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx28],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx30]);
            xmm1 = _mm_load_pd(&p_Ts[idx30]);
            xmm2 = _mm_load_pd(&p_M[idx30]);
            xmm3 = _mm_load_pd(&p_N0[idx30]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx30],xmm4);
       }
 
       for(; (i+23ull) < a_sz; i += 24ull)
       {
            const std::size_t idx2  = i+2ull;
            const std::size_t idx4  = i+4ull;
            const std::size_t idx6  = i+6ull;
            const std::size_t idx8  = i+8ull;
            const std::size_t idx10 = i+10ull;
            const std::size_t idx12 = i+12ull;
            const std::size_t idx14 = i+14ull;
            const std::size_t idx16 = i+16ull;
            const std::size_t idx18 = i+18ull;
            const std::size_t idx20 = i+20ull;
            const std::size_t idx22 = i+22ull;
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+0ull], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+0ull],_MM_HINT_T1);
#endif 
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx2]);
            xmm6 = _mm_load_pd(&p_Ts[idx2]);
            xmm7 = _mm_load_pd(&p_M[idx2]);
            xmm8 = _mm_load_pd(&p_N0[idx2]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx2],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx4]);
            xmm11 = _mm_load_pd(&p_Ts[idx4]);
            xmm12 = _mm_load_pd(&p_M[idx4]);
            xmm13 = _mm_load_pd(&p_N0[idx4]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx4],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx6]);
            xmm1 = _mm_load_pd(&p_Ts[idx6]);
            xmm2 = _mm_load_pd(&p_M[idx6]);
            xmm3 = _mm_load_pd(&p_N0[idx6]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx6],xmm4);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[idx8],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[idx8],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[idx8], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[idx8],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+idx8],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+idx8],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+idx8], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+idx8],_MM_HINT_T1);
#endif
            xmm5 = _mm_load_pd(&p_Ac[idx8]);
            xmm6 = _mm_load_pd(&p_Ts[idx8]);
            xmm7 = _mm_load_pd(&p_M[idx8]);
            xmm8 = _mm_load_pd(&p_N0[idx8]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx8],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx10]);
            xmm11 = _mm_load_pd(&p_Ts[idx10]);
            xmm12 = _mm_load_pd(&p_M[idx10]);
            xmm13 = _mm_load_pd(&p_N0[idx10]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx10],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx12]);
            xmm1 = _mm_load_pd(&p_Ts[idx12]);
            xmm2 = _mm_load_pd(&p_M[idx12]);
            xmm3 = _mm_load_pd(&p_N0[idx12]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx12],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx14]);
            xmm6 = _mm_load_pd(&p_Ts[idx14]);
            xmm7 = _mm_load_pd(&p_M[idx14]);
            xmm8 = _mm_load_pd(&p_N0[idx14]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx14],xmm9);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[idx16],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[idx16],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[idx16], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[idx16],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+idx16],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+idx16],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+idx16], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+idx16],_MM_HINT_T1);
#endif

            xmm10 = _mm_load_pd(&p_Ac[idx16]);
            xmm11 = _mm_load_pd(&p_Ts[idx16]);
            xmm12 = _mm_load_pd(&p_M[idx16]);
            xmm13 = _mm_load_pd(&p_N0[idx16]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx16],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx18]);
            xmm1 = _mm_load_pd(&p_Ts[idx18]);
            xmm2 = _mm_load_pd(&p_M[idx18]);
            xmm3 = _mm_load_pd(&p_N0[idx18]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx18],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx20]);
            xmm6 = _mm_load_pd(&p_Ts[idx20]);
            xmm7 = _mm_load_pd(&p_M[idx20]);
            xmm8 = _mm_load_pd(&p_N0[idx20]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx20],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx22]);
            xmm11 = _mm_load_pd(&p_Ts[idx22]);
            xmm12 = _mm_load_pd(&p_M[idx22]);
            xmm13 = _mm_load_pd(&p_N0[idx22]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx22],xmm14);
       }

       for(; (i+15ull) < a_sz; i += 16ull) 
       {
            const std::size_t idx2  = i+2ull;
            const std::size_t idx4  = i+4ull;
            const std::size_t idx6  = i+6ull;
            const std::size_t idx8  = i+8ull;
            const std::size_t idx10 = i+10ull;
            const std::size_t idx12 = i+12ull;
            const std::size_t idx14 = i+14ull;
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
             _mm_prefetch((const char*)&p_Ac[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+0ull], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+0ull],_MM_HINT_T1);
#endif 
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx2]);
            xmm6 = _mm_load_pd(&p_Ts[idx2]);
            xmm7 = _mm_load_pd(&p_M[idx2]);
            xmm8 = _mm_load_pd(&p_N0[idx2]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx2],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx4]);
            xmm11 = _mm_load_pd(&p_Ts[idx4]);
            xmm12 = _mm_load_pd(&p_M[idx4]);
            xmm13 = _mm_load_pd(&p_N0[idx4]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx4],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx6]);
            xmm1 = _mm_load_pd(&p_Ts[idx6]);
            xmm2 = _mm_load_pd(&p_M[idx6]);
            xmm3 = _mm_load_pd(&p_N0[idx6]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx6],xmm4);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[idx8],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[idx8],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[idx8], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[idx8],_MM_HINT_T0);
             _mm_prefetch((const char*)&p_Ac[i+L2_offset+idx8],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+idx8],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+idx8], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+idx8],_MM_HINT_T1);
#endif
            xmm5 = _mm_load_pd(&p_Ac[idx8]);
            xmm6 = _mm_load_pd(&p_Ts[idx8]);
            xmm7 = _mm_load_pd(&p_M[idx8]);
            xmm8 = _mm_load_pd(&p_N0[idx8]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx8],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx10]);
            xmm11 = _mm_load_pd(&p_Ts[idx10]);
            xmm12 = _mm_load_pd(&p_M[idx10]);
            xmm13 = _mm_load_pd(&p_N0[idx10]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx10],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx12]);
            xmm1 = _mm_load_pd(&p_Ts[idx12]);
            xmm2 = _mm_load_pd(&p_M[idx12]);
            xmm3 = _mm_load_pd(&p_N0[idx12]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx12],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx14]);
            xmm6 = _mm_load_pd(&p_Ts[idx14]);
            xmm7 = _mm_load_pd(&p_M[idx14]);
            xmm8 = _mm_load_pd(&p_N0[idx14]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx14],xmm9);
       }

       for(; (i+9ull) < a_sz; i += 10ull)
       {
            const std::size_t idx2  = i+2ull;
            const std::size_t idx4  = i+4ull;
            const std::size_t idx6  = i+6ull;
            const std::size_t idx8  = i+8ull;
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
             _mm_prefetch((const char*)&p_Ac[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+0ull], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+0ull],_MM_HINT_T1);
#endif 
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx2]);
            xmm6 = _mm_load_pd(&p_Ts[idx2]);
            xmm7 = _mm_load_pd(&p_M[idx2]);
            xmm8 = _mm_load_pd(&p_N0[idx2]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx2],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx4]);
            xmm11 = _mm_load_pd(&p_Ts[idx4]);
            xmm12 = _mm_load_pd(&p_M[idx4]);
            xmm13 = _mm_load_pd(&p_N0[idx4]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx4],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx6]);
            xmm1 = _mm_load_pd(&p_Ts[idx6]);
            xmm2 = _mm_load_pd(&p_M[idx6]);
            xmm3 = _mm_load_pd(&p_N0[idx6]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx6],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx8]);
            xmm6 = _mm_load_pd(&p_Ts[idx8]);
            xmm7 = _mm_load_pd(&p_M[idx8]);
            xmm8 = _mm_load_pd(&p_N0[idx8]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx8],xmm9);
       }

       for(; (i+3ull) < a_sz; i += 4ull)
       {
            const std::size_t idx2  = i+2ull;
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
             _mm_prefetch((const char*)&p_Ac[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+0ull], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+0ull],_MM_HINT_T1);
#endif 
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx2]);
            xmm6 = _mm_load_pd(&p_Ts[idx2]);
            xmm7 = _mm_load_pd(&p_M[idx2]);
            xmm8 = _mm_load_pd(&p_N0[idx2]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx2],xmm9);
       }

       for(; (i+1ull) < a_sz; i += 2ull)
       {
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);
       }

       for(; (i+0ull) < a_sz; i += 1ull)
       {
            const double valAc = p_Ac[i];
            const double valTs = p_Ts[i];
            const double valM  = p_M[i];
            const double valN0 = p_N0[i];
            p_sep_mam[i]      = scalar_kernel(valAc,valTs,valM,valN0,n);
       }
       return (0);
    }

    }
}

//////////////////////////////////////////////////////////////////////////////////

#if (ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F64_OVERRIDE_COMPILER_CMD_LINE) == 1
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
analytic_SEP_BEP_generic_2xf64_array1d_u12x(const double * __restrict__ Ac, 
                                        const double * __restrict__ Ts,
                                        const double * __restrict__ M,
                                        const double * __restrict__ N0,
                                        double       * __restrict__ sep_mam,
                                        const std::int32_t n,
                                        std::size_t sz,
                                        const callable_sse_kernel     & sse_kernel,
                                        const callable_scalar_kernel  & scalar_kernel) 
{
    if(__builtin_expect(sz<1ull,0)) {return (-1);}
    __m128d vAc;
    __m128d vTs;
    __m128d vM;
    __m128d vN0;
    const double * __restrict__ p_Ac      = Ac;
    const double * __restrict__ p_Ts      = Ts;
    const double * __restrict__ p_M       = M;
    const double * __restrict__ p_N0      = N0;
    double       * __restrict__ p_sep_mam = sep_mam;
    std::size_t             a_sz         = sz;
    if(a_sz>=1ull && a_sz<2ull)
    {
       for(std::uint32_t i = 0; i<2; ++i) 
       {
           const double valAc = p_Ac[i];
           const double valTs = p_Ts[i];
           const double valM  = p_M[i];
           const double valN0 = p_N0[i];
           p_sep_mam[i]      = scalar_kernel(valAc,valTs,valM,valN0,n);
       }
       return (0);
    }
    else if(a_sz==2ull)
    {
        vAc = _mm_loadu_pd(&p_Ac[0]);
        vTs = _mm_loadu_pd(&p_Ts[0]);
        vM  = _mm_loadu_pd(&p_M[0]);
        vN0 = _mm_loadu_pd(&p_N0[0]);
        _mm_storeu_pd(&p_sep_mam[0], sse_kernel(vAc,vTs,vM,vN0,n));
        return (0);
    }
    else if(a_sz>2ull && a_sz<=24ull)
    {
        std::size_t i,j;
        for(i = 0ull; i!=ROUND_DOWN(a_sz,1ull); i += 2ull)
        {
            vAc = _mm_loadu_pd(&p_Ac[i]);
            vTs = _mm_loadu_pd(&p_Ts[i]);
            vM  = _mm_loadu_pd(&p_M[i]);
            vN0 = _mm_loadu_pd(&p_N0[i]);
            _mm_storeu_pd(&p_sep_mam[i], sse_kernel(vAc,vTs,vM,vN0,n));
        }
        for(j = i; j<a_sz; ++j)  
        {
            const double valAc = p_Ac[j];
            const double valTs = p_Ts[j];
            const double valM  = p_M[j];
            const double valN0 = p_N0[j];
            p_sep_mam[j]      = scalar_kernel(valAc,valTs,valM,valN0,n);
        }
        return (0);
    }
    else if(a_sz>24ull)
    {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_USE_PEEL_LOOP) == 1
       while(((std::uintptr_t)&p_sep_mam & 15ull) && a_sz)
       {
           const double valAc = *p_Ac;
           const double valTs = *p_Ts;
           const double valM  = *p_M;
           const double valN0 = *p_N0;
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
       __m128d xmm0;
       __m128d xmm1;
       __m128d xmm2;
       __m128d xmm3;
       __m128d xmm4;
       __m128d xmm5;
       __m128d xmm6;
       __m128d xmm7;
       __m128d xmm8;
       __m128d xmm9;
       __m128d xmm10;
       __m128d xmm11;
       __m128d xmm12;
       __m128d xmm13;
       __m128d xmm14;
       std::size_t i,j;
       constexpr std::size_t L2_offset = 4000ull;
       constexpr std::size_t L1D_boundary = L2_offset;
    if(a_sz <= L1D_boundary)
    {
       for(i = 0ull; (i+23ull) < a_sz; i += 24ull) 
       {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
#endif 
            const std::size_t idx2  = i+2ull;
            const std::size_t idx4  = i+4ull;
            const std::size_t idx6  = i+6ull;
            const std::size_t idx8  = i+8ull;
            const std::size_t idx10 = i+10ull;
            const std::size_t idx12 = i+12ull;
            const std::size_t idx14 = i+14ull;
            const std::size_t idx16 = i+16ull;
            const std::size_t idx18 = i+18ull;
            const std::size_t idx20 = i+20ull;
            const std::size_t idx22 = i+22ull;

            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx2]);
            xmm6 = _mm_load_pd(&p_Ts[idx2]);
            xmm7 = _mm_load_pd(&p_M[idx2]);
            xmm8 = _mm_load_pd(&p_N0[idx2]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx2],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx4]);
            xmm11 = _mm_load_pd(&p_Ts[idx4]);
            xmm12 = _mm_load_pd(&p_M[idx4]);
            xmm13 = _mm_load_pd(&p_N0[idx4]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx4],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx6]);
            xmm1 = _mm_load_pd(&p_Ts[idx6]);
            xmm2 = _mm_load_pd(&p_M[idx6]);
            xmm3 = _mm_load_pd(&p_N0[idx6]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx6],xmm4);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[idx8],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[idx8],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[idx8], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[idx8],_MM_HINT_T0);
#endif
            xmm5 = _mm_load_pd(&p_Ac[idx8]);
            xmm6 = _mm_load_pd(&p_Ts[idx8]);
            xmm7 = _mm_load_pd(&p_M[idx8]);
            xmm8 = _mm_load_pd(&p_N0[idx8]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx8],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx10]);
            xmm11 = _mm_load_pd(&p_Ts[idx10]);
            xmm12 = _mm_load_pd(&p_M[idx10]);
            xmm13 = _mm_load_pd(&p_N0[idx10]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx10],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx12]);
            xmm1 = _mm_load_pd(&p_Ts[idx12]);
            xmm2 = _mm_load_pd(&p_M[idx12]);
            xmm3 = _mm_load_pd(&p_N0[idx12]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx12],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx14]);
            xmm6 = _mm_load_pd(&p_Ts[idx14]);
            xmm7 = _mm_load_pd(&p_M[idx14]);
            xmm8 = _mm_load_pd(&p_N0[idx14]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx14],xmm9);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[idx16],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[idx16],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[idx16], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[idx16],_MM_HINT_T0);
#endif

            xmm10 = _mm_load_pd(&p_Ac[idx16]);
            xmm11 = _mm_load_pd(&p_Ts[idx16]);
            xmm12 = _mm_load_pd(&p_M[idx16]);
            xmm13 = _mm_load_pd(&p_N0[idx16]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx16],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx18]);
            xmm1 = _mm_load_pd(&p_Ts[idx18]);
            xmm2 = _mm_load_pd(&p_M[idx18]);
            xmm3 = _mm_load_pd(&p_N0[idx18]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx18],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx20]);
            xmm6 = _mm_load_pd(&p_Ts[idx20]);
            xmm7 = _mm_load_pd(&p_M[idx20]);
            xmm8 = _mm_load_pd(&p_N0[idx20]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx20],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx22]);
            xmm11 = _mm_load_pd(&p_Ts[idx22]);
            xmm12 = _mm_load_pd(&p_M[idx22]);
            xmm13 = _mm_load_pd(&p_N0[idx22]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx22],xmm14);

       }
 
       for(; (i+15ull) < a_sz; i += 16ull) 
       {
            const std::size_t idx2  = i+2ull;
            const std::size_t idx4  = i+4ull;
            const std::size_t idx6  = i+6ull;
            const std::size_t idx8  = i+8ull;
            const std::size_t idx10 = i+10ull;
            const std::size_t idx12 = i+12ull;
            const std::size_t idx14 = i+14ull;
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
#endif 
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx2]);
            xmm6 = _mm_load_pd(&p_Ts[idx2]);
            xmm7 = _mm_load_pd(&p_M[idx2]);
            xmm8 = _mm_load_pd(&p_N0[idx2]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx2],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx4]);
            xmm11 = _mm_load_pd(&p_Ts[idx4]);
            xmm12 = _mm_load_pd(&p_M[idx4]);
            xmm13 = _mm_load_pd(&p_N0[idx4]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx4],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx6]);
            xmm1 = _mm_load_pd(&p_Ts[idx6]);
            xmm2 = _mm_load_pd(&p_M[idx6]);
            xmm3 = _mm_load_pd(&p_N0[idx6]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx6],xmm4);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[idx8],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[idx8],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[idx8], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[idx8],_MM_HINT_T0);
#endif
            xmm5 = _mm_load_pd(&p_Ac[idx8]);
            xmm6 = _mm_load_pd(&p_Ts[idx8]);
            xmm7 = _mm_load_pd(&p_M[idx8]);
            xmm8 = _mm_load_pd(&p_N0[idx8]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx8],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx10]);
            xmm11 = _mm_load_pd(&p_Ts[idx10]);
            xmm12 = _mm_load_pd(&p_M[idx10]);
            xmm13 = _mm_load_pd(&p_N0[idx10]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx10],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx12]);
            xmm1 = _mm_load_pd(&p_Ts[idx12]);
            xmm2 = _mm_load_pd(&p_M[idx12]);
            xmm3 = _mm_load_pd(&p_N0[idx12]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx12],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx14]);
            xmm6 = _mm_load_pd(&p_Ts[idx14]);
            xmm7 = _mm_load_pd(&p_M[idx14]);
            xmm8 = _mm_load_pd(&p_N0[idx14]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx14],xmm9);
       }

       for(; (i+9ull) < a_sz; i += 10ull)
       {
            const std::size_t idx2  = i+2ull;
            const std::size_t idx4  = i+4ull;
            const std::size_t idx6  = i+6ull;
            const std::size_t idx8  = i+8ull;
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
#endif 
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx2]);
            xmm6 = _mm_load_pd(&p_Ts[idx2]);
            xmm7 = _mm_load_pd(&p_M[idx2]);
            xmm8 = _mm_load_pd(&p_N0[idx2]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx2],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx4]);
            xmm11 = _mm_load_pd(&p_Ts[idx4]);
            xmm12 = _mm_load_pd(&p_M[idx4]);
            xmm13 = _mm_load_pd(&p_N0[idx4]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx4],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx6]);
            xmm1 = _mm_load_pd(&p_Ts[idx6]);
            xmm2 = _mm_load_pd(&p_M[idx6]);
            xmm3 = _mm_load_pd(&p_N0[idx6]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx6],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx8]);
            xmm6 = _mm_load_pd(&p_Ts[idx8]);
            xmm7 = _mm_load_pd(&p_M[idx8]);
            xmm8 = _mm_load_pd(&p_N0[idx8]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx8],xmm9);
       }

       for(; (i+3ull) < a_sz; i += 4ull)
       {
            const std::size_t idx2  = i+2ull;
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx2]);
            xmm6 = _mm_load_pd(&p_Ts[idx2]);
            xmm7 = _mm_load_pd(&p_M[idx2]);
            xmm8 = _mm_load_pd(&p_N0[idx2]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx2],xmm9);
       }

       for(; (i+1ull) < a_sz; i += 2ull)
       {
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);
       }

       for(; (i+0ull) < a_sz; i += 1ull)
       {
            const double valAc = p_Ac[i];
            const double valTs = p_Ts[i];
            const double valM  = p_M[i];
            const double valN0 = p_N0[i];
            p_sep_mam[i]      = scalar_kernel(valAc,valTs,valM,valN0,n);
       }
       return (0);
    }
    else 
    {
       for(i = 0ull; (i+23ull) < a_sz; i += 24ull) 
       {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
             _mm_prefetch((const char*)&p_Ac[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+0ull], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+0ull],_MM_HINT_T1);
#endif 
            const std::size_t idx2  = i+2ull;
            const std::size_t idx4  = i+4ull;
            const std::size_t idx6  = i+6ull;
            const std::size_t idx8  = i+8ull;
            const std::size_t idx10 = i+10ull;
            const std::size_t idx12 = i+12ull;
            const std::size_t idx14 = i+14ull;
            const std::size_t idx16 = i+16ull;
            const std::size_t idx18 = i+18ull;
            const std::size_t idx20 = i+20ull;
            const std::size_t idx22 = i+22ull;

            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx2]);
            xmm6 = _mm_load_pd(&p_Ts[idx2]);
            xmm7 = _mm_load_pd(&p_M[idx2]);
            xmm8 = _mm_load_pd(&p_N0[idx2]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx2],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx4]);
            xmm11 = _mm_load_pd(&p_Ts[idx4]);
            xmm12 = _mm_load_pd(&p_M[idx4]);
            xmm13 = _mm_load_pd(&p_N0[idx4]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx4],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx6]);
            xmm1 = _mm_load_pd(&p_Ts[idx6]);
            xmm2 = _mm_load_pd(&p_M[idx6]);
            xmm3 = _mm_load_pd(&p_N0[idx6]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx6],xmm4);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[idx8],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[idx8],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[idx8], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[idx8],_MM_HINT_T0);
             _mm_prefetch((const char*)&p_Ac[i+L2_offset+idx8],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+idx8],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+idx8], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+idx8],_MM_HINT_T1);
#endif
            xmm5 = _mm_load_pd(&p_Ac[idx8]);
            xmm6 = _mm_load_pd(&p_Ts[idx8]);
            xmm7 = _mm_load_pd(&p_M[idx8]);
            xmm8 = _mm_load_pd(&p_N0[idx8]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx8],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx10]);
            xmm11 = _mm_load_pd(&p_Ts[idx10]);
            xmm12 = _mm_load_pd(&p_M[idx10]);
            xmm13 = _mm_load_pd(&p_N0[idx10]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx10],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx12]);
            xmm1 = _mm_load_pd(&p_Ts[idx12]);
            xmm2 = _mm_load_pd(&p_M[idx12]);
            xmm3 = _mm_load_pd(&p_N0[idx12]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx12],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx14]);
            xmm6 = _mm_load_pd(&p_Ts[idx14]);
            xmm7 = _mm_load_pd(&p_M[idx14]);
            xmm8 = _mm_load_pd(&p_N0[idx14]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx14],xmm9);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[idx16],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[idx16],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[idx16], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[idx16],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+idx16],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+idx16],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+idx16], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+idx16],_MM_HINT_T1);
#endif

            xmm10 = _mm_load_pd(&p_Ac[idx16]);
            xmm11 = _mm_load_pd(&p_Ts[idx16]);
            xmm12 = _mm_load_pd(&p_M[idx16]);
            xmm13 = _mm_load_pd(&p_N0[idx16]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx16],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx18]);
            xmm1 = _mm_load_pd(&p_Ts[idx18]);
            xmm2 = _mm_load_pd(&p_M[idx18]);
            xmm3 = _mm_load_pd(&p_N0[idx18]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx18],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx20]);
            xmm6 = _mm_load_pd(&p_Ts[idx20]);
            xmm7 = _mm_load_pd(&p_M[idx20]);
            xmm8 = _mm_load_pd(&p_N0[idx20]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx20],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx22]);
            xmm11 = _mm_load_pd(&p_Ts[idx22]);
            xmm12 = _mm_load_pd(&p_M[idx22]);
            xmm13 = _mm_load_pd(&p_N0[idx22]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx22],xmm14);

       }
 
       for(; (i+15ull) < a_sz; i += 16ull) 
       {
            const std::size_t idx2  = i+2ull;
            const std::size_t idx4  = i+4ull;
            const std::size_t idx6  = i+6ull;
            const std::size_t idx8  = i+8ull;
            const std::size_t idx10 = i+10ull;
            const std::size_t idx12 = i+12ull;
            const std::size_t idx14 = i+14ull;
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
             _mm_prefetch((const char*)&p_Ac[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+0ull], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+0ull],_MM_HINT_T1);
#endif 
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx2]);
            xmm6 = _mm_load_pd(&p_Ts[idx2]);
            xmm7 = _mm_load_pd(&p_M[idx2]);
            xmm8 = _mm_load_pd(&p_N0[idx2]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx2],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx4]);
            xmm11 = _mm_load_pd(&p_Ts[idx4]);
            xmm12 = _mm_load_pd(&p_M[idx4]);
            xmm13 = _mm_load_pd(&p_N0[idx4]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx4],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx6]);
            xmm1 = _mm_load_pd(&p_Ts[idx6]);
            xmm2 = _mm_load_pd(&p_M[idx6]);
            xmm3 = _mm_load_pd(&p_N0[idx6]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx6],xmm4);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[idx8],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[idx8],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[idx8], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[idx8],_MM_HINT_T0);
             _mm_prefetch((const char*)&p_Ac[i+L2_offset+idx8],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+idx8],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+idx8], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+idx8],_MM_HINT_T1);
#endif
            xmm5 = _mm_load_pd(&p_Ac[idx8]);
            xmm6 = _mm_load_pd(&p_Ts[idx8]);
            xmm7 = _mm_load_pd(&p_M[idx8]);
            xmm8 = _mm_load_pd(&p_N0[idx8]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx8],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx10]);
            xmm11 = _mm_load_pd(&p_Ts[idx10]);
            xmm12 = _mm_load_pd(&p_M[idx10]);
            xmm13 = _mm_load_pd(&p_N0[idx10]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx10],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx12]);
            xmm1 = _mm_load_pd(&p_Ts[idx12]);
            xmm2 = _mm_load_pd(&p_M[idx12]);
            xmm3 = _mm_load_pd(&p_N0[idx12]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx12],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx14]);
            xmm6 = _mm_load_pd(&p_Ts[idx14]);
            xmm7 = _mm_load_pd(&p_M[idx14]);
            xmm8 = _mm_load_pd(&p_N0[idx14]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx14],xmm9);
       }

       for(; (i+9ull) < a_sz; i += 10ull)
       {
            const std::size_t idx2  = i+2ull;
            const std::size_t idx4  = i+4ull;
            const std::size_t idx6  = i+6ull;
            const std::size_t idx8  = i+8ull;
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
             _mm_prefetch((const char*)&p_Ac[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+0ull], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+0ull],_MM_HINT_T1);
#endif 
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx2]);
            xmm6 = _mm_load_pd(&p_Ts[idx2]);
            xmm7 = _mm_load_pd(&p_M[idx2]);
            xmm8 = _mm_load_pd(&p_N0[idx2]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx2],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx4]);
            xmm11 = _mm_load_pd(&p_Ts[idx4]);
            xmm12 = _mm_load_pd(&p_M[idx4]);
            xmm13 = _mm_load_pd(&p_N0[idx4]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx4],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx6]);
            xmm1 = _mm_load_pd(&p_Ts[idx6]);
            xmm2 = _mm_load_pd(&p_M[idx6]);
            xmm3 = _mm_load_pd(&p_N0[idx6]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx6],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx8]);
            xmm6 = _mm_load_pd(&p_Ts[idx8]);
            xmm7 = _mm_load_pd(&p_M[idx8]);
            xmm8 = _mm_load_pd(&p_N0[idx8]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx8],xmm9);
       }

       for(; (i+3ull) < a_sz; i += 4ull)
       {
            const std::size_t idx2  = i+2ull;
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx2]);
            xmm6 = _mm_load_pd(&p_Ts[idx2]);
            xmm7 = _mm_load_pd(&p_M[idx2]);
            xmm8 = _mm_load_pd(&p_N0[idx2]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx2],xmm9);
       }

       for(; (i+1ull) < a_sz; i += 2ull)
       {
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);
       }

       for(; (i+0ull) < a_sz; i += 1ull)
       {
            const double valAc = p_Ac[i];
            const double valTs = p_Ts[i];
            const double valM  = p_M[i];
            const double valN0 = p_N0[i];
            p_sep_mam[i]      = scalar_kernel(valAc,valTs,valM,valN0,n);
       }
       return (0);
    }
    }
}

///////////////////////////////////////////////////////////////////////////////////

#if (ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F64_OVERRIDE_COMPILER_CMD_LINE) == 1
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
analytic_SEP_BEP_generic_2xf64_array1d_u8x(const double * __restrict__ Ac, 
                                        const double * __restrict__ Ts,
                                        const double * __restrict__ M,
                                        const double * __restrict__ N0,
                                        double       * __restrict__ sep_mam,
                                        const std::int32_t n,
                                        std::size_t sz,
                                        const callable_sse_kernel     & sse_kernel,
                                        const callable_scalar_kernel  & scalar_kernel) 
{
    if(__builtin_expect(sz<1ull,0)) {return (-1);}
    __m128d vAc;
    __m128d vTs;
    __m128d vM;
    __m128d vN0;
    const double * __restrict__ p_Ac      = Ac;
    const double * __restrict__ p_Ts      = Ts;
    const double * __restrict__ p_M       = M;
    const double * __restrict__ p_N0      = N0;
    double       * __restrict__ p_sep_mam = sep_mam;
    std::size_t             a_sz         = sz;
    if(a_sz>=1ull && a_sz<2ull)
    {
       for(std::uint32_t i = 0; i<2; ++i) 
       {
           const double valAc = p_Ac[i];
           const double valTs = p_Ts[i];
           const double valM  = p_M[i];
           const double valN0 = p_N0[i];
           p_sep_mam[i]      = scalar_kernel(valAc,valTs,valM,valN0,n);
       }
       return (0);
    }
    else if(a_sz==2ull)
    {
        vAc = _mm_loadu_pd(&p_Ac[0]);
        vTs = _mm_loadu_pd(&p_Ts[0]);
        vM  = _mm_loadu_pd(&p_M[0]);
        vN0 = _mm_loadu_pd(&p_N0[0]);
        _mm_storeu_pd(&p_sep_mam[0], sse_kernel(vAc,vTs,vM,vN0,n));
        return (0);
    }
    else if(a_sz>2ull && a_sz<=16ull)
    {
        std::size_t i,j;
        for(i = 0ull; i!=ROUND_DOWN(a_sz,1ull); i += 2ull)
        {
            vAc = _mm_loadu_pd(&p_Ac[i]);
            vTs = _mm_loadu_pd(&p_Ts[i]);
            vM  = _mm_loadu_pd(&p_M[i]);
            vN0 = _mm_loadu_pd(&p_N0[i]);
            _mm_storeu_pd(&p_sep_mam[i], sse_kernel(vAc,vTs,vM,vN0,n));
        }
        for(j = i; j<a_sz; ++j)  
        {
            const double valAc = p_Ac[j];
            const double valTs = p_Ts[j];
            const double valM  = p_M[j];
            const double valN0 = p_N0[j];
            p_sep_mam[j]      = scalar_kernel(valAc,valTs,valM,valN0,n);
        }
        return (0);
    }
    else if(a_sz>16ull)
    {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_USE_PEEL_LOOP) == 1
       while(((std::uintptr_t)&p_sep_mam & 15ull) && a_sz)
       {
           const double valAc = *p_Ac;
           const double valTs = *p_Ts;
           const double valM  = *p_M;
           const double valN0 = *p_N0;
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
       __m128d xmm0;
       __m128d xmm1;
       __m128d xmm2;
       __m128d xmm3;
       __m128d xmm4;
       __m128d xmm5;
       __m128d xmm6;
       __m128d xmm7;
       __m128d xmm8;
       __m128d xmm9;
       __m128d xmm10;
       __m128d xmm11;
       __m128d xmm12;
       __m128d xmm13;
       __m128d xmm14;
       std::size_t i,j;
       constexpr std::size_t L2_offset = 8000ull;
       constexpr std::size_t L1D_boundary = L2_offset;
    if(a_sz <= L1D_boundary)
    {
        for(i = 0ull; (i+15ull) < a_sz; i += 16ull) 
       {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
#endif 
            const std::size_t idx2  = i+2ull;
            const std::size_t idx4  = i+4ull;
            const std::size_t idx6  = i+6ull;
            const std::size_t idx8  = i+8ull;
            const std::size_t idx10 = i+10ull;
            const std::size_t idx12 = i+12ull;
            const std::size_t idx14 = i+14ull;

            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx2]);
            xmm6 = _mm_load_pd(&p_Ts[idx2]);
            xmm7 = _mm_load_pd(&p_M[idx2]);
            xmm8 = _mm_load_pd(&p_N0[idx2]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx2],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx4]);
            xmm11 = _mm_load_pd(&p_Ts[idx4]);
            xmm12 = _mm_load_pd(&p_M[idx4]);
            xmm13 = _mm_load_pd(&p_N0[idx4]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx4],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx6]);
            xmm1 = _mm_load_pd(&p_Ts[idx6]);
            xmm2 = _mm_load_pd(&p_M[idx6]);
            xmm3 = _mm_load_pd(&p_N0[idx6]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx6],xmm4);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[idx8],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[idx8],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[idx8], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[idx8],_MM_HINT_T0);
#endif
            xmm5 = _mm_load_pd(&p_Ac[idx8]);
            xmm6 = _mm_load_pd(&p_Ts[idx8]);
            xmm7 = _mm_load_pd(&p_M[idx8]);
            xmm8 = _mm_load_pd(&p_N0[idx8]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx8],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx10]);
            xmm11 = _mm_load_pd(&p_Ts[idx10]);
            xmm12 = _mm_load_pd(&p_M[idx10]);
            xmm13 = _mm_load_pd(&p_N0[idx10]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx10],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx12]);
            xmm1 = _mm_load_pd(&p_Ts[idx12]);
            xmm2 = _mm_load_pd(&p_M[idx12]);
            xmm3 = _mm_load_pd(&p_N0[idx12]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx12],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx14]);
            xmm6 = _mm_load_pd(&p_Ts[idx14]);
            xmm7 = _mm_load_pd(&p_M[idx14]);
            xmm8 = _mm_load_pd(&p_N0[idx14]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx14],xmm9);
       }
 
       for(; (i+9ull) < a_sz; i += 10ull)
       {
            const std::size_t idx2  = i+2ull;
            const std::size_t idx4  = i+4ull;
            const std::size_t idx6  = i+6ull;
            const std::size_t idx8  = i+8ull;
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
#endif 
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx2]);
            xmm6 = _mm_load_pd(&p_Ts[idx2]);
            xmm7 = _mm_load_pd(&p_M[idx2]);
            xmm8 = _mm_load_pd(&p_N0[idx2]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx2],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx4]);
            xmm11 = _mm_load_pd(&p_Ts[idx4]);
            xmm12 = _mm_load_pd(&p_M[idx4]);
            xmm13 = _mm_load_pd(&p_N0[idx4]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx4],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx6]);
            xmm1 = _mm_load_pd(&p_Ts[idx6]);
            xmm2 = _mm_load_pd(&p_M[idx6]);
            xmm3 = _mm_load_pd(&p_N0[idx6]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx6],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx8]);
            xmm6 = _mm_load_pd(&p_Ts[idx8]);
            xmm7 = _mm_load_pd(&p_M[idx8]);
            xmm8 = _mm_load_pd(&p_N0[idx8]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx8],xmm9);
       }

       for(; (i+3ull) < a_sz; i += 4ull)
       {
            const std::size_t idx2  = i+2ull;
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx2]);
            xmm6 = _mm_load_pd(&p_Ts[idx2]);
            xmm7 = _mm_load_pd(&p_M[idx2]);
            xmm8 = _mm_load_pd(&p_N0[idx2]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx2],xmm9);
       }

       for(; (i+1ull) < a_sz; i += 2ull)
       {
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);
       }

       for(; (i+0ull) < a_sz; i += 1ull)
       {
            const double valAc = p_Ac[i];
            const double valTs = p_Ts[i];
            const double valM  = p_M[i];
            const double valN0 = p_N0[i];
            p_sep_mam[i]      = scalar_kernel(valAc,valTs,valM,valN0,n);
       }
       return (0);
    }
    else 
    {
       for(i = 0ull; (i+15ull) < a_sz; i += 16ull) 
       {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
             _mm_prefetch((const char*)&p_Ac[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+0ull], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+0ull],_MM_HINT_T1);
#endif 
            const std::size_t idx2  = i+2ull;
            const std::size_t idx4  = i+4ull;
            const std::size_t idx6  = i+6ull;
            const std::size_t idx8  = i+8ull;
            const std::size_t idx10 = i+10ull;
            const std::size_t idx12 = i+12ull;
            const std::size_t idx14 = i+14ull;
        
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx2]);
            xmm6 = _mm_load_pd(&p_Ts[idx2]);
            xmm7 = _mm_load_pd(&p_M[idx2]);
            xmm8 = _mm_load_pd(&p_N0[idx2]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx2],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx4]);
            xmm11 = _mm_load_pd(&p_Ts[idx4]);
            xmm12 = _mm_load_pd(&p_M[idx4]);
            xmm13 = _mm_load_pd(&p_N0[idx4]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx4],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx6]);
            xmm1 = _mm_load_pd(&p_Ts[idx6]);
            xmm2 = _mm_load_pd(&p_M[idx6]);
            xmm3 = _mm_load_pd(&p_N0[idx6]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx6],xmm4);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[idx8],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[idx8],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[idx8], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[idx8],_MM_HINT_T0);
             _mm_prefetch((const char*)&p_Ac[i+L2_offset+idx8],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+idx8],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+idx8], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+idx8],_MM_HINT_T1);
#endif
            xmm5 = _mm_load_pd(&p_Ac[idx8]);
            xmm6 = _mm_load_pd(&p_Ts[idx8]);
            xmm7 = _mm_load_pd(&p_M[idx8]);
            xmm8 = _mm_load_pd(&p_N0[idx8]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx8],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx10]);
            xmm11 = _mm_load_pd(&p_Ts[idx10]);
            xmm12 = _mm_load_pd(&p_M[idx10]);
            xmm13 = _mm_load_pd(&p_N0[idx10]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx10],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx12]);
            xmm1 = _mm_load_pd(&p_Ts[idx12]);
            xmm2 = _mm_load_pd(&p_M[idx12]);
            xmm3 = _mm_load_pd(&p_N0[idx12]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx12],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx14]);
            xmm6 = _mm_load_pd(&p_Ts[idx14]);
            xmm7 = _mm_load_pd(&p_M[idx14]);
            xmm8 = _mm_load_pd(&p_N0[idx14]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx14],xmm9);

#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[idx16],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[idx16],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[idx16], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[idx16],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+idx16],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+idx16],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+idx16], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+idx16],_MM_HINT_T1);
#endif

            xmm10 = _mm_load_pd(&p_Ac[idx16]);
            xmm11 = _mm_load_pd(&p_Ts[idx16]);
            xmm12 = _mm_load_pd(&p_M[idx16]);
            xmm13 = _mm_load_pd(&p_N0[idx16]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx16],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx18]);
            xmm1 = _mm_load_pd(&p_Ts[idx18]);
            xmm2 = _mm_load_pd(&p_M[idx18]);
            xmm3 = _mm_load_pd(&p_N0[idx18]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx18],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx20]);
            xmm6 = _mm_load_pd(&p_Ts[idx20]);
            xmm7 = _mm_load_pd(&p_M[idx20]);
            xmm8 = _mm_load_pd(&p_N0[idx20]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx20],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx22]);
            xmm11 = _mm_load_pd(&p_Ts[idx22]);
            xmm12 = _mm_load_pd(&p_M[idx22]);
            xmm13 = _mm_load_pd(&p_N0[idx22]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx22],xmm14);
       }
 
       for(; (i+9ull) < a_sz; i += 10ull)
       {
            const std::size_t idx2  = i+2ull;
            const std::size_t idx4  = i+4ull;
            const std::size_t idx6  = i+6ull;
            const std::size_t idx8  = i+8ull;
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
             _mm_prefetch((const char*)&p_Ac[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+0ull], _MM_HIN1_T1);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+0ull],_MM_HINT_T1);
#endif 
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx2]);
            xmm6 = _mm_load_pd(&p_Ts[idx2]);
            xmm7 = _mm_load_pd(&p_M[idx2]);
            xmm8 = _mm_load_pd(&p_N0[idx2]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx2],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx4]);
            xmm11 = _mm_load_pd(&p_Ts[idx4]);
            xmm12 = _mm_load_pd(&p_M[idx4]);
            xmm13 = _mm_load_pd(&p_N0[idx4]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx4],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx6]);
            xmm1 = _mm_load_pd(&p_Ts[idx6]);
            xmm2 = _mm_load_pd(&p_M[idx6]);
            xmm3 = _mm_load_pd(&p_N0[idx6]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx6],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx8]);
            xmm6 = _mm_load_pd(&p_Ts[idx8]);
            xmm7 = _mm_load_pd(&p_M[idx8]);
            xmm8 = _mm_load_pd(&p_N0[idx8]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx8],xmm9);
       }

       for(; (i+3ull) < a_sz; i += 4ull)
       {
            const std::size_t idx2  = i+2ull;
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx2]);
            xmm6 = _mm_load_pd(&p_Ts[idx2]);
            xmm7 = _mm_load_pd(&p_M[idx2]);
            xmm8 = _mm_load_pd(&p_N0[idx2]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx2],xmm9);
       }

       for(; (i+1ull) < a_sz; i += 2ull)
       {
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);
       }

       for(; (i+0ull) < a_sz; i += 1ull)
       {
            const double valAc = p_Ac[i];
            const double valTs = p_Ts[i];
            const double valM  = p_M[i];
            const double valN0 = p_N0[i];
            p_sep_mam[i]      = scalar_kernel(valAc,valTs,valM,valN0,n);
       }
       return (0);
    }

    }
}

/////////////////////////////////////////////////////////////////////////////////

#if (ANALYTIC_BEP_SEP_CH8_SSE_ARRAY1D_F64_OVERRIDE_COMPILER_CMD_LINE) == 1
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
analytic_SEP_BEP_generic_2xf64_array1d_u4x(const double * __restrict__ Ac, 
                                        const double * __restrict__ Ts,
                                        const double * __restrict__ M,
                                        const double * __restrict__ N0,
                                        double       * __restrict__ sep_mam,
                                        const std::int32_t n,
                                        std::size_t sz,
                                        const callable_sse_kernel     & sse_kernel,
                                        const callable_scalar_kernel  & scalar_kernel) 
{
    if(__builtin_expect(sz<1ull,0)) {return (-1);}
    __m128d vAc;
    __m128d vTs;
    __m128d vM;
    __m128d vN0;
    const double * __restrict__ p_Ac      = Ac;
    const double * __restrict__ p_Ts      = Ts;
    const double * __restrict__ p_M       = M;
    const double * __restrict__ p_N0      = N0;
    double       * __restrict__ p_sep_mam = sep_mam;
    std::size_t             a_sz         = sz;
    if(a_sz>=1ull && a_sz<2ull)
    {
       for(std::uint32_t i = 0; i<2; ++i) 
       {
           const double valAc = p_Ac[i];
           const double valTs = p_Ts[i];
           const double valM  = p_M[i];
           const double valN0 = p_N0[i];
           p_sep_mam[i]      = scalar_kernel(valAc,valTs,valM,valN0,n);
       }
       return (0);
    }
    else if(a_sz==2ull)
    {
        vAc = _mm_loadu_pd(&p_Ac[0]);
        vTs = _mm_loadu_pd(&p_Ts[0]);
        vM  = _mm_loadu_pd(&p_M[0]);
        vN0 = _mm_loadu_pd(&p_N0[0]);
        _mm_storeu_pd(&p_sep_mam[0], sse_kernel(vAc,vTs,vM,vN0,n));
        return (0);
    }
    else if(a_sz>2ull && a_sz<=8ull)
    {
        std::size_t i,j;
        for(i = 0ull; i!=ROUND_DOWN(a_sz,1ull); i += 2ull)
        {
            vAc = _mm_loadu_pd(&p_Ac[i]);
            vTs = _mm_loadu_pd(&p_Ts[i]);
            vM  = _mm_loadu_pd(&p_M[i]);
            vN0 = _mm_loadu_pd(&p_N0[i]);
            _mm_storeu_pd(&p_sep_mam[i], sse_kernel(vAc,vTs,vM,vN0,n));
        }
        for(j = i; j<a_sz; ++j)  
        {
            const double valAc = p_Ac[j];
            const double valTs = p_Ts[j];
            const double valM  = p_M[j];
            const double valN0 = p_N0[j];
            p_sep_mam[j]      = scalar_kernel(valAc,valTs,valM,valN0,n);
        }
        return (0);
    }
    else if(a_sz>8ull)
    {
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_USE_PEEL_LOOP) == 1
       while(((std::uintptr_t)&p_sep_mam & 15ull) && a_sz)
       {
           const double valAc = *p_Ac;
           const double valTs = *p_Ts;
           const double valM  = *p_M;
           const double valN0 = *p_N0;
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
       __m128d xmm0;
       __m128d xmm1;
       __m128d xmm2;
       __m128d xmm3;
       __m128d xmm4;
       __m128d xmm5;
       __m128d xmm6;
       __m128d xmm7;
       __m128d xmm8;
       __m128d xmm9;
       __m128d xmm10;
       __m128d xmm11;
       __m128d xmm12;
       __m128d xmm13;
       __m128d xmm14;
       std::size_t i,j;
       constexpr std::size_t L2_offset = 4000ull;
       constexpr std::size_t L1D_boundary = L2_offset;
    if(a_sz <= L1D_boundary)
    {
       for(i = 0ull; (i+7ull) < a_sz; i += 8ull) 
       {
            const std::size_t idx2  = i+2ull;
            const std::size_t idx4  = i+4ull;
            const std::size_t idx6  = i+6ull;
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
#endif 
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx2]);
            xmm6 = _mm_load_pd(&p_Ts[idx2]);
            xmm7 = _mm_load_pd(&p_M[idx2]);
            xmm8 = _mm_load_pd(&p_N0[idx2]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx2],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx4]);
            xmm11 = _mm_load_pd(&p_Ts[idx4]);
            xmm12 = _mm_load_pd(&p_M[idx4]);
            xmm13 = _mm_load_pd(&p_N0[idx4]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx4],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx6]);
            xmm1 = _mm_load_pd(&p_Ts[idx6]);
            xmm2 = _mm_load_pd(&p_M[idx6]);
            xmm3 = _mm_load_pd(&p_N0[idx6]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx6],xmm4);
       }

       for(; (i+3ull) < a_sz; i += 4ull)
       {
            const std::size_t idx2  = i+2ull;
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx2]);
            xmm6 = _mm_load_pd(&p_Ts[idx2]);
            xmm7 = _mm_load_pd(&p_M[idx2]);
            xmm8 = _mm_load_pd(&p_N0[idx2]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx2],xmm9);
       }
 
       for(; (i+1ull) < a_sz; i += 2ull)
       {
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);
       }

       for(; (i+0ull) < a_sz; i += 1ull)
       {
            const double valAc = p_Ac[i];
            const double valTs = p_Ts[i];
            const double valM  = p_M[i];
            const double valN0 = p_N0[i];
            p_sep_mam[i]      = scalar_kernel(valAc,valTs,valM,valN0,n);
       }
       return (0);
    }
    else 
    {
       for(i = 0ull; (i+7ull) < a_sz; i += 8ull) 
       {
            const std::size_t idx2  = i+2ull;
            const std::size_t idx4  = i+4ull;
            const std::size_t idx6  = i+6ull;
#if (ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_SOFTWARE_PREFETCH) == 1
            _mm_prefetch((const char*)&p_Ac[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ts[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_M[i+0ull], _MM_HINT_T0);
            _mm_prefetch((const char*)&p_N0[i+0ull],_MM_HINT_T0);
            _mm_prefetch((const char*)&p_Ac[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_Ts[i+L2_offset+0ull],_MM_HINT_T1);
            _mm_prefetch((const char*)&p_M[i+L2_offset+0ull], _MM_HIN1_T0);
            _mm_prefetch((const char*)&p_N0[i+L2_offset+0ull],_MM_HINT_T1);
#endif 
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx2]);
            xmm6 = _mm_load_pd(&p_Ts[idx2]);
            xmm7 = _mm_load_pd(&p_M[idx2]);
            xmm8 = _mm_load_pd(&p_N0[idx2]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx2],xmm9);

            xmm10 = _mm_load_pd(&p_Ac[idx4]);
            xmm11 = _mm_load_pd(&p_Ts[idx4]);
            xmm12 = _mm_load_pd(&p_M[idx4]);
            xmm13 = _mm_load_pd(&p_N0[idx4]);
            xmm14 = sse_kernel(xmm10,xmm11,xmm12,xmm13,n);
            _mm_store_pd(&p_sep_mam[idx4],xmm14);

            xmm0 = _mm_load_pd(&p_Ac[idx6]);
            xmm1 = _mm_load_pd(&p_Ts[idx6]);
            xmm2 = _mm_load_pd(&p_M[idx6]);
            xmm3 = _mm_load_pd(&p_N0[idx6]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[idx6],xmm4);

       }

       for(; (i+3ull) < a_sz; i += 4ull)
       {
            const std::size_t idx2  = i+2ull;
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);

            xmm5 = _mm_load_pd(&p_Ac[idx2]);
            xmm6 = _mm_load_pd(&p_Ts[idx2]);
            xmm7 = _mm_load_pd(&p_M[idx2]);
            xmm8 = _mm_load_pd(&p_N0[idx2]);
            xmm9 = sse_kernel(xmm5,xmm6,xmm7,xmm8,n);
            _mm_store_pd(&p_sep_mam[idx2],xmm9);
       }

       for(; (i+1ull) < a_sz; i += 2ull)
       {
            xmm0 = _mm_load_pd(&p_Ac[i+0ull]);
            xmm1 = _mm_load_pd(&p_Ts[i+0ull]);
            xmm2 = _mm_load_pd(&p_M[i+0ull]);
            xmm3 = _mm_load_pd(&p_N0[i+0ull]);
            xmm4 = sse_kernel(xmm0,xmm1,xmm2,xmm3,n);
            _mm_store_pd(&p_sep_mam[i+0ull],xmm4);
       }

       for(; (i+0ull) < a_sz; i += 1ull)
       {
            const double valAc = p_Ac[i];
            const double valTs = p_Ts[i];
            const double valM  = p_M[i];
            const double valN0 = p_N0[i];
            p_sep_mam[i]      = scalar_kernel(valAc,valTs,valM,valN0,n);
       }
       return (0);
    }

    }
}


}// fading_channel

}// gms
#endif /*__GMS_ANALYTIC_BEP_SEP_CH8_ARRAY1D_F64_H__*/