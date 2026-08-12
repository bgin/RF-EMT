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

#include <immintrin.h>
#include <cmath>
#include <algorithm>
#include <functional>
#include <random>
#include <cstdio>
#include "GMS_functionals_ch4_sse.h"
#include "GMS_integrands_func_ch4_sse.h"
#include "GMS_integrands_func_ch4.h"
#include "GMS_simpne_quad.h"

namespace 
{
static inline 
std::int32_t 
print_double(const char * val_name,double f, const std::int32_t offset) 
{
  (*((std::uint64_t *)&f)) += offset;
  std::int32_t printf_ret = std::printf("value: %s,%1.17g,%017llx,%llu\n",val_name,f, *(std::uint64_t *)&f, *(std::uint64_t *)&f);
  return (printf_ret);
}
}

#define FUNCTIONALS_CH4_SSE_COMPUTE_BODY_DBG(func_name,j,loop_incr,functional_table)\
{\
  for(std::int32_t i{0}; i < tot_elems; i += (loop_incr))\
  {\
        double * __restrict__ p_slice_out_buf{&p_out_buf[i]};\
        double * __restrict__ p_slice_in_buf{&p_in_buf[i]};\
        __asm__ volatile("int3");\
        (void)math::simpne(ntab,&p_slice_in_buf[0],&p_slice_out_buf[0],(functional_table)[(j)]);\
        [[maybe_unused]] std::int32_t printf_ret = print_double((func_name),(functional_table)[(j)],0);\
        ++(j);\
    }\
}

#define FUNCTIONALS_CH4_SSE_COMPUTE_BODY(j,loop_incr,functional_table)\
{\
  for(std::int32_t i{0}; i < tot_elems; i += (loop_incr))\
  {\
        double * __restrict__ p_slice_out_buf{&p_out_buf[i]};\
        double * __restrict__ p_slice_in_buf{&p_in_buf[i]};\
        (void)math::simpne(ntab,&p_slice_in_buf[0],&p_slice_out_buf[0],(functional_table)[(j)]);\
        ++(j);\
    }\
}

std::int32_t 
gms::fading_channel
::integrate_4_1_gauss_Q_func_sse(integrator_payload_sse_t * __restrict__ payload)
{
    using namespace math;
    if(__builtin_expect(nullptr==payload,0)) { return (-1);}
    double * __restrict__ p_in_buf      = payload->in_buf;
    double * __restrict__ p_out_buf     = payload->out_buf;
    double * __restrict__ p_functional  = payload->functional;
    double                lo1_x         = payload->lo1;
    double                hi1_x         = payload->hi2;
    double                lo2_x         = lo1_x;
    double                hi2_x         = payload->hi2;
    std::int32_t          nfunc_args    = payload->n_func_args;
    std::int32_t          ngauss_q_vals = payload->n_integrand_vals;
    const std::int32_t    tot_elems{nfunc_args*ngauss_q_vals};
    thread_local std::uniform_real_distribution<double> rv_func_arg;
    thread_local std::mt19937 rv_func_arg_gen;
    thread_local std::uint64_t seed_func_arg{};
    thread_local std::uniform_real_distribution<double> rv_gauss_q;
    thread_local std::mt19937 rv_gauss_q_gen;
    thread_local std::uint64_t seed_gauss_q{};
    rv_func_arg = std::uniform_real_distribution<double>(lo1_x,hi1_x);
    seed_func_arg = __rdtsc();
    rv_func_arg_gen = std::mt19937(seed_func_arg);
    rv_gauss_q = std::uniform_real_distribution<double>(lo2_x,hi2_x);
    seed_gauss_q = __rdtsc();
    rv_gauss_q_gen = std::mt19937(seed_gauss_q);
    for(std::int32_t i {0}; i < nfunc_args; ++i) 
    {
        const std::int32_t outer_idx = i*ngauss_q_vals;
        const double x = rv_func_arg.operator()(rv_func_arg_gen);
        for(std::int32_t j{0}; j < ngauss_q_vals; ++j)   
        {
            const std::int32_t inner_idx = outer_idx+j;
            const double y = rv_gauss_q.operator()(rv_gauss_q_gen);
            p_in_buf[inner_idx] = x+y;
        }
    }
    std::sort(&p_in_buf[0],&p_in_buf[tot_elems-1]);
    __m128d vx0;
    __m128d vx1;
    __m128d vx2;
    __m128d vx3;
    __m128d gauss_q_res0;
    __m128d gauss_q_res1;
    __m128d gauss_q_res2;
    __m128d gauss_q_res3;
    std::int32_t jj;
    for(jj = 0; (jj+7) < tot_elems; jj+=8)
    {
#if (FUNCTIONALS_CH4_SSE_INTERMIX_LOAD_COMPUTE) == 0
        vx0 = _mm_load_pd(&p_in_buf[jj+0]);
        vx1 = _mm_load_pd(&p_in_buf[jj+2]);
        vx2 = _mm_load_pd(&p_in_buf[jj+4]);
        vx3 = _mm_load_pd(&p_in_buf[jj+6]);
        gauss_q_res0 = integrand_4_1_gauss_Q_func_sse_pd(vx1);
        gauss_q_res1 = integrand_4_1_gauss_Q_func_sse_pd(vx2);
        gauss_q_res2 = integrand_4_1_gauss_Q_func_sse_pd(vx2);
        gauss_q_res3 = integrand_4_1_gauss_Q_func_sse_pd(vx3);
        _mm_store_pd(&p_out_buf[jj+0],gauss_q_res0);
        _mm_store_pd(&p_out_buf[jj+2],gauss_q_res1);
        _mm_store_pd(&p_out_buf[jj+4],gauss_q_res2);
        _mm_store_pd(&p_out_buf[jj+6],gauss_q_res3);
#else 
        vx0 = _mm_load_pd(&p_in_buf[jj+0]);
        gauss_q_res0 = integrand_4_1_gauss_Q_func_sse_pd(vx0);
        _mm_store_pd(&p_out_buf[jj+0],gauss_q_res0);
        vx1 = _mm_load_pd(&p_in_buf[jj+2]);
        gauss_q_res1 = integrand_4_1_gauss_Q_func_sse_pd(vx1);
        _mm_store_pd(&p_out_buf[jj+2],gauss_q_res1);
        vx2 = _mm_load_pd(&p_in_buf[jj+4]);
        gauss_q_res2 = integrand_4_1_gauss_Q_func_sse_pd(vx2);
        _mm_store_pd(&p_out_buf[jj+4],gauss_q_res2);
        vx3 = _mm_load_pd(&p_in_buf[jj+6]);
        gauss_q_res3 = integrand_4_1_gauss_Q_func_sse_pd(vx3);
        _mm_store_pd(&p_out_buf[jj+6],gauss_q_res3);
#endif
    }
    
    for(; (jj+3) < tot_elems; jj+=4)
    {
#if (FUNCTIONALS_CH4_SSE_INTERMIX_LOAD_COMPUTE) == 0
        vx0 = _mm_load_pd(&p_in_buf[jj+0]);
        vx1 = _mm_load_pd(&p_in_buf[jj+2]);
        gauss_q_res0 = integrand_4_1_gauss_Q_func_sse_pd(vx1);
        gauss_q_res1 = integrand_4_1_gauss_Q_func_sse_pd(vx2);
        _mm_store_pd(&p_out_buf[jj+0],gauss_q_res0);
        _mm_store_pd(&p_out_buf[jj+2],gauss_q_res1);
#else 
        vx0 = _mm_load_pd(&p_in_buf[jj+0]);
        gauss_q_res0 = integrand_4_1_gauss_Q_func_sse_pd(vx0);
        _mm_store_pd(&p_out_buf[jj+0],gauss_q_res0);
        vx1 = _mm_load_pd(&p_in_buf[jj+2]);
        gauss_q_res1 = integrand_4_1_gauss_Q_func_sse_pd(vx1);
        _mm_store_pd(&p_out_buf[jj+2],gauss_q_res1);
#endif
    }
    
    for(; (jj+1) < tot_elems; jj+=2)
    {
        vx0 = _mm_load_pd(&p_in_buf[jj+0]);
        gauss_q_res0 = integrand_4_1_gauss_Q_func_sse_pd(vx0);
        _mm_store_pd(&p_out_buf[jj+0],gauss_q_res0);
    }

    for(; (jj+0) < tot_elems; jj+=1)
    {
        const double x = p_in_buf[jj];
        const double gauss_q_res = integrand_4_1_gauss_Q_func(x);
        p_out_buf[jj] = gauss_q_res;
    }
#if 0
     for(std::int32_t i{0}; i<tot_elems; ++i) 
     {
        print_double("**in_buf**=",p_in_buf[i],0);
        print_double("**out_buf**=",p_out_buf[i],0);
     }
#endif
    const std::int32_t ntab = ngauss_q_vals;
    std::int32_t j{0};
    FUNCTIONALS_CH4_SSE_COMPUTE_BODY(j,ngauss_q_vals,p_functional);
    return (0);
}



