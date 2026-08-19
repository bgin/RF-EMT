

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
#include <immintrin.h> // __rdtsc
#include <functional>
#include <algorithm>
#include <random>
#include <cstdio>
#include "GMS_compute_functionals_ch4.h"
#include "GMS_cquadpack.h"


/*
    double (*integrand)(double,void * __restrict__);
    func_args_payload_t * __restrict__ func_args_payload{nullptr};
    double              * __restrict__ lo{nullptr}; // lower limit of integration
    double              * __restrict__ hi{nullptr}; // upper limit of integration
    double              * __restrict__ bound{nullptr}; //optional finite bound on integral.
    std::int32_t        * __restrict__ inf{nullptr}; // range of integration
    std::int32_t        * __restrict__ irule{nullptr};
    double              * __restrict__ epsabs{nullptr};
    double              * __restrict__ epsrel{nullptr};
    double              * __restrict__ abser{nullptr}; // integrator abosulute error
    double              * __restrict__ functional{nullptr}; // computed functional values
    std::int32_t        * __restrict__ neval{nullptr};      // number of evaluation
    std::int32_t        * __restrict__ ier{nullptr};        // integrator error indicator
    std::int32_t        * __restrict__ last{nullptr};
    std::int32_t                       n_func_vals;
    std::int32_t                       which_integrator; currently 1=dqage,2=dqagi
*/

std::int32_t 
gms 
::fading_channel
::compute_functional_gauss_Q_4_1(quadpack_integrator_payload_t * __restrict__ p_payload)
{
    if(__builtin_expect(nullptr==p_payload,0)) { return (-1);}
    double (*p_integrand)(const double,void * __restrict__) = p_payload->integrand;
    func_args_payload_t  * __restrict__ p_funcs_args_payload = p_payload->func_args_payload;
    double  * __restrict__      p_lo                         = p_payload->lo;
    double  * __restrict__      p_hi                         = p_payload->hi;
    double  * __restrict__      p_bound                      = p_payload->bound;
    std::int32_t * __restrict__ p_inf                        = p_payload->inf;
    std::int32_t * __restrict__ p_irule                      = p_payload->irule;
    double  *  __restrict__     p_epsabs                     = p_payload->epsabs;
    double  *  __restrict__     p_epsrel                     = p_payload->epsrel;
    double  *  __restrict__     p_abser                      = p_payload->abser;
    double  *  __restrict__     p_functional                 = p_payload->functional;
    std::int32_t * __restrict__ p_neval                      = p_payload->neval;
    std::int32_t * __restrict__ p_ier                        = p_payload->ier;
    std::int32_t * __restrict__ p_last                       = p_payload->last;
    const std::int32_t          nfunc_vals                   = p_payload->n_func_vals;
    const std::int32_t          integrator_type              = p_payload->which_integrator;
    thread_local std::uniform_real_distribution<double> rv_func_lo;
    thread_local std::mt19937 rv_func_lo_gen;
    thread_local std::uint64_t seed_func_lo{};
    if(integrator_type==1)
    {
        thread_local std::uniform_real_distribution<double> rv_func_hi;
        thread_local std::mt19937 rv_func_hi_gen;
        thread_local std::uint64_t seed_func_hi{};
        rv_func_lo = std::uniform_real_distribution<double>(-4.0,+4.0);
        seed_func_lo = __rdtsc();
        rv_func_lo_gen = std::mt19937(seed_func_lo);
        rv_func_hi = std::uniform_real_distribution<double>(+5.0,+100.0);
        seed_func_hi = __rdtsc();
        rv_func_hi_gen = std::mt19937(seed_func_hi);
        for(std::int32_t i{0}; i<nfunc_vals; ++i)
        {
            const double lo_limit{rv_func_lo.operator()(rv_func_lo_gen)};
            p_lo[i] = lo_limit;
            const double hi_limit{rv_func_hi.operator()(rv_func_hi_gen)};
            p_hi[i] = hi_limit;
        }
        std::sort(&p_lo[0],&p_lo[nfunc_vals-1],std::greater<double>());
        std::sort(&p_hi[0],&p_hi[nfunc_vals-1],std::less<double>());
        for(std::int32_t i{0}; i<nfunc_vals; ++i) 
        {   
            const double a{p_lo[i]};
            const double b{p_hi[i]};
            const double result = dqage(p_integrand,a,b,p_epsabs[0],p_epsrel[0],p_irule[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_last[i],nullptr);
            p_functional[i] = result;
        }
    }
    else if(integrator_type==2)
    {
        rv_func_lo = std::uniform_real_distribution<double>(-4.0,+10.0);
        seed_func_lo = __rdtsc();
        rv_func_lo_gen = std::mt19937(seed_func_lo);
        for(std::int32_t i{0}; i<nfunc_vals; ++i)
        {
            const double lo_limit{rv_func_lo.operator()(rv_func_lo_gen)};
            p_lo[i] = lo_limit;
        }
        std::sort(&p_lo[0],&p_lo[nfunc_vals-1],std::greater<double>());
        for(std::int32_t i{0}; i<nfunc_vals; ++i) 
        {
            const double bound{p_lo[i]};
            const double result = dqagi(p_integrand,bound,p_inf[0],p_epsabs[0],p_epsrel[0],
                                        &p_abser[i],&p_neval[i],&p_ier[i],nullptr);
            p_functional[i] = result;
        }
    }
    return (0);
}