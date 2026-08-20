

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
#include "GMS_machine_utils.h" // for RDTSCP wrapper


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
    std::int32_t                       which_integrator; currently 1=dqage,2=dqagi,3=dqags
*/

std::int32_t 
gms 
::fading_channel
::compute_functional_gauss_Q_4_1(quadpack_integrator_payload_t * __restrict__ p_payload)
{
    if(__builtin_expect(nullptr==p_payload,0)) { return (-1);}
    double (*p_integrand)(const double,void * __restrict__)  = p_payload->integrand;
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
    std::uint64_t * __restrict__ p_crude_tsc_meter           = p_payload->crude_tsc_measurement;
    double                      rand_low1                    = p_payload->rand_lo1;
    double                      rand_high1                   = p_payload->rand_hi1;
    double                      rand_low2                    = p_payload->rand_lo2;
    double                      rand_high2                   = p_payload->rand_hi2;
    double                      rand_low3                    = p_payload->rand_lo3;
    double                      rand_high3                   = p_payload->rand_hi3;
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
        rv_func_lo = std::uniform_real_distribution<double>(rand_low1,rand_high1);
        seed_func_lo = __rdtsc();
        rv_func_lo_gen = std::mt19937(seed_func_lo);
        rv_func_hi = std::uniform_real_distribution<double>(rand_low2,rand_high2);
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
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqage(p_integrand,a,b,p_epsabs[0],p_epsrel[0],p_irule[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_last[i],nullptr);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = result;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    else if(integrator_type==2)
    {
        rv_func_lo = std::uniform_real_distribution<double>(rand_low3,rand_high3);
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
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqagi(p_integrand,bound,p_inf[0],p_epsabs[0],p_epsrel[0],
                                        &p_abser[i],&p_neval[i],&p_ier[i],nullptr);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = result;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    else if(integrator_type==3)
    {
        thread_local std::uniform_real_distribution<double> rv_func_hi;
        thread_local std::mt19937 rv_func_hi_gen;
        thread_local std::uint64_t seed_func_hi{};
        rv_func_lo = std::uniform_real_distribution<double>(rand_low1,rand_high1);
        seed_func_lo = __rdtsc();
        rv_func_lo_gen = std::mt19937(seed_func_lo);
        rv_func_hi = std::uniform_real_distribution<double>(rand_low2,rand_high2);
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
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqags(p_integrand,a,b,p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],nullptr);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = result;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    return (0);
}


std::int32_t 
gms 
::fading_channel
::compute_functional_gauss_Q_4_2(quadpack_integrator_payload_t * __restrict__ p_payload)
{
    if(__builtin_expect(nullptr==p_payload,0)) { return (-1);}
    double (*p_integrand)(const double,void * __restrict__)  = p_payload->integrand;
    func_args_payload_t  * __restrict__ p_funcs_args_payload = p_payload->func_args_payload;
    double       * __restrict__ p_tmp_work                   = p_payload->tmp_work;
    std::int32_t * __restrict__ p_inf                        = p_payload->inf;
    std::int32_t * __restrict__ p_irule                      = p_payload->irule;
    double  *  __restrict__     p_epsabs                     = p_payload->epsabs;
    double  *  __restrict__     p_epsrel                     = p_payload->epsrel;
    double  *  __restrict__     p_abser                      = p_payload->abser;
    double  *  __restrict__     p_functional                 = p_payload->functional;
    std::int32_t * __restrict__ p_neval                      = p_payload->neval;
    std::int32_t * __restrict__ p_ier                        = p_payload->ier;
    std::int32_t * __restrict__ p_last                       = p_payload->last;
    std::uint64_t * __restrict__ p_crude_tsc_meter           = p_payload->crude_tsc_measurement;
    double                      rand_low1                    = p_payload->rand_lo1;
    double                      rand_high1                   = p_payload->rand_hi1;
    double                      rand_low2                    = p_payload->rand_lo2;
    double                      rand_high2                   = p_payload->rand_hi2;
    const std::int32_t          nfunc_vals                   = p_payload->n_func_vals;
    const std::int32_t          integrator_type              = p_payload->which_integrator;
    thread_local std::uniform_real_distribution<double> rv_func_x;
    thread_local std::mt19937 rv_func_x_gen;
    thread_local std::uint64_t seed_func_x{};
    rv_func_x = std::uniform_real_distribution<double>(rand_low2,rand_high2);
    seed_func_x = __rdtsc();
    rv_func_x_gen = std::mt19937(seed_func_x);
    for(std::int32_t i{0}; i<nfunc_vals; ++i)
    {
            const double x_val{rv_func_x.operator()(rv_func_x_gen)};
            p_tmp_work[i] = x_val;
    }
    std::sort(&p_tmp_work[0],&p_tmp_work[nfunc_vals-1],std::less<double>());
#if 0
    for(std::int32_t i{0}; i<nfunc_vals; ++i) { std::int32_t ret = std::printf("LOC=%d,arg1=%.17f\n",__LINE__,p_tmp_work[i]);}
#endif
    if(integrator_type==1)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i) 
        {   
            const double cpy_x{p_tmp_work[i]};
            p_funcs_args_payload[i].arg1 = cpy_x;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqage(p_integrand,rand_low1,rand_high1,p_epsabs[0],p_epsrel[0],p_irule[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_last[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = 0.318309886183790671537767527*result;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    else if(integrator_type==2)
    {
        return (-2);
    }
    else if(integrator_type==3)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i) 
        {   
            const double cpy_x{p_tmp_work[i]};
            p_funcs_args_payload[i].arg1 = cpy_x;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqags(p_integrand,rand_low1,rand_high1,p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = 0.318309886183790671537767527*result;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    return (0);
}


/*
double 
gms::fading_channel::test_gauss_Q_4_2(double * ptr) 
{
    return (3.14+ptr[0]);
}

std::int32_t 
compute_functional_gauss_Q_4_2(gms::fading_channel::quadpack_integrator_payload_t * __restrict__ p_payload)
{
    return (-1);
}
*/