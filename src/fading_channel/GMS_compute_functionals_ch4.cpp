

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
#include <cstdio> // for debugging
#include <cmath>
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
    else if(integrator_type==4)
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
            const double result = dqng(p_integrand,a,b,p_epsabs[0],p_epsrel[0],&p_abser[i],
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
    double       * __restrict__ p_tmp_work                   = p_payload->tmp_work1;
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
    else if(integrator_type==4)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i) 
        {   
            const double cpy_x{p_tmp_work[i]};
            p_funcs_args_payload[i].arg1 = cpy_x;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqng(p_integrand,rand_low1,rand_high1,p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = 0.318309886183790671537767527*result;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    return (0);
}

std::int32_t 
gms::fading_channel
::compute_functional_gauss_Q_4_6(quadpack_integrator_payload_t * __restrict__ p_payload)
{
    if(__builtin_expect(nullptr==p_payload,0)) { return (-1);}
    double (*p_integrand)(const double,void * __restrict__)  = p_payload->integrand;
    func_args_payload_t  * __restrict__ p_funcs_args_payload = p_payload->func_args_payload;
    double       * __restrict__ p_tmp_work1                  = p_payload->tmp_work1;
    double       * __restrict__ p_tmp_work2                  = p_payload->tmp_work2;
    double       * __restrict__ p_tmp_work3                  = p_payload->tmp_work3;
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
    thread_local std::uniform_real_distribution<double> rv_func_x1;
    thread_local std::mt19937 rv_func_x1_gen;
    thread_local std::uint64_t seed_func_x1{};
    thread_local std::uniform_real_distribution<double> rv_func_y1;
    thread_local std::mt19937 rv_func_y1_gen;
    thread_local std::uint64_t seed_func_y1{};
    thread_local std::uniform_real_distribution<double> rv_func_rho;
    thread_local std::mt19937 rv_func_rho_gen;
    thread_local std::uint64_t seed_func_rho{};
    rv_func_x1 = std::uniform_real_distribution<double>(rand_low1,rand_high1);
    seed_func_x1 = __rdtsc();
    rv_func_x1_gen = std::mt19937(seed_func_x1);
    rv_func_y1 = std::uniform_real_distribution<double>(rand_low2,rand_high2);
    seed_func_y1 = __rdtsc();
    rv_func_y1_gen = std::mt19937(seed_func_y1);
    rv_func_rho = std::uniform_real_distribution<double>(rand_low3,rand_high3);
    seed_func_rho = __rdtsc();
    rv_func_rho_gen = std::mt19937(seed_func_rho);
    for(std::int32_t i{0}; i<nfunc_vals; ++i)
    {
        const double val_x1{rv_func_x1.operator()(rv_func_x1_gen)};
        p_tmp_work1[i] = val_x1;
        const double val_y1{rv_func_y1.operator()(rv_func_y1_gen)};
        p_tmp_work2[i] = val_y1;
        const double val_rho{rv_func_rho.operator()(rv_func_rho_gen)};
        p_tmp_work3[i] = val_rho;

    }
    std::sort(&p_tmp_work1[0],&p_tmp_work1[nfunc_vals-1],std::less<double>());
    std::sort(&p_tmp_work2[0],&p_tmp_work2[nfunc_vals-1],std::less<double>());
    std::sort(&p_tmp_work3[0],&p_tmp_work3[nfunc_vals-1],std::less<double>());
    if(integrator_type==1)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i)
        {
            const double cpy_x1{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1 = cpy_x1;
            const double cpy_y1{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg2 = cpy_y1;
            const double phi_arg{cpy_y1/cpy_x1};
            const double cpy_rho{p_tmp_work3[i]};
            p_funcs_args_payload[i].arg3 = cpy_rho;
            const double phi_s{std::atan(phi_arg)};
            p_funcs_args_payload[i].arg5 = 1;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result1 = dqage(p_integrand,0.0,(1.570796326794896619231321692-phi_s),p_epsabs[0],p_epsrel[0],p_irule[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_last[i],&p_funcs_args_payload[i]);
            p_funcs_args_payload[i].arg5 = 2;
            const double result2 = dqage(p_integrand,0.0,phi_s,p_epsabs[0],p_epsrel[0],p_irule[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_last[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            const double full_result1{0.159154943091895335768883763*result1};
            const double full_result2{0.159154943091895335768883763*result2};
            p_functional[i] = full_result1+full_result2;
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
            const double cpy_x1{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1 = cpy_x1;
            const double cpy_y1{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg2 = cpy_y1;
            const double phi_arg{cpy_y1/cpy_x1};
            const double cpy_rho{p_tmp_work3[i]};
            p_funcs_args_payload[i].arg3 = cpy_rho;
            const double phi_s{std::atan(phi_arg)};
            p_funcs_args_payload[i].arg5 = 1;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result1 = dqags(p_integrand,0.0,(1.570796326794896619231321692-phi_s),p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            p_funcs_args_payload[i].arg5 = 2;
            const double result2 = dqags(p_integrand,0.0,phi_s,p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            const double full_result1{0.159154943091895335768883763*result1};
            const double full_result2{0.159154943091895335768883763*result2};
            p_functional[i] = full_result1+full_result2;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    else if(integrator_type==4)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i) 
        {
            const double cpy_x1{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1 = cpy_x1;
            const double cpy_y1{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg2 = cpy_y1;
            const double phi_arg{cpy_y1/cpy_x1};
            const double cpy_rho{p_tmp_work3[i]};
            p_funcs_args_payload[i].arg3 = cpy_rho;
            const double phi_s{std::atan(phi_arg)};
            p_funcs_args_payload[i].arg5 = 1;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result1 = dqng(p_integrand,0.0,(1.570796326794896619231321692-phi_s),p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            p_funcs_args_payload[i].arg5 = 2;
            const double result2 = dqng(p_integrand,0.0,phi_s,p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            const double full_result1{0.159154943091895335768883763*result1};
            const double full_result2{0.159154943091895335768883763*result2};
            p_functional[i] = full_result1+full_result2;
            p_crude_tsc_meter[i] = end-start;
        }
    }

    return (0);
}

std::int32_t 
gms::fading_channel
::compute_functional_gauss_Q_4_7(quadpack_integrator_payload_t * __restrict__ p_payload)
{
    if(__builtin_expect(nullptr==p_payload,0)) { return (-1);}
    double (*p_integrand)(const double,void * __restrict__)  = p_payload->integrand;
    func_args_payload_t  * __restrict__ p_funcs_args_payload = p_payload->func_args_payload;
    double       * __restrict__ p_tmp_work1                  = p_payload->tmp_work1;
    double       * __restrict__ p_tmp_work2                  = p_payload->tmp_work2;
    double       * __restrict__ p_tmp_work3                  = p_payload->tmp_work3;
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
    thread_local std::uniform_real_distribution<double> rv_func_x1;
    thread_local std::mt19937 rv_func_x1_gen;
    thread_local std::uint64_t seed_func_x1{};
    thread_local std::uniform_real_distribution<double> rv_func_y1;
    thread_local std::mt19937 rv_func_y1_gen;
    thread_local std::uint64_t seed_func_y1{};
    thread_local std::uniform_real_distribution<double> rv_func_rho;
    thread_local std::mt19937 rv_func_rho_gen;
    thread_local std::uint64_t seed_func_rho{};
    rv_func_x1 = std::uniform_real_distribution<double>(rand_low1,rand_high1);
    seed_func_x1 = __rdtsc();
    rv_func_x1_gen = std::mt19937(seed_func_x1);
    rv_func_y1 = std::uniform_real_distribution<double>(rand_low2,rand_high2);
    seed_func_y1 = __rdtsc();
    rv_func_y1_gen = std::mt19937(seed_func_y1);
    rv_func_rho = std::uniform_real_distribution<double>(rand_low3,rand_high3);
    seed_func_rho = __rdtsc();
    rv_func_rho_gen = std::mt19937(seed_func_rho);
    for(std::int32_t i{0}; i<nfunc_vals; ++i)
    {
        const double val_x1{rv_func_x1.operator()(rv_func_x1_gen)};
        p_tmp_work1[i] = val_x1;
        const double val_y1{rv_func_y1.operator()(rv_func_y1_gen)};
        p_tmp_work2[i] = val_y1;
        const double val_rho{rv_func_rho.operator()(rv_func_rho_gen)};
        p_tmp_work3[i] = val_rho;

    }
    std::sort(&p_tmp_work1[0],&p_tmp_work1[nfunc_vals-1],std::less<double>());
    std::sort(&p_tmp_work2[0],&p_tmp_work2[nfunc_vals-1],std::less<double>());
    std::sort(&p_tmp_work3[0],&p_tmp_work3[nfunc_vals-1],std::less<double>());
    if(integrator_type==1)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i)
        {
            const double cpy_x1{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1 = cpy_x1;
            const double cpy_y1{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg2 = cpy_y1;
            const double atan_arg{cpy_y1/cpy_x1};
            const double cpy_rho{p_tmp_work3[i]};
            p_funcs_args_payload[i].arg3 = cpy_rho;
            const double atan_val{std::atan(atan_arg)};
            p_funcs_args_payload[i].arg5 = 1;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result1 = dqage(p_integrand,0.0,(1.570796326794896619231321692-atan_val),p_epsabs[0],p_epsrel[0],p_irule[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_last[i],&p_funcs_args_payload[i]);
            p_funcs_args_payload[i].arg5 = 2;
            const double result2 = dqage(p_integrand,0.0,atan_val,p_epsabs[0],p_epsrel[0],p_irule[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_last[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            const double full_result1{0.159154943091895335768883763*result1};
            const double full_result2{0.159154943091895335768883763*result2};
            p_functional[i] = full_result1+full_result2;
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
            const double cpy_x1{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1 = cpy_x1;
            const double cpy_y1{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg2 = cpy_y1;
            const double atan_arg{cpy_y1/cpy_x1};
            const double cpy_rho{p_tmp_work3[i]};
            p_funcs_args_payload[i].arg3 = cpy_rho;
            const double atan_val{std::atan(atan_arg)};
            p_funcs_args_payload[i].arg5 = 1;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result1 = dqags(p_integrand,0.0,(1.570796326794896619231321692-atan_val),p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            p_funcs_args_payload[i].arg5 = 2;
            const double result2 = dqags(p_integrand,0.0,atan_val,p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            const double full_result1{0.159154943091895335768883763*result1};
            const double full_result2{0.159154943091895335768883763*result2};
            p_functional[i] = full_result1+full_result2;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    else if(integrator_type==4)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i) 
        {
            const double cpy_x1{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1 = cpy_x1;
            const double cpy_y1{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg2 = cpy_y1;
            const double atan_arg{cpy_y1/cpy_x1};
            const double cpy_rho{p_tmp_work3[i]};
            p_funcs_args_payload[i].arg3 = cpy_rho;
            const double atan_val{std::atan(atan_arg)};
            p_funcs_args_payload[i].arg5 = 1;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result1 = dqng(p_integrand,0.0,(1.570796326794896619231321692-atan_val),p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            p_funcs_args_payload[i].arg5 = 2;
            const double result2 = dqng(p_integrand,0.0,atan_val,p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            const double full_result1{0.159154943091895335768883763*result1};
            const double full_result2{0.159154943091895335768883763*result2};
            p_functional[i] = full_result1+full_result2;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    return (0);
}

std::int32_t 
gms::fading_channel
::compute_functional_gauss_Q_4_8(quadpack_integrator_payload_t * __restrict__ p_payload)
{
    if(__builtin_expect(nullptr==p_payload,0)) { return (-1);}
    double (*p_integrand)(const double,void * __restrict__)  = p_payload->integrand;
    func_args_payload_t  * __restrict__ p_funcs_args_payload = p_payload->func_args_payload;
    double       * __restrict__ p_tmp_work1                  = p_payload->tmp_work1;
    double       * __restrict__ p_tmp_work2                  = p_payload->tmp_work2;
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
    thread_local std::uniform_real_distribution<double> rv_func_x1;
    thread_local std::mt19937 rv_func_x1_gen;
    thread_local std::uint64_t seed_func_x1{};
    thread_local std::uniform_real_distribution<double> rv_func_y1;
    thread_local std::mt19937 rv_func_y1_gen;
    thread_local std::uint64_t seed_func_y1{};
    rv_func_x1 = std::uniform_real_distribution<double>(rand_low1,rand_high1);
    seed_func_x1 = __rdtsc();
    rv_func_x1_gen = std::mt19937(seed_func_x1);
    rv_func_y1 = std::uniform_real_distribution<double>(rand_low2,rand_high2);
    seed_func_y1 = __rdtsc();
    rv_func_y1_gen = std::mt19937(seed_func_y1);
    for(std::int32_t i{0}; i<nfunc_vals; ++i)
    {
        const double val_x1{rv_func_x1.operator()(rv_func_x1_gen)};
        p_tmp_work1[i] = val_x1;
        const double val_y1{rv_func_y1.operator()(rv_func_y1_gen)};
        p_tmp_work2[i] = val_y1;

    }
    std::sort(&p_tmp_work1[0],&p_tmp_work1[nfunc_vals-1],std::less<double>());
    std::sort(&p_tmp_work2[0],&p_tmp_work2[nfunc_vals-1],std::less<double>());
    if(integrator_type==1)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i)
        {
            const double cpy_x1{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1 = cpy_x1;
            const double cpy_y1{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg2 = cpy_y1;
            const double atan_arg{cpy_y1/cpy_x1};
            const double atan_val{std::atan(atan_arg)};
            p_funcs_args_payload[i].arg5 = 1;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result1 = dqage(p_integrand,0.0,(1.570796326794896619231321692-atan_val),p_epsabs[0],p_epsrel[0],p_irule[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_last[i],&p_funcs_args_payload[i]);
            p_funcs_args_payload[i].arg5 = 2;
            const double result2 = dqage(p_integrand,0.0,atan_val,p_epsabs[0],p_epsrel[0],p_irule[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_last[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            const double full_result1{0.159154943091895335768883763*result1};
            const double full_result2{0.159154943091895335768883763*result2};
            p_functional[i] = full_result1+full_result2;
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
            const double cpy_x1{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1 = cpy_x1;
            const double cpy_y1{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg2 = cpy_y1;
            const double atan_arg{cpy_y1/cpy_x1};
            const double atan_val{std::atan(atan_arg)};
            p_funcs_args_payload[i].arg5 = 1;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result1 = dqags(p_integrand,0.0,(1.570796326794896619231321692-atan_val),p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            p_funcs_args_payload[i].arg5 = 2;
            const double result2 = dqags(p_integrand,0.0,atan_val,p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            const double full_result1{0.159154943091895335768883763*result1};
            const double full_result2{0.159154943091895335768883763*result2};
            p_functional[i] = full_result1+full_result2;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    else if(integrator_type==4)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i) 
        {
            const double cpy_x1{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1 = cpy_x1;
            const double cpy_y1{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg2 = cpy_y1;
            const double atan_arg{cpy_y1/cpy_x1};
            const double atan_val{std::atan(atan_arg)};
            p_funcs_args_payload[i].arg5 = 1;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result1 = dqng(p_integrand,0.0,(1.570796326794896619231321692-atan_val),p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            p_funcs_args_payload[i].arg5 = 2;
            const double result2 = dqng(p_integrand,0.0,atan_val,p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            const double full_result1{0.159154943091895335768883763*result1};
            const double full_result2{0.159154943091895335768883763*result2};
            p_functional[i] = full_result1+full_result2;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    return (0);
}

std::int32_t 
gms::fading_channel
::compute_functional_marcum_Q_4_10(quadpack_integrator_payload_t * __restrict__ p_payload)
{
    if(__builtin_expect(nullptr==p_payload,0)) { return (-1);}
    double (*p_integrand)(const double,void * __restrict__)  = p_payload->integrand;
    func_args_payload_t  * __restrict__ p_funcs_args_payload = p_payload->func_args_payload;
    double       * __restrict__ p_tmp_work1                  = p_payload->tmp_work1;
    double       * __restrict__ p_tmp_work2                  = p_payload->tmp_work2;
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
    double                      up_lim                       = p_payload->rand_lo3; // upper limit of Marcum-Q first order integration
    const std::int32_t          nfunc_vals                   = p_payload->n_func_vals;
    const std::int32_t          integrator_type              = p_payload->which_integrator;
    thread_local std::uniform_real_distribution<double> rv_func_s;
    thread_local std::mt19937 rv_func_s_gen;
    thread_local std::uint64_t seed_func_s{};
    thread_local std::uniform_real_distribution<double> rv_func_y;
    thread_local std::mt19937 rv_func_y_gen;
    thread_local std::uint64_t seed_func_y{};
    rv_func_s = std::uniform_real_distribution<double>(rand_low1,rand_high1);
    seed_func_s = __rdtsc();
    rv_func_s_gen = std::mt19937(seed_func_s);
    rv_func_y = std::uniform_real_distribution<double>(rand_low2,rand_high2);
    seed_func_y = __rdtsc();
    rv_func_y_gen = std::mt19937(seed_func_y);
    for(std::int32_t i{0}; i<nfunc_vals; ++i)
    {
        const double val_s{rv_func_s.operator()(rv_func_s_gen)};
        p_tmp_work1[i] = val_s;
        const double val_y{rv_func_y.operator()(rv_func_y_gen)};
        p_tmp_work2[i] = val_y;

    }
    std::sort(&p_tmp_work1[0],&p_tmp_work1[nfunc_vals-1],std::less<double>());
    std::sort(&p_tmp_work2[0],&p_tmp_work2[nfunc_vals-1],std::less<double>());
    if(integrator_type==1)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i)
        {
            const double cpy_s{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1 = cpy_s;
            const double cpy_y{p_tmp_work2[i]};
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqage(p_integrand,cpy_y,up_lim,p_epsabs[0],p_epsrel[0],p_irule[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_last[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = result;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    else if(integrator_type==2)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i)
        {
            const double cpy_s{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1 = cpy_s;
            const double cpy_y{p_tmp_work2[i]};
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqagi(p_integrand,cpy_y,p_inf[0],p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = result;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    else if(integrator_type==3)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i)
        {
            const double cpy_s{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1 = cpy_s;
            const double cpy_y{p_tmp_work2[i]};
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqags(p_integrand,cpy_y,up_lim,p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = result;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    else if(integrator_type==4)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i)
        {
            const double cpy_s{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1 = cpy_s;
            const double cpy_y{p_tmp_work2[i]};
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqng(p_integrand,cpy_y,up_lim,p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = result;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    return (0);
}

std::int32_t 
gms::fading_channel
::compute_functional_marcum_Q_4_16(quadpack_integrator_payload_t * __restrict__ p_payload)
{
    if(__builtin_expect(nullptr==p_payload,0)) { return (-1);}
    double (*p_integrand)(const double,void * __restrict__)  = p_payload->integrand;
    func_args_payload_t  * __restrict__ p_funcs_args_payload = p_payload->func_args_payload;
    double       * __restrict__ p_tmp_work1                  = p_payload->tmp_work1;
    double       * __restrict__ p_tmp_work2                  = p_payload->tmp_work2;
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
    thread_local std::uniform_real_distribution<double> rv_func_psi;
    thread_local std::mt19937 rv_func_psi_gen;
    thread_local std::uint64_t seed_func_psi{};
    thread_local std::uniform_real_distribution<double> rv_func_beta;
    thread_local std::mt19937 rv_func_beta_gen;
    thread_local std::uint64_t seed_func_beta{};
    rv_func_psi = std::uniform_real_distribution<double>(rand_low1,rand_high1);
    seed_func_psi = __rdtsc();
    rv_func_psi_gen = std::mt19937(seed_func_psi);
    rv_func_beta = std::uniform_real_distribution<double>(rand_low2,rand_high2);
    seed_func_beta = __rdtsc();
    rv_func_beta_gen = std::mt19937(seed_func_beta);
    for(std::int32_t i{0}; i<nfunc_vals; ++i)
    {
        const double val_psi{rv_func_psi.operator()(rv_func_psi_gen)};
        p_tmp_work1[i] = val_psi;
        const double val_beta{rv_func_beta.operator()(rv_func_beta_gen)};
        p_tmp_work2[i] = val_beta;

    }
    std::sort(&p_tmp_work1[0],&p_tmp_work1[nfunc_vals-1],std::less<double>());// psi
    std::sort(&p_tmp_work2[0],&p_tmp_work2[nfunc_vals-1],std::less<double>());// beta
    if(integrator_type==1)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i)
        {
            const double cpy_psi{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg2 = cpy_psi;
            const double cpy_beta{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg1 = cpy_beta;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqage(p_integrand,-3.141592653589793238462643383,+3.141592653589793238462643383,p_epsabs[0],p_epsrel[0],p_irule[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_last[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = 0.159154943091895335768883763*result;
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
            const double cpy_psi{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg2 = cpy_psi;
            const double cpy_beta{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg1 = cpy_beta;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqags(p_integrand,-3.141592653589793238462643383,+3.141592653589793238462643383,p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = 0.159154943091895335768883763*result;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    else if(integrator_type==4)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i) 
        {
            const double cpy_psi{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg2 = cpy_psi;
            const double cpy_beta{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg1 = cpy_beta;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqng(p_integrand,-3.141592653589793238462643383,+3.141592653589793238462643383,p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = 0.159154943091895335768883763*result;
            p_crude_tsc_meter[i] = end-start;
        } 
    }
    return (0);
}

std::int32_t 
gms::fading_channel
::compute_functional_marcum_Q_4_20_lo(quadpack_integrator_payload_t * __restrict__ p_payload)                           
{
    if(__builtin_expect(nullptr==p_payload,0)) { return (-1);}
    double (*p_integrand)(const double,void * __restrict__)  = p_payload->integrand;
    func_args_payload_t  * __restrict__ p_funcs_args_payload = p_payload->func_args_payload;
    double       * __restrict__ p_tmp_work1                  = p_payload->tmp_work1;
    double       * __restrict__ p_tmp_work2                  = p_payload->tmp_work2;
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
    thread_local std::uniform_real_distribution<double> rv_func_psi;
    thread_local std::mt19937 rv_func_psi_gen;
    thread_local std::uint64_t seed_func_psi{};
    thread_local std::uniform_real_distribution<double> rv_func_beta;
    thread_local std::mt19937 rv_func_beta_gen;
    thread_local std::uint64_t seed_func_beta{};
    rv_func_psi = std::uniform_real_distribution<double>(rand_low1,rand_high1);
    seed_func_psi = __rdtsc();
    rv_func_psi_gen = std::mt19937(seed_func_psi);
    rv_func_beta = std::uniform_real_distribution<double>(rand_low2,rand_high2);
    seed_func_beta = __rdtsc();
    rv_func_beta_gen = std::mt19937(seed_func_beta);
    for(std::int32_t i{0}; i<nfunc_vals; ++i)
    {
        const double val_psi{rv_func_psi.operator()(rv_func_psi_gen)};
        p_tmp_work1[i] = val_psi;
        const double val_beta{rv_func_beta.operator()(rv_func_beta_gen)};
        p_tmp_work2[i] = val_beta;

    }
    std::sort(&p_tmp_work1[0],&p_tmp_work1[nfunc_vals-1],std::less<double>());// psi
    std::sort(&p_tmp_work2[0],&p_tmp_work2[nfunc_vals-1],std::less<double>());// beta
    if(integrator_type==1)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i)
        {
            const double cpy_psi{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg2 = cpy_psi;
            const double cpy_beta{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg1 = cpy_beta;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqage(p_integrand,-3.141592653589793238462643383,+0.0,p_epsabs[0],p_epsrel[0],p_irule[0],&p_abser[i],
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
            const double cpy_psi{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg2 = cpy_psi;
            const double cpy_beta{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg1 = cpy_beta;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqags(p_integrand,-3.141592653589793238462643383,+0.0,p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = 0.318309886183790671537767527*result;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    else if(integrator_type==4)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i) 
        {
            const double cpy_psi{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg2 = cpy_psi;
            const double cpy_beta{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg1 = cpy_beta;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqng(p_integrand,-3.141592653589793238462643383,+0.0,p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = 0.318309886183790671537767527*result;
            p_crude_tsc_meter[i] = end-start;
        } 
    }
    return (0);
}

