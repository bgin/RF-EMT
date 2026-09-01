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
#include "GMS_compute_functionals_ch5.h"
#include "GMS_cquadpack.h"
#include "GMS_machine_utils.h" // for RDTSCP wrapper
#include "GMS_tabulated_quadrature.h" // currently serves for the outer functional quadrature.
/*
    const std::string integrators_names[6] = {{"dqage"},{"dqagi"},{"dqags"},{"dqng"},{"dqagp"},{"dqaws"}};
    double (*integrand)(double,void * __restrict__);
    func_args_payload_t * __restrict__ func_args_payload{nullptr};
    double              * __restrict__ tmp_work1{nullptr}; //work storage (caller provided) used mainly for the functional a1st rguments sorting.
    double              * __restrict__ tmp_work2{nullptr};//work storage (caller provided) used mainly for the functional 2nd arguments sorting.
    double              * __restrict__ tmp_work3{nullptr};//work storage (caller provided) used mainly for the functional 3rd arguments sorting.
    double              * __restrict__ tmp_work4{nullptr};//work storage (caller provided) used mainly for the functional 4th arguments sorting.
    double              * __restrict__ tmp_work5{nullptr};//work storage (caller provided) used mainly for the functional 5th arguments sorting.
    double              * __restrict__ tmp_work6{nullptr};//work storage (caller provided) used mainly for the functional 6th arguments sorting.
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
    std::uint64_t       * __restrict__ crude_tsc_measurement{nullptr}; // crude RDTSCP measurement for approximate TSC assesment (shall not be used for the robust statistics)
    double                             rand_lo1{}; // set lower limit for random number generator
    double                             rand_hi1{}; // as above (1st argument pair)
    double                             rand_lo2{}; // set lower limit for random number generator
    double                             rand_hi2{}; // as above (2nd argument pair)
    double                             rand_lo3{}; // set lower limit for random number generator
    double                             rand_hi3{}; // as above (3rd argument pair)
    double                             rand_lo4{}; // set lower limit for random number generator
    double                             rand_hi4{}; // as above (4th argument pair)
    double                             rand_lo5{}; // set lower limit for random number generator
    double                             rand_hi5{}; // as above (5th argument pair)
    double                             rand_lo6{}; // set lower limit for random number generator
    double                             rand_hi6{}; // as above (6th argument pair)
    double                             rand_lo7{}; // set lower limit for random number generator
    double                             rand_hi7{}; // as above (7th argument pair)
    double                             rand_lo8{}; // set lower limit for random number generator
    double                             rand_hi8{}; // as above (8th argument pair)
    std::int32_t                       n_func_vals;
    std::int32_t                       which_integrator; //currently 1=dqage,2=dqagi,3=dqags
    bool                               randomly_generate_inputs; //as the name states: random input generation in use if true, otherwise provide deterministic inputs
*/

std::int32_t 
gms::fading_channel
::compute_functional_Rayleigh_chan_5_6(quadpack_integrator_payload_ch5_t * __restrict__ p_payload)
{
    if(__builtin_expect(nullptr==p_payload,0)) { return (-1);}
    double (*p_integrand)(const double,void * __restrict__)  = p_payload->integrand;
    func_args_ch5_payload_t  * __restrict__ p_funcs_args_payload = p_payload->func_args_payload;
    double  * __restrict__      p_tmp_work1                  = p_payload->tmp_work1;
    double  * __restrict__      p_tmp_work2                  = p_payload->tmp_work2;
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
    std::uint64_t * __restrict__ p_crude_tsc_start           = p_payload->crude_tsc_start;
    std::uint64_t * __restrict__ p_crude_tsc_end             = p_payload->crude_tsc_end;
    std::uint64_t * __restrict__ p_crude_tsc_meter           = p_payload->crude_tsc_measurement;
    double                      rand_low1                    = p_payload->rand_lo1;
    double                      rand_high1                   = p_payload->rand_hi1;
    double                      rand_low2                    = p_payload->rand_lo2;
    double                      rand_high2                   = p_payload->rand_hi2;
    const std::int32_t          nfunc_vals                   = p_payload->n_func_vals;
    const std::int32_t          integrator_type              = p_payload->which_integrator;
    const bool                  random_input_generation      = p_payload->randomly_generate_inputs;
    const bool                  rand_in_gen_eq_true          = random_input_generation==true;
    if(rand_in_gen_eq_true)
    {
       thread_local std::uniform_real_distribution<double> rv_func_a;
       thread_local std::mt19937 rv_func_a_gen;
       thread_local std::uint64_t seed_func_a{};
       thread_local std::uniform_real_distribution<double> rv_func_gamma;
       thread_local std::mt19937 rv_func_gamma_gen;
       thread_local std::uint64_t seed_func_gamma{};
       rv_func_a = std::uniform_real_distribution<double>(rand_low1,rand_high1);
       seed_func_a = __rdtsc();
       rv_func_a_gen = std::mt19937(seed_func_a);
       rv_func_gamma = std::uniform_real_distribution<double>(rand_low2,rand_high2);
       seed_func_gamma = __rdtsc();
       rv_func_gamma_gen = std::mt19937(seed_func_gamma);
       for(std::int32_t i{0}; i<nfunc_vals; ++i) 
       {
           const double val_a{rv_func_a.operator()(rv_func_a_gen)};
           p_tmp_work1[i] = val_a;
           const double val_gamma{rv_func_gamma.operator()(rv_func_gamma_gen)};
           p_tmp_work2[i] = val_gamma;
       }
       std::sort(&p_tmp_work1[0],&p_tmp_work1[nfunc_vals-1],std::less<double>());
       std::sort(&p_tmp_work2[0],&p_tmp_work2[nfunc_vals-1],std::less<double>());
    }
    if(integrator_type==1)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i)
        {
            const double cpy_a{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1d = cpy_a;
            const double cpy_gamma{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg2d = cpy_gamma;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqage(p_integrand,0.0,1.570796326794896619231321692,p_epsabs[0],p_epsrel[0],p_irule[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_last[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = 0.318309886183790671537767527*result;
            p_crude_tsc_start[i] = start;
            p_crude_tsc_end[i]   = end;
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
            const double cpy_a{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1d = cpy_a;
            const double cpy_gamma{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg2d = cpy_gamma;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqags(p_integrand,0.0,1.570796326794896619231321692,p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = 0.318309886183790671537767527*result;
            p_crude_tsc_start[i] = start;
            p_crude_tsc_end[i]   = end;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    else if(integrator_type==4)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i)
        {
            const double cpy_a{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1d = cpy_a;
            const double cpy_gamma{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg2d = cpy_gamma;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqng(p_integrand,0.0,1.570796326794896619231321692,p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = 0.318309886183790671537767527*result;
            p_crude_tsc_start[i] = start;
            p_crude_tsc_end[i]   = end;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    return (0);
}

std::int32_t 
gms::fading_channel
::compute_functional_Hoyt_chan_5_9(quadpack_integrator_payload_ch5_t * __restrict__ p_payload)
{
    if(__builtin_expect(nullptr==p_payload,0)) { return (-1);}
    double (*p_integrand)(const double,void * __restrict__)  = p_payload->integrand;
    func_args_ch5_payload_t  * __restrict__ p_funcs_args_payload = p_payload->func_args_payload;
    double  * __restrict__      p_tmp_work1                  = p_payload->tmp_work1;
    double  * __restrict__      p_tmp_work2                  = p_payload->tmp_work2;
    double  * __restrict__      p_tmp_work3                  = p_payload->tmp_work3; 
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
    std::uint64_t * __restrict__ p_crude_tsc_start           = p_payload->crude_tsc_start;
    std::uint64_t * __restrict__ p_crude_tsc_end             = p_payload->crude_tsc_end;
    std::uint64_t * __restrict__ p_crude_tsc_meter           = p_payload->crude_tsc_measurement;
    double                      rand_low1                    = p_payload->rand_lo1;
    double                      rand_high1                   = p_payload->rand_hi1;
    double                      rand_low2                    = p_payload->rand_lo2;
    double                      rand_high2                   = p_payload->rand_hi2;
    double                      rand_low3                    = p_payload->rand_lo3;
    double                      rand_high3                   = p_payload->rand_hi3;
    const std::int32_t          nfunc_vals                   = p_payload->n_func_vals;
    const std::int32_t          integrator_type              = p_payload->which_integrator;
    const bool                  random_input_generation      = p_payload->randomly_generate_inputs;
    const bool                  rand_in_gen_eq_true          = random_input_generation==true;
    if(rand_in_gen_eq_true)
    {
       thread_local std::uniform_real_distribution<double> rv_func_a;
       thread_local std::mt19937 rv_func_a_gen;
       thread_local std::uint64_t seed_func_a{};
       thread_local std::uniform_real_distribution<double> rv_func_gamma;
       thread_local std::mt19937 rv_func_gamma_gen;
       thread_local std::uint64_t seed_func_gamma{};
       thread_local std::uniform_real_distribution<double> rv_func_q;
       thread_local std::mt19937 rv_func_q_gen;
       thread_local std::uint64_t seed_func_q{};
       rv_func_a = std::uniform_real_distribution<double>(rand_low1,rand_high1);
       seed_func_a = __rdtsc();
       rv_func_a_gen = std::mt19937(seed_func_a);
       rv_func_gamma = std::uniform_real_distribution<double>(rand_low2,rand_high2);
       seed_func_gamma = __rdtsc();
       rv_func_gamma_gen = std::mt19937(seed_func_gamma);
       rv_func_q = std::uniform_real_distribution<double>(rand_low3,rand_high3);
       seed_func_q = __rdtsc();
       rv_func_q_gen = std::mt19937(seed_func_q);
       for(std::int32_t i{0}; i<nfunc_vals; ++i) 
       {
           const double val_a{rv_func_a.operator()(rv_func_a_gen)};
           p_tmp_work1[i] = val_a;
           const double val_gamma{rv_func_gamma.operator()(rv_func_gamma_gen)};
           p_tmp_work2[i] = val_gamma;
           const double val_q{rv_func_q.operator()(rv_func_q_gen)};
           p_tmp_work3[i] = val_q;
       }
       std::sort(&p_tmp_work1[0],&p_tmp_work1[nfunc_vals-1],std::less<double>());
       std::sort(&p_tmp_work2[0],&p_tmp_work2[nfunc_vals-1],std::less<double>());
       std::sort(&p_tmp_work3[0],&p_tmp_work3[nfunc_vals-1],std::less<double>());
    }
    if(integrator_type==1)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i)
        {
            const double cpy_a{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1d = cpy_a;
            const double cpy_gamma{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg2d = cpy_gamma;
            const double cpy_q{p_tmp_work3[i]};
            p_funcs_args_payload[i].arg3d = cpy_q;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqage(p_integrand,0.0,1.570796326794896619231321692,p_epsabs[0],p_epsrel[0],p_irule[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_last[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = 0.318309886183790671537767527*result;
            p_crude_tsc_start[i] = start;
            p_crude_tsc_end[i]   = end;
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
            const double cpy_a{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1d = cpy_a;
            const double cpy_gamma{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg2d = cpy_gamma;
            const double cpy_q{p_tmp_work3[i]};
            p_funcs_args_payload[i].arg3d = cpy_q;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqags(p_integrand,0.0,1.570796326794896619231321692,p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = 0.318309886183790671537767527*result;
            p_crude_tsc_start[i] = start;
            p_crude_tsc_end[i]   = end;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    else if(integrator_type==4)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i)
        {
            const double cpy_a{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1d = cpy_a;
            const double cpy_gamma{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg2d = cpy_gamma;
            const double cpy_q{p_tmp_work3[i]};
            p_funcs_args_payload[i].arg3d = cpy_q;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqng(p_integrand,0.0,1.570796326794896619231321692,p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = 0.318309886183790671537767527*result;
            p_crude_tsc_start[i] = start;
            p_crude_tsc_end[i]   = end;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    return (0);
}

std::int32_t 
gms::fading_channel
::compute_functional_Rice_chan_5_12(quadpack_integrator_payload_ch5_t * __restrict__ p_payload)
{
    if(__builtin_expect(nullptr==p_payload,0)) { return (-1);}
    double (*p_integrand)(const double,void * __restrict__)  = p_payload->integrand;
    func_args_ch5_payload_t  * __restrict__ p_funcs_args_payload = p_payload->func_args_payload;
    double  * __restrict__      p_tmp_work1                  = p_payload->tmp_work1;
    double  * __restrict__      p_tmp_work2                  = p_payload->tmp_work2;
    double  * __restrict__      p_tmp_work3                  = p_payload->tmp_work3; 
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
    std::uint64_t * __restrict__ p_crude_tsc_start           = p_payload->crude_tsc_start;
    std::uint64_t * __restrict__ p_crude_tsc_end             = p_payload->crude_tsc_end;
    std::uint64_t * __restrict__ p_crude_tsc_meter           = p_payload->crude_tsc_measurement;
    double                      rand_low1                    = p_payload->rand_lo1;
    double                      rand_high1                   = p_payload->rand_hi1;
    double                      rand_low2                    = p_payload->rand_lo2;
    double                      rand_high2                   = p_payload->rand_hi2;
    double                      rand_low3                    = p_payload->rand_lo3;
    double                      rand_high3                   = p_payload->rand_hi3;
    const std::int32_t          nfunc_vals                   = p_payload->n_func_vals;
    const std::int32_t          integrator_type              = p_payload->which_integrator;
    const bool                  random_input_generation      = p_payload->randomly_generate_inputs;
    const bool                  rand_in_gen_eq_true          = random_input_generation==true;
    if(rand_in_gen_eq_true)
    {
       thread_local std::uniform_real_distribution<double> rv_func_a;
       thread_local std::mt19937 rv_func_a_gen;
       thread_local std::uint64_t seed_func_a{};
       thread_local std::uniform_real_distribution<double> rv_func_gamma;
       thread_local std::mt19937 rv_func_gamma_gen;
       thread_local std::uint64_t seed_func_gamma{};
       thread_local std::uniform_real_distribution<double> rv_func_n;
       thread_local std::mt19937 rv_func_n_gen;
       thread_local std::uint64_t seed_func_n{};
       rv_func_a = std::uniform_real_distribution<double>(rand_low1,rand_high1);
       seed_func_a = __rdtsc();
       rv_func_a_gen = std::mt19937(seed_func_a);
       rv_func_gamma = std::uniform_real_distribution<double>(rand_low2,rand_high2);
       seed_func_gamma = __rdtsc();
       rv_func_gamma_gen = std::mt19937(seed_func_gamma);
       rv_func_n = std::uniform_real_distribution<double>(rand_low3,rand_high3);
       seed_func_n = __rdtsc();
       rv_func_n_gen = std::mt19937(seed_func_n);
       for(std::int32_t i{0}; i<nfunc_vals; ++i) 
       {
           const double val_a{rv_func_a.operator()(rv_func_a_gen)};
           p_tmp_work1[i] = val_a;
           const double val_gamma{rv_func_gamma.operator()(rv_func_gamma_gen)};
           p_tmp_work2[i] = val_gamma;
           const double val_n{rv_func_n.operator()(rv_func_n_gen)};
           p_tmp_work3[i] = val_n;
       }
       std::sort(&p_tmp_work1[0],&p_tmp_work1[nfunc_vals-1],std::less<double>());
       std::sort(&p_tmp_work2[0],&p_tmp_work2[nfunc_vals-1],std::less<double>());
       std::sort(&p_tmp_work3[0],&p_tmp_work3[nfunc_vals-1],std::less<double>());
    }
    if(integrator_type==1)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i)
        {
            const double cpy_a{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1d = cpy_a;
            const double cpy_gamma{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg3d = cpy_gamma;
            const double cpy_n{p_tmp_work3[i]};
            p_funcs_args_payload[i].arg2d = cpy_n;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqage(p_integrand,0.0,1.570796326794896619231321692,p_epsabs[0],p_epsrel[0],p_irule[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_last[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = 0.318309886183790671537767527*result;
            p_crude_tsc_start[i] = start;
            p_crude_tsc_end[i]   = end;
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
            const double cpy_a{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1d = cpy_a;
            const double cpy_gamma{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg3d = cpy_gamma;
            const double cpy_n{p_tmp_work3[i]};
            p_funcs_args_payload[i].arg2d = cpy_n;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqags(p_integrand,0.0,1.570796326794896619231321692,p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = 0.318309886183790671537767527*result;
            p_crude_tsc_start[i] = start;
            p_crude_tsc_end[i]   = end;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    else if(integrator_type==4)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i)
        {
            const double cpy_a{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1d = cpy_a;
            const double cpy_gamma{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg3d = cpy_gamma;
            const double cpy_n{p_tmp_work3[i]};
            p_funcs_args_payload[i].arg2d = cpy_n;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqng(p_integrand,0.0,1.570796326794896619231321692,p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = 0.318309886183790671537767527*result;
            p_crude_tsc_start[i] = start;
            p_crude_tsc_end[i]   = end;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    return (0);
}

std::int32_t 
gms::fading_channel
::compute_functional_Nakagami_m_chan_5_16(quadpack_integrator_payload_ch5_t * __restrict__ p_payload)
{
    if(__builtin_expect(nullptr==p_payload,0)) { return (-1);}
    double (*p_integrand)(const double,void * __restrict__)  = p_payload->integrand;
    func_args_ch5_payload_t  * __restrict__ p_funcs_args_payload = p_payload->func_args_payload;
    double  * __restrict__      p_tmp_work1                  = p_payload->tmp_work1;
    double  * __restrict__      p_tmp_work2                  = p_payload->tmp_work2;
    double  * __restrict__      p_tmp_work3                  = p_payload->tmp_work3; 
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
    std::uint64_t * __restrict__ p_crude_tsc_start           = p_payload->crude_tsc_start;
    std::uint64_t * __restrict__ p_crude_tsc_end             = p_payload->crude_tsc_end;
    std::uint64_t * __restrict__ p_crude_tsc_meter           = p_payload->crude_tsc_measurement;
    double                      rand_low1                    = p_payload->rand_lo1;
    double                      rand_high1                   = p_payload->rand_hi1;
    double                      rand_low2                    = p_payload->rand_lo2;
    double                      rand_high2                   = p_payload->rand_hi2;
    double                      rand_low3                    = p_payload->rand_lo3;
    double                      rand_high3                   = p_payload->rand_hi3;
    const std::int32_t          nfunc_vals                   = p_payload->n_func_vals;
    const std::int32_t          integrator_type              = p_payload->which_integrator;
    const bool                  random_input_generation      = p_payload->randomly_generate_inputs;
    const bool                  rand_in_gen_eq_true          = random_input_generation==true;
    if(rand_in_gen_eq_true)
    {
       thread_local std::uniform_real_distribution<double> rv_func_a;
       thread_local std::mt19937 rv_func_a_gen;
       thread_local std::uint64_t seed_func_a{};
       thread_local std::uniform_real_distribution<double> rv_func_gamma;
       thread_local std::mt19937 rv_func_gamma_gen;
       thread_local std::uint64_t seed_func_gamma{};
       thread_local std::uniform_real_distribution<double> rv_func_m;
       thread_local std::mt19937 rv_func_m_gen;
       thread_local std::uint64_t seed_func_m{};
       rv_func_a = std::uniform_real_distribution<double>(rand_low1,rand_high1);
       seed_func_a = __rdtsc();
       rv_func_a_gen = std::mt19937(seed_func_a);
       rv_func_gamma = std::uniform_real_distribution<double>(rand_low2,rand_high2);
       seed_func_gamma = __rdtsc();
       rv_func_gamma_gen = std::mt19937(seed_func_gamma);
       rv_func_m = std::uniform_real_distribution<double>(rand_low3,rand_high3);
       seed_func_m = __rdtsc();
       rv_func_m_gen = std::mt19937(seed_func_m);
       for(std::int32_t i{0}; i<nfunc_vals; ++i) 
       {
           const double val_a{rv_func_a.operator()(rv_func_a_gen)};
           p_tmp_work1[i] = val_a;
           const double val_gamma{rv_func_gamma.operator()(rv_func_gamma_gen)};
           p_tmp_work2[i] = val_gamma;
           const double val_m{rv_func_m.operator()(rv_func_m_gen)};
           p_tmp_work3[i] = val_m;
       }
       std::sort(&p_tmp_work1[0],&p_tmp_work1[nfunc_vals-1],std::less<double>());
       std::sort(&p_tmp_work2[0],&p_tmp_work2[nfunc_vals-1],std::less<double>());
       std::sort(&p_tmp_work3[0],&p_tmp_work3[nfunc_vals-1],std::less<double>());
    }
    if(integrator_type==1)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i)
        {
            const double cpy_a{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1d = cpy_a;
            const double cpy_gamma{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg3d = cpy_gamma;
            const double cpy_m{p_tmp_work3[i]};
            p_funcs_args_payload[i].arg2d = cpy_m;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqage(p_integrand,0.0,1.570796326794896619231321692,p_epsabs[0],p_epsrel[0],p_irule[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_last[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = 0.318309886183790671537767527*result;
            p_crude_tsc_start[i] = start;
            p_crude_tsc_end[i]   = end;
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
            const double cpy_a{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1d = cpy_a;
            const double cpy_gamma{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg3d = cpy_gamma;
            const double cpy_m{p_tmp_work3[i]};
            p_funcs_args_payload[i].arg2d = cpy_m;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqags(p_integrand,0.0,1.570796326794896619231321692,p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = 0.318309886183790671537767527*result;
            p_crude_tsc_start[i] = start;
            p_crude_tsc_end[i]   = end;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    else if(integrator_type==4)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i)
        {
            const double cpy_a{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1d = cpy_a;
            const double cpy_gamma{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg3d = cpy_gamma;
            const double cpy_m{p_tmp_work3[i]};
            p_funcs_args_payload[i].arg2d = cpy_m;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqng(p_integrand,0.0,1.570796326794896619231321692,p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = 0.318309886183790671537767527*result;
            p_crude_tsc_start[i] = start;
            p_crude_tsc_end[i]   = end;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    return (0);
}

std::int32_t 
gms::fading_channel
::compute_functional_LogNormShadow_chan_5_20(quadpack_integrator_payload_ch5_t * __restrict__ p_payload)
{
    if(__builtin_expect(nullptr==p_payload,0)) { return (-1);}
    double (*p_integrand)(const double,void * __restrict__)  = p_payload->integrand;
    func_args_ch5_payload_t  * __restrict__ p_funcs_args_payload = p_payload->func_args_payload;
    double  * __restrict__      p_tmp_work1                  = p_payload->tmp_work1;
    double  * __restrict__      p_tmp_work2                  = p_payload->tmp_work2;
    double  * __restrict__      p_tmp_work3                  = p_payload->tmp_work3; 
    double  * __restrict__      p_tmp_work4                  = p_payload->tmp_work4;
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
    std::uint64_t * __restrict__ p_crude_tsc_start           = p_payload->crude_tsc_start;
    std::uint64_t * __restrict__ p_crude_tsc_end             = p_payload->crude_tsc_end;
    std::uint64_t * __restrict__ p_crude_tsc_meter           = p_payload->crude_tsc_measurement;
    double                      rand_low1                    = p_payload->rand_lo1;
    double                      rand_high1                   = p_payload->rand_hi1;
    double                      rand_low2                    = p_payload->rand_lo2;
    double                      rand_high2                   = p_payload->rand_hi2;
    double                      rand_low3                    = p_payload->rand_lo3;
    double                      rand_high3                   = p_payload->rand_hi3;
    double                      rand_low4                    = p_payload->rand_lo4;
    double                      rand_high4                   = p_payload->rand_hi4;
    double                      rand_low5                    = p_payload->rand_lo5;
    double                      rand_high5                   = p_payload->rand_hi5;
    const std::int32_t          nfunc_vals                   = p_payload->n_func_vals;
    const std::int32_t          integrator_type              = p_payload->which_integrator;
    const bool                  random_input_generation      = p_payload->randomly_generate_inputs;
    const bool                  rand_in_gen_eq_true          = random_input_generation==true;
    if(rand_in_gen_eq_true)
    {
       thread_local std::uniform_real_distribution<double> rv_func_theta;
       thread_local std::mt19937 rv_func_theta_gen;
       thread_local std::uint64_t seed_func_theta{};
       thread_local std::uniform_real_distribution<double> rv_func_a;
       thread_local std::mt19937 rv_func_a_gen;
       thread_local std::uint64_t seed_func_a{};
       thread_local std::uniform_real_distribution<double> rv_func_mu;
       thread_local std::mt19937 rv_func_mu_gen;
       thread_local std::uint64_t seed_func_mu{};
       thread_local std::uniform_real_distribution<double> rv_func_sigma;
       thread_local std::mt19937 rv_func_sigma_gen;
       thread_local std::uint64_t seed_func_sigma{};
       rv_func_theta = std::uniform_real_distribution<double>(rand_low1,rand_high1);
       seed_func_theta = __rdtsc();
       rv_func_theta_gen = std::mt19937(seed_func_theta);
       rv_func_a = std::uniform_real_distribution<double>(rand_low2,rand_high2);
       seed_func_a = __rdtsc();
       rv_func_a_gen = std::mt19937(seed_func_a);
       rv_func_mu = std::uniform_real_distribution<double>(rand_low3,rand_high3);
       seed_func_mu = __rdtsc();
       rv_func_mu_gen = std::mt19937(seed_func_mu);
       rv_func_sigma = std::uniform_real_distribution<double>(rand_low4,rand_high4);
       seed_func_sigma = __rdtsc();
       rv_func_sigma_gen = std::mt19937(seed_func_sigma);

       for(std::int32_t i{0}; i<nfunc_vals; ++i) 
       {
           const double val_theta{rv_func_theta.operator()(rv_func_theta_gen)};
           p_tmp_work1[i] = val_theta;
           const double val_a{rv_func_a.operator()(rv_func_a_gen)};
           p_tmp_work2[i] = val_a;
           const double val_mu{rv_func_mu.operator()(rv_func_mu_gen)};
           p_tmp_work3[i] = val_mu;
           const double val_sigma{rv_func_sigma.operator()(rv_func_sigma_gen)};
           p_tmp_work4[i] = val_sigma;
       }
       std::sort(&p_tmp_work1[0],&p_tmp_work1[nfunc_vals-1],std::less<double>());
       std::sort(&p_tmp_work2[0],&p_tmp_work2[nfunc_vals-1],std::less<double>());
       std::sort(&p_tmp_work3[0],&p_tmp_work3[nfunc_vals-1],std::less<double>());
       std::sort(&p_tmp_work4[0],&p_tmp_work4[nfunc_vals-1],std::less<double>());
    }
    if(integrator_type==1)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i)
        {
            const double cpy_theta{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1d = cpy_theta;
            const double cpy_a{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg2d = cpy_a;
            const double cpy_mu{p_tmp_work3[i]};
            p_funcs_args_payload[i].arg3d = cpy_mu;
            const double cpy_sigma{p_tmp_work4[i]};
            p_funcs_args_payload[i].arg4d = cpy_sigma;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqage(p_integrand,rand_low5,rand_high5,p_epsabs[0],p_epsrel[0],p_irule[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_last[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = 0.564189583547756286948079452*result;
            p_crude_tsc_start[i] = start;
            p_crude_tsc_end[i]   = end;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    else if(integrator_type==2)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i)
        {
            const double cpy_theta{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1d = cpy_theta;
            const double cpy_a{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg2d = cpy_a;
            const double cpy_mu{p_tmp_work3[i]};
            p_funcs_args_payload[i].arg3d = cpy_mu;
            const double cpy_sigma{p_tmp_work4[i]};
            p_funcs_args_payload[i].arg4d = cpy_sigma;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqagi(p_integrand,rand_low5,p_inf[0],p_epsabs[0],p_epsrel[0],
                                        &p_abser[i],&p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = 0.564189583547756286948079452*result;
            p_crude_tsc_start[i] = start;
            p_crude_tsc_end[i]   = end;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    else if(integrator_type==3)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i)
        {
            const double cpy_theta{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1d = cpy_theta;
            const double cpy_a{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg2d = cpy_a;
            const double cpy_mu{p_tmp_work3[i]};
            p_funcs_args_payload[i].arg3d = cpy_mu;
            const double cpy_sigma{p_tmp_work4[i]};
            p_funcs_args_payload[i].arg4d = cpy_sigma;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqags(p_integrand,rand_low5,rand_high5,p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = 0.564189583547756286948079452*result;
            p_crude_tsc_start[i] = start;
            p_crude_tsc_end[i]   = end;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    else if(integrator_type==4)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i)
        {
            const double cpy_theta{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1d = cpy_theta;
            const double cpy_a{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg2d = cpy_a;
            const double cpy_mu{p_tmp_work3[i]};
            p_funcs_args_payload[i].arg3d = cpy_mu;
            const double cpy_sigma{p_tmp_work4[i]};
            p_funcs_args_payload[i].arg4d = cpy_sigma;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqng(p_integrand,rand_low5,rand_high5,p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = 0.564189583547756286948079452*result;
            p_crude_tsc_start[i] = start;
            p_crude_tsc_end[i]   = end;
            p_crude_tsc_meter[i] = end-start;
        }
    }
    return (0);
}

std::int32_t 
gms::fading_channel
::compute_outer_functional_LogNormShadow_chan_5_20(quadpack_integrator_payload_ch5_t * __restrict__ p_payload)
{
    if(__builtin_expect(nullptr==p_payload,0)) { return (-1);}
    double (*p_integrand)(const double,void * __restrict__)  = p_payload->integrand;
    func_args_ch5_payload_t  * __restrict__ p_funcs_args_payload = p_payload->func_args_payload;
    double  * __restrict__      p_tmp_work1                  = p_payload->tmp_work1;
    double  * __restrict__      p_tmp_work2                  = p_payload->tmp_work2;
    double  * __restrict__      p_tmp_work3                  = p_payload->tmp_work3; 
    double  * __restrict__      p_tmp_work4                  = p_payload->tmp_work4;
    double  * __restrict__      p_tmp_hiordq_work            = p_payload->tmp_work5;
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
    std::uint64_t * __restrict__ p_crude_tsc_start           = p_payload->crude_tsc_start;
    std::uint64_t * __restrict__ p_crude_tsc_end             = p_payload->crude_tsc_end;
    std::uint64_t * __restrict__ p_crude_tsc_meter           = p_payload->crude_tsc_measurement;
    std::uint64_t * __restrict__ p_crude_tsc_meas_outer      = p_payload->crude_tsc_meas_outer;
    double                      rand_low1                    = p_payload->rand_lo1;
    double                      rand_high1                   = p_payload->rand_hi1;
    double                      rand_low2                    = p_payload->rand_lo2;
    double                      rand_high2                   = p_payload->rand_hi2;
    double                      rand_low3                    = p_payload->rand_lo3;
    double                      rand_high3                   = p_payload->rand_hi3;
    double                      rand_low4                    = p_payload->rand_lo4;
    double                      rand_high4                   = p_payload->rand_hi4;
    double                      rand_low5                    = p_payload->rand_lo5;
    double                      rand_high5                   = p_payload->rand_hi5;
    const std::int32_t          nfunc_vals                   = p_payload->n_func_vals;
    const std::int32_t          integrator_type              = p_payload->which_integrator;
    const std::int32_t          tabular_integrator_type      = p_payload->which_tabulated_integrator;
    const bool                  random_input_generation      = p_payload->randomly_generate_inputs;
    const bool                  rand_in_gen_eq_true          = random_input_generation==true;
    if(rand_in_gen_eq_true)
    {
       thread_local std::uniform_real_distribution<double> rv_func_theta;
       thread_local std::mt19937 rv_func_theta_gen;
       thread_local std::uint64_t seed_func_theta{};
       thread_local std::uniform_real_distribution<double> rv_func_a;
       thread_local std::mt19937 rv_func_a_gen;
       thread_local std::uint64_t seed_func_a{};
       thread_local std::uniform_real_distribution<double> rv_func_mu;
       thread_local std::mt19937 rv_func_mu_gen;
       thread_local std::uint64_t seed_func_mu{};
       thread_local std::uniform_real_distribution<double> rv_func_sigma;
       thread_local std::mt19937 rv_func_sigma_gen;
       thread_local std::uint64_t seed_func_sigma{};
       rv_func_theta = std::uniform_real_distribution<double>(rand_low1,rand_high1);
       seed_func_theta = __rdtsc();
       rv_func_theta_gen = std::mt19937(seed_func_theta);
       rv_func_a = std::uniform_real_distribution<double>(rand_low2,rand_high2);
       seed_func_a = __rdtsc();
       rv_func_a_gen = std::mt19937(seed_func_a);
       rv_func_mu = std::uniform_real_distribution<double>(rand_low3,rand_high3);
       seed_func_mu = __rdtsc();
       rv_func_mu_gen = std::mt19937(seed_func_mu);
       rv_func_sigma = std::uniform_real_distribution<double>(rand_low4,rand_high4);
       seed_func_sigma = __rdtsc();
       rv_func_sigma_gen = std::mt19937(seed_func_sigma);

       for(std::int32_t i{0}; i<nfunc_vals; ++i) 
       {
           const double val_theta{rv_func_theta.operator()(rv_func_theta_gen)};
           p_tmp_work1[i] = val_theta;
           const double val_a{rv_func_a.operator()(rv_func_a_gen)};
           p_tmp_work2[i] = val_a;
           const double val_mu{rv_func_mu.operator()(rv_func_mu_gen)};
           p_tmp_work3[i] = val_mu;
           const double val_sigma{rv_func_sigma.operator()(rv_func_sigma_gen)};
           p_tmp_work4[i] = val_sigma;
       }
       std::sort(&p_tmp_work1[0],&p_tmp_work1[nfunc_vals-1],std::less<double>());
       std::sort(&p_tmp_work2[0],&p_tmp_work2[nfunc_vals-1],std::less<double>());
       std::sort(&p_tmp_work3[0],&p_tmp_work3[nfunc_vals-1],std::less<double>());
       std::sort(&p_tmp_work4[0],&p_tmp_work4[nfunc_vals-1],std::less<double>());
    }
    if(integrator_type==1)
    {
        for(std::int32_t i{0}; i<nfunc_vals; ++i)
        {
            const double cpy_theta{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1d = cpy_theta;
            const double cpy_a{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg2d = cpy_a;
            const double cpy_mu{p_tmp_work3[i]};
            p_funcs_args_payload[i].arg3d = cpy_mu;
            const double cpy_sigma{p_tmp_work4[i]};
            p_funcs_args_payload[i].arg4d = cpy_sigma;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqage(p_integrand,rand_low5,rand_high5,p_epsabs[0],p_epsrel[0],p_irule[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_last[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = 0.564189583547756286948079452*result;
            p_crude_tsc_start[i] = start;
            p_crude_tsc_end[i]   = end;
            p_crude_tsc_meter[i] = end-start;
        }
        const double abscissa_step = 1.0/static_cast<double>(nfunc_vals);
        static std::int32_t call_counter = 0;
#ifdef _OPENMP
#pragma threadprivate(call_counter)
#endif 
        if(tabular_integrator_type==1) 
        { 
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            (void)math::hiordq(nfunc_vals,abscissa_step,&p_functional[0],&p_tmp_hiordq_work[0],p_payload->outer_func_tmp_res); 
            p_payload->outer_func_tmp_res *= 0.318309886183790671537767527;
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            call_counter += 1;
            p_crude_tsc_meas_outer[call_counter] = end-start;
        }
        else if(tabular_integrator_type==2) 
        {   
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            (void)math::simpn(nfunc_vals,abscissa_step,&p_functional[0],p_payload->outer_func_tmp_res); 
            p_payload->outer_func_tmp_res *= 0.318309886183790671537767527;
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            call_counter += 1;
            p_crude_tsc_meas_outer[call_counter] = end-start;
        }
        else if(tabular_integrator_type==3) 
        {  
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            (void)math::wedint(nfunc_vals,abscissa_step,&p_functional[0],p_payload->outer_func_tmp_res);
            p_payload->outer_func_tmp_res *= 0.318309886183790671537767527;
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            call_counter += 1;
            p_crude_tsc_meas_outer[call_counter] = end-start;
        }
    }
    else if(integrator_type==2)
    {
         static std::int32_t call_counter = 0;
#ifdef _OPENMP
#pragma threadprivate(call_counter)
#endif 
        for(std::int32_t i{0}; i<nfunc_vals; ++i)
        {
            const double cpy_theta{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1d = cpy_theta;
            const double cpy_a{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg2d = cpy_a;
            const double cpy_mu{p_tmp_work3[i]};
            p_funcs_args_payload[i].arg3d = cpy_mu;
            const double cpy_sigma{p_tmp_work4[i]};
            p_funcs_args_payload[i].arg4d = cpy_sigma;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqagi(p_integrand,rand_low5,p_inf[0],p_epsabs[0],p_epsrel[0],
                                        &p_abser[i],&p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = 0.564189583547756286948079452*result;
            p_crude_tsc_start[i] = start;
            p_crude_tsc_end[i]   = end;
            p_crude_tsc_meter[i] = end-start;
        }
        const double abscissa_step = 1.0/static_cast<double>(nfunc_vals);
        if(tabular_integrator_type==1) 
        { 
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            (void)math::hiordq(nfunc_vals,abscissa_step,&p_functional[0],&p_tmp_hiordq_work[0],p_payload->outer_func_tmp_res); 
            p_payload->outer_func_tmp_res *= 0.318309886183790671537767527;
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            call_counter += 1;
            p_crude_tsc_meas_outer[call_counter] = end-start;
        }
        else if(tabular_integrator_type==2) 
        { 
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            (void)math::simpn(nfunc_vals,abscissa_step,&p_functional[0],p_payload->outer_func_tmp_res); 
            p_payload->outer_func_tmp_res *= 0.318309886183790671537767527;
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            call_counter += 1;
            p_crude_tsc_meas_outer[call_counter] = end-start;
        }
        else if(tabular_integrator_type==3) 
        {   
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            (void)math::wedint(nfunc_vals,abscissa_step,&p_functional[0],p_payload->outer_func_tmp_res);
            p_payload->outer_func_tmp_res *= 0.318309886183790671537767527;
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            call_counter += 1;
            p_crude_tsc_meas_outer[call_counter] = end-start;
        }
    }
    else if(integrator_type==3)
    {
        static std::int32_t call_counter = 0;
#ifdef _OPENMP
#pragma threadprivate(call_counter)
#endif 
        for(std::int32_t i{0}; i<nfunc_vals; ++i)
        {
            const double cpy_theta{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1d = cpy_theta;
            const double cpy_a{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg2d = cpy_a;
            const double cpy_mu{p_tmp_work3[i]};
            p_funcs_args_payload[i].arg3d = cpy_mu;
            const double cpy_sigma{p_tmp_work4[i]};
            p_funcs_args_payload[i].arg4d = cpy_sigma;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqags(p_integrand,rand_low5,rand_high5,p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = 0.564189583547756286948079452*result;
            p_crude_tsc_start[i] = start;
            p_crude_tsc_end[i]   = end;
            p_crude_tsc_meter[i] = end-start;
        }
        const double abscissa_step = 1.0/static_cast<double>(nfunc_vals);
        if(tabular_integrator_type==1) 
        {   
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            (void)math::hiordq(nfunc_vals,abscissa_step,&p_functional[0],&p_tmp_hiordq_work[0],p_payload->outer_func_tmp_res); 
            p_payload->outer_func_tmp_res *= 0.318309886183790671537767527;
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            call_counter += 1;
            p_crude_tsc_meas_outer[call_counter] = end-start;
        }
        else if(tabular_integrator_type==2) 
        { 
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            (void)math::simpn(nfunc_vals,abscissa_step,&p_functional[0],p_payload->outer_func_tmp_res); 
            p_payload->outer_func_tmp_res *= 0.318309886183790671537767527;
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            call_counter += 1;
            p_crude_tsc_meas_outer[call_counter] = end-start;
        }
        else if(tabular_integrator_type==3) 
        {   
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            (void)math::wedint(nfunc_vals,abscissa_step,&p_functional[0],p_payload->outer_func_tmp_res);
            p_payload->outer_func_tmp_res *= 0.318309886183790671537767527;
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            call_counter += 1;
            p_crude_tsc_meas_outer[call_counter] = end-start;
        }
    }
    else if(integrator_type==4)
    {
    static std::int32_t call_counter = 0;
#ifdef _OPENMP
#pragma threadprivate(call_counter)
#endif 
        for(std::int32_t i{0}; i<nfunc_vals; ++i)
        {
            const double cpy_theta{p_tmp_work1[i]};
            p_funcs_args_payload[i].arg1d = cpy_theta;
            const double cpy_a{p_tmp_work2[i]};
            p_funcs_args_payload[i].arg2d = cpy_a;
            const double cpy_mu{p_tmp_work3[i]};
            p_funcs_args_payload[i].arg3d = cpy_mu;
            const double cpy_sigma{p_tmp_work4[i]};
            p_funcs_args_payload[i].arg4d = cpy_sigma;
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            const double result = dqng(p_integrand,rand_low5,rand_high5,p_epsabs[0],p_epsrel[0],&p_abser[i],
                                        &p_neval[i],&p_ier[i],&p_funcs_args_payload[i]);
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            p_functional[i] = 0.564189583547756286948079452*result;
            p_crude_tsc_start[i] = start;
            p_crude_tsc_end[i]   = end;
            p_crude_tsc_meter[i] = end-start;
        }
        const double abscissa_step = 1.0/static_cast<double>(nfunc_vals);
        if(tabular_integrator_type==1) 
        { 
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            (void)math::hiordq(nfunc_vals,abscissa_step,&p_functional[0],&p_tmp_hiordq_work[0],p_payload->outer_func_tmp_res); 
            p_payload->outer_func_tmp_res *= 0.318309886183790671537767527;
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            call_counter += 1;
            p_crude_tsc_meas_outer[call_counter] = end-start;
        }
        else if(tabular_integrator_type==2) 
        {   
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            (void)math::simpn(nfunc_vals,abscissa_step,&p_functional[0],p_payload->outer_func_tmp_res); 
            p_payload->outer_func_tmp_res *= 0.318309886183790671537767527;
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            call_counter += 1;
            p_crude_tsc_meas_outer[call_counter] = end-start;
        }
        else if(tabular_integrator_type==3) 
        {   
            const std::uint64_t start{gms::common::rdtsc_serialized_start()};
            (void)math::wedint(nfunc_vals,abscissa_step,&p_functional[0],p_payload->outer_func_tmp_res);
            p_payload->outer_func_tmp_res *= 0.318309886183790671537767527;
            const std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            call_counter += 1;
            p_crude_tsc_meas_outer[call_counter] = end-start;
        }
    }
    return (0);
}