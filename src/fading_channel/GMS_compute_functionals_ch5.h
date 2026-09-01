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

#ifndef __GMS_COMPUTE_FUNCTIONALS_CH5_H__
#define __GMS_COMPUTE_FUNCTIONALS_CH5_H__ 280820260602

#include <cstdint>
#include <string>
#include <valarray>
#include <array>
#include "GMS_config.h"
#include "GMS_integrands_func_ch5.h"

namespace file_info 
{

     static const unsigned int GMS_COMPUTE_FUNCTIONALS_CH5_MAJOR = 1;
     static const unsigned int GMS_COMPUTE_FUNCTIONALS_CH5_MINOR = 1;
     static const unsigned int GMS_COMPUTE_FUNCTIONALS_CH5_MICRO = 0;
     static const unsigned int GMS_COMPUTE_FUNCTIONALS_CH5_FULLVER =
       1000U*GMS_COMPUTE_FUNCTIONALS_CH5_MAJOR+100U*GMS_COMPUTE_FUNCTIONALS_CH5_MINOR+
       10U*GMS_COMPUTE_FUNCTIONALS_CH5_MICRO;
     static const char GMS_COMPUTE_FUNCTIONALS_CH5_CREATION_DATE[] = "28-08-2026 06:03AM +00200 (FRI 28 AUG 2026 GMT+2)";
     static const char GMS_COMPUTE_FUNCTIONALS_CH5_BUILD_DATE[]    = __DATE__; 
     static const char GMS_COMPUTE_FUNCTIONALS_CH5_BUILD_TIME[]    = __TIME__;
     static const char GMS_COMPUTE_FUNCTIONALS_CH5_SYNOPSIS[]      = "Numerical integration by the QUADPACK of the integrands (ch.4) used for the fading channel PDF/CDF and BER,BEP,SER,SEP quadrature computation(QUADPACK).\
	                                                          Based on the M.K Simon, M.S. Alouini: Digital Communication over Fading Channels 1st ed\
															  ISBN-13 978-0471317791";

}

 
// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (COMPUTE_FUNCTIONALS_CH5_USE_PMC_INSTRUMENTATION)
#define COMPUTE_FUNCTIONALS_CH5_USE_PMC_INSTRUMENTATION 0
#endif 

#if (COMPUTE_FUNCTIONALS_CH5_USE_PMC_INSTRUMENTATION) == 1
#include "GMS_hw_perf_macros.h"

#define PMC_VARS                      \
uint64_t prog_counters_start[4] = {}; \
uint64_t prog_counters_stop[4]  = {}; \
uint64_t tsc_start,tsc_stop;          \
uint64_t act_cyc_start,act_cyc_stop;  \
uint64_t ref_cyc_start,ref_cyc_stop;  \
[[maybe_unused]] uint64_t dummy1;     \
[[maybe_unused]] uint64_t dummy2;     \
[[maybe_unused]] uint64_t dummy3;     \
int32_t core_counter_width;           \
double utilization,nom_ghz,avg_ghz;
#endif 

namespace gms 
{

namespace fading_channel
{

struct alignas(64) quadpack_integrator_payload_ch5_t 
{
    const std::string integrators_names[6] = {{"dqage"},{"dqagi"},{"dqags"},{"dqng"},{"dqagp"},{"dqaws"}};
    const std::string tabulated_integrators_names[4] = {{"hiordq"},{"simpn"},{"wedint"},{"avint"}};
    double (*integrand)(double,void * __restrict__);
    func_args_ch5_payload_t * __restrict__ func_args_payload{nullptr};
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
    double              * __restrict__ functional{nullptr}; // computed (inner if applicable) functional values
    double              * __restrict__ outer_functional{nullptr}; // computed outer functional values 
    std::int32_t        * __restrict__ neval{nullptr};      // number of evaluation
    std::int32_t        * __restrict__ ier{nullptr};        // integrator error indicator
    std::int32_t        * __restrict__ last{nullptr};
    std::uint64_t       * __restrict__ crude_tsc_start{nullptr};
    std::uint64_t       * __restrict__ crude_tsc_end{nullptr};
    std::uint64_t       * __restrict__ crude_tsc_measurement{nullptr}; // crude RDTSCP measurement for approximate TSC assesment (shall not be used for the robust statistics)
    std::uint64_t       * __restrict__ crude_tsc_meas_outer{nullptr}; // crude RDTSCP measurement for approximate TSC assesment if the outer functional (shall not be used for the robust statistics)
    double                            outer_func_tmp_res; // result for outer functional computation 
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
    std::int32_t                       n_func_vals{}; 
    std::int32_t                       n_outer_func_vals{}; // outer functional number of values.
    std::int32_t                       which_integrator{}; //currently 1=dqage,2=dqagi,3=dqags,4=dqng,5=dqagp,6=dqaws
    std::int32_t                       which_tabulated_integrator{};
    bool                               randomly_generate_inputs{}; //as the name states: random input generation in use if true, otherwise provide deterministic inputs
};

/* 
    Different implementation of the above-defined major data structure
    for possibly ADT implementation.
*/

struct alignas(64) quadpack_integrator_payload_ch5_v2_t
{
    const std::string integrators_names[6] = {{"dqage"},{"dqagi"},{"dqags"},{"dqng"},{"dqagp"},{"dqaws"}};
    const std::string tabulated_integrators_names[4] = {{"hiordq"},{"simpn"},{"wedint"},{"avint"}};
    double (*integrand)(double,void * __restrict__);
    std::valarray<func_args_ch5_payload_t> func_args_payload;
    std::valarray<double>              tmp_work1;
    std::valarray<double>              tmp_work2;
    std::valarray<double>              tmp_work3;
    std::valarray<double>              tmp_work4;
    std::valarray<double>              tmp_work5;
    std::valarray<double>              tmp_work6;
    std::valarray<double>              lo;
    std::valarray<double>              hi;
    std::valarray<double>              bound;
    std::valarray<std::int32_t>        inf;
    std::valarray<std::int32_t>        irule;
    std::valarray<double>              epsabs;
    std::valarray<double>              epsrel;
    std::valarray<double>              abser;
    std::valarray<double>              functional;
    std::valarray<double>              outer_functional;
    std::valarray<std::int32_t>        neval;
    std::valarray<std::int32_t>        ier;
    std::valarray<std::int32_t>        last;
    std::valarray<std::uint64_t>       crude_tsc_start;
    std::valarray<std::uint64_t>       crude_tsc_end;
    std::valarray<std::uint64_t>       crude_tsc_measurements;
    std::valarray<std::uint64_t>       crude_tsc_meas_outer;
    double                            outer_func_tmp_res; // result for outer functional computation 
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
    std::int32_t                       n_func_vals{};
    std::int32_t                       n_outer_func_vals{}; // outer functional values of nest 
    std::int32_t                       which_integrator{}; //currently 1=dqage,2=dqagi,3=dqags,4=dqng,5=dqagp,6=dqaws
    std::int32_t                       which_tabulated_integrator{};
    bool                               randomly_generate_inputs{}; //as the name states: random input generation in use if true, otherwise provide deterministic inputs
};

template<std::size_t N,std::size_t M>
struct alignas(64) quadpack_integrator_payload_ch5_v3_t
{
    static_assert(N>=50ull,"The number of Functional values <=50!!");
    static_assert(M>=50ull,"The number of outer Functional values <=50!!");
    const std::string integrators_names[6] = {{"dqage"},{"dqagi"},{"dqags"},{"dqng"},{"dqagp"},{"dqaws"}};
    const std::string tabulated_integrators_names[4] = {{"hiordq"},{"simpn"},{"wedint"},{"avint"}};
    double (*integrand)(double,void * __restrict__);
    std::array<func_args_ch5_payload_t,N> funcs_arg_payload;
    std::array<double,N>              tmp_work1;
    std::array<double,N>              tmp_work2;
    std::array<double,N>              tmp_work3;
    std::array<double,N>              tmp_work4;
    std::array<double,N>              tmp_work5;
    std::array<double,N>              tmp_work6;
    std::array<double,N>              lo;
    std::array<double,N>              hi;
    std::array<double,N>              bound;
    std::array<std::int32_t,N>        inf;
    std::array<std::int32_t,N>        irule;
    std::array<double,N>              epsabs;
    std::array<double,N>              epsrel;
    std::array<double,N>              abser;
    std::array<double,N>              functional;
    std::array<double,M>              outer_functional;
    std::array<std::int32_t,N>        neval;
    std::array<std::int32_t,N>        ier;
    std::array<std::int32_t,N>        last;
    std::array<std::uint64_t,N>       crude_tsc_start;
    std::array<std::uint64_t,N>       crude_tsc_end;
    std::array<std::uint64_t,N>       crude_tsc_measurememnts;
    std::array<std::uint64_t,M>       crude_tsc_meas_outer; 
    double                            outer_func_tmp_res; // result for outer functional computation 
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
    std::int32_t                       n_func_vals{static_cast<std::int32_t>(N)};
    std::int32_t                       n_outer_func_vals{static_cast<std::int32_t>(M)}; // outer functional values.
    std::int32_t                       which_integrator; //currently 1=dqage,2=dqagi,3=dqags,4=dqng,5=dqagp,6=dqaws
    std::int32_t                       which_tabulated_integrator{};
    bool                               randomly_generate_inputs; //as the name states: random input generation in use if true, otherwise provide deterministic inputs   
};

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
std::int32_t 
compute_functional_Rayleigh_chan_5_6(quadpack_integrator_payload_ch5_t * __restrict__);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
std::int32_t
compute_functional_Hoyt_chan_5_9(quadpack_integrator_payload_ch5_t * __restrict__);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
std::int32_t
compute_functional_Rice_chan_5_12(quadpack_integrator_payload_ch5_t * __restrict__);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
std::int32_t
compute_functional_Nakagami_m_chan_5_16(quadpack_integrator_payload_ch5_t * __restrict__);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
std::int32_t
compute_functional_LogNormShadow_chan_5_20(quadpack_integrator_payload_ch5_t * __restrict__);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
std::int32_t
compute_outer_functional_LogNormShadow_chan_5_20(quadpack_integrator_payload_ch5_t * __restrict__);


} // fading_channel

} //gms

#endif /*__GMS_COMPUTE_FUNCTIONALS_CH5_H__*/