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

#ifndef __GMS_COMPUTE_FUNCTIONALS_CH4_H__
#define __GMS_COMPUTE_FUNCTIONALS_CH4_H__ 180820261108

#include <cstdint>
#include <string>
#include "GMS_config.h"
#include "GMS_integrands_func_ch4.h"

namespace file_info 
{

     static const unsigned int GMS_COMPUTE_FUNCTIONALS_CH4_MAJOR = 1;
     static const unsigned int GMS_COMPUTE_FUNCTIONALS_CH4_MINOR = 1;
     static const unsigned int GMS_COMPUTE_FUNCTIONALS_CH4_MICRO = 0;
     static const unsigned int GMS_COMPUTE_FUNCTIONALS_CH4_FULLVER =
       1000U*GMS_COMPUTE_FUNCTIONALS_CH4_MAJOR+100U*GMS_COMPUTE_FUNCTIONALS_CH4_MINOR+
       10U*GMS_COMPUTE_FUNCTIONALS_CH4_MICRO;
     static const char GMS_COMPUTE_FUNCTIONALS_CH4_CREATION_DATE[] = "18-08-2026 11:07AM +00200 (TUE 18 AUG 2026 GMT+2)";
     static const char GMS_COMPUTE_FUNCTIONALS_CH4_BUILD_DATE[]    = __DATE__; 
     static const char GMS_COMPUTE_FUNCTIONALS_CH4_BUILD_TIME[]    = __TIME__;
     static const char GMS_COMPUTE_FUNCTIONALS_CH4_SYNOPSIS[]      = "Numerical integration by the QUADPACK of the integrands (ch.4) used for the fading channel PDF/CDF and BER,BEP,SER,SEP quadrature computation(QUADPACK).\
	                                                          Based on the M.K Simon, M.S. Alouini: Digital Communication over Fading Channels 1st ed\
															  ISBN-13 978-0471317791";

}

 
// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (COMPUTE_FUNCTIONALS_CH4_USE_PMC_INSTRUMENTATION)
#define COMPUTE_FUNCTIONALS_CH4_USE_PMC_INSTRUMENTATION 0
#endif 

#if (COMPUTE_FUNCTIONALS_CH4_USE_PMC_INSTRUMENTATION) == 1
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

struct alignas(64) quadpack_integrator_payload_t 
{
    const std::string integrators_names[3] = {{"dqage"},{"dqagi"},{"dqags"}};
    double (*integrand)(double,void * __restrict__);
    func_args_payload_t * __restrict__ func_args_payload{nullptr};
    double              * __restrict__ tmp_work{nullptr}; //work storage (caller provided) used mainly for the functional a1st rguments sorting.
    double              * __restrict__ tmp_work2{nullptr};//work storage (caller provided) used mainly for the functional 2nd arguments sorting.
    double              * __restrict__ tmp_work3{nullptr};//work storage (caller provided) used mainly for the functional 3rd arguments sorting.
    double              * __restrict__ tmp_work4{nullptr};//work storage (caller provided) used mainly for the functional 4th arguments sorting.
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
    std::int32_t                       n_func_vals;
    std::int32_t                       which_integrator; //currently 1=dqage,2=dqagi,3=dqags
};

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
std::int32_t 
compute_functional_gauss_Q_4_1(quadpack_integrator_payload_t * __restrict__);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
std::int32_t
compute_functional_gauss_Q_4_2(quadpack_integrator_payload_t * __restrict__);

//double test_gauss_Q_4_2(double * );


} // fading_channel

} //gms



#endif /*__GMS_COMPUTE_FUNCTIONALS_CH4_H__*/