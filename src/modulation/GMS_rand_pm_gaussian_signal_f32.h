/* Copyright (C) Bernard Gingold, 2020-2026 
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

#ifndef __GMS_RAND_PM_GAUSSIAN_SIGNAL_F32_H__
#define __GMS_RAND_PM_GAUSSIAN_SIGNAL_F32_H__ 051120261002


#include <random>
#include <cstdint>
#include <vector>
#include "GMS_config.h"
#include "GMS_dyn_array.h"

namespace file_info 
{

     const unsigned int GMS_RAND_PM_GAUSSIAN_SIGNAL_MAJOR = 1;
     const unsigned int GMS_RAND_PM_GAUSSIAN_SIGNAL_MINOR = 0;
     const unsigned int GMS_RAND_PM_GAUSSIAN_SIGNAL_MICRO = 0;
     const unsigned int GMS_RAND_PM_GAUSSIAN_SIGNAL_FULLVER =
       1000U*GMS_RAND_PM_GAUSSIAN_SIGNAL_MAJOR+100U*GMS_RAND_PM_GAUSSIAN_SIGNAL_MINOR+
       10U*GMS_RAND_PM_GAUSSIAN_SIGNAL_MICRO;
     const char * const GMS_RAND_PM_GAUSSIAN_SIGNAL_CREATION_DATE = "05-11-2026 10:03 AM +00200 (THR 11 NOV 2026 GMT+2)";
     const char * const GMS_RAND_PM_GAUSSIAN_SIGNAL_BUILD_DATE    = __DATE__ " " __TIME__;
     static const char  GMS_RAND_PM_GAUSSIAN_SIGNAL_AUTHOR[]      = "Programmer: Bernard Gingold, beniekg@gmail.com"; 
     const char * const GMS_RAND_PM_GAUSSIAN_SIGNAL_SYNOPSIS      = "Randomly phase modulated sum of Gaussians, suitable for a noise sources modeling/simulation.";

}

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
// To be added.
#if !defined (RAND_PM_GAUSSIAN_SIGNAL_F32_USE_NT_STORES)
#define RAND_PM_GAUSSIAN_SIGNAL_F32_USE_NT_STORES 0
#endif

#ifdef _OPENMP
// Default init a storage arrays for the first-touch (OpenMP) processing
#if !defined(RAND_PM_GAUSSIAN_SIGNAL_F32_INIT_STORAGE)
#define RAND_PM_GAUSSIAN_SIGNAL_F32_INIT_STORAGE 1
#endif 
#endif

#if (RAND_PM_GAUSSIAN_SIGNAL_F32_INIT_STORAGE) == 1
#define INIT_BY_STD_FILL 0
#endif 

// For inlining of trigo functions (asin,acos,sin,cos)
#if !defined(RAND_PM_GAUSSIAN_SIGNAL_F32_USE_CEPHES)
#define RAND_PM_GAUSSIAN_SIGNAL_F32_USE_CEPHES 1
#endif 

// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (RAND_PM_GAUSSIAN_SIGNAL_F32_USE_PMC_INSTRUMENTATION)
#define RAND_PM_GAUSSIAN_SIGNAL_F32_USE_PMC_INSTRUMENTATION 0
#endif 

#if (RAND_PM_GAUSSIAN_SIGNAL_F32_USE_PMC_INSTRUMENTATION) == 1
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

#if (RAND_PM_GAUSSIAN_SIGNAL_F32_USE_CEPHES) == 0
#include <cmath>
#else 
#include "GMS_cephes_double.h"
#endif

namespace gms
{

namespace radiolocation 
{

struct alignas(64) rand_pm_gaussian_signal_f32_t
{
    std::size_t m_nsamples;
    std::size_t m_num_gaussians;
    float       m_tstart;
    float       m_tend;
    float       m_min_sigma;
    float       m_max_sigma;
    float       m_max_amplitude;
    darray_r4_t m_signal_samples;

    rand_pm_gaussian_signal_f32_t() noexcept(true); 
    
    rand_pm_gaussian_signal_f32_t(const std::size_t,
                                  const std::size_t,
                                  const float,
                                  const float,
                                  const float,
                                  const float,
                                  const float) noexcept(false);

    rand_pm_gaussian_signal_f32_t(const rand_pm_gaussian_signal_f32_t &) noexcept(false);

    rand_pm_gaussian_signal_f32_t(rand_pm_gaussian_signal_f32_t &&) noexcept(true);

    ~rand_pm_gaussian_signal_f32_t() noexcept(true);

    rand_pm_gaussian_signal_f32_t &
                    operator=(const rand_pm_gaussian_signal_f32_t &) noexcept(false);

    rand_pm_gaussian_signal_f32_t &
                    operator=(rand_pm_gaussian_signal_f32_t &&) noexcept(true);

    auto operator==(const rand_pm_gaussian_signal_f32_t &)->std::vector<bool>;
                                                                                    
    auto operator!=(const rand_pm_gaussian_signal_f32_t &)->std::vector<bool>;
                                                                                    
    auto operator>=(const rand_pm_gaussian_signal_f32_t &)->std::vector<bool>;
                                                                                   
    auto operator<=(const rand_pm_gaussian_signal_f32_t &)->std::vector<bool>;

    __ATTR_ALWAYS_INLINE__ 
    inline
    bool init_samples_buffer(const float filler)
    {
        bool ret_status = false;
#if (INIT_BY_STD_FILL) == 0
        using namespace gms::common;
	    sse_memset_unroll8x_ps(&this->m_signal_samples.m_data[0],filler,this->m_nsamples);
#else 
        std::fill(this->m_signal_samples.m_data,this->m_sig_samples.m_data+this->m_nsamples,filler);
#endif 
        ret_status = true;
        return (ret_status);
    } 

    std::int32_t generate();
    
};

    auto operator<<(std::ostream &,
                    const rand_pm_gaussian_signal_f32_t &)->std::ostream &;

}

}

#endif /*__GMS_RAND_PM_GAUSSIAN_SIGNAL_F32_H__*/