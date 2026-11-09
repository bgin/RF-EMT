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

#ifndef __GMS_RAND_PM_GAUSSIAN_SIGNAL_F64_H__
#define __GMS_RAND_PM_GAUSSIAN_SIGNAL_F64_H__ 051120261002

#include <iostream>
#include <random>
#include <cstdint>
#include <vector>
#include "GMS_config.h"
#include "GMS_dyn_array.h"
#include "GMS_sse_memset.h"

namespace file_info 
{

     const unsigned int GMS_RAND_PM_GAUSSIAN_SIGNAL_F64_MAJOR = 1;
     const unsigned int GMS_RAND_PM_GAUSSIAN_SIGNAL_F64_MINOR = 0;
     const unsigned int GMS_RAND_PM_GAUSSIAN_SIGNAL_F64_MICRO = 0;
     const unsigned int GMS_RAND_PM_GAUSSIAN_SIGNAL_F64_FULLVER =
       1000U*GMS_RAND_PM_GAUSSIAN_SIGNAL_F64_MAJOR+100U*GMS_RAND_PM_GAUSSIAN_SIGNAL_F64_MINOR+
       10U*GMS_RAND_PM_GAUSSIAN_SIGNAL_F64_MICRO;
     const char * const GMS_RAND_PM_GAUSSIAN_SIGNAL_F64_CREATION_DATE = "05-11-2026 10:03 AM +00200 (THR 11 NOV 2026 GMT+2)";
     const char * const GMS_RAND_PM_GAUSSIAN_SIGNAL_F64_BUILD_DATE    = __DATE__ " " __TIME__;
     static const char  GMS_RAND_PM_GAUSSIAN_SIGNAL_F64_AUTHOR[]      = "Programmer: Bernard Gingold, beniekg@gmail.com"; 
     const char * const GMS_RAND_PM_GAUSSIAN_SIGNAL_F64_SYNOPSIS      = "Randomly phase modulated sum of Gaussians, suitable for a noise sources modeling/simulation.";

}

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
// To be added.
#if !defined (RAND_PM_GAUSSIAN_SIGNAL_F64_USE_NT_STORES)
#define RAND_PM_GAUSSIAN_SIGNAL_F64_USE_NT_STORES 0
#endif

#ifdef _OPENMP
// Default init a storage arrays for the first-touch (OpenMP) processing
#if !defined(RAND_PM_GAUSSIAN_SIGNAL_F64_INIT_STORAGE)
#define RAND_PM_GAUSSIAN_SIGNAL_F64_INIT_STORAGE 1
#endif 
#endif

#if (RAND_PM_GAUSSIAN_SIGNAL_F64_INIT_STORAGE) == 1
#define INIT_BY_STD_FILL 0
#endif 

// For inlining of trigo functions (asin,acos,sin,cos)
#if !defined(RAND_PM_GAUSSIAN_SIGNAL_F64_USE_CEPHES)
#define RAND_PM_GAUSSIAN_SIGNAL_F64_USE_CEPHES 1
#endif 

// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (RAND_PM_GAUSSIAN_SIGNAL_F64_USE_PMC_INSTRUMENTATION)
#define RAND_PM_GAUSSIAN_SIGNAL_F64_USE_PMC_INSTRUMENTATION 0
#endif 

#if (RAND_PM_GAUSSIAN_SIGNAL_F64_USE_PMC_INSTRUMENTATION) == 1
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

#if (RAND_PM_GAUSSIAN_SIGNAL_F64_USE_CEPHES) == 0
#include <cmath>
#else 
#include "GMS_cephes_double.h"
#endif

namespace gms
{

namespace radiolocation 
{

struct alignas(64) rand_pm_gaussian_signal_f64_t
{
    std::size_t m_nsamples;
    std::size_t m_num_gaussians;
    double       m_tstart;
    double       m_tend;
    double       m_min_sigma;
    double       m_max_sigma;
    double       m_max_amplitude;
    darray_r4_t m_signal_samples;

    rand_pm_gaussian_signal_f64_t() noexcept(true); 
    
    rand_pm_gaussian_signal_f64_t(const std::size_t,
                                  const std::size_t,
                                  const double,
                                  const double,
                                  const double,
                                  const double,
                                  const double) noexcept(false);

    rand_pm_gaussian_signal_f64_t(const rand_pm_gaussian_signal_f64_t &) noexcept(false);

    rand_pm_gaussian_signal_f64_t(rand_pm_gaussian_signal_f64_t &&) noexcept(true);

    ~rand_pm_gaussian_signal_f64_t() noexcept(true);

    rand_pm_gaussian_signal_f64_t &
                    operator=(const rand_pm_gaussian_signal_f64_t &) noexcept(false);

    rand_pm_gaussian_signal_f64_t &
                    operator=(rand_pm_gaussian_signal_f64_t &&) noexcept(true);

    auto operator==(const rand_pm_gaussian_signal_f64_t &)->std::vector<bool>;
                                                                                    
    auto operator!=(const rand_pm_gaussian_signal_f64_t &)->std::vector<bool>;
                                                                                    
    auto operator>=(const rand_pm_gaussian_signal_f64_t &)->std::vector<bool>;
                                                                                   
    auto operator<=(const rand_pm_gaussian_signal_f64_t &)->std::vector<bool>;

    __ATTR_ALWAYS_INLINE__ 
    inline
    bool init_samples_buffer(const double filler)
    {
        bool ret_status = false;
#if (INIT_BY_STD_FILL) == 0
        using namespace gms::common;
	    sse_memset_unroll8x_pd(&this->m_signal_samples.m_data[0],filler,this->m_nsamples);
#else 
        std::fill(this->m_signal_samples.m_data,this->m_sig_samples.m_data+this->m_nsamples,filler);
#endif 
        ret_status = true;
        return (ret_status);
    } 
    
    template<std::size_t N>
    std::int32_t generate()
    {
        std::int32_t ret_status = -9999;
        if(__builtin_expect(this->m_num_gaussians>N,0)) { ret_status = -1; return (ret_status);}
        struct alignas(16) gaussian_components_t 
        {
           double amplitude;
           double mu;
           double two_sigma_sq;
        };
        alignas(64) gaussian_components_t gc_store[N];
        thread_local std::uniform_real_distribution<double> dist_mu;
        thread_local std::uint64_t seed_mu{};
        thread_local std::mt19937  dist_mu_gen;
        thread_local std::uniform_real_distribution<double> dist_sigma;
        thread_local std::uint64_t seed_sigma{};
        thread_local std::mt19937 dist_sigma_gen;
        thread_local std::uniform_real_distribution<double> dist_amp;
        thread_local std::uint64_t seed_amp{};
        thread_local std::mt19937 dist_amp_gen;
        const double  padding = 3.0*this->m_max_sigma;
        seed_mu = __rdtsc();
        dist_mu = std::uniform_real_distribution<double>(this->m_tstart-padding,this->m_tend+padding);
        dist_mu_gen = std::mt19937(seed_mu);
        seed_sigma = __rdtsc();
        dist_sigma = std::uniform_real_distribution<double>(this->m_min_sigma,this->m_max_sigma);
        dist_sigma_gen = std::mt19937(seed_sigma);
        seed_amp   = __rdtsc();
        dist_amp   = std::uniform_real_distribution<double>(-this->m_max_amplitude,this->m_max_amplitude);
        dist_amp_gen   = std::mt19937(seed_amp);
        for(std::size_t i = 0; i<this->m_num_gaussians; ++i) 
        {
            gc_store[i].amplitude    = dist_amp.operator()(dist_amp_gen);
            gc_store[i].mu           = dist_mu.operator()(dist_mu_gen);
            double sigma_tmp          = dist_sigma.operator()(dist_sigma_gen);
            gc_store[i].two_sigma_sq = 2.0f*sigma_tmp*sigma_tmp; 
        }
        const double dt = (this->m_tend-this->m_tstart)/static_cast<double>(this->m_nsamples-1ull);
        constexpr double sqrt2PI = 2.506628274631000502415765284811;
        for(std::size_t j = 0ull; j<this->m_nsamples; ++j)   
        {
            const double t_j = static_cast<double>(j);
            const double t   = this->m_tstart+t_j*dt;
            double phase_accumulator = 0.0f;
            for(std::size_t k = 0; k<this->m_num_gaussians; ++k) 
            {   
                const double tmu           = gc_store[k].mu;
                double diff                = t-tmu;
                const double tamplitude    = gc_store[k].amplitude;
                const double tsigma        = gc_store[k].two_sigma_sq;
                const double normalization = 1.0f/(tsigma*sqrt2PI);
#if (RAND_PM_GAUSSIAN_SIGNAL_F64_USE_CEPHES) == 1
                const double exp_val    = normalization*gms::math::cephes_d::exp(-(diff*diff)/tsigma);
#else 
                const double exp_val    = normalization*std::exp(-(diff*diff)/tsigma);
#endif 
                phase_accumulator += (tamplitude*exp_val);
            }
            this->m_signal_samples.m_data[j] = phase_accumulator;
        }
        ret_status = 0;
        return (ret_status);
    }
    
};

    auto operator<<(std::ostream &,
                    const rand_pm_gaussian_signal_f64_t &)->std::ostream &;

}

}

#endif /*__GMS_RAND_PM_GAUSSIAN_SIGNAL_F64_H__*/