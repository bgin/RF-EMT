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

#ifndef __GMS_SKX_CPU_UTILIZATION_METRIC_H__
#define __GMS_SKX_CPU_UTILIZATION_METRIC_H__ 231020260604

#include <cstdint>
#include <string>
#include <valarray> // used instead of vector for non-runtime relocatable guarantee
#include "GMS_config.h"

namespace file_info 
{

     static const unsigned int GMS_SKX_CPU_UTILIZATION_METRIC_MAJOR = 1;
     static const unsigned int GMS_SKX_CPU_UTILIZATION_METRIC_MINOR = 1;
     static const unsigned int GMS_SKX_CPU_UTILIZATION_METRIC_MICRO = 0;
     static const unsigned int GMS_SKX_CPU_UTILIZATION_METRIC_FULLVER =
       1000U*GMS_SKX_CPU_UTILIZATION_METRIC_MAJOR+100U*GMS_SKX_CPU_UTILIZATION_METRIC_MINOR+
       10U*GMS_SKX_CPU_UTILIZATION_METRIC_MICRO;
     static const char GMS_SKX_CPU_UTILIZATION_METRIC_CREATION_DATE[]  = "23-10-2026 06:04AM +00200 (FRI 23 OCT 2026 GMT+2)";
     static const char GMS_SKX_CPU_UTILIZATION_METRIC_BUILD_DATE[]    = __DATE__; 
     static const char GMS_SKX_CPU_UTILIZATION_METRIC_BUILD_TIME[]    = __TIME__;
     static const char GMS_SKX_CPU_UTILIZATION_METRIC_AUTHOR[]        = "Programmer: Bernard Gingold, beniekg@gmail.com";
     static const char GMS_SKX_CPU_UTILIZATION_METRIC_SYNOPSIS[]      = "Preprocessing of Skylake-X (uarch:server) CPU CORE UTILIZATION metric.";
}

namespace gms
{

namespace hw_profiling
{

struct alignas(64) SKX_cpu_core_utilization_t 
{
    
    std::size_t           m_nruns;
    std::size_t           m_nsamples;
    std::size_t           m_totsamples;
    std::valarray<double> m_a; // a - event 
    std::valarray<double> m_b; // b - event
    std::valarray<double> m_c; // c - metric calculation result
    std::string           m_date; // date of processing
    std::string           m_time; // time of processing
    double                m_a_robust_mean; // Mathematical-Expectation
    double                m_b_robust_mean; // Mathematical-Expectation
    double                m_c_robust_mean; // Mathematical-Expectation
    double                m_a_normality; // Shapiro-Wilk normality test
    double                m_b_normality; // Shapiro-Wilk normality test
    double                m_c_normality; // Shapiro-Wilk normality test
    std::uint32_t         m_cpuid_basic[8]; //basic cpuid info
    struct alignas(64)    metric_metadata_t
    {
       const std::string     m_brief_descr = "Percentage of time spent in the active CPU power state C0";
       const std::string     m_metric_expr = "CPU_CLK_UNHALTED.REF_TSC / TSC";
       const std::string     m_metric_group= "";
       const std::string     m_metric_name = "cpu_utilization";
       const std::string     m_scale_unit  = "100%";
    } metadata_t;
    
    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    SKX_cpu_core_utilization_t() noexcept(true);
    
    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    SKX_cpu_core_utilization_t(const std::size_t nruns,
                               const std::size_t nsamples,
                               const std::string,
                               const std::string,
                               const std::uint32_t * __restrict__) noexcept(false);
    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    SKX_cpu_core_utilization_t(const SKX_cpu_core_utilization_t &) noexcept(false);
    
    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    SKX_cpu_core_utilization_t(SKX_cpu_core_utilization_t &&) noexcept(true);
    
    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    ~SKX_cpu_core_utilization_t();
    
    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    SKX_cpu_core_utilization_t &
                      operator=(const SKX_cpu_core_utilization_t &) noexcept(false);
    
    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    SKX_cpu_core_utilization_t &
                      operator=(SKX_cpu_core_utilization_t &&) noexcept(true);

    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    auto operator==(const SKX_cpu_core_utilization_t &)->std::valarray<bool>;
                    
    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    auto operator!=(const SKX_cpu_core_utilization_t &)->std::valarray<bool>;
                   
    __ATTR_ALWAYS_INLINE__
    inline const double * __restrict__ get_m_a_ptr() const 
    {
        return reinterpret_cast<const double* __restrict__>(&this->m_a.operator[](0ull));
    }

    __ATTR_ALWAYS_INLINE__ 
    inline const double * __restrict__ get_m_b_ptr() const 
    {
        return reinterpret_cast<const double* __restrict__>(&this->m_b.operator[](0ull));
    }

    __ATTR_ALWAYS_INLINE__
    inline const double * __restrict__ get_m_c_ptr() const 
    {
        return reinterpret_cast<const double* __restrict__>(&this->m_c.operator[](0ull));
    }

    
    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    std::int32_t 
         convert_m_a_to_fp64_u8x(const std::valarray<std::uint64_t>&) noexcept(true);

    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    std::int32_t 
         convert_m_b_to_fp64_u8x(const std::valarray<std::uint64_t>&) noexcept(true);

    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    std::int32_t 
         convert_m_a_to_fp64_u8x_avx512(const std::valarray<std::uint64_t>&) noexcept(true);

    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    std::int32_t 
         convert_m_b_to_fp64_u8x_avx512(const std::valarray<std::uint64_t>&) noexcept(true);

    __ATTR_HOT__ 
    __ATTR_ALIGN__(32)
    std::int32_t 
        parse_cpuid_basic_info();

    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    std::int32_t compute_metric_u4x() noexcept(true);
         
};

}

}

#endif /*__GMS_SKX_CPU_UTILIZATION_METRIC_H__*/