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

#ifndef __GMS_SKX_CPU_OPERATING_FREQ_METRIC_H__
#define __GMS_SKX_CPU_OPERATING_FREQ_METRIC_H__ 221020260552

#include <cstdint>
#include <string>
#include <valarray> // used instead of vector for non-runtime relocatable guarantee
#include "GMS_config.h"

namespace file_info 
{

     static const unsigned int GMS_SKX_CPU_OPERATING_FREQ_METRIC_MAJOR = 1;
     static const unsigned int GMS_SKX_CPU_OPERATING_FREQ_METRIC_MINOR = 1;
     static const unsigned int GMS_SKX_CPU_OPERATING_FREQ_METRIC_MICRO = 0;
     static const unsigned int GMS_SKX_CPU_OPERATING_FREQ_METRIC_FULLVER =
       1000U*GMS_SKX_CPU_OPERATING_FREQ_METRIC_MAJOR+100U*GMS_SKX_CPU_OPERATING_FREQ_METRIC_MINOR+
       10U*GMS_SKX_CPU_OPERATING_FREQ_METRIC_MICRO;
     static const char GMS_SKX_CPU_OPERATING_FREQ_METRIC_CREATION_DATE[]  = "22-10-2026 05:54AM +00200 (THR 22 OCT 2026 GMT+2)";
     static const char GMS_SKX_CPU_OPERATING_FREQ_METRIC_BUILD_DATE[]    = __DATE__; 
     static const char GMS_SKX_CPU_OPERATING_FREQ_METRIC_BUILD_TIME[]    = __TIME__;
     static const char GMS_SKX_CPU_OPERATING_FREQ_METRIC_AUTHOR[]        = "Programmer: Bernard Gingold, beniekg@gmail.com";
     static const char GMS_SKX_CPU_OPERATING_FREQ_METRIC_SYNOPSIS[]      = "Preprocessing of Skylake-X (uarch:server) CPU OPERATING FREQUENCY metric.";
}

namespace gms 
{

namespace hw_profiling
{

struct alignas(64) SKX_cpu_operating_frequency_t 
{
    
    std::size_t           m_nruns;
    std::size_t           m_nsamples;
    std::size_t           m_totsamples;
    std::valarray<double> m_a; // a - event 
    std::valarray<double> m_b; // b - event
    std::valarray<double> m_c; // c - event
    std::valarray<double> m_d; // d - metric calculation result
    std::string           m_cpu_arch;
    std::string           m_date; // date of processing
    std::string           m_time; // time of processing
    double                m_a_robust_mean; // Mathematical-Expectation
    double                m_b_robust_mean; // Mathematical-Expectation
    double                m_c_robust_mean; // Mathematical-Expectation
    double                m_d_robust_mean; // Mathematical-Expectation
    double                m_a_normality; // Shapiro-Wilk normality test
    double                m_b_normality; // Shapiro-Wilk normality test
    double                m_c_normality; // Shapiro-Wilk normality test
    double                m_d_normality; // Shapiro-Wilk normality test
    struct alignas(64)    metric_metadata_t
    {
       const std::string     m_brief_descr = "CPU operating frequency (in GHz)";
       const std::string     m_metric_expr = "((CPU_CLK_UNHALTED.THREAD / CPU_CLK_UNHALTED.REF_TSC * #SYSTEM_TSC_FREQ ) / 1000000000)";
       const std::string     m_metric_group= "";
       const std::string     m_metric_name = "cpu_operating_frequency";
       const std::string     m_scale_unit  = "1GHz";
    } metadata_t;
    
    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    SKX_cpu_operating_frequency_t() noexcept(true);
    
    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    SKX_cpu_operating_frequency_t(const std::size_t nruns,
                                  const std::size_t nsamples,
                                  const std::string &,
                                  const std::string,
                                  const std::string) noexcept(false);
    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    SKX_cpu_operating_frequency_t(const SKX_cpu_operating_frequency_t &) noexcept(false);
    
    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    SKX_cpu_operating_frequency_t(SKX_cpu_operating_frequency_t &&) noexcept(true);
    
    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    ~SKX_cpu_operating_frequency_t();
    
    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    SKX_cpu_operating_frequency_t &
                      operator=(const SKX_cpu_operating_frequency_t &) noexcept(false);
    
    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    SKX_cpu_operating_frequency_t &
                      operator=(SKX_cpu_operating_frequency_t &&) noexcept(true);

    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    auto operator==(const SKX_cpu_operating_frequency_t &)->std::valarray<bool>;
                    
    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    auto operator!=(const SKX_cpu_operating_frequency_t &)->std::valarray<bool>;
                   
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

    __ATTR_ALWAYS_INLINE__ 
    inline const double * __restrict__ get_m_d_ptr() const 
    {
        return reinterpret_cast<const double* __restrict__>(&this->m_d.operator[](0ull));
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
         convert_m_c_to_fp64_u8x(const std::valarray<std::uint64_t>&) noexcept(true);

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
         convert_m_c_to_fp64_u8x_avx512(const std::valarray<std::uint64_t>&) noexcept(true);

    
    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    std::int32_t compute_metric_u4x() noexcept(true);
         
};

}

}

#endif /*__GMS_SKX_CPU_OPERATING_FREQ_METRIC_H__*/