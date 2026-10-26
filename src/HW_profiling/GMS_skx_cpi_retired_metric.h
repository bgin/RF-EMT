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

#ifndef __GMS_SKX_CPI_RETIRED_METRIC_H__
#define __GMS_SKX_CPI_RETIRED_METRIC_H__ 23102021248

#include <cstdint>
#include <string>
#include <valarray> // used instead of vector for non-runtime relocatable guarantee
#if 0 /*Currently inactive due to compilation errors (g++,icpc)*/
#if !defined(__ICC) || !defined(__INTEL_COMPILER)
#include <boost/histogram.hpp>
#endif 
#endif 
#include "GMS_config.h"

namespace file_info 
{

     static const unsigned int GMS_SKX_CPI_RETIRED_METRIC_MAJOR = 1;
     static const unsigned int GMS_SKX_CPI_RETIRED_METRIC_MINOR = 1;
     static const unsigned int GMS_SKX_CPI_RETIRED_METRIC_MICRO = 0;
     static const unsigned int GMS_SKX_CPI_RETIRED_METRIC_FULLVER =
       1000U*GMS_SKX_CPI_RETIRED_METRIC_MAJOR+100U*GMS_SKX_CPI_RETIRED_METRIC_MINOR+
       10U*GMS_SKX_CPI_RETIRED_METRIC_MICRO;
     static const char GMS_SKX_CPI_RETIRED_METRIC_CREATION_DATE[]  = "23-10-2026 12:50PM +00200 (FRI 23 OCT 2026 GMT+2)";
     static const char GMS_SKX_CPI_RETIRED_METRIC_BUILD_DATE[]    = __DATE__; 
     static const char GMS_SKX_CPI_RETIRED_METRIC_BUILD_TIME[]    = __TIME__;
     static const char GMS_SKX_CPI_RETIRED_METRIC_AUTHOR[]        = "Programmer: Bernard Gingold, beniekg@gmail.com";
     static const char GMS_SKX_CPI_RETIRED_METRIC_SYNOPSIS[]      = "Preprocessing of Skylake-X (uarch:server) CPU Cycles per instruction retired metric.";
}

#if !defined(SKX_CPI_RETIRED_METRIC_VERBOSE)
#define SKX_CPI_RETIRED_METRIC_VERBOSE 0
#endif 

namespace gms
{

namespace hw_profiling
{

/*
//////////////// Begin: put this in header file ////////////////////

#include <algorithm>           // std::max_element
#include <boost/format.hpp>    // only needed for printing
#include <boost/histogram.hpp> // make_histogram, integer, indexed
#include <iostream>            // std::cout, std::endl
#include <sstream>             // std::ostringstream
#include <tuple>
#include <vector>

// use this when axis configuration is fix to get highest performance
struct HolderOfStaticHistogram {
  // put axis types here
  using axes_t = std::tuple<
    boost::histogram::axis::regular<>,
    boost::histogram::axis::integer<>
  >;
  using hist_t = boost::histogram::histogram<axes_t>;
  hist_t hist_;
};

// use this when axis configuration should be flexible
struct HolderOfDynamicHistogram {
  // put all axis types here that you are going to use
  using axis_t = boost::histogram::axis::variant<
    boost::histogram::axis::regular<>,
    boost::histogram::axis::variable<>,
    boost::histogram::axis::integer<>
  >;
  using axes_t = std::vector<axis_t>;
  using hist_t = boost::histogram::histogram<axes_t>;
  hist_t hist_;
};

//////////////// End: put this in header file ////////////////////

int main() {
  using namespace boost::histogram;

  HolderOfStaticHistogram hs;
  hs.hist_ = make_histogram(axis::regular<>(5, 0, 1), axis::integer<>(0, 3));
  // now assign a different histogram
  hs.hist_ = make_histogram(axis::regular<>(3, 1, 2), axis::integer<>(4, 6));
  // hs.hist_ = make_histogram(axis::regular<>(5, 0, 1)); does not work;
  // the static histogram cannot change the number or order of axis types

  HolderOfDynamicHistogram hd;
  hd.hist_ = make_histogram(axis::regular<>(5, 0, 1), axis::integer<>(0, 3));
  // now assign a different histogram
  hd.hist_ = make_histogram(axis::regular<>(3, -1, 2));
  // and assign another
  hd.hist_ = make_histogram(axis::integer<>(0, 5), axis::integer<>(3, 5));
}
*/



struct alignas(64) SKX_cpu_cpi_retired_t 
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
    std::uint32_t         m_cpuid_basic[8] = {0}; //basic cpuid info
#if 0 /*Boost Histogram compilation error!!*/
#if !defined(__ICC) && !defined(__INTEL_COMPILER)
    struct alignas(64) metric_histograms_t
    {
       using axis_reg_t = boost::histogram::axis::regular<>;
       boost::histogram::histogram<axis_reg_t> m_hist_a;
       boost::histogram::histogram<axis_reg_t> m_hist_b;
       boost::histogram::histogram<axis_reg_t> m_hist_c;
    } static m_histograms_t;
#endif 
#endif 
    struct alignas(64)    metric_metadata_t
    {
       const std::string     m_brief_descr = "Cycles per instruction retired; indicating how much time each executed instruction took; in units of cycles";
       const std::string     m_metric_expr = "CPU_CLK_UNHALTED.THREAD / INST_RETIRED.ANY";
       const std::string     m_metric_group= "";
       const std::string     m_metric_name = "cpi";
       const std::string     m_scale_unit  = "1-per instruction";
    } metadata_t;
    
    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    SKX_cpu_cpi_retired_t() noexcept(true);
    
    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    SKX_cpu_cpi_retired_t(     const std::size_t nruns,
                               const std::size_t nsamples,
                               const std::string,
                               const std::string,
                               const std::uint32_t * __restrict__) noexcept(false);
    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    SKX_cpu_cpi_retired_t(const SKX_cpu_cpi_retired_t &) noexcept(false);
    
    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    SKX_cpu_cpi_retired_t(SKX_cpu_cpi_retired_t &&) noexcept(true);
    
    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    ~SKX_cpu_cpi_retired_t();
    
    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    SKX_cpu_cpi_retired_t &
                      operator=(const SKX_cpu_cpi_retired_t &) noexcept(false);
    
    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    SKX_cpu_cpi_retired_t &
                      operator=(SKX_cpu_cpi_retired_t &&) noexcept(true);

    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    auto operator==(const SKX_cpu_cpi_retired_t &)->std::valarray<bool>;
                    
    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    auto operator!=(const SKX_cpu_cpi_retired_t &)->std::valarray<bool>;
                   
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
        parse_cpuid_basic_info(const bool);

    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    std::int32_t compute_metric_u4x() noexcept(true);
#if 0
#if !defined(__ICC) || !defined(__INTEL_COMPILER)
    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    std::int32_t process_histogram_data(const double,
                                        const double,
                                        const double,
                                        const double,
                                        const double,
                                        const double);
#endif 
#endif 
         
};

}

}

#endif /*__GMS_SKX_CPI_RETIRED_METRIC_H__*/