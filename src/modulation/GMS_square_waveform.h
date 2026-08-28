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

#ifndef __GMS_SQUARE_WAVEFORM_H__
#define __GMS_SQUARE_WAVEFORM_H__


namespace file_info 
{

     static const unsigned int GMS_SQUARE_WAVEFORM_MAJOR = 1;
     static const unsigned int GMS_SQUARE_WAVEFORM_MINOR = 1;
     static const unsigned int GMS_SQUARE_WAVEFORM_MICRO = 0;
     static const unsigned int GMS_SQUARE_WAVEFORM_FULLVER =
       1000U*GMS_SQUARE_WAVEFORM_MAJOR+100U*GMS_SQUARE_WAVEFORM_MINOR+
       10U*GMS_SQUARE_WAVEFORM_MICRO;
     static const char GMS_SQUARE_WAVEFORM_CREATION_DATE[] = "24-09-2025 08:43 +00200 (WED 24 SEP 2025 GMT+2)";
     static const char GMS_SQUARE_WAVEFORM_BUILD_DATE[]    = __DATE__; 
     static const char GMS_SQUARE_WAVEFORM_BUILD_TIME[]    = __TIME__;
     static const char GMS_SQUARE_WAVEFORM_SYNOPSIS[]      = "Square waveform generators.";

}

#include <cstdint>
#include <string>
#include <iostream>
#include "GMS_config.h"
#include "GMS_dyn_array.h"

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
// To be added.
#if !defined (SQUARE_WAVEFORM_USE_NT_STORES)
#define SQUARE_WAVEFORM_USE_NT_STORES 0
#endif

#ifdef _OPENMP
// Default init a storage arrays for the first-touch (OpenMP) processing
#if !defined(SQUARE_WAVEFORM_INIT_STORAGE)
#define SQUARE_WAVEFORM_INIT_STORAGE 1
#endif 
#endif

#if (SQUARE_WAVEFORM_INIT_STORAGE) == 1
#define INIT_BY_STD_FILL 0
#endif 

// For inlining of trigo functions (asin,acos,sin,cos)
#if !defined(SQUARE_WAVEFORM_USE_CEPHES)
#define SQUARE_WAVEFORM_USE_CEPHES 1
#endif 

// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (SQUARE_WAVEFORM_USE_PMC_INSTRUMENTATION)
#define SQUARE_WAVEFORM_USE_PMC_INSTRUMENTATION 1
#endif 

#if (SQUARE_WAVEFORM_USE_PMC_INSTRUMENTATION) == 1
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

namespace radiolocation 
{
          struct alignas(64) square_waveform_t final 
          {
                 std::size_t   __n_samples__;
                 std::uint32_t __n_waves__;
                 float         __A__; // amplitude 
                 float         __T__; // period
                 darray_r4_t   __sw_samples__;

                 square_waveform_t() = default;

                 square_waveform_t(const std::size_t,
                                   const std::uint32_t,
                                   const float,
                                   const float) noexcept(false);

                 square_waveform_t(square_waveform_t &&) noexcept(true);

                 ~square_waveform_t() noexcept(false);

                 square_waveform_t & operator=(const square_waveform_t &) = delete;

                 square_waveform_t & operator=(square_waveform_t &&);

                 void init_storage(const float);

                 static void create_signal_plot(const std::uint32_t,
                                                 const float * __restrict,
                                                 const float * __restrict,
                                                 const std::string &,
                                                 const std::string &,
                                                 const bool );

                void fourier_series_expansion();

        };

          auto 
          operator<<(std::ostream &,
                     square_waveform_t &)->std::ostream &;
}

}















#endif /*__GMS_SQUARE_WAVEFORM_H__*/
