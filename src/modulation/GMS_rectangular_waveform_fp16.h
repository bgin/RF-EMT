
/*MIT License
!Copyright (c) 2020 Bernard Gingold
!Permission is hereby granted, free of charge, to any person obtaining a copy
!of this software and associated documentation files (the "Software"), to deal
!in the Software without restriction, including without limitation the rights
!to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!copies of the Software, and to permit persons to whom the Software is
!furnished to do so, subject to the following conditions:
!The above copyright notice and this permission notice shall be included in all
!copies or substantial portions of the Software.
!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!SOFTWARE.
*/

#ifndef __GMS_RECTANGULAR_WAVEFORM_FP16_H__
#define __GMS_RECTANGULAR_WAVEFORM_FP16_H__ 221020250858

namespace file_info 
{

     static const unsigned int GMS_RECTANGULAR_WAVEFORM_FP16_MAJOR = 1;
     static const unsigned int GMS_RECTANGULAR_WAVEFORM_FP16_MINOR = 1;
     static const unsigned int GMS_RECTANGULAR_WAVEFORM_FP16_MICRO = 0;
     static const unsigned int GMS_RECTANGULAR_WAVEFORM_FP16_FULLVER =
       1000U*GMS_RECTANGULAR_WAVEFORM_FP16_MAJOR+100U*GMS_RECTANGULAR_WAVEFORM_FP16_MINOR+
       10U*GMS_RECTANGULAR_WAVEFORM_FP16_MICRO;
     static const char GMS_RECTANGULAR_WAVEFORM_FP16_CREATION_DATE[] = "22-10-2025 08:58 +00200 (WED 22 OCT 2025 GMT+2)";
     static const char GMS_RECTANGULAR_WAVEFORM_FP16_BUILD_DATE[]    = __DATE__; 
     static const char GMS_RECTANGULAR_WAVEFORM_FP16_BUILD_TIME[]    = __TIME__;
     static const char GMS_RECTANGULAR_WAVEFORM_FP16_SYNOPSIS[]      = "Rectangular waveform generators half_float stored.";

}

#include <cstdint>
#include <string>
#include "GMS_config.h"
#include "GMS_dyn_array.h"

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
// To be added.
#if !defined (RECTANGULAR_WAVEFORM_FP16_USE_NT_STORES)
#define RECTANGULAR_WAVEFORM_FP16_USE_NT_STORES 0
#endif

#ifdef _OPENMP
// Default init a storage arrays for the first-touch (OpenMP) processing
#if !defined(RECTANGULAR_WAVEFORM_FP16_INIT_STORAGE)
#define RECTANGULAR_WAVEFORM_FP16_INIT_STORAGE 1
#endif 
#endif

#if (RECTANGULAR_WAVEFORM_FP16_INIT_STORAGE) == 1
#define INIT_BY_STD_FILL 0
#endif 

// For inlining of trigo functions (asin,acos,sin,cos)
#if !defined(RECTANGULAR_WAVEFORM_FP16_USE_CEPHES)
#define RECTANGULAR_WAVEFORM_FP16_USE_CEPHES 1
#endif 

// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (RECTANGULAR_WAVEFORM_FP16_USE_PMC_INSTRUMENTATION)
#define RECTANGULAR_WAVEFORM_FP16_USE_PMC_INSTRUMENTATION 1
#endif 

#if (RECTANGULAR_WAVEFORM_FP16_USE_PMC_INSTRUMENTATION) == 1
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
           
          struct alignas(64) rectangular_waveform_fp16_t final 
          {
                 std::size_t   __n_samples__;
                 std::uint32_t __n_waves__;
                 float         __T__;   // period
                 float         __rho__; // pulse length
                 float         __A__;   // Amplitude 
                 darray_r2_t   __rw_samples__;

                 rectangular_waveform_fp16_t() = default;

                 rectangular_waveform_fp16_t(const std::size_t,
                                        const std::uint32_t,
                                        const float,
                                        const float,
                                        const float)  noexcept(false);

                  rectangular_waveform_fp16_t(rectangular_waveform_fp16_t &&) noexcept(true);

                  rectangular_waveform_fp16_t(const rectangular_waveform_fp16_t &) noexcept(false);

                  ~rectangular_waveform_fp16_t() noexcept(false);

                  rectangular_waveform_fp16_t & operator=(const rectangular_waveform_fp16_t &) = delete;

                  rectangular_waveform_fp16_t & operator=(rectangular_waveform_fp16_t &&);

                  void init_storage(const half_float::half);

                  static void create_signal_plot(const std::uint32_t,
                                                 const half_float::half * __restrict,
                                                 const half_float::half * __restrict,
                                                 const std::string &,
                                                 const std::string &,
                                                 const bool );

                  void fourier_series_expansion(const std::uint32_t);

                  void fourier_series_expansion_optim(const std::uint32_t);
          };

}

}




















#endif /*__GMS_RECTANGULAR_WAVEFORM_FP16_H__*/