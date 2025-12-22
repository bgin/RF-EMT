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

#ifndef __GMS_IQ_RECTWAVE_BITSTREAM_H__
#define __GMS_IQ_RECTWAVE_BITSTREAM_H__


#include <cstdint>
#include <string>
#include <iostream>
#include "GMS_config.h"
#include "GMS_dyn_array.h"
#include "GMS_sse_memset.h"


// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
// To be added.
#if !defined (IQ_RECTWAVE_BITSTREAM_USE_NT_STORES)
#define IQ_RECTWAVE_BITSTREAM_USE_NT_STORES 0
#endif

#ifdef _OPENMP
// Default init a storage arrays for the first-touch (OpenMP) processing
#if !defined(IQ_RECTWAVE_BITSTREAM_INIT_STORAGE)
#define IQ_RECTWAVE_BITSTREAM_INIT_STORAGE 1
#endif 
#endif

#if (IQ_RECTWAVE_BITSTREAM_INIT_STORAGE) == 1
#define INIT_BY_STD_FILL 0
#endif 

// For inlining of trigo functions (asin,acos,sin,cos)
#if !defined(IQ_RECTWAVE_BITSTREAM_USE_CEPHES)
#define IQ_RECTWAVE_BITSTREAM_USE_CEPHES 1
#endif 

// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (IQ_RECTWAVE_BITSTREAM_USE_PMC_INSTRUMENTATION)
#define IQ_RECTWAVE_BITSTREAM_USE_PMC_INSTRUMENTATION 0
#endif 

#if (IQ_RECTWAVE_BITSTREAM_USE_PMC_INSTRUMENTATION) == 1
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

#if (IQ_RECTWAVE_BITSTREAM_USE_CEPHES) == 0
#include <cmath>
#else 
#include "GMS_cephes_sin_cos.h"
#endif 

namespace gms 
{

namespace radiolocation 
{

              struct alignas(64) iq_rectwave_bitstream_t final 
              {
                     float        m_I_k;        // I-channel number of periods (k*T)
                     float        m_I_T;        // I-channel period
                     std::size_t  m_I_nsamples; // I-channel number of samples (k*T)
                     float        m_I_Ac;       // I-channel carrier Amplitude
                     float        m_I_fc;       // I-channel center frequency
                     float        m_I_ph0;      // I-channel initial phase shift
                     float        m_Q_k;        // Q-channel number of periods (k*T)
                     float        m_Q_T;        // Q-channel period
                     std::size_t  m_Q_nsamples; // Q-channel number of samples (k*T)
                     float        m_Q_Ac;       // Q-channel carrier Amplitude
                     float        m_Q_fc;       // Q-channel center frequency
                     float        m_Q_ph0;      // Q-channel initial phase shift
                     darray_r4_t  m_I_ch_bitstream;
                     darray_r4_t  m_Q_ch_bitstream;

                     iq_rectwave_bitstream_t() noexcept(true);

                     iq_rectwave_bitstream_t(const float, 
                                             const float,
                                             const float,
                                             const float,
                                             const float,
                                             const float,
                                             const float,
                                             const float,
                                             const float,
                                             const float) noexcept(true);

                     iq_rectwave_bitstream_t(const iq_rectwave_bitstream_t &) noexcept(true);

                     iq_rectwave_bitstream_t(iq_rectwave_bitstream_t &&) noexcept(true);

                     ~iq_rectwave_bitstream_t() noexcept(true);

                     iq_rectwave_bitstream_t & operator=(const iq_rectwave_bitstream_t &) noexcept(true);

                     iq_rectwave_bitstream_t & operator=(iq_rectwave_bitstream_t &&) noexcept(true);

                     static void create_signal_plot(   const std::uint32_t,
                                                   const float * __restrict,
                                                   const float * __restrict,
                                                   const std::string &,
                                                   const std::string &,
                                                   const bool ) noexcept(false);

                    __ATTR_ALWAYS_INLINE__
                    inline 
                    void clear_I_bitstream() noexcept(true)
                    {
#if (IQ_RECTWAVE_BITSTREAM_INIT_STORAGE) == 0
                       using namespace gms::common;
	                   sse_memset_unroll8x_ps(this->m_I_ch_bitstream.m_data,0.0f,this->m_I_nsamples);
#else 
                      std::fill(this->m_I_ch_bitstream.m_data,this->m_I_ch_bitstream.m_data+this->m_I_nsamples,0.0f);
#endif
                    }


                    __ATTR_ALWAYS_INLINE__
                    inline
                    void clear_Q_bitstream() noexcept(true)
                    {
#if (IQ_RECTWAVE_BITSTREAM_INIT_STORAGE) == 0
                      using namespace gms::common;
	                  sse_memset_unroll8x_ps(this->m_Q_ch_bitstream.m_data,0.0f,this->m_Q_nsamples);
#else 
                      std::fill(this->m_Q_ch_bitstream.m_data,this->m_Q_ch_bitstream.m_data+this->m_Q_nsamples,0.0f);
#endif                      
                    }

                    __ATTR_HOT__
                    __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
                    __ATTR_OPTIMIZE_03__
#endif     
                     std::int32_t 
                     generate_I_channel_bitstream_scalar();

                    __ATTR_HOT__
                    __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
                    __ATTR_OPTIMIZE_03__
#endif     
                     std::int32_t 
                     generate_Q_channel_bitstream_scalar();

                    __ATTR_HOT__
                    __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
                    __ATTR_OPTIMIZE_03__
#endif     
                     std::int32_t
                     generate_I_channel_bitstream_sse();

                     __ATTR_HOT__
                    __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
                    __ATTR_OPTIMIZE_03__
#endif     
                     std::int32_t
                     generate_Q_channel_bitstream_sse();

                     __ATTR_HOT__
                    __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
                    __ATTR_OPTIMIZE_03__
#endif     
                     std::int32_t
                     generate_I_channel_bitstream_sse_u4x();

                     __ATTR_HOT__
                    __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
                    __ATTR_OPTIMIZE_03__
#endif     
                     std::int32_t
                     generate_Q_channel_bitstream_sse_u4x();

                     __ATTR_HOT__
                    __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
                    __ATTR_OPTIMIZE_03__
#endif     
                     std::int32_t
                     generate_I_channel_bitstream_avx();

                     __ATTR_HOT__
                    __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
                    __ATTR_OPTIMIZE_03__
#endif     
                     std::int32_t
                     generate_Q_channel_bitstream_avx();

                     __ATTR_HOT__
                    __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
                    __ATTR_OPTIMIZE_03__
#endif     
                     std::int32_t
                     generate_I_channel_bitstream_avx_u4x();

                     __ATTR_HOT__
                    __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
                    __ATTR_OPTIMIZE_03__
#endif     
                     std::int32_t
                     generate_Q_channel_bitstream_avx_u4x();

                     __ATTR_HOT__
                    __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
                    __ATTR_OPTIMIZE_03__
#endif     
                     std::int32_t
                     generate_I_channel_bitstream_avx512();

                     __ATTR_HOT__
                    __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
                    __ATTR_OPTIMIZE_03__
#endif     
                     std::int32_t
                     generate_Q_channel_bitstream_avx512();

                     __ATTR_HOT__
                    __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
                    __ATTR_OPTIMIZE_03__
#endif     
                     std::int32_t
                     generate_I_channel_bitstream_avx512_u4x();

                     __ATTR_HOT__
                    __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
                    __ATTR_OPTIMIZE_03__
#endif     
                     std::int32_t
                     generate_Q_channel_bitstream_avx512_u4x();
                                          
              };

              auto operator<<(std::ostream &,
                          const iq_rectwave_bitstream_t &)->std::ostream &;

} //radiolocation

} //gms



















#endif /*__GMS_IQ_RECTWAVE_BITSTREAM_H__*/