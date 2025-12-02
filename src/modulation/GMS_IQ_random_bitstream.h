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

#ifndef __GMS_IQ_RANDOM_BITSTREAM_H__
#define __GMS_IQ_RANDOM_BITSTREAM_H__

namespace file_info 
{

     static const unsigned int GMS_IQ_RANDOM_BITSTREAM_MAJOR = 1;
     static const unsigned int GMS_IQ_RANDOM_BITSTREAM_MINOR = 1;
     static const unsigned int GMS_IQ_RANDOM_BITSTREAM_MICRO = 0;
     static const unsigned int GMS_IQ_RANDOM_BITSTREAM_FULLVER =
       1000U*GMS_IQ_RANDOM_BITSTREAM_MAJOR+100U*GMS_IQ_RANDOM_BITSTREAM_MINOR+
       10U*GMS_IQ_RANDOM_BITSTREAM_MICRO;
     static const char GMS_IQ_RANDOM_BITSTREAM_CREATION_DATE[] = "25-11-2025 06:38AM +00200 (TUE 25 NOV 2025 GMT+2)";
     static const char GMS_IQ_RANDOM_BITSTREAM_BUILD_DATE[]    = __DATE__; 
     static const char GMS_IQ_RANDOM_BITSTREAM_BUILD_TIME[]    = __TIME__;
     static const char GMS_IQ_RANDOM_BITSTREAM_SYNOPSIS[]      = "I/Q Random (unifromely distributed) bitstream (PWM rectangular wave).";

}

#include <cstdint>
#include <string>
#include <iostream>
#include <thread>
#include <valarray>
#include <random>
#include "GMS_config.h"
#include "GMS_dyn_array.h"
#include "GMS_sse_memset.h"

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
// To be added.
#if !defined (IQ_RANDOM_BITSTREAM_USE_NT_STORES)
#define IQ_RANDOM_BITSTREAM_USE_NT_STORES 0
#endif

#ifdef _OPENMP
// Default init a storage arrays for the first-touch (OpenMP) processing
#if !defined(IQ_RANDOM_BITSTREAM_INIT_STORAGE)
#define IQ_RANDOM_BITSTREAM_INIT_STORAGE 1
#endif 
#endif

#if (IQ_RANDOM_BITSTREAM_INIT_STORAGE) == 1
#define INIT_BY_STD_FILL 0
#endif 

// For inlining of trigo functions (asin,acos,sin,cos)
#if !defined(IQ_RANDOM_BITSTREAM_USE_CEPHES)
#define IQ_RANDOM_BITSTREAM_USE_CEPHES 1
#endif 

// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (IQ_RANDOM_BITSTREAM_USE_PMC_INSTRUMENTATION)
#define IQ_RANDOM_BITSTREAM_USE_PMC_INSTRUMENTATION 0
#endif 

#if (IQ_RANDOM_BITSTREAM_USE_PMC_INSTRUMENTATION) == 1
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

#if (IQ_RANDOM_BITSTREAM_USE_CEPHES) == 0
#include <cmath>
#else 
#include "GMS_cephes_sin_cos.h"
#endif 


namespace gms 
{

namespace radiolocation 
{

          struct alignas(64) iq_random_bitstream_t final 
          {
                 float m_I_duration;
                 float m_I_sample_rate;
                 float m_I_w0; // center frequency of underlying sinusoid
                 float m_I_ph0;// initial phase shift
                 std::size_t m_I_nsamples; // inefficient structure packing -- hope that alignas will properly align.
                 float m_Q_duration;
                 float m_Q_sample_rate;
                 float m_Q_w0;
                 float m_Q_ph0;
                 std::size_t m_Q_nsamples;
                 darray_r4_t m_I_bitstream;
                 darray_r4_t m_Q_bitstream;

                 iq_random_bitstream_t() = default;

                 iq_random_bitstream_t(const float,
                                       const float,
                                       const float,
                                       const float,
                                       const float,
                                       const float,
                                       const float,
                                       const float) noexcept(false);

                 iq_random_bitstream_t(const iq_random_bitstream_t &) noexcept(false);

                 iq_random_bitstream_t(iq_random_bitstream_t &&) noexcept(true);

                 ~iq_random_bitstream_t() noexcept(false);

                 iq_random_bitstream_t & 
                                    operator=(const iq_random_bitstream_t &) noexcept(false);

                 iq_random_bitstream_t &
                                    operator=(iq_random_bitstream_t &&) noexcept(true);

                 static void create_signal_plot( const std::uint32_t,
                                                   const float * __restrict,
                                                   const float * __restrict,
                                                   const std::string &,
                                                   const std::string &,
                                                   const bool ) noexcept(false);

                  __ATTR_ALWAYS_INLINE__
                  inline 
                  void clear_I_bitstream() noexcept(true)
                  {
#if (IQ_RANDOM_BITSTREAM_INIT_STORAGE) == 0
                      using namespace gms::common;
	                  sse_memset_unroll8x_ps(this->m_I_bitstream.m_data,0.0f,this->m_I_nsamples);
#else 
                      std::fill(this->m_I_bitstream.m_data,this->m_I_bitstream.m_data+this->m_I_nsamples,0.0f);
#endif
                  }


                  __ATTR_ALWAYS_INLINE__
                  inline
                  void clear_Q_bitstream() noexcept(true)
                  {
#if (IQ_RANDOM_BITSTREAM_INIT_STORAGE) == 0
                      using namespace gms::common;
	                  sse_memset_unroll8x_ps(this->m_Q_bitstream.m_data,0.0f,this->m_Q_nsamples);
#else 
                      std::fill(this->m_Q_bitstream.m_data,this->m_Q_bitstream.m_data+this->m_Q_nsamples,0.0f);
#endif                      
                  }

                   __ATTR_ALWAYS_INLINE__
                  inline          
                  float 
                  channel_I_uniform_sample_n1_0_1(const float lo, // either -1 or 0
                                                  const float hi) // always 1
                  {
                      thread_local unsigned long long I_seed{0ull};
                      thread_local std::uniform_real_distribution<float> rand_I_bits;
                      I_seed      = __rdtsc();
                      rand_I_bits = std::uniform_real_distribution<float>(lo,hi);
                      auto rseed_I{std::mt19937(I_seed)};
                      float uni_noise_sample = rand_I_bits.operator()(rseed_I);
                      return (uni_noise_sample);
                  }

                  __ATTR_ALWAYS_INLINE__
                  inline          
                  float 
                  channel_Q_uniform_sample_n1_0_1(const float lo, // either -1 or 0
                                                  const float hi) // always 1
                  {
                      thread_local unsigned long long Q_seed{0ull};
                      thread_local std::uniform_real_distribution<float> rand_Q_bits;
                      Q_seed      = __rdtsc();
                      rand_Q_bits = std::uniform_real_distribution<float>(lo,hi);
                      auto rseed_Q{std::mt19937(Q_seed)};
                      float uni_noise_sample = rand_Q_bits.operator()(rseed_Q);
                      return (uni_noise_sample);
                  }

                 
                  std::int32_t generate_channel_I_random_bitstream(const float lo,
                                                                   const float hi,
                                                                   const std::int32_t generate_random) 
                  {
                       using namespace gms::math;
                       constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
                       const float inv_sample_rate{1.0f/this->m_I_sample_rate};
                       switch (generate_random)
                       {
                          case 0 : // non-random
                              for(std::size_t i = 0ull; i != this->m_I_nsamples; ++i) 
                              {
                                    
                                    const float t_i{static_cast<float>(i*inv_sample_rate)};
#if (IQ_RANDOM_BITSTREAM_USE_CEPHES) == 1 
                                    const float cos_val = ceph_cosf(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_w0*t_i));
                                    this->m_I_bitstream.m_data[i] = (cos_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
                                    const float cos_val = std::cos(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_w0*t_i));
                                    this->m_I_bitstream.m_data[i] = (cos_val>0.0f)?+1.0f:-1.0f;
#endif
                              }
                           break;
                           case 1 : // random
                              for(std::size_t i = 0ull; i != this->m_I_nsamples; ++i) 
                              {
                                    const float cos_r_t{channel_Q_uniform_sample_n1_0_1(lo,hi)};
                                    const float cos_r_f{channel_Q_uniform_sample_n1_0_1(lo,hi)};
                                    const float t_i{static_cast<float>(i*inv_sample_rate)};
#if (IQ_RANDOM_BITSTREAM_USE_CEPHES) == 1 
                                    const float cos_val = ceph_cosf(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_w0*t_i));
                                    this->m_I_bitstream.m_data[i] = (cos_val>=0.0f)?cos_r_t:cos_r_f;
                                    
#else 
                                    const float cos_val = std::cos(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_w0*t_i));
                                    this->m_I_bitstream.m_data[i] = (cos_val>0.0f)?cos_r_t:cos_r_f;
#endif
                              }
                           break;
                           default :
                              return (-1);
                       }
                                          
                       return (0);
                  } 

                                   
                  std::int32_t generate_channel_Q_random_bitstream(const float lo, 
                                                                   const float hi,
                                                                   const std::int32_t generate_random) 
                  {
                       using namespace gms::math;
                       constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
                       const float inv_sample_rate{1.0f/this->m_Q_sample_rate};  
                       switch (generate_random)
                       {     
                          case 0 : // non-random
                              for(std::size_t i = 0ull; i != this->m_Q_nsamples; ++i) 
                              {
                                    const float t_i{static_cast<float>(i*inv_sample_rate)};
#if (IQ_RANDOM_BITSTREAM_USE_CEPHES) == 1 
                                    const float sin_val = ceph_sinf(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_w0*t_i));
                                    this->m_Q_bitstream.m_data[i] = (sin_val>0.0f)?+1.0f:-1.0f;
                                    
#else 
                                    const float sin_val = std::sin(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_w0*t_i));
                                    this->m_Q_bitstream.m_data[i] = (sin_val>0.0f)?+1.0f:-1.0f;
#endif
                              }
                          break;
                          case 1 : // random 
                              for(std::size_t i = 0ull; i != this->m_Q_nsamples; ++i) 
                              {
                                    const float sin_r_t{channel_Q_uniform_sample_n1_0_1(lo,hi)};
                                    const float sin_r_f{channel_Q_uniform_sample_n1_0_1(lo,hi)};
                                    const float t_i{static_cast<float>(i*inv_sample_rate)};
#if (IQ_RANDOM_BITSTREAM_USE_CEPHES) == 1 
                                    const float sin_val = ceph_sinf(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_w0*t_i));
                                    this->m_Q_bitstream.m_data[i] = (sin_val>0.0f)?sin_r_t:sin_r_f;
                                    
#else 
                                    const float sin_val = std::sin(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_w0*t_i));
                                    this->m_Q_bitstream.m_data[i] = (sin_val>0.0f)?sin_r_t:sin_r_f;
#endif
                              }
                          break;
                          default : 
                              return (-1);
                       }
                            return (0);
                       }

                       
                    

          };

          auto operator<<(std::ostream &,
                          const iq_random_bitstream_t &)->std::ostream &;

}

}



#endif /*__GMS_IQ_RANDOM_BITSTREAM_H__*/