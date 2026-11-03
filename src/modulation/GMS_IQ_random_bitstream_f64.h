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

#ifndef __GMS_IQ_RANDOM_BITSTREAM_F64_H__
#define __GMS_IQ_RANDOM_BITSTREAM_F64_H__ 031120260719

namespace file_info 
{

     static const unsigned int GMS_IQ_RANDOM_BITSTREAM_F64_MAJOR = 1;
     static const unsigned int GMS_IQ_RANDOM_BITSTREAM_F64_MINOR = 0;
     static const unsigned int GMS_IQ_RANDOM_BITSTREAM_F64_MICRO = 0;
     static const unsigned int GMS_IQ_RANDOM_BITSTREAM_F64_FULLVER =
       1000U*GMS_IQ_RANDOM_BITSTREAM_F64_MAJOR+100U*GMS_IQ_RANDOM_BITSTREAM_F64_MINOR+
       10U*GMS_IQ_RANDOM_BITSTREAM_F64_MICRO;
     static const char GMS_IQ_RANDOM_BITSTREAM_F64_CREATION_DATE[] = "03-11-2026 07:19AM +00200 (TUE 03 NOV 2026 GMT+2)";
     static const char GMS_IQ_RANDOM_BITSTREAM_F64_BUILD_DATE[]    = __DATE__; 
     static const char GMS_IQ_RANDOM_BITSTREAM_F64_BUILD_TIME[]    = __TIME__;
     static const char GMS_IQ_RANDOM_BITSTREAM_F64_AUTHOR[]        = "Programmer: Bernard Gingold, beniekg@gmail.com"; 
     static const char GMS_IQ_RANDOM_BITSTREAM_F64_SYNOPSIS[]      = "I/Q Random (uniformely distributed) bitstream (PWM rectangular wave) double precision.";

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
#if !defined (IQ_RANDOM_BITSTREAM_F64_USE_NT_STORES)
#define IQ_RANDOM_BITSTREAM_F64_USE_NT_STORES 0
#endif

#ifdef _OPENMP
// Default init a storage arrays for the first-touch (OpenMP) processing
#if !defined(IQ_RANDOM_BITSTREAM_F64_INIT_STORAGE)
#define IQ_RANDOM_BITSTREAM_F64_INIT_STORAGE 1
#endif 
#endif

#if (IQ_RANDOM_BITSTREAM_F64_INIT_STORAGE) == 1
#define INIT_BY_STD_FILL 0
#endif 

// For inlining of trigo functions (asin,acos,sin,cos)
#if !defined(IQ_RANDOM_BITSTREAM_F64_USE_CEPHES)
#define IQ_RANDOM_BITSTREAM_F64_USE_CEPHES 1
#endif 

// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (IQ_RANDOM_BITSTREAM_F64_USE_PMC_INSTRUMENTATION)
#define IQ_RANDOM_BITSTREAM_F64_USE_PMC_INSTRUMENTATION 0
#endif 

#if (IQ_RANDOM_BITSTREAM_F64_USE_PMC_INSTRUMENTATION) == 1
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

#if (IQ_RANDOM_BITSTREAM_F64_USE_CEPHES) == 0
#include <cmath>
#else 
#include "GMS_cephes_double.h"
#endif 


namespace gms 
{

namespace radiolocation 
{

          struct alignas(64) iq_random_bitstream_f64_t final 
          {
                 double m_I_duration;
                 double m_I_sample_rate;
                 double m_I_w0; // center frequency of underlying sinusoid
                 double m_I_ph0;// initial phase shift
                 std::size_t m_I_nsamples; // inefficient structure packing -- hope that alignas will properly align.
                 double m_Q_duration;
                 double m_Q_sample_rate;
                 double m_Q_w0;
                 double m_Q_ph0;
                 std::size_t m_Q_nsamples;
                 darray_r8_t m_I_bitstream;
                 darray_r8_t m_Q_bitstream;

                 iq_random_bitstream_f64_t() = default;

                 iq_random_bitstream_f64_t(const double,
                                       const double,
                                       const double,
                                       const double,
                                       const double,
                                       const double,
                                       const double,
                                       const double) noexcept(false);

                 iq_random_bitstream_f64_t(const iq_random_bitstream_f64_t &) noexcept(false);

                 iq_random_bitstream_f64_t(iq_random_bitstream_f64_t &&) noexcept(true);

                 ~iq_random_bitstream_f64_t() noexcept(false);

                 iq_random_bitstream_f64_t & 
                                    operator=(const iq_random_bitstream_f64_t &) noexcept(false);

                 iq_random_bitstream_f64_t &
                                    operator=(iq_random_bitstream_f64_t &&) noexcept(true);

                 static void create_signal_plot(   const std::uint32_t,
                                                   const double * __restrict,
                                                   const double * __restrict,
                                                   const std::string &,
                                                   const std::string &,
                                                   const bool ) noexcept(false);

                  __ATTR_ALWAYS_INLINE__
                  inline 
                  void clear_I_bitstream() noexcept(true)
                  {
#if (IQ_RANDOM_BITSTREAM_F64_INIT_STORAGE) == 0
                      using namespace gms::common;
	                  sse_memset_unroll8x_pd(this->m_I_bitstream.m_data,0.0f,this->m_I_nsamples);
#else 
                      std::fill(this->m_I_bitstream.m_data,this->m_I_bitstream.m_data+this->m_I_nsamples,0.0f);
#endif
                  }


                  __ATTR_ALWAYS_INLINE__
                  inline
                  void clear_Q_bitstream() noexcept(true)
                  {
#if (IQ_RANDOM_BITSTREAM_F64_INIT_STORAGE) == 0
                      using namespace gms::common;
	                  sse_memset_unroll8x_pd(this->m_Q_bitstream.m_data,0.0f,this->m_Q_nsamples);
#else 
                      std::fill(this->m_Q_bitstream.m_data,this->m_Q_bitstream.m_data+this->m_Q_nsamples,0.0f);
#endif                      
                  }

                   __ATTR_ALWAYS_INLINE__
                  inline          
                  double 
                  channel_I_uniform_sample_n1_0_1(const double lo, // either -1 or 0
                                                  const double hi) // always 1
                  {
                      thread_local unsigned long long I_seed{0ull};
                      thread_local std::uniform_real_distribution<double> rand_I_bits;
                      I_seed      = __rdtsc();
                      rand_I_bits = std::uniform_real_distribution<double>(lo,hi);
                      auto rseed_I{std::mt19937(I_seed)};
                      double uni_noise_sample = rand_I_bits.operator()(rseed_I);
                      return (uni_noise_sample);
                  }

                  __ATTR_ALWAYS_INLINE__
                  inline          
                  double 
                  channel_Q_uniform_sample_n1_0_1(const double lo, // either -1 or 0
                                                  const double hi) // always 1
                  {
                      thread_local unsigned long long Q_seed{0ull};
                      thread_local std::uniform_real_distribution<double> rand_Q_bits;
                      Q_seed      = __rdtsc();
                      rand_Q_bits = std::uniform_real_distribution<double>(lo,hi);
                      auto rseed_Q{std::mt19937(Q_seed)};
                      double uni_noise_sample = rand_Q_bits.operator()(rseed_Q);
                      return (uni_noise_sample);
                  }

                 
                  std::int32_t generate_channel_I_random_bitstream(const double lo,
                                                                   const double hi,
                                                                   const std::int32_t generate_random) 
                  {
                       using namespace gms::math;
                       constexpr double C6283185307179586476925286766559{6.283185307179586476925286766559};
                       const double inv_sample_rate{1.0f/this->m_I_sample_rate};
                       switch (generate_random)
                       {
                          case 0 : // non-random
                              for(std::size_t i = 0ull; i != this->m_I_nsamples; ++i) 
                              {
                                    
                                    const double t_i{static_cast<double>(i*inv_sample_rate)};
#if (IQ_RANDOM_BITSTREAM_F64_USE_CEPHES) == 1 
                                    const double cos_val = cephes_d::cos(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_w0*t_i));
                                    this->m_I_bitstream.m_data[i] = (cos_val>=0.0)?+1.0:-1.0;
                                    
#else 
                                    const double cos_val = std::cos(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_w0*t_i));
                                    this->m_I_bitstream.m_data[i] = (cos_val>0.0)?+1.0:-1.0;
#endif
                              }
                           break;
                           case 1 : // random
                              for(std::size_t i = 0ull; i != this->m_I_nsamples; ++i) 
                              {
                                    const double cos_r_t{channel_Q_uniform_sample_n1_0_1(lo,hi)};
                                    const double cos_r_f{channel_Q_uniform_sample_n1_0_1(lo,hi)};
                                    const double t_i{static_cast<double>(i*inv_sample_rate)};
#if (IQ_RANDOM_BITSTREAM_F64_USE_CEPHES) == 1 
                                    const double cos_val = cephes_d::cos(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_w0*t_i));
                                    this->m_I_bitstream.m_data[i] = (cos_val>=0.0)?cos_r_t:cos_r_f;
                                    
#else 
                                    const double cos_val = std::cos(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_w0*t_i));
                                    this->m_I_bitstream.m_data[i] = (cos_val>0.0)?cos_r_t:cos_r_f;
#endif
                              }
                           break;
                           default :
                              return (-1);
                       }
                                          
                       return (0);
                  } 

                                   
                  std::int32_t generate_channel_Q_random_bitstream(const double lo, 
                                                                   const double hi,
                                                                   const std::int32_t generate_random) 
                  {
                       using namespace gms::math;
                       constexpr double C6283185307179586476925286766559{6.283185307179586476925286766559};
                       const double inv_sample_rate{1.0f/this->m_Q_sample_rate};  
                       switch (generate_random)
                       {     
                          case 0 : // non-random
                              for(std::size_t i = 0ull; i != this->m_Q_nsamples; ++i) 
                              {
                                    const double t_i{static_cast<double>(i*inv_sample_rate)};
#if (IQ_RANDOM_BITSTREAM_F64_USE_CEPHES) == 1 
                                    const double sin_val = cephes_d::sin(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_w0*t_i));
                                    this->m_Q_bitstream.m_data[i] = (sin_val>0.0)?+1.0:-1.0;
                                    
#else 
                                    const double sin_val = std::sin(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_w0*t_i));
                                    this->m_Q_bitstream.m_data[i] = (sin_val>0.0f)?+1.0f:-1.0f;
#endif
                              }
                          break;
                          case 1 : // random 
                              for(std::size_t i = 0ull; i != this->m_Q_nsamples; ++i) 
                              {
                                    const double sin_r_t{channel_Q_uniform_sample_n1_0_1(lo,hi)};
                                    const double sin_r_f{channel_Q_uniform_sample_n1_0_1(lo,hi)};
                                    const double t_i{static_cast<double>(i*inv_sample_rate)};
#if (IQ_RANDOM_BITSTREAM_F64_USE_CEPHES) == 1 
                                    const double sin_val = cephes_d::sin(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_w0*t_i));
                                    this->m_Q_bitstream.m_data[i] = (sin_val>0.0)?sin_r_t:sin_r_f;
                                    
#else 
                                    const double sin_val = std::sin(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_w0*t_i));
                                    this->m_Q_bitstream.m_data[i] = (sin_val>0.0)?sin_r_t:sin_r_f;
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
                          const iq_random_bitstream_f64_t &)->std::ostream &;

}

}

#endif /*__GMS_IQ_RANDOM_BITSTREAM_F64_H__*/