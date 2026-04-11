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

#ifndef __GMS_CPFSK_SIGNAL_H__
#define __GMS_CPFSK_SIGNAL_H__ 040420260943

namespace file_info 
{

     static const unsigned int GMS_CPFSK_SIGNAL_MAJOR = 1;
     static const unsigned int GMS_CPFSK_SIGNAL_MINOR = 1;
     static const unsigned int GMS_CPFSK_SIGNAL_MICRO = 0;
     static const unsigned int GMS_CPFSK_SIGNAL_FULLVER =
       1000U*GMS_CPFSK_SIGNAL_MAJOR+100U*GMS_CPFSK_SIGNAL_MINOR+
       10U*GMS_CPFSK_SIGNAL_MICRO;
     static const char GMS_CPFSK_SIGNAL_CREATION_DATE[] = "04-04-2026 09:42PM +00200 (SAT 04 APR 2026 GMT+2)";
     static const char GMS_CPFSK_SIGNAL_BUILD_DATE[]    = __DATE__; 
     static const char GMS_CPFSK_SIGNAL_BUILD_TIME[]    = __TIME__;
     static const char GMS_CPFSK_SIGNAL_AUTHOR[]        =  "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
     static const char GMS_CPFSK_SIGNAL_SYNOPSIS[]      =  "Continous Phase Frequency Shift Keying (CPFSK) signal.";

}

#include <cstdint>
#include <string>
#include <iostream>
#include "GMS_config.h"
#include "GMS_dyn_array.h"
#include "GMS_sse_memset.h"

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
// To be added.
#if !defined (CPFSK_SIGNAL_USE_NT_STORES)
#define CPFSK_SIGNAL_USE_NT_STORES 0
#endif

#ifdef _OPENMP
// Default init a storage arrays for the first-touch (OpenMP) processing
#if !defined(CPFSK_SIGNAL_INIT_STORAGE)
#define CPFSK_SIGNAL_INIT_STORAGE 1
#endif 
#endif

#if (CPFSK_SIGNAL_INIT_STORAGE) == 1
#define INIT_BY_STD_FILL 0
#endif 

// For inlining of trigo functions (asin,acos,sin,cos)
#if !defined(CPFSK_SIGNAL_USE_CEPHES)
#define CPFSK_SIGNAL_USE_CEPHES 1
#endif 

// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (CPFSK_SIGNAL_USE_PMC_INSTRUMENTATION)
#define CPFSK_SIGNAL_USE_PMC_INSTRUMENTATION 0
#endif 

#if (CPFSK_SIGNAL_USE_PMC_INSTRUMENTATION) == 1
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

#if (CPFSK_SIGNAL_USE_CEPHES) == 0
#include <cmath>
#else 
#include "GMS_cephes_sin_cos.h"
#endif 

namespace gms
{

namespace radiolocation
{

           struct alignas(64) cpfsk_signal_t final 
           {

                  std::size_t   m_T; // T samples
                  std::size_t   m_nfreqs; // number of PSD frequencies
                  std::uint32_t m_k; // k-th pulse
                  float         m_A; // An amplitude
                  float         m_h; // h-modulation index.
                  float         m_fc; // carrier frequency
                  float         m_ph0; // initial phase offset
                  float         m_Phik;// phase offset coefficient
                  darray_r4_t   m_cpfsk_signal;
                  darray_r4_t   m_cpfsk_psd; // CPFSK Power Spectral Density (PSD) only M-ary 1REC type.

                  enum class User_data_bitstream_type : std::int32_t 
                  {
                       IQ_random_bitstream_type,
                       IQ_rectw_bitstream_type,
                       IQ_rectw_bitstream_vseq_type
                  };

                  inline cpfsk_signal_t() noexcept(true);

                  inline cpfsk_signal_t(const std::size_t,
                                        const std::size_t,
                                        const std::uint32_t,
                                        const float,
                                        const float,
                                        const float,
                                        const float) noexcept(false);

                  inline cpfsk_signal_t(const cpfsk_signal_t &) noexcept(false);

                  inline cpfsk_signal_t(cpfsk_signal_t &&) noexcept(true);

                  inline ~cpfsk_signal_t() noexcept(true);

                  inline cpfsk_signal_t & 
                  operator=(const cpfsk_signal_t &) noexcept(false);

                  inline cpfsk_signal_t 
                  operator=(cpfsk_signal_t &&) noexcept(true);

                  static void create_sequence_plot(const std::uint32_t,
                                                   const float * __restrict,
                                                   const float * __restrict,
                                                   const std::string &,
                                                   const std::string &,
                                                   const bool );

                 
                  __ATTR_ALWAYS_INLINE__ 
                  inline void clear_cpfsk_signal(const float filler) noexcept(true)
                  {
#if (INIT_BY_STD_FILL) == 0
                      using namespace gms::common;
	                 sse_memset_unroll8x_ps(this->m_cpfsk_signal.m_data,filler,this->m_T);
#else 
                      std::fill(this->m_cpfsk_signal.m_data,this->m_cpfsk_signal.m_data+this->m_T,filler);
#endif   
                   }

                   __ATTR_ALWAYS_INLINE__ 
                  inline void clear_cpfsk_psd(const float filler) noexcept(true)
                  {
#if (INIT_BY_STD_FILL) == 0
                      using namespace gms::common;
	                 sse_memset_unroll8x_ps(this->m_cpfsk_psd.m_data,filler,this->m_nfreqs);
#else 
                      std::fill(this->m_cpfsk_pds.m_data,this->m_cpfsk_pds.m_data+this->m_nfreqs,filler);
#endif   
                   }

                   /*Must be called before the CPFSK generators!!*/
                   __ATTR_ALWAYS_INLINE__ 
                   inline std::int32_t compute_Phik(const darray_r4_t &user_data) noexcept(true) 
                   {
                          if(__builtin_expect(static_cast<std::size_t>(this->m_k)!=user_data.mnx,0)) 
                          { 
                             return (-1);
                          }
                          constexpr float PI{3.14159265358979323846264338328f};
                          const float * __restrict__ p_user_data{user_data.m_data};
                          const float ak(p_user_data[this->m_k-1ull]);
                          constexpr std::int32_t mod_2PI{360};
                          float sum_f{0.0f};
                          std::int32_t sum_i{};
                          std::int32_t mod_term{};
                          for(std::size_t i{0ull}; i != user_data.mnx; ++i) 
                          {
                              const float t_k(static_cast<float>(i));
                              const float ai(user_data.m_data[i]);
                              const float diff(ai-t_k*ak);
                              sum_f += diff;
                          }
                          sum_i = static_cast<std::int32_t>(sum_f);
                          mod_term = sum_i%mod_2PI;
                          this->m_Phik = PI*this->m_h*static_cast<float>(mod_term);
                          return (0);
                   }


__ATTR_HOT__
__ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
__ATTR_OPTIMIZE_03__
#endif                    
                  std::int32_t 
                  compute_cpfsk_psd_analytic(const darray_r4_t &,
                                             const std::int32_t) noexcept(true);
    
 __ATTR_HOT__
__ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
__ATTR_OPTIMIZE_03__
#endif                                  
                   std::int32_t 
                   generate_cpfsk_signal_scalar_u8x(const float) noexcept(true);

 __ATTR_HOT__
__ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
__ATTR_OPTIMIZE_03__
#endif                                   
                   std::int32_t 
                   generate_cpfsk_signal_sse(const float) noexcept(true);

__ATTR_HOT__
__ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
__ATTR_OPTIMIZE_03__
#endif                                   
                   std::int32_t 
                   generate_cpfsk_signal_sse_u8x(const float) noexcept(true);
                                                 

__ATTR_HOT__
__ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
__ATTR_OPTIMIZE_03__
#endif                                   
                   std::int32_t 
                   generate_cpfsk_signal_avx(const float,
                                             const bool) noexcept(true);

__ATTR_HOT__
__ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
__ATTR_OPTIMIZE_03__
#endif                                   
                   std::int32_t 
                   generate_cpfsk_signal_avx_u8x(const float,
                                                 const bool) noexcept(true);

__ATTR_HOT__
__ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
__ATTR_OPTIMIZE_03__
#endif                                   
                   std::int32_t 
                   generate_cpfsk_signal_avx512(const float,
                                                const bool) noexcept(true);

__ATTR_HOT__
__ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
__ATTR_OPTIMIZE_03__
#endif                                   
                   std::int32_t 
                   generate_cpfsk_signal_avx512_u8x(const float,
                                                    const bool) noexcept(true);

           };

          auto 
          operator<<(std::ostream &, const cpfsk_signal_t &)->std::ostream &;

} //radiolocation

} //gms































#endif /*__GMS_CPFSK_SIGNAL_H__*/