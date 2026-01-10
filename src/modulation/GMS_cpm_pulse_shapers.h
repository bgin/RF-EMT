
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

#ifndef __GMS_CPM_PULSE_SHAPERS_H__
#define __GMS_CPM_PULSE_SHAPERS_H__ 100120260718


namespace file_info 
{

     static const unsigned int GMS_CPM_PULSE_SHAPERS_MAJOR = 1;
     static const unsigned int GMS_CPM_PULSE_SHAPERS_MINOR = 1;
     static const unsigned int GMS_CPM_PULSE_SHAPERS_MICRO = 0;
     static const unsigned int GMS_CPM_PULSE_SHAPERS_FULLVER =
       1000U*GMS_CPM_PULSE_SHAPERS_MAJOR+100U*GMS_CPM_PULSE_SHAPERS_MINOR+
       10U*GMS_CPM_PULSE_SHAPERS_MICRO;
     static const char GMS_CPM_PULSE_SHAPERS_CREATION_DATE[] = "10-01-2025 07:18PM +00200 (SAT 10 JAN 2026 GMT+2)";
     static const char GMS_CPM_PULSE_SHAPERS_BUILD_DATE[]    = __DATE__; 
     static const char GMS_CPM_PULSE_SHAPERS_BUILD_TIME[]    = __TIME__;
     static const char GMS_CPM_PULSE_SHAPERS_AUTHOR[]        =  "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
     static const char GMS_CPM_PULSE_SHAPERS_SYNOPSIS[]      = "Continous Phase Modulation Pulse Shapers.";

}

#include <cstdint>
#include <string>
#include <iostream>
#include "GMS_config.h"
#include "GMS_dyn_array.h"
#include "GMS_sse_memset.h"
#include "GMS_avx_memset.h"
#include "GMS_avx512_memset.h"

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
// To be added.
#if !defined (CPM_PULSE_SHAPERS_NT_STORES)
#define CPM_PULSE_SHAPERS_USE_NT_STORES 0
#endif

#ifdef _OPENMP
// Default init a storage arrays for the first-touch (OpenMP) processing
#if !defined(CPM_PULSE_SHAPERS_INIT_STORAGE)
#define CPM_PULSE_SHAPERS_INIT_STORAGE 1
#endif 
#endif

#if (CPM_PULSE_SHAPERS_INIT_STORAGE) == 1
#define INIT_BY_STD_FILL 0
#endif 

// For inlining of trigo functions (asin,acos,sin,cos)
#if !defined(CPM_PULSE_SHAPERS_USE_CEPHES)
#define CPM_PULSE_SHAPERS_USE_CEPHES 1
#endif 

// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (CPM_PULSE_SHAPERS_USE_PMC_INSTRUMENTATION)
#define CPM_PULSE_SHAPERS_USE_PMC_INSTRUMENTATION 0
#endif 

#if (CPM_PULSE_SHAPERS_USE_PMC_INSTRUMENTATION) == 1
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

#if (CPM_PULSE_SHAPERS_USE_CEPHES) == 0
#include <cmath>
#else 
#include "GMS_cephes_sin_cos.h"
#endif 

namespace gms 
{

namespace radiolocation 
{

          struct alignas(64) cpm_pulse_shapers_t final 
          {
                 std::size_t    m_nTLsamples; //T*L samples
                 std::size_t    m_nTsamples; //T samples
                 float          m_T; //sample period length
                 float          m_L; //symbols length
                 float          m_beta; // LSRC argument
                 float          m_BbT;   // GMSK argument
                 darray_r4_t    m_lrec_pulse; //LREC samples
                 darray_r4_t    m_lrc_pulse;  //L(Raised-Cosine)
                 darray_r4_t    m_lsrc_pulse; //L(Spectrally-Raised-Cosine)
                 darray_r4_t    m_tfm_pulse;  //Tamed-Frequency-Modulation
                 darray_r4_t    m_gmsk_pulse; // Gaussian-MSK

                 cpm_pulse_shapers_t() noexcept(true);

                 cpm_pulse_shapers_t(const std::size_t,
                                     const std::size_t,
                                     const float,
                                     const float,
                                     const float,
                                     const float) noexcept(false);

                 cpm_pulse_shapers_t(const cpm_pulse_shapers_t &) noexcept(false);

                 cpm_pulse_shapers_t(cpm_pulse_shapers_t &&) noexcept(true);

                 ~cpm_pulse_shapers_t() noexcept(true);

                 cpm_pulse_shapers_t & operator=(const cpm_pulse_shapers_t &) noexcept(false);

                 cpm_pulse_shapers_t & operator=(cpm_pulse_shapers_t &&) noexcept(true);

                  __ATTR_ALWAYS_INLINE__ 
                  inline void clear_lrec_pulse(const float filler) noexcept(true)
                  {
#if (INIT_BY_STD_FILL) == 0
                      using namespace gms::common;
	                  sse_memset_unroll8x_ps(this->m_lrec_pulse.m_data,filler,this->m_nTLsamples);
#else 
                      std::fill(this->m_lrec_pulse.m_data,this->m_lrec_pulse.m_data+this->m_nTLsamples,filler);
#endif   
                   }

                    __ATTR_ALWAYS_INLINE__ 
                  inline void clear_lrc_pulse(const float filler) noexcept(true)
                  {
#if (INIT_BY_STD_FILL) == 0
                      using namespace gms::common;
	                  sse_memset_unroll8x_ps(this->m_lrc_pulse.m_data,filler,this->m_nTLsamples);
#else 
                      std::fill(this->m_lrc_pulse.m_data,this->m_lrc_pulse.m_data+this->m_nTLsamples,filler);
#endif   
                   }

                    __ATTR_ALWAYS_INLINE__ 
                  inline void clear_lsrc_pulse(const float filler) noexcept(true)
                  {
#if (INIT_BY_STD_FILL) == 0
                      using namespace gms::common;
	                  sse_memset_unroll8x_ps(this->m_lsrc_pulse.m_data,filler,this->m_nTLsamples);
#else 
                      std::fill(this->m_lsrc_pulse.m_data,this->m_lsrc_pulse.m_data+this->m_nTLsamples,filler);
#endif   
                   }

                    __ATTR_ALWAYS_INLINE__ 
                  inline void clear_tfm_pulse(const float filler) noexcept(true)
                  {
#if (INIT_BY_STD_FILL) == 0
                      using namespace gms::common;
	                  sse_memset_unroll8x_ps(this->m_tfm_pulse.m_data,filler,this->m_nTsamples);
#else 
                      std::fill(this->m_tfm_pulse.m_data,this->m_tfm_pulse.m_data+this->m_nTsamples,filler);
#endif   
                   }

                    __ATTR_ALWAYS_INLINE__ 
                  inline void clear_gms_pulse(const float filler) noexcept(true)
                  {
#if (INIT_BY_STD_FILL) == 0
                      using namespace gms::common;
	                  sse_memset_unroll8x_ps(this->m_gmsk_pulse.m_data,filler,this->m_nTsamples);
#else 
                      std::fill(this->m_gmsk_pulse.m_data,this->m_gmsk_pulse.m_data+this->m_nTsamples,filler);
#endif   
                   }

                   static void create_signal_plot( const std::uint32_t,
                                                   const float * __restrict,
                                                   const float * __restrict,
                                                   const std::string &,
                                                   const std::string &,
                                                   const bool );

                   // pulse shapers
                   
                    __ATTR_ALWAYS_INLINE__ 
                   inline std::int32_t 
                   generate_lrec_pulse_scalar_u8x() 
                   {
                       std::size_t i,j;
                       const float twoLT{2.0f*this->m_L*this->m_T};
                       const float value{1.0f/twoLT};
                       for(i = 0ull; i != ROUND_TO_EIGHT(this->m_nTLsamples,8ull)l i += 8ull) 
                       {
                           this->m_lrec_pulse.m_data[i+0ull] = value;
                           this->m_lrec_pulse.m_data[i+1ull] = value;
                           this->m_lrec_pulse.m_data[i+2ull] = value;
                           this->m_lrec_pulse.m_data[i+3ull] = value;
                           this->m_lrec_pulse.m_data[i+4ull] = value;
                           this->m_lrec_pulse.m_data[i+5ull] = value;
                           this->m_lrec_pulse.m_data[i+6ull] = value;
                           this->m_lrec_pulse.m_data[i+7ull] = value;
                       }

                       for(j = i; j != this->m_nTLsamples; ++j)  
                       {
                           this->m_lrec_pulse.m_data[j] = value;
                       }
                       return (0);
                   }

                   __ATTR_ALWAYS_INLINE__ 
                   inline std::int32_t 
                   generate_lrec_pulse_sse_u8x() 
                   {
                       using namespace gms::common;
                       const float twoLT{2.0f*this->m_L*this->m_T};
                       const float value{1.0f/twoLT};
                       sse_memset_unroll8x_ps(&this->m_lrec_pulse.m_data[0],value,this->m_nTLsamples);
                       return (0);
                   }

                   __ATTR_ALWAYS_INLINE__ 
                   inline std::int32_t 
                   generate_lrec_pulse_avx_u8x() 
                   {
                       using namespace gms::common;
                       const float twoLT{2.0f*this->m_L*this->m_T};
                       const float value{1.0f/twoLT};
                       avx_memset_unroll8x_ps(&this->m_lrec_pulse.m_data[0],value,this->m_nTLsamples);
                       return (0);
                   }

                   __ATTR_ALWAYS_INLINE__ 
                   inline std::int32_t 
                   generate_lrec_pulse_avx512_u8x() 
                   {
                       using namespace gms::common;
                       const float twoLT{2.0f*this->m_L*this->m_T};
                       const float value{1.0f/twoLT};
                       avx512_memset_unroll8x_ps(&this->m_lrec_pulse.m_data[0],value,this->m_nTLsamples);
                       return (0);
                   }

                   
                   
          };    

          auto 
          operator<<(std::ostream &, const cpm_pulse_shapers_t &)->std::ostream &;     

} //radiolocation

}// gms

















#endif /*__GMS_CPM_PULSE_SHAPERS_H__*/