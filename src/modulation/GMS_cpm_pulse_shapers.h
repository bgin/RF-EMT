
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
#include "GMS_loop_indices_LUT.h"

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

#if !defined(CPM_PULSE_SHAPERS_USE_VEC_RCP14_PRECISION)
#define CPM_PULSE_SHAPERS_USE_VEC_RCP14_PRECISION 1
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

                   __ATTR_ALWAYS_INLINE__
                   inline std::int32_t 
                   generate_lrc_pulse_scalar_u8x()
                   {
                       using namespace gms::math;
                       constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
                       std::size_t i,j;
                       const float twoLT{2.0f*this->m_L*this->m_T};
                       const float inv2LT{1.0f/twoLT};
                       const float invLT{1.0f/(this->m_L*this->m_T)};
                       for(i = 0ull; i != ROUND_TO_EIGHT(this->m_nTLsamples,8ull); i += 8ull) 
                       {
                           const float t_i_0{static_cast<float>(i+0ull)};
#if (CPM_PULSE_SHAPERS_USE_CEPHES) == 1
                           const float cos_term_0{ceph_cosf(C6283185307179586476925286766559*t_i_0*invLT)};
#else
                           const float cos_term_0{std::cos(C6283185307179586476925286766559*t_i_0*invLT)};
#endif
                           const float lrc_val_0{inv2LT*(1.0f-cos_term_0)};
                           this->m_lrc_pulse.m_data[i+0ull] = lrc_val_0;
                           const float t_i_1{static_cast<float>(i+1ull)};
#if (CPM_PULSE_SHAPERS_USE_CEPHES) == 1
                           const float cos_term_1{ceph_cosf(C6283185307179586476925286766559*t_i_1*invLT)};
#else
                           const float cos_term_1{std::cos(C6283185307179586476925286766559*t_i_1*invLT)};
#endif
                           const float lrc_val_1{inv2LT*(1.0f-cos_term_1)};
                           this->m_lrc_pulse.m_data[i+1ull] = lrc_val_1;
                           const float t_i_2{static_cast<float>(i+2ull)};
#if (CPM_PULSE_SHAPERS_USE_CEPHES) == 1
                           const float cos_term_2{ceph_cosf(C6283185307179586476925286766559*t_i_2*invLT)};
#else
                           const float cos_term_2{std::cos(C6283185307179586476925286766559*t_i_2*invLT)};
#endif
                           const float lrc_val_2{inv2LT*(1.0f-cos_term_2)};
                           this->m_lrc_pulse.m_data[i+2ull] = lrc_val_2;
                           const float t_i_3{static_cast<float>(i+3ull)};
#if (CPM_PULSE_SHAPERS_USE_CEPHES) == 1
                           const float cos_term_3{ceph_cosf(C6283185307179586476925286766559*t_i_3*invLT)};
#else
                           const float cos_term_3{std::cos(C6283185307179586476925286766559*t_i_3*invLT)};
#endif
                           const float lrc_val_3{inv2LT*(1.0f-cos_term_3)};
                           this->m_lrc_pulse.m_data[i+3ull] = lrc_val_3;
                           const float t_i_4{static_cast<float>(i+4ull)};
#if (CPM_PULSE_SHAPERS_USE_CEPHES) == 1
                           const float cos_term_4{ceph_cosf(C6283185307179586476925286766559*t_i_4*invLT)};
#else
                           const float cos_term_4{std::cos(C6283185307179586476925286766559*t_i_4*invLT)};
#endif
                           const float lrc_val_4{inv2LT*(1.0f-cos_term_4)};
                           this->m_lrc_pulse.m_data[i+4ull] = lrc_val_4;
                           const float t_i_5{static_cast<float>(i+5ull)};
#if (CPM_PULSE_SHAPERS_USE_CEPHES) == 1 
                           const float cos_term_5{ceph_cosf(C6283185307179586476925286766559*t_i_5*invLT)};
#else
                           const float cos_term_5{std::cos(C6283185307179586476925286766559*t_i_5*invLT)};
#endif
                           const float lrc_val_5{inv2LT*(1.0f-cos_term_5)};
                           this->m_lrc_pulse.m_data[i+5ull] = lrc_val_5;
                           const float t_i_6{static_cast<float>(i+6ull)};
#if (CPM_PULSE_SHAPERS_USE_CEPHES) == 1
                           const float cos_term_6{ceph_cosf(C6283185307179586476925286766559*t_i_6*invLT)};
#else
                           const float cos_term_6{std::cos(C6283185307179586476925286766559*t_i_6*invLT)};
#endif
                           const float lrc_val_6{inv2LT*(1.0f-cos_term_6)};
                           this->m_lrc_pulse.m_data[i+6ull] = lrc_val_6;
                           const float t_i_7{static_cast<float>(i+7ull)};
#if (CPM_PULSE_SHAPERS_USE_CEPHES) == 1
                           const float cos_term_7{ceph_cosf(C6283185307179586476925286766559*t_i_7*invLT)};
#else
                           const float cos_term_7{std::cos(C6283185307179586476925286766559*t_i_7*invLT)};
#endif
                           const float lrc_val_7{inv2LT*(1.0f-cos_term_7)};
                           this->m_lrc_pulse.m_data[i+7ull] = lrc_val_7;

                       }

                       for(j = i; j != this->m_nTLsamples; ++j)   
                       {
                            const float t_j{static_cast<float>(j)};
#if (CPM_PULSE_SHAPERS_USE_CEPHES) == 1
                           const float cos_term_0{ceph_cosf(C6283185307179586476925286766559*t_j*invLT)};
#else
                           const float cos_term_0{std::cos(C6283185307179586476925286766559*t_j*invLT)};
#endif
                           const float lrc_val_0{inv2LT*(1.0f-cos_term_0)};
                           this->m_lrc_pulse.m_data[j] = lrc_val_0;
                       }
                       return (0);
                   }

                  __ATTR_ALWAYS_INLINE__
                  inline std::int32_t 
                  generate_lrc_pulse_sse_rolled()
                  {
                       using namespace gms::math;
                       constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
                       const float twoLT{2.0f*this->m_L*this->m_T};
                       const float inv2LT{1.0f/twoLT};
                       const float invLT{1.0f/(this->m_L*this->m_T)};
                       const __m128 vC6283185307179586476925286766559{_mm_set1_ps(C6283185307179586476925286766559)};
                       const __m128 vinv2LT{_mm_set1_ps(inv2LT)};
                       const __m128 vinvLT{_mm_set1_ps(invLT)};
                       const __m128 vone{_mm_set1_ps(1.0f)};
                       std::size_t i,j;
                       float jj;
                       constexpr std::size_t LUT_loop_idx_threshold{2257ull};
                       if(__builtin_expect(this->m_nTLsamples>LUT_loop_idx_threshold,0)) 
                       {
                           for(i=0ull,jj=0.0f;i != ROUND_TO_FOUR(this->m_nTLsamples,4ull); i+=4ull,jj+=4.0f)
                           {
                               const __m128 vt_i{_mm_setr_ps(jj*invLT,(jj+1.0f)*invLT,(jj+2.0f)*invLT,(jj+3.0f)*invLT)};
                               const __m128 vcos_term{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i))};
                               const __m128 vlrc_value{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term))};
                               _mm_store_ps(&this->m_lrc_pulse.m_data[i],vlrc_value);
                           }

                           for(j = i;j != this->m_nTLsamples; ++j)   
                           {
                               const float t_j{static_cast<float>(j)};
#if (CPM_PULSE_SHAPERS_USE_CEPHES) == 1
                               const float cos_term_0{ceph_cosf(C6283185307179586476925286766559*t_j*invLT)};
#else
                               const float cos_term_0{std::cos(C6283185307179586476925286766559*t_j*invLT)};
#endif
                               const float lrc_val_0{inv2LT*(1.0f-cos_term_0)};
                               this->m_lrc_pulse.m_data[j] = lrc_val_0;
                           }
                       }
                       else 
                       {
                           for(i = 0ull; i != ROUND_TO_FOUR(this->m_nTLsamples,4ull); i += 4ull) 
                           {
                               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i],_MM_HINT_T0);
                               const __m128 vt_i{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]),vinvLT)};
                               const __m128 vcos_term{_mm_cos_ps(_mm_mul_ps(vC6283185307179586476925286766559,vt_i))};
                               const __m128 vlrc_value{_mm_mul_ps(vinv2LT,_mm_sub_ps(vone,vcos_term))};
                               _mm_store_ps(&this->m_lrc_pulse.m_data[i],vlrc_value);
                           }
                           
                           for(j = i;j != this->m_nTLsamples; ++j)   
                           {
                               const float t_j{gms::math::LUT_loop_indices_2257_align16[j]};
#if (CPM_PULSE_SHAPERS_USE_CEPHES) == 1
                               const float cos_term_0{ceph_cosf(C6283185307179586476925286766559*t_j*invLT)};
#else
                               const float cos_term_0{std::cos(C6283185307179586476925286766559*t_j*invLT)};
#endif
                               const float lrc_val_0{inv2LT*(1.0f-cos_term_0)};
                               this->m_lrc_pulse.m_data[j] = lrc_val_0;
                           }
                       }
                       return (0);
                  }  

                  __ATTR_HOT__
                  __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
                  __ATTR_OPTIMIZE_03__
#endif                
                  std::int32_t 
                  generate_lrc_pulse_sse_u8x() noexcept(true);

                  __ATTR_ALWAYS_INLINE__
                  inline std::int32_t 
                  generate_lsrc_pulse_scalar_u8x() 
                  {
                       using namespace gms::math;
                       if(__builtin_expect(this->m_beta<0.0f,0) || 
                          __builtin_expect(this->m_beta>1.0f,0)) { return (-1);}                    
                       constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
                       const float invLT{1.0f/(this->m_L*this->m_T)};
                       std::size_t i,j;
                       for(i = 0ull; i != ROUND_TO_EIGHT(this->m_nTLsamples,8ull); i += 8ull) 
                       {
                           const float t_i_0{static_cast<float>(i+0ull)};
                           const float arg_term_0{C6283185307179586476925286766559*t_i_0*invLT};
                           const float cosinc_denom_0{4.0f*this->m_beta*t_i_0*invLT};
#if (CPM_PULSE_SHAPERS_USE_CEPHES) == 1
                           const float sin_term_0{ceph_sinf(arg_term_0)};
                           const float cos_term_0{ceph_cosf(this->m_beta*arg_term_0)};
#else 
                           const float sin_term_0{std::sin(arg_term_0)};
                           const float cos_term_0{std::cos(this->m_beta*arg_term_0)};
#endif
                           const float sinc_term_0{sin_term_0/arg_term_0};
                           const float cosinc_term_0{cos_term_0/(1.0f-(cosinc_denom_0*cosinc_denom_0))};
                           const float lsrc_sample_0{invLT*sinc_term_0*cosinc_term_0};
                           this->m_lsrc_pulse.m_data[i+0ull] = lsrc_sample_0;
                           const float t_i_1{static_cast<float>(i+1ull)};
                           const float arg_term_1{C6283185307179586476925286766559*t_i_1*invLT};
                           const float cosinc_denom_1{4.0f*this->m_beta*t_i_1*invLT};
#if (CPM_PULSE_SHAPERS_USE_CEPHES) == 1
                           const float sin_term_1{ceph_sinf(arg_term_1)};
                           const float cos_term_1{ceph_cosf(this->m_beta*arg_term_1)};
#else 
                           const float sin_term_1{std::sin(arg_term_1)};
                           const float cos_term_1{std::cos(this->m_beta*arg_term_1)};
#endif
                           const float sinc_term_1{sin_term_1/arg_term_1};
                           const float cosinc_term_1{cos_term_1/(1.0f-(cosinc_denom_1*cosinc_denom_1))};
                           const float lsrc_sample_1{invLT*sinc_term_1*cosinc_term_1};
                           this->m_lsrc_pulse.m_data[i+1ull] = lsrc_sample_1;
                           const float t_i_2{static_cast<float>(i+2ull)};
                           const float arg_term_2{C6283185307179586476925286766559*t_i_2*invLT};
                           const float cosinc_denom_2{4.0f*this->m_beta*t_i_2*invLT};
#if (CPM_PULSE_SHAPERS_USE_CEPHES) == 1
                           const float sin_term_2{ceph_sinf(arg_term_2)};
                           const float cos_term_2{ceph_cosf(this->m_beta*arg_term_2)};
#else 
                           const float sin_term_2{std::sin(arg_term_2)};
                           const float cos_term_2{std::cos(this->m_beta*arg_term_2)};
#endif
                           const float sinc_term_2{sin_term_2/arg_term_2};
                           const float cosinc_term_2{cos_term_2/(1.0f-(cosinc_denom_2*cosinc_denom_2))};
                           const float lsrc_sample_2{invLT*sinc_term_2*cosinc_term_2};
                           this->m_lsrc_pulse.m_data[i+2ull] = lsrc_sample_2;
                           const float t_i_3{static_cast<float>(i+3ull)};
                           const float arg_term_3{C6283185307179586476925286766559*t_i_3*invLT};
                           const float cosinc_denom_3{4.0f*this->m_beta*t_i_3*invLT};
#if (CPM_PULSE_SHAPERS_USE_CEPHES) == 1
                           const float sin_term_3{ceph_sinf(arg_term_3)};
                           const float cos_term_3{ceph_cosf(this->m_beta*arg_term_3)};
#else 
                           const float sin_term_3{std::sin(arg_term_3)};
                           const float cos_term_3{std::cos(this->m_beta*arg_term_3)};
#endif
                           const float sinc_term_3{sin_term_3/arg_term_3};
                           const float cosinc_term_3{cos_term_3/(1.0f-(cosinc_denom_3*cosinc_denom_3))};
                           const float lsrc_sample_3{invLT*sinc_term_3*cosinc_term_3};
                           this->m_lsrc_pulse.m_data[i+3ull] = lsrc_sample_3;
                           const float t_i_4{static_cast<float>(i+4ull)};
                           const float arg_term_4{C6283185307179586476925286766559*t_i_4*invLT};
                           const float cosinc_denom_4{4.0f*this->m_beta*t_i_4*invLT};
#if (CPM_PULSE_SHAPERS_USE_CEPHES) == 1
                           const float sin_term_4{ceph_sinf(arg_term_4)};
                           const float cos_term_4{ceph_cosf(this->m_beta*arg_term_4)};
#else 
                           const float sin_term_4{std::sin(arg_term_4)};
                           const float cos_term_4{std::cos(this->m_beta*arg_term_4)};
#endif
                           const float sinc_term_4{sin_term_4/arg_term_4};
                           const float cosinc_term_4{cos_term_4/(1.0f-(cosinc_denom_4*cosinc_denom_4))};
                           const float lsrc_sample_4{invLT*sinc_term_4*cosinc_term_4};
                           this->m_lsrc_pulse.m_data[i+4ull] = lsrc_sample_4;
                           const float t_i_5{static_cast<float>(i+5ull)};
                           const float arg_term_5{C6283185307179586476925286766559*t_i_5*invLT};
                           const float cosinc_denom_5{4.0f*this->m_beta*t_i_5*invLT};
#if (CPM_PULSE_SHAPERS_USE_CEPHES) == 1
                           const float sin_term_5{ceph_sinf(arg_term_5)};
                           const float cos_term_5{ceph_cosf(this->m_beta*arg_term_5)};
#else 
                           const float sin_term_5{std::sin(arg_term_5)};
                           const float cos_term_5{std::cos(this->m_beta*arg_term_5)};
#endif
                           const float sinc_term_5{sin_term_5/arg_term_5};
                           const float cosinc_term_5{cos_term_5/(1.0f-(cosinc_denom_5*cosinc_denom_5))};
                           const float lsrc_sample_5{invLT*sinc_term_5*cosinc_term_5};
                           this->m_lsrc_pulse.m_data[i+5ull] = lsrc_sample_5;
                           const float t_i_6{static_cast<float>(i+6ull)};
                           const float arg_term_6{C6283185307179586476925286766559*t_i_6*invLT};
                           const float cosinc_denom_6{4.0f*this->m_beta*t_i_6*invLT};
#if (CPM_PULSE_SHAPERS_USE_CEPHES) == 1
                           const float sin_term_6{ceph_sinf(arg_term_6)};
                           const float cos_term_6{ceph_cosf(this->m_beta*arg_term_6)};
#else 
                           const float sin_term_6{std::sin(arg_term_6)};
                           const float cos_term_6{std::cos(this->m_beta*arg_term_6)};
#endif
                           const float sinc_term_6{sin_term_6/arg_term_6};
                           const float cosinc_term_6{cos_term_6/(1.0f-(cosinc_denom_6*cosinc_denom_6))};
                           const float lsrc_sample_6{invLT*sinc_term_6*cosinc_term_6};
                           this->m_lsrc_pulse.m_data[i+6ull] = lsrc_sample_6;
                           const float t_i_7{static_cast<float>(i+7ull)};
                           const float arg_term_7{C6283185307179586476925286766559*t_i_7*invLT};
                           const float cosinc_denom_7{4.0f*this->m_beta*t_i_7*invLT};
#if (CPM_PULSE_SHAPERS_USE_CEPHES) == 1
                           const float sin_term_7{ceph_sinf(arg_term_7)};
                           const float cos_term_7{ceph_cosf(this->m_beta*arg_term_7)};
#else 
                           const float sin_term_7{std::sin(arg_term_7)};
                           const float cos_term_7{std::cos(this->m_beta*arg_term_7)};
#endif
                           const float sinc_term_7{sin_term_7/arg_term_7};
                           const float cosinc_term_7{cos_term_7/(1.0f-(cosinc_denom_7*cosinc_denom_7))};
                           const float lsrc_sample_7{invLT*sinc_term_7*cosinc_term_7};
                           this->m_lsrc_pulse.m_data[i+7ull] = lsrc_sample_7;
                       }

                       for(j = i;j != this->m_nTLsamples; ++j)  
                       {
                           const float t_j{static_cast<float>(j)};
                           const float arg_term_0{C6283185307179586476925286766559*t_j*invLT};
                           const float cosinc_denom_0{4.0f*this->m_beta*t_j*invLT};
#if (CPM_PULSE_SHAPERS_USE_CEPHES) == 1
                           const float sin_term_0{ceph_sinf(arg_term_0)};
                           const float cos_term_0{ceph_cosf(this->m_beta*arg_term_0)};
#else 
                           const float sin_term_0{std::sin(arg_term_0)};
                           const float cos_term_0{std::cos(this->m_beta*arg_term_0)};
#endif
                           const float sinc_term_0{sin_term_0/arg_term_0};
                           const float cosinc_term_0{cos_term_0/(1.0f-(cosinc_denom_0*cosinc_denom_0))};
                           const float lsrc_sample_0{invLT*sinc_term_0*cosinc_term_0};
                           this->m_lsrc_pulse.m_data[j] = lsrc_sample_0;
                       }
                       return (0);
                  }

                  __ATTR_ALWAYS_INLINE__
                  inline std::int32_t 
                  generate_lsrc_pulse_sse_rolled()
                  {
                      using namespace gms::math;
                      constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
                      constexpr std::size_t LUT_loop_idx_threshold{2257ull};
                      std::size_t i,j;
                      float       jj;
                      const float invLT{1.0f/(this->m_L*this->m_T)};
                      const __m128 vinvLT{_mm_set1_ps(invLT)};
                      const __m128 vC6283185307179586476925286766559{_mm_set1_ps(C6283185307179586476925286766559)};
                      const __m128 vone{_mm_set1_ps(1.0f)};
                      const __m128 vfour{_mm_set1_ps(4.0f)};
                      const __m128 vbeta{_mm_set1_ps(this->m_beta)};
                      if(__builtin_expect(this->m_nTLsamples>LUT_loop_indices_2257_align16,0)) 
                      {
                           for(i = 0ull,jj = 0.0f; i != ROUND_TO_FOUR(this->m_nTLsamples,4ull); i += 4ull,jj += 4.0f) 
                           {
                               const __m128 vt_i{_mm_setr_ps(jj*vinvLT,(jj+1.0f)*vinvLT,(jj+2.0f)*vinvLT,(jj+3.0f)*vinvLT)};
                               const __m128 varg_term{_mm_mul_ps(C6283185307179586476925286766559,vt_i)};
                               const __m128 vcosinc_denom{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i))};
                               const __m128 vsin_term{_mm_sin_ps(varg_term)};
                               const __m128 vsinc_term{_mm_div_ps(vsin_term,varg_term)};
                               const __m128 vcos_term{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term))};
                               const __m128 vcosinc_term{_mm_div_ps(vcos_term,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom,vcosinc_denom)))};
                               const __m128 vlsrc_sample{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term,vcosinc_term))};
                               _mm_store_ps(&this->m_lsrc_pulse.m_data[i], vlsrc_sample);
                           }

                           for(j = i; j != this->m_nTLsamples; ++j)  
                           {
                                const float t_j{static_cast<float>(j)};
                                const float arg_term_0{C6283185307179586476925286766559*t_j*invLT};
                                const float cosinc_denom_0{4.0f*this->m_beta*t_j*invLT};
#if (CPM_PULSE_SHAPERS_USE_CEPHES) == 1
                                const float sin_term_0{ceph_sinf(arg_term_0)};
                                const float cos_term_0{ceph_cosf(this->m_beta*arg_term_0)};
#else 
                                const float sin_term_0{std::sin(arg_term_0)};
                                const float cos_term_0{std::cos(this->m_beta*arg_term_0)};
#endif
                                const float sinc_term_0{sin_term_0/arg_term_0};
                                const float cosinc_term_0{cos_term_0/(1.0f-(cosinc_denom_0*cosinc_denom_0))};
                                const float lsrc_sample_0{invLT*sinc_term_0*cosinc_term_0};
                                this->m_lsrc_pulse.m_data[j] = lsrc_sample_0;
                           }
                      }
                      else 
                      {
                           for(i = 0ull; i != ROUND_TO_FOUR(this->m_nTLsamples,4ull); i += 4ull) 
                           { 
                               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i],_MM_HINT_T0);
                               const __m128 vt_i{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]),vinvLT)};
                               const __m128 varg_term{_mm_mul_ps(C6283185307179586476925286766559,vt_i)};
                               const __m128 vcosinc_denom{_mm_mul_ps(vfour,_mm_mul_ps(vbeta,vt_i))};
                               const __m128 vsin_term{_mm_sin_ps(varg_term)};
                               const __m128 vsinc_term{_mm_div_ps(vsin_term,varg_term)};
                               const __m128 vcos_term{_mm_cos_ps(_mm_mul_ps(vbeta,varg_term))};
                               const __m128 vcosinc_term{_mm_div_ps(vcos_term,_mm_sub_ps(vone,_mm_mul_ps(vcosinc_denom,vcosinc_denom)))};
                               const __m128 vlsrc_sample{_mm_mul_ps(vinvLT,_mm_mul_ps(vsinc_term,vcosinc_term))};
                               _mm_store_ps(&this->m_lsrc_pulse.m_data[i], vlsrc_sample);
                           }

                            for(j = i; j != this->m_nTLsamples; ++j)  
                           {
                                const float t_j{gms::math::LUT_loop_indices_2257_align16[j]*invLT};
                                const float arg_term_0{C6283185307179586476925286766559*t_j};
                                const float cosinc_denom_0{4.0f*this->m_beta*t_j};
#if (CPM_PULSE_SHAPERS_USE_CEPHES) == 1
                                const float sin_term_0{ceph_sinf(arg_term_0)};
                                const float cos_term_0{ceph_cosf(this->m_beta*arg_term_0)};
#else 
                                const float sin_term_0{std::sin(arg_term_0)};
                                const float cos_term_0{std::cos(this->m_beta*arg_term_0)};
#endif
                                const float sinc_term_0{sin_term_0/arg_term_0};
                                const float cosinc_term_0{cos_term_0/(1.0f-(cosinc_denom_0*cosinc_denom_0))};
                                const float lsrc_sample_0{invLT*sinc_term_0*cosinc_term_0};
                                this->m_lsrc_pulse.m_data[j] = lsrc_sample_0;
                           }
                      }
                      return (0);
                  }

                  __ATTR_HOT__
                  __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
                  __ATTR_OPTIMIZE_03__
#endif                
                  std::int32_t 
                  generate_lsrc_pulse_sse_u8x() noexcept(true);

                  __ATTR_ALWAYS_INLINE__ 
                  inline float 
                  tfm_g0_value(const float t,
                               const float invT,
                               const float invTT)
                  {
                       using namespace gms::math;
                       constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
                       constexpr float C314159265358979323846264338328{3.14159265358979323846264338328f};
                       constexpr float C157079632679489661923132169164{1.57079632679489661923132169164f};
                       constexpr float C0125{1.0f/8.0f};
                       const float invPIt{1.0f/(C314159265358979323846264338328*t)};
                       const float PItT{C314159265358979323846264338328*t*invT};
                       const float tpow3{t*t*t};
                       const float twoPItT{C6283185307179586476925286766559*t*invT};
                       const float denom{24.0f*C314159265358979323846264338328*tpow3*invTT};
#if (CPM_PULSE_SHAPERS_USE_CEPHES) == 1
                       const float sin_term{ceph_sinf(PItT)};
#else 
                       const float sin_term{std::sin(PItT)};
#endif 
                       const float cot_term{std::tan(C157079632679489661923132169164-PItT)};
                       const float sqrPItT{PItT*PItT};
                       const float numerator{2.0f-twoPItT*cot_term-sqrPItT};
                       const float ratio{invPIt-(numerator/denom)};
                       const float g0_value{sin_term*ratio};
                       return (g0_value);
                  }

                  __ATTR_ALWAYS_INLINE__ 
                  inline std::int32_t 
                  generate_tfm_pulse_scalar_u8x() noexcept(true) 
                  {
                       using namespace gms::math;
                       constexpr float one_ov_eight{1.0f/8.0f};
                       std::size_t i,j;  
                       const float invT{1.0f/this->m_T};
                       const float invTT{1.0f/(this->m_T*this->m_T)};
                       for(i = 0ull; i != ROUND_TO_EIGHT(this->m_nTsamples,8ull); i += 8ull) 
                       {
                            const float t_i_0{static_cast<float>(i+0ull)};
                            const float g0_left_0{tfm_g0_value((t_i_0-this->m_T),invT,invTT)};
                            const float g0_mid_0{2.0f*tfm_g0_value(t_i_0,invT,invTT)};
                            const float g0_right_0{tfm_g0_value((t_i_0+this->m_T,invT,invTT))};
                            const float tfm_sample_0{one_ov_eight*(g0_left_0+g0_mid_0+g0_right_0)};
                            this->m_tfm_pulse.m_data[i+0ull] = tfm_sample_0;
                            const float t_i_1{static_cast<float>(i+1ull)};
                            const float g0_left_1{tfm_g0_value((t_i_1-this->m_T),invT,invTT)};
                            const float g0_mid_1{2.0f*tfm_g0_value(t_i_1,invT,invTT)};
                            const float g0_right_1{tfm_g0_value((t_i_1+this->m_T,invT,invTT))};
                            const float tfm_sample_1{one_ov_eight*(g0_left_1+g0_mid_1+g0_right_1)};
                            this->m_tfm_pulse.m_data[i+1ull] = tfm_sample_1;
                            const float t_i_2{static_cast<float>(i+2ull)};
                            const float g0_left_2{tfm_g0_value((t_i_2-this->m_T),invT,invTT)};
                            const float g0_mid_2{2.0f*tfm_g0_value(t_i_2,invT,invTT)};
                            const float g0_right_2{tfm_g0_value((t_i_2+this->m_T,invT,invTT))};
                            const float tfm_sample_2{one_ov_eight*(g0_left_2+g0_mid_2+g0_right_2)};
                            this->m_tfm_pulse.m_data[i+2ull] = tfm_sample_2;
                            const float t_i_3{static_cast<float>(i+3ull)};
                            const float g0_left_3{tfm_g0_value((t_i_3-this->m_T),invT,invTT)};
                            const float g0_mid_3{2.0f*tfm_g0_value(t_i_3,invT,invTT)};
                            const float g0_right_3{tfm_g0_value((t_i_3+this->m_T,invT,invTT))};
                            const float tfm_sample_3{one_ov_eight*(g0_left_3+g0_mid_3+g0_right_3)};
                            this->m_tfm_pulse.m_data[i+3ull] = tfm_sample_3;
                            const float t_i_4{static_cast<float>(i+4ull)};
                            const float g0_left_4{tfm_g0_value((t_i_4-this->m_T),invT,invTT)};
                            const float g0_mid_4{2.0f*tfm_g0_value(t_i_4,invT,invTT)};
                            const float g0_right_4{tfm_g0_value((t_i_4+this->m_T,invT,invTT))};
                            const float tfm_sample_4{one_ov_eight*(g0_left_4+g0_mid_4+g0_right_4)};
                            this->m_tfm_pulse.m_data[i+4ull] = tfm_sample_4;
                            const float t_i_5{static_cast<float>(5+1ull)};
                            const float g0_left_5{tfm_g0_value((t_i_5-this->m_T),invT,invTT)};
                            const float g0_mid_5{2.0f*tfm_g0_value(t_i_5,invT,invTT)};
                            const float g0_right_5{tfm_g0_value((t_i_5+this->m_T,invT,invTT))};
                            const float tfm_sample_5{one_ov_eight*(g0_left_5+g0_mid_5+g0_right_5)};
                            this->m_tfm_pulse.m_data[i+5ull] = tfm_sample_5;
                            const float t_i_6{static_cast<float>(6+1ull)};
                            const float g0_left_6{tfm_g0_value((t_i_6-this->m_T),invT,invTT)};
                            const float g0_mid_6{2.0f*tfm_g0_value(t_i_6,invT,invTT)};
                            const float g0_right_6{tfm_g0_value((t_i_6+this->m_T,invT,invTT))};
                            const float tfm_sample_6{one_ov_eight*(g0_left_6+g0_mid_6+g0_right_6)};
                            this->m_tfm_pulse.m_data[i+6ull] = tfm_sample_6;
                            const float t_i_7{static_cast<float>(i+7ull)};
                            const float g0_left_7{tfm_g0_value((t_i_7-this->m_T),invT,invTT)};
                            const float g0_mid_7{2.0f*tfm_g0_value(t_i_7,invT,invTT)};
                            const float g0_right_7{tfm_g0_value((t_i_7+this->m_T,invT,invTT))};
                            const float tfm_sample_7{one_ov_eight*(g0_left_7+g0_mid_7+g0_right_7)};
                            this->m_tfm_pulse.m_data[i+7ull] = tfm_sample_7;
                       }

                       for(j = i; j != this->m_nTsamples; ++j)  
                       {
                            const float t_j{static_cast<float>(j)};
                            const float g0_left_0{tfm_g0_value((t_j-this->m_T),invT,invTT)};
                            const float g0_mid_0{2.0f*tfm_g0_value(t_j,invT,invTT)};
                            const float g0_right_0{tfm_g0_value((t_j+this->m_T,invT,invTT))};
                            const float tfm_sample_0{one_ov_eight*(g0_left_0+g0_mid_0+g0_right_0)};
                            this->m_tfm_pulse.m_data[j] = tfm_sample_0;
                       }
                       return (0);
                  }

                  __ATTR_ALWAYS_INLINE__
                  inline __m128 
                  _mm_tfm_g0_vec_ps(const __m128 vt,
                                    const __m128 vinvT,
                                    const __m128 vinvTT) 
                  {
                        
                       const __m128 vC6283185307179586476925286766559{_mm_set1_ps(6.283185307179586476925286766559f)};
                       const __m128 vC314159265358979323846264338328{_mm_set1_ps(3.14159265358979323846264338328f)};
                       const __m128 vC157079632679489661923132169164{_mm_set1_ps(1.57079632679489661923132169164f)};
                       const __m128 v24{_mm_set1_ps(24.0f)};
                       const __m128 v2{_mm_set1_ps(2.0f)};
#if (CPM_PULSE_SHAPERS_USE_VEC_RCP14_PRECISION) == 1 
                       const __m128 vinvPIt{_mm_rcp14_ps(_mm_mul_ps(vC314159265358979323846264338328,vt))};
#else 
                       const __m128 vONE{_mm_set1_ps(1.0f)};
                       const __m128 vinvPIt{_mm_div_ps(vONE,_mm_mul_ps(vC314159265358979323846264338328,vt))};
#endif 
                       const __m128 vPItT{_mm_mul_ps(vC314159265358979323846264338328,_mm_mul_ps(vt,vinvT))};
                       const __m128 vtpow3{_mm_mul_ps(vt,_mm_mul_ps(vt,vt))};
                       const __m128 v2PItT{_mm_mul_ps(vC6283185307179586476925286766559,_mm_mul_ps(vt,vinvT))};
                       const __m128 vdenom{_mm_mul_ps(_mm_mul_ps(v24,vC314159265358979323846264338328),_mm_mul_ps(vtpow3,vinvTT))};
                       const __m128 vsin_term{_mm_sin_ps(vPItT)};
                       const __m128 vcot_term{_mm_tan_ps(_mm_sub_ps(vC157079632679489661923132169164,vPItT))};
                       const __m128 vsqrPItT{_mm_mul_ps(vPItT,vPItT)};
                       const __m128 vnumerator{_mm_sub_ps(v2,_mm_fmsub_ps(v2PItT,vcot_term,vsqrPItT))};
                       const __m128 vratio{_mm_sub_ps(vinvPIt,_mm_div_ps(vnumerator,vdenom))};
                       const __m128 vg0_vec{_mm_mul_ps(vsin_term,vratio)};
                       return (vg0_vec);
                  }

                  __ATTR_ALWAYS_INLINE__ 
                  inline std::int32_t 
                  generate_tfm_pulse_sse_rolled() noexcept(true) 
                  {
                       using namespace gms::math;
                       constexpr std::size_t LUT_loop_idx_threshold{2257ull};
                       constexpr float one_ov_eight{1.0f/8.0f};
                       std::size_t i,j;  
                       float jj;
                       const float invT{1.0f/this->m_T};
                       const float invTT{1.0f/(this->m_T*this->m_T)};
                       const __m128 vinvT{_mm_set1_ps(invT)};
                       const __m128 vinvTT{_mm_set1_ps(invTT)};
                       const __m128 v2{_mm_set1_ps(2.0f)};
                       if(__builtin_expect(this->m_nTsamples>LUT_loop_indices_2257_align16,0))
                       {
                             for(i = 0ull,jj = 0.0f;i != ROUND_TO_FOUR(this->m_nTsamples,4ull);i += 4ull,j += 4.0f) 
                             {
                                   const __m128 vt_i_left{_mm_setr_ps(jj-this->m_T,jj-this->m_T+1.0f,jj-this->m_T+2.0f,jj-this->m_T+3.0f)};
                                   const __m128 g0_vec_left{_mm_tfm_g0_vec_ps(vt_i_left,vinvT,vinvTT)};
                                   const __m128 vt_i_mid{_mm_setr_ps(jj,jj+1.0f,jj+2.0f,jj+3.0f)};
                                   const __m128 g0_vec_mid{_mm_tfm_g0_vec_ps(vt_i_mid,vinvT,vinvTT)};
                                   const __m128 vt_i_right{_mm_setr_ps(jj+this->m_T,jj+this->m_T+1.0f,jj+this->m_T+2.0f,jj+this->m_T+3.0f)};
                                   const __m128 g0_vec_right{_mm_tfm_g0_vec_ps(vt_i_right,vinvT,vinvTT)};
                                   const __m128 tfm_g_vec_sample{_mm_mul_ps(one_ov_eight,
                                                _mm_add_ps(g0_vec_left,_mm_add_ps(_mm_mul_ps(v2,g0_vec_mid),g0_vec_right)))};
                                   _mm_store_ps(&this->m_tfm_pulse.m_data[i],tfm_g_vec_sample);
                             }

                             for(j = i; j != this->m_nTsamples; ++j)  
                             {
                                   const float t_j{static_cast<float>(j)};
                                   const float g0_left_0{tfm_g0_value((t_j-this->m_T),invT,invTT)};
                                   const float g0_mid_0{2.0f*tfm_g0_value(t_j,invT,invTT)};
                                   const float g0_right_0{tfm_g0_value((t_j+this->m_T,invT,invTT))};
                                   const float tfm_sample_0{one_ov_eight*(g0_left_0+g0_mid_0+g0_right_0)};
                                   this->m_tfm_pulse.m_data[j] = tfm_sample_0;
                             }
                       }
                       else 
                       {
                             const __m128 vT{_mm_set1_ps(this->m_T)};
                             for(i = 0ull; i != ROUND_TO_FOUR(this->m_nTsamples,4ull); i += 4ull) 
                             {
                                     _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i],_MM_HINT_T0);
                                     const __m128 vidx{_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i])};
                                     const __m128 vt_i_left{_mm_sub_ps(vidx,vT)};
                                     const __m128 g0_vec_left{_mm_tfm_g0_vec_ps(vt_i_left,vinvT,vinvTT)};
                                     const __m128 vt_i_mid{vidx};
                                     const __m128 g0_vec_mid{_mm_tfm_g0_vec_ps(vt_i_mid,vinvT,vinvTT)};
                                     const __m128 vt_i_right{_mm_add_ps(vidx,vT)};
                                     const __m128 g0_vec_right{_mm_tfm_g0_vec_ps(vt_i_right,vinvT,vinvTT)};
                                     const __m128 tfm_g_vec_sample{_mm_mul_ps(one_ov_eight,
                                                _mm_add_ps(g0_vec_left,_mm_add_ps(_mm_mul_ps(v2,g0_vec_mid),g0_vec_right)))};
                                     _mm_store_ps(&this->m_tfm_pulse.m_data[i],tfm_g_vec_sample);
                             }

                             for(j = i; j != this->m_nTsamples; ++j)  
                             {
                                   const float t_j{gms::math::LUT_loop_indices_2257_align16[j]};
                                   const float g0_left_0{tfm_g0_value((t_j-this->m_T),invT,invTT)};
                                   const float g0_mid_0{2.0f*tfm_g0_value(t_j,invT,invTT)};
                                   const float g0_right_0{tfm_g0_value((t_j+this->m_T,invT,invTT))};
                                   const float tfm_sample_0{one_ov_eight*(g0_left_0+g0_mid_0+g0_right_0)};
                                   this->m_tfm_pulse.m_data[j] = tfm_sample_0;
                             }
                       }
                       return (0);
                  }

                   __ATTR_HOT__
                  __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
                  __ATTR_OPTIMIZE_03__
#endif                
                  std::int32_t 
                  generate_tfm_pulse_sse_u8x() noexcept(true);

          };    

          auto 
          operator<<(std::ostream &, const cpm_pulse_shapers_t &)->std::ostream &;     

} //radiolocation

}// gms

















#endif /*__GMS_CPM_PULSE_SHAPERS_H__*/