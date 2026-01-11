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

#ifndef __GMS_SINUSOIDAL_FSK_H__
#define __GMS_SINUSOIDAL_FSK_H__ 101220250143

namespace file_info 
{

     static const unsigned int GMS_SINUSOIDAL_FSK_MAJOR = 1;
     static const unsigned int GMS_SINUSOIDAL_FSK_MINOR = 1;
     static const unsigned int GMS_SINUSOIDAL_FSK_MICRO = 0;
     static const unsigned int GMS_SINUSOIDAL_FSK_FULLVER =
       1000U*GMS_SINUSOIDAL_FSK_MAJOR+100U*GMS_SINUSOIDAL_FSK_MINOR+
       10U*GMS_SINUSOIDAL_FSK_MICRO;
     static const char GMS_SINUSOIDAL_FSK_CREATION_DATE[] = "10-12-2025 01:38PM +00200 (WED 10 DEC 2025 GMT+2)";
     static const char GMS_SINUSOIDAL_FSK_BUILD_DATE[]    = __DATE__; 
     static const char GMS_SINUSOIDAL_FSK_BUILD_TIME[]    = __TIME__;
     static const char GMS_SINUSOIDAL_FSK_SYNOPSIS[]      = "Snusoidal Frequency Shift Keyeing signal.";

}

#include <cstdint>
#include <string>
#include <iostream>
#include <thread>
#include <immintrin.h>
#include <random>
#include "GMS_config.h"
#include "GMS_dyn_array.h"
#include "GMS_sse_memset.h"
#include "GMS_iq_rectw_bitstream_vsequence.h"
#include "GMS_loop_indices_LUT.h"

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
// To be added.
#if !defined (SINUSOIDAL_FSK_USE_NT_STORES)
#define SINUSOIDAL_FSK_USE_NT_STORES 0
#endif

#ifdef _OPENMP
// Default init a storage arrays for the first-touch (OpenMP) processing
#if !defined(SINUSOIDAL_FSK_INIT_STORAGE)
#define SINUSOIDAL_FSK_INIT_STORAGE 1
#endif 
#endif

#if (SINUSOIDAL_FSK_INIT_STORAGE) == 1
#define INIT_BY_STD_FILL 0
#endif 

// For inlining of trigo functions (asin,acos,sin,cos)
#if !defined(SINUSOIDAL_FSK_USE_CEPHES)
#define SINUSOIDAL_FSK_USE_CEPHES 1
#endif 

// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (SINUSOIDAL_FSK_USE_PMC_INSTRUMENTATION)
#define SINUSOIDAL_FSK_USE_PMC_INSTRUMENTATION 0
#endif 

#if (SINUSOIDAL_FSK_USE_PMC_INSTRUMENTATION) == 1
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

#if (SINUSOIDAL_FSK_USE_CEPHES) == 0
#include <cmath>
#else 
#include "GMS_cephes_sin_cos.h"
#endif 

namespace gms
{

namespace radiolocation 
{

          

              /*
                Useful for noise added modulation
           */
             struct alignas(64) sinusoidal_fsk_pdf_params_t 
             {
                   float uni_real_a_r;
                   float uni_real_a_i;
                   float uni_real_b_r;
                   float uni_real_b_i;
                   float norm_mu_r;
                   float norm_mu_i;
                   float norm_sigma_r;
                   float norm_sigma_i;
                   float cauchy_a_r;
                   float cauchy_a_i;
                   float cauchy_b_r;
                   float cauchy_b_i;
                   float log_norm_m_r;
                   float log_norm_m_i;
                   float log_norm_s_r;
                   float log_norm_s_i;
                   float expo_gamma_r;
                   float expo_gamma_i;
                   float weibull_a_r;
                   float weibull_a_i;
                   float weibull_b_r;
                   float weibull_b_i;
                   float gamma_alph_r;
                   float gamma_alph_i;
                   float gamma_bet_r;
                   float gamma_bet_i;
                   int32_t poisson_mu;

            };

            

            struct alignas(64) sinusoidal_fsk_t final 
            {
                   std::size_t   m_sfsk_nsamples; // sfsk signal number of samples //function of k*T
                   std::size_t   m_I_ch_nsamples; // I-channel i.e. rect wave number of samples
                   std::size_t   m_Q_ch_nsamples; // Q-channel i.e. rect wave number of samples
                   std::size_t   m_psfunc_nsamples; // pulse shaping function number of samples
                   float         m_k; // number of k periods (k*T)
                   float         m_T; // a period T
                   float         m_Ac; // carrier amplitude
                   float         m_fc; // carrier frequency (central)
                   float         m_ph0; // carrier initial phase shift
                   float         m_I_fc; // bitstream sine wave I-channel center frequency
                   float         m_Q_fc; // bitstream sine wave Q-channel center frequency
                   float         m_I_ph0;// bitstream I-channel initial phase shift
                   float         m_Q_ph0;// bitstream Q-channel initial phase shift
                   darray_r4_t   m_I_ch_bitstream; //user data, i.e. [-1,1] currently regularly alternating pattern
                   darray_r4_t   m_Q_ch_bitstream; //user data, i.e. [-1,1] currently regularly alternating pattern
                   darray_r4_t   m_psfunc; // pulse shaping function (phi(t) = cos(pi*t/2*T*0.25*sin(2*pi*t/T)))
                   darray_r4_t   m_sfsk_signal;

                     /*
                Useful for noise added modulation
           */
                 enum class sinusoidal_fsk_rand_gens : int32_t 
                 {
                        rg_minstd_rand0,
                        rg_minstd_rand,
                        rg_mt19937,
                        rg_mt19937_64,
                        rg_ranlux24_base,
                        rg_ranlux48_base,
                        rg_ranlux24,
                        rg_ranlux48,
                        rg_knuth_b
                  };

                 enum class sinusoidal_fsk_rand_distr : int32_t 
                 {
                        uniform,
                        normal,
                        cauchy,
                        log_norm,
                        expo_gamma,
                        weibull,
                        gamma
                 };

                 enum class I_channel_bitstream_optimized_path : int32_t 
                 {
                        default_scalar_path,
                        vector_sse_path,
                        vector_sse_u4x_path,
                        vector_avx_path,
                        vector_avx_u4x_path,
                        vector_avx512_path,
                        vector_avx512_u4x_path
                 };

                 enum class Q_channel_bitstream_optimized_path : int32_t 
                 {
                        default_scalar_path,
                        vector_sse_path,
                        vector_sse_u4x_path,
                        vector_avx_path,
                        vector_avx_u4x_path,
                        vector_avx512_path,
                        vector_avx512_u4x_path
                 };

                 enum class pulse_shaping_function_optimized_path : int32_t 
                 {
                        default_scalar_path,
                        vector_sse_path,
                        vector_avx_path,
                        vector_avx512_path
                 };    

                   sinusoidal_fsk_t() noexcept(true);

                   sinusoidal_fsk_t(const std::size_t,
                                    const std::size_t,
                                    const std::size_t,
                                    const std::size_t,
                                    const float,
                                    const float,
                                    const float,
                                    const float,
                                    const float,
                                    const float,
                                    const float,
                                    const float,
                                    const float) noexcept(false);

                   sinusoidal_fsk_t(const sinusoidal_fsk_t &) noexcept(false);

                   sinusoidal_fsk_t(sinusoidal_fsk_t &&) noexcept(true);

                   ~sinusoidal_fsk_t() noexcept(true);

                   sinusoidal_fsk_t & operator=(const sinusoidal_fsk_t &) noexcept(false);

                   sinusoidal_fsk_t & operator=(sinusoidal_fsk_t &&) noexcept(true);

                   __ATTR_ALWAYS_INLINE__ 
                   inline void clear_I_ch_bitstream() noexcept(true)
                   {
#if (INIT_BY_STD_FILL) == 0
                      using namespace gms::common;
	                    sse_memset_unroll8x_ps(this->m_I_ch_bitstream.m_data,0.0f,this->m_I_ch_nsamples);
#else 
                      std::fill(this->m_I_ch_bitstream.m_data,this->m_I_ch_bitstream.m_data+this->m_I_ch_nsamples,0.0f);
#endif   
                   }

                   __ATTR_ALWAYS_INLINE__ 
                   inline void clear_Q_ch_bitstream() noexcept(true)
                   {
#if (INIT_BY_STD_FILL) == 0
                      using namespace gms::common;
	                    sse_memset_unroll8x_ps(this->m_Q_ch_bitstream.m_data,0.0f,this->m_Q_ch_nsamples);
#else 
                      std::fill(this->m_Q_ch_bitstream.m_data,this->m_Q_ch_bitstream.m_data+this->m_Q_ch_nsamples,0.0f);
#endif   
                   }

                   __ATTR_ALWAYS_INLINE__ 
                   inline void clear_psfunc() noexcept(true)
                   {
#if (INIT_BY_STD_FILL) == 0
                      using namespace gms::common;
	                    sse_memset_unroll8x_ps(this->m_psfunc.m_data,0.0f,this->m_psfunc_nsamples);
#else 
                      std::fill(this->m_psfunc.m_data,this->m_psfunc.m_data+this->m_psfunc_nsamples,0.0f);
#endif                    
                   }

                   __ATTR_ALWAYS_INLINE__
                   inline void clear_sfsk_signal() noexcept(true) 
                   {
#if (INIT_BY_STD_FILL) == 0
                      using namespace gms::common;
	                    sse_memset_unroll8x_ps(this->m_sfsk_signal.m_data,0.0f,this->m_sfsk_nsamples);
#else 
                      std::fill(this->m_sfsk_signal.m_data,this->m_sfsk_signal.m_data+this->m_sfsk_nsamples,0.0f);
#endif                       
                   }

                   static void create_signal_plot( const std::uint32_t,
                                                   const float * __restrict,
                                                   const float * __restrict,
                                                   const std::string &,
                                                   const std::string &,
                                                   const bool );

                  __ATTR_ALWAYS_INLINE__
                  inline          
                  float 
                  channel_I_sample_noise_n1_0_1(const float lo, // either -1 or 0
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
                  channel_Q_sample_noise_n1_0_1(const float lo, // either -1 or 0
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

                  __ATTR_ALWAYS_INLINE__
                  inline 
                  std::int32_t 
                  generate_I_channel_bitstream(const std::int32_t which_samples_len)
                  {
                      using namespace gms::math;
                      constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
                      switch (which_samples_len) 
                      {
                           case 1 : // over T duration
                           {
                                 const float invT{1.0f/this->m_T};
                                 for(std::size_t i{0ull}; i != this->m_I_ch_nsamples; ++i) 
                                 {
                                      const float t_i{static_cast<float>(i*invT)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
                                      const float cos_val{ceph_cosf(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_i))};
                                      this->m_I_ch_bitstream.m_data[i] = (cos_val>=0.0f)?+1.0f:-1.0f;
#else 
                                      const float cos_val{std::cos(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_i))};
                                      this->m_I_ch_bitstream.m_data[i] = (cos_val>=0.0f)?+1.0f:-1.0f;
#endif 
                                 }
                           }
                           break;
                           case 2 : // over k*T durations
                           {
                                 const float invT{1.0f/this->m_I_ch_nsamples};
                                 for(std::size_t i{0ull}; i != this->m_I_ch_nsamples; ++i) 
                                 {
                                      const float t_i{static_cast<float>(i*invT)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
                                      const float cos_val{ceph_cosf(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_i))};
                                      this->m_I_ch_bitstream.m_data[i] = (cos_val>=0.0f)?+1.0f:-1.0f;
#else 
                                      const float cos_val{std::cos(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_i))};
                                      this->m_I_ch_bitstream.m_data[i] = (cos_val>=0.0f)?+1.0f:-1.0f;
#endif 
                                 }
                           }
                           break;
                           default : 
                               return (-1);
                      }

                      return (0);
                     
                  }


                  __ATTR_ALWAYS_INLINE__
                  inline 
                  std::int32_t 
                  generate_Q_channel_bitstream(const std::int32_t which_samples_len)
                  {
                      using namespace gms::math;
                      constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
                      switch (which_samples_len) 
                      {
                           case 1 : // over T duration
                           {
                                 const float invT{1.0f/this->m_T};
                                 for(std::size_t i{0ull}; i != this->m_Q_ch_nsamples; ++i) 
                                 {
                                      const float t_i{static_cast<float>(i*invT)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
                                      const float cos_val{ceph_cosf(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_i))};
                                      this->m_Q_ch_bitstream.m_data[i] = (cos_val>=0.0f)?+1.0f:-1.0f;
#else 
                                      const float cos_val{std::cos(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_i))};
                                      this->m_Q_ch_bitstream.m_data[i] = (cos_val>=0.0f)?+1.0f:-1.0f;
#endif 
                                 }
                           }
                           break;
                           case 2 : // over k*T durations
                           {
                                 const float invT{1.0f/this->m_Q_ch_nsamples};
                                 for(std::size_t i{0ull}; i != this->m_Q_ch_nsamples; ++i) 
                                 {
                                      const float t_i{static_cast<float>(i*invT)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
                                      const float cos_val{ceph_cosf(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_i))};
                                      this->m_Q_ch_bitstream.m_data[i] = (cos_val>=0.0f)?+1.0f:-1.0f;
#else 
                                      const float cos_val{std::cos(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_i))};
                                      this->m_Q_ch_bitstream.m_data[i] = (cos_val>=0.0f)?+1.0f:-1.0f;
#endif 
                                 }
                           }
                           break;
                           default : 
                               return (-1);
                      }

                      return (0);
                     
                  }

                  __ATTR_ALWAYS_INLINE__ 
                  inline 
                  std::int32_t 
                  generate_I_channel_bitstream_sse()
                  {
                       using namespace gms::math;
                       std::size_t i,j;
                       const float inv_sr{1.0f/this->m_I_ch_nsamples};
                       const __m128 vinv_sr{_mm_set1_ps(inv_sr)};
                       const __m128 vfc{_mm_set1_ps(this->m_I_fc)};
                       const __m128 vph0{_mm_set1_ps(this->m_I_ph0)};
                       const __m128 v2pi{_mm_set1_ps(6.283185307179586476925286766559f)};
                       const __m128 vzero{_mm_setzero_ps()};
                       const __m128 vpone{_mm_set1_ps(+1.0f)};
                       const __m128 vnone{_mm_set1_ps(-1.0f)};
                       float jj;
                       constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
                       constexpr std::size_t LUT_loop_idx_threshold{2257ull};
                       if(this->m_I_ch_nsamples>LUT_loop_idx_threshold)
                       {
                            for(i = 0ull,jj = 0.0f; i != ROUND_TO_FOUR(this->m_I_ch_nsamples,4ull); i += 4ull,jj += 4.0f) 
                            {
                               const __m128 vt_i{_mm_setr_ps(jj*inv_sr,(jj+1.0f)*inv_sr,(jj+2.0f)*inv_sr,(jj+3.0f)*inv_sr)};
                               const __m128 vcos_val{_mm_cos_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vfc,vt_i),vph0))};
                               const __mmask8 vcos_ge_0{_mm_cmp_ps_mask(vcos_val,vzero,_CMP_GE_OQ)};
                               _mm_store_ps(&this->m_I_ch_bitstream.m_data[i], _mm_mask_blend_ps(vcos_ge_0,vnone,vpone));   
                           }
                       
                           for(j = i; j != this->m_I_ch_nsamples; ++j) 
                           {
                                const float t_j{static_cast<float>(j*inv_sr)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
                                const float cos_val = ceph_cosf(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_j));
                                this->m_I_ch_bitstream.m_data[j] = (cos_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
                                const float cos_val = std::cos(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_j));
                                this->m_I_ch_bitstream.m_data[j] = (cos_val>=0.0f)?+1.0f:-1.0f;
#endif
                           }
                       }
                       else 
                       {
                           for(i = 0ull; i != ROUND_TO_FOUR(this->m_I_ch_nsamples,4ull); i += 4ull) 
                           {
                               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i],_MM_HINT_T0);
                               const __m128 vt_i{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]),vinv_sr)};
                               const __m128 vcos_val{_mm_cos_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vfc,vt_i),vph0))};
                               const __mmask8 vcos_ge_0{_mm_cmp_ps_mask(vcos_val,vzero,_CMP_GE_OQ)};
                               _mm_store_ps(&this->m_I_ch_bitstream.m_data[i], _mm_mask_blend_ps(vcos_ge_0,vnone,vpone));   
                           }
                       
                           for(j = i; j != this->m_I_ch_nsamples; ++j) 
                           {
                                const float t_j{gms::math::LUT_loop_indices_2257_align16[j]*inv_sr};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
                                const float cos_val = ceph_cosf(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_j));
                                this->m_I_ch_bitstream.m_data[j] = (cos_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
                                const float cos_val = std::cos(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_j));
                                this->m_I_ch_bitstream.m_data[j] = (cos_val>=0.0f)?+1.0f:-1.0f;
#endif
                           }
                       }

                       return (0);
                  }

                  __ATTR_ALWAYS_INLINE__
                  inline 
                  std::int32_t 
                  generate_I_channel_bitstream_sse_u4x()
                  {
                       using namespace gms::math;
                       const float inv_sr{1.0f/this->m_I_ch_nsamples};
                       const __m128 vinv_sr{_mm_set1_ps(inv_sr)};
                       const __m128 vfc{_mm_set1_ps(this->m_I_fc)};
                       const __m128 vph0{_mm_set1_ps(this->m_I_ph0)};
                       const __m128 v2pi{_mm_set1_ps(6.283185307179586476925286766559f)};
                       const __m128 vzero{_mm_setzero_ps()};
                       const __m128 vpone{_mm_set1_ps(+1.0f)};
                       const __m128 vnone{_mm_set1_ps(-1.0f)};
                       std::size_t i;
                       float j;
                       constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
                       constexpr std::size_t LUT_loop_idx_threshold{2257ull};
                       if(this->m_I_ch_nsamples>LUT_loop_idx_threshold)
                       {
                             for(i = 0ull,j = 0.0f; (i+15ull) < this->m_I_ch_nsamples; i += 16ull,j += 16.0f) 
                             {
                                  
                                   const __m128 vt_i_0{_mm_setr_ps(j*inv_sr,(j+1.0f)*inv_sr,(j+2.0f)*inv_sr,(j+3.0f)*inv_sr)};
                                   const __m128 vcos_val_0    = _mm_cos_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vfc,vt_i_0),vph0));
                                   const __mmask8 vcos_ge_0_0 = _mm_cmp_ps_mask(vcos_val_0,vzero,_CMP_GE_OQ);
                                   _mm_store_ps(&this->m_I_ch_bitstream.m_data[i+0ull], _mm_mask_blend_ps(vcos_ge_0_0,vnone,vpone));
                                   const __m128 vt_i_1{_mm_setr_ps((j+4.0f)*inv_sr,(j+5.0f)*inv_sr,(j+6.0f)*inv_sr,(j+7.0f)*inv_sr)};
                                   const __m128 vcos_val_1    = _mm_cos_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vfc,vt_i_1),vph0));
                                   const __mmask8 vcos_ge_0_1 = _mm_cmp_ps_mask(vcos_val_1,vzero,_CMP_GE_OQ);
                                   _mm_store_ps(&this->m_I_ch_bitstream.m_data[i+4ull], _mm_mask_blend_ps(vcos_ge_0_1,vnone,vpone));
                                   const __m128 vt_i_2{_mm_setr_ps((j+8.0f)*inv_sr,(j+9.0f)*inv_sr,(j+10.0f)*inv_sr,(j+11.0f)*inv_sr)};
                                   const __m128 vcos_val_2    = _mm_cos_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vfc,vt_i_2),vph0));
                                   const __mmask8 vcos_ge_0_2 = _mm_cmp_ps_mask(vcos_val_2,vzero,_CMP_GE_OQ);
                                   _mm_store_ps(&this->m_I_ch_bitstream.m_data[i+8ull], _mm_mask_blend_ps(vcos_ge_0_2,vnone,vpone));
                                   const __m128 vt_i_3{_mm_setr_ps((j+12.0f)*inv_sr,(j+13.0f)*inv_sr,(j+14.0f)*inv_sr,(j+15.0f)*inv_sr)};
                                   const __m128 vcos_val_3    = _mm_cos_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vfc,vt_i_3),vph0));
                                   const __mmask8 vcos_ge_0_3 = _mm_cmp_ps_mask(vcos_val_3,vzero,_CMP_GE_OQ);
                                   _mm_store_ps(&this->m_I_ch_bitstream.m_data[i+12ull], _mm_mask_blend_ps(vcos_ge_0_3,vnone,vpone));
                              }

                              for(; (i+7ull) < this->m_I_ch_nsamples; i += 8ull,j += 8.0f) 
                              {
                                    const __m128 vt_i_0{_mm_setr_ps(j*inv_sr,(j+1.0f)*inv_sr,(j+2.0f)*inv_sr,(j+3.0f)*inv_sr)};
                                    const __m128 vcos_val_0    = _mm_cos_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vfc,vt_i_0),vph0));
                                    const __mmask8 vcos_ge_0_0 = _mm_cmp_ps_mask(vcos_val_0,vzero,_CMP_GE_OQ);
                                    _mm_store_ps(&this->m_I_ch_bitstream.m_data[i+0ull], _mm_mask_blend_ps(vcos_ge_0_0,vnone,vpone));
                                    const __m128 vt_i_1{_mm_setr_ps((j+4.0f)*inv_sr,(j+5.0f)*inv_sr,(j+6.0f)*inv_sr,(j+7.0f)*inv_sr)};
                                    const __m128 vcos_val_1    = _mm_cos_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vfc,vt_i_1),vph0));
                                    const __mmask8 vcos_ge_0_1 = _mm_cmp_ps_mask(vcos_val_1,vzero,_CMP_GE_OQ);
                                    _mm_store_ps(&this->m_I_ch_bitstream.m_data[i+4ull], _mm_mask_blend_ps(vcos_ge_0_1,vnone,vpone));
                          
                              }

                              for(; (i+3ull) < this->m_I_ch_nsamples; i += 4ull,j += 4.0f) 
                              {
                                    const __m128 vt_i_0{_mm_setr_ps(j*inv_sr,(j+1.0f)*inv_sr,(j+2.0f)*inv_sr,(j+3.0f)*inv_sr)};
                                    const __m128 vcos_val_0    = _mm_cos_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vfc,vt_i_0),vph0));
                                    const __mmask8 vcos_ge_0_0 = _mm_cmp_ps_mask(vcos_val_0,vzero,_CMP_GE_OQ);
                                    _mm_store_ps(&this->m_I_ch_bitstream.m_data[i+0ull], _mm_mask_blend_ps(vcos_ge_0_0,vnone,vpone));
                               }
                                            
                              for(; (i+0ull) < this->m_I_ch_nsamples; i += 1ull,j += 1.0f) 
                              {
                                      const float t_i{j*inv_sr};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
                                      const float cos_val = ceph_cosf(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_i));
                                      this->m_I_ch_bitstream.m_data[i] = (cos_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
                                      const float cos_val = std::cos(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_i));
                                      this->m_I_ch_bitstream.m_data[i] = (cos_val>=0.0f)?+1.0f:-1.0f;
#endif
                               }
                       }
                       else 
                       {
                       
                            for(i = 0ull; (i+15ull) < this->m_I_ch_nsamples; i += 16ull) 
                            {
                                 
                                 _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+0ull],_MM_HINT_T0);
                                 const __m128 vt_i_0{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]),vinv_sr)};
                                 const __m128 vcos_val_0    = _mm_cos_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vfc,vt_i_0),vph0));
                                 const __mmask8 vcos_ge_0_0 = _mm_cmp_ps_mask(vcos_val_0,vzero,_CMP_GE_OQ);
                                 _mm_store_ps(&this->m_I_ch_bitstream.m_data[i+0ull], _mm_mask_blend_ps(vcos_ge_0_0,vnone,vpone));
                                                                 
                                 const __m128 vt_i_1{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+4ull]),vinv_sr)};
                                 const __m128 vcos_val_1    = _mm_cos_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vfc,vt_i_1),vph0));
                                 const __mmask8 vcos_ge_0_1 = _mm_cmp_ps_mask(vcos_val_1,vzero,_CMP_GE_OQ);
                                 _mm_store_ps(&this->m_I_ch_bitstream.m_data[i+4ull], _mm_mask_blend_ps(vcos_ge_0_1,vnone,vpone));
                                                               
                                 const __m128 vt_i_2{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+8ull]),vinv_sr)};
                                 const __m128 vcos_val_2    = _mm_cos_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vfc,vt_i_2),vph0));
                                 const __mmask8 vcos_ge_0_2 = _mm_cmp_ps_mask(vcos_val_2,vzero,_CMP_GE_OQ);
                                 _mm_store_ps(&this->m_I_ch_bitstream.m_data[i+8ull], _mm_mask_blend_ps(vcos_ge_0_2,vnone,vpone));
                                                                 
                                 const __m128 vt_i_3{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+12ull]),vinv_sr)};
                                 const __m128 vcos_val_3    = _mm_cos_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vfc,vt_i_3),vph0));
                                 const __mmask8 vcos_ge_0_3 = _mm_cmp_ps_mask(vcos_val_3,vzero,_CMP_GE_OQ);
                                 _mm_store_ps(&this->m_I_ch_bitstream.m_data[i+12ull], _mm_mask_blend_ps(vcos_ge_0_3,vnone,vpone));
                           }

                           for(; (i+7ull) < this->m_I_ch_nsamples; i += 8ull) 
                           {
                          
                                 const __m128 vt_i_0{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]),vinv_sr)};
                                 const __m128 vcos_val_0    = _mm_cos_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vfc,vt_i_0),vph0));
                                 const __mmask8 vcos_ge_0_0 = _mm_cmp_ps_mask(vcos_val_0,vzero,_CMP_GE_OQ);
                                 _mm_store_ps(&this->m_I_ch_bitstream.m_data[i+0ull], _mm_mask_blend_ps(vcos_ge_0_0,vnone,vpone));
                          
                                 const __m128 vt_i_1{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+4ull]),vinv_sr)};
                                 const __m128 vcos_val_1    = _mm_cos_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vfc,vt_i_1),vph0));
                                 const __mmask8 vcos_ge_0_1 = _mm_cmp_ps_mask(vcos_val_1,vzero,_CMP_GE_OQ);
                                 _mm_store_ps(&this->m_I_ch_bitstream.m_data[i+4ull], _mm_mask_blend_ps(vcos_ge_0_1,vnone,vpone));
                          
                           }

                           for(; (i+3ull) < this->m_I_ch_nsamples; i += 4ull) 
                           {
                               const __m128 vt_i_0{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]),vinv_sr)};
                               const __m128 vcos_val_0    = _mm_cos_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vfc,vt_i_0),vph0));
                               const __mmask8 vcos_ge_0_0 = _mm_cmp_ps_mask(vcos_val_0,vzero,_CMP_GE_OQ);
                               _mm_store_ps(&this->m_I_ch_bitstream.m_data[i+0ull], _mm_mask_blend_ps(vcos_ge_0_0,vnone,vpone));
                           }
                      
                       
                           for(; (i+0ull) < this->m_I_ch_nsamples; i += 1ull) 
                           {
                            
                                    const float t_i{gms::math::LUT_loop_indices_2257_align16[i]*inv_sr};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
                                    const float cos_val = ceph_cosf(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_i));
                                    this->m_I_ch_bitstream.m_data[i] = (cos_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
                                    const float cos_val = std::cos(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_i));
                                    this->m_I_ch_bitstream.m_data[i] = (cos_val>=0.0f)?+1.0f:-1.0f;
#endif
                           }
                       }
                       return (0);
                  }

                   __ATTR_ALWAYS_INLINE__
                  inline 
                  std::int32_t 
                  generate_Q_channel_bitstream_sse()                                                   
                  {
                       using namespace gms::math;
                       const float inv_sr{1.0f/this->m_Q_ch_nsamples};
                       const __m128 vinv_sr{_mm_set1_ps(inv_sr)};
                       const __m128 vw0{_mm_set1_ps(this->m_Q_fc)};
                       const __m128 vph0{_mm_set1_ps(this->m_Q_ph0)};
                       const __m128 v2pi{_mm_set1_ps(6.283185307179586476925286766559f)};
                       const __m128 vzero{_mm_setzero_ps()};
                       const __m128 vpone{_mm_set1_ps(+1.0f)};
                       const __m128 vnone{_mm_set1_ps(-1.0f)};
                       std::size_t i,j;
                       float jj;
                       constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
                       constexpr std::size_t LUT_loop_idx_threshold{2257ull};
                       if(this->m_Q_ch_nsamples>LUT_loop_idx_threshold)
                       {
                            for(i = 0ull,jj = 0.0f; i != ROUND_TO_FOUR(this->m_Q_ch_nsamples,4ull); i += 4ull,jj += 4.0f) 
                            {
                                const __m128 vt_i{_mm_setr_ps(jj*inv_sr,(jj+1.0f)*inv_sr,(jj+2.0f)*inv_sr,(jj+3.0f)*inv_sr)};
                                const __m128 vsin_val    = _mm_sin_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vw0,vt_i),vph0));
                                const __mmask8 vsin_ge_0 = _mm_cmp_ps_mask(vsin_val,vzero,_CMP_GE_OQ);
                                _mm_store_ps(&this->m_Q_ch_bitstream.m_data[i], _mm_mask_blend_ps(vsin_ge_0,vnone,vpone));   
                            }
                       
                            for(j = i; j != this->m_Q_ch_nsamples; ++j) 
                            {
                                  const float t_j{static_cast<float>(j*inv_sr)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
                                  const float sin_val = ceph_sinf(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_j));
                                  this->m_Q_ch_bitstream.m_data[j] = (sin_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
                                  const float sin_val = std::sin(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_j));
                                  this->m_Q_ch_bitstream.m_data[j] = (sin_val>=0.0f)?+1.0f:-1.0f;
#endif
                            }
                       }
                       else 
                       {
                             for(i = 0ull; i != ROUND_TO_FOUR(this->m_Q_ch_nsamples,4ull); i += 4ull) 
                             {
                                  
                                 _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i],_MM_HINT_T0);
                                 const __m128 vt_i{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i]),vinv_sr)};
                                 const __m128 vsin_val    = _mm_sin_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vw0,vt_i),vph0));
                                 const __mmask8 vsin_ge_0 = _mm_cmp_ps_mask(vsin_val,vzero,_CMP_GE_OQ);
                                 _mm_store_ps(&this->m_Q_ch_bitstream.m_data[i], _mm_mask_blend_ps(vsin_ge_0,vnone,vpone));   
                             }
                       
                             for(j = i; j != this->m_Q_ch_nsamples; ++j) 
                             {
                                   //const float t_j{static_cast<float>(j*inv_sr)};
                                   const float t_j{gms::math::LUT_loop_indices_2257_align16[j]*inv_sr};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
                                   const float sin_val = ceph_sinf(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_j));
                                   this->m_Q_ch_bitstream.m_data[j] = (sin_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
                                   const float sin_val = std::sin(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_j));
                                   this->m_Q_ch_bitstream.m_data[j] = (sin_val>=0.0f)?+1.0f:-1.0f;
#endif
                             }
                       }
                      
                      
                       return (0);
                  }

                  __ATTR_ALWAYS_INLINE__
                  inline 
                  std::int32_t 
                  generate_Q_channel_bitstream_sse_u4x()                                                   
                  {
                       using namespace gms::math;
                       const float inv_sr{1.0f/this->m_Q_ch_nsamples};
                       const __m128 vinv_sr{_mm_set1_ps(inv_sr)};
                       const __m128 vw0{_mm_set1_ps(this->m_Q_fc)};
                       const __m128 vph0{_mm_set1_ps(this->m_Q_ph0)};
                       const __m128 v2pi{_mm_set1_ps(6.283185307179586476925286766559f)};
                       const __m128 vzero{_mm_setzero_ps()};
                       const __m128 vpone{_mm_set1_ps(+1.0f)};
                       const __m128 vnone{_mm_set1_ps(-1.0f)};
                       std::size_t i;
                       float j;
                       constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
                       constexpr std::size_t LUT_loop_idx_threshold{2257ull};
                       if(this->m_Q_ch_nsamples>LUT_loop_idx_threshold)
                       {
                             for(i = 0ull,j = 0.0f; (i+15ull) < this->m_Q_ch_nsamples; i += 16ull,j += 16.0f) 
                             {
                                   const __m128 vt_i_0{_mm_setr_ps(j*inv_sr,(j+1.0f)*inv_sr,(j+2.0f)*inv_sr,(j+3.0f)*inv_sr)};
                                   const __m128 vsin_val_0    = _mm_sin_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vw0,vt_i_0),vph0));
                                   const __mmask8 vsin_ge_0_0 = _mm_cmp_ps_mask(vsin_val_0,vzero,_CMP_GE_OQ);
                                   _mm_store_ps(&this->m_Q_ch_bitstream.m_data[i+0ull], _mm_mask_blend_ps(vsin_ge_0_0,vnone,vpone));
                                   const __m128 vt_i_1{_mm_setr_ps((j+4.0f)*inv_sr,(j+5.0f)*inv_sr,(j+6.0f)*inv_sr,(j+7.0f)*inv_sr)};
                                   const __m128 vsin_val_1    = _mm_sin_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vw0,vt_i_1),vph0));
                                   const __mmask8 vsin_ge_0_1 = _mm_cmp_ps_mask(vsin_val_1,vzero,_CMP_GE_OQ);
                                   _mm_store_ps(&this->m_Q_ch_bitstream.m_data[i+4ull], _mm_mask_blend_ps(vsin_ge_0_1,vnone,vpone));
                                   const __m128 vt_i_2{_mm_setr_ps((j+8.0f)*inv_sr,(j+9.0f)*inv_sr,(j+10.0f)*inv_sr,(j+11.0f)*inv_sr)};
                                   const __m128 vsin_val_2    = _mm_sin_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vw0,vt_i_2),vph0));
                                   const __mmask8 vsin_ge_0_2 = _mm_cmp_ps_mask(vsin_val_2,vzero,_CMP_GE_OQ);
                                   _mm_store_ps(&this->m_Q_ch_bitstream.m_data[i+8ull], _mm_mask_blend_ps(vsin_ge_0_2,vnone,vpone));
                                   const __m128 vt_i_3{_mm_setr_ps((j+12.0f)*inv_sr,(j+13.0f)*inv_sr,(j+14.0f)*inv_sr,(j+15.0f)*inv_sr)};
                                   const __m128 vsin_val_3    = _mm_sin_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vw0,vt_i_3),vph0));
                                   const __mmask8 vsin_ge_0_3 = _mm_cmp_ps_mask(vsin_val_3,vzero,_CMP_GE_OQ);
                                   _mm_store_ps(&this->m_Q_ch_bitstream.m_data[i+12ull], _mm_mask_blend_ps(vsin_ge_0_3,vnone,vpone));
                             }

                             for(; (i+7ul) < this->m_Q_ch_nsamples; i += 8ull, j += 8.0f) 
                             {
                                   const __m128 vt_i_0{_mm_setr_ps(j*inv_sr,(j+1.0f)*inv_sr,(j+2.0f)*inv_sr,(j+3.0f)*inv_sr)};
                                   const __m128 vsin_val_0    = _mm_sin_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vw0,vt_i_0),vph0));
                                   const __mmask8 vsin_ge_0_0 = _mm_cmp_ps_mask(vsin_val_0,vzero,_CMP_GE_OQ);
                                   _mm_store_ps(&this->m_Q_ch_bitstream.m_data[i+0ull], _mm_mask_blend_ps(vsin_ge_0_0,vnone,vpone));
                                   const __m128 vt_i_1{_mm_setr_ps((j+4.0f)*inv_sr,(j+5.0f)*inv_sr,(j+6.0f)*inv_sr,(j+7.0f)*inv_sr)};
                                   const __m128 vsin_val_1    = _mm_sin_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vw0,vt_i_1),vph0));
                                   const __mmask8 vsin_ge_0_1 = _mm_cmp_ps_mask(vsin_val_1,vzero,_CMP_GE_OQ);
                                   _mm_store_ps(&this->m_Q_ch_bitstream.m_data[i+4ull], _mm_mask_blend_ps(vsin_ge_0_1,vnone,vpone));
                          
                             }

                             for(; (i+3ull) < this->m_Q_ch_nsamples; i += 4ull, j += 4.0f) 
                             {
                                   const __m128 vt_i_0{_mm_setr_ps(j*inv_sr,(j+1.0f)*inv_sr,(j+2.0f)*inv_sr,(j+3.0f)*inv_sr)};
                                   const __m128 vsin_val_0    = _mm_sin_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vw0,vt_i_0),vph0));
                                   const __mmask8 vsin_ge_0_0 = _mm_cmp_ps_mask(vsin_val_0,vzero,_CMP_GE_OQ);
                                   _mm_store_ps(&this->m_Q_ch_bitstream.m_data[i+0ull], _mm_mask_blend_ps(vsin_ge_0_0,vnone,vpone));
                             }
                      
                       
                             for(; (i+0ull) < this->m_Q_ch_nsamples; i += 1ull, j += 1.0f) 
                             {
                                   const float t_i{j*inv_sr};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
                                   const float sin_val = ceph_sinf(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_i));
                                   this->m_Q_ch_bitstream.m_data[i] = (sin_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
                                   const float sin_val = std::sin(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_i));
                                   this->m_Q_ch_bitstream.m_data[i] = (sin_val>=0.0f)?+1.0f:-1.0f;
#endif
                             }
                       }
                       else 
                       {
                             for(i = 0ull; (i+15ull) < this->m_Q_ch_nsamples; i += 16ull) 
                             {
                                  
                                   _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i],_MM_HINT_T0);
                                   const __m128 vt_i_0{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]),vinv_sr)};
                                   const __m128 vsin_val_0    = _mm_sin_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vw0,vt_i_0),vph0));
                                   const __mmask8 vsin_ge_0_0 = _mm_cmp_ps_mask(vsin_val_0,vzero,_CMP_GE_OQ);
                                   _mm_store_ps(&this->m_Q_ch_bitstream.m_data[i+0ull], _mm_mask_blend_ps(vsin_ge_0_0,vnone,vpone));
                               
                                   const __m128 vt_i_1{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+4ull]),vinv_sr)}; 
                                   const __m128 vsin_val_1    = _mm_sin_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vw0,vt_i_1),vph0));
                                   const __mmask8 vsin_ge_0_1 = _mm_cmp_ps_mask(vsin_val_1,vzero,_CMP_GE_OQ);
                                   _mm_store_ps(&this->m_Q_ch_bitstream.m_data[i+4ull], _mm_mask_blend_ps(vsin_ge_0_1,vnone,vpone));
                                                                      
                                   const __m128 vt_i_2{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+8ull]),vinv_sr)};
                                   const __m128 vsin_val_2    = _mm_sin_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vw0,vt_i_2),vph0));
                                   const __mmask8 vsin_ge_0_2 = _mm_cmp_ps_mask(vsin_val_2,vzero,_CMP_GE_OQ);
                                   _mm_store_ps(&this->m_Q_ch_bitstream.m_data[i+8ull], _mm_mask_blend_ps(vsin_ge_0_2,vnone,vpone));
                                   
                                   const __m128 vt_i_3{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+12ull]),vinv_sr)};
                                   const __m128 vsin_val_3    = _mm_sin_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vw0,vt_i_3),vph0));
                                   const __mmask8 vsin_ge_0_3 = _mm_cmp_ps_mask(vsin_val_3,vzero,_CMP_GE_OQ);
                                   _mm_store_ps(&this->m_Q_ch_bitstream.m_data[i+12ull], _mm_mask_blend_ps(vsin_ge_0_3,vnone,vpone));
                              }

                              for(; (i+7ul) < this->m_Q_ch_nsamples; i += 8ull) 
                              {
                                     _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i],_MM_HINT_T0);
                                   const __m128 vt_i_0{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]),vinv_sr)};
                                   const __m128 vsin_val_0    = _mm_sin_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vw0,vt_i_0),vph0));
                                   const __mmask8 vsin_ge_0_0 = _mm_cmp_ps_mask(vsin_val_0,vzero,_CMP_GE_OQ);
                                   _mm_store_ps(&this->m_Q_ch_bitstream.m_data[i+0ull], _mm_mask_blend_ps(vsin_ge_0_0,vnone,vpone));
                                   
                                   const __m128 vt_i_1{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+4ull]),vinv_sr)}; 
                                   const __m128 vsin_val_1    = _mm_sin_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vw0,vt_i_1),vph0));
                                   const __mmask8 vsin_ge_0_1 = _mm_cmp_ps_mask(vsin_val_1,vzero,_CMP_GE_OQ);
                                   _mm_store_ps(&this->m_Q_ch_bitstream.m_data[i+4ull], _mm_mask_blend_ps(vsin_ge_0_1,vnone,vpone));
                                   
                              }

                              for(; (i+3ull) < this->m_Q_ch_nsamples; i += 4ull) 
                              {
                                   _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align16[i+0ull],_MM_HINT_T0);
                                   const __m128 vt_i_0{_mm_mul_ps(_mm_load_ps(&gms::math::LUT_loop_indices_2257_align16[i+0ull]),vinv_sr)};
                                   const __m128 vsin_val_0    = _mm_sin_ps(_mm_fmadd_ps(v2pi,_mm_mul_ps(vw0,vt_i_0),vph0));
                                   const __mmask8 vsin_ge_0_0 = _mm_cmp_ps_mask(vsin_val_0,vzero,_CMP_GE_OQ);
                                   _mm_store_ps(&this->m_Q_ch_bitstream.m_data[i+0ull], _mm_mask_blend_ps(vsin_ge_0_0,vnone,vpone));
                              }
                      
                              for(; (i+0ull) < this->m_Q_ch_nsamples; i += 1ull) 
                              {
                                  
                                   const float t_i{gms::math::LUT_loop_indices_2257_align16[i]*inv_sr};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
                                   const float sin_val = ceph_sinf(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_i));
                                   this->m_Q_ch_bitstream.m_data[i] = (sin_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
                                   const float sin_val = std::sin(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_i));
                                   this->m_Q_ch_bitstream.m_data[i] = (sin_val>=0.0f)?+1.0f:-1.0f;
#endif
                              }
                       }
                      
                      
                       return (0);
                  }

                   __ATTR_ALWAYS_INLINE__
                  inline 
                  std::int32_t 
                  generate_I_channel_bitstream_avx(const bool do_constants_prefetch)                                                    
                  {
                       using namespace gms::math;
                       const float inv_sr{1.0f/this->m_I_ch_nsamples};
                       const __m256 vw0{_mm256_set1_ps(this->m_I_fc)};
                       const __m256 vph0{_mm256_set1_ps(this->m_I_ph0)};
                       const __m256 vinv_sr{_mm256_set1_ps(inv_sr)};
                       __m256 c0_8,v2pi;
                       __m256 vzero,vpone;
                       __m256 vnone;
                       if(!do_constants_prefetch)
                       {
                          c0_8 = _mm256_setr_ps(0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f);
                          v2pi = _mm256_set1_ps(6.283185307179586476925286766559f);
                          vzero= _mm256_setzero_ps();
                          vpone= _mm256_set1_ps(+1.0f);
                          vnone= _mm256_set1_ps(-1.0f);
                       }
                       else 
                       {
                          __ATTR_ALIGN__(32) const float prefetched_constants[48] = {
                                                     0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f,
                                                     6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                     6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                     6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                     6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                     0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,
                                                     +1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,
                                                     -1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,
                                                     0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f}; // the last eight floats serve as a padding
                              _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[32],_MM_HINT_T0);
                              c0_8 = _mm256_load_ps(&prefetched_constants[0]);
                              v2pi = _mm256_load_ps(&prefetched_constants[8]);
                              vzero= _mm256_load_ps(&prefetched_constants[16]);
                              vpone= _mm256_load_ps(&prefetched_constants[24]);
                              vnone= _mm256_load_ps(&prefetched_constants[32]);
                       }
                      
                       std::size_t i,j;
                       float jj;
                       constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
                       constexpr std::size_t LUT_loop_idx_threshold{2257ull};
                       if(this->m_I_ch_nsamples>LUT_loop_idx_threshold)
                       {
                            for(i = 0ull,jj = 0.0f; i != ROUND_TO_EIGHT(this->m_I_ch_nsamples,8ull); i += 8ull,jj += 8.0f) 
                            {
                               const __m256 vt_i{_mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(jj),c0_8),vinv_sr)};
                               const __m256 vcos_val    = _mm256_cos_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i),vph0));
                               const __mmask8 vcos_ge_0 = _mm256_cmp_ps_mask(vcos_val,vzero,_CMP_GE_OQ);
                               _mm256_store_ps(&this->m_I_ch_bitstream.m_data[i], _mm256_mask_blend_ps(vcos_ge_0,vnone,vpone));   
                           }
                       
                           for(j = i; j != this->m_I_ch_nsamples; ++j) 
                           {
                                const float t_j{static_cast<float>(j*inv_sr)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
                                const float cos_val = ceph_cosf(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_j));
                                this->m_I_ch_bitstream.m_data[j] = (cos_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
                                const float cos_val = std::cos(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_j));
                                this->m_I_ch_bitstream.m_data[j] = (cos_val>=0.0f)?+1.0f:-1.0f;
#endif
                           }
                       }
                       else 
                       {
                       
                           for(i = 0ull; i != ROUND_TO_EIGHT(this->m_I_ch_nsamples,8ull); i += 8ull) 
                           {
                               
                               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align32[i],_MM_HINT_T0);
                               const __m256 vt_i{_mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i]),vinv_sr)};
                               const __m256 vcos_val    = _mm256_cos_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i),vph0));
                               const __mmask8 vcos_ge_0 = _mm256_cmp_ps_mask(vcos_val,vzero,_CMP_GE_OQ);
                               _mm256_store_ps(&this->m_I_ch_bitstream.m_data[i], _mm256_mask_blend_ps(vcos_ge_0,vnone,vpone));   
                           }
                       
                           for(j = i; j != this->m_I_ch_nsamples; ++j) 
                           {
                                const float t_j{gms::math::LUT_loop_indices_2257_align32[j]*inv_sr};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
                                const float cos_val = ceph_cosf(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_j));
                                this->m_I_ch_bitstream.m_data[j] = (cos_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
                                const float cos_val = std::cos(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_j));
                                this->m_I_ch_bitstream.m_data[j] = (cos_val>=0.0f)?+1.0f:-1.0f;
#endif
                           }

                       }
                       return (0);
                  }

                   __ATTR_ALWAYS_INLINE__
                  inline 
                  std::int32_t 
                  generate_I_channel_bitstream_avx_u4x(const bool do_constants_prefetch)                                                  
                  {
                       using namespace gms::math;
                       
                       const float inv_sr{1.0f/this->m_I_ch_nsamples};
                       const __m256 vinv_sr{_mm256_set1_ps(inv_sr)};
                       const __m256 vw0{_mm256_set1_ps(this->m_I_fc)};
                       const __m256 vph0{_mm256_set1_ps(this->m_I_ph0)};
                       __m256 c0_7,c8_15,c16_23,c24_31;
                       __m256 v2pi,vzero,vpone,vnone;
                       if(!do_constants_prefetch)
                       {
                           c0_7   = _mm256_setr_ps(0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f);
                           c8_15  = _mm256_setr_ps(8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f);
                           c16_23 = _mm256_setr_ps(16.0f,17.0f,18.0f,19.0f,20.0f,21.0f,22.0f,23.0f);
                           c24_31 = _mm256_setr_ps(24.0f,25.0f,26.0f,27.0f,28.0f,29.0f,30.0f,31.0f);
                           v2pi   = _mm256_set1_ps(6.283185307179586476925286766559f);
                           vzero  = _mm256_setzero_ps();
                           vpone  = _mm256_set1_ps(+1.0f);
                           vnone  = _mm256_set1_ps(-1.0f);  
                       }
                       else 
                       {
                              __ATTR_ALIGN__(32) const float prefetched_constants[64] = {
                                                      0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f,          
                                                      8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f,
                                                      16.0f,17.0f,18.0f,19.0f,20.0f,21.0f,22.0f,23.0f,
                                                      24.0f,25.0f,26.0f,27.0f,28.0f,29.0f,30.0f,31.0f,
                                                      6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                      6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                      6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                      6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                      0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,
                                                      +1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,
                                                      -1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f};
                              _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[48],_MM_HINT_T0);
                              c0_7   = _mm256_load_ps(&prefetched_constants[0]);
                              c8_15  = _mm256_load_ps(&prefetched_constants[8]);
                              c16_23 = _mm256_load_ps(&prefetched_constants[16]);
                              c24_31 = _mm256_load_ps(&prefetched_constants[24]);
                              v2pi   = _mm256_load_ps(&prefetched_constants[32]);
                              vzero  = _mm256_load_ps(&prefetched_constants[40]);
                              vpone  = _mm256_load_ps(&prefetched_constants[48]);
                              vnone  = _mm256_load_ps(&prefetched_constants[56]);
                       }
                       std::size_t i;
                       float j;
                       constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
                       constexpr std::size_t LUT_loop_idx_threshold{2257ull};
                       if(this->m_I_ch_nsamples>LUT_loop_idx_threshold)
                       {
                             for(i = 0ull,j = 0.0f; (i+31ull) < this->m_I_ch_nsamples; i += 32ull,j += 32.0f) 
                             {
                                   const __m256 vt_i_0{_mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c0_7),vinv_sr)};
                                   const __m256 vcos_val_0    = _mm256_cos_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i_0),vph0));
                                   const __mmask8 vcos_ge_0_0 = _mm256_cmp_ps_mask(vcos_val_0,vzero,_CMP_GE_OQ);
                                   _mm256_store_ps(&this->m_I_ch_bitstream.m_data[i+0ull], _mm256_mask_blend_ps(vcos_ge_0_0,vnone,vpone));
                                   const __m256 vt_i_1{_mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c8_15),vinv_sr)};
                                   const __m256 vcos_val_1    = _mm256_cos_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i_1),vph0));
                                   const __mmask8 vcos_ge_0_1 = _mm256_cmp_ps_mask(vcos_val_1,vzero,_CMP_GE_OQ);
                                   _mm256_store_ps(&this->m_I_ch_bitstream.m_data[i+8ull], _mm256_mask_blend_ps(vcos_ge_0_1,vnone,vpone));
                                   const __m256 vt_i_2{_mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c16_23),vinv_sr)};
                                   const __m256 vcos_val_2    = _mm256_cos_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i_2),vph0));
                                   const __mmask8 vcos_ge_0_2 = _mm256_cmp_ps_mask(vcos_val_2,vzero,_CMP_GE_OQ);
                                   _mm256_store_ps(&this->m_I_ch_bitstream.m_data[i+16ull], _mm256_mask_blend_ps(vcos_ge_0_2,vnone,vpone));
                                   const __m256 vt_i_3{_mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c24_31),vinv_sr)};
                                   const __m256 vcos_val_3    = _mm256_cos_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i_3),vph0));
                                   const __mmask8 vcos_ge_0_3 = _mm256_cmp_ps_mask(vcos_val_3,vzero,_CMP_GE_OQ);
                                   _mm256_store_ps(&this->m_I_ch_bitstream.m_data[i+24ull], _mm256_mask_blend_ps(vcos_ge_0_3,vnone,vpone));
                              }

                              for(; (i+15ull) < this->m_I_ch_nsamples; i += 16ull,j += 16.0f) 
                              {  
                                   const __m256 vt_i_0{_mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c0_7),vinv_sr)};
                                   const __m256 vcos_val_0    = _mm256_cos_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i_0),vph0));
                                   const __mmask8 vcos_ge_0_0 = _mm256_cmp_ps_mask(vcos_val_0,vzero,_CMP_GE_OQ);
                                   _mm256_store_ps(&this->m_I_ch_bitstream.m_data[i+0ull], _mm256_mask_blend_ps(vcos_ge_0_0,vnone,vpone));
                                   const __m256 vt_i_1{_mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c8_15),vinv_sr)};
                                   const __m256 vcos_val_1    = _mm256_cos_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i_1),vph0));
                                   const __mmask8 vcos_ge_0_1 = _mm256_cmp_ps_mask(vcos_val_1,vzero,_CMP_GE_OQ);
                                   _mm256_store_ps(&this->m_I_ch_bitstream.m_data[i+8ull], _mm256_mask_blend_ps(vcos_ge_0_1,vnone,vpone));
                              }

                              for(; (i+7ull) < this->m_I_ch_nsamples; i += 8ull,j += 8.0f) 
                              {
                                   const __m256 vt_i_0{_mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c0_7),vinv_sr)};
                                   const __m256 vcos_val_0    = _mm256_cos_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i_0),vph0));
                                   const __mmask8 vcos_ge_0_0 = _mm256_cmp_ps_mask(vcos_val_0,vzero,_CMP_GE_OQ);
                                   _mm256_store_ps(&this->m_I_ch_bitstream.m_data[i+0ull], _mm256_mask_blend_ps(vcos_ge_0_0,vnone,vpone));
                              }
                                            
                              for(; (i+0ull) < this->m_I_ch_nsamples; i += 1ull,j += 1.0f) 
                              {
                                      const float t_i{j*inv_sr};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
                                      const float cos_val = ceph_cosf(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_i));
                                      this->m_I_ch_bitstream.m_data[i] = (cos_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
                                      const float cos_val = std::cos(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_i));
                                      this->m_I_ch_bitstream.m_data[i] = (cos_val>=0.0f)?+1.0f:-1.0f;
#endif
                               }
                       }
                       else 
                       {
                       
                            for(i = 0ull; (i+31ull) < this->m_I_ch_nsamples; i += 32ull) 
                            {
                                
                                 _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align32[i+0ull],_MM_HINT_T0);
                                 const __m256 vt_i_0{_mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i+0ull]),vinv_sr)};
                                 const __m256 vcos_val_0    = _mm256_cos_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i_0),vph0));
                                 const __mmask8 vcos_ge_0_0 = _mm256_cmp_ps_mask(vcos_val_0,vzero,_CMP_GE_OQ);
                                 _mm256_store_ps(&this->m_I_ch_bitstream.m_data[i+0ull], _mm256_mask_blend_ps(vcos_ge_0_0,vnone,vpone));
                                
                                const __m256 vt_i_1{_mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i+8ull]),vinv_sr)};
                                 const __m256 vcos_val_1    = _mm256_cos_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i_1),vph0));
                                 const __mmask8 vcos_ge_0_1 = _mm256_cmp_ps_mask(vcos_val_1,vzero,_CMP_GE_OQ);
                                 _mm256_store_ps(&this->m_I_ch_bitstream.m_data[i+8ull], _mm256_mask_blend_ps(vcos_ge_0_1,vnone,vpone));
                                
                                 _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align32[i+16ull],_MM_HINT_T0);
                                 const __m256 vt_i_2{_mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i+8ull]),vinv_sr)};
                                 const __m256 vcos_val_2    = _mm256_cos_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i_2),vph0));
                                 const __mmask8 vcos_ge_0_2 = _mm256_cmp_ps_mask(vcos_val_2,vzero,_CMP_GE_OQ);
                                 _mm256_store_ps(&this->m_I_ch_bitstream.m_data[i+16ull], _mm256_mask_blend_ps(vcos_ge_0_2,vnone,vpone));
                               
                                 const __m256 vt_i_3{_mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i+24ull]),vinv_sr)};
                                 const __m256 vcos_val_3    = _mm256_cos_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i_3),vph0));
                                 const __mmask8 vcos_ge_0_3 = _mm256_cmp_ps_mask(vcos_val_3,vzero,_CMP_GE_OQ);
                                 _mm256_store_ps(&this->m_I_ch_bitstream.m_data[i+24ull], _mm256_mask_blend_ps(vcos_ge_0_3,vnone,vpone));
                           }

                           for(; (i+15ull) < this->m_I_ch_nsamples; i += 16ull) 
                           {         
                                _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align32[i+0ull],_MM_HINT_T0);
                                 const __m256 vt_i_0{_mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i+0ull]),vinv_sr)};
                                 const __m256 vcos_val_0    = _mm256_cos_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i_0),vph0));
                                 const __mmask8 vcos_ge_0_0 = _mm256_cmp_ps_mask(vcos_val_0,vzero,_CMP_GE_OQ);
                                 _mm256_store_ps(&this->m_I_ch_bitstream.m_data[i+0ull], _mm256_mask_blend_ps(vcos_ge_0_0,vnone,vpone));
                                
                                 const __m256 vt_i_1{_mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i+8ull]),vinv_sr)};
                                 const __m256 vcos_val_1    = _mm256_cos_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i_1),vph0));
                                 const __mmask8 vcos_ge_0_1 = _mm256_cmp_ps_mask(vcos_val_1,vzero,_CMP_GE_OQ);
                                 _mm256_store_ps(&this->m_I_ch_bitstream.m_data[i+8ull], _mm256_mask_blend_ps(vcos_ge_0_1,vnone,vpone));
                          
                           }

                           for(; (i+7ull) < this->m_I_ch_nsamples; i += 8ull) 
                           {
                                 _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align32[i],_MM_HINT_T0);
                                 const __m256 vt_i_0{_mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i+0ull]),vinv_sr)};
                                 const __m256 vcos_val_0    = _mm256_cos_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i_0),vph0));
                                 const __mmask8 vcos_ge_0_0 = _mm256_cmp_ps_mask(vcos_val_0,vzero,_CMP_GE_OQ);
                                 _mm256_store_ps(&this->m_I_ch_bitstream.m_data[i+0ull], _mm256_mask_blend_ps(vcos_ge_0_0,vnone,vpone));
                           }
                      
                       
                           for(; (i+0ull) < this->m_I_ch_nsamples; i += 1ull) 
                           {
                            
                                    const float t_i{gms::math::LUT_loop_indices_2257_align32[i]*inv_sr};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
                                    const float cos_val = ceph_cosf(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_i));
                                    this->m_I_ch_bitstream.m_data[i] = (cos_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
                                    const float cos_val = std::cos(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_i));
                                    this->m_I_ch_bitstream.m_data[i] = (cos_val>=0.0f)?+1.0f:-1.0f;
#endif
                           }
                       }
                       return (0);
                  }

                   __ATTR_ALWAYS_INLINE__
                  inline 
                  std::int32_t 
                  generate_Q_channel_bitstream_avx(const bool do_constants_prefetch)                                                    
                  {
                       using namespace gms::math;
                       const float inv_sr{1.0f/this->m_Q_ch_nsamples};
                       const __m256 vw0{_mm256_set1_ps(this->m_Q_fc)};
                       const __m256 vph0{_mm256_set1_ps(this->m_Q_ph0)};
                       const __m256 vinv_sr{_mm256_set1_ps(inv_sr)};
                       __m256 c0_8,v2pi;
                       __m256 vzero,vpone;
                       __m256 vnone;
                       if(!do_constants_prefetch)
                       {
                          c0_8 = _mm256_setr_ps(0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f);
                          v2pi = _mm256_set1_ps(6.283185307179586476925286766559f);
                          vzero= _mm256_setzero_ps();
                          vpone= _mm256_set1_ps(+1.0f);
                          vnone= _mm256_set1_ps(-1.0f);
                       }
                       else 
                       {
                          __ATTR_ALIGN__(32) const float prefetched_constants[48] = {
                                                     0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f,
                                                     6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                     6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                     6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                     6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                     0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,
                                                     +1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,
                                                     -1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,
                                                     0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f}; // the last eight floats serve as a padding
                              _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[32],_MM_HINT_T0);
                              c0_8 = _mm256_load_ps(&prefetched_constants[0]);
                              v2pi = _mm256_load_ps(&prefetched_constants[8]);
                              vzero= _mm256_load_ps(&prefetched_constants[16]);
                              vpone= _mm256_load_ps(&prefetched_constants[24]);
                              vnone= _mm256_load_ps(&prefetched_constants[32]);
                       }
                      
                       std::size_t i,j;
                       float jj;
                       constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
                       constexpr std::size_t LUT_loop_idx_threshold{2257ull};
                       if(this->m_Q_ch_nsamples>LUT_loop_idx_threshold)
                       {
                            for(i = 0ull,jj = 0.0f; i != ROUND_TO_EIGHT(this->m_Q_ch_nsamples,8ull); i += 8ull,jj += 8.0f) 
                            {
                               const __m256 vt_i{_mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(jj),c0_8),vinv_sr)};
                               const __m256 vsin_val    = _mm256_sin_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i),vph0));
                               const __mmask8 vsin_ge_0 = _mm256_cmp_ps_mask(vsin_val,vzero,_CMP_GE_OQ);
                               _mm256_store_ps(&this->m_Q_ch_bitstream.m_data[i], _mm256_mask_blend_ps(vsin_ge_0,vnone,vpone));   
                           }
                       
                           for(j = i; j != this->m_Q_ch_nsamples; ++j) 
                           {
                                const float t_j{static_cast<float>(j*inv_sr)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
                                const float sin_val = ceph_sinf(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_j));
                                this->m_Q_ch_bitstream.m_data[j] = (sin_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
                                const float sin_val = std::sin(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_j));
                                this->m_Q_ch_bitstream.m_data[j] = (sin_val>=0.0f)?+1.0f:-1.0f;
#endif
                           }
                       }
                       else 
                       {
                       
                           for(i = 0ull; i != ROUND_TO_EIGHT(this->m_Q_ch_nsamples,8ull); i += 8ull) 
                           {
                               
                               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align32[i],_MM_HINT_T0);
                               const __m256 vt_i{_mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i]),vinv_sr)};
                               const __m256 vsin_val    = _mm256_sin_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i),vph0));
                               const __mmask8 vsin_ge_0 = _mm256_cmp_ps_mask(vsin_val,vzero,_CMP_GE_OQ);
                               _mm256_store_ps(&this->m_Q_ch_bitstream.m_data[i], _mm256_mask_blend_ps(vsin_ge_0,vnone,vpone));   
                           }
                       
                           for(j = i; j != this->m_Q_ch_nsamples; ++j) 
                           {
                                const float t_j{gms::math::LUT_loop_indices_2257_align32[j]*inv_sr};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
                                const float sin_val = ceph_sinf(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_j));
                                this->m_Q_ch_bitstream.m_data[j] = (sin_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
                                const float sin_val = std::sin(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_j));
                                this->m_Q_ch_bitstream.m_data[j] = (sin_val>=0.0f)?+1.0f:-1.0f;
#endif
                           }

                       }
                       return (0);
                  }

                   __ATTR_ALWAYS_INLINE__
                  inline 
                  std::int32_t 
                  generate_Q_channel_bitstream_avx_u4x(const bool do_constants_prefetch)                                                  
                  {
                       using namespace gms::math;
                       
                       const float inv_sr{1.0f/this->m_Q_ch_nsamples};
                       const __m256 vinv_sr{_mm256_set1_ps(inv_sr)};
                       const __m256 vw0{_mm256_set1_ps(this->m_Q_fc)};
                       const __m256 vph0{_mm256_set1_ps(this->m_Q_ph0)};
                       __m256 c0_7,c8_15,c16_23,c24_31;
                       __m256 v2pi,vzero,vpone,vnone;
                       if(!do_constants_prefetch)
                       {
                           c0_7   = _mm256_setr_ps(0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f);
                           c8_15  = _mm256_setr_ps(8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f);
                           c16_23 = _mm256_setr_ps(16.0f,17.0f,18.0f,19.0f,20.0f,21.0f,22.0f,23.0f);
                           c24_31 = _mm256_setr_ps(24.0f,25.0f,26.0f,27.0f,28.0f,29.0f,30.0f,31.0f);
                           v2pi   = _mm256_set1_ps(6.283185307179586476925286766559f);
                           vzero  = _mm256_setzero_ps();
                           vpone  = _mm256_set1_ps(+1.0f);
                           vnone  = _mm256_set1_ps(-1.0f);  
                       }
                       else 
                       {
                              __ATTR_ALIGN__(32) const float prefetched_constants[64] = {
                                                      0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f,          
                                                      8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f,
                                                      16.0f,17.0f,18.0f,19.0f,20.0f,21.0f,22.0f,23.0f,
                                                      24.0f,25.0f,26.0f,27.0f,28.0f,29.0f,30.0f,31.0f,
                                                      6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                      6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                      6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                      6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                      0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,
                                                      +1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,
                                                      -1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f};
                              _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[48],_MM_HINT_T0);
                              c0_7   = _mm256_load_ps(&prefetched_constants[0]);
                              c8_15  = _mm256_load_ps(&prefetched_constants[8]);
                              c16_23 = _mm256_load_ps(&prefetched_constants[16]);
                              c24_31 = _mm256_load_ps(&prefetched_constants[24]);
                              v2pi   = _mm256_load_ps(&prefetched_constants[32]);
                              vzero  = _mm256_load_ps(&prefetched_constants[40]);
                              vpone  = _mm256_load_ps(&prefetched_constants[48]);
                              vnone  = _mm256_load_ps(&prefetched_constants[56]);
                       }
                       std::size_t i;
                       float j;
                       constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
                       constexpr std::size_t LUT_loop_idx_threshold{2257ull};
                       if(this->m_Q_ch_nsamples>LUT_loop_idx_threshold)
                       {
                             for(i = 0ull,j = 0.0f; (i+31ull) < this->m_Q_ch_nsamples; i += 32ull,j += 32.0f) 
                             {
                                   const __m256 vt_i_0{_mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c0_7),vinv_sr)};
                                   const __m256 vsin_val_0    = _mm256_sin_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i_0),vph0));
                                   const __mmask8 vsin_ge_0_0 = _mm256_cmp_ps_mask(vsin_val_0,vzero,_CMP_GE_OQ);
                                   _mm256_store_ps(&this->m_Q_ch_bitstream.m_data[i+0ull], _mm256_mask_blend_ps(vsin_ge_0_0,vnone,vpone));
                                   const __m256 vt_i_1{_mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c8_15),vinv_sr)};
                                   const __m256 vsin_val_1    = _mm256_sin_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i_1),vph0));
                                   const __mmask8 vsin_ge_0_1 = _mm256_cmp_ps_mask(vsin_val_1,vzero,_CMP_GE_OQ);
                                   _mm256_store_ps(&this->m_Q_ch_bitstream.m_data[i+8ull], _mm256_mask_blend_ps(vsin_ge_0_1,vnone,vpone));
                                   const __m256 vt_i_2{_mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c16_23),vinv_sr)};
                                   const __m256 vsin_val_2    = _mm256_sin_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i_2),vph0));
                                   const __mmask8 vsin_ge_0_2 = _mm256_cmp_ps_mask(vsin_val_2,vzero,_CMP_GE_OQ);
                                   _mm256_store_ps(&this->m_Q_ch_bitstream.m_data[i+16ull], _mm256_mask_blend_ps(vsin_ge_0_2,vnone,vpone));
                                   const __m256 vt_i_3{_mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c24_31),vinv_sr)};
                                   const __m256 vsin_val_3    = _mm256_sin_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i_3),vph0));
                                   const __mmask8 vsin_ge_0_3 = _mm256_cmp_ps_mask(vsin_val_3,vzero,_CMP_GE_OQ);
                                   _mm256_store_ps(&this->m_Q_ch_bitstream.m_data[i+24ull], _mm256_mask_blend_ps(vsin_ge_0_3,vnone,vpone));
                              }

                              for(; (i+15ull) < this->m_Q_ch_nsamples; i += 16ull,j += 16.0f) 
                              {  
                                   const __m256 vt_i_0{_mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c0_7),vinv_sr)};
                                   const __m256 vsin_val_0    = _mm256_sin_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i_0),vph0));
                                   const __mmask8 vsin_ge_0_0 = _mm256_cmp_ps_mask(vsin_val_0,vzero,_CMP_GE_OQ);
                                   _mm256_store_ps(&this->m_Q_ch_bitstream.m_data[i+0ull], _mm256_mask_blend_ps(vsin_ge_0_0,vnone,vpone));
                                   const __m256 vt_i_1{_mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c8_15),vinv_sr)};
                                   const __m256 vsin_val_1    = _mm256_sin_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i_1),vph0));
                                   const __mmask8 vsin_ge_0_1 = _mm256_cmp_ps_mask(vsin_val_1,vzero,_CMP_GE_OQ);
                                   _mm256_store_ps(&this->m_Q_ch_bitstream.m_data[i+8ull], _mm256_mask_blend_ps(vsin_ge_0_1,vnone,vpone));
                              }

                              for(; (i+7ull) < this->m_Q_ch_nsamples; i += 8ull,j += 8.0f) 
                              {
                                   const __m256 vt_i_0{_mm256_mul_ps(_mm256_add_ps(_mm256_set1_ps(j),c0_7),vinv_sr)};
                                   const __m256 vsin_val_0    = _mm256_sin_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i_0),vph0));
                                   const __mmask8 vsin_ge_0_0 = _mm256_cmp_ps_mask(vsin_val_0,vzero,_CMP_GE_OQ);
                                   _mm256_store_ps(&this->m_Q_ch_bitstream.m_data[i+0ull], _mm256_mask_blend_ps(vsin_ge_0_0,vnone,vpone));
                              }
                                            
                              for(; (i+0ull) < this->m_Q_ch_nsamples; i += 1ull,j += 1.0f) 
                              {
                                      const float t_i{j*inv_sr};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
                                      const float sin_val = ceph_sinf(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_i));
                                      this->m_Q_ch_bitstream.m_data[i] = (sin_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
                                      const float sin_val = std::sin(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_i));
                                      this->m_Q_ch_bitstream.m_data[i] = (sin_val>=0.0f)?+1.0f:-1.0f;
#endif
                               }
                       }
                       else 
                       {
                       
                            for(i = 0ull; (i+31ull) < this->m_Q_ch_nsamples; i += 32ull) 
                            {
                                
                                 _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align32[i+0ull],_MM_HINT_T0);
                                 const __m256 vt_i_0{_mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i+0ull]),vinv_sr)};
                                 const __m256 vsin_val_0    = _mm256_sin_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i_0),vph0));
                                 const __mmask8 vsin_ge_0_0 = _mm256_cmp_ps_mask(vsin_val_0,vzero,_CMP_GE_OQ);
                                 _mm256_store_ps(&this->m_Q_ch_bitstream.m_data[i+0ull], _mm256_mask_blend_ps(vsin_ge_0_0,vnone,vpone));
                                
                                 const __m256 vt_i_1{_mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i+8ull]),vinv_sr)};
                                 const __m256 vsin_val_1    = _mm256_sin_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i_1),vph0));
                                 const __mmask8 vsin_ge_0_1 = _mm256_cmp_ps_mask(vsin_val_1,vzero,_CMP_GE_OQ);
                                 _mm256_store_ps(&this->m_Q_ch_bitstream.m_data[i+8ull], _mm256_mask_blend_ps(vsin_ge_0_1,vnone,vpone));
                                
                                 _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align32[i+16ull],_MM_HINT_T0);
                                 const __m256 vt_i_2{_mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i+16ull]),vinv_sr)};
                                 const __m256 vsin_val_2    = _mm256_sin_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i_2),vph0));
                                 const __mmask8 vsin_ge_0_2 = _mm256_cmp_ps_mask(vsin_val_2,vzero,_CMP_GE_OQ);
                                 _mm256_store_ps(&this->m_Q_ch_bitstream.m_data[i+16ull], _mm256_mask_blend_ps(vsin_ge_0_2,vnone,vpone));
                               
                                 const __m256 vt_i_3{_mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i+24ull]),vinv_sr)};
                                 const __m256 vsin_val_3    = _mm256_sin_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i_3),vph0));
                                 const __mmask8 vsin_ge_0_3 = _mm256_cmp_ps_mask(vsin_val_3,vzero,_CMP_GE_OQ);
                                 _mm256_store_ps(&this->m_Q_ch_bitstream.m_data[i+24ull], _mm256_mask_blend_ps(vsin_ge_0_3,vnone,vpone));
                           }

                           for(; (i+15ull) < this->m_Q_ch_nsamples; i += 16ull) 
                           {         
                                _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align32[i+0ull],_MM_HINT_T0);
                                 const __m256 vt_i_0{_mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i+0ull]),vinv_sr)};
                                 const __m256 vsin_val_0    = _mm256_sin_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i_0),vph0));
                                 const __mmask8 vsin_ge_0_0 = _mm256_cmp_ps_mask(vsin_val_0,vzero,_CMP_GE_OQ);
                                 _mm256_store_ps(&this->m_Q_ch_bitstream.m_data[i+0ull], _mm256_mask_blend_ps(vsin_ge_0_0,vnone,vpone));
                                
                                 const __m256 vt_i_1{_mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i+8ull]),vinv_sr)};
                                 const __m256 vsin_val_1    = _mm256_sin_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i_1),vph0));
                                 const __mmask8 vsin_ge_0_1 = _mm256_cmp_ps_mask(vsin_val_1,vzero,_CMP_GE_OQ);
                                 _mm256_store_ps(&this->m_Q_ch_bitstream.m_data[i+8ull], _mm256_mask_blend_ps(vsin_ge_0_1,vnone,vpone));
                          
                           }

                           for(; (i+7ull) < this->m_Q_ch_nsamples; i += 8ull) 
                           {
                                 const __m256 vt_i_0{_mm256_mul_ps(_mm256_load_ps(&gms::math::LUT_loop_indices_2257_align32[i+0ull]),vinv_sr)};
                                 const __m256 vsin_val_0    = _mm256_sin_ps(_mm256_fmadd_ps(v2pi,_mm256_mul_ps(vw0,vt_i_0),vph0));
                                 const __mmask8 vsin_ge_0_0 = _mm256_cmp_ps_mask(vsin_val_0,vzero,_CMP_GE_OQ);
                                 _mm256_store_ps(&this->m_Q_ch_bitstream.m_data[i+0ull], _mm256_mask_blend_ps(vsin_ge_0_0,vnone,vpone));
                           }
                      
                       
                           for(; (i+0ull) < this->m_Q_ch_nsamples; i += 1ull) 
                           {
                            
                                    const float t_i{gms::math::LUT_loop_indices_2257_align32[i]*inv_sr};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
                                    const float sin_val = ceph_sinf(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_i));
                                    this->m_Q_ch_bitstream.m_data[i] = (sin_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
                                    const float sin_val = std::sin(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_i));
                                    this->m_Q_ch_bitstream.m_data[i] = (sin_val>=0.0f)?+1.0f:-1.0f;
#endif
                           }
                       }
                       return (0);
                  }

                   __ATTR_ALWAYS_INLINE__
                  inline 
                  std::int32_t 
                  generate_I_channel_bitstream_avx512(const bool do_constants_prefetch)                                                   
                  {
                       using namespace gms::math;
                      
                       const float inv_sr{1.0f/this->m_I_ch_nsamples};
                       const __m512 vinv_sr{_mm512_set1_ps(inv_sr)};
                       const __m512 vw0{_mm512_set1_ps(this->m_I_fc)};
                       const __m512 vph0{_mm512_set1_ps(this->m_I_ph0)};
                       __m512 c0_15,v2pi;
                       __m512 vzero,vpone;
                       __m512 vnone;
                       if(!do_constants_prefetch)
                       {
                            c0_15 = _mm512_setr_ps(0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f,8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f);
                            v2pi  = _mm512_set1_ps(6.283185307179586476925286766559f);
                            vzero = _mm512_setzero_ps();
                            vpone = _mm512_set1_ps(+1.0f);
                            vnone = _mm512_set1_ps(-1.0f);
                       }
                       else 
                       {
                            __ATTR_ALIGN__(64) const float prefetched_constants[80] = {
                                                0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f,8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f,
                                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,
                                                +1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,
                                                -1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f}; 

                              _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[48],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[64],_MM_HINT_T0);
                              c0_15 = _mm512_load_ps(&prefetched_constants[0]); 
                              v2pi  = _mm512_load_ps(&prefetched_constants[16]);
                              vzero = _mm512_load_ps(&prefetched_constants[32]);
                              vpone = _mm512_load_ps(&prefetched_constants[48]);
                              vnone = _mm512_load_ps(&prefetched_constants[64]);
                       }
                       
                       std::size_t i,j;
                       float jj;
                       constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
                       constexpr std::size_t LUT_loop_idx_threshold{2257ull};
                       if(this->m_I_ch_nsamples>LUT_loop_idx_threshold)
                       {
                            for(i = 0ull,jj = 0.0f; i != ROUND_TO_SIXTEEN(this->m_I_ch_nsamples,16ull); i += 16ull,jj += 16.0f) 
                            {
                               const __m512 vt_i{_mm512_mul_ps(_mm512_add_ps(_mm512_set1_ps(jj),c0_15),vinv_sr)};
                               const __m512 vcos_val    = _mm512_cos_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i),vph0));
                               const __mmask16 vcos_ge_0 = _mm512_cmp_ps_mask(vcos_val,vzero,_CMP_GE_OQ);
                               _mm512_store_ps(&this->m_I_ch_bitstream.m_data[i], _mm512_mask_blend_ps(vcos_ge_0,vnone,vpone));   
                           }
                       
                           for(j = i; j != this->m_I_ch_nsamples; ++j) 
                           {
                                const float t_j{static_cast<float>(j*inv_sr)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
                                const float cos_val = ceph_cosf(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_j));
                                this->m_I_ch_bitstream.m_data[j] = (cos_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
                                const float cos_val = std::cos(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_j));
                                this->m_I_ch_bitstream.m_data[j] = (cos_val>=0.0f)?+1.0f:-1.0f;
#endif
                           }
                       }
                       else 
                       {
                       
                           for(i = 0ull; i != ROUND_TO_SIXTEEN(this->m_I_ch_nsamples,16ull); i += 16ull) 
                           {
                               
                               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i],_MM_HINT_T0);
                               const __m512 vt_i{_mm512_mul_ps(_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i]),vinv_sr)};
                               const __m512 vcos_val    = _mm512_cos_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i),vph0));
                               const __mmask16 vcos_ge_0 = _mm512_cmp_ps_mask(vcos_val,vzero,_CMP_GE_OQ);
                               _mm512_store_ps(&this->m_I_ch_bitstream.m_data[i], _mm512_mask_blend_ps(vcos_ge_0,vnone,vpone));   
                           }
                       
                           for(j = i; j != this->m_I_ch_nsamples; ++j) 
                           {
                                
                                const float t_j{gms::math::LUT_loop_indices_2257_align64[j]*inv_sr};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
                                const float cos_val = ceph_cosf(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_j));
                                this->m_I_ch_bitstream.m_data[j] = (cos_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
                                const float cos_val = std::cos(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_j));
                                this->m_I_ch_bitstream.m_data[j] = (cos_val>=0.0f)?+1.0f:-1.0f;
#endif
                           }

                       }
                       return (0);
                  }

                   __ATTR_ALWAYS_INLINE__
                  inline 
                  std::int32_t 
                  generate_I_channel_bitstream_avx512_u4x(const bool do_prefetch_constants)                                                 
                  {
                       using namespace gms::math;                                           
                       const float inv_sr{1.0f/this->m_I_ch_nsamples};
                       const __m512 vinv_sr{_mm512_set1_ps(inv_sr)};
                       const __m512 vw0{_mm512_set1_ps(this->m_I_fc)};
                       const __m512 vph0{_mm512_set1_ps(this->m_I_ph0)};
                       __m512 c0_15,c16_31;
                       __m512 c32_47,c48_63;
                       __m512 v2pi,vzero;
                       __m512 vpone,vnone;
                       if(!do_prefetch_constants)
                       {
                             c0_15 = _mm512_setr_ps(0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f,8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f);
                             c16_31= _mm512_setr_ps(16.0f,17.0f,18.0f,19.0f,20.0f,21.0f,22.0f,23.0f,24.0f,25.0f,26.0f,27.0f,28.0f,29.0f,30.0f,31.0f);
                             c32_47= _mm512_setr_ps(32.0f,33.0f,34.0f,35.0f,36.0f,37.0f,38.0f,39.0f,40.0f,41.0f,42.0f,43.0f,44.0f,45.0f,46.0f,47.0f);
                             c48_63= _mm512_setr_ps(48.0f,49.0f,50.0f,51.0f,52.0f,53.0f,54.0f,55.0f,56.0f,57.0f,58.0f,59.0f,60.0f,61.0f,62.0f,63.0f);
                             v2pi  = _mm512_set1_ps(6.283185307179586476925286766559f);
                             vzero = _mm512_setzero_ps();
                             vpone = _mm512_set1_ps(+1.0f);
                             vnone = _mm512_set1_ps(-1.0f);
                       }
                       else 
                       {
                             __ATTR_ALIGN__(64) const float prefetched_constants[128] = {
                                                            0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f,8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f,
                                                            16.0f,17.0f,18.0f,19.0f,20.0f,21.0f,22.0f,23.0f,24.0f,25.0f,26.0f,27.0f,28.0f,29.0f,30.0f,31.0f,
                                                            32.0f,33.0f,34.0f,35.0f,36.0f,37.0f,38.0f,39.0f,40.0f,41.0f,42.0f,43.0f,44.0f,45.0f,46.0f,47.0f,
                                                            48.0f,49.0f,50.0f,51.0f,52.0f,53.0f,54.0f,55.0f,56.0f,57.0f,58.0f,59.0f,60.0f,61.0f,62.0f,63.0f,
                                                            6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                            6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                            6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                            6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                            6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                            6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                            6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                            6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                            0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,
                                                            +1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,
                                                           -1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f}; 
                              
                              _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[48],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[80],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[112],_MM_HINT_T0);
                              c0_15 = _mm512_load_ps(&prefetched_constants[0]);
                              c16_31= _mm512_load_ps(&prefetched_constants[16]);
                              c32_47= _mm512_load_ps(&prefetched_constants[32]);
                              c48_63= _mm512_load_ps(&prefetched_constants[48]);
                              v2pi  = _mm512_load_ps(&prefetched_constants[64]);
                              vzero = _mm512_load_ps(&prefetched_constants[80]);
                              vpone = _mm512_load_ps(&prefetched_constants[96]);
                              vnone = _mm512_load_ps(&prefetched_constants[112]);
                       }
                       
                       std::size_t i;
                       float j;
                       constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
                       constexpr std::size_t LUT_loop_idx_threshold{2257ull};
                       if(this->m_I_ch_nsamples>LUT_loop_idx_threshold)
                       {
                             for(i = 0ull,j = 0.0f; (i+63ull) < this->m_I_ch_nsamples; i += 64ull,j += 64.0f) 
                             {
                                   const __m512 vt_i_0{_mm512_mul_ps(_mm512_add_ps(_mm512_set1_ps(j),c0_15),vinv_sr)};
                                   const __m512 vcos_val_0    = _mm512_cos_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i_0),vph0));
                                   const __mmask16 vcos_ge_0_0 = _mm512_cmp_ps_mask(vcos_val_0,vzero,_CMP_GE_OQ);
                                   _mm512_store_ps(&this->m_I_ch_bitstream.m_data[i+0ull], _mm512_mask_blend_ps(vcos_ge_0_0,vnone,vpone));
                                   const __m512 vt_i_1{_mm512_mul_ps(_mm512_add_ps(_mm512_set1_ps(j),c16_31),vinv_sr)};
                                   const __m512 vcos_val_1    = _mm512_cos_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i_1),vph0));
                                   const __mmask16 vcos_ge_0_1 = _mm512_cmp_ps_mask(vcos_val_1,vzero,_CMP_GE_OQ);
                                   _mm512_store_ps(&this->m_I_ch_bitstream.m_data[i+16ull], _mm512_mask_blend_ps(vcos_ge_0_1,vnone,vpone));
                                   const __m512 vt_i_2{_mm512_mul_ps(_mm512_add_ps(_mm512_set1_ps(j),c32_47),vinv_sr)};
                                   const __m512 vcos_val_2    = _mm512_cos_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i_2),vph0));
                                   const __mmask16 vcos_ge_0_2 = _mm512_cmp_ps_mask(vcos_val_2,vzero,_CMP_GE_OQ);
                                   _mm512_store_ps(&this->m_I_ch_bitstream.m_data[i+32ull], _mm512_mask_blend_ps(vcos_ge_0_2,vnone,vpone));
                                   const __m512 vt_i_3{_mm512_mul_ps(_mm512_add_ps(_mm512_set1_ps(j),c48_63),vinv_sr)};
                                   const __m512 vcos_val_3    = _mm512_cos_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i_3),vph0));
                                   const __mmask16 vcos_ge_0_3 = _mm512_cmp_ps_mask(vcos_val_3,vzero,_CMP_GE_OQ);
                                   _mm512_store_ps(&this->m_I_ch_bitstream.m_data[i+48ull], _mm512_mask_blend_ps(vcos_ge_0_3,vnone,vpone));
                              }

                              for(; (i+31ull) < this->m_I_ch_nsamples; i += 32ull,j += 32.0f) 
                              {  
                                   const __m512 vt_i_0{_mm512_mul_ps(_mm512_add_ps(_mm512_set1_ps(j),c0_15),vinv_sr)};
                                   const __m512 vcos_val_0    = _mm512_cos_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i_0),vph0));
                                   const __mmask16 vcos_ge_0_0 = _mm512_cmp_ps_mask(vcos_val_0,vzero,_CMP_GE_OQ);
                                   _mm512_store_ps(&this->m_I_ch_bitstream.m_data[i+0ull], _mm512_mask_blend_ps(vcos_ge_0_0,vnone,vpone));
                                    const __m512 vt_i_1{_mm512_mul_ps(_mm512_add_ps(_mm512_set1_ps(j),c16_31),vinv_sr)};
                                   const __m512 vcos_val_1    = _mm512_cos_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i_1),vph0));
                                   const __mmask16 vcos_ge_0_1 = _mm512_cmp_ps_mask(vcos_val_1,vzero,_CMP_GE_OQ);
                                   _mm512_store_ps(&this->m_I_ch_bitstream.m_data[i+16ull], _mm512_mask_blend_ps(vcos_ge_0_1,vnone,vpone));
                              }

                              for(; (i+15ull) < this->m_I_ch_nsamples; i += 16ull,j += 16.0f) 
                              {
                                   const __m512 vt_i_0{_mm512_mul_ps(_mm512_add_ps(_mm512_set1_ps(j),c0_15),vinv_sr)};
                                   const __m512 vcos_val_0    = _mm512_cos_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i_0),vph0));
                                   const __mmask16 vcos_ge_0_0 = _mm512_cmp_ps_mask(vcos_val_0,vzero,_CMP_GE_OQ);
                                   _mm512_store_ps(&this->m_I_ch_bitstream.m_data[i+0ull], _mm512_mask_blend_ps(vcos_ge_0_0,vnone,vpone));
                              }
                                            
                              for(; (i+0ull) < this->m_I_ch_nsamples; i += 1ull,j += 1.0f) 
                              {
                                      const float t_i{j*inv_sr};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
                                      const float cos_val = ceph_cosf(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_i));
                                      this->m_I_ch_bitstream.m_data[i] = (cos_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
                                      const float cos_val = std::cos(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_i));
                                      this->m_I_ch_bitstream.m_data[i] = (cos_val>=0.0f)?+1.0f:-1.0f;
#endif
                               }
                       }
                       else 
                       {
                       
                            for(i = 0ull; (i+63ull) < this->m_I_ch_nsamples; i += 64ull) 
                            {
                                 
                                 _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+0ull],_MM_HINT_T0);
                                 const __m512 vt_i_0{_mm512_mul_ps(_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i+0ull]),vinv_sr)};
                                 const __m512 vcos_val_0    = _mm512_cos_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i_0),vph0));
                                 const __mmask16 vcos_ge_0_0 = _mm512_cmp_ps_mask(vcos_val_0,vzero,_CMP_GE_OQ);
                                 _mm512_store_ps(&this->m_I_ch_bitstream.m_data[i+0ull], _mm512_mask_blend_ps(vcos_ge_0_0,vnone,vpone));
                                
                                 _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+16ull],_MM_HINT_T0);
                                 const __m512 vt_i_1{_mm512_mul_ps(_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i+16ull]),vinv_sr)};
                                 const __m512 vcos_val_1    = _mm512_cos_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i_1),vph0));
                                 const __mmask16 vcos_ge_0_1 = _mm512_cmp_ps_mask(vcos_val_1,vzero,_CMP_GE_OQ);
                                 _mm512_store_ps(&this->m_I_ch_bitstream.m_data[i+16ull], _mm512_mask_blend_ps(vcos_ge_0_1,vnone,vpone));
                                 
                                 _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+32ull],_MM_HINT_T0);
                                 const __m512 vt_i_2{_mm512_mul_ps(_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i+32ull]),vinv_sr)};
                                 const __m512 vcos_val_2    = _mm512_cos_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i_2),vph0));
                                 const __mmask16 vcos_ge_0_2 = _mm512_cmp_ps_mask(vcos_val_2,vzero,_CMP_GE_OQ);
                                 _mm512_store_ps(&this->m_I_ch_bitstream.m_data[i+32ull], _mm512_mask_blend_ps(vcos_ge_0_2,vnone,vpone));
                                 
                                 _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+48ull],_MM_HINT_T0);
                                 const __m512 vt_i_3{_mm512_mul_ps(_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i+48ull]),vinv_sr)};
                                 const __m512 vcos_val_3    = _mm512_cos_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i_3),vph0));
                                 const __mmask16 vcos_ge_0_3 = _mm512_cmp_ps_mask(vcos_val_3,vzero,_CMP_GE_OQ);
                                 _mm512_store_ps(&this->m_I_ch_bitstream.m_data[i+48ull], _mm512_mask_blend_ps(vcos_ge_0_3,vnone,vpone));
                           }

                           for(; (i+31ull) < this->m_I_ch_nsamples; i += 32ull) 
                           {         
                                _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+0ull],_MM_HINT_T0);
                                 const __m512 vt_i_0{_mm512_mul_ps(_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i+0ull]),vinv_sr)};
                                 const __m512 vcos_val_0    = _mm512_cos_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i_0),vph0));
                                 const __mmask16 vcos_ge_0_0 = _mm512_cmp_ps_mask(vcos_val_0,vzero,_CMP_GE_OQ);
                                 _mm512_store_ps(&this->m_I_ch_bitstream.m_data[i+0ull], _mm512_mask_blend_ps(vcos_ge_0_0,vnone,vpone));
                                
                                 _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+16ull],_MM_HINT_T0);
                                 const __m512 vt_i_1{_mm512_mul_ps(_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i+16ull]),vinv_sr)};
                                 const __m512 vcos_val_1    = _mm512_cos_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i_1),vph0));
                                 const __mmask16 vcos_ge_0_1 = _mm512_cmp_ps_mask(vcos_val_1,vzero,_CMP_GE_OQ);
                                 _mm512_store_ps(&this->m_I_ch_bitstream.m_data[i+16ull], _mm512_mask_blend_ps(vcos_ge_0_1,vnone,vpone));
                          
                           }

                           for(; (i+15ull) < this->m_I_ch_nsamples; i += 16ull) 
                           {
                                _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+0ull],_MM_HINT_T0);
                                 const __m512 vt_i_0{_mm512_mul_ps(_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i+0ull]),vinv_sr)};
                                 const __m512 vcos_val_0    = _mm512_cos_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i_0),vph0));
                                 const __mmask16 vcos_ge_0_0 = _mm512_cmp_ps_mask(vcos_val_0,vzero,_CMP_GE_OQ);
                                 _mm512_store_ps(&this->m_I_ch_bitstream.m_data[i+0ull], _mm512_mask_blend_ps(vcos_ge_0_0,vnone,vpone));
                           }
                      
                       
                           for(; (i+0ull) < this->m_I_ch_nsamples; i += 1ull) 
                           {
                            
                                    const float t_i{gms::math::LUT_loop_indices_2257_align64[i]*inv_sr};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
                                    const float cos_val = ceph_cosf(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_i));
                                    this->m_I_ch_bitstream.m_data[i] = (cos_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
                                    const float cos_val = std::cos(this->m_I_ph0+(C6283185307179586476925286766559*this->m_I_fc*t_i));
                                    this->m_I_ch_bitstream.m_data[i] = (cos_val>=0.0f)?+1.0f:-1.0f;
#endif
                           }
                       }
                       return (0);
                  }


                  __ATTR_ALWAYS_INLINE__
                  inline 
                  std::int32_t 
                  generate_Q_channel_bitstream_avx512(const bool do_constants_prefetch)                                                   
                  {
                       using namespace gms::math;
                      
                       const float inv_sr{1.0f/this->m_Q_ch_nsamples};
                       const __m512 vinv_sr{_mm512_set1_ps(inv_sr)};
                       const __m512 vw0{_mm512_set1_ps(this->m_Q_fc)};
                       const __m512 vph0{_mm512_set1_ps(this->m_Q_ph0)};
                       __m512 c0_15,v2pi;
                       __m512 vzero,vpone;
                       __m512 vnone;
                       if(!do_constants_prefetch)
                       {
                            c0_15 = _mm512_setr_ps(0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f,8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f);
                            v2pi  = _mm512_set1_ps(6.283185307179586476925286766559f);
                            vzero = _mm512_setzero_ps();
                            vpone = _mm512_set1_ps(+1.0f);
                            vnone = _mm512_set1_ps(-1.0f);
                       }
                       else 
                       {
                            __ATTR_ALIGN__(64) const float prefetched_constants[80] = {
                                                0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f,8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f,
                                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,
                                                +1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,
                                                -1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f}; 

                              _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[48],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[64],_MM_HINT_T0);
                              c0_15 = _mm512_load_ps(&prefetched_constants[0]); 
                              v2pi  = _mm512_load_ps(&prefetched_constants[16]);
                              vzero = _mm512_load_ps(&prefetched_constants[32]);
                              vpone = _mm512_load_ps(&prefetched_constants[48]);
                              vnone = _mm512_load_ps(&prefetched_constants[64]);
                       }
                       
                       std::size_t i,j;
                       float jj;
                       constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
                       constexpr std::size_t LUT_loop_idx_threshold{2257ull};
                       if(this->m_Q_ch_nsamples>LUT_loop_idx_threshold)
                       {
                            for(i = 0ull,jj = 0.0f; i != ROUND_TO_SIXTEEN(this->m_Q_ch_nsamples,16ull); i += 16ull,jj += 16.0f) 
                            {
                               const __m512 vt_i{_mm512_mul_ps(_mm512_add_ps(_mm512_set1_ps(jj),c0_15),vinv_sr)};
                               const __m512 vsin_val    = _mm512_sin_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i),vph0));
                               const __mmask16 vsin_ge_0 = _mm512_cmp_ps_mask(vsin_val,vzero,_CMP_GE_OQ);
                               _mm512_store_ps(&this->m_Q_ch_bitstream.m_data[i], _mm512_mask_blend_ps(vsin_ge_0,vnone,vpone));   
                           }
                       
                           for(j = i; j != this->m_Q_ch_nsamples; ++j) 
                           {
                                const float t_j{static_cast<float>(j*inv_sr)};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
                                const float sin_val = ceph_sinf(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_j));
                                this->m_Q_ch_bitstream.m_data[j] = (sin_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
                                const float sin_val = std::sin(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_j));
                                this->m_Q_ch_bitstream.m_data[j] = (sin_val>=0.0f)?+1.0f:-1.0f;
#endif
                           }
                       }
                       else 
                       {
                       
                           for(i = 0ull; i != ROUND_TO_SIXTEEN(this->m_Q_ch_nsamples,16ull); i += 16ull) 
                           {
                               
                               _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i],_MM_HINT_T0);
                               const __m512 vt_i{_mm512_mul_ps(_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i]),vinv_sr)};
                               const __m512 vsin_val    = _mm512_sin_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i),vph0));
                               const __mmask16 vsin_ge_0 = _mm512_cmp_ps_mask(vsin_val,vzero,_CMP_GE_OQ);
                               _mm512_store_ps(&this->m_Q_ch_bitstream.m_data[i], _mm512_mask_blend_ps(vsin_ge_0,vnone,vpone));   
                           }
                       
                           for(j = i; j != this->m_Q_ch_nsamples; ++j) 
                           {
                                
                                const float t_j{gms::math::LUT_loop_indices_2257_align64[j]*inv_sr};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
                                const float sin_val = ceph_sinf(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_j));
                                this->m_Q_ch_bitstream.m_data[j] = (sin_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
                                const float sin_val = std::sin(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_j));
                                this->m_Q_ch_bitstream.m_data[j] = (sin_val>=0.0f)?+1.0f:-1.0f;
#endif
                           }

                       }
                       return (0);
                  }


                  __ATTR_ALWAYS_INLINE__
                  inline 
                  std::int32_t 
                  generate_Q_channel_bitstream_avx512_u4x(const bool do_prefetch_constants)                                                 
                  {
                       using namespace gms::math;                                           
                       const float inv_sr{1.0f/this->m_Q_ch_nsamples};
                       const __m512 vinv_sr{_mm512_set1_ps(inv_sr)};
                       const __m512 vw0{_mm512_set1_ps(this->m_Q_fc)};
                       const __m512 vph0{_mm512_set1_ps(this->m_Q_ph0)};
                       __m512 c0_15,c16_31;
                       __m512 c32_47,c48_63;
                       __m512 v2pi,vzero;
                       __m512 vpone,vnone;
                       if(!do_prefetch_constants)
                       {
                             c0_15 = _mm512_setr_ps(0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f,8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f);
                             c16_31= _mm512_setr_ps(16.0f,17.0f,18.0f,19.0f,20.0f,21.0f,22.0f,23.0f,24.0f,25.0f,26.0f,27.0f,28.0f,29.0f,30.0f,31.0f);
                             c32_47= _mm512_setr_ps(32.0f,33.0f,34.0f,35.0f,36.0f,37.0f,38.0f,39.0f,40.0f,41.0f,42.0f,43.0f,44.0f,45.0f,46.0f,47.0f);
                             c48_63= _mm512_setr_ps(48.0f,49.0f,50.0f,51.0f,52.0f,53.0f,54.0f,55.0f,56.0f,57.0f,58.0f,59.0f,60.0f,61.0f,62.0f,63.0f);
                             v2pi  = _mm512_set1_ps(6.283185307179586476925286766559f);
                             vzero = _mm512_setzero_ps();
                             vpone = _mm512_set1_ps(+1.0f);
                             vnone = _mm512_set1_ps(-1.0f);
                       }
                       else 
                       {
                             __ATTR_ALIGN__(64) const float prefetched_constants[128] = {
                                                            0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f,8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f,
                                                            16.0f,17.0f,18.0f,19.0f,20.0f,21.0f,22.0f,23.0f,24.0f,25.0f,26.0f,27.0f,28.0f,29.0f,30.0f,31.0f,
                                                            32.0f,33.0f,34.0f,35.0f,36.0f,37.0f,38.0f,39.0f,40.0f,41.0f,42.0f,43.0f,44.0f,45.0f,46.0f,47.0f,
                                                            48.0f,49.0f,50.0f,51.0f,52.0f,53.0f,54.0f,55.0f,56.0f,57.0f,58.0f,59.0f,60.0f,61.0f,62.0f,63.0f,
                                                            6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                            6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                            6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                            6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                            6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                            6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                            6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                            6.283185307179586476925286766559f,6.283185307179586476925286766559f,
                                                            0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,
                                                            +1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,+1.0f,
                                                           -1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f}; 
                              
                              _mm_prefetch((const char*)&prefetched_constants[0],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[16],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[48],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[80],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&prefetched_constants[112],_MM_HINT_T0);
                              c0_15 = _mm512_load_ps(&prefetched_constants[0]);
                              c16_31= _mm512_load_ps(&prefetched_constants[16]);
                              c32_47= _mm512_load_ps(&prefetched_constants[32]);
                              c48_63= _mm512_load_ps(&prefetched_constants[48]);
                              v2pi  = _mm512_load_ps(&prefetched_constants[64]);
                              vzero = _mm512_load_ps(&prefetched_constants[80]);
                              vpone = _mm512_load_ps(&prefetched_constants[96]);
                              vnone = _mm512_load_ps(&prefetched_constants[112]);
                       }
                       
                       std::size_t i;
                       float j;
                       constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
                       constexpr std::size_t LUT_loop_idx_threshold{2257ull};
                       if(this->m_Q_ch_nsamples>LUT_loop_idx_threshold)
                       {
                             for(i = 0ull,j = 0.0f; (i+63ull) < this->m_Q_ch_nsamples; i += 64ull,j += 64.0f) 
                             {
                                   const __m512 vt_i_0{_mm512_mul_ps(_mm512_add_ps(_mm512_set1_ps(j),c0_15),vinv_sr)};
                                   const __m512 vsin_val_0    = _mm512_sin_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i_0),vph0));
                                   const __mmask16 vsin_ge_0_0 = _mm512_cmp_ps_mask(vsin_val_0,vzero,_CMP_GE_OQ);
                                   _mm512_store_ps(&this->m_Q_ch_bitstream.m_data[i+0ull], _mm512_mask_blend_ps(vsin_ge_0_0,vnone,vpone));
                                   const __m512 vt_i_1{_mm512_mul_ps(_mm512_add_ps(_mm512_set1_ps(j),c16_31),vinv_sr)};
                                   const __m512 vsin_val_1    = _mm512_sin_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i_1),vph0));
                                   const __mmask16 vsin_ge_0_1 = _mm512_cmp_ps_mask(vsin_val_1,vzero,_CMP_GE_OQ);
                                   _mm512_store_ps(&this->m_Q_ch_bitstream.m_data[i+16ull], _mm512_mask_blend_ps(vsin_ge_0_1,vnone,vpone));
                                   const __m512 vt_i_2{_mm512_mul_ps(_mm512_add_ps(_mm512_set1_ps(j),c32_47),vinv_sr)};
                                   const __m512 vsin_val_2    = _mm512_sin_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i_2),vph0));
                                   const __mmask16 vsin_ge_0_2 = _mm512_cmp_ps_mask(vsin_val_2,vzero,_CMP_GE_OQ);
                                   _mm512_store_ps(&this->m_Q_ch_bitstream.m_data[i+32ull], _mm512_mask_blend_ps(vsin_ge_0_2,vnone,vpone));
                                   const __m512 vt_i_3{_mm512_mul_ps(_mm512_add_ps(_mm512_set1_ps(j),c48_63),vinv_sr)};
                                   const __m512 vsin_val_3    = _mm512_sin_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i_3),vph0));
                                   const __mmask16 vsin_ge_0_3 = _mm512_cmp_ps_mask(vsin_val_3,vzero,_CMP_GE_OQ);
                                   _mm512_store_ps(&this->m_Q_ch_bitstream.m_data[i+48ull], _mm512_mask_blend_ps(vsin_ge_0_3,vnone,vpone));
                              }

                              for(; (i+31ull) < this->m_Q_ch_nsamples; i += 32ull,j += 32.0f) 
                              {  
                                   const __m512 vt_i_0{_mm512_mul_ps(_mm512_add_ps(_mm512_set1_ps(j),c0_15),vinv_sr)};
                                   const __m512 vsin_val_0    = _mm512_sin_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i_0),vph0));
                                   const __mmask16 vsin_ge_0_0 = _mm512_cmp_ps_mask(vsin_val_0,vzero,_CMP_GE_OQ);
                                   _mm512_store_ps(&this->m_Q_ch_bitstream.m_data[i+0ull], _mm512_mask_blend_ps(vsin_ge_0_0,vnone,vpone));
                                    const __m512 vt_i_1{_mm512_mul_ps(_mm512_add_ps(_mm512_set1_ps(j),c16_31),vinv_sr)};
                                   const __m512 vsin_val_1    = _mm512_sin_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i_1),vph0));
                                   const __mmask16 vsin_ge_0_1 = _mm512_cmp_ps_mask(vsin_val_1,vzero,_CMP_GE_OQ);
                                   _mm512_store_ps(&this->m_Q_ch_bitstream.m_data[i+16ull], _mm512_mask_blend_ps(vsin_ge_0_1,vnone,vpone));
                              }

                              for(; (i+15ull) < this->m_Q_ch_nsamples; i += 16ull,j += 16.0f) 
                              {
                                   const __m512 vt_i_0{_mm512_mul_ps(_mm512_add_ps(_mm512_set1_ps(j),c0_15),vinv_sr)};
                                   const __m512 vsin_val_0    = _mm512_sin_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i_0),vph0));
                                   const __mmask16 vsin_ge_0_0 = _mm512_cmp_ps_mask(vsin_val_0,vzero,_CMP_GE_OQ);
                                   _mm512_store_ps(&this->m_Q_ch_bitstream.m_data[i+0ull], _mm512_mask_blend_ps(vsin_ge_0_0,vnone,vpone));
                              }
                                            
                              for(; (i+0ull) < this->m_Q_ch_nsamples; i += 1ull,j += 1.0f) 
                              {
                                      const float t_i{j*inv_sr};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
                                      const float sin_val = ceph_sinf(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_i));
                                      this->m_Q_ch_bitstream.m_data[i] = (sin_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
                                      const float sin_val = std::sin(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_i));
                                      this->m_Q_ch_bitstream.m_data[i] = (sin_val>=0.0f)?+1.0f:-1.0f;
#endif
                               }
                       }
                       else 
                       {
                       
                            for(i = 0ull; (i+63ull) < this->m_Q_ch_nsamples; i += 64ull) 
                            {
                                 
                                 _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+0ull],_MM_HINT_T0);
                                 const __m512 vt_i_0{_mm512_mul_ps(_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i+0ull]),vinv_sr)};
                                 const __m512 vsin_val_0    = _mm512_sin_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i_0),vph0));
                                 const __mmask16 vsin_ge_0_0 = _mm512_cmp_ps_mask(vsin_val_0,vzero,_CMP_GE_OQ);
                                 _mm512_store_ps(&this->m_Q_ch_bitstream.m_data[i+0ull], _mm512_mask_blend_ps(vsin_ge_0_0,vnone,vpone));
                                
                                 _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+16ull],_MM_HINT_T0);
                                 const __m512 vt_i_1{_mm512_mul_ps(_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i+16ull]),vinv_sr)};
                                 const __m512 vsin_val_1    = _mm512_sin_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i_1),vph0));
                                 const __mmask16 vsin_ge_0_1 = _mm512_cmp_ps_mask(vsin_val_1,vzero,_CMP_GE_OQ);
                                 _mm512_store_ps(&this->m_Q_ch_bitstream.m_data[i+16ull], _mm512_mask_blend_ps(vsin_ge_0_1,vnone,vpone));
                                 
                                 _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+32ull],_MM_HINT_T0);
                                 const __m512 vt_i_2{_mm512_mul_ps(_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i+32ull]),vinv_sr)};
                                 const __m512 vsin_val_2    = _mm512_sin_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i_2),vph0));
                                 const __mmask16 vsin_ge_0_2 = _mm512_cmp_ps_mask(vsin_val_2,vzero,_CMP_GE_OQ);
                                 _mm512_store_ps(&this->m_Q_ch_bitstream.m_data[i+32ull], _mm512_mask_blend_ps(vsin_ge_0_2,vnone,vpone));
                                 
                                 _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+48ull],_MM_HINT_T0);
                                 const __m512 vt_i_3{_mm512_mul_ps(_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i+48ull]),vinv_sr)};
                                 const __m512 vsin_val_3    = _mm512_sin_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i_3),vph0));
                                 const __mmask16 vsin_ge_0_3 = _mm512_cmp_ps_mask(vsin_val_3,vzero,_CMP_GE_OQ);
                                 _mm512_store_ps(&this->m_Q_ch_bitstream.m_data[i+48ull], _mm512_mask_blend_ps(vsin_ge_0_3,vnone,vpone));
                           }

                           for(; (i+31ull) < this->m_Q_ch_nsamples; i += 32ull) 
                           {         
                                _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+0ull],_MM_HINT_T0);
                                 const __m512 vt_i_0{_mm512_mul_ps(_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i+0ull]),vinv_sr)};
                                 const __m512 vsin_val_0    = _mm512_sin_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i_0),vph0));
                                 const __mmask16 vsin_ge_0_0 = _mm512_cmp_ps_mask(vsin_val_0,vzero,_CMP_GE_OQ);
                                 _mm512_store_ps(&this->m_Q_ch_bitstream.m_data[i+0ull], _mm512_mask_blend_ps(vsin_ge_0_0,vnone,vpone));
                                
                                 _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+16ull],_MM_HINT_T0);
                                 const __m512 vt_i_1{_mm512_mul_ps(_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i+16ull]),vinv_sr)};
                                 const __m512 vsin_val_1    = _mm512_sin_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i_1),vph0));
                                 const __mmask16 vsin_ge_0_1 = _mm512_cmp_ps_mask(vsin_val_1,vzero,_CMP_GE_OQ);
                                 _mm512_store_ps(&this->m_Q_ch_bitstream.m_data[i+16ull], _mm512_mask_blend_ps(vsin_ge_0_1,vnone,vpone));
                          
                           }

                           for(; (i+15ull) < this->m_Q_ch_nsamples; i += 16ull) 
                           {
                                _mm_prefetch((const char*)&gms::math::LUT_loop_indices_2257_align64[i+0ull],_MM_HINT_T0);
                                 const __m512 vt_i_0{_mm512_mul_ps(_mm512_load_ps(&gms::math::LUT_loop_indices_2257_align64[i+0ull]),vinv_sr)};
                                 const __m512 vsin_val_0    = _mm512_sin_ps(_mm512_fmadd_ps(v2pi,_mm512_mul_ps(vw0,vt_i_0),vph0));
                                 const __mmask16 vsin_ge_0_0 = _mm512_cmp_ps_mask(vsin_val_0,vzero,_CMP_GE_OQ);
                                 _mm512_store_ps(&this->m_Q_ch_bitstream.m_data[i+0ull], _mm512_mask_blend_ps(vsin_ge_0_0,vnone,vpone));
                           }
                      
                       
                           for(; (i+0ull) < this->m_Q_ch_nsamples; i += 1ull) 
                           {
                            
                                    const float t_i{gms::math::LUT_loop_indices_2257_align64[i]*inv_sr};
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1 
                                    const float sin_val = ceph_sinf(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_i));
                                    this->m_Q_ch_bitstream.m_data[i] = (sin_val>=0.0f)?+1.0f:-1.0f;
                                    
#else 
                                    const float sin_val = std::sin(this->m_Q_ph0+(C6283185307179586476925286766559*this->m_Q_fc*t_i));
                                    this->m_Q_ch_bitstream.m_data[i] = (sin_val>=0.0f)?+1.0f:-1.0f;
#endif
                           }
                       }
                       return (0);
                  }

                  __ATTR_ALWAYS_INLINE__
                  inline 
                  float pulse_shaping_function_value(const float t) const 
                  {
                        
                        using namespace gms::math; 
                        constexpr float C314159265358979323846264338328{3.14159265358979323846264338328f};
                        constexpr float C6283185307179586476925286766559{6.283185307179586476925286766559f};
                        const     float inv2T{1.0f/(2.0f*this->m_T)};
                        const     float invT{1.0/this->m_T};
                        float           psf_value{0.0f};
                        const float     cos_arg = C314159265358979323846264338328*t*inv2T;
#if (SINUSOIDAL_FSK_USE_CEPHES) == 1
                        const float sin_arg = 0.25f*ceph_sinf(C6283185307179586476925286766559*t*invT);
                        psf_value           = ceph_cosf(cos_arg-sin_arg);
#else 
                        const float sin_arg = 0.25f*std:sin(C6283185307179586476925286766559*t*invT);
                        psf_value           = ceph_cosf(cos_arg-sin_arg);
#endif 
                        return (psf_value);
                  }


          __ATTR_HOT__
          __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
          __ATTR_OPTIMIZE_03__
#endif                   
                  std::int32_t 
                  generate_pulse_shaping_function_u8x();

          __ATTR_HOT__
          __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
          __ATTR_OPTIMIZE_03__
#endif     
                  std::int32_t 
                  generate_pulse_shaping_function_sse_u4x(const bool);

          __ATTR_HOT__
          __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
          __ATTR_OPTIMIZE_03__
#endif     
                  std::int32_t 
                  generate_pulse_shaping_function_avx_u4x(const bool);   


          __ATTR_HOT__
          __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
          __ATTR_OPTIMIZE_03__
#endif     
                  std::int32_t 
                  generate_pulse_shaping_function_avx512_u4x(const bool); 


                  template<I_channel_bitstream_optimized_path I_ch_path,
                           Q_channel_bitstream_optimized_path Q_ch_path,
                           pulse_shaping_function_optimized_path psf_path>
                  std::int32_t  generate_fsk_signal_u8x(const iq_rectw_bitstream_vsequence_t &,
                                                        const std::int32_t,
                                                        const std::int32_t,
                                                        const bool,
                                                        const bool,
                                                        const bool,
                                                        const bool);
                                          
          
                 template<I_channel_bitstream_optimized_path I_ch_path,
                           Q_channel_bitstream_optimized_path Q_ch_path,
                           pulse_shaping_function_optimized_path psf_path> 
                std::int32_t generate_fsk_signal_sse_u4x(const iq_rectw_bitstream_vsequence_t &,
                                                         const std::int32_t,
                                                         const std::int32_t,
                                                         const bool,
                                                         const bool,
                                                         const bool,
                                                         const bool,
                                                         const bool); 

           
                 template<I_channel_bitstream_optimized_path I_ch_path,
                           Q_channel_bitstream_optimized_path Q_ch_path,
                           pulse_shaping_function_optimized_path psf_path>   
                std::int32_t generate_fsk_signal_avx_u4x(const iq_rectw_bitstream_vsequence_t &,
                                                         const std::int32_t,
                                                         const std::int32_t,
                                                         const bool,
                                                         const bool,
                                                         const bool,
                                                         const bool,
                                                         const bool);   

                                     
                 template<I_channel_bitstream_optimized_path I_ch_path,
                          Q_channel_bitstream_optimized_path Q_ch_path,
                          pulse_shaping_function_optimized_path psf_path>
                 std::int32_t generate_fsk_signal_avx512_u4x(const iq_rectw_bitstream_vsequence_t &,
                                                             const std::int32_t,
                                                             const std::int32_t,
                                                             const bool,
                                                             const bool,
                                                             const bool,
                                                             const bool,
                                                             const bool);

                 


            };

           auto 
           operator<<(std::ostream &,
                      const sinusoidal_fsk_t &)->std::ostream &;

} //radiolocation

}// gms
























#endif /*__GMS_SINUSOIDAL_FSK_H__*/