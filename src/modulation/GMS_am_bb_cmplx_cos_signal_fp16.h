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

#ifndef __GMS_AM_BB_CMPLX_COS_SIGNAL_FP16_H__
#define __GMS_AM_BB_CMPLX_COS_SIGNAL_FP16_H__ 201020251218

namespace file_info 
{

     static const unsigned int GMS_AM_BB_CMPLX_COS_SIGNAL_FP16_MAJOR = 1;
     static const unsigned int GMS_AM_BB_CMPLX_COS_SIGNAL_FP16_MINOR = 1;
     static const unsigned int GMS_AM_BB_CMPLX_COS_SIGNAL_FP16_MICRO = 0;
     static const unsigned int GMS_AM_BB_CMPLX_COS_SIGNAL_FP16_FULLVER =
       1000U*GMS_AM_BB_CMPLX_COS_SIGNAL_FP16_MAJOR+100U*GMS_AM_BB_CMPLX_COS_SIGNAL_FP16_MINOR+
       10U*GMS_AM_BB_CMPLX_COS_SIGNAL_FP16_MICRO;
     static const char GMS_AM_BB_CMPLX_COS_SIGNAL_FP16_CREATION_DATE[] = "20-10-2025 12:18 +00200 (MON 07 OCT 2025 GMT+2)";
     static const char GMS_AM_BB_CMPLX_COS_SIGNAL_FP16_BUILD_DATE[]    = __DATE__; 
     static const char GMS_AM_BB_CMPLX_COS_SIGNAL_FP16_BUILD_TIME[]    = __TIME__;
     static const char GMS_AM_BB_CMPLX_COS_SIGNAL_FP16_SYNOPSIS[]      = "AM baseband complex cosine wave signal half-float stored.";

}

#include <cstdint>
#include <string>
#include <iostream>
#include "GMS_config.h"
#include "GMS_dyn_array.h"

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
// To be added.
#if !defined (AM_BB_CMPLX_COS_SIGNAL_FP16_USE_NT_STORES)
#define AM_BB_CMPLX_COS_SIGNAL_FP16_USE_NT_STORES 0
#endif

#ifdef _OPENMP
// Default init a storage arrays for the first-touch (OpenMP) processing
#if !defined(AM_BB_CMPLX_COS_SIGNAL_FP16_INIT_STORAGE)
#define AM_BB_CMPLX_COS_SIGNAL_FP16_INIT_STORAGE 1
#endif 
#endif

#if (AM_BB_CMPLX_COS_SIGNAL_FP16_INIT_STORAGE) == 1
#define INIT_BY_STD_FILL 0
#endif 

// For inlining of trigo functions (asin,acos,sin,cos)
#if !defined(AM_BB_CMPLX_COS_SIGNAL_FP16_USE_CEPHES)
#define AM_BB_CMPLX_COS_SIGNAL_FP16_USE_CEPHES 1
#endif 

// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (AM_BB_CMPLX_COS_SIGNAL_FP16_USE_PMC_INSTRUMENTATION)
#define AM_BB_CMPLX_COS_SIGNAL_FP16_USE_PMC_INSTRUMENTATION 1
#endif 

#if (AM_BB_CMPLX_COS_SIGNAL_FP16_USE_PMC_INSTRUMENTATION) == 1
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

#if (AM_BB_CMPLX_COS_SIGNAL_FP16_USE_CEPHES) == 0
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
           enum class am_bb_cmplx_cos_signal_fp16_rand_gens : int32_t 
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

             enum class am_bb_cmplx_cos_signal_fp16_rand_distr : int32_t 
             {
                     uniform,
                     normal,
                     cauchy,
                     log_norm,
                     expo_gamma,
                     weibull,
                     gamma
             };

              /*
                Useful for noise added modulation
           */
             struct alignas(64) am_bb_cmplx_cos_signal_fp16_pdf_params_t 
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

            struct alignas(64) am_bb_cmplx_cos_signal_fp16_t final 
            {
                   
                   std::size_t       m_nsamples; 
                   std::uint32_t     m_nK; // number of inner loop wave components
                   float             m_n;  // number of 2PI cycles
                   float             m_A;  // An amplitude
                   float             m_P;  // A period 
                   darray_c2_t       m_sig_samples;

                   am_bb_cmplx_cos_signal_fp16_t() = default;

                   am_bb_cmplx_cos_signal_fp16_t(  const std::size_t,
                                                   const std::uint32_t,
                                                   const float,
                                                   const float,
                                                   const float) noexcept(false);

                   am_bb_cmplx_cos_signal_fp16_t(am_bb_cmplx_cos_signal_fp16_t &&) noexcept(true);

                   am_bb_cmplx_cos_signal_fp16_t(const am_bb_cmplx_cos_signal_fp16_t &) noexcept(false);

                   ~am_bb_cmplx_cos_signal_fp16_t() noexcept(false);

                   am_bb_cmplx_cos_signal_fp16_t & operator=(const am_bb_cmplx_cos_signal_fp16_t &) = delete;

                   am_bb_cmplx_cos_signal_fp16_t & operator=(am_bb_cmplx_cos_signal_fp16_t &&) noexcept(true);

                   void init_storage(const half_float::half);
                                     
                   
                   static void create_signal_plot( const std::uint32_t,
                                                   const half_float::half * __restrict,
                                                   const half_float::half * __restrict,
                                                   const std::string &,
                                                   const std::string &,
                                                   const bool );

                    __ATTR_ALWAYS_INLINE__
                    inline float cos_sample(const float t,
                                            const float invT) 
                    {
                          using namespace gms::math;
                          constexpr float C628318530717958647692528676656{6.28318530717958647692528676656f};
                          const float arg{this->m_n*C628318530717958647692528676656*t*invT};
                          float cos_val{};
#if (AM_BB_CMPLX_COS_SIGNAL_USE_CEPHES) == 1
                          cos_val = 0.5f*this->m_A*(1.0f-ceph_cosf(arg));
                          return (cos_val);
#else                     
                          cos_val = 0.5f*this->m_A*(1.0f-std::cos(arg));
                          return (cos_val);
#endif 
                    }

                    template<class Functor> float 
                    cos_sample_noisy(const float t,
                                     const float invT,
                                     const float scale,
                                     Functor &f)
                    {
                          using namespace gms::math;
                          constexpr float C628318530717958647692528676656{6.28318530717958647692528676656f};
                          const float arg{this->m_n*C628318530717958647692528676656*t*invT};
                          const float r_noise{f()*scale};
                          float cos_val{};
#if (AM_BB_CMPLX_COS_SIGNAL_FP16_USE_CEPHES) == 1
                          cos_val = 0.5f*this->m_A*(1.0f-ceph_cosf(arg));
                          cos_val += r_noise;
                          return (cos_val);
#else                     
                          cos_val = 0.5f*this->m_A*(1.0f-std::cos(arg));
                          cos_val += r_noise;
                          return (cos_val);
#endif 
                    }

                     __ATTR_ALWAYS_INLINE__ 
                    inline const half_float::half * get_I_channel() const noexcept
                    {
                         const half_float::half * __restrict__ p_I{&reinterpret_cast<half_float::half(&)[2]>(this->m_sig_samples.m_data[0])[0]};
                         return (p_I);
                    }

                    __ATTR_ALWAYS_INLINE__
                    inline const half_float::half * get_Q_channel() const noexcept
                    {
                         const half_float::half * __restrict__ p_Q{&reinterpret_cast<half_float::half(&)[2]>(this->m_sig_samples.m_data[0])[1]};
                         return (p_Q);
                    }

                    __ATTR_ALWAYS_INLINE__ 
                    inline half_float::half * get_I_channel() noexcept
                    {
                         half_float::half * __restrict__ p_I{&reinterpret_cast<half_float::half(&)[2]>(this->m_sig_samples.m_data[0])[0]};
                         return (p_I);
                    }

                    __ATTR_ALWAYS_INLINE__
                    inline  half_float::half * get_Q_channel() noexcept
                    {
                         half_float::half * __restrict__ p_Q{&reinterpret_cast<half_float::half(&)[2]>(this->m_sig_samples.m_data[0])[1]};
                         return (p_Q);
                    }


                      /*  Data symbol-transmitted*/
                    __ATTR_ALIGN__(32)
                    __ATTR_HOT__
#if defined (__INTEL_COMPILER) || defined (__ICC)
                    __ATTR_OPTIMIZE_03__
#endif
                    std::int32_t 
                    create_signal_user_data(const float * __restrict__,
                                            const std::uint32_t,
                                            const std::uint32_t);
                    
                    
                     /*  Data symbol-transmitted, outer loop unrolled four times.*/
                    __ATTR_ALIGN__(32)
                    __ATTR_HOT__
#if defined (__INTEL_COMPILER) || defined (__ICC)
                    __ATTR_OPTIMIZE_03__
#endif
                    std::int32_t 
                    create_signal_user_data_u4x(const float * __restrict__,
                                                const std::uint32_t,
                                                const std::uint32_t);

                      /*  Data symbol-transmitted*/
                    __ATTR_ALIGN__(32)
                    __ATTR_HOT__
#if defined (__INTEL_COMPILER) || defined (__ICC)
                    __ATTR_OPTIMIZE_03__
#endif
                    std::int32_t 
                    create_noisy_signal_user_data(am_bb_cmplx_cos_signal_fp16_pdf_params_t &,
                                                  const float,
                                                  const float * __restrict__,
                                                  const std::uint32_t,
                                                  const std::uint32_t);
                    __ATTR_ALIGN__(32)
                    __ATTR_HOT__
#if defined (__INTEL_COMPILER) || defined (__ICC)
                    __ATTR_OPTIMIZE_03__
#endif
                    std::int32_t 
                    create_noisy_signal_user_data(am_bb_cmplx_cos_signal_fp16_pdf_params_t &,
                                                  am_bb_cmplx_cos_signal_fp16_rand_distr,
                                                  const float,
                                                  const float * __restrict__,
                                                  const std::uint32_t,
                                                  const std::uint32_t);


                     

                   
            };

            auto operator<<(std::ostream &,
                            am_bb_cmplx_cos_signal_fp16_t &)->std::ostream &;
            

}

}

















#endif /*__GMS_AM_BB_CMPLX_COS_SIGNAL_FP16_H__*/
