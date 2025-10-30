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

#ifndef __GMS_PM_BB_CMPLX_SINE_SIGNAL_H__
#define __GMS_PM_BB_CMPLX_SINE_SIGNAL_H__ 301020250928

namespace file_info 
{

     static const unsigned int GMS_PM_BB_CMPLX_SINE_SIGNAL_MAJOR = 1;
     static const unsigned int GMS_PM_BB_CMPLX_SINE_SIGNAL_MINOR = 1;
     static const unsigned int GMS_PM_BB_CMPLX_SINE_SIGNAL_MICRO = 0;
     static const unsigned int GMS_PM_BB_CMPLX_SINE_SIGNAL_FULLVER =
       1000U*GMS_PM_BB_CMPLX_SINE_SIGNAL_MAJOR+100U*GMS_PM_BB_CMPLX_SINE_SIGNAL_MINOR+
       10U*GMS_PM_BB_CMPLX_SINE_SIGNAL_MICRO;
     static const char GMS_PM_BB_CMPLX_SINE_SIGNAL_CREATION_DATE[] = "30-10-2025 09:29AM +00200 (THR 30 OCT 2025 GMT+2)";
     static const char GMS_PM_BB_CMPLX_SINE_SIGNAL_BUILD_DATE[]    = __DATE__; 
     static const char GMS_PM_BB_CMPLX_SINE_SIGNAL_BUILD_TIME[]    = __TIME__;
     static const char GMS_PM_BB_CMPLX_SINE_SIGNAL_SYNOPSIS[]      = "PM baseband complex sine wave signal.";

}

#include <cstdint>
#include <string>
#include <iostream>
#include "GMS_config.h"
#include "GMS_dyn_array.h"

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
// To be added.
#if !defined (PM_BB_CMPLX_SINE_SIGNAL_USE_NT_STORES)
#define PM_BB_CMPLX_SINE_SIGNAL_USE_NT_STORES 0
#endif

#ifdef _OPENMP
// Default init a storage arrays for the first-touch (OpenMP) processing
#if !defined(PM_BB_CMPLX_SINE_SIGNAL_INIT_STORAGE)
#define PM_BB_CMPLX_SINE_SIGNAL_INIT_STORAGE 1
#endif 
#endif

#if (PM_BB_CMPLX_SINE_SIGNAL_INIT_STORAGE) == 1
#define INIT_BY_STD_FILL 0
#endif 

// For inlining of trigo functions (asin,acos,sin,cos)
#if !defined(PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES)
#define PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES 1
#endif 

// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (PM_BB_CMPLX_SINE_SIGNAL_USE_PMC_INSTRUMENTATION)
#define PM_BB_CMPLX_SINE_SIGNAL_USE_PMC_INSTRUMENTATION 1
#endif 

#if (PM_BB_CMPLX_SINE_SIGNAL_USE_PMC_INSTRUMENTATION) == 1
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

#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 0
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
           enum class pm_bb_cmplx_sin_signal_rand_gens : int32_t 
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

             enum class pm_bb_cmplx_sin_signal_rand_distr : int32_t 
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
             struct alignas(64) pm_bb_cmplx_sin_signal_pdf_params_t 
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

            struct alignas(64) pm_bb_cmplx_sin_signal_t final 
            {
                   std::size_t       m_nsamples; 
                   std::uint32_t     m_nK; // number of inner loop wave components
                    // subcarrier center frequency
                   float             m_w0;
                   float             m_ph; // initial phase offset.
                   float             m_n;  // number of 2PI cycles
                   float             m_A;  // An amplitude
                   float             m_P;  // A period 
                   darray_c4_t       m_sig_samples;

                   pm_bb_cmplx_sin_signal_t() = default;

                   pm_bb_cmplx_sin_signal_t(  const std::size_t,
                                              const std::uint32_t,
                                              const float,
                                              const float,
                                              const float,
                                              const float,
                                              const float) noexcept(false);

                   pm_bb_cmplx_sin_signal_t(pm_bb_cmplx_sin_signal_t &&) noexcept(true);

                   pm_bb_cmplx_sin_signal_t(const pm_bb_cmplx_sin_signal_t &) = delete;

                   ~pm_bb_cmplx_sin_signal_t() noexcept(false);

                   pm_bb_cmplx_sin_signal_t & operator=(const pm_bb_cmplx_sin_signal_t &) = delete;

                   pm_bb_cmplx_sin_signal_t & operator=(pm_bb_cmplx_sin_signal_t &&) noexcept(true);

                   void init_storage(const float);
                                     

                   static void create_signal_plot( const std::uint32_t,
                                                   const float * __restrict,
                                                   const float * __restrict,
                                                   const std::string &,
                                                   const std::string &,
                                                   const bool );

                    __ATTR_ALWAYS_INLINE__
                    inline float sin_sample(const float t,
                                            const float invT) 
                    {
                          using namespace gms::math;
                          constexpr float C314159265358979323846264338328{3.14159265358979323846264338328f};
                          const float arg{this->m_n*C314159265358979323846264338328*t*invT};
                          float sin_val{};
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1
                          sin_val = this->m_A*ceph_sinf(arg);
                          return (sin_val);
#else                     
                          sin_val = this->m_A*std::sin(arg);
                          return (sin_val);
#endif 
                    }

                    template<class Functor> float 
                    sin_sample_noisy(const float t,
                                     const float invT,
                                     const float scale,
                                     Functor &f)
                    {
                          using namespace gms::math;
                          constexpr float C314159265358979323846264338328{3.14159265358979323846264338328f};
                          const float arg{this->m_n*C314159265358979323846264338328*t*invT};
                          const float r_noise{f()*scale};
                          float sin_val{};
                          this->m_A += r_noise;
#if (PM_BB_CMPLX_SINE_SIGNAL_USE_CEPHES) == 1
                          sin_val = this->m_A*ceph_sinf(arg);
                          return (sin_val);
#else                     
                          sin_val = this->m_A*std::sin(arg);
                          return (sin_val);
#endif 
                    }

                     __ATTR_ALWAYS_INLINE__ 
                    inline const float * get_I_channel() const noexcept
                    {
                         const float * __restrict__ p_I{&reinterpret_cast<float(&)[2]>(this->m_sig_samples.m_data[0])[0]};
                         return (p_I);
                    }

                    __ATTR_ALWAYS_INLINE__
                    inline const float * get_Q_channel() const noexcept
                    {
                         const float * __restrict__ p_Q{&reinterpret_cast<float(&)[2]>(this->m_sig_samples.m_data[0])[1]};
                         return (p_Q);
                    }

                    __ATTR_ALWAYS_INLINE__ 
                    inline float * get_I_channel() noexcept
                    {
                         float * __restrict__ p_I{&reinterpret_cast<float(&)[2]>(this->m_sig_samples.m_data[0])[0]};
                         return (p_I);
                    }

                    __ATTR_ALWAYS_INLINE__
                    inline  float * get_Q_channel() noexcept
                    {
                         float * __restrict__ p_Q{&reinterpret_cast<float(&)[2]>(this->m_sig_samples.m_data[0])[1]};
                         return (p_Q);
                    }


                      /*  Data symbol-transmitted*/
                    std::int32_t 
                    create_signal_user_data(const float * __restrict__,
                                            const std::uint32_t,
                                            const std::uint32_t);
                    
                    
                     /*  Data symbol-transmitted, outer loop unrolled four times.*/
                    std::int32_t 
                    create_signal_user_data_u4x(const float * __restrict__,
                                                const std::uint32_t,
                                                const std::uint32_t);

                      /*  Data symbol-transmitted*/
                    std::int32_t 
                    create_noisy_signal_user_data(pm_bb_cmplx_sin_signal_pdf_params_t &,
                                                  const float,
                                                  const float * __restrict__,
                                                  const std::uint32_t,
                                                  const std::uint32_t);

                    std::int32_t 
                    create_noisy_signal_user_data(pm_bb_cmplx_sin_signal_pdf_params_t &,
                                                  pm_bb_cmplx_sin_signal_rand_distr,
                                                  const float,
                                                  const float * __restrict__,
                                                  const std::uint32_t,
                                                  const std::uint32_t);


                     

                   
            };

            auto operator<<(std::ostream &,
                            pm_bb_cmplx_sin_signal_t &)->std::ostream &;
            

}

}

















#endif /*__GMS_PM_BB_CMPLX_SINE_SIGNAL_H__*/