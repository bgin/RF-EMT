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

#ifndef __GMS_SINUSOIDAL_WEIGHTED_OQPSK_H__
#define __GMS_SINUSOIDAL_WEIGHTED_OQPSK_H__ 251120250638

namespace file_info 
{

     static const unsigned int GMS_SINUSOIDAL_WEIGHTED_OQPSK_MAJOR = 1;
     static const unsigned int GMS_SINUSOIDAL_WEIGHTED_OQPSK_MINOR = 1;
     static const unsigned int GMS_SINUSOIDAL_WEIGHTED_OQPSK_MICRO = 0;
     static const unsigned int GMS_SINUSOIDAL_WEIGHTED_OQPSK_FULLVER =
       1000U*GMS_SINUSOIDAL_WEIGHTED_OQPSK_MAJOR+100U*GMS_SINUSOIDAL_WEIGHTED_OQPSK_MINOR+
       10U*GMS_SINUSOIDAL_WEIGHTED_OQPSK_MICRO;
     static const char GMS_SINUSOIDAL_WEIGHTED_OQPSK_CREATION_DATE[] = "25-11-2025 06:38AM +00200 (TUE 25 NOV 2025 GMT+2)";
     static const char GMS_SINUSOIDAL_WEIGHTED_OQPSK_BUILD_DATE[]    = __DATE__; 
     static const char GMS_SINUSOIDAL_WEIGHTED_OQPSK_BUILD_TIME[]    = __TIME__;
     static const char GMS_SINUSOIDAL_WEIGHTED_OQPSK_SYNOPSIS[]      = "SInusoidally weighted OQPSK signal representation.";

}

#include <cstdint>
#include <string>
#include <iostream>
#include <thread>
#include "GMS_config.h"
#include "GMS_dyn_array.h"

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
// To be added.
#if !defined (SINUSOIDAL_WEIGHTED_OQPSK_USE_NT_STORES)
#define SINUSOIDAL_WEIGHTED_OQPSK_USE_NT_STORES 0
#endif

#ifdef _OPENMP
// Default init a storage arrays for the first-touch (OpenMP) processing
#if !defined(SINUSOIDAL_WEIGHTED_OQPSK_INIT_STORAGE)
#define SINUSOIDAL_WEIGHTED_OQPSK_INIT_STORAGE 1
#endif 
#endif

#if (SINUSOIDAL_WEIGHTED_OQPSK_INIT_STORAGE) == 1
#define INIT_BY_STD_FILL 0
#endif 

// For inlining of trigo functions (asin,acos,sin,cos)
#if !defined(SINUSOIDAL_WEIGHTED_OQPSK_USE_CEPHES)
#define SINUSOIDAL_WEIGHTED_OQPSK_USE_CEPHES 1
#endif 

// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (SINUSOIDAL_WEIGHTED_OQPSK_USE_PMC_INSTRUMENTATION)
#define SINUSOIDAL_WEIGHTED_OQPSK_USE_PMC_INSTRUMENTATION 0
#endif 

#if (SINUSOIDAL_WEIGHTED_OQPSK_USE_PMC_INSTRUMENTATION) == 1
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

#if (SINUSOIDAL_WEIGHTED_OQPSK_USE_CEPHES) == 0
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
           enum class sinusoidal_weighted_oqpsk_rand_gens : int32_t 
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

             enum class sinusoidal_weighted_oqpsk_rand_distr : int32_t 
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
             struct alignas(64) sinusoidal_weighted_oqpsk_pdf_params_t 
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

           struct alignas(64) sinusoidal_weighted_oqpsk_t final 
           {
                  std::size_t     m_nsamples;
                  std::size_t     m_I_ch_nsamples; // I-channel i.e. rect wave number of samples
                  std::size_t     m_Q_ch_nsamples; // Q-channel i.e. rec wave number of samples
                  float           m_A_I; // I channel amplitude
                  float           m_A_Q; // Q channel amplitude
                  float           m_T;   // period
                  float           m_cw0; // cosine carrier central frequency
                  float           m_sw0; // sine carrier central frequency
                  float           m_cph0;// cosine carrier phase offset
                  float           m_sph0;// sine carrier phase offset
                  darray_r4_t     m_I_bitstream; //user data, i.e. [-1,1]
                  darray_r4_t     m_Q_bitstream; //user data, i.e. [-1,1]
                  darray_r4_t     m_I_channel;
                  darray_r4_t     m_Q_channel;
                  darray_r4_t     m_msk_signal; // summed waveform of I/Q channels
                  
                  sinusoidal_weighted_oqpsk_t() = delete;

                  sinusoidal_weighted_oqpsk_t(const std::size_t,
                                              const std::size_t,
                                              const std::size_t,
                                              const float,
                                              const float,
                                              const float,
                                              const float,
                                              const float,
                                              const float,
                                              const float) noexcept(false);

                  sinusoidal_weighted_oqpsk_t(const sinusoidal_weighted_oqpsk_t &) noexcept(false);

                  sinusoidal_weighted_oqpsk_t(sinusoidal_weighted_oqpsk_t &&) noexcept(true);

                  ~sinusoidal_weighted_oqpsk_t() noexcept(true);

                  sinusoidal_weighted_oqpsk_t &
                             operator=(const sinusoidal_weighted_oqpsk_t &) noexcept(false);

                  sinusoidal_weighted_oqpsk_t &
                             operator=(sinusoidal_weighted_oqpsk_t &&) noexcept(true);

                                               

                  static void create_signal_plot( const std::uint32_t,
                                                   const float * __restrict,
                                                   const float * __restrict,
                                                   const std::string &,
                                                   const std::string &,
                                                   const bool );

                  private:
                  
                  __ATTR_ALWAYS_INLINE__
                  inline 
                  void create_random_bitstream() noexcept(false)
                  {
                       std::size_t i,j;
                       thread_local unsigned long long I_seed{0ull};
                       thread_local unsigned long long Q_seed{0ull};
                       thread_local std::uniform_real_distribution<float> rand_I_bits;
                       thread_local std::uniform_real_distribution<float> rand_Q_bits;
                       float r_I_0{0.0f};
                       float r_I_1{r_I_0};
                       float r_I_2{r_I_0};
                       float r_I_3{r_I_0};
                       float r_Q_0{r_I_0};
                       float r_Q_1{r_I_0};
                       float r_Q_2{r_I_0};
                       float r_Q_3{r_I_0};
                       
                       I_seed = __rdtsc();
                       Q_seed = __rdtsc();
                       rand_I_bits = std::uniform_real_distribution<float>(-1.0f,1.0f);
                       rand_Q_bits = std::uniform_real_distribution<float>(-1.0f,1.0f);
                       auto rseed_I{std::mt19937(I_seed)};
                       auto rseed_Q{std::mt19937(Q_seed)};
                       if(__builtin_expect(this->m_I_ch_nsamples==this->m_Q_ch_nsamples,1))
                       {
                           for(i = 0ull; i != ROUND_TO_FOUR(this->m_I_ch_nsamples,4ull); i += 4ull) 
                           {
                               r_I_0 = rand_I_bits.operator()(rseed_I);
                               this->m_I_bitstream.m_data[i+0ull] = r_I_0;
                               r_Q_0 = rand_Q_bits.operator()(rseed_Q);
                               this->m_Q_bitstream.m_data[i+0ull] = r_Q_0;
                               r_I_1 = rand_I_bits.operator()(rseed_I);
                               this->m_I_bitstream.m_data[i+1ull] = r_I_1;
                               r_Q_1 = rand_Q_bits.operator()(rseed_Q);
                               this->m_Q_bitstream.m_data[i+1ull] = r_Q_1;
                               r_I_2 = rand_I_bits.operator()(rseed_I);
                               this->m_I_bitstream.m_data[i+2ull] = r_I_2;
                               r_Q_2 = rand_Q_bits.operator()(rseed_Q);
                               this->m_Q_bitstream.m_data[i+2ull] = r_Q_2;
                               r_I_3 = rand_I_bits.operator()(rseed_I);
                               this->m_I_bitstream.m_data[i+3ull] = r_I_3;
                               r_Q_3 = rand_Q_bits.operator()(rseed_Q);
                               this->m_Q_bitstream.m_data[i+3ull] = r_Q_3;
                           }

                           for(j = i; j != this->m_I_ch_nsamples; ++j) 
                           {
                               r_I_0 = rand_I_bits.operator()(rseed_I);
                               this->m_I_bitstream.m_data[j] = r_I_0;
                               r_Q_0 = rand_Q_bits.operator()(rseed_Q);
                               this->m_Q_bitstream.m_data[j] = r_Q_0;
                           }
                       }
                       else 
                       {
                           for(i = 0ull; i != ROUND_TO_FOUR(this->m_I_ch_nsamples,4ull); i += 4ull) 
                           {
                               r_I_0 = rand_I_bits.operator()(rseed_I);
                               this->m_I_bitstream.m_data[i+0ull] = r_I_0;
                               r_I_1 = rand_I_bits.operator()(rseed_I);
                               this->m_I_bitstream.m_data[i+1ull] = r_I_1;
                               r_I_2 = rand_I_bits.operator()(rseed_I);
                               this->m_I_bitstream.m_data[i+2ull] = r_I_2;
                               r_I_3 = rand_I_bits.operator()(rseed_I);
                               this->m_I_bitstream.m_data[i+3ull] = r_I_3;
                           }

                           for(j = i; j != this->m_I_ch_nsamples; ++j) 
                           {
                               r_I_0 = rand_I_bits.operator()(rseed_I);
                               this->m_I_bitstream.m_data[j] = r_I_0;
                           }

                           for(i = 0ull; i != ROUND_TO_FOUR(this->m_Q_ch_nsamples,4ull); i += 4ull) 
                           {
                               r_Q_0 = rand_Q_bits.operator()(rseed_Q);
                               this->m_Q_bitstream.m_data[i+0ull] = r_I_0;
                               r_Q_1 = rand_Q_bits.operator()(rseed_Q);
                               this->m_Q_bitstream.m_data[i+1ull] = r_I_1;
                               r_Q_2 = rand_Q_bits.operator()(rseed_Q);
                               this->m_Q_bitstream.m_data[i+2ull] = r_I_2;
                               r_I_3 = rand_Q_bits.operator()(rseed_Q);
                               this->m_Q_bitstream.m_data[i+3ull] = r_I_3;
                           }

                           for(j = i; j != this->m_Q_ch_nsamples; ++j) 
                           {
                               r_I_0 = rand_Q_bits.operator()(rseed_Q);
                               this->m_Q_bitstream.m_data[j] = r_I_0;
                           }
                       }
                  }

                  public:
          __ATTR_HOT__
          __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
          __ATTR_OPTIMIZE_03__
#endif 
                  std::int32_t oqpsk_signal_user_data(const darray_r4_t &,
                                                      const darray_r4_t &);

          __ATTR_HOT__
          __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
          __ATTR_OPTIMIZE_03__
#endif 
                  std::int32_t oqpsk_signal_noisy_user_data(const darray_r4_t &,
                                                            const darray_r4_t &,
                                                            sinusoidal_weighted_oqpsk_pdf_params_t &,
                                                            sinusoidal_weighted_oqpsk_rand_distr);


          __ATTR_HOT__
          __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
          __ATTR_OPTIMIZE_03__
#endif 
                  std::int32_t oqpsk_signal_random_data();

          __ATTR_HOT__
          __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
          __ATTR_OPTIMIZE_03__
#endif 
                  std::int32_t oqpsk_signal_noisy_random_data(sinusoidal_weighted_oqpsk_pdf_params_t &,
                                                              sinusoidal_weighted_oqpsk_rand_distr);                                 
           };

           auto 
           operator<<(std::ostream &,
                      const sinusoidal_weighted_oqpsk_t &)->std::ostream &;

} // radiolocation


} // gms

















#endif /*__GMS_SINUSOIDAL_WEIGHTED_OQPSK_H__*/