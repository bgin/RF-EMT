

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

#ifndef __GMS_TRAPEZOID_WAVEFORM_H__
#define __GMS_TRAPEZOID_WAVEFORM_H__ 100920250803


namespace file_info 
{

     static const unsigned int GMS_TRAPEZOID_WAVEFORM_MAJOR = 1;
     static const unsigned int GMS_TRAPEZOID_WAVEFORM_MINOR = 1;
     static const unsigned int GMS_TRAPEZOID_WAVEFORM_MICRO = 0;
     static const unsigned int GMS_TRAPEZOID_WAVEFORM_FULLVER =
       1000U*GMS_TRAPEZOID_WAVEFORM_MAJOR+100U*GMS_TRAPEZOID_WAVEFORM_MINOR+
       10U*GMS_TRAPEZOID_WAVEFORM_MICRO;
     static const char GMS_TRAPEZOID_WAVEFORM_CREATION_DATE[] = "10-09-2025 08:03 +00200 (WED  10 SEP 2025 GMT+2)";
     static const char GMS_TRAPEZOID_WAVEFORM_BUILD_DATE[]    = __DATE__; 
     static const char GMS_TRAPEZOID_WAVEFORM_BUILD_TIME[]    = __TIME__;
     static const char GMS_TRAPEZOID_WAVEFORM_SYNOPSIS[]      = "Trapezoid waveform generators.";

}

#include <cstdint>
#include <string>
#include "GMS_config.h"
#include "GMS_dyn_array.h"

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
// To be added.
#if !defined (TRAPEZOID_WAVEFORM_USE_NT_STORES)
#define TRAPEZOID_WAVEFORM_USE_NT_STORES 0
#endif

#ifdef _OPENMP
// Default init a storage arrays for the first-touch (OpenMP) processing
#if !defined(TRAPEZOID_WAVEFORM_INIT_STORAGE)
#define TRAPEZOID_WAVEFORM_INIT_STORAGE 1
#endif 
#endif

#if (TRAPEZOID_WAVEFORM_INIT_STORAGE) == 1
#define INIT_BY_STD_FILL 0
#endif 

// For inlining of trigo functions (asin,acos,sin,cos)
#if !defined(TRAPEZOID_WAVEFORM_USE_CEPHES)
#define TRAPEZOID_WAVEFORM_USE_CEPHES 1
#endif 

// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (TRAPEZOID_WAVEFORM_USE_PMC_INSTRUMENTATION)
#define TRAPEZOID_WAVEFORM_USE_PMC_INSTRUMENTATION 1
#endif 

#if (TRAPEZOID_WAVEFORM_USE_PMC_INSTRUMENTATION) == 1
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


namespace gms 
{

namespace radiolocation
{
            

             enum class trapezw_rand_gens : int32_t 
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

             struct alignas(64) trapezw_pdf_params_t 
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

            struct alignas(64) trapezoid_waveform_t final 
            {
                   
                   std::size_t                      __n_samples__;
                   std::uint32_t                    __n_waves__;
                   std::uint32_t                    __n_param_a__;
                   std::uint32_t                    __n_param_l__;
                   std::uint32_t                    __n_param_c__;
                   std::uint32_t                    __n_param_m__;
                   darray_r4_t                      __trapezw_samples__;

                   trapezoid_waveform_t() = delete;

                   trapezoid_waveform_t(const std::size_t,
                                        const std::uint32_t,
                                        const std::uint32_t,
                                        const std::uint32_t,
                                        const std::uint32_t,
                                        const std::uint32_t) noexcept(false);

                    trapezoid_waveform_t(trapezoid_waveform_t &&) noexcept(true);

                   ~trapezoid_waveform_t() noexcept(false);

                    trapezoid_waveform_t & operator=(const trapezoid_waveform_t &) = delete;

                    trapezoid_waveform_t & operator=(trapezoid_waveform_t &&);

                    void init_storage(const float);

                    static void create_signal_plot(const std::uint32_t,
                                                   const float * __restrict,
                                                   const float * __restrict,
                                                   const std::string &,
                                                   const std::string &,
                                                   const bool );

                    /*Create a single sample of trapezoid wave*/
                    float sample_of_trapezoid_wave(const float,
                                                   const float,
                                                   const float,
                                                   const float,
                                                   const float);

                    /* Create single trapezoid waveform*/
                    void single_trapezoid_wave(const float,
                                               const float,
                                               const float,
                                               const float);

                  
                     /* Create series of trapezoid waves (shaping the curve)*/
                    void series_of_trapezoid_waves(const float,
                                                   const float,
                                                   const float,
                                                   const float,
                                                   const std::uint32_t);
                     
                       /* Create series of trapezoid waves, unroll the outer loop 2 times*/
                    void series_of_trapezoid_waves_u2x(const float,
                                                       const float,
                                                       const float,
                                                       const float);
                    
                      /* Create series of trapezoid waves, unroll the outer loop 4 times*/
                    void series_of_trapezoid_waves_u4x(const float,
                                                       const float,
                                                       const float,
                                                       const float);
                   
                    /*
                      Create single trapezoid wave with a,l,c,m parameters
                      monotonically growing. and (curve shaped)
                    */
                   void single_trapezoid_wave_v2(const float,
                                                 const float,
                                                 const float,
                                                 const float,
                                                 const std::uint32_t);

                    /*
                      Create single trapezoid wave with a,l,c,m parameters
                      caller passed (curve shaped)
                    */
                   void single_trapezoid_wave_v2(const float * __restrict,
                                                 const float * __restrict,
                                                 const float * __restrict,
                                                 const float * __restrict,
                                                 const std::uint32_t);

                                                    

                    
            };
}

}






























#endif /*__GMS_TRAPEZOID_WAVEFORM_H__*/