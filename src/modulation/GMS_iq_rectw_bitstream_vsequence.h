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

#ifndef __GMS_IQ_RECTW_BITSTREAM_VSEQUENCE_H__
#define __GMS_IQ_RECTW_BITSTREAM_VSEQUENCE_H__ 241220250123


namespace file_info 
{

     static const unsigned int GMS_IQ_RECTW_BITSTREAM_VSEQUENCE_MAJOR = 1;
     static const unsigned int GMS_IQ_RECTW_BITSTREAM_VSEQUENCE_MINOR = 1;
     static const unsigned int GMS_IQ_RECTW_BITSTREAM_VSEQUENCE_MICRO = 0;
     static const unsigned int GMS_IQ_RECTW_BITSTREAM_VSEQUENCE_FULLVER =
       1000U*GMS_IQ_RECTW_BITSTREAM_VSEQUENCE_MAJOR+100U*GMS_IQ_RECTW_BITSTREAM_VSEQUENCE_MINOR+
       10U*GMS_IQ_RECTW_BITSTREAM_VSEQUENCE_MICRO;
     static const char GMS_IQ_RECTW_BITSTREAM_VSEQUENCE_CREATION_DATE[] = "24-12-2025 01:23PM +00200 (WED 24 NOV 2025 GMT+2)";
     static const char GMS_IQ_RECTW_BITSTREAM_VSEQUENCE_BUILD_DATE[]    = __DATE__; 
     static const char GMS_IQ_RECTW_BITSTREAM_VSEQUENCE_BUILD_TIME[]    = __TIME__;
     static const char GMS_IQ_RECTW_BITSTREAM_VSEQUENCE_SYNOPSIS[]      = "I/Q Rectwave bitstream signal sequence (vertically summed).";

}

#include <cstdint>
#include <string>
#include <vector>
#include "GMS_config.h"
#include "GMS_iq_rectwave_bitstream.h"


#if !defined(IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH)
#define IQ_RECTW_BITSTREAM_VSEQUENCE_SOFT_PREFETCH 1 
#endif

// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (IQ_RECTW_BITSTREAM_VSEQUENCE_USE_PMC_INSTRUMENTATION)
#define IQ_RECTW_BITSTREAM_VSEQUENCE_USE_PMC_INSTRUMENTATION 0
#endif 

#if (IQ_RECTW_BITSTREAM_VSEQUENCE_USE_PMC_INSTRUMENTATION) == 1
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

          struct alignas(64) I_rectw_bitstream_pdf_params_t 
          {
                  float uni_real_a_r;  
                  float uni_real_b_r;  
                  float norm_mu_r;    
                  float norm_sigma_r;                
                  float cauchy_a_r;                      
                  float cauchy_b_r;                      
                  float log_norm_m_r;                      
                  float log_norm_s_r;                     
                  float expo_gamma_r;                      
                  float weibull_a_r;                     
                  float weibull_b_r;                      
                  float gamma_alph_r;                      
                  float gamma_bet_r;                     
                  int32_t poisson_mu;
          };

          struct alignas(32) Q_rectw_bitstream_pdf_params_t 
          {
                 float uni_real_a_i;
                 float uni_real_b_i;
                 float norm_mu_i;
                 float norm_sigma_i;
                 float cauchy_a_i;
                 float cauchy_b_i;
                 float log_norm_m_i;
                 float log_norm_s_i;
                 float expo_gamma_i;
                 float weibull_a_i;
                 float weibull_b_i;
                 float gamma_alph_i;
                 float gamma_bet_i;
                 int32_t poisson_mu;
          };

          struct alignas(64) iq_rectw_bitstream_vsequence_t final 
          {
                 std::size_t                                   m_I_nsamples;
                 std::size_t                                   m_IQ_nsignals;                
                 darray_r4_t                                   m_I_k_values;
                 darray_r4_t                                   m_I_T_values;
                 darray_r4_t                                   m_I_Ac_values;
                 darray_r4_t                                   m_I_fc_values;
                 darray_r4_t                                   m_I_ph0_values;   
                 std::size_t                                   m_Q_nsamples;  
                 darray_r4_t                                   m_Q_k_values;
                 darray_r4_t                                   m_Q_T_values;
                 darray_r4_t                                   m_Q_Ac_values;
                 darray_r4_t                                   m_Q_fc_values;
                 darray_r4_t                                   m_Q_ph0_values;
                 std::vector<iq_rectwave_bitstream_t>          m_iq_rectw_bitstreams;
                 std::vector<I_rectw_bitstream_pdf_params_t>   m_I_rw_bitstream_pdf_params;
                 std::vector<Q_rectw_bitstream_pdf_params_t>   m_Q_rw_bitstream_pdf_params;
                 darray_r4_t                                   m_I_vsequence;
                 darray_r4_t                                   m_Q_vsequence;
                 
                  
                  enum class IQ_rectw_bitstreams_pdf_distros : int32_t 
                  {
                       UNIFORM,
                       NORMAL,
                       CAUCHY,
                       LOG_NORMAL,
                       EXPONENTIAL_GAMMA,
                       WEIBULL,
                       GAMMA,
                       POISSON
                  };

                  enum class IQ_rectw_bitstreams_optim_path : int32_t 
                  {
                       SCALAR_PATH,
                       SSE_VEC_PATH,
                       SSE_U4X_VEC_PATH,
                       AVX_VEC_PATH,
                       AVX_U4X_VEC_PATH,
                       AVX512_VEC_PATH,
                       AVX512_U4X_VEC_PATH
                  };

                  //Optional
                  enum class IQ_rectw_bitstream_max_signal_components : std::size_t 
                  {
                       MAX_SIGNAL_COMPONENTS_25  = 25ull,
                       MAX_SIGNAL_COMPONENTS_50  = 50ull,
                       MAX_SIGNAL_COMPONENTS_100 = 100ull,
                       MAX_SIGNAL_COMPONENTS_150 = 150ull,
                       MAX_SIGNAL_COMPONENTS_200 = 200ull,
                       MAX_SIGNAL_COMPONENTS_250 = 250ull
                  };

                  iq_rectw_bitstream_vsequence_t() noexcept(true);

                  iq_rectw_bitstream_vsequence_t(const std::size_t,
                                                 const std::size_t,
                                                 const darray_r4_t &,
                                                 const darray_r4_t &,
                                                 const darray_r4_t &,
                                                 const darray_r4_t &,
                                                 const darray_r4_t &,
                                                 const std::size_t,
                                                 const darray_r4_t &,
                                                 const darray_r4_t &,
                                                 const darray_r4_t &,
                                                 const darray_r4_t &,
                                                 const darray_r4_t &,
                                                 const std::vector<I_rectw_bitstream_pdf_params_t> &,
                                                 const std::vector<Q_rectw_bitstream_pdf_params_t> &) noexcept(false);

                  iq_rectw_bitstream_vsequence_t(const iq_rectw_bitstream_vsequence_t &) noexcept(false);

                  iq_rectw_bitstream_vsequence_t(iq_rectw_bitstream_vsequence_t &&) noexcept(true);

                  ~iq_rectw_bitstream_vsequence_t();

                  iq_rectw_bitstream_vsequence_t &
                                     operator=(const iq_rectw_bitstream_vsequence_t &) noexcept(false);

                  iq_rectw_bitstream_vsequence_t &
                                     operator=(iq_rectw_bitstream_vsequence_t &&) noexcept(true);

                  static void create_sequence_plot(  const std::uint32_t,
                                                   const float * __restrict,
                                                   const float * __restrict,
                                                   const std::string &,
                                                   const std::string &,
                                                   const bool );
                  
                  /*Important!!: -- The size of two input data arrays must be equal to m_IQ_signals and must not exceed the 
                                  template argument max_signals*/

                  template<IQ_rectw_bitstreams_optim_path  optim_path,
                           IQ_rectw_bitstream_max_signal_components max_signals>         
                  std::int32_t 
                  generate_iq_rectw_bitstream_sequence_sse_u16x(std::int32_t * __restrict__,
                                                                std::int32_t * __restrict__);
                                                                
                  template<IQ_rectw_bitstreams_optim_path  optim_path,
                           IQ_rectw_bitstream_max_signal_components max_signals>                          
                  std::int32_t 
                  generate_iq_rectw_bitstream_sequence_sse_u10x(std::int32_t * __restrict__,
                                                                std::int32_t * __restrict__);
                                                                

                  template<IQ_rectw_bitstreams_optim_path  optim_path,
                           IQ_rectw_bitstream_max_signal_components max_signals>                          
                  std::int32_t 
                  generate_iq_rectw_bitstream_sequence_sse_u6x(std::int32_t * __restrict__,
                                                               std::int32_t * __restrict__);
                                                              

                  template<IQ_rectw_bitstreams_optim_path  optim_path,
                           IQ_rectw_bitstream_max_signal_components max_signals>                          
                  std::int32_t 
                  generate_iq_rectw_bitstream_sequence_avx_u16x(std::int32_t * __restrict__,
                                                                std::int32_t * __restrict__);
                                                                
                  template<IQ_rectw_bitstreams_optim_path  optim_path,
                           IQ_rectw_bitstream_max_signal_components max_signals>                          
                  std::int32_t 
                  generate_iq_rectw_bitstream_sequence_avx_u10x(std::int32_t * __restrict__,
                                                                std::int32_t * __restrict__);
                                                                
                  template<IQ_rectw_bitstreams_optim_path  optim_path,
                           IQ_rectw_bitstream_max_signal_components max_signals>                         
                  std::int32_t 
                  generate_iq_rectw_bitstream_sequence_avx_u6x(std::int32_t * __restrict__,
                                                               std::int32_t * __restrict__);
                                                               
                  template<IQ_rectw_bitstreams_optim_path  optim_path,
                           IQ_rectw_bitstream_max_signal_components max_signals>                          
                  std::int32_t 
                  generate_iq_rectw_bitstream_sequence_avx512_u16x(std::int32_t * __restrict__,
                                                                   std::int32_t * __restrict__);
                                                                   
                  template<IQ_rectw_bitstreams_optim_path  optim_path,
                           IQ_rectw_bitstream_max_signal_components max_signals>                           
                  std::int32_t 
                  generate_iq_rectw_bitstream_sequence_avx512_u10x(std::int32_t * __restrict__,
                                                                   std::int32_t * __restrict__);
                                                                   
                  template<IQ_rectw_bitstreams_optim_path  optim_path,
                           IQ_rectw_bitstream_max_signal_components max_signals>                          
                  std::int32_t 
                  generate_iq_rectw_bitstream_sequence_avx512_u6x(std::int32_t * __restrict__,
                                                                  std::int32_t * __restrict__);
                                                                  

          };

          

}// radiolocation

}// gms












































#endif /*__GMS_IQ_RECTW_BITSTREAM_VSEQUENCE_H__*/
