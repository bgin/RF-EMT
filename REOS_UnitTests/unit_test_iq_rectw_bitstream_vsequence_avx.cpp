
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <vector>
#include <thread>
#include "GMS_dyn_array.h"
#include "GMS_iq_rectw_bitstream_vsequence.h"

/*
   icpc -o unit_test_iq_rectwave_bitstream_vsequence_avx -O3 -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_loop_indices_LUT.h GMS_loop_indices_LUT.cpp GMS_iq_rectwave_bitstream.h GMS_iq_rectwave_bitstream.cpp GMS_iq_rectw_bitstream_vsequence.h GMS_iq_rectw_bitstream_vsequence.cpp unit_test_iq_rectw_bitstream_vsequence_avx.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_loop_indices_LUT.h GMS_loop_indices_LUT.cpp GMS_iq_rectwave_bitstream.h GMS_iq_rectwave_bitstream.cpp GMS_iq_rectw_bitstream_vsequence.h GMS_iq_rectw_bitstream_vsequence.cpp unit_test_iq_rectw_bitstream_vsequence_avx.cpp

*/

__attribute__((hot))
__attribute__((noinline))
void 
unit_test_iq_rectw_bitstream_vsequence_avx_u16x();

void 
unit_test_iq_rectw_bitstream_vsequence_avx_u16x()
{
     using namespace gms::radiolocation;
     using namespace gms;
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     constexpr std::size_t IQ_nsignals{20ull};
     constexpr std::size_t I_nsamples{10000ull};
     constexpr std::size_t Q_nsamples{10000ull};
     darray_r4_t           I_k_values   = darray_r4_t(IQ_nsignals);     
     darray_r4_t           I_T_values   = darray_r4_t(IQ_nsignals);    
     darray_r4_t           I_Ac_values  = darray_r4_t(IQ_nsignals);    
     darray_r4_t           I_fc_values  = darray_r4_t(IQ_nsignals);    
     darray_r4_t           I_ph0_values = darray_r4_t(IQ_nsignals);     
     darray_r4_t           Q_k_values   = darray_r4_t(IQ_nsignals);     
     darray_r4_t           Q_T_values   = darray_r4_t(IQ_nsignals);    
     darray_r4_t           Q_Ac_values  = darray_r4_t(IQ_nsignals);     
     darray_r4_t           Q_fc_values  = darray_r4_t(IQ_nsignals);     
     darray_r4_t           Q_ph0_values = darray_r4_t(IQ_nsignals);     
     std::vector<I_rectw_bitstream_pdf_params_t>   I_pdf_params = std::vector<I_rectw_bitstream_pdf_params_t>(IQ_nsignals);
     std::vector<Q_rectw_bitstream_pdf_params_t>   Q_pdf_params = std::vector<Q_rectw_bitstream_pdf_params_t>(IQ_nsignals);
     const float           I_k[IQ_nsignals] = {1.0f,2.0f,4.0f,5.0f,10.0f,20.0f,25.0f,40.0f,50.0f,80.0f,100.0f,
                                               125.0f,200.0f,250.0f,400.0f,500.0f,625.0f,1000.0f,2000.0f,2500.0f};  
     const float           I_T[IQ_nsignals] = {10000.0f,5000.0f,2500.0f,2000.0f,1000.0f,500.0f,400.0f,250.0f,200.0f,
                                               125.0f,100.0f,80.0f,50.0f,40.0f,25.0f,20.0f,15.0f,10.0f,5.0f,4.0f};
     const float           Q_k[IQ_nsignals] = {1.0f,2.0f,4.0f,5.0f,10.0f,20.0f,25.0f,40.0f,50.0f,80.0f,100.0f,
                                               125.0f,200.0f,250.0f,400.0f,500.0f,625.0f,1000.0f,2000.0f,2500.0f};  
     const float           Q_T[IQ_nsignals] = {10000.0f,5000.0f,2500.0f,2000.0f,1000.0f,500.0f,400.0f,250.0f,200.0f,
                                               125.0f,100.0f,80.0f,50.0f,40.0f,25.0f,20.0f,15.0f,10.0f,5.0f,4.0f};   
     std::int32_t          I_ch_ret_codes[IQ_nsignals]    = {};
     std::int32_t          Q_ch_ret_codes[IQ_nsignals]    = {};
     std::normal_distribution<float>       gauss_pdf_I_Ac;
     std::uniform_real_distribution<float> uni_pdf_I_fc;
     std::uniform_real_distribution<float> uni_pdf_I_ph0;
     std::normal_distribution<float>       gauss_pdf_Q_Ac;
     std::uniform_real_distribution<float> uni_pdf_Q_fc;
     std::uniform_real_distribution<float> uni_pdf_Q_ph0;
     unsigned long long    seed_I_Ac;
     unsigned long long    seed_I_fc;
     unsigned long long    seed_I_ph0;
     unsigned long long    seed_Q_Ac;
     unsigned long long    seed_Q_fc;
     unsigned long long    seed_Q_ph0;
     const        float                 pdf_I_Ac_incr{0.25f};
     const        float                 pdf_Q_Ac_incr{0.5f};
     float                              pdf_I_Ac_start{0.0f};
     float                              pdf_Q_Ac_start{0.25f};
     seed_I_Ac = __rdtsc();
     seed_Q_Ac = __rdtsc();
     seed_I_fc = __rdtsc();
     seed_I_ph0= __rdtsc();
     seed_Q_fc = __rdtsc();
     seed_Q_ph0= __rdtsc();
     uni_pdf_I_fc = std::uniform_real_distribution<float>(1000.0f,9500.0f);
     auto r_seed_I_fc{std::mt19937(seed_I_fc)};
     uni_pdf_I_ph0= std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,+3.14159265358979323846264338328f);
     auto r_seed_I_ph0{std::mt19937(seed_I_ph0)};
     uni_pdf_Q_fc = std::uniform_real_distribution<float>(4500.0f,20000.0f);
     auto r_seed_Q_fc{std::mt19937(seed_Q_fc)};
     uni_pdf_Q_ph0= std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,+3.14159265358979323846264338328f);
     auto r_seed_Q_ph0{std::mt19937(seed_Q_ph0)};
     for(std::size_t i{0ull}; i != IQ_nsignals; ++i) 
     {
         pdf_I_Ac_start += pdf_I_Ac_incr;         
         I_pdf_params.operator[](i).uni_real_a_r = 0.0f;
         I_pdf_params.operator[](i).uni_real_b_r = 0.0f;
         I_pdf_params.operator[](i).norm_mu_r    = 0.0f;
         I_pdf_params.operator[](i).norm_sigma_r = pdf_I_Ac_start;
         I_pdf_params.operator[](i).cauchy_a_r   = 0.0f;
         I_pdf_params.operator[](i).cauchy_b_r   = 0.0f;
         I_pdf_params.operator[](i).log_norm_m_r = 0.0f;
         I_pdf_params.operator[](i).log_norm_s_r = 0.0f;
         I_pdf_params.operator[](i).expo_gamma_r = 0.0f;
         I_pdf_params.operator[](i).weibull_a_r  = 0.0f;
         I_pdf_params.operator[](i).weibull_b_r  = 0.0f;
         I_pdf_params.operator[](i).gamma_alph_r = 0.0f;
         I_pdf_params.operator[](i).gamma_bet_r  = 0.0f;
         I_pdf_params.operator[](i).poisson_mu   = 0;
         gauss_pdf_I_Ac = std::normal_distribution<float>(0.0f,I_pdf_params.operator[](i).norm_sigma_r);
         auto r_seed_I_Ac{std::mt19937(seed_I_Ac)};
         pdf_Q_Ac_start += pdf_Q_Ac_incr;
         Q_pdf_params.operator[](i).uni_real_a_i = 0.0f;
         Q_pdf_params.operator[](i).uni_real_b_i = 0.0f;
         Q_pdf_params.operator[](i).norm_mu_i    = 0.0f;
         Q_pdf_params.operator[](i).norm_sigma_i = pdf_Q_Ac_start;
         Q_pdf_params.operator[](i).cauchy_a_i   = 0.0f;
         Q_pdf_params.operator[](i).cauchy_b_i   = 0.0f;
         Q_pdf_params.operator[](i).log_norm_m_i = 0.0f;
         Q_pdf_params.operator[](i).log_norm_s_i = 0.0f;
         Q_pdf_params.operator[](i).expo_gamma_i = 0.0f;
         Q_pdf_params.operator[](i).weibull_a_i  = 0.0f;
         Q_pdf_params.operator[](i).weibull_b_i  = 0.0f;
         Q_pdf_params.operator[](i).gamma_alph_i = 0.0f;
         Q_pdf_params.operator[](i).gamma_bet_i  = 0.0f;
         Q_pdf_params.operator[](i).poisson_mu   = 0;
         gauss_pdf_Q_Ac = std::normal_distribution<float>(0.0f,Q_pdf_params.operator[](i).norm_sigma_i);
         auto r_seed_Q_Ac{std::mt19937(seed_Q_Ac)};
         I_k_values.m_data[i]   = I_k[i];
         I_T_values.m_data[i]   = I_T[i];
         Q_k_values.m_data[i]   = Q_k[i];
         Q_T_values.m_data[i]   = Q_T[i];
         I_Ac_values.m_data[i]  = gauss_pdf_I_Ac.operator()(r_seed_I_Ac);
         if(I_Ac_values.m_data[i]==0.0f) {I_Ac_values.m_data[i] = 1.0f;}
         Q_Ac_values.m_data[i]  = gauss_pdf_Q_Ac.operator()(r_seed_Q_Ac);
         if(Q_Ac_values.m_data[i]==0.0f) {Q_Ac_values.m_data[i] = 1.0f;}
         I_fc_values.m_data[i]  = uni_pdf_I_fc.operator()(r_seed_I_fc);
         I_ph0_values.m_data[i] = uni_pdf_I_ph0.operator()(r_seed_I_ph0);
         Q_fc_values.m_data[i]  = uni_pdf_Q_fc.operator()(r_seed_Q_fc);
         Q_ph0_values.m_data[i] = uni_pdf_Q_ph0.operator()(r_seed_Q_ph0);

     }

      int32_t status{};
     
      // For breaking into gdb
      //__asm__ ("int3");
      iq_rectw_bitstream_vsequence_t IQ_rectw_bitstream_vseq = iq_rectw_bitstream_vsequence_t(I_nsamples,IQ_nsignals,I_k_values,I_T_values,
                                                                                              I_Ac_values,I_fc_values,I_ph0_values,
                                                                                              Q_nsamples,Q_k_values,Q_T_values,
                                                                                              Q_Ac_values,Q_fc_values,Q_ph0_values,
                                                                                              I_pdf_params,Q_pdf_params);
/*    
      for(std::size_t i{0ull}; i < IQ_rectw_bitstream_vseq.m_iq_rectw_bitstreams.size(); ++i) 
      {
           std::printf("==========================I-Channel=======================\n");
           std::printf("signal:%llu, m_I_nsamples=%llu\n",i,  IQ_rectw_bitstream_vseq.m_iq_rectw_bitstreams.operator[](i).m_I_nsamples);
           std::printf("signal:%llu, m_I_k=%.7f\n",       i,  IQ_rectw_bitstream_vseq.m_iq_rectw_bitstreams.operator[](i).m_I_k);
           std::printf("signal:%llu, m_I_T=%.7f\n",i,         IQ_rectw_bitstream_vseq.m_iq_rectw_bitstreams.operator[](i).m_I_T);
           std::printf("signal:%llu, m_I_Ac=%.7f\n",i,        IQ_rectw_bitstream_vseq.m_iq_rectw_bitstreams.operator[](i).m_I_Ac);
           std::printf("signal:%llu, m_I_fc=%.7f\n",i,        IQ_rectw_bitstream_vseq.m_iq_rectw_bitstreams.operator[](i).m_I_fc);
           std::printf("signal:%llu, m_I_ph0=%.7f\n", i,      IQ_rectw_bitstream_vseq.m_iq_rectw_bitstreams.operator[](i).m_I_ph0);
           std::printf("==========================Q-Channel=======================\n");
           std::printf("signal:%llu, m_Q_nsamples=%llu\n",i,  IQ_rectw_bitstream_vseq.m_iq_rectw_bitstreams.operator[](i).m_Q_nsamples);
           std::printf("signal:%llu, m_Q_k=%.7f\n",       i,  IQ_rectw_bitstream_vseq.m_iq_rectw_bitstreams.operator[](i).m_Q_k);
           std::printf("signal:%llu, m_Q_T=%.7f\n",i,         IQ_rectw_bitstream_vseq.m_iq_rectw_bitstreams.operator[](i).m_Q_T);
           std::printf("signal:%llu, m_Q_Ac=%.7f\n",i,        IQ_rectw_bitstream_vseq.m_iq_rectw_bitstreams.operator[](i).m_Q_Ac);
           std::printf("signal:%llu, m_Q_fc=%.7f\n",i,        IQ_rectw_bitstream_vseq.m_iq_rectw_bitstreams.operator[](i).m_Q_fc);
           std::printf("signal:%llu, m_Q_ph0=%.7f\n", i,      IQ_rectw_bitstream_vseq.m_iq_rectw_bitstreams.operator[](i).m_Q_ph0);
           
      }
*/      
      char * ctor_name{gms::common::demangle(typeid(IQ_rectw_bitstream_vseq).name(),status)};
      if(status==0 && ctor_name != NULL)
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
      }
      else
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(IQ_rectw_bitstream_vseq).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
      }   
    
      std::int32_t ret_status = IQ_rectw_bitstream_vseq.generate_iq_rectw_bitstream_sequence_avx_u16x<iq_rectw_bitstream_vsequence_t
                                                                                                      ::IQ_rectw_bitstreams_optim_path
                                                                                                      ::AVX_U4X_VEC_PATH,
                                                                                                      iq_rectw_bitstream_vsequence_t
                                                                                                      ::IQ_rectw_bitstream_max_signal_components
                                                                                                      ::MAX_SIGNAL_COMPONENTS_25>(I_ch_ret_codes,
                                                                                                                                  Q_ch_ret_codes);
      if(ret_status<0)
      {
             std::printf("[UNIT-TEST]: -- generate_iq_rectw_bitstream_sequence_avx_u16x -- **FAILED**, status=%d\n",ret_status);
             return;
      }
      else 
      {
             for(std::int32_t i{0};i != 20; ++i) 
             {
                 const std::int32_t I_ch_ret_code_status = I_ch_ret_codes[i];
                 if(I_ch_ret_code_status<0)
                 {
                    std::printf("[FAIL] -- I-channel -- I_ch_ret_codes[%d]=%d\n",i,I_ch_ret_codes[i]);
                 }
                 const std::int32_t Q_ch_ret_code_status = Q_ch_ret_codes[i];
                 if(Q_ch_ret_code_status<0)
                 {
                    std::printf("[FAIL] -- Q-channel -- Q_ch_ret_codes[%d]=%d\n",i,Q_ch_ret_codes[i]);
                 }
             }
      }
      
      std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
    
      iq_rectw_bitstream_vsequence_t::create_sequence_plot(static_cast<std::int32_t>(IQ_rectw_bitstream_vseq.m_I_nsamples),
                                                           IQ_rectw_bitstream_vseq.m_I_vsequence.m_data,nullptr,
                                                           "generate_iq_rectw_bitstream_I_sequence_avx_u16x_test_1_",
                                                           "I-channel vertically summed",false);
      iq_rectw_bitstream_vsequence_t::create_sequence_plot(static_cast<std::int32_t>(IQ_rectw_bitstream_vseq.m_Q_nsamples),
                                                           IQ_rectw_bitstream_vseq.m_Q_vsequence.m_data,nullptr,
                                                           "generate_iq_rectw_bitstream_Q_sequence_avx_u16x_test_1_",
                                                           "Q-channel vertically summed",false); 
                                                          
      printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
      
}

__attribute__((hot))
__attribute__((noinline))
void 
unit_test_iq_rectw_bitstream_vsequence_avx_u10x();

void 
unit_test_iq_rectw_bitstream_vsequence_avx_u10x()
{
     using namespace gms::radiolocation;
     using namespace gms;
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     constexpr std::size_t IQ_nsignals{20ull};
     constexpr std::size_t I_nsamples{10000ull};
     constexpr std::size_t Q_nsamples{10000ull};
     darray_r4_t           I_k_values   = darray_r4_t(IQ_nsignals);     
     darray_r4_t           I_T_values   = darray_r4_t(IQ_nsignals);    
     darray_r4_t           I_Ac_values  = darray_r4_t(IQ_nsignals);    
     darray_r4_t           I_fc_values  = darray_r4_t(IQ_nsignals);    
     darray_r4_t           I_ph0_values = darray_r4_t(IQ_nsignals);     
     darray_r4_t           Q_k_values   = darray_r4_t(IQ_nsignals);     
     darray_r4_t           Q_T_values   = darray_r4_t(IQ_nsignals);    
     darray_r4_t           Q_Ac_values  = darray_r4_t(IQ_nsignals);     
     darray_r4_t           Q_fc_values  = darray_r4_t(IQ_nsignals);     
     darray_r4_t           Q_ph0_values = darray_r4_t(IQ_nsignals);     
     std::vector<I_rectw_bitstream_pdf_params_t>   I_pdf_params = std::vector<I_rectw_bitstream_pdf_params_t>(IQ_nsignals);
     std::vector<Q_rectw_bitstream_pdf_params_t>   Q_pdf_params = std::vector<Q_rectw_bitstream_pdf_params_t>(IQ_nsignals);
     const float           I_k[IQ_nsignals] = {1.0f,2.0f,4.0f,5.0f,10.0f,20.0f,25.0f,40.0f,50.0f,80.0f,100.0f,
                                               125.0f,200.0f,250.0f,400.0f,500.0f,625.0f,1000.0f,2000.0f,2500.0f};  
     const float           I_T[IQ_nsignals] = {10000.0f,5000.0f,2500.0f,2000.0f,1000.0f,500.0f,400.0f,250.0f,200.0f,
                                               125.0f,100.0f,80.0f,50.0f,40.0f,25.0f,20.0f,15.0f,10.0f,5.0f,4.0f};
     const float           Q_k[IQ_nsignals] = {1.0f,2.0f,4.0f,5.0f,10.0f,20.0f,25.0f,40.0f,50.0f,80.0f,100.0f,
                                               125.0f,200.0f,250.0f,400.0f,500.0f,625.0f,1000.0f,2000.0f,2500.0f};  
     const float           Q_T[IQ_nsignals] = {10000.0f,5000.0f,2500.0f,2000.0f,1000.0f,500.0f,400.0f,250.0f,200.0f,
                                               125.0f,100.0f,80.0f,50.0f,40.0f,25.0f,20.0f,15.0f,10.0f,5.0f,4.0f};   
     std::int32_t          I_ch_ret_codes[IQ_nsignals]    = {};
     std::int32_t          Q_ch_ret_codes[IQ_nsignals]    = {};
     std::normal_distribution<float>       gauss_pdf_I_Ac;
     std::uniform_real_distribution<float> uni_pdf_I_fc;
     std::uniform_real_distribution<float> uni_pdf_I_ph0;
     std::normal_distribution<float>       gauss_pdf_Q_Ac;
     std::uniform_real_distribution<float> uni_pdf_Q_fc;
     std::uniform_real_distribution<float> uni_pdf_Q_ph0;
     unsigned long long    seed_I_Ac;
     unsigned long long    seed_I_fc;
     unsigned long long    seed_I_ph0;
     unsigned long long    seed_Q_Ac;
     unsigned long long    seed_Q_fc;
     unsigned long long    seed_Q_ph0;
     const        float                 pdf_I_Ac_incr{0.25f};
     const        float                 pdf_Q_Ac_incr{0.5f};
     float                              pdf_I_Ac_start{0.0f};
     float                              pdf_Q_Ac_start{0.25f};
     seed_I_Ac = __rdtsc();
     seed_Q_Ac = __rdtsc();
     seed_I_fc = __rdtsc();
     seed_I_ph0= __rdtsc();
     seed_Q_fc = __rdtsc();
     seed_Q_ph0= __rdtsc();
     uni_pdf_I_fc = std::uniform_real_distribution<float>(1000.0f,9500.0f);
     auto r_seed_I_fc{std::mt19937(seed_I_fc)};
     uni_pdf_I_ph0= std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,+3.14159265358979323846264338328f);
     auto r_seed_I_ph0{std::mt19937(seed_I_ph0)};
     uni_pdf_Q_fc = std::uniform_real_distribution<float>(4500.0f,20000.0f);
     auto r_seed_Q_fc{std::mt19937(seed_Q_fc)};
     uni_pdf_Q_ph0= std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,+3.14159265358979323846264338328f);
     auto r_seed_Q_ph0{std::mt19937(seed_Q_ph0)};
     for(std::size_t i{0ull}; i != IQ_nsignals; ++i) 
     {
         pdf_I_Ac_start += pdf_I_Ac_incr;         
         I_pdf_params.operator[](i).uni_real_a_r = 0.0f;
         I_pdf_params.operator[](i).uni_real_b_r = 0.0f;
         I_pdf_params.operator[](i).norm_mu_r    = 0.0f;
         I_pdf_params.operator[](i).norm_sigma_r = pdf_I_Ac_start;
         I_pdf_params.operator[](i).cauchy_a_r   = 0.0f;
         I_pdf_params.operator[](i).cauchy_b_r   = 0.0f;
         I_pdf_params.operator[](i).log_norm_m_r = 0.0f;
         I_pdf_params.operator[](i).log_norm_s_r = 0.0f;
         I_pdf_params.operator[](i).expo_gamma_r = 0.0f;
         I_pdf_params.operator[](i).weibull_a_r  = 0.0f;
         I_pdf_params.operator[](i).weibull_b_r  = 0.0f;
         I_pdf_params.operator[](i).gamma_alph_r = 0.0f;
         I_pdf_params.operator[](i).gamma_bet_r  = 0.0f;
         I_pdf_params.operator[](i).poisson_mu   = 0;
         gauss_pdf_I_Ac = std::normal_distribution<float>(0.0f,I_pdf_params.operator[](i).norm_sigma_r);
         auto r_seed_I_Ac{std::mt19937(seed_I_Ac)};
         pdf_Q_Ac_start += pdf_Q_Ac_incr;
         Q_pdf_params.operator[](i).uni_real_a_i = 0.0f;
         Q_pdf_params.operator[](i).uni_real_b_i = 0.0f;
         Q_pdf_params.operator[](i).norm_mu_i    = 0.0f;
         Q_pdf_params.operator[](i).norm_sigma_i = pdf_Q_Ac_start;
         Q_pdf_params.operator[](i).cauchy_a_i   = 0.0f;
         Q_pdf_params.operator[](i).cauchy_b_i   = 0.0f;
         Q_pdf_params.operator[](i).log_norm_m_i = 0.0f;
         Q_pdf_params.operator[](i).log_norm_s_i = 0.0f;
         Q_pdf_params.operator[](i).expo_gamma_i = 0.0f;
         Q_pdf_params.operator[](i).weibull_a_i  = 0.0f;
         Q_pdf_params.operator[](i).weibull_b_i  = 0.0f;
         Q_pdf_params.operator[](i).gamma_alph_i = 0.0f;
         Q_pdf_params.operator[](i).gamma_bet_i  = 0.0f;
         Q_pdf_params.operator[](i).poisson_mu   = 0;
         gauss_pdf_Q_Ac = std::normal_distribution<float>(0.0f,Q_pdf_params.operator[](i).norm_sigma_i);
         auto r_seed_Q_Ac{std::mt19937(seed_Q_Ac)};
         I_k_values.m_data[i]   = I_k[i];
         I_T_values.m_data[i]   = I_T[i];
         Q_k_values.m_data[i]   = Q_k[i];
         Q_T_values.m_data[i]   = Q_T[i];
         I_Ac_values.m_data[i]  = gauss_pdf_I_Ac.operator()(r_seed_I_Ac);
         if(I_Ac_values.m_data[i]==0.0f) {I_Ac_values.m_data[i] = 1.0f;}
         Q_Ac_values.m_data[i]  = gauss_pdf_Q_Ac.operator()(r_seed_Q_Ac);
         if(Q_Ac_values.m_data[i]==0.0f) {Q_Ac_values.m_data[i] = 1.0f;}
         I_fc_values.m_data[i]  = uni_pdf_I_fc.operator()(r_seed_I_fc);
         I_ph0_values.m_data[i] = uni_pdf_I_ph0.operator()(r_seed_I_ph0);
         Q_fc_values.m_data[i]  = uni_pdf_Q_fc.operator()(r_seed_Q_fc);
         Q_ph0_values.m_data[i] = uni_pdf_Q_ph0.operator()(r_seed_Q_ph0);

     }

      int32_t status{};
     
      // For breaking into gdb
      //__asm__ ("int3");
      iq_rectw_bitstream_vsequence_t IQ_rectw_bitstream_vseq = iq_rectw_bitstream_vsequence_t(I_nsamples,IQ_nsignals,I_k_values,I_T_values,
                                                                                              I_Ac_values,I_fc_values,I_ph0_values,
                                                                                              Q_nsamples,Q_k_values,Q_T_values,
                                                                                              Q_Ac_values,Q_fc_values,Q_ph0_values,
                                                                                              I_pdf_params,Q_pdf_params);

      char * ctor_name{gms::common::demangle(typeid(IQ_rectw_bitstream_vseq).name(),status)};
      if(status==0 && ctor_name != NULL)
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
      }
      else
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(IQ_rectw_bitstream_vseq).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
      }   
    
      std::int32_t ret_status = IQ_rectw_bitstream_vseq.generate_iq_rectw_bitstream_sequence_avx_u10x<iq_rectw_bitstream_vsequence_t
                                                                                                      ::IQ_rectw_bitstreams_optim_path
                                                                                                      ::AVX_U4X_VEC_PATH,
                                                                                                      iq_rectw_bitstream_vsequence_t
                                                                                                      ::IQ_rectw_bitstream_max_signal_components
                                                                                                      ::MAX_SIGNAL_COMPONENTS_25>(I_ch_ret_codes,
                                                                                                                                  Q_ch_ret_codes);
      if(ret_status<0)
      {
             std::printf("[UNIT-TEST]: -- generate_iq_rectw_bitstream_sequence_avx_u10x -- **FAILED**, status=%d\n",ret_status);
             return;
      }
      else 
      {
             for(std::int32_t i{0};i != 20; ++i) 
             {
                 const std::int32_t I_ch_ret_code_status = I_ch_ret_codes[i];
                 if(I_ch_ret_code_status<0)
                 {
                    std::printf("[FAIL] -- I-channel -- I_ch_ret_codes[%d]=%d\n",i,I_ch_ret_codes[i]);
                 }
                 const std::int32_t Q_ch_ret_code_status = Q_ch_ret_codes[i];
                 if(Q_ch_ret_code_status<0)
                 {
                    std::printf("[FAIL] -- Q-channel -- Q_ch_ret_codes[%d]=%d\n",i,Q_ch_ret_codes[i]);
                 }
             }
      }
      
      std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
    
      iq_rectw_bitstream_vsequence_t::create_sequence_plot(static_cast<std::int32_t>(IQ_rectw_bitstream_vseq.m_I_nsamples),
                                                           IQ_rectw_bitstream_vseq.m_I_vsequence.m_data,nullptr,
                                                           "generate_iq_rectw_bitstream_I_sequence_avx_u10x_test_1_",
                                                           "I-channel vertically summed",false);
      iq_rectw_bitstream_vsequence_t::create_sequence_plot(static_cast<std::int32_t>(IQ_rectw_bitstream_vseq.m_Q_nsamples),
                                                           IQ_rectw_bitstream_vseq.m_Q_vsequence.m_data,nullptr,
                                                           "generate_iq_rectw_bitstream_Q_sequence_avx_u10x_test_1_",
                                                           "Q-channel vertically summed",false); 
                                                          
      printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
      
}

__attribute__((hot))
__attribute__((noinline))
void 
unit_test_iq_rectw_bitstream_vsequence_avx_u6x();

void 
unit_test_iq_rectw_bitstream_vsequence_avx_u6x()
{
     using namespace gms::radiolocation;
     using namespace gms;
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     constexpr std::size_t IQ_nsignals{20ull};
     constexpr std::size_t I_nsamples{10000ull};
     constexpr std::size_t Q_nsamples{10000ull};
     darray_r4_t           I_k_values   = darray_r4_t(IQ_nsignals);     
     darray_r4_t           I_T_values   = darray_r4_t(IQ_nsignals);    
     darray_r4_t           I_Ac_values  = darray_r4_t(IQ_nsignals);    
     darray_r4_t           I_fc_values  = darray_r4_t(IQ_nsignals);    
     darray_r4_t           I_ph0_values = darray_r4_t(IQ_nsignals);     
     darray_r4_t           Q_k_values   = darray_r4_t(IQ_nsignals);     
     darray_r4_t           Q_T_values   = darray_r4_t(IQ_nsignals);    
     darray_r4_t           Q_Ac_values  = darray_r4_t(IQ_nsignals);     
     darray_r4_t           Q_fc_values  = darray_r4_t(IQ_nsignals);     
     darray_r4_t           Q_ph0_values = darray_r4_t(IQ_nsignals);     
     std::vector<I_rectw_bitstream_pdf_params_t>   I_pdf_params = std::vector<I_rectw_bitstream_pdf_params_t>(IQ_nsignals);
     std::vector<Q_rectw_bitstream_pdf_params_t>   Q_pdf_params = std::vector<Q_rectw_bitstream_pdf_params_t>(IQ_nsignals);
     const float           I_k[IQ_nsignals] = {1.0f,2.0f,4.0f,5.0f,10.0f,20.0f,25.0f,40.0f,50.0f,80.0f,100.0f,
                                               125.0f,200.0f,250.0f,400.0f,500.0f,625.0f,1000.0f,2000.0f,2500.0f};  
     const float           I_T[IQ_nsignals] = {10000.0f,5000.0f,2500.0f,2000.0f,1000.0f,500.0f,400.0f,250.0f,200.0f,
                                               125.0f,100.0f,80.0f,50.0f,40.0f,25.0f,20.0f,15.0f,10.0f,5.0f,4.0f};
     const float           Q_k[IQ_nsignals] = {1.0f,2.0f,4.0f,5.0f,10.0f,20.0f,25.0f,40.0f,50.0f,80.0f,100.0f,
                                               125.0f,200.0f,250.0f,400.0f,500.0f,625.0f,1000.0f,2000.0f,2500.0f};  
     const float           Q_T[IQ_nsignals] = {10000.0f,5000.0f,2500.0f,2000.0f,1000.0f,500.0f,400.0f,250.0f,200.0f,
                                               125.0f,100.0f,80.0f,50.0f,40.0f,25.0f,20.0f,15.0f,10.0f,5.0f,4.0f};   
     std::int32_t          I_ch_ret_codes[IQ_nsignals]    = {};
     std::int32_t          Q_ch_ret_codes[IQ_nsignals]    = {};
     std::normal_distribution<float>       gauss_pdf_I_Ac;
     std::uniform_real_distribution<float> uni_pdf_I_fc;
     std::uniform_real_distribution<float> uni_pdf_I_ph0;
     std::normal_distribution<float>       gauss_pdf_Q_Ac;
     std::uniform_real_distribution<float> uni_pdf_Q_fc;
     std::uniform_real_distribution<float> uni_pdf_Q_ph0;
     unsigned long long    seed_I_Ac;
     unsigned long long    seed_I_fc;
     unsigned long long    seed_I_ph0;
     unsigned long long    seed_Q_Ac;
     unsigned long long    seed_Q_fc;
     unsigned long long    seed_Q_ph0;
     const        float                 pdf_I_Ac_incr{0.25f};
     const        float                 pdf_Q_Ac_incr{0.5f};
     float                              pdf_I_Ac_start{0.0f};
     float                              pdf_Q_Ac_start{0.25f};
     seed_I_Ac = __rdtsc();
     seed_Q_Ac = __rdtsc();
     seed_I_fc = __rdtsc();
     seed_I_ph0= __rdtsc();
     seed_Q_fc = __rdtsc();
     seed_Q_ph0= __rdtsc();
     uni_pdf_I_fc = std::uniform_real_distribution<float>(1000.0f,9500.0f);
     auto r_seed_I_fc{std::mt19937(seed_I_fc)};
     uni_pdf_I_ph0= std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,+3.14159265358979323846264338328f);
     auto r_seed_I_ph0{std::mt19937(seed_I_ph0)};
     uni_pdf_Q_fc = std::uniform_real_distribution<float>(4500.0f,20000.0f);
     auto r_seed_Q_fc{std::mt19937(seed_Q_fc)};
     uni_pdf_Q_ph0= std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,+3.14159265358979323846264338328f);
     auto r_seed_Q_ph0{std::mt19937(seed_Q_ph0)};
     for(std::size_t i{0ull}; i != IQ_nsignals; ++i) 
     {
         pdf_I_Ac_start += pdf_I_Ac_incr;         
         I_pdf_params.operator[](i).uni_real_a_r = 0.0f;
         I_pdf_params.operator[](i).uni_real_b_r = 0.0f;
         I_pdf_params.operator[](i).norm_mu_r    = 0.0f;
         I_pdf_params.operator[](i).norm_sigma_r = pdf_I_Ac_start;
         I_pdf_params.operator[](i).cauchy_a_r   = 0.0f;
         I_pdf_params.operator[](i).cauchy_b_r   = 0.0f;
         I_pdf_params.operator[](i).log_norm_m_r = 0.0f;
         I_pdf_params.operator[](i).log_norm_s_r = 0.0f;
         I_pdf_params.operator[](i).expo_gamma_r = 0.0f;
         I_pdf_params.operator[](i).weibull_a_r  = 0.0f;
         I_pdf_params.operator[](i).weibull_b_r  = 0.0f;
         I_pdf_params.operator[](i).gamma_alph_r = 0.0f;
         I_pdf_params.operator[](i).gamma_bet_r  = 0.0f;
         I_pdf_params.operator[](i).poisson_mu   = 0;
         gauss_pdf_I_Ac = std::normal_distribution<float>(0.0f,I_pdf_params.operator[](i).norm_sigma_r);
         auto r_seed_I_Ac{std::mt19937(seed_I_Ac)};
         pdf_Q_Ac_start += pdf_Q_Ac_incr;
         Q_pdf_params.operator[](i).uni_real_a_i = 0.0f;
         Q_pdf_params.operator[](i).uni_real_b_i = 0.0f;
         Q_pdf_params.operator[](i).norm_mu_i    = 0.0f;
         Q_pdf_params.operator[](i).norm_sigma_i = pdf_Q_Ac_start;
         Q_pdf_params.operator[](i).cauchy_a_i   = 0.0f;
         Q_pdf_params.operator[](i).cauchy_b_i   = 0.0f;
         Q_pdf_params.operator[](i).log_norm_m_i = 0.0f;
         Q_pdf_params.operator[](i).log_norm_s_i = 0.0f;
         Q_pdf_params.operator[](i).expo_gamma_i = 0.0f;
         Q_pdf_params.operator[](i).weibull_a_i  = 0.0f;
         Q_pdf_params.operator[](i).weibull_b_i  = 0.0f;
         Q_pdf_params.operator[](i).gamma_alph_i = 0.0f;
         Q_pdf_params.operator[](i).gamma_bet_i  = 0.0f;
         Q_pdf_params.operator[](i).poisson_mu   = 0;
         gauss_pdf_Q_Ac = std::normal_distribution<float>(0.0f,Q_pdf_params.operator[](i).norm_sigma_i);
         auto r_seed_Q_Ac{std::mt19937(seed_Q_Ac)};
         I_k_values.m_data[i]   = I_k[i];
         I_T_values.m_data[i]   = I_T[i];
         Q_k_values.m_data[i]   = Q_k[i];
         Q_T_values.m_data[i]   = Q_T[i];
         I_Ac_values.m_data[i]  = gauss_pdf_I_Ac.operator()(r_seed_I_Ac);
         if(I_Ac_values.m_data[i]==0.0f) {I_Ac_values.m_data[i] = 1.0f;}
         Q_Ac_values.m_data[i]  = gauss_pdf_Q_Ac.operator()(r_seed_Q_Ac);
         if(Q_Ac_values.m_data[i]==0.0f) {Q_Ac_values.m_data[i] = 1.0f;}
         I_fc_values.m_data[i]  = uni_pdf_I_fc.operator()(r_seed_I_fc);
         I_ph0_values.m_data[i] = uni_pdf_I_ph0.operator()(r_seed_I_ph0);
         Q_fc_values.m_data[i]  = uni_pdf_Q_fc.operator()(r_seed_Q_fc);
         Q_ph0_values.m_data[i] = uni_pdf_Q_ph0.operator()(r_seed_Q_ph0);

     }

      int32_t status{};
     
      // For breaking into gdb
      //__asm__ ("int3");
      iq_rectw_bitstream_vsequence_t IQ_rectw_bitstream_vseq = iq_rectw_bitstream_vsequence_t(I_nsamples,IQ_nsignals,I_k_values,I_T_values,
                                                                                              I_Ac_values,I_fc_values,I_ph0_values,
                                                                                              Q_nsamples,Q_k_values,Q_T_values,
                                                                                              Q_Ac_values,Q_fc_values,Q_ph0_values,
                                                                                              I_pdf_params,Q_pdf_params);

      char * ctor_name{gms::common::demangle(typeid(IQ_rectw_bitstream_vseq).name(),status)};
      if(status==0 && ctor_name != NULL)
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
      }
      else
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(IQ_rectw_bitstream_vseq).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
      }   
    
      std::int32_t ret_status = IQ_rectw_bitstream_vseq.generate_iq_rectw_bitstream_sequence_avx_u6x<iq_rectw_bitstream_vsequence_t
                                                                                                      ::IQ_rectw_bitstreams_optim_path
                                                                                                      ::AVX_U4X_VEC_PATH,
                                                                                                      iq_rectw_bitstream_vsequence_t
                                                                                                      ::IQ_rectw_bitstream_max_signal_components
                                                                                                      ::MAX_SIGNAL_COMPONENTS_25>(I_ch_ret_codes,
                                                                                                                                  Q_ch_ret_codes);
      if(ret_status<0)
      {
             std::printf("[UNIT-TEST]: -- generate_iq_rectw_bitstream_sequence_avx_u6x -- **FAILED**, status=%d\n",ret_status);
             return;
      }
      else 
      {
             for(std::int32_t i{0};i != 20; ++i) 
             {
                 const std::int32_t I_ch_ret_code_status = I_ch_ret_codes[i];
                 if(I_ch_ret_code_status<0)
                 {
                    std::printf("[FAIL] -- I-channel -- I_ch_ret_codes[%d]=%d\n",i,I_ch_ret_codes[i]);
                 }
                 const std::int32_t Q_ch_ret_code_status = Q_ch_ret_codes[i];
                 if(Q_ch_ret_code_status<0)
                 {
                    std::printf("[FAIL] -- Q-channel -- Q_ch_ret_codes[%d]=%d\n",i,Q_ch_ret_codes[i]);
                 }
             }
      }
      
      std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
    
      iq_rectw_bitstream_vsequence_t::create_sequence_plot(static_cast<std::int32_t>(IQ_rectw_bitstream_vseq.m_I_nsamples),
                                                           IQ_rectw_bitstream_vseq.m_I_vsequence.m_data,nullptr,
                                                           "generate_iq_rectw_bitstream_I_sequence_avx_u6x_test_1_",
                                                           "I-channel vertically summed",false);
      iq_rectw_bitstream_vsequence_t::create_sequence_plot(static_cast<std::int32_t>(IQ_rectw_bitstream_vseq.m_Q_nsamples),
                                                           IQ_rectw_bitstream_vseq.m_Q_vsequence.m_data,nullptr,
                                                           "generate_iq_rectw_bitstream_Q_sequence_avx_u6x_test_1_",
                                                           "Q-channel vertically summed",false); 
                                                          
      printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
      
}



int main()
{
    unit_test_iq_rectw_bitstream_vsequence_avx_u16x();
    unit_test_iq_rectw_bitstream_vsequence_avx_u10x();
    unit_test_iq_rectw_bitstream_vsequence_avx_u6x();
    return 0;
}