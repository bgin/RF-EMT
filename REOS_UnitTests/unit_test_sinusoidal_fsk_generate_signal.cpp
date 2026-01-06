#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_dyn_array.h"
#include "GMS_sinusoidal_fsk.h"

/*
   icpc -o unit_test_sinusoidal_fsk_generate_signal -O3 -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_sse_memcpy.h GMS_sse_memcpy.cpp GMS_avx_memcpy.h GMS_avx_memcpy.cpp GMS_avx512_memcpy.h GMS_avx512_memcpu.cpp GMS_cephes_sin_cos.h GMS_loop_indices_LUT.h GMS_loop_indices_LUT.cpp GMS_iq_rectwave_bitstream.h GMS_iq_rectwave_bitstream.cpp GMS_iq_rectw_bitstream_vsequence.h GMS_iq_rectw_bitstream_vsequence.cpp GMS_sinusoidal_fsk.h GMS_sinusoidal_fsk.cpp unit_test_sinusoidal_fsk_generate_signal.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_sse_memcpy.h GMS_sse_memcpy.cpp GMS_avx_memcpy.h GMS_avx_memcpy.cpp GMS_avx512_memcpy.h GMS_avx512_memcpu.cpp GMS_cephes_sin_cos.h GMS_loop_indices_LUT.h GMS_loop_indices_LUT.cpp GMS_iq_rectwave_bitstream.h GMS_iq_rectwave_bitstream.cpp GMS_iq_rectw_bitstream_vsequence.h GMS_iq_rectw_bitstream_vsequence.cpp GMS_sinusoidal_fsk.h GMS_sinusoidal_fsk.cpp unit_test_sinusoidal_fsk_generate_signal.cpp

*/

__attribute__((hot))
__attribute__((noinline))
void 
unit_test_generate_fsk_signal_u8x();

void  
unit_test_generate_fsk_signal_u8x()
{
     using namespace gms::radiolocation;
     using namespace gms;
     constexpr float k{50.0f};
     constexpr float T{200.0f};
     constexpr std::size_t sfsk_nsamples{static_cast<std::size_t>(k*T)}; // 10 periods of 20 samples each
     constexpr std::size_t I_nsamples{sfsk_nsamples};
     constexpr std::size_t Q_nsamples{sfsk_nsamples};
     constexpr std::size_t psfunc_nsamples{sfsk_nsamples};      
     constexpr float Ac{1.0f};
     constexpr float fc{9000.0f};
     constexpr float I_fc{9000.0f};
     constexpr float Q_fc{I_fc};
     constexpr float I_ph0{0.78539816339744830961566084582f};
     constexpr float Q_ph0{0.78539816339744830961566084582f};
     constexpr float ph0{0.0f};
     constexpr std::size_t IQ_nsignals{20ull};
     constexpr std::size_t iq_I_nsamples{10000ull};
     constexpr std::size_t iq_Q_nsamples{10000ull};
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

      int32_t status_1{};
     
      // For breaking into gdb
      //__asm__ ("int3");
      iq_rectw_bitstream_vsequence_t IQ_rectw_bitstream_vseq = iq_rectw_bitstream_vsequence_t(iq_I_nsamples,IQ_nsignals,I_k_values,I_T_values,
                                                                                              I_Ac_values,I_fc_values,I_ph0_values,
                                                                                              iq_Q_nsamples,Q_k_values,Q_T_values,
                                                                                              Q_Ac_values,Q_fc_values,Q_ph0_values,
                                                                                              I_pdf_params,Q_pdf_params);

      char * ctor_name_1{gms::common::demangle(typeid(IQ_rectw_bitstream_vseq).name(),status_1)};
      if(status_1==0 && ctor_name_1 != NULL)
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name_1);
          gms::common::gms_mm_free(ctor_name_1);
      }
      else
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(IQ_rectw_bitstream_vseq).name());
          if(ctor_name_1 != NULL) gms::common::gms_mm_free(ctor_name_1);
      }   
    
      std::int32_t ret_status_1 = IQ_rectw_bitstream_vseq.generate_iq_rectw_bitstream_sequence_avx512_u16x<iq_rectw_bitstream_vsequence_t
                                                                                                      ::IQ_rectw_bitstreams_optim_path
                                                                                                      ::AVX512_U4X_VEC_PATH,
                                                                                                      iq_rectw_bitstream_vsequence_t
                                                                                                      ::IQ_rectw_bitstream_max_signal_components
                                                                                                      ::MAX_SIGNAL_COMPONENTS_25>(I_ch_ret_codes,
                                                                                                                                  Q_ch_ret_codes);
      if(ret_status_1<0)
      {
             std::printf("[UNIT-TEST]: -- generate_iq_rectw_bitstream_sequence_avx512_u16x -- **FAILED**, status=%d\n",ret_status_1);
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

      int32_t status_2{};
      // For breaking into gdb
      //__asm__ ("int3");
      sinusoidal_fsk_t SFSK_signal(sfsk_nsamples,I_nsamples,Q_nsamples,
                                   psfunc_nsamples,k,T,Ac,fc,ph0,
                                   I_fc,Q_fc,I_ph0,Q_ph0);
      char * ctor_name_2{gms::common::demangle(typeid(SFSK_signal).name(),status_2)};
      if(status_2==0 && ctor_name_2 != NULL)
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name_2);
          gms::common::gms_mm_free(ctor_name_2);
      }
      else
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(SFSK_signal).name());
          if(ctor_name_2 != NULL) gms::common::gms_mm_free(ctor_name_2);
      }   

      
     std::printf("[UNIT-TEST]: -- Calling: generate_fsk_signal_u8x\n");
     const std::int32_t ret_status_2 = 
     SFSK_signal.generate_fsk_signal_u8x<sinusoidal_fsk_t::I_channel_bitstream_optimized_path::vector_avx512_u4x_path,
                                                     sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::vector_avx512_u4x_path, 
                                                     sinusoidal_fsk_t::pulse_shaping_function_optimized_path::vector_avx512_path>(IQ_rectw_bitstream_vseq,1,1,true,true,true,true);
          if(ret_status_2<0)
          {
             std::printf("[UNIT-TEST]: -- generate_fsk_signal_u8x -- **FAILED**, status=%d\n",ret_status_2);
             return;
          }
          std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
          sinusoidal_fsk_t::create_signal_plot(static_cast<std::uint32_t>(SFSK_signal.m_sfsk_nsamples),SFSK_signal.m_sfsk_signal.m_data,
                                                                          nullptr,"generate_fsk_signal_bitstream_vseq_u8x_test_1_",
                                                                          "Sinusoidal-FSK Signal",false);   
         printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
     
}

__attribute__((hot))
__attribute__((noinline))
void 
unit_test_generate_fsk_signal_sse_u4x();

void  
unit_test_generate_fsk_signal_sse_u4x()
{
     using namespace gms::radiolocation;
     using namespace gms;
     constexpr float k{50.0f};
     constexpr float T{200.0f};
     constexpr std::size_t sfsk_nsamples{static_cast<std::size_t>(k*T)}; // 10 periods of 20 samples each
     constexpr std::size_t I_nsamples{sfsk_nsamples};
     constexpr std::size_t Q_nsamples{sfsk_nsamples};
     constexpr std::size_t psfunc_nsamples{sfsk_nsamples};      
     constexpr float Ac{1.0f};
     constexpr float fc{9000.0f};
     constexpr float I_fc{9000.0f};
     constexpr float Q_fc{I_fc};
     constexpr float I_ph0{0.78539816339744830961566084582f};
     constexpr float Q_ph0{0.78539816339744830961566084582f};
     constexpr float ph0{0.0f};
     constexpr std::size_t IQ_nsignals{20ull};
     constexpr std::size_t iq_I_nsamples{10000ull};
     constexpr std::size_t iq_Q_nsamples{10000ull};
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

      int32_t status_1{};
     
      // For breaking into gdb
      //__asm__ ("int3");
      iq_rectw_bitstream_vsequence_t IQ_rectw_bitstream_vseq = iq_rectw_bitstream_vsequence_t(iq_I_nsamples,IQ_nsignals,I_k_values,I_T_values,
                                                                                              I_Ac_values,I_fc_values,I_ph0_values,
                                                                                              iq_Q_nsamples,Q_k_values,Q_T_values,
                                                                                              Q_Ac_values,Q_fc_values,Q_ph0_values,
                                                                                              I_pdf_params,Q_pdf_params);

      char * ctor_name_1{gms::common::demangle(typeid(IQ_rectw_bitstream_vseq).name(),status_1)};
      if(status_1==0 && ctor_name_1 != NULL)
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name_1);
          gms::common::gms_mm_free(ctor_name_1);
      }
      else
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(IQ_rectw_bitstream_vseq).name());
          if(ctor_name_1 != NULL) gms::common::gms_mm_free(ctor_name_1);
      }   
    
      std::int32_t ret_status_1 = IQ_rectw_bitstream_vseq.generate_iq_rectw_bitstream_sequence_avx512_u16x<iq_rectw_bitstream_vsequence_t
                                                                                                      ::IQ_rectw_bitstreams_optim_path
                                                                                                      ::AVX512_U4X_VEC_PATH,
                                                                                                      iq_rectw_bitstream_vsequence_t
                                                                                                      ::IQ_rectw_bitstream_max_signal_components
                                                                                                      ::MAX_SIGNAL_COMPONENTS_25>(I_ch_ret_codes,
                                                                                                                                  Q_ch_ret_codes);
      if(ret_status_1<0)
      {
             std::printf("[UNIT-TEST]: -- generate_iq_rectw_bitstream_sequence_avx512_u16x -- **FAILED**, status=%d\n",ret_status_1);
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

      int32_t status_2{};
      // For breaking into gdb
      //__asm__ ("int3");
      sinusoidal_fsk_t SFSK_signal(sfsk_nsamples,I_nsamples,Q_nsamples,
                                   psfunc_nsamples,k,T,Ac,fc,ph0,
                                   I_fc,Q_fc,I_ph0,Q_ph0);
      char * ctor_name_2{gms::common::demangle(typeid(SFSK_signal).name(),status_2)};
      if(status_2==0 && ctor_name_2 != NULL)
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name_2);
          gms::common::gms_mm_free(ctor_name_2);
      }
      else
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(SFSK_signal).name());
          if(ctor_name_2 != NULL) gms::common::gms_mm_free(ctor_name_2);
      }   

      
     std::printf("[UNIT-TEST]: -- Calling: generate_fsk_signal_sse_u4x\n");
     const std::int32_t ret_status_2 = 
     SFSK_signal.generate_fsk_signal_sse_u4x<sinusoidal_fsk_t::I_channel_bitstream_optimized_path::vector_avx512_u4x_path,
                                                     sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::vector_avx512_u4x_path, 
                                                     sinusoidal_fsk_t::pulse_shaping_function_optimized_path::vector_avx512_path>(IQ_rectw_bitstream_vseq,1,1,true,true,true,true,true);
          if(ret_status_2<0)
          {
             std::printf("[UNIT-TEST]: -- generate_fsk_signal_sse_u4x -- **FAILED**, status=%d\n",ret_status_2);
             return;
          }
          std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
          sinusoidal_fsk_t::create_signal_plot(static_cast<std::uint32_t>(SFSK_signal.m_sfsk_nsamples),SFSK_signal.m_sfsk_signal.m_data,
                                                                          nullptr,"generate_fsk_signal_bitstream_vseq_sse_u4x_test_1_",
                                                                          "Sinusoidal-FSK Signal",false);   
         printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
     
}

__attribute__((hot))
__attribute__((noinline))
void 
unit_test_generate_fsk_signal_avx_u4x();

void  
unit_test_generate_fsk_signal_avx_u4x()
{
     using namespace gms::radiolocation;
     using namespace gms;
     constexpr float k{50.0f};
     constexpr float T{200.0f};
     constexpr std::size_t sfsk_nsamples{static_cast<std::size_t>(k*T)}; 
     constexpr std::size_t I_nsamples{sfsk_nsamples};
     constexpr std::size_t Q_nsamples{sfsk_nsamples};
     constexpr std::size_t psfunc_nsamples{sfsk_nsamples};      
     constexpr float Ac{1.0f};
     constexpr float fc{9000.0f};
     constexpr float I_fc{9000.0f};
     constexpr float Q_fc{I_fc};
     constexpr float I_ph0{0.78539816339744830961566084582f};
     constexpr float Q_ph0{0.78539816339744830961566084582f};
     constexpr float ph0{0.0f};
     constexpr std::size_t IQ_nsignals{20ull};
     constexpr std::size_t iq_I_nsamples{10000ull};
     constexpr std::size_t iq_Q_nsamples{10000ull};
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

      int32_t status_1{};
     
      // For breaking into gdb
      //__asm__ ("int3");
      iq_rectw_bitstream_vsequence_t IQ_rectw_bitstream_vseq = iq_rectw_bitstream_vsequence_t(iq_I_nsamples,IQ_nsignals,I_k_values,I_T_values,
                                                                                              I_Ac_values,I_fc_values,I_ph0_values,
                                                                                              iq_Q_nsamples,Q_k_values,Q_T_values,
                                                                                              Q_Ac_values,Q_fc_values,Q_ph0_values,
                                                                                              I_pdf_params,Q_pdf_params);

      char * ctor_name_1{gms::common::demangle(typeid(IQ_rectw_bitstream_vseq).name(),status_1)};
      if(status_1==0 && ctor_name_1 != NULL)
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name_1);
          gms::common::gms_mm_free(ctor_name_1);
      }
      else
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(IQ_rectw_bitstream_vseq).name());
          if(ctor_name_1 != NULL) gms::common::gms_mm_free(ctor_name_1);
      }   
    
      std::int32_t ret_status_1 = IQ_rectw_bitstream_vseq.generate_iq_rectw_bitstream_sequence_avx512_u16x<iq_rectw_bitstream_vsequence_t
                                                                                                      ::IQ_rectw_bitstreams_optim_path
                                                                                                      ::AVX512_U4X_VEC_PATH,
                                                                                                      iq_rectw_bitstream_vsequence_t
                                                                                                      ::IQ_rectw_bitstream_max_signal_components
                                                                                                      ::MAX_SIGNAL_COMPONENTS_25>(I_ch_ret_codes,
                                                                                                                                  Q_ch_ret_codes);
      if(ret_status_1<0)
      {
             std::printf("[UNIT-TEST]: -- generate_iq_rectw_bitstream_sequence_avx512_u16x -- **FAILED**, status=%d\n",ret_status_1);
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

      int32_t status_2{};
      // For breaking into gdb
      //__asm__ ("int3");
      sinusoidal_fsk_t SFSK_signal(sfsk_nsamples,I_nsamples,Q_nsamples,
                                   psfunc_nsamples,k,T,Ac,fc,ph0,
                                   I_fc,Q_fc,I_ph0,Q_ph0);
      char * ctor_name_2{gms::common::demangle(typeid(SFSK_signal).name(),status_2)};
      if(status_2==0 && ctor_name_2 != NULL)
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name_2);
          gms::common::gms_mm_free(ctor_name_2);
      }
      else
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(SFSK_signal).name());
          if(ctor_name_2 != NULL) gms::common::gms_mm_free(ctor_name_2);
      }   

      
     std::printf("[UNIT-TEST]: -- Calling: generate_fsk_signal_avx_u4x\n");
     const std::int32_t ret_status_2 = 
     SFSK_signal.generate_fsk_signal_avx_u4x<sinusoidal_fsk_t::I_channel_bitstream_optimized_path::vector_avx512_u4x_path,
                                                     sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::vector_avx512_u4x_path, 
                                                     sinusoidal_fsk_t::pulse_shaping_function_optimized_path::vector_avx512_path>(IQ_rectw_bitstream_vseq,1,1,true,true,true,true,true);
          if(ret_status_2<0)
          {
             std::printf("[UNIT-TEST]: -- generate_fsk_signal_avx_u4x -- **FAILED**, status=%d\n",ret_status_2);
             return;
          }
          std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
          sinusoidal_fsk_t::create_signal_plot(static_cast<std::uint32_t>(SFSK_signal.m_sfsk_nsamples),SFSK_signal.m_sfsk_signal.m_data,
                                                                          nullptr,"generate_fsk_signal_bitstream_vseq_avx_u4x_test_1_",
                                                                          "Sinusoidal-FSK Signal",false);   
         printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
     
}

__attribute__((hot))
__attribute__((noinline))
void 
unit_test_generate_fsk_signal_avx512_u4x();

void  
unit_test_generate_fsk_signal_avx512_u4x()
{
     using namespace gms::radiolocation;
     using namespace gms;
     constexpr float k{50.0f};
     constexpr float T{200.0f};
     constexpr std::size_t sfsk_nsamples{static_cast<std::size_t>(k*T)}; // 10 periods of 20 samples each
     constexpr std::size_t I_nsamples{sfsk_nsamples};
     constexpr std::size_t Q_nsamples{sfsk_nsamples};
     constexpr std::size_t psfunc_nsamples{sfsk_nsamples};      
     constexpr float Ac{1.0f};
     constexpr float fc{9000.0f};
     constexpr float I_fc{9000.0f};
     constexpr float Q_fc{I_fc};
     constexpr float I_ph0{0.78539816339744830961566084582f};
     constexpr float Q_ph0{0.78539816339744830961566084582f};
     constexpr float ph0{0.0f};
     constexpr std::size_t IQ_nsignals{20ull};
     constexpr std::size_t iq_I_nsamples{10000ull};
     constexpr std::size_t iq_Q_nsamples{10000ull};
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

      int32_t status_1{};
     
      // For breaking into gdb
      //__asm__ ("int3");
      iq_rectw_bitstream_vsequence_t IQ_rectw_bitstream_vseq = iq_rectw_bitstream_vsequence_t(iq_I_nsamples,IQ_nsignals,I_k_values,I_T_values,
                                                                                              I_Ac_values,I_fc_values,I_ph0_values,
                                                                                              iq_Q_nsamples,Q_k_values,Q_T_values,
                                                                                              Q_Ac_values,Q_fc_values,Q_ph0_values,
                                                                                              I_pdf_params,Q_pdf_params);

      char * ctor_name_1{gms::common::demangle(typeid(IQ_rectw_bitstream_vseq).name(),status_1)};
      if(status_1==0 && ctor_name_1 != NULL)
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name_1);
          gms::common::gms_mm_free(ctor_name_1);
      }
      else
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(IQ_rectw_bitstream_vseq).name());
          if(ctor_name_1 != NULL) gms::common::gms_mm_free(ctor_name_1);
      }   
    
      std::int32_t ret_status_1 = IQ_rectw_bitstream_vseq.generate_iq_rectw_bitstream_sequence_avx512_u16x<iq_rectw_bitstream_vsequence_t
                                                                                                      ::IQ_rectw_bitstreams_optim_path
                                                                                                      ::AVX512_U4X_VEC_PATH,
                                                                                                      iq_rectw_bitstream_vsequence_t
                                                                                                      ::IQ_rectw_bitstream_max_signal_components
                                                                                                      ::MAX_SIGNAL_COMPONENTS_25>(I_ch_ret_codes,
                                                                                                                                  Q_ch_ret_codes);
      if(ret_status_1<0)
      {
             std::printf("[UNIT-TEST]: -- generate_iq_rectw_bitstream_sequence_avx512_u16x -- **FAILED**, status=%d\n",ret_status_1);
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

      int32_t status_2{};
      // For breaking into gdb
      //__asm__ ("int3");
      sinusoidal_fsk_t SFSK_signal(sfsk_nsamples,I_nsamples,Q_nsamples,
                                   psfunc_nsamples,k,T,Ac,fc,ph0,
                                   I_fc,Q_fc,I_ph0,Q_ph0);
      char * ctor_name_2{gms::common::demangle(typeid(SFSK_signal).name(),status_2)};
      if(status_2==0 && ctor_name_2 != NULL)
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name_2);
          gms::common::gms_mm_free(ctor_name_2);
      }
      else
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(SFSK_signal).name());
          if(ctor_name_2 != NULL) gms::common::gms_mm_free(ctor_name_2);
      }   

      
     std::printf("[UNIT-TEST]: -- Calling: generate_fsk_signal_avx512_u4x\n");
     const std::int32_t ret_status_2 = 
     SFSK_signal.generate_fsk_signal_avx512_u4x<sinusoidal_fsk_t::I_channel_bitstream_optimized_path::vector_avx512_u4x_path,
                                                     sinusoidal_fsk_t::Q_channel_bitstream_optimized_path::vector_avx512_u4x_path, 
                                                     sinusoidal_fsk_t::pulse_shaping_function_optimized_path::vector_avx512_path>(IQ_rectw_bitstream_vseq,1,1,true,true,true,true,true);
          if(ret_status_2<0)
          {
             std::printf("[UNIT-TEST]: -- generate_fsk_signal_sse_u4x -- **FAILED**, status=%d\n",ret_status_2);
             return;
          }
          std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
          sinusoidal_fsk_t::create_signal_plot(static_cast<std::uint32_t>(SFSK_signal.m_sfsk_nsamples),SFSK_signal.m_sfsk_signal.m_data,
                                                                          nullptr,"generate_fsk_signal_bitstream_vseq_avx512_u4x_test_1_",
                                                                          "Sinusoidal-FSK Signal",false);   
         printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
     
}


int main()
{
    unit_test_generate_fsk_signal_u8x();
    unit_test_generate_fsk_signal_sse_u4x();
    unit_test_generate_fsk_signal_avx_u4x();
    unit_test_generate_fsk_signal_avx512_u4x();
    return 0;
}