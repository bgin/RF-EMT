#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_dyn_array.h"
#include "GMS_am_bb_cmplx_trapez_signal_fp16.h"

/*
   icpc -o unit_test_am_bb_cmplx_trapez_signal_fp16 -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_half.h GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_indices.h GMS_am_bb_cmplx_trapez_signal_fp16.h GMS_am_bb_cmplx_trapez_signal_fp16.cpp unit_test_am_bb_cmplx_trapez_signal_fp16.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_half.h GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_indices.h GMS_am_bb_cmplx_trapez_signal_fp16.h GMS_am_bb_cmplx_trapez_signal_fp16.cpp unit_test_am_bb_cmplx_trapez_signal_fp16.cpp

*/

__attribute__((hot))
__attribute__((noinline))
void unit_test_am_bb_cmplx_trapez_signal_fp16();

void unit_test_am_bb_cmplx_trapez_signal_fp16()
{
     using namespace gms::radiolocation;
     using namespace gms;
     constexpr std::size_t   I_n_samples{1024ull};
     constexpr std::uint32_t I_n_K{8};
     constexpr std::size_t   Re_n_samples{I_n_samples/2ull};
     constexpr float I_a{10.0f};
     constexpr float I_m{5.0f};
     constexpr float I_l{5.0f};
     constexpr float I_c{2.0f};
     //float re_buf[Re_n_samples];
     //float im_buf[Re_n_samples];
     darray_r4_t I_sym; 
     std::clock_t seed_I;
     
     constexpr std::size_t I_sym_len{I_n_samples*static_cast<std::size_t>(I_n_K)};
    
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     am_bb_cmplx_trapez_signal_fp16_t symbol_1 = am_bb_cmplx_trapez_signal_fp16_t(I_n_samples,I_n_K,I_a,I_m,I_l,I_c);
     I_sym = darray_r4_t(I_sym_len);
     char * ctor_name{gms::common::demangle(typeid(symbol_1).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(symbol_1).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
     //std::printf("[UNIT-TEST:] -- I_n_samples=%llu,I_n_K=%d,Q_n_samples=%llu,Q_n_K=%d\n\n",__cmplx_trapw_env__.__I_n_samples__,
     //            __cmplx_trapw_env__.__I_n_K__,__cmplx_trapw_env__.__Q_n_samples__,__cmplx_trapw_env__.__Q_n_K__);

     std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     symbol_1.m_sig_samples.info_size_alignment();
     std::printf("[UNIT-TEST:] -- Generating random symbol data: I_channel(0,1)\n");
     seed_I = std::clock();
     auto rand_sym_I{std::bind(std::uniform_int_distribution<int32_t>(0,1),std::mt19937(seed_I))};
     
     for(std::size_t __i{0ull}; __i != I_sym_len; ++__i) 
     {
          I_sym.m_data[__i] =static_cast<float>(rand_sym_I());
          
          //printf("%.7f\n",I_sym.m_data[__i]);
     }
     std::printf("[UNIT-TEST:] -- Calling  trapezoid wave components modulated by the user data.\n");
     std::int32_t stat_I = symbol_1.create_signal_user_data(&I_sym.m_data[0],static_cast<std::uint32_t>(I_n_samples),I_n_K);
     if(stat_I!=0)
     {
          std::printf("[UNIT-TEST:] -- chan_I_data_symbol failed -- error=%d\n",stat_I);
          return;
     }
     
    /*for(std::size_t __i{0ull}; __i != Re_n_samples; ++__i) 
     {
         re_buf[__i] = symbol_1.m_sig_samples.m_data[__i].real();
         im_buf[__i] = symbol_1.m_sig_samples.m_data[__i].imag();
     }
     */
     const half_float::half * __restrict__ p_re{&reinterpret_cast<half_float::half(&)[2]>(symbol_1.m_sig_samples.m_data[0])[0]};
     const half_float::half * __restrict__ p_im{&reinterpret_cast<half_float::half(&)[2]>(symbol_1.m_sig_samples.m_data[0])[1]};
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     am_bb_cmplx_trapez_signal_fp16_t::create_signal_plot(Re_n_samples,p_re,nullptr,
                                              "am_bb_cmplx_trapez_signal_fp16_I_test_1_","I-channel",false);
     am_bb_cmplx_trapez_signal_fp16_t::create_signal_plot(Re_n_samples,p_im,nullptr,
                                              "am_bb_cmplx_trapez_signal_fp16_Q_test_1_","Q-channel",false);                                         
     
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_am_bb_cmplx_trapez_noisy_signal_fp16();

void unit_test_am_bb_cmplx_trapez_noisy_signal_fp16()
{
     using namespace gms::radiolocation;
     using namespace gms;
     constexpr std::size_t   I_n_samples{1024ull};
     constexpr std::uint32_t I_n_K{8};
     constexpr std::size_t   Re_n_samples{I_n_samples/2ull};
     constexpr float I_a{10.0f};
     constexpr float I_m{5.0f};
     constexpr float I_l{5.0f};
     constexpr float I_c{2.0f};
     //float re_buf[Re_n_samples];
     //float im_buf[Re_n_samples];
     am_bb_cmplx_trapez_signal_fp16_pdf_params_t pdf_params;
     darray_r4_t I_sym; 
     std::clock_t seed_I;
     float scale{0.001f};
     constexpr std::size_t I_sym_len{I_n_samples*static_cast<std::size_t>(I_n_K)};
    
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     am_bb_cmplx_trapez_signal_fp16_t symbol_1 = am_bb_cmplx_trapez_signal_fp16_t(I_n_samples,I_n_K,I_a,I_m,I_l,I_c);
     I_sym = darray_r4_t(I_sym_len);
     pdf_params.uni_real_a_r = 0.0f;
     pdf_params.uni_real_b_r = 1.0f;
     char * ctor_name{gms::common::demangle(typeid(symbol_1).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(symbol_1).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
     //std::printf("[UNIT-TEST:] -- I_n_samples=%llu,I_n_K=%d,Q_n_samples=%llu,Q_n_K=%d\n\n",__cmplx_trapw_env__.__I_n_samples__,
     //            __cmplx_trapw_env__.__I_n_K__,__cmplx_trapw_env__.__Q_n_samples__,__cmplx_trapw_env__.__Q_n_K__);

     std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     symbol_1.m_sig_samples.info_size_alignment();
     std::printf("[UNIT-TEST:] -- Generating random symbol data: I_channel(0,1)\n");
     seed_I = std::clock();
     auto rand_sym_I{std::bind(std::uniform_int_distribution<int32_t>(0,1),std::mt19937(seed_I))};
     
     for(std::size_t __i{0ull}; __i != I_sym_len; ++__i) 
     {
          I_sym.m_data[__i] =static_cast<float>(rand_sym_I());
          
          //printf("%.7f\n",I_sym.m_data[__i]);
     }
     std::printf("[UNIT-TEST:] -- Calling  trapezoid wave components modulated by the user data and noise corrupted.\n");
     //std::int32_t stat_I = symbol_1.create_signal_user_data(&I_sym.m_data[0],static_cast<std::uint32_t>(I_n_samples),I_n_K);
     std::int32_t stat_I{symbol_1.create_noisy_signal_user_data(pdf_params,scale,&I_sym.m_data[0],static_cast<std::uint32_t>(I_n_samples),I_n_K)};
     if(stat_I!=0)
     {
          std::printf("[UNIT-TEST:] -- chan_I_data_symbol failed -- error=%d\n",stat_I);
          return;
     }
     
    /*for(std::size_t __i{0ull}; __i != Re_n_samples; ++__i) 
     {
         re_buf[__i] = symbol_1.m_sig_samples.m_data[__i].real();
         im_buf[__i] = symbol_1.m_sig_samples.m_data[__i].imag();
     }
     */
     const half_float::half * __restrict__ p_re{&reinterpret_cast<half_float::half(&)[2]>(symbol_1.m_sig_samples.m_data[0])[0]};
     const half_float::half * __restrict__ p_im{&reinterpret_cast<half_float::half(&)[2]>(symbol_1.m_sig_samples.m_data[0])[1]};
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     am_bb_cmplx_trapez_signal_fp16_t::create_signal_plot(Re_n_samples,p_re,nullptr,
                                              "am_bb_cmplx_trapez_noisy_signal_fp16_I_test_1_","I-channel",false);
     am_bb_cmplx_trapez_signal_fp16_t::create_signal_plot(Re_n_samples,p_im,nullptr,
                                              "am_bb_cmplx_trapez_noisy_signal_fp16_Q_test_1_","Q-channel",false);                                         
     
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_am_bb_cmplx_trapez_noisy_signal__fp16_2();

void unit_test_am_bb_cmplx_trapez_noisy_signal_fp16_2()
{
     using namespace gms::radiolocation;
     using namespace gms;
     constexpr std::size_t   I_n_samples{256ull};
     constexpr std::uint32_t I_n_K{16};
     constexpr std::size_t   Re_n_samples{I_n_samples/2ull};
     constexpr float I_a{10.0f};
     constexpr float I_m{5.0f};
     constexpr float I_l{5.0f};
     constexpr float I_c{2.0f};
     //float re_buf[Re_n_samples];
     //float im_buf[Re_n_samples];
     am_bb_cmplx_trapez_signal_fp16_pdf_params_t pdf_params;
     darray_r4_t I_sym; 
     std::clock_t seed_I;
     float scale{0.01f};
     constexpr std::size_t I_sym_len{I_n_samples*static_cast<std::size_t>(I_n_K)};
    
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     am_bb_cmplx_trapez_signal_fp16_t symbol_1 = am_bb_cmplx_trapez_signal_fp16_t(I_n_samples,I_n_K,I_a,I_m,I_l,I_c);
     I_sym = darray_r4_t(I_sym_len);
     pdf_params.uni_real_a_r = 0.0f;
     pdf_params.uni_real_b_r = 1.0f;
     pdf_params.norm_mu_r    = 0.0f;
     pdf_params.norm_sigma_r = 1.0f;
     pdf_params.cauchy_a_r   = 0.0f;
     pdf_params.cauchy_b_r   = 1.0f;
     char * ctor_name{gms::common::demangle(typeid(symbol_1).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(symbol_1).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
     //std::printf("[UNIT-TEST:] -- I_n_samples=%llu,I_n_K=%d,Q_n_samples=%llu,Q_n_K=%d\n\n",__cmplx_trapw_env__.__I_n_samples__,
     //            __cmplx_trapw_env__.__I_n_K__,__cmplx_trapw_env__.__Q_n_samples__,__cmplx_trapw_env__.__Q_n_K__);

     std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     symbol_1.m_sig_samples.info_size_alignment();
     std::printf("[UNIT-TEST:] -- Generating random symbol data: I_channel(0,1)\n");
     seed_I = std::clock();
     auto rand_sym_I{std::bind(std::uniform_int_distribution<int32_t>(0,1),std::mt19937(seed_I))};
     
     for(std::size_t __i{0ull}; __i != I_sym_len; ++__i) 
     {
          I_sym.m_data[__i] =static_cast<float>(rand_sym_I());
          
          //printf("%.7f\n",I_sym.m_data[__i]);
     }
     std::printf("[UNIT-TEST:] -- Calling  trapezoid wave components modulated by the user data and noise corrupted.\n");
     //std::int32_t stat_I = symbol_1.create_signal_user_data(&I_sym.m_data[0],static_cast<std::uint32_t>(I_n_samples),I_n_K);
     std::int32_t stat_I{symbol_1.create_noisy_signal_user_data(pdf_params,am_bb_cmplx_trapez_signal_fp16_rand_distr::cauchy,scale,&I_sym.m_data[0],static_cast<std::uint32_t>(I_n_samples),I_n_K)};
     if(stat_I!=0)
     {
          std::printf("[UNIT-TEST:] -- chan_I_data_symbol failed -- error=%d\n",stat_I);
          return;
     }
     
    /*for(std::size_t __i{0ull}; __i != Re_n_samples; ++__i) 
     {
         re_buf[__i] = symbol_1.m_sig_samples.m_data[__i].real();
         im_buf[__i] = symbol_1.m_sig_samples.m_data[__i].imag();
     }
     */
     const half_float::half * __restrict__ p_re{&reinterpret_cast<half_float::half(&)[2]>(symbol_1.m_sig_samples.m_data[0])[0]};
     const half_float::half * __restrict__ p_im{&reinterpret_cast<half_float::half(&)[2]>(symbol_1.m_sig_samples.m_data[0])[1]};
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     am_bb_cmplx_trapez_signal_fp16_t::create_signal_plot(Re_n_samples,p_re,nullptr,
                                              "am_bb_cmplx_trapez_noisy_cauchy_signal_fp16_I_test_1_","I-channel",false);
     am_bb_cmplx_trapez_signal_fp16_t::create_signal_plot(Re_n_samples,p_im,nullptr,
                                              "am_bb_cmplx_trapez_noisy_cauchy_signal_fp16_Q_test_1_","Q-channel",false);                                         
     
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}



int main()
{
     unit_test_am_bb_cmplx_trapez_signal_fp16();
     unit_test_am_bb_cmplx_trapez_noisy_signal_fp16();
     unit_test_am_bb_cmplx_trapez_noisy_signal_fp16_2();
     return 0;
}
