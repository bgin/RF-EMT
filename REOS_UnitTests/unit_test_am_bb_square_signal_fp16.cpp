#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_dyn_array.h"
#include "GMS_am_bb_square_signal_fp16.h"

/*
   icpc -o unit_test_am_bb_square_signal_fp16 -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h GMS_half.h GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_indices.h GMS_am_bb_square_signal_fp16.h GMS_am_bb_square_signal_fp16.cpp unit_test_am_bb_square_signal_fp16.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_half.h GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_indices.h GMS_am_bb_square_signal_fp16.h GMS_am_bb_square_signal_fp16.cpp unit_test_am_bb_square_signal_fp16.cpp

*/

__attribute__((hot))
__attribute__((noinline))
void unit_test_am_bb_square_signal_fp16_rand_data_1();

void unit_test_am_bb_square_signal_fp16_rand_data_1()
{
     using namespace gms::radiolocation;
     using namespace gms;
     constexpr std::size_t   I_n_samples{1024ull};
     constexpr std::uint32_t I_n_K{8};
     constexpr float A{1.0f};
     constexpr float n{16.0f};
     constexpr float P{8.0f};
     darray_r4_t I_sym; 
     std::clock_t seed_I;
     constexpr std::size_t I_sym_len{I_n_samples*static_cast<std::size_t>(I_n_K)};
    
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     am_bb_square_signal_fp16_t symbol_1 = am_bb_square_signal_fp16_t(I_n_samples,I_n_K,n,A,P);
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
     std::printf("[UNIT-TEST:] -- Generating random symbol data: I_channel(-1,1)\n");
     seed_I = std::clock();
     auto rand_sym_I{std::bind(std::uniform_real_distribution<float>(-1.0f,1.0f),std::mt19937(seed_I))};
     
     for(std::size_t __i{0ull}; __i != I_sym_len; ++__i) 
     {
          I_sym.m_data[__i] = rand_sym_I();
         
     }
     std::printf("[UNIT-TEST:] -- Calling square wave components modulated by the user data.\n");
     std::int32_t stat_I = symbol_1.create_signal_user_data(&I_sym.m_data[0],1u,static_cast<std::uint32_t>(I_n_samples),I_n_K);
     if(stat_I!=0)
     {
          std::printf("[UNIT-TEST:] -- create_signal_user_data: **FAILED** -- error=%d\n",stat_I);
          return;
     }
     
     //std::cout << __cmplx_trapw_env__ << std::endl;
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     am_bb_square_signal_fp16_t::create_signal_plot(symbol_1.m_nsamples,symbol_1.m_sig_samples.m_data,nullptr,
                                              "am_bb_square_signal_fp16_user_data_test_1_","AM_Baseband_SquareWave_Signal",false);
     
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_am_bb_square_signal_fp16_rand_data_2();

void unit_test_am_bb_square_signal_fp16_rand_data_2()
{
     using namespace gms::radiolocation;
     using namespace gms;
     constexpr std::size_t   I_n_samples{128ull};
     constexpr std::uint32_t I_n_K{16};
     constexpr float A{1.0f};
     constexpr float n{16.0f};
     constexpr float P{8.0f};
     
     darray_r4_t I_sym; 
     std::clock_t seed_I;
     constexpr std::size_t I_sym_len{I_n_samples*static_cast<std::size_t>(I_n_K)};
    
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     am_bb_square_signal_fp16_t symbol_1 = am_bb_square_signal_fp16_t(I_n_samples,I_n_K,n,A,P);
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
     std::printf("[UNIT-TEST:] -- Generating random symbol data: I_channel(-1,1)\n");
     seed_I = std::clock();
     auto rand_sym_I{std::bind(std::uniform_real_distribution<float>(-1.0f,1.0f),std::mt19937(seed_I))};
     
     for(std::size_t __i{0ull}; __i != I_sym_len; ++__i) 
     {
          I_sym.m_data[__i] = rand_sym_I();
         
     }
     std::printf("[UNIT-TEST:] -- Calling  square wave components modulated by the user data.\n");
     std::int32_t stat_I = symbol_1.create_signal_user_data(&I_sym.m_data[0],0u,static_cast<std::uint32_t>(I_n_samples),I_n_K);
     if(stat_I!=0)
     {
          std::printf("[UNIT-TEST:] -- create_signal_user_data: **FAILED** -- error=%d\n",stat_I);
          return;
     }
     
     //std::cout << __cmplx_trapw_env__ << std::endl;
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     am_bb_square_signal_fp16_t::create_signal_plot(symbol_1.m_nsamples,symbol_1.m_sig_samples.m_data,nullptr,
                                              "am_bb_square_signal_fp16_user_data_test_2_","AM_Baseband_SquareWave_Signal",false);
     
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}


__attribute__((hot))
__attribute__((noinline))
void unit_test_am_bb_square_signal_fp16_rand_data_3();

void unit_test_am_bb_square_signal_fp16_rand_data_3()
{
     using namespace gms::radiolocation;
     using namespace gms;
     constexpr std::size_t   I_n_samples{1024ull};
     constexpr std::uint32_t I_n_K{64};
     constexpr float A{1.0f};
     constexpr float n{16.0f};
     constexpr float P{8.0f};
     darray_r4_t I_sym; 
     std::clock_t seed_I;
     constexpr std::size_t I_sym_len{I_n_samples*static_cast<std::size_t>(I_n_K)};
    
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     am_bb_square_signal_fp16_t symbol_1 = am_bb_square_signal_fp16_t(I_n_samples,I_n_K,n,A,P);
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
     std::printf("[UNIT-TEST:] -- Generating random symbol data: I_channel(-1,1)\n");
     seed_I = std::clock();
     auto rand_sym_I{std::bind(std::uniform_real_distribution<float>(-1.0f,1.0f),std::mt19937(seed_I))};
     
     for(std::size_t __i{0ull}; __i != I_sym_len; ++__i) 
     {
          I_sym.m_data[__i] = rand_sym_I();
         
     }
     
     std::printf("[UNIT-TEST:] -- Calling  square wave components modulated by the user data.\n");
     std::int32_t stat_I = symbol_1.create_signal_user_data_u4x(&I_sym.m_data[0],0u,static_cast<std::uint32_t>(I_n_samples),I_n_K);
     if(stat_I!=0)
     {
          std::printf("[UNIT-TEST:] -- create_signal_user_data: **FAILED** -- error=%d\n",stat_I);
          return;
     }
     
     //std::cout << __cmplx_trapw_env__ << std::endl;
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     am_bb_square_signal_fp16_t::create_signal_plot(symbol_1.m_nsamples,symbol_1.m_sig_samples.m_data,nullptr,
                                              "am_bb_square_signal_fp16_user_data_test_3_","AM_Baseband_SquareWave_Signal",false);
     
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}

int main()
{
    unit_test_am_bb_square_signal_fp16_rand_data_1();
    unit_test_am_bb_square_signal_fp16_rand_data_2();
    unit_test_am_bb_square_signal_fp16_rand_data_3();
    return 0;
}
