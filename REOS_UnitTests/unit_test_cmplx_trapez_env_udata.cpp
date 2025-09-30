#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_dyn_array.h"
#include "GMS_cmplx_trapezw_env.h"

/*
   icpc -o unit_test_cmplx_trapezw_env_udata -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_indices.h GMS_cmplx_trapezw_env.h GMS_cmplx_trapezw_env.cpp unit_test_cmplx_trapezw_env_udata.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_indices.h GMS_cmplx_trapezw_env.h GMS_cmplx_trapezw_env.cpp unit_test_cmplx_trapezw_env_udata.cpp

*/

__attribute__((hot))
__attribute__((noinline))
void unit_test_cmplx_trapezw_env_rand_data();

void unit_test_cmplx_trapezw_env_rand_data()
{
     using namespace gms::radiolocation;
     using namespace gms;
     constexpr std::size_t   I_n_samples{1024ull};
     constexpr std::uint32_t I_n_K{8};
     constexpr std::size_t   Q_n_samples{1024ull};
     constexpr std::uint32_t Q_n_K{8};
     constexpr float I_a{10.0f};
     constexpr float I_m{5.0f};
     constexpr float I_l{5.0f};
     constexpr float I_c{2.0f};
     constexpr float Q_a{10.f};
     constexpr float Q_m{5.0f};
     constexpr float Q_l{5.0f};
     constexpr float Q_c{2.0f};
     darray_r4_t I_sym;
     darray_r4_t Q_sym;
     std::clock_t seed_I;
     std::clock_t seed_Q;
     constexpr std::size_t I_sym_len{I_n_samples*static_cast<std::size_t>(I_n_K)};
     constexpr std::size_t Q_sym_len{Q_n_samples*static_cast<std::size_t>(Q_n_K)};
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     cmplx_trapezw_env_t __cmplx_trapw_env__ = cmplx_trapezw_env_t(I_n_samples,Q_n_samples,
                                                                   I_n_K,Q_n_K,
                                                                   I_a,I_m,I_l,I_c,
                                                                   Q_a,Q_m,Q_l,Q_c);
     I_sym = darray_r4_t(I_sym_len);
     Q_sym = darray_r4_t(Q_sym_len);

     char * ctor_name{gms::common::demangle(typeid(__cmplx_trapw_env__).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(__cmplx_trapw_env__).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
     //std::printf("[UNIT-TEST:] -- I_n_samples=%llu,I_n_K=%d,Q_n_samples=%llu,Q_n_K=%d\n\n",__cmplx_trapw_env__.__I_n_samples__,
     //            __cmplx_trapw_env__.__I_n_K__,__cmplx_trapw_env__.__Q_n_samples__,__cmplx_trapw_env__.__Q_n_K__);

     std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     __cmplx_trapw_env__.__I_chan__.info_size_alignment();
     __cmplx_trapw_env__.__Q_chan__.info_size_alignment();
     std::printf("[UNIT-TEST:] -- Generating random symbol data: I_channel(-1,1), Q_channel(1,0)\n");
     seed_I = std::clock();
     auto rand_sym_I{std::bind(std::uniform_real_distribution<float>(-1.0f,1.0f),std::mt19937(seed_I))};
     seed_Q = std::clock();
     auto rand_sym_Q{std::bind(std::uniform_real_distribution<float>(0.0f,1.0f),std::mt19937(seed_Q))};
     for(std::size_t __i{0ull}; __i != I_sym_len; ++__i) 
     {
          I_sym.m_data[__i] = rand_sym_I();
          Q_sym.m_data[__i] = rand_sym_Q();
     }
     std::printf("[UNIT-TEST:] -- Calling I/Q trapezoid wave components modulated by the user data.\n");
     std::int32_t stat_I = __cmplx_trapw_env__.chan_I_data_symbol(&I_sym.m_data[0],static_cast<std::uint32_t>(I_n_samples),I_n_K);
     if(stat_I!=0)
     {
          std::printf("[UNIT-TEST:] -- chan_I_data_symbol failed -- error=%d\n",stat_I);
          return;
     }
     std::int32_t stat_Q = __cmplx_trapw_env__.chan_Q_data_symbol(&Q_sym.m_data[0],static_cast<std::uint32_t>(Q_n_samples),Q_n_K);
     if(stat_Q!=0)
     {
          std::printf("[UNIT-TEST:] -- chan_Q_data_symbol failed -- error=%d\n",stat_Q);
          return;
     }

     //std::cout << __cmplx_trapw_env__ << std::endl;
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     cmplx_trapezw_env_t::create_signal_plot(__cmplx_trapw_env__.__I_n_samples__,__cmplx_trapw_env__.__I_chan__.m_data,nullptr,
                                           "chan_I_user_data_test_1_","Complex_Trapezoid_Envelope_I_channel",false);
     cmplx_trapezw_env_t::create_signal_plot(__cmplx_trapw_env__.__Q_n_samples__,__cmplx_trapw_env__.__Q_chan__.m_data,nullptr,
                                           "chan_Q_user_data_test_1_","Complex_Trapezoid_Envelope_Q_channel",false);
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_cmplx_trapezw_env_rand_data_u4x();

void unit_test_cmplx_trapezw_env_rand_data_u4x()
{
     using namespace gms::radiolocation;
     using namespace gms;
     constexpr std::size_t   I_n_samples{1024ull};
     constexpr std::uint32_t I_n_K{32};
     constexpr std::size_t   Q_n_samples{1024ull};
     constexpr std::uint32_t Q_n_K{32};
     constexpr float I_a{10.0f};
     constexpr float I_m{5.0f};
     constexpr float I_l{5.0f};
     constexpr float I_c{2.0f};
     constexpr float Q_a{10.f};
     constexpr float Q_m{5.0f};
     constexpr float Q_l{5.0f};
     constexpr float Q_c{2.0f};
     darray_r4_t I_sym;
     darray_r4_t Q_sym;
     std::clock_t seed_I;
     std::clock_t seed_Q;
     constexpr std::size_t I_sym_len{I_n_samples*static_cast<std::size_t>(I_n_K)};
     constexpr std::size_t Q_sym_len{Q_n_samples*static_cast<std::size_t>(Q_n_K)};
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     cmplx_trapezw_env_t __cmplx_trapw_env__ = cmplx_trapezw_env_t(I_n_samples,Q_n_samples,
                                                                   I_n_K,Q_n_K,
                                                                   I_a,I_m,I_l,I_c,
                                                                   Q_a,Q_m,Q_l,Q_c);
     I_sym = darray_r4_t(I_sym_len);
     Q_sym = darray_r4_t(Q_sym_len);

     char * ctor_name{gms::common::demangle(typeid(__cmplx_trapw_env__).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(__cmplx_trapw_env__).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
     //std::printf("[UNIT-TEST:] -- I_n_samples=%llu,I_n_K=%d,Q_n_samples=%llu,Q_n_K=%d\n\n",__cmplx_trapw_env__.__I_n_samples__,
     //            __cmplx_trapw_env__.__I_n_K__,__cmplx_trapw_env__.__Q_n_samples__,__cmplx_trapw_env__.__Q_n_K__);

     std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     __cmplx_trapw_env__.__I_chan__.info_size_alignment();
     __cmplx_trapw_env__.__Q_chan__.info_size_alignment();
     std::printf("[UNIT-TEST:] -- Generating random symbol data: I_channel(-1,1), Q_channel(1,0)\n");
     seed_I = std::clock();
     auto rand_sym_I{std::bind(std::uniform_real_distribution<float>(-1.0f,1.0f),std::mt19937(seed_I))};
     seed_Q = std::clock();
     auto rand_sym_Q{std::bind(std::uniform_real_distribution<float>(0.0f,1.0f),std::mt19937(seed_Q))};
     for(std::size_t __i{0ull}; __i != I_sym_len; ++__i) 
     {
          I_sym.m_data[__i] = rand_sym_I();
          Q_sym.m_data[__i] = rand_sym_Q();
     }
     std::printf("[UNIT-TEST:] -- Calling I/Q trapezoid wave components modulated by the user data.\n");
     std::int32_t stat_I = __cmplx_trapw_env__.chan_I_data_symbol_u4x(&I_sym.m_data[0],static_cast<std::uint32_t>(I_n_samples),I_n_K);
     if(stat_I!=0)
     {
          std::printf("[UNIT-TEST:] -- chan_I_data_symbol failed -- error=%d\n",stat_I);
          return;
     }
     std::int32_t stat_Q = __cmplx_trapw_env__.chan_Q_data_symbol_u4x(&Q_sym.m_data[0],static_cast<std::uint32_t>(Q_n_samples),Q_n_K);
     if(stat_Q!=0)
     {
          std::printf("[UNIT-TEST:] -- chan_Q_data_symbol failed -- error=%d\n",stat_Q);
          return;
     }

     //std::cout << __cmplx_trapw_env__ << std::endl;
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     cmplx_trapezw_env_t::create_signal_plot(__cmplx_trapw_env__.__I_n_samples__,__cmplx_trapw_env__.__I_chan__.m_data,nullptr,
                                           "chan_I_user_data_test_2_","Complex_Trapezoid_Envelope_I_channel",false);
     cmplx_trapezw_env_t::create_signal_plot(__cmplx_trapw_env__.__Q_n_samples__,__cmplx_trapw_env__.__Q_chan__.m_data,nullptr,
                                           "chan_Q_user_data_test_2_","Complex_Trapezoid_Envelope_Q_channel",false);
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_cmplx_trapezw_env_rand_data_2();

void unit_test_cmplx_trapezw_env_rand_data_2()
{
     using namespace gms::radiolocation;
     using namespace gms;
     constexpr std::size_t   I_n_samples{1024ull};
     constexpr std::uint32_t I_n_K{8};
     constexpr std::size_t   Q_n_samples{1024ull};
     constexpr std::uint32_t Q_n_K{8};
     constexpr float I_a{10.0f};
     constexpr float I_m{5.0f};
     constexpr float I_l{5.0f};
     constexpr float I_c{2.0f};
     constexpr float Q_a{10.f};
     constexpr float Q_m{5.0f};
     constexpr float Q_l{5.0f};
     constexpr float Q_c{2.0f};
     std::clock_t seed_I;
     std::clock_t seed_Q;
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     cmplx_trapezw_env_t __cmplx_trapw_env__ = cmplx_trapezw_env_t(I_n_samples,Q_n_samples,
                                                                   I_n_K,Q_n_K,
                                                                   I_a,I_m,I_l,I_c,
                                                                   Q_a,Q_m,Q_l,Q_c);
    
     char * ctor_name{gms::common::demangle(typeid(__cmplx_trapw_env__).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(__cmplx_trapw_env__).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
     //std::printf("[UNIT-TEST:] -- I_n_samples=%llu,I_n_K=%d,Q_n_samples=%llu,Q_n_K=%d\n\n",__cmplx_trapw_env__.__I_n_samples__,
     //            __cmplx_trapw_env__.__I_n_K__,__cmplx_trapw_env__.__Q_n_samples__,__cmplx_trapw_env__.__Q_n_K__);

     std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     __cmplx_trapw_env__.__I_chan__.info_size_alignment();
     __cmplx_trapw_env__.__Q_chan__.info_size_alignment();
     std::printf("[UNIT-TEST:] -- Generating random symbol data: I_channel(-1,1), Q_channel(1,0)\n");
     seed_I = std::clock();
     auto rand_sym_I{std::bind(std::uniform_real_distribution<float>(-1.0f,1.0f),std::mt19937(seed_I))};
     seed_Q = std::clock();
     auto rand_sym_Q{std::bind(std::uniform_real_distribution<float>(0.0f,1.0f),std::mt19937(seed_Q))};
     std::printf("[UNIT-TEST:] -- Calling I/Q trapezoid wave components modulated by the user data.\n");
     std::int32_t stat_I = __cmplx_trapw_env__.chan_I_rand_data_symbol(rand_sym_I);
     if(stat_I!=0)
     {
          std::printf("[UNIT-TEST:] -- chan_I_data_symbol failed -- error=%d\n",stat_I);
          return;
     }
     std::int32_t stat_Q = __cmplx_trapw_env__.chan_Q_rand_data_symbol(rand_sym_Q);
     if(stat_Q!=0)
     {
          std::printf("[UNIT-TEST:] -- chan_Q_data_symbol failed -- error=%d\n",stat_Q);
          return;
     }

     //std::cout << __cmplx_trapw_env__ << std::endl;
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     cmplx_trapezw_env_t::create_signal_plot(__cmplx_trapw_env__.__I_n_samples__,__cmplx_trapw_env__.__I_chan__.m_data,nullptr,
                                           "chan_I_user_data_test_3_","Complex_Trapezoid_Envelope_I_channel",false);
     cmplx_trapezw_env_t::create_signal_plot(__cmplx_trapw_env__.__Q_n_samples__,__cmplx_trapw_env__.__Q_chan__.m_data,nullptr,
                                           "chan_Q_user_data_test_3_","Complex_Trapezoid_Envelope_Q_channel",false);
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}



int main()
{   
    unit_test_cmplx_trapezw_env_rand_data_2();
    unit_test_cmplx_trapezw_env_rand_data_u4x();
    return 0;
}