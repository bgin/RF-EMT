#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_dyn_array.h"
#include "GMS_cmplx_trapezw_env.h"

/*
   icpc -o unit_test_cmplx_trapezw_env -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_cmplx_trapezw_env.h GMS_cmplx_trapezw_env.cpp unit_test_cmplx_trapezw_env.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cmplx_trapezw_env.h GMS_cmplx_trapezw_env.cpp unit_test_cmplx_trapezw_env.cpp

*/

__attribute__((hot))
__attribute__((noinline))
void unit_test_cmplx_trapezw_env_s1024_uncode();

void unit_test_cmplx_trapezw_env_s1024_uncode()
{
     using namespace gms::radiolocation;
     using namespace gms;
     constexpr std::size_t   I_n_samples{1024ull};
     constexpr std::uint32_t I_n_K{32ull};
     constexpr std::size_t   Q_n_samples{1024ull};
     constexpr std::uint32_t Q_n_K{32ull};
     constexpr float I_a{10.0f};
     constexpr float I_m{5.0f};
     constexpr float I_l{5.0f};
     constexpr float I_c{2.0f};
     constexpr float Q_a{2.0f};
     constexpr float Q_m{10.0f};
     constexpr float Q_l{10.0f};
     constexpr float Q_c{5.0f};
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     cmplx_trapezw_env_t __cmplx_trapw_env__ = cmplx_trapezw_env_t(I_n_samples,I_n_K,
                                                                   Q_n_samples,Q_n_K,
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
     std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     __cmplx_trapw_env__.__I_chan__.info_size_alignment();
     __cmplx_trapw_env__.__Q_chan__.info_size_alignment();
     std::printf("[UNIT-TEST:] -- Calling I/Q uncoded trapezoid wave components.\n");
     __cmplx_trapw_env__.cmplx_env_I_chan_uncoded();
     __cmplx_trapw_env__.cmplx_env_Q_chan_uncoded();
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     cmplx_trapezw_env_t::create_signal_plot(__cmplx_trapw_env__.__I_n_samples__,__cmplx_trapw_env__.__I_chan__.m_data,nullptr,
                                             "cmplx_env_I_chan_uncode_1024_test_1_","Complex_Trapezoid_Envelope_I_channel",false);
     cmplx_trapezw_env_t::create_signal_plot(__cmplx_trapw_env__.__Q_n_samples__,__cmplx_trapw_env__.__Q_chan__.m_data,nullptr,
                                             "cmplx_env_Q_chan_uncode_1024_test_1_","Complex_Trapezoid_Envelope_Q_channel",false);
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}

int main()
{   
    unit_test_cmplx_trapezw_env_s1024_uncode();
    return 0;
}