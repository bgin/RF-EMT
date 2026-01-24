
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <vector>
#include <thread>
#include "GMS_cpm_pulse_shapers.h"


/*
   icpc -o unit_test_cpm_lsrc_pulse_shapers -O3 -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_avx_memset.h GMS_avx_memset.cpp GMS_avx512_memset.h GMS_avx512_memset.cpp GMS_cephes_sin_cos.h GMS_loop_indices_LUT.h GMS_loop_indices_LUT.cpp GMS_cpm_pulse_shapers.h GMS_cpm_pulse_shapers.cpp unit_test_cpm_lrc_pulse_shapers.cpp
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_avx_memset.h GMS_avx_memset.cpp GMS_avx512_memset.h GMS_avx512_memset.cpp GMS_cephes_sin_cos.h GMS_loop_indices_LUT.h GMS_loop_indices_LUT.cpp GMS_cpm_pulse_shapers.h GMS_cpm_pulse_shapers.cpp unit_test_cpm_lsrc_pulse_shapers.cpp

*/


__attribute__((hot))
__attribute__((noinline))
void unit_test_cpm_lsrc_pulse_shapers();

void unit_test_cpm_lsrc_pulse_shapers()
{
   using namespace gms::radiolocation;
   using namespace gms;
   constexpr std::size_t nTLsamples{1500ull};
   constexpr std::size_t nTsamples{500ull};
   constexpr float       T{500.0f};
   constexpr float       L{3.0f};
   constexpr float       beta{0.854f};
   constexpr float       BbT{0.245f};
   cpm_pulse_shapers_t   lsrc_pulse = cpm_pulse_shapers_t(nTLsamples,nTsamples,
                                                            T,L,beta,BbT);
   int32_t status{};
     //__asm__ ("int3");
   printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);    
   char * ctor_name{gms::common::demangle(typeid(lsrc_pulse).name(),status)};
   if(status==0 && ctor_name != NULL)
   {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
   }
   else
   {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(lsrc_pulse).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
   }  
   printf("[UNIT-TEST]: #1 -- generate_lsrc_pulse_scalar_u8x()\n");
   lsrc_pulse.generate_lsrc_pulse_scalar_u8x();
   std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
   cpm_pulse_shapers_t::create_signal_plot(static_cast<std::int32_t>(lsrc_pulse.m_nTLsamples),lsrc_pulse.m_lsrc_pulse.m_data,nullptr,
                                              "cpm_lsrc_pulse_scalar_u8x_test_1_","L-Spectrally Raised Cosine Pulse",false);
   lsrc_pulse.clear_lsrc_pulse(0.0f);
   printf("[UNIT-TEST]: #2 -- generate_lsrc_pulse_sse_rolled()\n");
   lsrc_pulse.generate_lsrc_pulse_sse_rolled();
   std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
   cpm_pulse_shapers_t::create_signal_plot(static_cast<std::int32_t>(lsrc_pulse.m_nTLsamples),lsrc_pulse.m_lsrc_pulse.m_data,nullptr,
                                              "cpm_lsrc_pulse_sse_rolled_test_1_","L-Spectrally Raised Cosine Pulse",false);
   lsrc_pulse.clear_lsrc_pulse(0.0f);
   printf("[UNIT-TEST]: #2 -- generate_lsrc_pulse_sse_u8x()\n");
   lsrc_pulse.generate_lsrc_pulse_sse_u8x();
   std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
   cpm_pulse_shapers_t::create_signal_plot(static_cast<std::int32_t>(lsrc_pulse.m_nTLsamples),lsrc_pulse.m_lsrc_pulse.m_data,nullptr,
                                              "cpm_lsrc_pulse_sse_u8x_test_1_","L-Spectrally Raised Cosine Pulse",false);
   printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__); 
}

int main()
{
   unit_test_cpm_lsrc_pulse_shapers();
   return (0);
}