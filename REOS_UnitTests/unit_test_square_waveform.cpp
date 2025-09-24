#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_square_waveform.h"

/*
   icpc -o unit_test_square_wave_series -fp-model -std=c++17 fast=2 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_square_waveform.h GMS_square_waveform.cpp unit_test_square_waveform.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_square_waveform.h GMS_square_waveform.cpp unit_test_square_waveform.cpp

*/

__attribute__((hot))
__attribute__((noinline))
void unit_test_squarew_series_1();

void unit_test_squarew_series_1()
{
     using namespace gms::radiolocation;
     constexpr std::size_t   n_samples{256ull};
     constexpr std::uint32_t n_waves{70ull};
     constexpr float A{1.0f};
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     square_waveform_t __squarew_1__ = square_waveform_t(n_samples,n_waves,A);
     char * ctor_name{gms::common::demangle(typeid(__squarew_1__).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(__squarew_1__).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
     std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     __squarew_1__.__sw_samples__.info_size_alignment();
     std::printf("[UNIT-TEST:] -- Creating square waveform.\n");
     __squarew_1__.fourier_series_expansion();
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     square_waveform_t::create_signal_plot(__squarew_1__.__n_samples__,__squarew_1__.__sw_samples__.m_data,nullptr,
                                              "square_wave_expanded_test_4_","Square_Waveform_Expansion",false);
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}

int main()
{
    unit_test_squarew_series_1();
    return 0;
}