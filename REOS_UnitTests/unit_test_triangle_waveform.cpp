#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_triangle_waveform.h"

/*
   icpc -o unit_test_triangle_wave_series -fp-model -std=c++17 fast=2 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_triangle_waveform.h GMS_triangle_waveform.cpp unit_test_triangle_waveform.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_triangle_waveform.h GMS_triangle_waveform.cpp unit_test_triangle_waveform.cpp

*/

__attribute__((hot))
__attribute__((noinline))
void unit_test_trianglew_series_1();

void unit_test_trianglew_series_1()
{
     using namespace gms::radiolocation;
     constexpr std::size_t   n_samples{256ull};
     constexpr std::uint32_t n_waves{70ull};
     constexpr float A{1.0f};
     constexpr float T{8.0f};
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     triangle_waveform_t __trianglew_1__ = triangle_waveform_t(n_samples,n_waves,A,T);
     char * ctor_name{gms::common::demangle(typeid(__trianglew_1__).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(__trianglew_1__).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
     std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     __trianglew_1__.__tw_samples__.info_size_alignment();
     std::printf("[UNIT-TEST:] -- Creating square waveform.\n");
     __trianglew_1__.fourier_series_expansion();
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     triangle_waveform_t::create_signal_plot(__trianglew_1__.__n_samples__,__trianglew_1__.__tw_samples__.m_data,nullptr,
                                              "Triangle_wave_expanded_test_1_","Triangle_Waveform_Expansion",false);
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}

int main()
{
    unit_test_trianglew_series_1();
    return 0;
}