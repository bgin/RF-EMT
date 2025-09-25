#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_sawtooth_waveform.h"

/*
   icpc -o unit_test_sawtooth_wave_series -fp-model -std=c++17 fast=2 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_sawtooth_waveform.h GMS_sawtooth_waveform.cpp unit_test_sawtooth_waveform.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_sawtooth_waveform.h GMS_sawtooth_waveform.cpp unit_test_sawtooth_waveform.cpp

*/

__attribute__((hot))
__attribute__((noinline))
void unit_test_sawtoothw_series_1();

void unit_test_sawtoothw_series_1()
{
     using namespace gms::radiolocation;
     constexpr std::size_t   n_samples{256ull};
     constexpr std::uint32_t n_waves{64ull};
     constexpr float A{1.0f};
     constexpr float T{1.0f};
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     sawtooth_waveform_t __sawtoothw_1__ = sawtooth_waveform_t(n_samples,n_waves,A,T);
     char * ctor_name{gms::common::demangle(typeid(__sawtoothw_1__).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(__sawtoothw_1__).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
     std::cout << __sawtoothw_1__ << std::endl;
     std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     __sawtoothw_1__.__sw_samples__.info_size_alignment();
     std::printf("[UNIT-TEST:] -- Creating sawtooth waveform.\n");
     std::uint32_t res = __sawtoothw_1__.fourier_series_expansion(0);
     if(res < 0)
     {
          std::printf("[UNIT-TEST:] -- fourier_series_expansion: FAILED with an error=%d\n",res);
          return;
     }
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     sawtooth_waveform_t::create_signal_plot(__sawtoothw_1__.__n_samples__,__sawtoothw_1__.__sw_samples__.m_data,nullptr,
                                              "sawtooth_wave_expanded_test_7_","Sawtooth_Waveform_Expansion",false);
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_sawtoothw_series_2();

void unit_test_sawtoothw_series_2()
{
     using namespace gms::radiolocation;
     constexpr std::size_t   n_samples{8192ull};
     constexpr std::uint32_t n_waves{64ull};
     constexpr float A{1.0f};
     constexpr float T{16.0f};
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     sawtooth_waveform_t __sawtoothw_1__ = sawtooth_waveform_t(n_samples,n_waves,A,T);
     char * ctor_name{gms::common::demangle(typeid(__sawtoothw_1__).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(__sawtoothw_1__).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
     std::cout << __sawtoothw_1__ << std::endl;
     std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     __sawtoothw_1__.__sw_samples__.info_size_alignment();
     std::printf("[UNIT-TEST:] -- Creating sawtooth waveform.\n");
     std::uint32_t res = __sawtoothw_1__.fourier_series_expansion_u4x(0);
     if(res < 0)
     {
          std::printf("[UNIT-TEST:] -- fourier_series_expansion: FAILED with an error=%d\n",res);
          return;
     }
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     sawtooth_waveform_t::create_signal_plot(__sawtoothw_1__.__n_samples__,__sawtoothw_1__.__sw_samples__.m_data,nullptr,
                                              "fourier_series_expansion_u4x_test_2_","Sawtooth_Waveform_Expansion",false);
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}

int main()
{
    //unit_test_sawtoothw_series_1();
    unit_test_sawtoothw_series_2();
    return 0;
}