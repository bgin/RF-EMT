#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_rectangular_waveform_fp16.h"

/*
   icpc -o unit_test_rect_wave_fp16 -fp-model -std=c++17 fast=2 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h GMS_half.h GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_rectangular_waveform_fp16.h GMS_rectangular_waveform_fp16.cpp unit_test_rectangular_wave__fp16.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h GMS_half.h GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_rectangular_waveform_fp16.h GMS_rectangular_waveform_fp16.cpp unit_test_rectangular_wave_fp16.cpp

*/

__attribute__((hot))
__attribute__((noinline))
void unit_test_rectangw_fp16_series_1();

void unit_test_rectangw_fp16_series_1()
{
     using namespace gms::radiolocation;
     constexpr std::size_t   n_samples{1024ull};
     constexpr std::uint32_t n_waves{32ull};
     constexpr float T{512.0f};
     constexpr float A{1.0f};
     constexpr float rho{256.0f};
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     rectangular_waveform_fp16_t __rectw_1__  = rectangular_waveform_fp16_t(n_samples,n_waves,
                                                      T,rho,A);
     char * ctor_name{gms::common::demangle(typeid(__rectw_1__).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(__rectw_1__).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
     std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     __rectw_1__.__rw_samples__.info_size_alignment();
     std::printf("[UNIT-TEST:] -- Creating single trapezoid waveform.\n");
     __rectw_1__.fourier_series_expansion(0);
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     rectangular_waveform_fp16_t::create_signal_plot(__rectw_1__.__n_samples__,__rectw_1__.__rw_samples__.m_data,nullptr,
                                              "rectangular_wave_fp16_expanded_test_4_","Rectangular_Waveform_FP16_Expansion",false);
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_rectangw_fp16_series_2();

void unit_test_rectangw_fp16_series_2()
{
     using namespace gms::radiolocation;
     constexpr std::size_t   n_samples{1024ull};
     constexpr std::uint32_t n_waves{32ull};
     constexpr float T{2.0f};
     constexpr float A{1.0f};
     constexpr float rho{1.0f};
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     rectangular_waveform_fp16_t __rectw_1__  = rectangular_waveform_fp16_t(n_samples,n_waves,
                                                      T,rho,A);
     char * ctor_name{gms::common::demangle(typeid(__rectw_1__).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(__rectw_1__).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
     std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     __rectw_1__.__rw_samples__.info_size_alignment();
     std::printf("[UNIT-TEST:] -- Creating single trapezoid waveform.\n");
     __rectw_1__.fourier_series_expansion(1);
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     rectangular_waveform_fp16_t::create_signal_plot(__rectw_1__.__n_samples__,__rectw_1__.__rw_samples__.m_data,nullptr,
                                              "rectangular_wave_fp16_expanded_test_v2_1_","Rectangular_Waveform_FP16_Expansion",false);
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_rectangw_fp16_series_optimized();

void unit_test_rectangw_fp16_series_optimized()
{
     using namespace gms::radiolocation;
     constexpr std::size_t   n_samples{1024ull};
     constexpr std::uint32_t n_waves{32ull};
     constexpr float T{2.0f};
     constexpr float A{1.0f};
     constexpr float rho{1.0f};
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     rectangular_waveform_fp16_t __rectw_1__  = rectangular_waveform_fp16_t(n_samples,n_waves,
                                                      T,rho,A);
     char * ctor_name{gms::common::demangle(typeid(__rectw_1__).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(__rectw_1__).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
     std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     __rectw_1__.__rw_samples__.info_size_alignment();
     std::printf("[UNIT-TEST:] -- Creating single trapezoid waveform.\n");
     __rectw_1__.fourier_series_expansion_optim(1);
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     rectangular_waveform_fp16_t::create_signal_plot(__rectw_1__.__n_samples__,__rectw_1__.__rw_samples__.m_data,nullptr,
                                              "rectangular_wave_fp16_expanded_optim_test_v2_1_","Rectangular_Waveform_FP16_Expansion",false);
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}

int main()
{
     unit_test_rectangw_fp16_series_1();
     unit_test_rectangw_fp16_series_2();
     unit_test_rectangw_fp16_series_optimized();
     return 0;
}