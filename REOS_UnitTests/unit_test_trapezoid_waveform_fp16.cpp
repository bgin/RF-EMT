#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_dyn_array.h"
#include "GMS_trapezoid_waveform_fp16.h"

/*
   icpc -o unit_test_create_trapezw_fp16 -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_half.h GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_trapezoid_waveform_fp16.h GMS_trapezoid_waveform_fp16.cpp unit_test_trapezoid_waveform_fp16.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_half.h GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_trapezoid_waveform_fp16.h GMS_trapezoid_waveform_fp16.cpp unit_test_trapezoid_waveform_fp16.cpp

*/

__attribute__((hot))
__attribute__((noinline))
void unit_test_create_trapezw_fp16_coded_1024();

void unit_test_create_trapezw_fp16_coded_1024()
{
     using namespace gms::radiolocation;
     using namespace gms;
     using namespace half_float;
     constexpr std::size_t   n_samples{1024ull};
     constexpr std::uint32_t n_waves{32ull};
     constexpr std::uint32_t n_a_param{32ull};
     constexpr std::uint32_t n_l_param{32ull};
     constexpr std::uint32_t n_c_param{32ull};
     constexpr std::uint32_t n_m_param{32ull};
     constexpr float a{10.0f};
     constexpr float m{5.0f};
     constexpr float l{5.0f};
     constexpr float c{2.0f};
     float coded_seq[32];
     float sum;
     std::clock_t seed;
     
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     seed = std::clock();
     auto rand_x{std::bind(std::uniform_real_distribution<float>(-1.0f,1.0f),
                                                       std::mt19937(seed))}; 
     std::fill_n(coded_seq,32,rand_x());
     trapezoid_waveform_fp16_t __trapezw_1__ = trapezoid_waveform_fp16_t(n_samples,n_waves,
                                                               n_a_param,n_l_param,
                                                               n_c_param,n_m_param);
     char * ctor_name{gms::common::demangle(typeid(__trapezw_1__).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(__trapezw_1__).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
     std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     __trapezw_1__.__trapezw_samples__.info_size_alignment();
     std::printf("[UNIT-TEST:] -- Creating horizontal summation of trapezoid waves coded.\n");
     const float T{static_cast<float>(__trapezw_1__.__n_samples__)};
     for(std::size_t __t{0}; __t != __trapezw_1__.__n_samples__; ++__t) 
     {
         const float t{static_cast<float>(__t)};
         sum = 0.0f;
         for(std::uint32_t __k{0}; __k != __trapezw_1__.__n_waves__; ++__k)
         {
             const float k{static_cast<float>(__k)};
             sum += (__trapezw_1__.sample_of_trapezoid_wave(t-k*T,a,m,l,c)*coded_seq[__k]);
         }
         __trapezw_1__.__trapezw_samples__.m_data[__t] = half_cast<half>(sum);
     }
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     trapezoid_waveform_fp16_t::create_signal_plot(__trapezw_1__.__n_samples__,__trapezw_1__.__trapezw_samples__.m_data,nullptr,
                                              "series_of_trapezoid_waves_fp16_1024_test_1","HSummation_Trapezoid_Waveform",false);
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_create_trapezw_fp16_coded_128();

void unit_test_create_trapezw_fp16_coded_128()
{
     using namespace gms::radiolocation;
     using namespace gms;
     using namespace half_float;
     constexpr std::size_t   n_samples{128ull};
     constexpr std::uint32_t n_waves{8ull};
     constexpr std::uint32_t n_a_param{8ull};
     constexpr std::uint32_t n_l_param{8ull};
     constexpr std::uint32_t n_c_param{8ull};
     constexpr std::uint32_t n_m_param{8ull};
     constexpr float a{10.0f};
     constexpr float m{5.0f};
     constexpr float l{5.0f};
     constexpr float c{2.0f};
     float coded_seq[8];
     float sum;
     std::clock_t seed;
     
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     seed = std::clock();
     auto rand_x{std::bind(std::uniform_real_distribution<float>(-1.0f,1.0f),
                                                       std::mt19937(seed))}; 
     std::fill_n(coded_seq,8,rand_x());
     trapezoid_waveform_fp16_t __trapezw_1__ = trapezoid_waveform_fp16_t(n_samples,n_waves,
                                                               n_a_param,n_l_param,
                                                               n_c_param,n_m_param);
     char * ctor_name{gms::common::demangle(typeid(__trapezw_1__).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(__trapezw_1__).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
     std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     __trapezw_1__.__trapezw_samples__.info_size_alignment();
     std::printf("[UNIT-TEST:] -- Creating horizontal summation of trapezoid waves coded.\n");
     const float T{static_cast<float>(__trapezw_1__.__n_samples__)};
     for(std::size_t __t{0}; __t != __trapezw_1__.__n_samples__; ++__t) 
     {
         const float t{static_cast<float>(__t)};
         sum = 0.0f;
         for(std::uint32_t __k{0}; __k != __trapezw_1__.__n_waves__; ++__k)
         {
             const float k{static_cast<float>(__k)};
             sum += (__trapezw_1__.sample_of_trapezoid_wave(t-k*T,a,m,l,c)*coded_seq[__k]);
         }
         __trapezw_1__.__trapezw_samples__.m_data[__t] = half_cast<half>(sum);
     }
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     trapezoid_waveform_fp16_t::create_signal_plot(__trapezw_1__.__n_samples__,__trapezw_1__.__trapezw_samples__.m_data,nullptr,
                                              "series_of_trapezoid_waves_fp16_128_test_1","HSummation_Trapezoid_Waveform",false);
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_create_trapezw_fp16_coded_256();

void unit_test_create_trapezw_fp16_coded_256()
{
     using namespace gms::radiolocation;
     using namespace gms;
     using namespace half_float;
     constexpr std::size_t   n_samples{256ull};
     constexpr std::uint32_t n_waves{128ull};
     constexpr std::uint32_t n_a_param{8ull};
     constexpr std::uint32_t n_l_param{8ull};
     constexpr std::uint32_t n_c_param{8ull};
     constexpr std::uint32_t n_m_param{8ull};
     constexpr float a{1.0f};
     constexpr float m{2.0f};
     constexpr float l{1.0f};
     constexpr float c{0.5f};
     float coded_seq[128];
     float sum;
     std::clock_t seed;
     
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     seed = std::clock();
     auto rand_x{std::bind(std::uniform_real_distribution<float>(-1.0f,1.0f),
                                                       std::mt19937(seed))}; 
     std::fill_n(coded_seq,128,rand_x());
     trapezoid_waveform_fp16_t __trapezw_1__ = trapezoid_waveform_fp16_t(n_samples,n_waves,
                                                               n_a_param,n_l_param,
                                                               n_c_param,n_m_param);
     char * ctor_name{gms::common::demangle(typeid(__trapezw_1__).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(__trapezw_1__).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
     std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     __trapezw_1__.__trapezw_samples__.info_size_alignment();
     std::printf("[UNIT-TEST:] -- Creating horizontal summation of trapezoid waves coded.\n");
     const float T{static_cast<float>(__trapezw_1__.__n_samples__)};
     for(std::size_t __t{0}; __t != __trapezw_1__.__n_samples__; ++__t) 
     {
         const float t{static_cast<float>(__t)};
         sum = 0.0f;
         for(std::uint32_t __k{0}; __k != __trapezw_1__.__n_waves__; ++__k)
         {
             const float k{static_cast<float>(__k)};
             sum += (__trapezw_1__.sample_of_trapezoid_wave(t-k*T,a,m,l,c)*coded_seq[__k]);
         }
         __trapezw_1__.__trapezw_samples__.m_data[__t] = half_cast<half>(sum);
     }
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     trapezoid_waveform_fp16_t::create_signal_plot(__trapezw_1__.__n_samples__,__trapezw_1__.__trapezw_samples__.m_data,nullptr,
                                              "series_of_trapezoid_waves_fp16_256_test_1","HSummation_Trapezoid_Waveform",false);
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_create_trapezw_fp16_coded_32_8();

void unit_test_create_trapezw_fp16_coded_32_8()
{
     using namespace gms::radiolocation;
     using namespace gms;
     using namespace half_float;
     constexpr std::size_t   n_samples{32ull};
     constexpr std::uint32_t n_waves{8ull};
     constexpr std::uint32_t n_a_param{8ull};
     constexpr std::uint32_t n_l_param{8ull};
     constexpr std::uint32_t n_c_param{8ull};
     constexpr std::uint32_t n_m_param{8ull};
     constexpr float a{1.0f};
     constexpr float m{10.0f};
     constexpr float l{1.0f};
     constexpr float c{4.0f};
     float coded_seq[8];
     float sum;
     std::clock_t seed;
     
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     seed = std::clock();
     auto rand_x{std::bind(std::uniform_real_distribution<float>(-1.0f,1.0f),
                                                       std::mt19937(seed))}; 
     std::fill_n(coded_seq,8,rand_x());
     trapezoid_waveform_fp16_t __trapezw_1__ = trapezoid_waveform_fp16_t(n_samples,n_waves,
                                                               n_a_param,n_l_param,
                                                               n_c_param,n_m_param);
     char * ctor_name{gms::common::demangle(typeid(__trapezw_1__).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(__trapezw_1__).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
     std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     __trapezw_1__.__trapezw_samples__.info_size_alignment();
     std::printf("[UNIT-TEST:] -- Creating horizontal summation of trapezoid waves coded.\n");
     const float T{static_cast<float>(__trapezw_1__.__n_samples__)};
     for(std::size_t __t{0}; __t != __trapezw_1__.__n_samples__; ++__t) 
     {
         const float t{static_cast<float>(__t)};
         sum = 0.0f;
         for(std::uint32_t __k{0}; __k != __trapezw_1__.__n_waves__; ++__k)
         {
             const float k{static_cast<float>(__k)};
             sum += (__trapezw_1__.sample_of_trapezoid_wave(t-k*T,a,m,l,c)*coded_seq[__k]);
         }
         __trapezw_1__.__trapezw_samples__.m_data[__t] = half_cast<half>(sum);
     }
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     trapezoid_waveform_fp16_t::create_signal_plot(__trapezw_1__.__n_samples__,__trapezw_1__.__trapezw_samples__.m_data,nullptr,
                                              "series_of_trapezoid_waves_fp16_32_test_2","HSummation_Trapezoid_Waveform",false);
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_create_trapezw_fp16_coded_256_8_v2();

void unit_test_create_trapezw_fp16_coded_256_8_v2()
{
     using namespace gms::radiolocation;
     using namespace gms;
     using namespace half_float;
     constexpr std::size_t   n_samples{256ull};
     constexpr std::uint32_t n_waves{8ull};
     constexpr std::uint32_t n_a_param{8ull};
     constexpr std::uint32_t n_l_param{8ull};
     constexpr std::uint32_t n_c_param{8ull};
     constexpr std::uint32_t n_m_param{8ull};
     constexpr float a{10.0f};
     constexpr float m{10.0f};
     constexpr float l{1.0f};
     constexpr float c{4.0f};
     float coded_seq[8];
     float sum;
     std::clock_t seed;
     
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     seed = std::clock();
     auto rand_x{std::bind(std::uniform_real_distribution<float>(-1.0f,1.0f),
                                                       std::mt19937(seed))}; 
     std::fill_n(coded_seq,8,rand_x());
     trapezoid_waveform_fp16_t __trapezw_1__ = trapezoid_waveform_fp16_t(n_samples,n_waves,
                                                               n_a_param,n_l_param,
                                                               n_c_param,n_m_param);
     char * ctor_name{gms::common::demangle(typeid(__trapezw_1__).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(__trapezw_1__).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
     std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     __trapezw_1__.__trapezw_samples__.info_size_alignment();
     std::printf("[UNIT-TEST:] -- Creating horizontal summation of trapezoid waves coded.\n");
     const float T{static_cast<float>(__trapezw_1__.__n_samples__)};
     for(std::size_t __t{0}; __t != __trapezw_1__.__n_samples__; ++__t) 
     {
         const float t{static_cast<float>(__t)};
         sum = 0.0f;
         for(std::uint32_t __k{0}; __k != __trapezw_1__.__n_waves__; ++__k)
         {
             const float k{static_cast<float>(__k)};
             sum += (__trapezw_1__.sample_of_trapezoid_wave(t-k*T,a,m,l,c)*coded_seq[__k]);
         }
         __trapezw_1__.__trapezw_samples__.m_data[__t] = half_cast<half>(sum);
     }
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     trapezoid_waveform_fp16_t::create_signal_plot(__trapezw_1__.__n_samples__,__trapezw_1__.__trapezw_samples__.m_data,nullptr,
                                              "series_of_trapezoid_waves_fp16_256_test_2","HSummation_Trapezoid_Waveform",false);
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_create_trapezw_fp16_coded_256_16();

void unit_test_create_trapezw_fp16_coded_256_16()
{
     using namespace gms::radiolocation;
     using namespace gms;
     using namespace half_float;
     constexpr std::size_t   n_samples{256ull};
     constexpr std::uint32_t n_waves{8ull};
     constexpr std::uint32_t n_a_param{8ull};
     constexpr std::uint32_t n_l_param{8ull};
     constexpr std::uint32_t n_c_param{8ull};
     constexpr std::uint32_t n_m_param{8ull};
     constexpr float a{10.0f};
     constexpr float m{5.0f};
     constexpr float l{2.0f};
     constexpr float c{2.0f};
     float coded_seq[8];
     float sum;
     std::clock_t seed;
     
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     seed = std::clock();
     auto rand_x{std::bind(std::uniform_real_distribution<float>(1.0f,0.0f),
                                                       std::mt19937(seed))}; 
     std::fill_n(coded_seq,8,rand_x());
     trapezoid_waveform_fp16_t __trapezw_1__ = trapezoid_waveform_fp16_t(n_samples,n_waves,
                                                               n_a_param,n_l_param,
                                                               n_c_param,n_m_param);
     char * ctor_name{gms::common::demangle(typeid(__trapezw_1__).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(__trapezw_1__).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
     std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     __trapezw_1__.__trapezw_samples__.info_size_alignment();
     std::printf("[UNIT-TEST:] -- Creating horizontal summation of trapezoid waves coded.\n");
     const float T{static_cast<float>(__trapezw_1__.__n_samples__)};
     for(std::size_t __t{0}; __t != __trapezw_1__.__n_samples__; ++__t) 
     {
         const float t{static_cast<float>(__t)};
         sum = 0.0f;
         for(std::uint32_t __k{0}; __k != __trapezw_1__.__n_waves__; ++__k)
         {
             const float k{static_cast<float>(__k)};
             sum += (__trapezw_1__.sample_of_trapezoid_wave(t-k*T,a,m,l,c)*coded_seq[__k]);
         }
         __trapezw_1__.__trapezw_samples__.m_data[__t] = half_cast<half>(sum);
     }
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     trapezoid_waveform_fp16_t::create_signal_plot(__trapezw_1__.__n_samples__,__trapezw_1__.__trapezw_samples__.m_data,nullptr,
                                              "series_of_trapezoid_waves_fp16_256_test_9","HSummation_Trapezoid_Waveform",false);
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}


/*
   %equation
a = 10  %Amplitude
m = 5   %Time Period
l = 5   %Horizontal Spread
c = 2   %Vertical Spread
x = 0:.1:10 %Sample Points
Trapezoidal_Wave = a/pi*(asin(sin((pi/m)*x+l))+acos(cos((pi/m)*x+l)))-a/2+c;
*/


int main()
{
     unit_test_create_trapezw_fp16_coded_1024();
     unit_test_create_trapezw_fp16_coded_128();
     unit_test_create_trapezw_fp16_coded_256();
     unit_test_create_trapezw_fp16_coded_32_8();
     unit_test_create_trapezw_fp16_coded_256_8_v2();
     unit_test_create_trapezw_fp16_coded_256_16();
     return 0;
}