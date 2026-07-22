/* Copyright (C) Bernard Gingold, 2020-2026 
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_dyn_array.h"
#include "GMS_trapezoid_waveform.h"

/*
   icpc -o unit_test_2_trapezw_single_v2 -fp-model -std=c++17 fast=2 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_trapezoid_waveform.h GMS_trapezoid_waveform.cpp unit_test_2_single_trapezoid_wave_v2.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_trapezoid_waveform.h GMS_trapezoid_waveform.cpp unit_test_2_single_trapezoid_wave_v2.cpp

*/

__attribute__((hot))
__attribute__((noinline))
void unit_test_2_create_trapezw_single_v2_1024();

void unit_test_2_create_trapezw_single_v2_1024()
{
     using namespace gms::radiolocation;
     using namespace gms;
     constexpr std::size_t   n_samples{1024ull};
     constexpr std::uint32_t n_waves{32ull};
     constexpr std::uint32_t n_a_param{32ull};
     constexpr std::uint32_t n_l_param{32ull};
     constexpr std::uint32_t n_c_param{32ull};
     constexpr std::uint32_t n_m_param{32ull};
     float a_params[n_a_param];
     float m_params[n_m_param];
     float l_params[n_m_param];
     float c_params[n_m_param];
     std::clock_t seed_a{};
     std::clock_t seed_m{};
     std::clock_t seed_l{};
     std::clock_t seed_c{};
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     seed_a = std::clock();
     auto rand_a{std::bind(std::uniform_real_distribution<float>(0.01f,0.1f),std::mt19937(seed_a))};
     std::fill_n(a_params,n_a_param,rand_a());
     seed_m = std::clock();
     auto rand_m{std::bind(std::uniform_real_distribution<float>(1.0f,2.0f),std::mt19937(seed_m))};
     std::fill_n(m_params,n_m_param,rand_m());
     seed_l = std::clock();
     auto rand_l{std::bind(std::uniform_real_distribution<float>(0.0f,0.9f),std::mt19937(seed_l))};
     std::fill_n(l_params,n_l_param,rand_l());
     seed_c = std::clock();
     auto rand_c{std::bind(std::uniform_real_distribution<float>(1.0f,2.0f),std::mt19937(seed_c))};
     std::fill_n(c_params,n_c_param,rand_c());
     trapezoid_waveform_t __trapezw_1__ = trapezoid_waveform_t(n_samples,n_waves,
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
     std::printf("[UNIT-TEST:] -- Creating series of trapezoid waves coded.\n");
     __trapezw_1__.single_trapezoid_wave_v2(&a_params[0],&m_params[0],&l_params[0],&c_params[0],2);
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     trapezoid_waveform_t::create_signal_plot(__trapezw_1__.__n_samples__,__trapezw_1__.__trapezw_samples__.m_data,nullptr,
                                              "single_trapezoid_wave_v2_1024_test_2","Series_Trapezoid_Waveform",false);
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_2_create_trapezw_single_v2_128();

void unit_test_2_create_trapezw_single_v2_128()
{
     using namespace gms::radiolocation;
     using namespace gms;
     constexpr std::size_t   n_samples{128ull};
     constexpr std::uint32_t n_waves{16ull};
     constexpr std::uint32_t n_a_param{32ull}; // not used here
     constexpr std::uint32_t n_l_param{32ull}; // same
     constexpr std::uint32_t n_c_param{32ull}; // same 
     constexpr std::uint32_t n_m_param{32ull}; // same 
     float a_params[n_a_param];
     float m_params[n_m_param];
     float l_params[n_m_param];
     float c_params[n_m_param];
     std::clock_t seed_a{};
     std::clock_t seed_m{};
     std::clock_t seed_l{};
     std::clock_t seed_c{};
          
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     seed_a = std::clock();
     auto rand_a{std::bind(std::uniform_real_distribution<float>(0.01f,0.1f),std::mt19937(seed_a))};
     std::fill_n(a_params,n_a_param,rand_a());
     seed_m = std::clock();
     auto rand_m{std::bind(std::uniform_real_distribution<float>(1.0f,2.0f),std::mt19937(seed_m))};
     std::fill_n(m_params,n_m_param,rand_m());
     seed_l = std::clock();
     auto rand_l{std::bind(std::uniform_real_distribution<float>(1.0f,2.0f),std::mt19937(seed_l))};
     std::fill_n(l_params,n_l_param,rand_l());
     seed_c = std::clock();
     auto rand_c{std::bind(std::uniform_real_distribution<float>(1.0f,2.0f),std::mt19937(seed_c))};
     std::fill_n(c_params,n_c_param,rand_c());
     trapezoid_waveform_t __trapezw_1__ = trapezoid_waveform_t(n_samples,n_waves,
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
     std::printf("[UNIT-TEST:] -- Creating series of trapezoid waves coded.\n");
     __trapezw_1__.single_trapezoid_wave_v2(&a_params[0],&m_params[0],&l_params[0],&c_params[0],2);
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     trapezoid_waveform_t::create_signal_plot(__trapezw_1__.__n_samples__,__trapezw_1__.__trapezw_samples__.m_data,nullptr,
                                              "single_trapezoid_wave_v2_128_test_2","Single_Trapezoid_Waveform",false);
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}


int main()
{
     std::clock_t seed = std::clock();
     int32_t which{-1};
     auto rand_d{std::bind(std::uniform_int_distribution<int32_t>(0,1),std::mt19937(seed))};
     which = rand_d();
     switch (which)
     {
         case 0 : 
            unit_test_2_create_trapezw_single_v2_128();
         break;
         case 1 : 
            unit_test_2_create_trapezw_single_v2_1024();
         break;
         default : 
            std::printf("[UNIT-TEST:] -- Invalid switch variable = %d\n",which);
            std::terminate();
     }
     
     
     return 0;
}
