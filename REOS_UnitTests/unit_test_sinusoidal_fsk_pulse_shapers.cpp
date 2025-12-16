#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_dyn_array.h"
#include "GMS_sinusoidal_fsk.h"

/*
   icpc -o unit_test_sinusoidal_fsk_pulse_shapers -O3 -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_loop_indices_LUT.h GMS_loop_indices_LUT.cpp GMS_sinusoidal_fsk.h GMS_sinusoidal_fsk.cpp unit_test_sinusoidal_fsk_pulse_shapers.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_loop_indices_LUT.h GMS_loop_indices_LUT.cpp GMS_sinusoidal_fsk.h GMS_sinusoidal_fsk.cpp unit_test_sinusoidal_fsk_pulse_shapers.cpp

*/

__attribute__((hot))
__attribute__((noinline))
void 
unit_test_pulse_shaping_function_u8x();

void  
unit_test_pulse_shaping_function_u8x()
{
      using namespace gms::radiolocation;
      using namespace gms;
      constexpr float k{10.0f};
      constexpr float T{200.0f};
      constexpr std::size_t sfsk_nsamples{static_cast<std::size_t>(k*T)}; // 10 periods of 200 samples each
      constexpr std::size_t I_nsamples{sfsk_nsamples};
      constexpr std::size_t Q_nsamples{sfsk_nsamples};
      constexpr std::size_t psfunc_nsamples{sfsk_nsamples};
      
      constexpr float Ac{1.0f};
      constexpr float fc{9000.0f};
      constexpr float I_fc{9000.0f};
      constexpr float Q_fc{I_fc};
      constexpr float I_ph0{0.78539816339744830961566084582f};
      constexpr float Q_ph0{0.78539816339744830961566084582f};
      constexpr float ph0{0.0f};
      int32_t status{};
      printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
      // For breaking into gdb
      //__asm__ ("int3");
      sinusoidal_fsk_t SFSK_signal(sfsk_nsamples,I_nsamples,Q_nsamples,
                                   psfunc_nsamples,k,T,Ac,fc,ph0,
                                   I_fc,Q_fc,I_ph0,Q_ph0);
      char * ctor_name{gms::common::demangle(typeid(SFSK_signal).name(),status)};
      if(status==0 && ctor_name != NULL)
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
      }
      else
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(SFSK_signal).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
      }   

      
          std::printf("[UNIT-TEST]: -- Calling: generate_pulse_shaping_function_u8x\n");
          const std::int32_t ret_status = SFSK_signal.generate_pulse_shaping_function_u8x();
          if(ret_status<0)
          {
             std::printf("[UNIT-TEST]: -- generate_pulse_shaping_function_u8x -- **FAILED**, status=%d\n",ret_status);
             return;
          }
          std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
          sinusoidal_fsk_t::create_signal_plot(static_cast<std::uint32_t>(SFSK_signal.m_psfunc_nsamples),SFSK_signal.m_psfunc.m_data,
                                                                          nullptr,"generate_pulse_shaping_function_u8x_test_2_",
                                                                          "Pulse_Shaping_Function",false);   
         printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
     
}

__attribute__((hot))
__attribute__((noinline))
void 
unit_test_pulse_shaping_function_sse_u4x();

void  
unit_test_pulse_shaping_function_sse_u4x()
{
      using namespace gms::radiolocation;
      using namespace gms;
      constexpr float k{10.0f};
      constexpr float T{200.0f};
      constexpr std::size_t sfsk_nsamples{static_cast<std::size_t>(k*T)}; // 10 periods of 200 samples each
      constexpr std::size_t I_nsamples{sfsk_nsamples};
      constexpr std::size_t Q_nsamples{sfsk_nsamples};
      constexpr std::size_t psfunc_nsamples{sfsk_nsamples};
      
      constexpr float Ac{1.0f};
      constexpr float fc{9000.0f};
      constexpr float I_fc{9000.0f};
      constexpr float Q_fc{I_fc};
      constexpr float I_ph0{0.78539816339744830961566084582f};
      constexpr float Q_ph0{0.78539816339744830961566084582f};
      constexpr float ph0{0.0f};
      int32_t status{};
      printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
      // For breaking into gdb
      //__asm__ ("int3");
      sinusoidal_fsk_t SFSK_signal(sfsk_nsamples,I_nsamples,Q_nsamples,
                                   psfunc_nsamples,k,T,Ac,fc,ph0,
                                   I_fc,Q_fc,I_ph0,Q_ph0);
      char * ctor_name{gms::common::demangle(typeid(SFSK_signal).name(),status)};
      if(status==0 && ctor_name != NULL)
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
      }
      else
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(SFSK_signal).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
      }   

      
          std::printf("[UNIT-TEST]: -- Calling: generate_pulse_shaping_function_sse_u4x\n");
          const std::int32_t ret_status = SFSK_signal.generate_pulse_shaping_function_sse_u4x(true);
          if(ret_status<0)
          {
             std::printf("[UNIT-TEST]: -- generate_pulse_shaping_function_sse_u4x -- **FAILED**, status=%d\n",ret_status);
             return;
          }
          std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
          sinusoidal_fsk_t::create_signal_plot(static_cast<std::uint32_t>(SFSK_signal.m_psfunc_nsamples),SFSK_signal.m_psfunc.m_data,
                                                                          nullptr,"generate_pulse_shaping_function_sse_u4x_test_2_",
                                                                          "Pulse_Shaping_Function_SSE_v",false);   
         printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
     
}

__attribute__((hot))
__attribute__((noinline))
void 
unit_test_pulse_shaping_function_avx_u4x();

void  
unit_test_pulse_shaping_function_avx_u4x()
{
      using namespace gms::radiolocation;
      using namespace gms;
      constexpr float k{10.0f};
      constexpr float T{200.0f};
      constexpr std::size_t sfsk_nsamples{static_cast<std::size_t>(k*T)}; // 10 periods of 200 samples each
      constexpr std::size_t I_nsamples{sfsk_nsamples};
      constexpr std::size_t Q_nsamples{sfsk_nsamples};
      constexpr std::size_t psfunc_nsamples{sfsk_nsamples};
      
      constexpr float Ac{1.0f};
      constexpr float fc{9000.0f};
      constexpr float I_fc{9000.0f};
      constexpr float Q_fc{I_fc};
      constexpr float I_ph0{0.78539816339744830961566084582f};
      constexpr float Q_ph0{0.78539816339744830961566084582f};
      constexpr float ph0{0.0f};
      int32_t status{};
      printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
      // For breaking into gdb
      //__asm__ ("int3");
      sinusoidal_fsk_t SFSK_signal(sfsk_nsamples,I_nsamples,Q_nsamples,
                                   psfunc_nsamples,k,T,Ac,fc,ph0,
                                   I_fc,Q_fc,I_ph0,Q_ph0);
      char * ctor_name{gms::common::demangle(typeid(SFSK_signal).name(),status)};
      if(status==0 && ctor_name != NULL)
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
      }
      else
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(SFSK_signal).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
      }   

      
          std::printf("[UNIT-TEST]: -- Calling: generate_pulse_shaping_function_avx_u4x\n");
          const std::int32_t ret_status = SFSK_signal.generate_pulse_shaping_function_avx_u4x(true);
          if(ret_status<0)
          {
             std::printf("[UNIT-TEST]: -- generate_pulse_shaping_function_avx_u4x -- **FAILED**, status=%d\n",ret_status);
             return;
          }
          std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
          sinusoidal_fsk_t::create_signal_plot(static_cast<std::uint32_t>(SFSK_signal.m_psfunc_nsamples),SFSK_signal.m_psfunc.m_data,
                                                                          nullptr,"generate_pulse_shaping_function_avx_u4x_test_2_",
                                                                          "Pulse_Shaping_Function_AVX_v",false);   
         printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
     
}

__attribute__((hot))
__attribute__((noinline))
void 
unit_test_pulse_shaping_function_avx512_u4x();

void  
unit_test_pulse_shaping_function_avx512_u4x()
{
      using namespace gms::radiolocation;
      using namespace gms;
      constexpr float k{10.0f};
      constexpr float T{200.0f};
      constexpr std::size_t sfsk_nsamples{static_cast<std::size_t>(k*T)}; // 10 periods of 200 samples each
      constexpr std::size_t I_nsamples{sfsk_nsamples};
      constexpr std::size_t Q_nsamples{sfsk_nsamples};
      constexpr std::size_t psfunc_nsamples{sfsk_nsamples};
      
      constexpr float Ac{1.0f};
      constexpr float fc{9000.0f};
      constexpr float I_fc{9000.0f};
      constexpr float Q_fc{I_fc};
      constexpr float I_ph0{0.78539816339744830961566084582f};
      constexpr float Q_ph0{0.78539816339744830961566084582f};
      constexpr float ph0{0.0f};
      int32_t status{};
      printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
      // For breaking into gdb
      //__asm__ ("int3");
      sinusoidal_fsk_t SFSK_signal(sfsk_nsamples,I_nsamples,Q_nsamples,
                                   psfunc_nsamples,k,T,Ac,fc,ph0,
                                   I_fc,Q_fc,I_ph0,Q_ph0);
      char * ctor_name{gms::common::demangle(typeid(SFSK_signal).name(),status)};
      if(status==0 && ctor_name != NULL)
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
      }
      else
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(SFSK_signal).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
      }   

      
          std::printf("[UNIT-TEST]: -- Calling: generate_pulse_shaping_function_avx512_u4x\n");
          const std::int32_t ret_status = SFSK_signal.generate_pulse_shaping_function_avx512_u4x(true);
          if(ret_status<0)
          {
             std::printf("[UNIT-TEST]: -- generate_pulse_shaping_function_avx512_u4x -- **FAILED**, status=%d\n",ret_status);
             return;
          }
          std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
          sinusoidal_fsk_t::create_signal_plot(static_cast<std::uint32_t>(SFSK_signal.m_psfunc_nsamples),SFSK_signal.m_psfunc.m_data,
                                                                          nullptr,"generate_pulse_shaping_function_avx512_u4x_test_2_",
                                                                          "Pulse_Shaping_Function_AVX512_v",false);   
         printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
     
}



int main()
{
    //unit_test_pulse_shaping_function_u8x();
    unit_test_pulse_shaping_function_sse_u4x();
    unit_test_pulse_shaping_function_avx_u4x();
    unit_test_pulse_shaping_function_avx512_u4x();
    return 0;
}