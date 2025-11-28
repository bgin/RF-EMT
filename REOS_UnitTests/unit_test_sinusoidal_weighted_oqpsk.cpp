#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_dyn_array.h"
#include "GMS_sinusoidal_weighted_oqpsk.h"

/*
   icpc -o unit_test_sinusoid_weighted_oqpsk -O3 -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_sinusoidal_weighted_oqpsk.h GMS_sinusoidal_weighted_oqpsk.cpp unit_test_sinusoidal_weighted_oqpsk.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_sinusoidal_weighted_oqpsk.h GMS_sinusoidal_weighted_oqpsk.cpp unit_test_sinusoidal_weighted_oqpsk.cpp

*/

__attribute__((hot))
__attribute__((noinline))
std::int32_t 
unit_test_generate_I_channel_random_bitstream();

std::int32_t 
unit_test_generate_I_channel_random_bitstream()
{
      using namespace gms::radiolocation;
      using namespace gms;
      constexpr std::size_t nsamples{10ull*200ull}; // 10 periods of 200 samples each
      constexpr std::size_t I_nsamples{nsamples};
      constexpr std::size_t Q_nsamples{nsamples};
      constexpr float A_I{1.0f};
      constexpr float A_Q{1.0f};
      constexpr float T{10.0f};
      constexpr float sw0{9000.0f};
      constexpr float cw0{9000.0f};
      constexpr float cph0{0.78539816339744830961566084582f};
      constexpr float sph0{1.57079632679489661923132169164f};
      // same for both I/Q channels
      constexpr float I_duration{10.0f};
      constexpr float I_sample_rate{200.0f};
      constexpr float I_w0{10.0f};
      constexpr float I_ph0{0.78539816339744830961566084582f};
      int32_t status{};
      printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
      // For breaking into gdb
      //__asm__ ("int3");
      sinusoidal_weighted_oqpsk_t OQPSK_signal = sinusoidal_weighted_oqpsk_t(nsamples,I_nsamples,Q_nsamples,
                                                                             A_I,A_Q,T,cw0,sw0,cph0,sph0);
      char * ctor_name{gms::common::demangle(typeid(OQPSK_signal).name(),status)};
      if(status==0 && ctor_name != NULL)
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
      }
      else
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(OQPSK_signal).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
      }   

      std::printf("[UNIT-TEST]: -- Calling: I-channel user data random bitstream generator\n");
      const std::int32_t ret_status = OQPSK_signal.generate_I_channel_random_bitstream(I_duration,I_w0,I_ph0,I_sample_rate);
      if(ret_status<0)
      {
          std::printf("[UNIT-TEST]: -- generate_I_channel_random_bitstream -- **FAILED**, status=%d\n",ret_status);
          return (ret_status);
      }
      std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
      sinusoidal_weighted_oqpsk_t::create_signal_plot(static_cast<std::uint32_t>(OQPSK_signal.m_I_ch_nsamples),OQPSK_signal.m_I_bitstream.m_data,
                                                                                 nullptr,"generate_I_channel_random_bitstream_test_4_",
                                                                                 "I-Channel random bitstream",false);
      printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
      return (0);
}


int main()
{
    std::int32_t status_generate_I_channel_random_bitstream = unit_test_generate_I_channel_random_bitstream();
    if(status_generate_I_channel_random_bitstream==0)
    {
         //std::printf( ANSI_COLOR_GREEN "[UNIT-TEST-- generate_I_channel_random_bitstream]: PASSED!!" ANSI_RESET_ALL "\n");
         std::printf("[UNIT-TEST-- generate_I_channel_random_bitstream]: PASSED!! \n");
    }
    else 
    {
         //std::printf( ANSI_COLOR_RED "[UNIT-TEST-- generate_I_channel_random_bitstream]: FAILED!!" ANSI_RESET_ALL "\n");
         std::printf("[UNIT-TEST-- generate_I_channel_random_bitstream]: FAILED!! \n");
    }
    return 0;
}