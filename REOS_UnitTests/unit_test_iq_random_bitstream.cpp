#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_dyn_array.h"
#include "GMS_IQ_random_bitstream.h"

/*
   icpc -o unit_test_iq_random_bitstream -O3 -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_IQ_random_bitstream.h GMS_IQ_random_bitstream.cpp unit_test_iq_random_bitstream.cpp
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_IQ_random_bitstream.h GMS_IQ_random_bitstream.cpp unit_test_iq_random_bitstream.cpp

*/

__attribute__((hot))
__attribute__((noinline))
void
unit_test_iq_random_bitstream();

void
unit_test_iq_random_bitstream()
{
      using namespace gms::radiolocation;
      using namespace gms;
      // same for both I/Q channels
      constexpr float I_duration{10.0f};
      constexpr float I_sample_rate{200.0f};
      constexpr float I_w0{10.0f};
      constexpr float I_ph0{0.78539816339744830961566084582f};
      constexpr float Q_duration{I_duration};
      constexpr float Q_sample_rate{I_sample_rate};
      constexpr float Q_w0{I_w0};
      constexpr float Q_ph0{I_ph0};
      int32_t status{};
      printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
      // For breaking into gdb
      //__asm__ ("int3");
      iq_random_bitstream_t IQ_rand_bitstream = iq_random_bitstream_t(I_duration,I_sample_rate,I_w0,I_ph0,
                                                                      Q_duration,Q_sample_rate,Q_w0,Q_ph0);
      char * ctor_name{gms::common::demangle(typeid(IQ_rand_bitstream).name(),status)};
      if(status==0 && ctor_name != NULL)
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
      }
      else
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(IQ_rand_bitstream).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
      }  

      {
          std::printf("[UNIT-TEST]: -- Calling: I-channel user data random bitstream generator\n");
          const std::int32_t ret_status = IQ_rand_bitstream.generate_channel_I_random_bitstream<10ull>();
          if(ret_status<0)
          {
             std::printf("[UNIT-TEST]: -- generate_I_channel_random_bitstream -- **FAILED**, status=%d\n",ret_status);
            
          }
          std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
          iq_random_bitstream_t::create_signal_plot(static_cast<std::uint32_t>(IQ_rand_bitstream.m_I_nsamples),IQ_rand_bitstream.m_I_bitstream.m_data,
                                                                                 nullptr,"generate_I_channel_random_bitstream_test_1_",
                                                                                 "I-Channel random bitstream",false);
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      {
         std::printf("[UNIT-TEST]: -- Calling: Q-channel user data random bitstream generator\n");
         const std::int32_t ret_status = IQ_rand_bitstream.generate_channel_Q_random_bitstream<10ull>();
         if(ret_status<0)
         {
             std::printf("[UNIT-TEST]: -- generate_Q_channel_random_bitstream -- **FAILED**, status=%d\n",ret_status);
             
         }
         std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
         iq_random_bitstream_t::create_signal_plot(static_cast<std::uint32_t>(IQ_rand_bitstream.m_Q_nsamples),IQ_rand_bitstream.m_Q_bitstream.m_data,
                                                                                 nullptr,"generate_Q_channel_random_bitstream_test_1_",
                                                                                 "Q-Channel random bitstream",false);
      }
      printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}


int main() 
{
    unit_test_iq_random_bitstream();
    return 0;
}