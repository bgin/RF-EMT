#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_dyn_array.h"
#include "GMS_sinusoidal_weighted_oqpsk.h"

/*
   icpc -o unit_test_sinusoid_weighted_oqpsk -O3 -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_loop_indices_LUT.h GMS_loop_indices_LUT.cpp GMS_sinusoidal_weighted_oqpsk.h GMS_sinusoidal_weighted_oqpsk.cpp unit_test_sinusoidal_weighted_oqpsk.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_loop_indices_LUT.h GMS_loop_indices_LUT.cpp GMS_sinusoidal_weighted_oqpsk.h GMS_sinusoidal_weighted_oqpsk.cpp unit_test_sinusoidal_weighted_oqpsk.cpp

*/

__attribute__((hot))
__attribute__((noinline))
void 
unit_test_generate_I_Q_channel_bitstream();

void  
unit_test_generate_I_Q_channel_bitstream()
{
      using namespace gms::radiolocation;
      using namespace gms;
      constexpr std::size_t nsamples{10ull*200ull}; // 10 periods of 200 samples each
      constexpr std::size_t I_nsamples{nsamples};
      constexpr std::size_t Q_nsamples{nsamples};
      constexpr std::size_t nfreqs{120ull};
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
      constexpr float Q_duration{I_duration};
      constexpr float Q_sample_rate{I_sample_rate};
      constexpr float Q_w0{I_w0};
      constexpr float Q_ph0{I_ph0};
      int32_t status{};
      printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
      // For breaking into gdb
      //__asm__ ("int3");
      sinusoidal_weighted_oqpsk_t OQPSK_signal = sinusoidal_weighted_oqpsk_t(nsamples,I_nsamples,Q_nsamples,
                                                                             nfreqs,A_Q,T,cw0,sw0,cph0,sph0);
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

      {
          std::printf("[UNIT-TEST]: -- Calling: I-channel user data bitstream generator\n");
          const std::int32_t ret_status = OQPSK_signal.generate_I_channel_bitstream(I_duration,I_w0,I_ph0,I_sample_rate);
          if(ret_status<0)
          {
             std::printf("[UNIT-TEST]: -- generate_I_channel_bitstream -- **FAILED**, status=%d\n",ret_status);
             return;
          }
          std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
          sinusoidal_weighted_oqpsk_t::create_signal_plot(static_cast<std::uint32_t>(OQPSK_signal.m_I_ch_nsamples),OQPSK_signal.m_I_bitstream.m_data,
                                                                                 nullptr,"generate_I_channel_bitstream_test_5_",
                                                                                 "I-Channel bitstream",false);
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      {
         std::printf("[UNIT-TEST]: -- Calling: Q-channel user data bitstream generator\n");
         const std::int32_t ret_status = OQPSK_signal.generate_Q_channel_bitstream(Q_duration,Q_w0,Q_ph0,Q_sample_rate);
         if(ret_status<0)
         {
             std::printf("[UNIT-TEST]: -- generate_Q_channel_bitstream -- **FAILED**, status=%d\n",ret_status);
             return;
         }
         std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
         sinusoidal_weighted_oqpsk_t::create_signal_plot(static_cast<std::uint32_t>(OQPSK_signal.m_Q_ch_nsamples),OQPSK_signal.m_Q_bitstream.m_data,
                                                                                 nullptr,"generate_Q_channel_bitstream_test_1_",
                                                                                 "Q-Channel bitstream",false);
      }
      printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
     
}


__attribute__((hot))
__attribute__((noinline))
void 
unit_test_generate_I_Q_channel_bitstream_sse();

void  
unit_test_generate_I_Q_channel_bitstream_sse()
{
      using namespace gms::radiolocation;
      using namespace gms;
      constexpr std::size_t nsamples{10ull*200ull}; // 10 periods of 200 samples each
      constexpr std::size_t I_nsamples{nsamples};
      constexpr std::size_t Q_nsamples{nsamples};
      constexpr std::size_t nfreqs{120ull};
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
      constexpr float Q_duration{I_duration};
      constexpr float Q_sample_rate{I_sample_rate};
      constexpr float Q_w0{I_w0};
      constexpr float Q_ph0{I_ph0};
      int32_t status{};
      printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
      // For breaking into gdb
      //__asm__ ("int3");
      sinusoidal_weighted_oqpsk_t OQPSK_signal = sinusoidal_weighted_oqpsk_t(nsamples,I_nsamples,Q_nsamples,
                                                                             nfreqs,A_Q,T,cw0,sw0,cph0,sph0);
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

      {
          std::printf("[UNIT-TEST]: -- Calling: I-channel user data bitstream generator\n");
          const std::int32_t ret_status = OQPSK_signal.generate_I_channel_bitstream_sse(I_duration,I_w0,I_ph0,I_sample_rate);
          if(ret_status<0)
          {
             std::printf("[UNIT-TEST]: -- generate_I_channel_bitstream_sse -- **FAILED**, status=%d\n",ret_status);
             return;
          }
          std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
          sinusoidal_weighted_oqpsk_t::create_signal_plot(static_cast<std::uint32_t>(OQPSK_signal.m_I_ch_nsamples),OQPSK_signal.m_I_bitstream.m_data,
                                                                                 nullptr,"generate_I_channel_bitstream_sse_test_1_",
                                                                                 "I-Channel rectwave bitstream",false);
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      {
         std::printf("[UNIT-TEST]: -- Calling: Q-channel user data  bitstream generator\n");
         const std::int32_t ret_status = OQPSK_signal.generate_Q_channel_bitstream_sse(Q_duration,Q_w0,Q_ph0,Q_sample_rate);
         if(ret_status<0)
         {
             std::printf("[UNIT-TEST]: -- generate_Q_channel_bitstream_sse -- **FAILED**, status=%d\n",ret_status);
             return;
         }
         std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
         sinusoidal_weighted_oqpsk_t::create_signal_plot(static_cast<std::uint32_t>(OQPSK_signal.m_Q_ch_nsamples),OQPSK_signal.m_Q_bitstream.m_data,
                                                                                 nullptr,"generate_Q_channel_bitstream_sse_test_1_",
                                                                                 "Q-Channel rectwave bitstream",false);
      }
      printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
     
}

__attribute__((hot))
__attribute__((noinline))
void 
unit_test_generate_I_Q_channel_bitstream_avx();

void  
unit_test_generate_I_Q_channel_bitstream_avx()
{
      using namespace gms::radiolocation;
      using namespace gms;
      constexpr std::size_t nsamples{10ull*200ull}; // 10 periods of 200 samples each
      constexpr std::size_t I_nsamples{nsamples};
      constexpr std::size_t Q_nsamples{nsamples};
      constexpr std::size_t nfreqs{120ull};
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
      constexpr float Q_duration{I_duration};
      constexpr float Q_sample_rate{I_sample_rate};
      constexpr float Q_w0{I_w0};
      constexpr float Q_ph0{I_ph0};
      int32_t status{};
      printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
      // For breaking into gdb
      //__asm__ ("int3");
      sinusoidal_weighted_oqpsk_t OQPSK_signal = sinusoidal_weighted_oqpsk_t(nsamples,I_nsamples,Q_nsamples,
                                                                             nfreqs,A_I,A_Q,T,cw0,sw0,cph0,sph0);
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

      {
          std::printf("[UNIT-TEST]: -- Calling: generate_I_channel_bitstream_avx\n");
          const std::int32_t ret_status = OQPSK_signal.generate_I_channel_bitstream_avx(I_duration,I_w0,I_ph0,I_sample_rate);
          if(ret_status<0)
          {
             std::printf("[UNIT-TEST]: -- generate_I_channel_bitstream_avx -- **FAILED**, status=%d\n",ret_status);
             return;
          }
          std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
          sinusoidal_weighted_oqpsk_t::create_signal_plot(static_cast<std::uint32_t>(OQPSK_signal.m_I_ch_nsamples),OQPSK_signal.m_I_bitstream.m_data,
                                                                                 nullptr,"generate_I_channel_bitstream_avx_test_1_",
                                                                                 "I-Channel rectwave bitstream",false);
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      {
         std::printf("[UNIT-TEST]: -- Calling: generate_Q_channel_bitstream_avx\n");
         const std::int32_t ret_status = OQPSK_signal.generate_Q_channel_bitstream_avx(Q_duration,Q_w0,Q_ph0,Q_sample_rate);
         if(ret_status<0)
         {
             std::printf("[UNIT-TEST]: -- generate_Q_channel_bitstream_avx -- **FAILED**, status=%d\n",ret_status);
             return;
         }
         std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
         sinusoidal_weighted_oqpsk_t::create_signal_plot(static_cast<std::uint32_t>(OQPSK_signal.m_Q_ch_nsamples),OQPSK_signal.m_Q_bitstream.m_data,
                                                                                 nullptr,"generate_Q_channel_bitstream_avx_test_1_",
                                                                                 "Q-Channel rectwave bitstream",false);
      }
      printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
     
}


__attribute__((hot))
__attribute__((noinline))
void 
unit_test_generate_I_Q_channel_bitstream_sse_u4x();

void  
unit_test_generate_I_Q_channel_bitstream_sse_u4x()
{
      using namespace gms::radiolocation;
      using namespace gms;
      //constexpr std::size_t nsamples{10ull*200ull}; // 10 periods of 200 samples each
      constexpr std::size_t nsamples{10ull*300ull}; // exceeds the LUT indices storage.
      constexpr std::size_t I_nsamples{nsamples};
      constexpr std::size_t Q_nsamples{nsamples};
      constexpr std::size_t nfreqs{120ull};
      constexpr float A_I{1.0f};
      constexpr float A_Q{1.0f};
      constexpr float T{10.0f};
      constexpr float sw0{9000.0f};
      constexpr float cw0{9000.0f};
      constexpr float cph0{0.78539816339744830961566084582f};
      constexpr float sph0{1.57079632679489661923132169164f};
      // same for both I/Q channels
      constexpr float I_duration{10.0f};
      constexpr float I_sample_rate{300.0f};
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
      sinusoidal_weighted_oqpsk_t OQPSK_signal = sinusoidal_weighted_oqpsk_t(nsamples,I_nsamples,Q_nsamples,
                                                                             nfreqs,A_I,A_Q,T,cw0,sw0,cph0,sph0);
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

      {
          std::printf("[UNIT-TEST]: -- Calling: generate_I_channel_bitstream_sse_u4x.\n");
          const std::int32_t ret_status = OQPSK_signal.generate_I_channel_bitstream_sse_u4x(I_duration,I_w0,I_ph0,I_sample_rate);
          if(ret_status<0)
          {
             std::printf("[UNIT-TEST]: -- generate_I_channel_bitstream_sse_u4x -- **FAILED**, status=%d\n",ret_status);
             return;
          }
          std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
          sinusoidal_weighted_oqpsk_t::create_signal_plot(static_cast<std::uint32_t>(OQPSK_signal.m_I_ch_nsamples),OQPSK_signal.m_I_bitstream.m_data,
                                                                                 nullptr,"generate_I_channel_bitstream_sse_u4x_test_2_",
                                                                                 "I-Channel rectwave bitstream",false);
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      {
         std::printf("[UNIT-TEST]: -- Calling: generate_Q_channel_bitstream_sse_u4x\n");
         const std::int32_t ret_status = OQPSK_signal.generate_Q_channel_bitstream_sse_u4x(Q_duration,Q_w0,Q_ph0,Q_sample_rate);
         if(ret_status<0)
         {
             std::printf("[UNIT-TEST]: -- generate_Q_channel_bitstream_sse_u4x -- **FAILED**, status=%d\n",ret_status);
             return;
         }
         std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
         sinusoidal_weighted_oqpsk_t::create_signal_plot(static_cast<std::uint32_t>(OQPSK_signal.m_Q_ch_nsamples),OQPSK_signal.m_Q_bitstream.m_data,
                                                                                 nullptr,"generate_Q_channel_bitstream_sse_u4x_test_2_",
                                                                                 "Q-Channel rectwave bitstream",false);
      }
      printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
     
}


int main()
{
    unit_test_generate_I_Q_channel_random_bitstream();
    unit_test_generate_I_Q_channel_bitstream_sse();
    unit_test_generate_I_Q_channel_bitstream_sse_u4x();
    unit_test_generate_I_Q_channel_bitstream_avx();
    
    return 0;
}