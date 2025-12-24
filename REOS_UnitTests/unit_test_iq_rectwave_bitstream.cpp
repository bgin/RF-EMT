#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_dyn_array.h"
#include "GMS_iq_rectwave_bitstream.h"

/*
   icpc -o unit_test_rectwave_bitstream -O3 -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_loop_indices_LUT.h GMS_loop_indices_LUT.cpp GMS_iq_rectwave_bitstream.h GMS_iq_rectwave_bitstream.cpp unit_test_rectwave_bitstream.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_loop_indices_LUT.h GMS_loop_indices_LUT.cpp GMS_iq_rectwave_bitstream.h GMS_iq_rectwave_bitstream.cpp unit_test_rectwave_bitstream.cpp

*/

__attribute__((hot))
__attribute__((noinline))
void 
unit_test_generate_I_Q_rectwave_bitstream_scalar();

void  
unit_test_generate_I_Q_rectwave_bitstream_scalar()
{
      using namespace gms::radiolocation;
      using namespace gms;
      constexpr float I_k{10.0f};
      constexpr float I_T{200.0f};
      constexpr float I_Ac{1.54f};
      constexpr float I_fc{9500.0f};
      constexpr float I_ph0{0.78539816339744830961566084582f};
      constexpr float Q_k{40.0f};
      constexpr float Q_T{50.0f};
      constexpr float Q_Ac{4.54f};
      constexpr float Q_fc{10500.0f};
      constexpr float Q_ph0{-0.78539816339744830961566084582f};
      int32_t status{};
      printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
      // For breaking into gdb
      //__asm__ ("int3");
      iq_rectwave_bitstream_t IQ_rectw_bitstream(I_k,I_T,I_Ac,I_fc,I_ph0,
                                                 Q_k,Q_T,Q_Ac,Q_fc,Q_ph0);
      char * ctor_name{gms::common::demangle(typeid(IQ_rectw_bitstream).name(),status)};
      if(status==0 && ctor_name != NULL)
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
      }
      else
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(IQ_rectw_bitstream).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
      }   

      {
          std::printf("[UNIT-TEST]: -- Calling: I-channel user data bitstream generator\n");
          const std::int32_t ret_status = IQ_rectw_bitstream.generate_I_channel_bitstream_scalar();
          if(ret_status<0)
          {
             std::printf("[UNIT-TEST]: -- generate_I_channel_bitstream_scalar -- **FAILED**, status=%d\n",ret_status);
             return;
          }
          std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
          iq_rectwave_bitstream_t::create_signal_plot(static_cast<std::uint32_t>(IQ_rectw_bitstream.m_I_nsamples),IQ_rectw_bitstream.m_I_ch_bitstream.m_data,
                                                                                 nullptr,"generate_I_channel_bitstream_scalar_test_1_",
                                                                                 "I-Channel rectwave bitstream",false);
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      {
          std::printf("[UNIT-TEST]: -- Calling: Q-channel user data bitstream generator\n");
          const std::int32_t ret_status = IQ_rectw_bitstream.generate_Q_channel_bitstream_scalar();
          if(ret_status<0)
          {
             std::printf("[UNIT-TEST]: -- generate_Q_channel_bitstream_scalar -- **FAILED**, status=%d\n",ret_status);
             return;
          }
          std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
          iq_rectwave_bitstream_t::create_signal_plot(static_cast<std::uint32_t>(IQ_rectw_bitstream.m_Q_nsamples),IQ_rectw_bitstream.m_Q_ch_bitstream.m_data,
                                                                                 nullptr,"generate_Q_channel_bitstream_scalar_test_1_",
                                                                                 "Q-Channel rectwave bitstream",false);
      }
      printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
     
}

__attribute__((hot))
__attribute__((noinline))
void 
unit_test_generate_I_Q_rectwave_bitstream_sse();

void  
unit_test_generate_I_Q_rectwave_bitstream_sse()
{
      using namespace gms::radiolocation;
      using namespace gms;
      constexpr float I_k{10.0f};
      constexpr float I_T{200.0f};
      constexpr float I_Ac{1.54f};
      constexpr float I_fc{9500.0f};
      constexpr float I_ph0{0.78539816339744830961566084582f};
      constexpr float Q_k{40.0f};
      constexpr float Q_T{50.0f};
      constexpr float Q_Ac{4.54f};
      constexpr float Q_fc{10500.0f};
      constexpr float Q_ph0{-0.78539816339744830961566084582f};
      int32_t status{};
      printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
      // For breaking into gdb
      //__asm__ ("int3");
      iq_rectwave_bitstream_t IQ_rectw_bitstream(I_k,I_T,I_Ac,I_fc,I_ph0,
                                                 Q_k,Q_T,Q_Ac,Q_fc,Q_ph0);
      char * ctor_name{gms::common::demangle(typeid(IQ_rectw_bitstream).name(),status)};
      if(status==0 && ctor_name != NULL)
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
      }
      else
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(IQ_rectw_bitstream).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
      }   

      {
          std::printf("[UNIT-TEST]: -- Calling: generate_I_channel_bitstream_sse.\n");
          const std::int32_t ret_status = IQ_rectw_bitstream.generate_I_channel_bitstream_sse();
          if(ret_status<0)
          {
             std::printf("[UNIT-TEST]: -- generate_I_channel_bitstream_sse -- **FAILED**, status=%d\n",ret_status);
             return;
          }
          std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
          iq_rectwave_bitstream_t::create_signal_plot(static_cast<std::uint32_t>(IQ_rectw_bitstream.m_I_nsamples),IQ_rectw_bitstream.m_I_ch_bitstream.m_data,
                                                                                 nullptr,"generate_I_channel_bitstream_sse_test_1_",
                                                                                 "I-Channel rectwave bitstream",false);
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      {
          std::printf("[UNIT-TEST]: -- Calling: generate_Q_channel_bitstream_sse.\n");
          const std::int32_t ret_status = IQ_rectw_bitstream.generate_Q_channel_bitstream_sse();
          if(ret_status<0)
          {
             std::printf("[UNIT-TEST]: -- generate_Q_channel_bitstream_sse -- **FAILED**, status=%d\n",ret_status);
             return;
          }
          std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
          iq_rectwave_bitstream_t::create_signal_plot(static_cast<std::uint32_t>(IQ_rectw_bitstream.m_Q_nsamples),IQ_rectw_bitstream.m_Q_ch_bitstream.m_data,
                                                                                 nullptr,"generate_Q_channel_bitstream_sse_test_1_",
                                                                                 "Q-Channel rectwave bitstream",false);
      }
      printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
     
}

__attribute__((hot))
__attribute__((noinline))
void 
unit_test_generate_I_Q_rectwave_bitstream_sse_u4x();

void  
unit_test_generate_I_Q_rectwave_bitstream_sse_u4x()
{
      using namespace gms::radiolocation;
      using namespace gms;
      constexpr float I_k{10.0f};
      constexpr float I_T{200.0f};
      constexpr float I_Ac{1.54f};
      constexpr float I_fc{9500.0f};
      constexpr float I_ph0{0.78539816339744830961566084582f};
      constexpr float Q_k{40.0f};
      constexpr float Q_T{50.0f};
      constexpr float Q_Ac{4.54f};
      constexpr float Q_fc{10500.0f};
      constexpr float Q_ph0{-0.78539816339744830961566084582f};
      int32_t status{};
      printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
      // For breaking into gdb
      //__asm__ ("int3");
      iq_rectwave_bitstream_t IQ_rectw_bitstream(I_k,I_T,I_Ac,I_fc,I_ph0,
                                                 Q_k,Q_T,Q_Ac,Q_fc,Q_ph0);
      char * ctor_name{gms::common::demangle(typeid(IQ_rectw_bitstream).name(),status)};
      if(status==0 && ctor_name != NULL)
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
      }
      else
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(IQ_rectw_bitstream).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
      }   

      {
          std::printf("[UNIT-TEST]: -- Calling: generate_I_channel_bitstream_sse_u4x.\n");
          const std::int32_t ret_status = IQ_rectw_bitstream.generate_I_channel_bitstream_sse_u4x();
          if(ret_status<0)
          {
             std::printf("[UNIT-TEST]: -- generate_I_channel_bitstream_sse_u4x -- **FAILED**, status=%d\n",ret_status);
             return;
          }
          std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
          iq_rectwave_bitstream_t::create_signal_plot(static_cast<std::uint32_t>(IQ_rectw_bitstream.m_I_nsamples),IQ_rectw_bitstream.m_I_ch_bitstream.m_data,
                                                                                 nullptr,"generate_I_channel_bitstream_sse_u4x_test_1_",
                                                                                 "I-Channel rectwave bitstream",false);
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      {
          std::printf("[UNIT-TEST]: -- Calling: generate_Q_channel_bitstream_sse_u4x.\n");
          const std::int32_t ret_status = IQ_rectw_bitstream.generate_Q_channel_bitstream_sse_u4x();
          if(ret_status<0)
          {
             std::printf("[UNIT-TEST]: -- generate_Q_channel_bitstream_sse_u4x -- **FAILED**, status=%d\n",ret_status);
             return;
          }
          std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
          iq_rectwave_bitstream_t::create_signal_plot(static_cast<std::uint32_t>(IQ_rectw_bitstream.m_Q_nsamples),IQ_rectw_bitstream.m_Q_ch_bitstream.m_data,
                                                                                 nullptr,"generate_Q_channel_bitstream_sse_u4x_test_1_",
                                                                                 "Q-Channel rectwave bitstream",false);
      }
      printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
     
}

__attribute__((hot))
__attribute__((noinline))
void 
unit_test_generate_I_Q_rectwave_bitstream_avx();

void  
unit_test_generate_I_Q_rectwave_bitstream_avx()
{
      using namespace gms::radiolocation;
      using namespace gms;
      constexpr float I_k{10.0f};
      constexpr float I_T{200.0f};
      constexpr float I_Ac{1.54f};
      constexpr float I_fc{9500.0f};
      constexpr float I_ph0{0.78539816339744830961566084582f};
      constexpr float Q_k{40.0f};
      constexpr float Q_T{50.0f};
      constexpr float Q_Ac{4.54f};
      constexpr float Q_fc{10500.0f};
      constexpr float Q_ph0{-0.78539816339744830961566084582f};
      int32_t status{};
      printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
      // For breaking into gdb
      //__asm__ ("int3");
      iq_rectwave_bitstream_t IQ_rectw_bitstream(I_k,I_T,I_Ac,I_fc,I_ph0,
                                                 Q_k,Q_T,Q_Ac,Q_fc,Q_ph0);
      char * ctor_name{gms::common::demangle(typeid(IQ_rectw_bitstream).name(),status)};
      if(status==0 && ctor_name != NULL)
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
      }
      else
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(IQ_rectw_bitstream).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
      }   

      {
          std::printf("[UNIT-TEST]: -- Calling: generate_I_channel_bitstream_avx.\n");
          const std::int32_t ret_status = IQ_rectw_bitstream.generate_I_channel_bitstream_avx(true);
          if(ret_status<0)
          {
             std::printf("[UNIT-TEST]: -- generate_I_channel_bitstream_avx -- **FAILED**, status=%d\n",ret_status);
             return;
          }
          std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
          iq_rectwave_bitstream_t::create_signal_plot(static_cast<std::uint32_t>(IQ_rectw_bitstream.m_I_nsamples),IQ_rectw_bitstream.m_I_ch_bitstream.m_data,
                                                                                 nullptr,"generate_I_channel_bitstream_avx_test_1_",
                                                                                 "I-Channel rectwave bitstream",false);
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      {
          std::printf("[UNIT-TEST]: -- Calling: generate_Q_channel_bitstream_avx.\n");
          const std::int32_t ret_status = IQ_rectw_bitstream.generate_Q_channel_bitstream_avx(true);
          if(ret_status<0)
          {
             std::printf("[UNIT-TEST]: -- generate_Q_channel_bitstream_avx -- **FAILED**, status=%d\n",ret_status);
             return;
          }
          std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
          iq_rectwave_bitstream_t::create_signal_plot(static_cast<std::uint32_t>(IQ_rectw_bitstream.m_Q_nsamples),IQ_rectw_bitstream.m_Q_ch_bitstream.m_data,
                                                                                 nullptr,"generate_Q_channel_bitstream_avx_test_1_",
                                                                                 "Q-Channel rectwave bitstream",false);
      }
      printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
     
}

__attribute__((hot))
__attribute__((noinline))
void 
unit_test_generate_I_Q_rectwave_bitstream_avx_u4x();

void  
unit_test_generate_I_Q_rectwave_bitstream_avx_u4x()
{
      using namespace gms::radiolocation;
      using namespace gms;
      constexpr float I_k{10.0f};
      constexpr float I_T{200.0f};
      constexpr float I_Ac{1.54f};
      constexpr float I_fc{9500.0f};
      constexpr float I_ph0{0.78539816339744830961566084582f};
      constexpr float Q_k{40.0f};
      constexpr float Q_T{50.0f};
      constexpr float Q_Ac{4.54f};
      constexpr float Q_fc{10500.0f};
      constexpr float Q_ph0{-0.78539816339744830961566084582f};
      int32_t status{};
      printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
      // For breaking into gdb
      //__asm__ ("int3");
      iq_rectwave_bitstream_t IQ_rectw_bitstream(I_k,I_T,I_Ac,I_fc,I_ph0,
                                                 Q_k,Q_T,Q_Ac,Q_fc,Q_ph0);
      char * ctor_name{gms::common::demangle(typeid(IQ_rectw_bitstream).name(),status)};
      if(status==0 && ctor_name != NULL)
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
      }
      else
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(IQ_rectw_bitstream).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
      }   

      {
          std::printf("[UNIT-TEST]: -- Calling: generate_I_channel_bitstream_avx_u4x.\n");
          const std::int32_t ret_status = IQ_rectw_bitstream.generate_I_channel_bitstream_avx_u4x(true);
          if(ret_status<0)
          {
             std::printf("[UNIT-TEST]: -- generate_I_channel_bitstream_avx_u4x -- **FAILED**, status=%d\n",ret_status);
             return;
          }
          std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
          iq_rectwave_bitstream_t::create_signal_plot(static_cast<std::uint32_t>(IQ_rectw_bitstream.m_I_nsamples),IQ_rectw_bitstream.m_I_ch_bitstream.m_data,
                                                                                 nullptr,"generate_I_channel_bitstream_avx_u4x_test_1_",
                                                                                 "I-Channel rectwave bitstream",false);
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      {
          std::printf("[UNIT-TEST]: -- Calling: generate_Q_channel_bitstream_avx_u4x.\n");
          const std::int32_t ret_status = IQ_rectw_bitstream.generate_Q_channel_bitstream_avx_u4x(true);
          if(ret_status<0)
          {
             std::printf("[UNIT-TEST]: -- generate_Q_channel_bitstream_avx_u4x -- **FAILED**, status=%d\n",ret_status);
             return;
          }
          std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
          iq_rectwave_bitstream_t::create_signal_plot(static_cast<std::uint32_t>(IQ_rectw_bitstream.m_Q_nsamples),IQ_rectw_bitstream.m_Q_ch_bitstream.m_data,
                                                                                 nullptr,"generate_Q_channel_bitstream_avx_u4x_test_1_",
                                                                                 "Q-Channel rectwave bitstream",false);
      }
      printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
     
}

__attribute__((hot))
__attribute__((noinline))
void 
unit_test_generate_I_Q_rectwave_bitstream_avx512();

void  
unit_test_generate_I_Q_rectwave_bitstream_avx512()
{
      using namespace gms::radiolocation;
      using namespace gms;
      constexpr float I_k{10.0f};
      constexpr float I_T{200.0f};
      constexpr float I_Ac{1.54f};
      constexpr float I_fc{9500.0f};
      constexpr float I_ph0{0.78539816339744830961566084582f};
      constexpr float Q_k{40.0f};
      constexpr float Q_T{50.0f};
      constexpr float Q_Ac{4.54f};
      constexpr float Q_fc{10500.0f};
      constexpr float Q_ph0{-0.78539816339744830961566084582f};
      int32_t status{};
      printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
      // For breaking into gdb
      //__asm__ ("int3");
      iq_rectwave_bitstream_t IQ_rectw_bitstream(I_k,I_T,I_Ac,I_fc,I_ph0,
                                                 Q_k,Q_T,Q_Ac,Q_fc,Q_ph0);
      char * ctor_name{gms::common::demangle(typeid(IQ_rectw_bitstream).name(),status)};
      if(status==0 && ctor_name != NULL)
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
      }
      else
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(IQ_rectw_bitstream).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
      }   

      {
          std::printf("[UNIT-TEST]: -- Calling: generate_I_channel_bitstream_avx512.\n");
          const std::int32_t ret_status = IQ_rectw_bitstream.generate_I_channel_bitstream_avx512(true);
          if(ret_status<0)
          {
             std::printf("[UNIT-TEST]: -- generate_I_channel_bitstream_avx512 -- **FAILED**, status=%d\n",ret_status);
             return;
          }
          std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
          iq_rectwave_bitstream_t::create_signal_plot(static_cast<std::uint32_t>(IQ_rectw_bitstream.m_I_nsamples),IQ_rectw_bitstream.m_I_ch_bitstream.m_data,
                                                                                 nullptr,"generate_I_channel_bitstream_avx512_test_1_",
                                                                                 "I-Channel rectwave bitstream",false);
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      {
          std::printf("[UNIT-TEST]: -- Calling: generate_Q_channel_bitstream_avx512.\n");
          const std::int32_t ret_status = IQ_rectw_bitstream.generate_Q_channel_bitstream_avx512(true);
          if(ret_status<0)
          {
             std::printf("[UNIT-TEST]: -- generate_Q_channel_bitstream_avx512 -- **FAILED**, status=%d\n",ret_status);
             return;
          }
          std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
          iq_rectwave_bitstream_t::create_signal_plot(static_cast<std::uint32_t>(IQ_rectw_bitstream.m_Q_nsamples),IQ_rectw_bitstream.m_Q_ch_bitstream.m_data,
                                                                                 nullptr,"generate_Q_channel_bitstream_avx512_test_1_",
                                                                                 "Q-Channel rectwave bitstream",false);
      }
      printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
     
}

__attribute__((hot))
__attribute__((noinline))
void 
unit_test_generate_I_Q_rectwave_bitstream_avx512_u4x();

void  
unit_test_generate_I_Q_rectwave_bitstream_avx512_u4x()
{
      using namespace gms::radiolocation;
      using namespace gms;
      constexpr float I_k{10.0f};
      constexpr float I_T{200.0f};
      constexpr float I_Ac{1.54f};
      constexpr float I_fc{95.0f};
      constexpr float I_ph0{0.78539816339744830961566084582f};
      constexpr float Q_k{40.0f};
      constexpr float Q_T{50.0f};
      constexpr float Q_Ac{4.54f};
      constexpr float Q_fc{105.0f};
      constexpr float Q_ph0{-0.78539816339744830961566084582f};
      int32_t status{};
      printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
      // For breaking into gdb
      //__asm__ ("int3");
      iq_rectwave_bitstream_t IQ_rectw_bitstream(I_k,I_T,I_Ac,I_fc,I_ph0,
                                                 Q_k,Q_T,Q_Ac,Q_fc,Q_ph0);
      char * ctor_name{gms::common::demangle(typeid(IQ_rectw_bitstream).name(),status)};
      if(status==0 && ctor_name != NULL)
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
      }
      else
      {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(IQ_rectw_bitstream).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
      }   

      {
          std::printf("[UNIT-TEST]: -- Calling: generate_I_channel_bitstream_avx512_u4x.\n");
          const std::int32_t ret_status = IQ_rectw_bitstream.generate_I_channel_bitstream_avx512_u4x(true);
          if(ret_status<0)
          {
             std::printf("[UNIT-TEST]: -- generate_I_channel_bitstream_avx512_u4x -- **FAILED**, status=%d\n",ret_status);
             return;
          }
          std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
          iq_rectwave_bitstream_t::create_signal_plot(static_cast<std::uint32_t>(IQ_rectw_bitstream.m_I_nsamples),IQ_rectw_bitstream.m_I_ch_bitstream.m_data,
                                                                                 nullptr,"generate_I_channel_bitstream_avx512_u4x_test_2_",
                                                                                 "I-Channel rectwave bitstream",false);
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      {
          std::printf("[UNIT-TEST]: -- Calling: generate_Q_channel_bitstream_avx512_u4x.\n");
          const std::int32_t ret_status = IQ_rectw_bitstream.generate_Q_channel_bitstream_avx512_u4x(true);
          if(ret_status<0)
          {
             std::printf("[UNIT-TEST]: -- generate_Q_channel_bitstream_avx512_u4x -- **FAILED**, status=%d\n",ret_status);
             return;
          }
          std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
          iq_rectwave_bitstream_t::create_signal_plot(static_cast<std::uint32_t>(IQ_rectw_bitstream.m_Q_nsamples),IQ_rectw_bitstream.m_Q_ch_bitstream.m_data,
                                                                                 nullptr,"generate_Q_channel_bitstream_avx512_u4x_test_2_",
                                                                                 "Q-Channel rectwave bitstream",false);
      }
      printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
     
}


int main()
{
    unit_test_generate_I_Q_rectwave_bitstream_scalar();
    unit_test_generate_I_Q_rectwave_bitstream_sse();
    unit_test_generate_I_Q_rectwave_bitstream_sse_u4x();
    unit_test_generate_I_Q_rectwave_bitstream_avx();
    unit_test_generate_I_Q_rectwave_bitstream_avx_u4x();
    unit_test_generate_I_Q_rectwave_bitstream_avx512();
    unit_test_generate_I_Q_rectwave_bitstream_avx512_u4x();
    return 0;
}