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
      constexpr float I_sample_rate{10.0f};
      constexpr float I_w0{4.0f};
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
          const std::int32_t ret_status = IQ_rand_bitstream.generate_channel_I_random_bitstream(0.0f,0.0f,0);
          if(ret_status<0)
          {
             std::printf("[UNIT-TEST]: -- generate_I_channel_random_bitstream -- **FAILED**, status=%d\n",ret_status);
            
          }
          std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
          iq_random_bitstream_t::create_signal_plot(static_cast<std::uint32_t>(IQ_rand_bitstream.m_I_nsamples),IQ_rand_bitstream.m_I_bitstream.m_data,
                                                                                 nullptr,"generate_I_channel_random_bitstream_test_5_",
                                                                                 "I-Channel random bitstream",false);
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      {
         std::printf("[UNIT-TEST]: -- Calling: Q-channel user data random bitstream generator\n");
         const std::int32_t ret_status = IQ_rand_bitstream.generate_channel_Q_random_bitstream(0.0f,0.0f,0);
         if(ret_status<0)
         {
             std::printf("[UNIT-TEST]: -- generate_Q_channel_random_bitstream -- **FAILED**, status=%d\n",ret_status);
             
         }
         std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
         iq_random_bitstream_t::create_signal_plot(static_cast<std::uint32_t>(IQ_rand_bitstream.m_Q_nsamples),IQ_rand_bitstream.m_Q_bitstream.m_data,
                                                                                 nullptr,"generate_Q_channel_random_bitstream_test_5_",
                                                                                 "Q-Channel random bitstream",false);
      }
      printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void
unit_test_iq_random_bitstream_randomized();

void unit_test_iq_random_bitstream_randomized()
{
      using namespace gms::radiolocation;
      using namespace gms;
      std::vector<iq_random_bitstream_t>                       iq_rand_bitstreams_vec;
      thread_local std::uniform_int_distribution<std::int32_t> rand_niters;
      thread_local std::uniform_real_distribution<float>       rand_I_duration;
      thread_local std::uniform_real_distribution<float>       rand_I_sample_rate;
      thread_local std::uniform_real_distribution<float>       rand_I_w0;
      thread_local std::uniform_real_distribution<float>       rand_I_ph0;
      thread_local std::uniform_real_distribution<float>       rand_Q_duration;
      thread_local std::uniform_real_distribution<float>       rand_Q_sample_rate;
      thread_local std::uniform_real_distribution<float>       rand_Q_w0;
      thread_local std::uniform_real_distribution<float>       rand_Q_ph0;
      thread_local unsigned long long seed_niters{0ull};
      thread_local unsigned long long seed_I_duration{0ull};
      thread_local unsigned long long seed_I_sample_rate{0ull};
      thread_local unsigned long long seed_I_w0{0ull};
      thread_local unsigned long long seed_I_ph0{0ull};
      thread_local unsigned long long seed_Q_duration{0ull};
      thread_local unsigned long long seed_Q_sample_rate{0ull};
      thread_local unsigned long long seed_Q_w0{0ull};
      thread_local unsigned long long seed_Q_ph0{0ull};
      std::int32_t n_iters;
      int32_t status{-9999};
      constexpr std::int32_t lo_niters{5};
      constexpr std::int32_t hi_niters{50};
      printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
      // Initalization of random generators.
      seed_niters        = __rdtsc();
      rand_niters        = std::uniform_int_distribution<std::int32_t>(lo_niters,hi_niters);
      auto rgen_niters{std::mt19937(seed_niters)};
      seed_I_duration    = __rdtsc();
      rand_I_duration    = std::uniform_real_distribution<float>(1.0f,50.0f);
      auto rgen_I_duration{std::mt19937(seed_I_duration)};
      seed_I_sample_rate = __rdtsc();
      rand_I_sample_rate = std::uniform_real_distribution<float>(20.0f,100.0f);
      auto rgen_I_sample_rate{std::mt19937(seed_I_sample_rate)};
      seed_I_w0          = __rdtsc();
      rand_I_w0          = std::uniform_real_distribution<float>(2.0f,1000.0f);
      auto rgen_I_w0{std::mt19937(seed_I_w0)};
      seed_I_ph0         = __rdtsc();
      rand_I_ph0         = std::uniform_real_distribution<float>(-1.57079632679489661923132169164f,+1.57079632679489661923132169164f);
      auto rgen_I_ph0{std::mt19937(seed_I_ph0)};
      seed_Q_duration    = __rdtsc();
      rand_Q_duration    = std::uniform_real_distribution<float>(10.0f,70.0f);
      auto rgen_Q_duration{std::mt19937(seed_Q_duration)};
      seed_Q_sample_rate = __rdtsc();
      rand_Q_sample_rate = std::uniform_real_distribution<float>(10.0f,30.0f);
      auto rgen_Q_sample_rate{std::mt19937(seed_Q_sample_rate)};
      seed_Q_w0          = __rdtsc();
      rand_Q_w0          = std::uniform_real_distribution<float>(2.0f,10.0f);
      auto rgen_Q_w0{std::mt19937(seed_Q_w0)};
      seed_Q_ph0         = __rdtsc();
      rand_Q_ph0         = std::uniform_real_distribution<float>(0.0f,+1.57079632679489661923132169164f);
      auto rgen_Q_ph0{std::mt19937(seed_Q_ph0)};

      n_iters = rand_niters.operator()(rgen_niters);
      // sanity checks;
      std::printf("[UNIT-TEST]: Sanity-Check: -- n_iters=%d\n",n_iters);
      if(n_iters < lo_niters) {n_iters = lo_niters;} 
      if(n_iters > hi_niters) {n_iters = hi_niters;}
      //iq_rand_bitstreams_vec.reserve(static_cast<std::size_t>(n_iters));
      iq_rand_bitstreams_vec = std::vector<iq_random_bitstream_t>(static_cast<std::size_t>(n_iters));
      for(std::size_t i{0}; i != static_cast<std::size_t>(n_iters); ++i) 
      {
          const float rval_I_duration{rand_I_duration.operator()(rgen_I_duration)};
          const float rval_I_sample_rate{rand_I_sample_rate.operator()(rgen_I_sample_rate)};
          const float rval_I_w0{rand_I_w0.operator()(rgen_I_w0)};
          const float rval_I_ph0{rand_I_ph0.operator()(rgen_I_ph0)};

          const float rval_Q_duration{rand_Q_duration.operator()(rgen_Q_duration)};
          const float rval_Q_sample_rate{rand_Q_sample_rate.operator()(rgen_Q_sample_rate)};
          const float rval_Q_w0{rand_Q_w0.operator()(rgen_Q_w0)};
          const float rval_Q_ph0{rand_Q_ph0.operator()(rgen_Q_ph0)};

          iq_rand_bitstreams_vec.operator[](i).operator=(iq_random_bitstream_t(rval_I_duration,rval_I_sample_rate,
                                                                               rval_I_w0,rval_I_ph0,
                                                                               rval_Q_duration,rval_Q_sample_rate,
                                                                               rval_Q_w0,rval_Q_ph0));
          
      }
      
      for(std::size_t i{0}; i != iq_rand_bitstreams_vec.size(); ++i) 
      {
           const iq_random_bitstream_t IQ_rand_bitstream = iq_rand_bitstreams_vec.operator[](i);
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

           std::printf("[UNIT-TEST]: -- Calling: #%llu, I-channel user data random bitstream generator.\n",i);
           const std::int32_t I_rand_ch_ret_stat = iq_rand_bitstreams_vec.operator[](i).generate_channel_I_random_bitstream(0.0f,0.0f,0);
           if(I_rand_ch_ret_stat<0)
           {
             std::printf("[UNIT-TEST]: -- generate_I_channel_random_bitstream -- **FAILED**, bitstream=%llu, status=%d\n",i,I_rand_ch_ret_stat);
             return;
           }
           
           std::printf("[UNIT-TEST]: -- Calling: #%llu, Q-channel user data random bitstream generator.\n",i);
           const std::int32_t Q_rand_ch_ret_stat = iq_rand_bitstreams_vec.operator[](i).generate_channel_Q_random_bitstream(0.0f,0.0f,0);
           if(Q_rand_ch_ret_stat<0)
           {
             std::printf("[UNIT-TEST]: -- generate_Q_channel_random_bitstream -- **FAILED**, bitstream=%llu, status=%d\n",i,Q_rand_ch_ret_stat);
             return;
           }
           std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command files.\n");
           std::string I_channel_descr = "generate_I_channel_random_bitstream_test_1_iter_"+std::to_string(i);
           iq_random_bitstream_t::create_signal_plot(static_cast<std::int32_t>(iq_rand_bitstreams_vec.operator[](i).m_I_nsamples),
                                                     iq_rand_bitstreams_vec.operator[](i).m_I_bitstream.m_data,nullptr,
                                                     I_channel_descr.c_str(),"I-Channel random bitstream",false);

           std::string Q_channel_descr = "generate_Q_channel_random_bitstream_test_1_iter_"+std::to_string(i);
           iq_random_bitstream_t::create_signal_plot(static_cast<std::int32_t>(iq_rand_bitstreams_vec.operator[](i).m_Q_nsamples),
                                                     iq_rand_bitstreams_vec.operator[](i).m_Q_bitstream.m_data,nullptr,
                                                     Q_channel_descr.c_str(),"Q-Channel random bitstream",false);
      }
       printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);

}


int main() 
{
    unit_test_iq_random_bitstream();
    unit_test_iq_random_bitstream_randomized();
    return 0;
}