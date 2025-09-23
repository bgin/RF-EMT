#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_dyn_array.h"
#include "GMS_cmplx_trapezw_env.h"

/*
   icpc -o unit_test_cmplx_trapezw_env -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_cmplx_trapezw_env.h GMS_cmplx_trapezw_env.cpp unit_test_cmplx_trapezw_env.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cmplx_trapezw_env.h GMS_cmplx_trapezw_env.cpp unit_test_cmplx_trapezw_env.cpp

*/

__attribute__((hot))
__attribute__((noinline))
void unit_test_cmplx_trapezw_env_s256();

void unit_test_cmplx_trapezw_env_s256()
{
     using namespace gms::radiolocation;
     using namespace gms;
     constexpr std::size_t   I_n_samples{256ull};
     constexpr std::uint32_t I_n_K{8ull};
     constexpr std::size_t   Q_n_samples{256ull};
     constexpr std::uint32_t Q_n_K{8ull};
     constexpr float I_a{0.1f};
     constexpr float I_m{50.0f};
     constexpr float I_l{5.0f};
     constexpr float I_c{2.0f};
     constexpr float Q_a{0.1f};
     constexpr float Q_m{50.0f};
     constexpr float Q_l{5.0f};
     constexpr float Q_c{2.0f};
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     cmplx_trapezw_env_t __cmplx_trapw_env__ = cmplx_trapezw_env_t(I_n_samples,Q_n_samples,
                                                                   I_n_K,Q_n_K,
                                                                   I_a,I_m,I_l,I_c,
                                                                   Q_a,Q_m,Q_l,Q_c);

     char * ctor_name{gms::common::demangle(typeid(__cmplx_trapw_env__).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(__cmplx_trapw_env__).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
     std::printf("[UNIT-TEST:] -- I_n_samples=%llu,I_n_K=%d,Q_n_samples=%llu,Q_n_K=%d\n\n",__cmplx_trapw_env__.__I_n_samples__,
                 __cmplx_trapw_env__.__I_n_K__,__cmplx_trapw_env__.__Q_n_samples__,__cmplx_trapw_env__.__Q_n_K__);
     std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     __cmplx_trapw_env__.__I_chan__.info_size_alignment();
     __cmplx_trapw_env__.__Q_chan__.info_size_alignment();
     std::printf("[UNIT-TEST:] -- Calling I/Q trapezoid wave components.\n");
     __cmplx_trapw_env__.chan_I_squarew_modulated(0);
     __cmplx_trapw_env__.chan_Q_squarew_modulated(0);
     std::cout << __cmplx_trapw_env__ << std::endl;
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     cmplx_trapezw_env_t::create_signal_plot(__cmplx_trapw_env__.__I_n_samples__,__cmplx_trapw_env__.__I_chan__.m_data,nullptr,
                                           "chan_I_modulated_256_test_4_","Complex_Trapezoid_Envelope_I_channel",false);
     cmplx_trapezw_env_t::create_signal_plot(__cmplx_trapw_env__.__Q_n_samples__,__cmplx_trapw_env__.__Q_chan__.m_data,nullptr,
                                           "chan_Q_modulated_256_test_4_","Complex_Trapezoid_Envelope_Q_channel",false);
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_sin_cos_square_wave();

void unit_test_sin_cos_square_wave()
{
     using namespace gms::radiolocation;
     using namespace gms;
     constexpr std::size_t   I_n_samples{256ull};
     constexpr std::uint32_t I_n_K{8ull};
     constexpr std::size_t   Q_n_samples{256ull};
     constexpr std::uint32_t Q_n_K{8ull};
     constexpr float I_a{0.1f};
     constexpr float I_m{50.0f};
     constexpr float I_l{5.0f};
     constexpr float I_c{2.0f};
     constexpr float Q_a{0.1f};
     constexpr float Q_m{50.0f};
     constexpr float Q_l{5.0f};
     constexpr float Q_c{2.0f};
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     cmplx_trapezw_env_t __cmplx_trapw_env__ = cmplx_trapezw_env_t(I_n_samples,Q_n_samples,
                                                                   I_n_K,Q_n_K,
                                                                   I_a,I_m,I_l,I_c,
                                                                   Q_a,Q_m,Q_l,Q_c);
     //std::cout << __cmplx_trapw_env__ << std::endl;
     char * ctor_name{gms::common::demangle(typeid(__cmplx_trapw_env__).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(__cmplx_trapw_env__).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
     std::printf("[UNIT-TEST:] -- I_n_samples=%llu,I_n_K=%d,Q_n_samples=%llu,Q_n_K=%d\n\n",__cmplx_trapw_env__.__I_n_samples__,
                 __cmplx_trapw_env__.__I_n_K__,__cmplx_trapw_env__.__Q_n_samples__,__cmplx_trapw_env__.__Q_n_K__);
     std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     __cmplx_trapw_env__.__I_chan__.info_size_alignment();
     __cmplx_trapw_env__.__Q_chan__.info_size_alignment();
     std::printf("[UNIT-TEST:] -- Calling square-wave signum approximation.\n");
     const float invT{1.0f/__cmplx_trapw_env__.__I_n_samples__};
     for(std::size_t __i{0}; __i != __cmplx_trapw_env__.__I_n_samples__; ++__i)
     {
         const float t{static_cast<float>(__i)};
         __cmplx_trapw_env__.__I_chan__.m_data[__i] = __cmplx_trapw_env__.sin_squarew_I_sample(t,invT);
         __cmplx_trapw_env__.__Q_chan__.m_data[__i] = __cmplx_trapw_env__.cos_squarew_Q_sample(t,invT);
     }
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     cmplx_trapezw_env_t::create_signal_plot(__cmplx_trapw_env__.__I_n_samples__,__cmplx_trapw_env__.__I_chan__.m_data,nullptr,
                                             "sin_squarew_256_test_2_","Sin_squarew_test_2",false);
     cmplx_trapezw_env_t::create_signal_plot(__cmplx_trapw_env__.__Q_n_samples__,__cmplx_trapw_env__.__Q_chan__.m_data,nullptr,
                                             "cos_squarew_256_test_2_","Cos_squarew_test_2",false);
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);


}





int main()
{   
    unit_test_cmplx_trapezw_env_s256();
    unit_test_sin_cos_square_wave();
    return 0;
}