#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_dyn_array.h"
#include "GMS_am_bb_cosine_sequence.h"

/*
   icpc -o unit_test_am_bb_cosine_sequence_avx -O3 -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h GMS_half.h GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_indices.h GMS_am_bb_cosine_signal.h GMS_am_bb_cosine_signal.cpp GMS_am_bb_cosine_sequence.h GMS_am_bb_cosine_sequence.cpp unit_test_am_bb_cosine_sequence_avx.cpp
   ASM: 
   icpc -S --O3 fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_half.h GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_indices.h GMS_am_bb_cosine_sequence.h GMS_am_bb_cosine_sequence.cpp unit_test_am_bb_cosine_sequence_avx.cpp

*/

__attribute__((hot))
__attribute__((noinline))
void unit_test_am_bb_cosine_sequence_rand_data_avx_u16x();

void unit_test_am_bb_cosine_sequence_rand_data_avx_u16x()
{
     using namespace gms::radiolocation;
     using namespace gms;
     constexpr std::size_t   I_n_samples{1024ull};
     constexpr std::uint32_t I_n_K{8};
     constexpr std::size_t   nsignals{4ull};
     constexpr float A{1.0f};
     constexpr float n{16.0f};
     constexpr float P{8.0f};
     darray_r4_t I_sym; 
     std::vector<float> A_values          = std::vector<float>(nsignals);
     std::vector<float> n_values          = std::vector<float>(nsignals);
     std::vector<float> P_values          = std::vector<float>(nsignals);
     std::vector<std::uint32_t> nK_values = std::vector<std::uint32_t>(nsignals);
     std::vector<std::int32_t>  ret_codes = std::vector<std::int32_t>(nsignals);
     for(std::size_t i{0}; i != nsignals; ++i) 
     {
         A_values.operator[](i)  = A;
         n_values.operator[](i)  = n;
         P_values.operator[](i)  = P;
         nK_values.operator[](i) = I_n_K; 
     }
     std::clock_t seed_I;
     constexpr std::size_t I_sym_len{I_n_samples*static_cast<std::size_t>(I_n_K)};
    
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     am_bb_cos_signal_sequence_t cos_sequence = am_bb_cos_signal_sequence_t(I_n_samples,
                                                                             nsignals,nK_values,
                                                                             n_values,A_values,P_values);
     I_sym = darray_r4_t(I_sym_len);
     char * ctor_name{gms::common::demangle(typeid(cos_sequence).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(cos_sequence).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
     //std::printf("[UNIT-TEST:] -- I_n_samples=%llu,I_n_K=%d,Q_n_samples=%llu,Q_n_K=%d\n\n",__cmplx_trapw_env__.__I_n_samples__,
     //            __cmplx_trapw_env__.__I_n_K__,__cmplx_trapw_env__.__Q_n_samples__,__cmplx_trapw_env__.__Q_n_K__);

     //std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     //symbol_1.m_sig_samples.info_size_alignment();
     std::printf("[UNIT-TEST:] -- Generating random symbol data: I_channel(-1,1)\n");
     seed_I = std::clock();
     auto rand_sym_I{std::bind(std::uniform_real_distribution<float>(-1.0f,1.0f),std::mt19937(seed_I))};
     
     for(std::size_t __i{0ull}; __i != I_sym_len; ++__i) 
     {
          I_sym.m_data[__i] = rand_sym_I();
          //std::printf("%.7f\n",I_sym.m_data[__i]);
     }
     std::printf("[UNIT-TEST:] -- Calling  cosine wave components modulated by the user data -- vertically summed sequence.\n");
     cos_sequence.signal_sequence_avx_u16x_u(&I_sym.m_data[0],&ret_codes.operator[](0));     
     //std::cout << __cmplx_trapw_env__ << std::endl;
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     am_bb_cos_signal_sequence_t::create_sequence_plot(cos_sequence.m_nsamples,&cos_sequence.m_vsummed_sequence.operator[](0),nullptr,
                                              "am_bb_cosine_sequence_avx_u16x_test_1_","AM-BB_Cosine_Signals_Sequence",false);
     
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_am_bb_cosine_sequence_varying_params_avx_u16x();

void unit_test_am_bb_cosine_sequence_varying_params_avx_u16x()
{
     using namespace gms::radiolocation;
     using namespace gms;
     constexpr std::size_t   I_n_samples{1024ull};
     constexpr std::uint32_t I_n_K{64};
     constexpr std::size_t   nsignals{128ull};
     darray_r4_t I_sym; 
     std::vector<float> A_values          = std::vector<float>(nsignals);
     std::vector<float> n_values          = std::vector<float>(nsignals);
     std::vector<float> P_values          = std::vector<float>(nsignals);
     std::vector<std::uint32_t> nK_values = std::vector<std::uint32_t>(nsignals);
     std::vector<std::int32_t>  ret_codes = std::vector<std::int32_t>(nsignals);
     std::clock_t seed_A;
     std::clock_t seed_n;
     std::clock_t seed_P;
     seed_A = std::clock();
     auto rand_A{std::bind(std::uniform_real_distribution<float>(0.01f,1.0f),std::mt19937(seed_A))};
     seed_n = std::clock();
     auto rand_n{std::bind(std::uniform_real_distribution<float>(1.0f,128.0f),std::mt19937(seed_n))};
     seed_P = std::clock();
     auto rand_P{std::bind(std::uniform_real_distribution<float>(1.0f,128.0f),std::mt19937(seed_P))};
     for(std::size_t i{0}; i != nsignals; ++i) 
     {
         A_values.operator[](i)  = rand_A.operator()();
         n_values.operator[](i)  = rand_n.operator()();
         P_values.operator[](i)  = rand_P.operator()();
         nK_values.operator[](i) = I_n_K; 
     }
     std::clock_t seed_I;
     constexpr std::size_t I_sym_len{I_n_samples*static_cast<std::size_t>(I_n_K)};
    
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     am_bb_cos_signal_sequence_t cos_sequence = am_bb_cos_signal_sequence_t(I_n_samples,
                                                                             nsignals,nK_values,
                                                                             n_values,A_values,P_values);
     I_sym = darray_r4_t(I_sym_len);
     char * ctor_name{gms::common::demangle(typeid(cos_sequence).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(cos_sequence).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
     //std::printf("[UNIT-TEST:] -- I_n_samples=%llu,I_n_K=%d,Q_n_samples=%llu,Q_n_K=%d\n\n",__cmplx_trapw_env__.__I_n_samples__,
     //            __cmplx_trapw_env__.__I_n_K__,__cmplx_trapw_env__.__Q_n_samples__,__cmplx_trapw_env__.__Q_n_K__);

     //std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     //symbol_1.m_sig_samples.info_size_alignment();
     std::printf("[UNIT-TEST:] -- Generating random symbol data: I_channel(-1,1)\n");
     seed_I = std::clock();
     auto rand_sym_I{std::bind(std::uniform_real_distribution<float>(-1.0f,1.0f),std::mt19937(seed_I))};
     
     for(std::size_t __i{0ull}; __i != I_sym_len; ++__i) 
     {
          I_sym.m_data[__i] = rand_sym_I();
          //std::printf("%.7f\n",I_sym.m_data[__i]);
     }
     std::printf("[UNIT-TEST:] -- Calling  cosine wave components modulated by the user data -- vertically summed sequence.\n");
     cos_sequence.signal_sequence_avx_u16x_u(&I_sym.m_data[0],&ret_codes.operator[](0));     
     for(std::size_t i{0}; i != cos_sequence.m_nsignals; ++i) 
     {
          const std::int32_t ret_code{ret_codes.operator[](i)};
          if(ret_code < 0ull) 
          {
               std::printf("[**ERROR**]: -- Found a failed call at index=%llu, of value=%d\n",i,ret_code);
          }
     }
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     am_bb_cos_signal_sequence_t::create_sequence_plot(cos_sequence.m_nsamples,&cos_sequence.m_vsummed_sequence.operator[](0),nullptr,
                                              "am_bb_cosine_sequence_avx_u16x_test_2_","AM-BB_Cosine_Signals_Sequence",false);
     
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}


__attribute__((hot))
__attribute__((noinline))
void unit_test_am_bb_cosine_sequence_rand_data_avx_u10x();

void unit_test_am_bb_cosine_sequence_rand_data_avx_u10x()
{
     using namespace gms::radiolocation;
     using namespace gms;
     constexpr std::size_t   I_n_samples{1024ull};
     constexpr std::uint32_t I_n_K{8};
     constexpr std::size_t   nsignals{128ull};
     constexpr float A{1.0f};
     constexpr float n{16.0f};
     constexpr float P{8.0f};
     darray_r4_t I_sym; 
     std::vector<float> A_values          = std::vector<float>(nsignals);
     std::vector<float> n_values          = std::vector<float>(nsignals);
     std::vector<float> P_values          = std::vector<float>(nsignals);
     std::vector<std::uint32_t> nK_values = std::vector<std::uint32_t>(nsignals);
     std::vector<std::int32_t>  ret_codes = std::vector<std::int32_t>(nsignals);
     for(std::size_t i{0}; i != nsignals; ++i) 
     {
         A_values.operator[](i)  = A;
         n_values.operator[](i)  = n;
         P_values.operator[](i)  = P;
         nK_values.operator[](i) = I_n_K; 
         
     }
     std::clock_t seed_I;
     constexpr std::size_t I_sym_len{I_n_samples*static_cast<std::size_t>(I_n_K)};
    
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     am_bb_cos_signal_sequence_t cos_sequence = am_bb_cos_signal_sequence_t(I_n_samples,
                                                                             nsignals,nK_values,
                                                                             n_values,A_values,P_values);
     I_sym = darray_r4_t(I_sym_len);
     char * ctor_name{gms::common::demangle(typeid(cos_sequence).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(cos_sequence).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
     //std::printf("[UNIT-TEST:] -- I_n_samples=%llu,I_n_K=%d,Q_n_samples=%llu,Q_n_K=%d\n\n",__cmplx_trapw_env__.__I_n_samples__,
     //            __cmplx_trapw_env__.__I_n_K__,__cmplx_trapw_env__.__Q_n_samples__,__cmplx_trapw_env__.__Q_n_K__);

     //std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     //symbol_1.m_sig_samples.info_size_alignment();
     std::printf("[UNIT-TEST:] -- Generating random symbol data: I_channel(-1,1)\n");
     seed_I = std::clock();
     auto rand_sym_I{std::bind(std::uniform_real_distribution<float>(-1.0f,1.0f),std::mt19937(seed_I))};
     
     for(std::size_t __i{0ull}; __i != I_sym_len; ++__i) 
     {
          I_sym.m_data[__i] = rand_sym_I();
          //std::printf("%.7f\n",I_sym.m_data[__i]);
     }
     std::printf("[UNIT-TEST:] -- Calling  cosine wave components modulated by the user data -- vertically summed sequence.\n");
     cos_sequence.signal_sequence_avx_u10x_u(&I_sym.m_data[0],&ret_codes.operator[](0));     
     //std::cout << __cmplx_trapw_env__ << std::endl;
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     am_bb_cos_signal_sequence_t::create_sequence_plot(cos_sequence.m_nsamples,&cos_sequence.m_vsummed_sequence.operator[](0),nullptr,
                                              "am_bb_cosine_sequence_avx_u10x_test_1_","AM-BB_Cosine_Signals_Sequence",false);
     
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_am_bb_cosine_sequence_rand_data_avx_u6x();

void unit_test_am_bb_cosine_sequence_rand_data_avx_u6x()
{
     using namespace gms::radiolocation;
     using namespace gms;
     constexpr std::size_t   I_n_samples{1024ull};
     constexpr std::uint32_t I_n_K{8};
     constexpr std::size_t   nsignals{128ull};
     constexpr float A{1.0f};
     constexpr float n{16.0f};
     constexpr float P{8.0f};
     darray_r4_t I_sym; 
     std::vector<float> A_values          = std::vector<float>(nsignals);
     std::vector<float> n_values          = std::vector<float>(nsignals);
     std::vector<float> P_values          = std::vector<float>(nsignals);
     std::vector<std::uint32_t> nK_values = std::vector<std::uint32_t>(nsignals);
     std::vector<std::int32_t>  ret_codes = std::vector<std::int32_t>(nsignals);
     for(std::size_t i{0}; i != nsignals; ++i) 
     {
         A_values.operator[](i)  = A;
         n_values.operator[](i)  = n;
         P_values.operator[](i)  = P;
         nK_values.operator[](i) = I_n_K; 
     }
     std::clock_t seed_I;
     constexpr std::size_t I_sym_len{I_n_samples*static_cast<std::size_t>(I_n_K)};
    
     int32_t status{};
     //__asm__ ("int3");
     //printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     am_bb_cos_signal_sequence_t cos_sequence = am_bb_cos_signal_sequence_t(I_n_samples,
                                                                             nsignals,nK_values,
                                                                             n_values,A_values,P_values);
     I_sym = darray_r4_t(I_sym_len);
     char * ctor_name{gms::common::demangle(typeid(cos_sequence).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(cos_sequence).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
     //std::printf("[UNIT-TEST:] -- I_n_samples=%llu,I_n_K=%d,Q_n_samples=%llu,Q_n_K=%d\n\n",__cmplx_trapw_env__.__I_n_samples__,
     //            __cmplx_trapw_env__.__I_n_K__,__cmplx_trapw_env__.__Q_n_samples__,__cmplx_trapw_env__.__Q_n_K__);

     //std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     //symbol_1.m_sig_samples.info_size_alignment();
     std::printf("[UNIT-TEST:] -- Generating random symbol data: I_channel(-1,1)\n");
     seed_I = std::clock();
     auto rand_sym_I{std::bind(std::uniform_real_distribution<float>(-1.0f,1.0f),std::mt19937(seed_I))};
     
     for(std::size_t __i{0ull}; __i != I_sym_len; ++__i) 
     {
          I_sym.m_data[__i] = rand_sym_I();
          //std::printf("%.7f\n",I_sym.m_data[__i]);
     }
     std::printf("[UNIT-TEST:] -- Calling  cosine wave components modulated by the user data -- vertically summed sequence.\n");
     cos_sequence.signal_sequence_avx_u6x_u(&I_sym.m_data[0],&ret_codes.operator[](0));     
     //std::cout << __cmplx_trapw_env__ << std::endl;
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     am_bb_cos_signal_sequence_t::create_sequence_plot(cos_sequence.m_nsamples,&cos_sequence.m_vsummed_sequence.operator[](0),nullptr,
                                              "am_bb_cosine_sequence_avx_u6x_test_1_","AM-BB_Cosine_Signals_Sequence",false);
     
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}


int main()
{
    unit_test_am_bb_cosine_sequence_rand_data_avx_u16x();
    unit_test_am_bb_cosine_sequence_rand_data_avx_u10x();
    unit_test_am_bb_cosine_sequence_rand_data_avx_u6x();
    unit_test_am_bb_cosine_sequence_varying_params_avx_u16x();
    return 0;
}
