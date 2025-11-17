#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <cmath>
#include "GMS_dyn_array.h"
#include "GMS_am_bb_sine_sequence.h"

/*
   icpc -o unit_test_am_bb_sine_sequence_avx512 -O3 -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h GMS_half.h GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_indices.h GMS_am_bb_cosine_signal.h GMS_am_bb_cosine_signal.cpp GMS_am_bb_sine_sequence.h GMS_am_bb_sine_sequence.cpp unit_test_am_bb_sine_sequence_avx512.cpp
   ASM: 
   icpc -S --O3 fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_half.h GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_indices.h GMS_am_bb_sine_sequence.h GMS_am_bb_sine_sequence.cpp unit_test_am_bb_sine_sequence_avx512.cpp

*/

__attribute__((hot))
__attribute__((noinline))
void unit_test_am_bb_sine_sequence_rand_data_avx512_u16x();

void unit_test_am_bb_sine_sequence_rand_data_avx512_u16x()
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
     am_bb_sine_signal_sequence_t sin_sequence = am_bb_sine_signal_sequence_t(I_n_samples,
                                                                             nsignals,nK_values,
                                                                             n_values,A_values,P_values);
     I_sym = darray_r4_t(I_sym_len);
     char * ctor_name{gms::common::demangle(typeid(sin_sequence).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(sin_sequence).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
    
     std::printf("[UNIT-TEST:] -- Generating random symbol data: I_channel(-1,1)\n");
     seed_I = std::clock();
     auto rand_sym_I{std::bind(std::uniform_real_distribution<float>(-1.0f,1.0f),std::mt19937(seed_I))};
     
     for(std::size_t __i{0ull}; __i != I_sym_len; ++__i) 
     {
          I_sym.m_data[__i] = rand_sym_I();
          //std::printf("%.7f\n",I_sym.m_data[__i]);
     }
     std::printf("[UNIT-TEST:] -- Calling  sine wave components modulated by the user data -- vertically summed sequence.\n");
     std::int32_t ret_stat = sin_sequence.signal_sequence_avx512_u16x(&I_sym.m_data[0],&ret_codes.operator[](0),1); 
     if(ret_stat<0) { 
        std::printf("[UNIT-TEST:] -- **FAIL** -> signal_sequence_avx512_u16x stat=%d\n",ret_stat);
        return;
     }    
     //std::cout << __cmplx_trapw_env__ << std::endl;
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     am_bb_sine_signal_sequence_t::create_sequence_plot(sin_sequence.m_nsamples,&sin_sequence.m_vsummed_sequence.operator[](0),nullptr,
                                              "am_bb_sine_sequence_avx512_u16x_rand_data_test_1_","AM-BB_sine_Signals_Sequence",false);
     
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void signal_sequence_vertical_summation(const float * __restrict__,
                                        const float * __restrict__,
                                        const float * __restrict__,
                                        const float * __restrict__,
                                        float * __restrict__,
                                        const std::size_t);

void signal_sequence_vertical_summation(const float * __restrict__ sig_1,
                                        const float * __restrict__ sig_2,
                                        const float * __restrict__ sig_3,
                                        const float * __restrict__ sig_4,
                                        float * __restrict__       sig_vsumm,
                                        const std::size_t          sig_len) 
{
     std::size_t i,j;
     for(i = 0ull; i != ROUND_TO_SIXTEEN(sig_len,16ull); i += 16ull) 
     {
         __m512 xmm1 = _mm512_load_ps(&sig_1[i]);
         __m512 xmm2 = _mm512_load_ps(&sig_2[i]);
         __m512 sum0 = _mm512_add_ps(xmm1,xmm2);
         __m512 xmm3 = _mm512_load_ps(&sig_3[i]);
         __m512 xmm4 = _mm512_load_ps(&sig_4[i]);
         __m512 sum1 = _mm512_add_ps(xmm3,xmm4);
         __m512 sum2 = _mm512_add_ps(sum0,sum1);
         _mm512_store_ps(&sig_vsumm[i],sum2);
     }

     for(j = i; j != sig_len; j += 1ull) 
     {
          const float x1{sig_1[j]};
          const float x2{sig_2[j]};
          const float sum0{x1+x2};
          const float x3{sig_3[j]};
          const float x4{sig_4[j]};
          const float sum1{x3+x4};
          const float sum2{sum0+sum1};
          sig_vsumm[j] = sum2;
     }
}    

static inline 
bool approximatelyEqual(const float a,
		              const float b,
				    const float epsilon) 
{
	const float fabsa = std::fabs(a);
	const float fabsb = std::fabs(b);
     return std::fabs(a - b) <=
            ((fabsa < fabsb ? fabsb : fabsa) * epsilon);
}

__attribute__((hot))
__attribute__((noinline))
void check_fp_approx_equality(const float * __restrict__,
                              const float * __restrict__,
                              const std::size_t);

void check_fp_approx_equality(const float * __restrict__ sig_orig,
                              const float * __restrict__ sig_ref,
                              const std::size_t sig_len) 
{
     for(std::size_t i{0ull}; i != sig_len; ++i) 
     {
          const float so{sig_orig[i]};
          const float sr{sig_ref[i]};
          const bool approx_equal{approximatelyEqual(so,sr,std::numeric_limits<double>::epsilon())};
          if(!approx_equal)
          {    
               const float delta{so-sr};
               std::printf("[**TEST-FAIL**] -- Found inequality at loc=%llu, orig=%.7f,ref=%.7f,delta=%.7f\n",
                           i,so,sr,delta );
          }
     }
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_am_bb_sine_sequence_correctness_avx512_u16x();

void unit_test_am_bb_sine_sequence_correctness_avx512_u16x()
{
     using namespace gms::radiolocation;
     using namespace gms;
     constexpr std::size_t   I_n_samples{4096ull};
     constexpr std::uint32_t I_n_K{8};
     constexpr std::size_t   nsignals{4ull};
     constexpr float A{1.0f};
     constexpr float n{16.0f};
     constexpr float P{8.0f};
    // std::clock_t seed_I_n_samples;
    // std::clock_t seed_I_n_K;
    // std::clock_t seed_nsignals;
    // seed_I_n_samples = std::clock();
    // auto rand_I_n_samples{std::bind(std::uniform_int_distribution<std::size_t>(124ull,8192ull),std::mt19937(seed_I_n_samples))};
    // seed_I_n_K       = std::clock();
   //  auto rand_I_n_K{std::bind(std::uniform_int_distribution<std::uint32_t>(4,128),std::mt19937(seed_I_n_K))};
   //  seed_nsignals    = std::clock();
    // auto rand_nsignals{std::bind(std::uniform_int_distribution<std::size_t>(2ull,256ull),std::mt19937(seed_nsignals))};
   //  const std::size_t   I_n_samples{rand_I_n_samples.operator()()};
   //  const std::uint32_t I_n_K{rand_I_n_K.operator()()};
   //  const std::size_t   nsignals{rand_nsignals.operator()()};
     darray_r4_t I_sym; 
     darray_r4_t sig_vsum_reference;
     std::vector<float> A_values          = std::vector<float>(nsignals);
     std::vector<float> n_values          = std::vector<float>(nsignals);
     std::vector<float> P_values          = std::vector<float>(nsignals);
     std::vector<std::uint32_t> nK_values = std::vector<std::uint32_t>(nsignals);
     std::vector<std::int32_t>  ret_codes = std::vector<std::int32_t>(nsignals);
     std::clock_t seed_func{};
     for(std::size_t i{0}; i != nsignals; ++i) 
     {
         A_values.operator[](i)  = A;
         n_values.operator[](i)  = n;
         P_values.operator[](i)  = P;
         nK_values.operator[](i) = I_n_K; 
     }
     std::clock_t seed_I;
     const std::size_t I_sym_len{I_n_samples*static_cast<std::size_t>(I_n_K)};
    
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n\n", __PRETTY_FUNCTION__);
     am_bb_sine_signal_sequence_t sin_sequence = am_bb_sine_signal_sequence_t(I_n_samples,
                                                                             nsignals,nK_values,
                                                                             n_values,A_values,P_values);
     am_bb_sine_signal_t symbol_1 = am_bb_sine_signal_t(I_n_samples,I_n_K,n,A,P);
     am_bb_sine_signal_t symbol_2 = am_bb_sine_signal_t(I_n_samples,I_n_K,n,A,P);
     am_bb_sine_signal_t symbol_3 = am_bb_sine_signal_t(I_n_samples,I_n_K,n,A,P);
     am_bb_sine_signal_t symbol_4 = am_bb_sine_signal_t(I_n_samples,I_n_K,n,A,P);
          
     I_sym = darray_r4_t(I_sym_len);
     sig_vsum_reference = darray_r4_t(I_n_samples);
     std::printf("[UNIT-TEST] -- Number of samples=%llu,number of narrowband signals=%d,number of baseband signals=%llu\n\n",
                            I_n_samples,I_n_K,nsignals);
     char * ctor_name{gms::common::demangle(typeid(sin_sequence).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(sin_sequence).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
    
     std::printf("[UNIT-TEST:] -- Generating random symbol data: I_channel(-1,1)\n\n");
     seed_I = std::clock();
     auto rand_sym_I{std::bind(std::uniform_real_distribution<float>(-1.0f,1.0f),std::mt19937(seed_I))};
     
     for(std::size_t __i{0ull}; __i != I_sym_len; ++__i) 
     {
          I_sym.m_data[__i] = rand_sym_I();
          //std::printf("%.7f\n",I_sym.m_data[__i]);
     }
     seed_func = std::clock();
     auto rand_func{std::bind(std::uniform_int_distribution<std::uint32_t>(0,1),std::mt19937(seed_func))};
     const std::uint32_t which_func{rand_func()};
     //cos_sequence.signal_sequence_avx_u16x_u(&I_sym.m_data[0],&ret_codes.operator[](0)); 
     std::int32_t ret_stat = sin_sequence.signal_sequence_avx512_u16x(&I_sym.m_data[0],&ret_codes.operator[](0),which_func);
     if(ret_stat<0) { 
        std::printf("[UNIT-TEST:] -- **FAIL** -> signal_sequence_avx512_u16x stat=%d\n",ret_stat);
        return;
     }
     std::printf("[UNIT-TEST:] -- Calling  sine wave components modulated by the user data -- vertically summed sequence (REFERENCE!!).\n\n");
     switch (which_func) 
     {
     case 0 : 
         symbol_1.create_signal_user_data(&I_sym.m_data[0],I_n_samples,I_n_K);
         symbol_2.create_signal_user_data(&I_sym.m_data[0],I_n_samples,I_n_K); 
         symbol_3.create_signal_user_data(&I_sym.m_data[0],I_n_samples,I_n_K);
         symbol_4.create_signal_user_data(&I_sym.m_data[0],I_n_samples,I_n_K);
     break;
     case 1 : 
         symbol_1.create_signal_user_data_u4x(&I_sym.m_data[0],I_n_samples,I_n_K);
         symbol_2.create_signal_user_data_u4x(&I_sym.m_data[0],I_n_samples,I_n_K); 
         symbol_3.create_signal_user_data_u4x(&I_sym.m_data[0],I_n_samples,I_n_K);
         symbol_4.create_signal_user_data_u4x(&I_sym.m_data[0],I_n_samples,I_n_K);
     break;
     default : 
         std::printf("wrong case argument=%d\n",which_func);
         return;
     }
     
     signal_sequence_vertical_summation(&symbol_1.m_sig_samples.m_data[0],
                                        &symbol_2.m_sig_samples.m_data[0],
                                        &symbol_3.m_sig_samples.m_data[0],
                                        &symbol_4.m_sig_samples.m_data[0],
                                        &sig_vsum_reference.m_data[0],
                                        I_n_samples);
     check_fp_approx_equality(&sin_sequence.m_vsummed_sequence.operator[](0),
                             &sig_vsum_reference.m_data[0],I_n_samples);
     //std::cout << __cmplx_trapw_env__ << std::endl;
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n\n");
     am_bb_sine_signal_sequence_t::create_sequence_plot(sin_sequence.m_nsamples,&sin_sequence.m_vsummed_sequence.operator[](0),nullptr,
                                              "am_bb_sine_sequence_avx512_u16x_correctness_test_1_","AM-BB_sine_Signals_Sequence",false);
     am_bb_sine_signal_sequence_t::create_sequence_plot(I_n_samples,&sig_vsum_reference.m_data[0],nullptr,
                                               "am_bb_sine_signal_reference_vsumm_u16x_test_1_","AM-BB_sine_Signals_Sequence",false);
     
     printf("[UNIT_TEST]: function=%s -- **END**\n\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_am_bb_sine_sequence_correctness_avx512_u10x();

void unit_test_am_bb_sine_sequence_correctness_avx512_u10x()
{
     using namespace gms::radiolocation;
     using namespace gms;
     constexpr std::size_t   I_n_samples{4096ull};
     constexpr std::uint32_t I_n_K{8};
     constexpr std::size_t   nsignals{4ull};
     constexpr float A{1.0f};
     constexpr float n{16.0f};
     constexpr float P{8.0f};
     /*
     std::clock_t seed_I_n_samples;
     std::clock_t seed_I_n_K;
     std::clock_t seed_nsignals;
     seed_I_n_samples = std::clock();
     auto rand_I_n_samples{std::bind(std::uniform_int_distribution<std::size_t>(124ull,8192ull),std::mt19937(seed_I_n_samples))};
     seed_I_n_K       = std::clock();
     auto rand_I_n_K{std::bind(std::uniform_int_distribution<std::uint32_t>(4,128),std::mt19937(seed_I_n_K))};
     seed_nsignals    = std::clock();
     auto rand_nsignals{std::bind(std::uniform_int_distribution<std::size_t>(2ull,256ull),std::mt19937(seed_nsignals))};
     const std::size_t   I_n_samples{rand_I_n_samples.operator()()};
     const std::uint32_t I_n_K{rand_I_n_K.operator()()};
     const std::size_t   nsignals{rand_nsignals.operator()()};
     */
     darray_r4_t I_sym; 
     darray_r4_t sig_vsum_reference;
     std::vector<float> A_values          = std::vector<float>(nsignals);
     std::vector<float> n_values          = std::vector<float>(nsignals);
     std::vector<float> P_values          = std::vector<float>(nsignals);
     std::vector<std::uint32_t> nK_values = std::vector<std::uint32_t>(nsignals);
     std::vector<std::int32_t>  ret_codes = std::vector<std::int32_t>(nsignals);
     std::clock_t seed_func{};
     for(std::size_t i{0}; i != nsignals; ++i) 
     {
         A_values.operator[](i)  = A;
         n_values.operator[](i)  = n;
         P_values.operator[](i)  = P;
         nK_values.operator[](i) = I_n_K; 
     }
     std::clock_t seed_I;
     const std::size_t I_sym_len{I_n_samples*static_cast<std::size_t>(I_n_K)};
    
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n\n", __PRETTY_FUNCTION__);
     am_bb_sine_signal_sequence_t sin_sequence = am_bb_sine_signal_sequence_t(I_n_samples,
                                                                             nsignals,nK_values,
                                                                             n_values,A_values,P_values);
     am_bb_sine_signal_t symbol_1 = am_bb_sine_signal_t(I_n_samples,I_n_K,n,A,P);
     am_bb_sine_signal_t symbol_2 = am_bb_sine_signal_t(I_n_samples,I_n_K,n,A,P);
     am_bb_sine_signal_t symbol_3 = am_bb_sine_signal_t(I_n_samples,I_n_K,n,A,P);
     am_bb_sine_signal_t symbol_4 = am_bb_sine_signal_t(I_n_samples,I_n_K,n,A,P);
          
     I_sym = darray_r4_t(I_sym_len);
     sig_vsum_reference = darray_r4_t(I_n_samples);
     std::printf("[UNIT-TEST] -- Number of samples=%llu,number of narrowband signals=%d,number of baseband signals=%llu\n\n",
                            I_n_samples,I_n_K,nsignals);
     char * ctor_name{gms::common::demangle(typeid(sin_sequence).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(sin_sequence).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
    
     std::printf("[UNIT-TEST:] -- Generating random symbol data: I_channel(-1,1)\n\n");
     seed_I = std::clock();
     auto rand_sym_I{std::bind(std::uniform_real_distribution<float>(-1.0f,1.0f),std::mt19937(seed_I))};
     
     for(std::size_t __i{0ull}; __i != I_sym_len; ++__i) 
     {
          I_sym.m_data[__i] = rand_sym_I();
          //std::printf("%.7f\n",I_sym.m_data[__i]);
     }
     //cos_sequence.signal_sequence_avx_u16x_u(&I_sym.m_data[0],&ret_codes.operator[](0)); 
     seed_func = std::clock();
     auto rand_func{std::bind(std::uniform_int_distribution<std::uint32_t>(0,1),std::mt19937(seed_func))};
     const std::uint32_t which_func{rand_func()};
     std::int32_t ret_stat = sin_sequence.signal_sequence_avx512_u10x(&I_sym.m_data[0],&ret_codes.operator[](0),which_func);
     if(ret_stat<0) { 
        std::printf("[UNIT-TEST:] -- **FAIL** -> signal_sequence_avx512_u10x stat=%d\n",ret_stat);
        return;
     } 
     std::printf("[UNIT-TEST:] -- Calling  sine wave components modulated by the user data -- vertically summed sequence (REFERENCE!!).\n\n");
     switch (which_func)
     {
       case 0 :
           symbol_1.create_signal_user_data(&I_sym.m_data[0],I_n_samples,I_n_K);
           symbol_2.create_signal_user_data(&I_sym.m_data[0],I_n_samples,I_n_K); 
           symbol_3.create_signal_user_data(&I_sym.m_data[0],I_n_samples,I_n_K);
           symbol_4.create_signal_user_data(&I_sym.m_data[0],I_n_samples,I_n_K);
       break;
       case 1 :
           symbol_1.create_signal_user_data_u4x(&I_sym.m_data[0],I_n_samples,I_n_K);
           symbol_2.create_signal_user_data_u4x(&I_sym.m_data[0],I_n_samples,I_n_K); 
           symbol_3.create_signal_user_data_u4x(&I_sym.m_data[0],I_n_samples,I_n_K);
           symbol_4.create_signal_user_data_u4x(&I_sym.m_data[0],I_n_samples,I_n_K);
       break;
       default : 
           std::printf("wrong case argument=%d\n",which_func);
           return;
     }
     signal_sequence_vertical_summation(&symbol_1.m_sig_samples.m_data[0],
                                        &symbol_2.m_sig_samples.m_data[0],
                                        &symbol_3.m_sig_samples.m_data[0],
                                        &symbol_4.m_sig_samples.m_data[0],
                                        &sig_vsum_reference.m_data[0],
                                        I_n_samples);
     check_fp_approx_equality(&sin_sequence.m_vsummed_sequence.operator[](0),
                             &sig_vsum_reference.m_data[0],I_n_samples);
     //std::cout << __cmplx_trapw_env__ << std::endl;
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n\n");
     am_bb_sine_signal_sequence_t::create_sequence_plot(sin_sequence.m_nsamples,&sin_sequence.m_vsummed_sequence.operator[](0),nullptr,
                                              "am_bb_cosine_sequence_avx512_u10x_correctness_test_1_","AM-BB_sine_Signals_Sequence",false);
     am_bb_sine_signal_sequence_t::create_sequence_plot(I_n_samples,&sig_vsum_reference.m_data[0],nullptr,
                                               "am_bb_sine_signal_reference_vsumm_u10x_test_1_","AM-BB_sine_Signals_Sequence",false);
     
     printf("[UNIT_TEST]: function=%s -- **END**\n\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_am_bb_sine_sequence_correctness_avx512_u6x();

void unit_test_am_bb_sine_sequence_correctness_avx512_u6x()
{
     using namespace gms::radiolocation;
     using namespace gms;
     constexpr std::size_t   I_n_samples{4096ull};
     constexpr std::uint32_t I_n_K{8};
     constexpr std::size_t   nsignals{4ull};
     constexpr float A{1.0f};
     constexpr float n{16.0f};
     constexpr float P{8.0f};
     /*
     std::clock_t seed_I_n_samples;
     std::clock_t seed_I_n_K;
     std::clock_t seed_nsignals;
     seed_I_n_samples = std::clock();
     auto rand_I_n_samples{std::bind(std::uniform_int_distribution<std::size_t>(124ull,8192ull),std::mt19937(seed_I_n_samples))};
     seed_I_n_K       = std::clock();
     auto rand_I_n_K{std::bind(std::uniform_int_distribution<std::uint32_t>(4,128),std::mt19937(seed_I_n_K))};
     seed_nsignals    = std::clock();
     auto rand_nsignals{std::bind(std::uniform_int_distribution<std::size_t>(2ull,256ull),std::mt19937(seed_nsignals))};
     const std::size_t   I_n_samples{rand_I_n_samples.operator()()};
     const std::uint32_t I_n_K{rand_I_n_K.operator()()};
     const std::size_t   nsignals{rand_nsignals.operator()()};
     */
     darray_r4_t I_sym; 
     darray_r4_t sig_vsum_reference;
     std::vector<float> A_values          = std::vector<float>(nsignals);
     std::vector<float> n_values          = std::vector<float>(nsignals);
     std::vector<float> P_values          = std::vector<float>(nsignals);
     std::vector<std::uint32_t> nK_values = std::vector<std::uint32_t>(nsignals);
     std::vector<std::int32_t>  ret_codes = std::vector<std::int32_t>(nsignals);
     std::clock_t seed_func{};
     for(std::size_t i{0}; i != nsignals; ++i) 
     {
         A_values.operator[](i)  = A;
         n_values.operator[](i)  = n;
         P_values.operator[](i)  = P;
         nK_values.operator[](i) = I_n_K; 
     }
     std::clock_t seed_I;
     const std::size_t I_sym_len{I_n_samples*static_cast<std::size_t>(I_n_K)};
    
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n\n", __PRETTY_FUNCTION__);
     am_bb_sine_signal_sequence_t sin_sequence = am_bb_sine_signal_sequence_t(I_n_samples,
                                                                             nsignals,nK_values,
                                                                             n_values,A_values,P_values);
     am_bb_sine_signal_t symbol_1 = am_bb_sine_signal_t(I_n_samples,I_n_K,n,A,P);
     am_bb_sine_signal_t symbol_2 = am_bb_sine_signal_t(I_n_samples,I_n_K,n,A,P);
     am_bb_sine_signal_t symbol_3 = am_bb_sine_signal_t(I_n_samples,I_n_K,n,A,P);
     am_bb_sine_signal_t symbol_4 = am_bb_sine_signal_t(I_n_samples,I_n_K,n,A,P);
          
     I_sym = darray_r4_t(I_sym_len);
     sig_vsum_reference = darray_r4_t(I_n_samples);
     std::printf("[UNIT-TEST] -- Number of samples=%llu,number of narrowband signals=%d,number of baseband signals=%llu\n\n",
                            I_n_samples,I_n_K,nsignals);
     char * ctor_name{gms::common::demangle(typeid(sin_sequence).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(sin_sequence).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
    
     std::printf("[UNIT-TEST:] -- Generating random symbol data: I_channel(-1,1)\n\n");
     seed_I = std::clock();
     auto rand_sym_I{std::bind(std::uniform_real_distribution<float>(-1.0f,1.0f),std::mt19937(seed_I))};
     
     for(std::size_t __i{0ull}; __i != I_sym_len; ++__i) 
     {
          I_sym.m_data[__i] = rand_sym_I();
          //std::printf("%.7f\n",I_sym.m_data[__i]);
     }
     seed_func = std::clock();
     auto rand_func{std::bind(std::uniform_int_distribution<std::uint32_t>(0,1),std::mt19937(seed_func))};
     const std::uint32_t which_func{rand_func()}; 
     std::int32_t ret_stat = sin_sequence.signal_sequence_avx512_u6x(&I_sym.m_data[0],&ret_codes.operator[](0),which_func);
     if(ret_stat<0) { 
        std::printf("[UNIT-TEST:] -- **FAIL** -> signal_sequence_avx512_u6x stat=%d\n",ret_stat);
        return;
     } 
     std::printf("[UNIT-TEST:] -- Calling  cosine wave components modulated by the user data -- vertically summed sequence (REFERENCE!!).\n\n");
     switch (which_func) 
     {
     case 0 :
          symbol_1.create_signal_user_data(&I_sym.m_data[0],I_n_samples,I_n_K);
          symbol_2.create_signal_user_data(&I_sym.m_data[0],I_n_samples,I_n_K); 
          symbol_3.create_signal_user_data(&I_sym.m_data[0],I_n_samples,I_n_K);
          symbol_4.create_signal_user_data(&I_sym.m_data[0],I_n_samples,I_n_K);
     break;
     case 1 : 
          symbol_1.create_signal_user_data_u4x(&I_sym.m_data[0],I_n_samples,I_n_K);
          symbol_2.create_signal_user_data_u4x(&I_sym.m_data[0],I_n_samples,I_n_K); 
          symbol_3.create_signal_user_data_u4x(&I_sym.m_data[0],I_n_samples,I_n_K);
          symbol_4.create_signal_user_data_u4x(&I_sym.m_data[0],I_n_samples,I_n_K);
     default : 
           std::printf("wrong case argument=%d\n",which_func);
           return;
     }    
     signal_sequence_vertical_summation(&symbol_1.m_sig_samples.m_data[0],
                                        &symbol_2.m_sig_samples.m_data[0],
                                        &symbol_3.m_sig_samples.m_data[0],
                                        &symbol_4.m_sig_samples.m_data[0],
                                        &sig_vsum_reference.m_data[0],
                                        I_n_samples);
     check_fp_approx_equality(&sin_sequence.m_vsummed_sequence.operator[](0),
                             &sig_vsum_reference.m_data[0],I_n_samples);
     //std::cout << __cmplx_trapw_env__ << std::endl;
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n\n");
     am_bb_sine_signal_sequence_t::create_sequence_plot(sin_sequence.m_nsamples,&sin_sequence.m_vsummed_sequence.operator[](0),nullptr,
                                              "am_bb_sine_sequence_avx512_u6x_correctness_test_1_","AM-BB_sine_Signals_Sequence",false);
     am_bb_sine_signal_sequence_t::create_sequence_plot(I_n_samples,&sig_vsum_reference.m_data[0],nullptr,
                                               "am_bb_sine_signal_reference_vsumm_u6x_test_1_","AM-BB_sine_Signals_Sequence",false);
     
     printf("[UNIT_TEST]: function=%s -- **END**\n\n", __PRETTY_FUNCTION__);
}



__attribute__((hot))
__attribute__((noinline))
void unit_test_am_bb_sine_sequence_varying_params_avx512_u16x();

void unit_test_am_bb_sine_sequence_varying_params_avx512_u16x()
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
     std::clock_t seed_func{};
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
     am_bb_sine_signal_sequence_t sin_sequence = am_bb_sine_signal_sequence_t(I_n_samples,
                                                                             nsignals,nK_values,
                                                                             n_values,A_values,P_values);
     I_sym = darray_r4_t(I_sym_len);
     char * ctor_name{gms::common::demangle(typeid(sin_sequence).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(sin_sequence).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
   
     std::printf("[UNIT-TEST:] -- Generating random symbol data: I_channel(-1,1)\n");
     seed_I = std::clock();
     auto rand_sym_I{std::bind(std::uniform_real_distribution<float>(-1.0f,1.0f),std::mt19937(seed_I))};
     
     for(std::size_t __i{0ull}; __i != I_sym_len; ++__i) 
     {
          I_sym.m_data[__i] = rand_sym_I();
          //std::printf("%.7f\n",I_sym.m_data[__i]);
     }
     std::printf("[UNIT-TEST:] -- Calling  cosine wave components modulated by the user data -- vertically summed sequence.\n");
     seed_func = std::clock();
     auto rand_func{std::bind(std::uniform_int_distribution<std::uint32_t>(0,1),std::mt19937(seed_func))};
     const std::uint32_t which_func{rand_func()};
     std::int32_t ret_stat = sin_sequence.signal_sequence_avx512_u16x(&I_sym.m_data[0],&ret_codes.operator[](0),which_func);  
     if(ret_stat<0) { 
        std::printf("[UNIT-TEST:] -- **FAIL** -> signal_sequence_avx512_u16x stat=%d\n",ret_stat);
        return;
     }   
     for(std::size_t i{0}; i != sin_sequence.m_nsignals; ++i) 
     {
          const std::int32_t ret_code{ret_codes.operator[](i)};
          if(ret_code < 0ull) 
          {
               std::printf("[**ERROR**]: -- Found a failed call at index=%llu, of value=%d\n",i,ret_code);
          }
     }
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     am_bb_sine_signal_sequence_t::create_sequence_plot(sin_sequence.m_nsamples,&sin_sequence.m_vsummed_sequence.operator[](0),nullptr,
                                              "am_bb_sine_sequence_avx512_u16x_var_params_test_1_","AM-BB_sine_Signals_Sequence",false);
     
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}


__attribute__((hot))
__attribute__((noinline))
void unit_test_am_bb_sine_sequence_rand_data_avx512_u10x();

void unit_test_am_bb_sine_sequence_rand_data_avx512_u10x()
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
     std::clock_t seed_func{};
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
     am_bb_sine_signal_sequence_t sin_sequence = am_bb_sine_signal_sequence_t(I_n_samples,
                                                                             nsignals,nK_values,
                                                                             n_values,A_values,P_values);
     I_sym = darray_r4_t(I_sym_len);
     char * ctor_name{gms::common::demangle(typeid(sin_sequence).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(sin_sequence).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
     
     std::printf("[UNIT-TEST:] -- Generating random symbol data: I_channel(-1,1)\n");
     seed_I = std::clock();
     auto rand_sym_I{std::bind(std::uniform_real_distribution<float>(-1.0f,1.0f),std::mt19937(seed_I))};
     
     for(std::size_t __i{0ull}; __i != I_sym_len; ++__i) 
     {
          I_sym.m_data[__i] = rand_sym_I();
          //std::printf("%.7f\n",I_sym.m_data[__i]);
     }
     std::printf("[UNIT-TEST:] -- Calling  cosine wave components modulated by the user data -- vertically summed sequence.\n");
     seed_func = std::clock();
     auto rand_func{std::bind(std::uniform_int_distribution<std::uint32_t>(0,1),std::mt19937(seed_func))};
     const std::uint32_t which_func{rand_func()};
     std::int32_t ret_stat = sin_sequence.signal_sequence_avx512_u10x(&I_sym.m_data[0],&ret_codes.operator[](0),which_func); 
     if(ret_stat<0) { 
        std::printf("[UNIT-TEST:] -- **FAIL** -> signal_sequence_avx512_u10x stat=%d\n",ret_stat);
        return;
     }    
     //std::cout << __cmplx_trapw_env__ << std::endl;
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     am_bb_sine_signal_sequence_t::create_sequence_plot(sin_sequence.m_nsamples,&sin_sequence.m_vsummed_sequence.operator[](0),nullptr,
                                              "am_bb_sine_sequence_avx512_u10x_rand_data_test_1_","AM-BB_sine_Signals_Sequence",false);
     
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_am_bb_sine_sequence_rand_data_avx512_u6x();

void unit_test_am_bb_sine_sequence_rand_data_avx512_u6x()
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
     std::clock_t seed_func{};
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
     am_bb_sine_signal_sequence_t sin_sequence = am_bb_sine_signal_sequence_t(I_n_samples,
                                                                             nsignals,nK_values,
                                                                             n_values,A_values,P_values);
     I_sym = darray_r4_t(I_sym_len);
     char * ctor_name{gms::common::demangle(typeid(sin_sequence).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(sin_sequence).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
   
     std::printf("[UNIT-TEST:] -- Generating random symbol data: I_channel(-1,1)\n");
     seed_I = std::clock();
     auto rand_sym_I{std::bind(std::uniform_real_distribution<float>(-1.0f,1.0f),std::mt19937(seed_I))};
     
     for(std::size_t __i{0ull}; __i != I_sym_len; ++__i) 
     {
          I_sym.m_data[__i] = rand_sym_I();
          //std::printf("%.7f\n",I_sym.m_data[__i]);
     }
     std::printf("[UNIT-TEST:] -- Calling  sine wave components modulated by the user data -- vertically summed sequence.\n");
     seed_func = std::clock();
     auto rand_func{std::bind(std::uniform_int_distribution<std::uint32_t>(0,1),std::mt19937(seed_func))};
     const std::uint32_t which_func{rand_func()};
     std::int32_t ret_stat = sin_sequence.signal_sequence_avx512_u6x(&I_sym.m_data[0],&ret_codes.operator[](0),which_func);  
     if(ret_stat<0) { 
        std::printf("[UNIT-TEST:] -- **FAIL** -> signal_sequence_avx512_u6x stat=%d\n",ret_stat);
        return;
     }   
     //std::cout << __cmplx_trapw_env__ << std::endl;
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     am_bb_sine_signal_sequence_t::create_sequence_plot(sin_sequence.m_nsamples,&sin_sequence.m_vsummed_sequence.operator[](0),nullptr,
                                              "am_bb_sine_sequence_avx512_u6x_trand_data_test_1_","AM-BB_sine_Signals_Sequence",false);
     
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}


int main()
{
    unit_test_am_bb_sine_sequence_rand_data_avx512_u16x();
    unit_test_am_bb_sine_sequence_rand_data_avx512_u10x();
    unit_test_am_bb_sine_sequence_rand_data_avx512_u6x();
    unit_test_am_bb_sine_sequence_varying_params_avx512_u16x();
    unit_test_am_bb_sine_sequence_correctness_avx512_u16x();
    unit_test_am_bb_sine_sequence_correctness_avx512_u10x();
    unit_test_am_bb_sine_sequence_correctness_avx512_u6x();
    return 0;
}
