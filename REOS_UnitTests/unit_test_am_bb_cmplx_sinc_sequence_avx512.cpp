#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <cmath>
#include "GMS_dyn_array.h"
#include "GMS_sse_memset.h"
#include "GMS_am_bb_cmplx_sinc_sequence.h"

/*
   icpc -o unit_test_am_bb_cmplx_sinc_sequence_avx512 -O3 -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h GMS_half.h GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_indices.h GMS_am_bb_cmplx_sinc_signal.h GMS_am_bb_cmplx_sinc_signal.cpp GMS_am_bb_cmplx_sinc_sequence.h GMS_am_bb_cmplx_sinc_sequence.cpp unit_test_am_bb_cmplx_sinc_sequence_avx512.cpp
   ASM: 
   icpc -S --O3 fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h GMS_half.h GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_indices.h GMS_am_bb_cmplx_sinc_signal.h GMS_am_bb_cmplx_sinc_signal.cpp GMS_am_bb_cmplx_sinc_sequence.h GMS_am_bb_cmplx_sinc_sequence.cpp unit_test_am_bb_cmplx_sinc_sequence_avx512.cpp

*/



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
void unit_test_am_bb_cmplx_sinc_sequence_correctness_avx512_u16x();

void unit_test_am_bb_cmplx_sinc_sequence_correctness_avx512_u16x()
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
     darray_r4_t sig_vsum_reference_re;
     darray_r4_t sig_vsum_reference_im;
     std::vector<float> A_values          = std::vector<float>(nsignals);
     std::vector<float> n_values          = std::vector<float>(nsignals);
     std::vector<float> P_values          = std::vector<float>(nsignals);
     std::vector<std::uint32_t> nK_values = std::vector<std::uint32_t>(nsignals);
     std::vector<std::int32_t>  ret_codes = std::vector<std::int32_t>(nsignals);
     float rand_scales[4] = {0.1f,0.5f,0.7f,1.0f};
     am_bb_cmplx_sinc_signal_pdf_params_t pdf_params;
     pdf_params.cauchy_a_r   = 0.0f;
     pdf_params.cauchy_b_r   = 1.0f;
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
     am_bb_cmplx_sinc_signal_sequence_t sinc_sequence = am_bb_cmplx_sinc_signal_sequence_t(I_n_samples,
                                                                                         nsignals,nK_values,
                                                                                        n_values,A_values,P_values);
     am_bb_cmplx_sinc_signal_t symbol_1 = am_bb_cmplx_sinc_signal_t(I_n_samples,I_n_K,n,A,P);
     am_bb_cmplx_sinc_signal_t symbol_2 = am_bb_cmplx_sinc_signal_t(I_n_samples,I_n_K,n,A,P);
     am_bb_cmplx_sinc_signal_t symbol_3 = am_bb_cmplx_sinc_signal_t(I_n_samples,I_n_K,n,A,P);
     am_bb_cmplx_sinc_signal_t symbol_4 = am_bb_cmplx_sinc_signal_t(I_n_samples,I_n_K,n,A,P);
          
     I_sym = darray_r4_t(I_sym_len);
     sig_vsum_reference_re = darray_r4_t(I_n_samples/2ull);
     sig_vsum_reference_im = darray_r4_t(I_n_samples/2ull);
     
     std::printf("[UNIT-TEST] -- Number of samples=%llu,number of narrowband signals=%d,number of baseband signals=%llu\n\n",
                            I_n_samples,I_n_K,nsignals);
     char * ctor_name{gms::common::demangle(typeid(sinc_sequence).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(sinc_sequence).name());
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
     for(std::uint32_t i_case{0u}; i_case != 4u; ++i_case) 
     {
         std::printf("[UNIT-TEST:] -- Calling: signal_sequence_avx512_u16x -- case=%d\n\n",i_case);  
         std::int32_t ret_stat = sinc_sequence.signal_sequence_avx512_u16x(&I_sym.m_data[0],pdf_params,am_bb_cmplx_sinc_signal_rand_distr::cauchy,&rand_scales[0],&ret_codes.operator[](0),i_case);
         if(ret_stat<0) { 
            std::printf("[UNIT-TEST:] -- **FAIL** -> signal_sequence_avx512_u16x stat=%d\n",ret_stat);
            return;
          }
         std::printf("[UNIT-TEST:] -- Calling  sinc wave components modulated by the user data -- vertically summed sequence (REFERENCE!!).\n\n");
         symbol_1.create_signal_user_data(&I_sym.m_data[0],I_n_samples,I_n_K);
         symbol_2.create_signal_user_data(&I_sym.m_data[0],I_n_samples,I_n_K); 
         symbol_3.create_signal_user_data(&I_sym.m_data[0],I_n_samples,I_n_K);
         symbol_4.create_signal_user_data(&I_sym.m_data[0],I_n_samples,I_n_K);
         if(i_case==0u)
         {
            std::printf("[UNIT-TEST:] -- case=%d, sanity-check!!\n",i_case);
            signal_sequence_vertical_summation(&symbol_1.get_I_channel()[0],
                                        &symbol_2.get_I_channel()[0],
                                        &symbol_3.get_I_channel()[0],
                                        &symbol_4.get_I_channel()[0],
                                        &sig_vsum_reference_re.m_data[0],
                                        I_n_samples/2ull);
            signal_sequence_vertical_summation(&symbol_1.get_Q_channel()[0],
                                        &symbol_2.get_Q_channel()[0],
                                        &symbol_3.get_Q_channel()[0],
                                        &symbol_4.get_Q_channel()[0],
                                        &sig_vsum_reference_im.m_data[0],
                                        I_n_samples/2ull);                                   
            check_fp_approx_equality(&sinc_sequence.m_re_sequence.m_data[0],
                             &sig_vsum_reference_re.m_data[0],I_n_samples/2ull);
            check_fp_approx_equality(&sinc_sequence.m_im_sequence.m_data[0],
                             &sig_vsum_reference_im.m_data[0],I_n_samples/2ull);                        
         }
         if(i_case>0u && i_case<4)
         { 
            std::printf("[UNIT-TEST:] -- case=%d, plotting results\n",i_case);
            std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n\n");
            std::string s1 = "am_bb_cmplx_sinc_seq_avx512_u16x_re_correct_test_1_case_"+std::to_string(i_case);
            std::string s2 = "am_bb_cmplx_sinc_signal_reference_vsumm_re_u16x_test_1_"+std::to_string(i_case);
         
            am_bb_cmplx_sinc_signal_sequence_t::create_sequence_plot(sinc_sequence.m_nsamples,&sinc_sequence.m_re_sequence.m_data[0],nullptr,
                                              s1.c_str(),"AM-BB_sinc_signals_sequence",false);
            am_bb_cmplx_sinc_signal_sequence_t::create_sequence_plot(I_n_samples/2ull,&sig_vsum_reference_re.m_data[0],nullptr,
                                               s2.c_str(),"AM-BB_sinc_signals_sequence",false);
            std::string s3 = "am_bb_cmplx_sinc_seq_avx512_u16x_im_correct_test_1_"+std::to_string(i_case);
            std::string s4 = "am_bb_cmplx_sinc_signal_reference_vsumm_im_u16x_test_1_"+std::to_string(i_case);
            am_bb_cmplx_sinc_signal_sequence_t::create_sequence_plot(sinc_sequence.m_nsamples,&sinc_sequence.m_im_sequence.m_data[0],nullptr,
                                              s3.c_str(),"AM-BB_sinc_signals_sequence",false);
            am_bb_cmplx_sinc_signal_sequence_t::create_sequence_plot(I_n_samples/2ull,&sig_vsum_reference_im.m_data[0],nullptr,
                                              s4.c_str() ,"AM-BB_sinc_signals_sequence",false);  
         }
         gms::common::sse_memset_unroll8x_ps(&sinc_sequence.m_re_sequence.m_data[0],0.0f,I_n_samples/2ull);
         gms::common::sse_memset_unroll8x_ps(&sig_vsum_reference_re.m_data[0],0.0f,I_n_samples/2ull);
         gms::common::sse_memset_unroll8x_ps(&sinc_sequence.m_im_sequence.m_data[0],0.0f,I_n_samples/2ull);
         gms::common::sse_memset_unroll8x_ps(&sig_vsum_reference_im.m_data[0],0.0f,I_n_samples/2ull);  
         std::memset(&ret_codes[0],0,sinc_sequence.m_nsignals);
     }
     printf("[UNIT_TEST]: function=%s -- **END**\n\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_am_bb_cmplx_sinc_sequence_correctness_avx512_u10x();

void unit_test_am_bb_cmplx_sinc_sequence_correctness_avx512_u10x()
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
     darray_r4_t sig_vsum_reference_re;
     darray_r4_t sig_vsum_reference_im;
     std::vector<float> A_values          = std::vector<float>(nsignals);
     std::vector<float> n_values          = std::vector<float>(nsignals);
     std::vector<float> P_values          = std::vector<float>(nsignals);
     std::vector<std::uint32_t> nK_values = std::vector<std::uint32_t>(nsignals);
     std::vector<std::int32_t>  ret_codes = std::vector<std::int32_t>(nsignals);
     float rand_scales[4] = {0.1f,0.5f,0.7f,1.0f};
     am_bb_cmplx_sinc_signal_pdf_params_t pdf_params;
     pdf_params.cauchy_a_r   = 0.0f;
     pdf_params.cauchy_b_r   = 1.0f;
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
     am_bb_cmplx_sinc_signal_sequence_t sinc_sequence = am_bb_cmplx_sinc_signal_sequence_t(I_n_samples,
                                                                                         nsignals,nK_values,
                                                                                        n_values,A_values,P_values);
     am_bb_cmplx_sinc_signal_t symbol_1 = am_bb_cmplx_sinc_signal_t(I_n_samples,I_n_K,n,A,P);
     am_bb_cmplx_sinc_signal_t symbol_2 = am_bb_cmplx_sinc_signal_t(I_n_samples,I_n_K,n,A,P);
     am_bb_cmplx_sinc_signal_t symbol_3 = am_bb_cmplx_sinc_signal_t(I_n_samples,I_n_K,n,A,P);
     am_bb_cmplx_sinc_signal_t symbol_4 = am_bb_cmplx_sinc_signal_t(I_n_samples,I_n_K,n,A,P);
          
     I_sym = darray_r4_t(I_sym_len);
     sig_vsum_reference_re = darray_r4_t(I_n_samples/2ull);
     sig_vsum_reference_im = darray_r4_t(I_n_samples/2ull);
     
     std::printf("[UNIT-TEST] -- Number of samples=%llu,number of narrowband signals=%d,number of baseband signals=%llu\n\n",
                            I_n_samples,I_n_K,nsignals);
     char * ctor_name{gms::common::demangle(typeid(sinc_sequence).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(sinc_sequence).name());
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
     for(std::uint32_t i_case{0u}; i_case != 4u; ++i_case) 
     {
         std::printf("[UNIT-TEST:] -- Calling: signal_sequence_avx512_u10x -- case=%d\n\n",i_case);  
         std::int32_t ret_stat = sinc_sequence.signal_sequence_avx512_u10x(&I_sym.m_data[0],pdf_params,am_bb_cmplx_sinc_signal_rand_distr::cauchy,&rand_scales[0],&ret_codes.operator[](0),i_case);
         if(ret_stat<0) { 
            std::printf("[UNIT-TEST:] -- **FAIL** -> signal_sequence_avx_u10x stat=%d\n",ret_stat);
            return;
          }
         std::printf("[UNIT-TEST:] -- Calling  sinc wave components modulated by the user data -- vertically summed sequence (REFERENCE!!).\n\n");
         symbol_1.create_signal_user_data(&I_sym.m_data[0],I_n_samples,I_n_K);
         symbol_2.create_signal_user_data(&I_sym.m_data[0],I_n_samples,I_n_K); 
         symbol_3.create_signal_user_data(&I_sym.m_data[0],I_n_samples,I_n_K);
         symbol_4.create_signal_user_data(&I_sym.m_data[0],I_n_samples,I_n_K);
         if(i_case==0u)
         {
            std::printf("[UNIT-TEST:] -- case=%d, sanity-check!!\n",i_case);
            signal_sequence_vertical_summation(&symbol_1.get_I_channel()[0],
                                        &symbol_2.get_I_channel()[0],
                                        &symbol_3.get_I_channel()[0],
                                        &symbol_4.get_I_channel()[0],
                                        &sig_vsum_reference_re.m_data[0],
                                        I_n_samples/2ull);
            signal_sequence_vertical_summation(&symbol_1.get_Q_channel()[0],
                                        &symbol_2.get_Q_channel()[0],
                                        &symbol_3.get_Q_channel()[0],
                                        &symbol_4.get_Q_channel()[0],
                                        &sig_vsum_reference_im.m_data[0],
                                        I_n_samples/2ull);                                   
            check_fp_approx_equality(&sinc_sequence.m_re_sequence.m_data[0],
                             &sig_vsum_reference_re.m_data[0],I_n_samples/2ull);
            check_fp_approx_equality(&sinc_sequence.m_im_sequence.m_data[0],
                             &sig_vsum_reference_im.m_data[0],I_n_samples/2ull);                        
         }
         if(i_case>0u && i_case<4)
         { 
            std::printf("[UNIT-TEST:] -- case=%d, plotting results\n",i_case);
            std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n\n");
            std::string s1 = "am_bb_cmplx_sinc_seq_avx512_u10x_re_correct_test_1_case_"+std::to_string(i_case);
            std::string s2 = "am_bb_cmplx_sinc_signal_reference_vsumm_re_u10x_test_1_"+std::to_string(i_case);
         
            am_bb_cmplx_sinc_signal_sequence_t::create_sequence_plot(sinc_sequence.m_nsamples,&sinc_sequence.m_re_sequence.m_data[0],nullptr,
                                              s1.c_str(),"AM-BB_sinc_signals_sequence",false);
            am_bb_cmplx_sinc_signal_sequence_t::create_sequence_plot(I_n_samples/2ull,&sig_vsum_reference_re.m_data[0],nullptr,
                                               s2.c_str(),"AM-BB_sinc_signals_sequence",false);
            std::string s3 = "am_bb_cmplx_sinc_seq_avx512_u10x_im_correct_test_1_"+std::to_string(i_case);
            std::string s4 = "am_bb_cmplx_sinc_signal_reference_vsumm_im_u10x_test_1_"+std::to_string(i_case);
            am_bb_cmplx_sinc_signal_sequence_t::create_sequence_plot(sinc_sequence.m_nsamples,&sinc_sequence.m_im_sequence.m_data[0],nullptr,
                                              s3.c_str(),"AM-BB_sinc_signals_sequence",false);
            am_bb_cmplx_sinc_signal_sequence_t::create_sequence_plot(I_n_samples/2ull,&sig_vsum_reference_im.m_data[0],nullptr,
                                              s4.c_str() ,"AM-BB_sinc_signals_sequence",false);  
         }
         gms::common::sse_memset_unroll8x_ps(&sinc_sequence.m_re_sequence.m_data[0],0.0f,I_n_samples/2ull);
         gms::common::sse_memset_unroll8x_ps(&sig_vsum_reference_re.m_data[0],0.0f,I_n_samples/2ull);
         gms::common::sse_memset_unroll8x_ps(&sinc_sequence.m_im_sequence.m_data[0],0.0f,I_n_samples/2ull);
         gms::common::sse_memset_unroll8x_ps(&sig_vsum_reference_im.m_data[0],0.0f,I_n_samples/2ull);  
         std::memset(&ret_codes[0],0,sinc_sequence.m_nsignals);
     }
     printf("[UNIT_TEST]: function=%s -- **END**\n\n", __PRETTY_FUNCTION__);
}


__attribute__((hot))
__attribute__((noinline))
void unit_test_am_bb_cmplx_sinc_sequence_correctness_avx512_u6x();

void unit_test_am_bb_cmplx_sinc_sequence_correctness_avx512_u6x()
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
     darray_r4_t sig_vsum_reference_re;
     darray_r4_t sig_vsum_reference_im;
     std::vector<float> A_values          = std::vector<float>(nsignals);
     std::vector<float> n_values          = std::vector<float>(nsignals);
     std::vector<float> P_values          = std::vector<float>(nsignals);
     std::vector<std::uint32_t> nK_values = std::vector<std::uint32_t>(nsignals);
     std::vector<std::int32_t>  ret_codes = std::vector<std::int32_t>(nsignals);
     float rand_scales[4] = {0.1f,0.5f,0.7f,1.0f};
     am_bb_cmplx_sinc_signal_pdf_params_t pdf_params;
     pdf_params.cauchy_a_r   = 0.0f;
     pdf_params.cauchy_b_r   = 1.0f;
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
     am_bb_cmplx_sinc_signal_sequence_t sinc_sequence = am_bb_cmplx_sinc_signal_sequence_t(I_n_samples,
                                                                                         nsignals,nK_values,
                                                                                        n_values,A_values,P_values);
     am_bb_cmplx_sinc_signal_t symbol_1 = am_bb_cmplx_sinc_signal_t(I_n_samples,I_n_K,n,A,P);
     am_bb_cmplx_sinc_signal_t symbol_2 = am_bb_cmplx_sinc_signal_t(I_n_samples,I_n_K,n,A,P);
     am_bb_cmplx_sinc_signal_t symbol_3 = am_bb_cmplx_sinc_signal_t(I_n_samples,I_n_K,n,A,P);
     am_bb_cmplx_sinc_signal_t symbol_4 = am_bb_cmplx_sinc_signal_t(I_n_samples,I_n_K,n,A,P);
          
     I_sym = darray_r4_t(I_sym_len);
     sig_vsum_reference_re = darray_r4_t(I_n_samples/2ull);
     sig_vsum_reference_im = darray_r4_t(I_n_samples/2ull);
     
     std::printf("[UNIT-TEST] -- Number of samples=%llu,number of narrowband signals=%d,number of baseband signals=%llu\n\n",
                            I_n_samples,I_n_K,nsignals);
     char * ctor_name{gms::common::demangle(typeid(sinc_sequence).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(sinc_sequence).name());
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
     for(std::uint32_t i_case{0u}; i_case != 4u; ++i_case) 
     {
         std::printf("[UNIT-TEST:] -- Calling: signal_sequence_avx512_u6x -- case=%d\n\n",i_case);  
         std::int32_t ret_stat = sinc_sequence.signal_sequence_avx512_u6x(&I_sym.m_data[0],pdf_params,am_bb_cmplx_sinc_signal_rand_distr::cauchy,&rand_scales[0],&ret_codes.operator[](0),i_case);
         if(ret_stat<0) { 
            std::printf("[UNIT-TEST:] -- **FAIL** -> signal_sequence_avx512_u6x stat=%d\n",ret_stat);
            return;
          }
         std::printf("[UNIT-TEST:] -- Calling  sinc wave components modulated by the user data -- vertically summed sequence (REFERENCE!!).\n\n");
         symbol_1.create_signal_user_data(&I_sym.m_data[0],I_n_samples,I_n_K);
         symbol_2.create_signal_user_data(&I_sym.m_data[0],I_n_samples,I_n_K); 
         symbol_3.create_signal_user_data(&I_sym.m_data[0],I_n_samples,I_n_K);
         symbol_4.create_signal_user_data(&I_sym.m_data[0],I_n_samples,I_n_K);
         if(i_case==0u)
         {
            std::printf("[UNIT-TEST:] -- case=%d, sanity-check!!\n",i_case);
            signal_sequence_vertical_summation(&symbol_1.get_I_channel()[0],
                                        &symbol_2.get_I_channel()[0],
                                        &symbol_3.get_I_channel()[0],
                                        &symbol_4.get_I_channel()[0],
                                        &sig_vsum_reference_re.m_data[0],
                                        I_n_samples/2ull);
            signal_sequence_vertical_summation(&symbol_1.get_Q_channel()[0],
                                        &symbol_2.get_Q_channel()[0],
                                        &symbol_3.get_Q_channel()[0],
                                        &symbol_4.get_Q_channel()[0],
                                        &sig_vsum_reference_im.m_data[0],
                                        I_n_samples/2ull);                                   
            check_fp_approx_equality(&sinc_sequence.m_re_sequence.m_data[0],
                             &sig_vsum_reference_re.m_data[0],I_n_samples/2ull);
            check_fp_approx_equality(&sinc_sequence.m_im_sequence.m_data[0],
                             &sig_vsum_reference_im.m_data[0],I_n_samples/2ull);                        
         }
         if(i_case>0u && i_case<4)
         { 
            std::printf("[UNIT-TEST:] -- case=%d, plotting results\n",i_case);
            std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n\n");
            std::string s1 = "am_bb_cmplx_sinc_seq_avx512_u6x_re_correct_test_1_case_"+std::to_string(i_case);
            std::string s2 = "am_bb_cmplx_sinc_signal_reference_vsumm_re_u6x_test_1_"+std::to_string(i_case);
         
            am_bb_cmplx_sinc_signal_sequence_t::create_sequence_plot(sinc_sequence.m_nsamples,&sinc_sequence.m_re_sequence.m_data[0],nullptr,
                                              s1.c_str(),"AM-BB_sinc_signals_sequence",false);
            am_bb_cmplx_sinc_signal_sequence_t::create_sequence_plot(I_n_samples/2ull,&sig_vsum_reference_re.m_data[0],nullptr,
                                               s2.c_str(),"AM-BB_sinc_signals_sequence",false);
            std::string s3 = "am_bb_cmplx_sinc_seq_avx512_u6x_im_correct_test_1_"+std::to_string(i_case);
            std::string s4 = "am_bb_cmplx_sinc_signal_reference_vsumm_im_u6x_test_1_"+std::to_string(i_case);
            am_bb_cmplx_sinc_signal_sequence_t::create_sequence_plot(sinc_sequence.m_nsamples,&sinc_sequence.m_im_sequence.m_data[0],nullptr,
                                              s3.c_str(),"AM-BB_sinc_signals_sequence",false);
            am_bb_cmplx_sinc_signal_sequence_t::create_sequence_plot(I_n_samples/2ull,&sig_vsum_reference_im.m_data[0],nullptr,
                                              s4.c_str() ,"AM-BB_sinc_signals_sequence",false);  
         }
         gms::common::sse_memset_unroll8x_ps(&sinc_sequence.m_re_sequence.m_data[0],0.0f,I_n_samples/2ull);
         gms::common::sse_memset_unroll8x_ps(&sig_vsum_reference_re.m_data[0],0.0f,I_n_samples/2ull);
         gms::common::sse_memset_unroll8x_ps(&sinc_sequence.m_im_sequence.m_data[0],0.0f,I_n_samples/2ull);
         gms::common::sse_memset_unroll8x_ps(&sig_vsum_reference_im.m_data[0],0.0f,I_n_samples/2ull);  
         std::memset(&ret_codes[0],0,sinc_sequence.m_nsignals);
     }
     printf("[UNIT_TEST]: function=%s -- **END**\n\n", __PRETTY_FUNCTION__);
}



int main()
{
     unit_test_am_bb_cmplx_sinc_sequence_correctness_avx512_u16x();
     unit_test_am_bb_cmplx_sinc_sequence_correctness_avx512_u10x();
     unit_test_am_bb_cmplx_sinc_sequence_correctness_avx512_u6x();
     return 0;
}