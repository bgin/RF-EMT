/* Copyright (C) Bernard Gingold, 2020-2026 
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <sched.h>
#include <cerrno>
#include <sys/types.h>
#include <sys/resource.h> //getpid 
#include <signal.h>
#include <unistd.h>
#include <immintrin.h>
#include <omp.h> 
#include <algorithm>
#include <random>
#include <ctime>
#include <functional>
#include "GMS_malloc.h"
#include "GMS_ansi_colors.h"
#include "GMS_simd_cephes_cos.h"
#include "GMS_vec_cephes_cos_ps.h"

/*
   icpc -o perf_test_select_vs_masked_move_32kiB -O3 -fp-model fast=2 -fno-exceptions -fopenmp  -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_ansi_colors.h GMS_malloc.h GMS_vec_cephes_cos_ps.h GMS_simd_cephes_cos.h perf_test_select_vs_masked_move_32kiB.cpp
   ASM: 
   icpc -S -O3 -fverbose-asm -masm=intel  -fno-exceptions -fopenmp  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_ansi_colors.h GMS_malloc.h GMS_vec_cephes_cos_ps.h GMS_simd_cephes_cos.h perf_test_select_vs_masked_move_32kiB.cpp

*/

#define BUFFER_STORE_SIZE 80
#define FORMAT_STORE_SIZE 80 

// More accurate for the floating-point comparison.
static inline
bool 
almostEqual2sComplement(float A,float B, std::int32_t maxUlps) 
{
  // Make sure maxUlps is non-negative and small enough that the
  // default NAN won't compare as equal to anything.
  // This check disabled for now so that my tests run properly.
  // This assert should be enabled for normal use.
  // assert(maxUlps > 0 && maxUlps < 4 * 1024 * 1024);
  int aInt =  *reinterpret_cast<std::int32_t*>(&A);   //*(int *)&A;
  // Make aInt lexicographically ordered as a twos-complement int
  if (aInt < 0)
    aInt = 0x80000000 - aInt;
  // Make bInt lexicographically ordered as a twos-complement int
  int bInt =   *reinterpret_cast<std::int32_t*>(&B);  //*(int *)&B;
  if (bInt < 0)
    bInt = 0x80000000 - bInt;
  int intDiff = std::abs(aInt - bInt);
  if (intDiff <= maxUlps)
    return true;
  return false;
}

static inline 
void 
printNumber(const char * val_name,float f, const std::int32_t offset) 
{
  (*((int *)&f)) += offset;
  std::printf("value: %s,%+1.11g,0x%08X,%d\n",val_name,f, *(int *)&f, *(int *)&f);
}

__attribute__((hot))
int32_t set_affinity_and_priority(const int32_t,const int32_t);
// 0 -- success, 1 and  2 failure.
int32_t set_affinity_and_priority(const int32_t cpu, const int32_t priority)
{
    cpu_set_t cpu_set;
    sched_param sp;
    int32_t status{-1};
    CPU_ZERO(&cpu_set);
    CPU_SET(cpu,&cpu_set);
    if(sched_setaffinity(0,sizeof(cpu_set), &cpu_set) < 0) 
    {
         status = 1;
         return status;
    }
    printf("Affinity set to cpu: %d\n",cpu);
    __builtin_memset(&sp,0,sizeof(sp));
    sp.sched_priority = priority; //99
    if((sched_setscheduler(0,SCHED_FIFO,&sp)) == -1) 
    {
        status = 2;
        return status;
    }
    status = 0;
    return status;
}

__attribute__((hot))
void print_thread_affinity();

void print_thread_affinity()
{
     char default_format[FORMAT_STORE_SIZE];
     char format_specifier[] = "host=%20H tid=%0.4n binds_to=%A";
     char buffer[BUFFER_STORE_SIZE];
     std::size_t nchars{};
     std::size_t diff;
     nchars = omp_get_affinity_format(default_format,(std::size_t)FORMAT_STORE_SIZE);
     diff   = nchars-(std::size_t)FORMAT_STORE_SIZE;
     if(diff>0ull)
        nchars += diff;
     omp_set_affinity_format(format_specifier);
     nchars = omp_capture_affinity(&buffer[0],(std::size_t)BUFFER_STORE_SIZE,NULL);
     printf("tid=%d affinity:%s\n",omp_get_thread_num(),buffer);
}

__attribute__((hot))
__attribute__((noinline))
void perf_test_mm512_cephes_cos_ps_v3_32kiB(const float * __restrict__,
                                        float * __restrict__,
                                        unsigned __int64 * __restrict__ ,
                                        unsigned __int64 * __restrict__ ,
                                        unsigned __int64 * __restrict__ ,
                                        std::uint32_t    * __restrict__ ,
                                        std::uint32_t    * __restrict__ ,
                                        const int32_t ,
                                        const int32_t,
                                        uint32_t &);

__attribute__((hot))
__attribute__((noinline))
void perf_test_mm512_cephes_cos_ps_v3_32kiB(const float * __restrict__ cos_args,
                                        float * __restrict__ cos_values,
                                        unsigned __int64 * __restrict__ cephes_cos_s,
                                        unsigned __int64 * __restrict__ cephes_cos_e,
                                        unsigned __int64 * __restrict__ cephes_cos_d,
                                        std::uint32_t    * __restrict__ cephes_cos_tsc_aux_s,
                                        std::uint32_t    * __restrict__ cephes_cos_tsc_aux_e,
                                        const int32_t n_runs,
                                        const int32_t n_samples,
                                        uint32_t & tid) 
{
    tid = 9999;
/*
    Latency and Throughput
    Architecture	Latency	Throughput (CPI)
    Skylake	          42	     5.5
*/
    constexpr std::size_t RDTSCP_LAT{42ull};
    constexpr std::int32_t zmm_size{16};
    [[maybe_unused]]
    volatile unsigned __int64 warmup_start;
    [[maybe_unused]]
    volatile unsigned __int64 warmup_end;
    [[maybe_unused]]
    volatile __m512 warmup_cos_values;
    //register __m512 result_cos_values;
    const float * __restrict__ p_cos_args   = cos_args;
    float       * __restrict__ p_cos_values = cos_values; 
    uint32_t mem_start{0};
    uint32_t mem_end{0}; 
    // warmup of RDTSCP
    warmup_start = __rdtscp(&mem_start);
    warmup_end   = __rdtscp(&mem_end);
    // warmup call of _mm512_cephes_cosf_ps_v2
    __m512 a_cos_args = _mm512_load_ps((const float*)&p_cos_args[0]); // load first 64-bytes
    warmup_cos_values = gms::math::_mm512_ceph_cosf_ps_v3<true>(a_cos_args);

    // Main testing loop begins.
    for(int32_t __i{0}; __i != n_runs; ++__i)
    {
        for(int32_t __j{0}; __j != n_samples; ++__j) 
        {   
            register __m512 data_in = _mm512_load_ps((const float*)&p_cos_args[(__i*n_samples+__j)*zmm_size]);
            __asm__ __volatile__ ("lfence");
            unsigned __int64 start_curr{__rdtscp(&mem_start)};       
            register __m512 result_cos_values = gms::math::_mm512_ceph_cosf_ps_v3<false>(data_in);
            unsigned __int64 end_curr{__rdtscp(&mem_end)};
            __asm__ __volatile__ ("lfence");
            //remove latency
            unsigned __int64 start_corrected{start_curr-RDTSCP_LAT};
            cephes_cos_s[__i*n_samples+__j] = start_corrected;
            unsigned __int64 end_corrected{end_curr-RDTSCP_LAT};
            cephes_cos_e[__i*n_samples+__j]   = end_corrected;
            cephes_cos_d[__i*n_samples+__j] = end_corrected-start_corrected;
            cephes_cos_tsc_aux_s[__i*n_samples+__j] = mem_start;
            cephes_cos_tsc_aux_e[__i*n_samples+__j] = mem_end;
            _mm512_store_ps((float*)&p_cos_values[(__i*n_samples+__j)*zmm_size],result_cos_values);
        }
    }
    tid = omp_get_thread_num();
}


__attribute__((hot))
__attribute__((noinline))
void perf_test_simd_ceph_cosf_16xf32_32kiB(const float * __restrict__ ,
                                    float * __restrict__ ,
                                    unsigned __int64 * __restrict__ ,
                                    unsigned __int64 * __restrict__ ,
                                    unsigned __int64 * __restrict__ ,
                                    std::uint32_t    * __restrict__ ,
                                    std::uint32_t    * __restrict__ ,
                                    const int32_t ,
                                    const int32_t,
                                    uint32_t &);

__attribute__((hot))
__attribute__((noinline))
void perf_test_simd_ceph_cosf_16xf32_32kiB(const float * __restrict__ cos_args,
                                    float * __restrict__ cos_values,
                                    unsigned __int64 * __restrict__ svml_cos_s,
                                    unsigned __int64 * __restrict__ svml_cos_e,
                                    unsigned __int64 * __restrict__ svml_cos_d,
                                    std::uint32_t    * __restrict__ svml_cos_tsc_aux_s,
                                    std::uint32_t    * __restrict__ svml_cos_tsc_aux_e,
                                    const int32_t n_runs,
                                    const int32_t n_samples,
                                    uint32_t & tid) 
{
    tid = 9999;
/*
    Latency and Throughput
    Architecture	Latency	Throughput (CPI)
    Skylake	          42	     5.5
*/
    constexpr std::size_t RDTSCP_LAT{42ull};
    constexpr std::int32_t zmm_size{16};
    [[maybe_unused]]
    volatile unsigned __int64 warmup_start;
    [[maybe_unused]]
    volatile unsigned __int64 warmup_end;
    [[maybe_unused]]
    volatile __m512 warmup_cos_values;
    const float * __restrict__ p_cos_args   = cos_args;
    float       * __restrict__ p_cos_values = cos_values;
    uint32_t mem_start{0};
    uint32_t mem_end{0}; 
    // warmup of RDTSCP
    warmup_start = __rdtscp(&mem_start);
    warmup_end   = __rdtscp(&mem_end);
    // warmup call of _mm512_cos_ps
    __m512 a_cos_args = _mm512_load_ps((const float*)&p_cos_args[0]);
    warmup_cos_values = gms::math::simd_ceph_cosf_16xf32(a_cos_args);

    // Main testing loop begins.
    for(int32_t __i{0}; __i != n_runs; ++__i)
    {
        for(int32_t __j{0}; __j != n_samples; ++__j) 
        {   
            register __m512 data_in = _mm512_load_ps((const float*)&p_cos_args[(__i*n_samples+__j)*zmm_size]);
            __asm__ __volatile__ ("lfence");
            unsigned __int64 start_curr{__rdtscp(&mem_start)};
            register __m512 result_cos_values = gms::math::simd_ceph_cosf_16xf32(data_in);
            unsigned __int64 end_curr{__rdtscp(&mem_end)};
            __asm__ __volatile__ ("lfence");
            //remove latency
            unsigned __int64 start_corrected{start_curr-RDTSCP_LAT};
            svml_cos_s[__i*n_samples+__j] = start_corrected;
            unsigned __int64 end_corrected{end_curr-RDTSCP_LAT};
            svml_cos_e[__i*n_samples+__j]   = end_corrected;
            svml_cos_d[__i*n_samples+__j] = end_corrected-start_corrected;
            svml_cos_tsc_aux_s[__i*n_samples+__j] = mem_start;
            svml_cos_tsc_aux_e[__i*n_samples+__j] = mem_end;
            _mm512_store_ps((float*)&p_cos_values[(__i*n_samples+__j)*zmm_size],result_cos_values);
        }
    }
    tid = omp_get_thread_num(); 
}


__attribute__((hot))
void test_runner_omp_sections_1st_seq();

void test_runner_omp_sections_1st_seq()
{
    using namespace gms::common;
    constexpr int32_t n_runs{10};
    constexpr int32_t n_samples{50};
    constexpr int32_t n_total{n_runs*n_samples};
    constexpr int32_t sz{n_total*16};
    constexpr std::size_t sz_bytes{static_cast<std::size_t>(sz)*sizeof(float)};
    unsigned __int64 ceph_cosf_16xf32_s[n_total]   = {UINT64_MAX};
    unsigned __int64 ceph_cosf_16xf32_e[n_total]   = {UINT64_MAX}; 
    unsigned __int64 ceph_cosf_16xf32_d[n_total]   = {UINT64_MAX};
    unsigned __int64 cephes_cos_ps_s[n_total]      = {UINT64_MAX};
    unsigned __int64 cephes_cos_ps_e[n_total]      = {UINT64_MAX};
    unsigned __int64 cephes_cos_ps_d[n_total]      = {UINT64_MAX};
    std::uint32_t    ceph_cosf_16xf32_tsc_aux_s[n_total] = {UINT32_MAX};
    std::uint32_t    ceph_cosf_16xf32_tsc_aux_e[n_total]   = {UINT32_MAX};
    std::uint32_t    cephes_cos_tsc_aux_s[n_total] = {UINT32_MAX};
    std::uint32_t    cephes_cos_tsc_aux_e[n_total] = {UINT32_MAX};
    float       * __restrict__ cos_data_in         = NULL;
    float       * __restrict__ cephes_cos_results  = NULL;
    float       * __restrict__ ceph_cosf_16xf32_results   = NULL;
    cos_data_in        = reinterpret_cast<float * __restrict__>(gms_mm_malloc(sz_bytes,64ull));
    cephes_cos_results = reinterpret_cast<float * __restrict__>(gms_mm_malloc(sz_bytes,64ull));
    ceph_cosf_16xf32_results   = reinterpret_cast<float * __restrict__>(gms_mm_malloc(sz_bytes,64ull));
    thread_local std::uniform_real_distribution<float> cos_rand_in;
    thread_local std::uniform_int_distribution<std::int32_t> rand_rad_range;
    thread_local std::uint64_t seed{};
    thread_local std::uint64_t seed_range{};
    float lo;
    float hi;
    seed_range = __rdtsc();
    rand_rad_range = std::uniform_int_distribution<std::int32_t>(0,1);
    auto rand_rad_range_gen{std::mt19937(seed_range)};
    const std::int32_t which_range{rand_rad_range.operator()(rand_rad_range_gen)};
    if(which_range==1)
    {
        lo = -3.14159265358979324f;
        hi = +3.14159265358979324f;
        std::printf("[UNIT_TEST]: Input range: lo=%.7f(rad),hi=%.7f(rad)\n",lo,hi);
    }
    else if(which_range==0)
    {
        lo = -3.14159265358979324f*8192.0f;
        hi = +3.14159265358979324f*8192.0f;
        std::printf("[UNIT_TEST]: Input range: lo=%.7f(rad),hi=%.7f(rad)\n",lo,hi);
    }
    seed = __rdtsc();
    cos_rand_in = std::uniform_real_distribution<float>(lo,hi);
    auto rand_uni_gen{std::mt19937(seed)};
    for(std::size_t i = 0ull; i != ROUND_TO_EIGHT(sz,7); i += 8ull)
    {
         const float in_cos_arg_1{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+0ull] = in_cos_arg_1;
         const float in_cos_arg_2{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+1ull] = in_cos_arg_2;
         const float in_cos_arg_3{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+2ull] = in_cos_arg_3;
         const float in_cos_arg_4{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+3ull] = in_cos_arg_4;
         const float in_cos_arg_5{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+4ull] = in_cos_arg_5;
         const float in_cos_arg_6{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+5ull] = in_cos_arg_6;
         const float in_cos_arg_7{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+6ull] = in_cos_arg_7;
         const float in_cos_arg_8{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+7ull] = in_cos_arg_8;
    }
    
    uint32_t ceph_cosf_16xf32_tid;
    uint32_t cephes_cos_tid;
    int32_t setenv_ret;
    setenv_ret = setenv("OMP_PROC_BIND","true",1);
    if(setenv_ret==-1)
    {
        printf("[**ERROR**]: -- setenv reported an error=%d\n",setenv_ret);
    }
    setenv_ret = setenv("OMP_PROC_BIND","spread",1);
    if(setenv_ret==-1)
    {
        printf("[**ERROR**]: -- setenv reported an error=%d\n",setenv_ret);
    }
#pragma omp parallel sections 
{
    #pragma omp section 
    {
        perf_test_simd_ceph_cosf_16xf32_32kiB(&cos_data_in[0],
                                        &ceph_cosf_16xf32_results[0],
                                        &ceph_cosf_16xf32_s[0],
                                        &ceph_cosf_16xf32_e[0],
                                        &ceph_cosf_16xf32_d[0],
                                        &ceph_cosf_16xf32_tsc_aux_s[0],
                                        &ceph_cosf_16xf32_tsc_aux_e[0],
                                        n_runs,
                                        n_samples,
                                        ceph_cosf_16xf32_tid);
        print_thread_affinity();
    }

    #pragma omp section 
    {
        perf_test_mm512_cephes_cos_ps_v3_32kiB(&cos_data_in[0],
                                          &cephes_cos_results[0],
                                          &cephes_cos_ps_s[0],
                                          &cephes_cos_ps_e[0],
                                          &cephes_cos_ps_d[0],
                                          &cephes_cos_tsc_aux_s[0],
                                          &cephes_cos_tsc_aux_e[0],
                                          n_runs,
                                          n_samples,
                                          cephes_cos_tid);
        print_thread_affinity();
    }
}
     
     printf(ANSI_COLOR_GREEN "[PERF-TEST]:  simd_ceph_cosf_16xf32: Started!!\n");
     printf(ANSI_COLOR_GREEN "[PERF-TEST] -- Executed by Core=%d \n",ceph_cosf_16xf32_tid);
     for(int32_t __i{0}; __i != n_runs; ++__i)
     {
         for(int32_t __j{0}; __j != n_samples; ++__j) 
         {  
            unsigned __int64 s{ceph_cosf_16xf32_s[__i*n_samples+__j]};
            unsigned __int64 e{ceph_cosf_16xf32_e[__i*n_samples+__j]};
            unsigned __int64 d{ceph_cosf_16xf32_d[__i*n_samples+__j]};
            std::uint32_t tas{ceph_cosf_16xf32_tsc_aux_s[__i*n_samples+__j]}; // tsc_aux_s
            std::uint32_t tac{ceph_cosf_16xf32_tsc_aux_e[__i*n_samples+__j]}; // tsc_aux_e
            printf(ANSI_COLOR_GREEN "[PMC: RDTSCP] -- Run=%d, start=%llu,end=%llu,delta=%llu,tas=%d,tac=%d\n",__i,s,e,d,tas,tac);
         }
     }
     printf(ANSI_COLOR_GREEN "[PERF-TEST]: simd_ceph_cosf_16xf32: Finished -- dumping-results" ANSI_RESET_ALL"\n\n");

     printf(ANSI_COLOR_WHITE "[PERF-TEST]:  _mm512_cephes_cosf_ps_v3: Started!!\n");
     printf(ANSI_COLOR_WHITE "[PERF-TEST]: -- Executed by Core=%d \n",cephes_cos_tid);
     for(int32_t __i{0}; __i != n_runs; ++__i)
     {
         for(int32_t __j{0}; __j != n_samples; ++__j) 
         {  
            unsigned __int64 s{cephes_cos_ps_s[__i*n_samples+__j]};
            unsigned __int64 e{cephes_cos_ps_e[__i*n_samples+__j]};
            unsigned __int64 d{cephes_cos_ps_d[__i*n_samples+__j]};
            std::uint32_t tas{cephes_cos_tsc_aux_s[__i*n_samples+__j]}; // tsc_aux_s
            std::uint32_t tac{cephes_cos_tsc_aux_e[__i*n_samples+__j]}; // tsc_aux_e
            printf(ANSI_COLOR_WHITE "[PMC: RDTSCP] -- Run=%d, start=%llu,end=%llu,delta=%llu,tas=%d,tac=%d\n",__i,s,e,d,tas,tac);
         }
     }
     printf(ANSI_COLOR_BLUE "[PERF-TEST]: _mm512_cephes_cosf_ps_v3: Finished -- dumping-results" ANSI_RESET_ALL"\n\n");
#if 1
     printf(ANSI_COLOR_MAGENTA "[PERF-TEST]: Result Comparison (using safe floating-point compare as an integer)\n");
     float ceph_cosf_16xf32_result{};
     float cephes_result{};
     bool fcomp_result{};
     for(std::size_t i = 0;i != ROUND_TO_FOUR(sz,3ull); i += 4ull) 
     {
        ceph_cosf_16xf32_result   = ceph_cosf_16xf32_results[i+0ull];
        cephes_result = cephes_cos_results[i+0ull];
        fcomp_result    = almostEqual2sComplement(ceph_cosf_16xf32_result,cephes_result,7);
        if(!fcomp_result) 
        {
            printf(ANSI_COLOR_RED "[PERF-TEST]: ***FAILED***, ceph_cosf_16xf32=%.7f,cephes_cos=%.7f\n",ceph_cosf_16xf32_result,cephes_result);
            printNumber("simd_ceph_cosf_16xf32",ceph_cosf_16xf32_result,0);
            printNumber("_mm512_cephes_cosf_ps_v3",cephes_result,0);
        }
        else
        {
            printf(ANSI_COLOR_CYAN "[PERF-TEST]: ***CORRECT***, ceph_cosf_16xf32=%.7f,cephes_cos=%.7f\n",ceph_cosf_16xf32_result,cephes_result);
            printNumber("simd_ceph_cosf_16xf32",ceph_cosf_16xf32_result,0);
            printNumber("_mm512_cephes_cosf_ps_v3",cephes_result,0);
        }
        ceph_cosf_16xf32_result   = ceph_cosf_16xf32_results[i+1ull];
        cephes_result = cephes_cos_results[i+1ull];
        fcomp_result    = almostEqual2sComplement(ceph_cosf_16xf32_result,cephes_result,7);
        if(!fcomp_result) 
        {
            printf(ANSI_COLOR_RED "[PERF-TEST]: ***FAILED***, ceph_cosf_16xf32=%.7f,cephes_cos=%.7f\n",ceph_cosf_16xf32_result,cephes_result);
            printNumber("simd_ceph_cosf_16xf32",ceph_cosf_16xf32_result,0);
            printNumber("_mm512_cephes_cosf_ps_v3",cephes_result,0);
        }
        else
        {
            printf(ANSI_COLOR_CYAN "[PERF-TEST]: ***CORRECT***, ceph_cosf_16xf32=%.7f,cephes_cos=%.7f\n",ceph_cosf_16xf32_result,cephes_result);
            printNumber("simd_ceph_cosf_16xf32",ceph_cosf_16xf32_result,0);
            printNumber("_mm512_cephes_cosf_ps_v3",cephes_result,0);
        }
        ceph_cosf_16xf32_result   = ceph_cosf_16xf32_results[i+2ull];
        cephes_result = cephes_cos_results[i+2ull];
        fcomp_result    = almostEqual2sComplement(ceph_cosf_16xf32_result,cephes_result,7);
        if(!fcomp_result) 
        {
            printf(ANSI_COLOR_RED "[PERF-TEST]: ***FAILED***, ceph_cosf_16xf32=%.7f,cephes_cos=%.7f\n",ceph_cosf_16xf32_result,cephes_result);
            printNumber("simd_ceph_cosf_16xf32",ceph_cosf_16xf32_result,0);
            printNumber("_mm512_cephes_cosf_ps_v3",cephes_result,0);
        }
        else
        {
            printf(ANSI_COLOR_CYAN "[PERF-TEST]: ***CORRECT***, ceph_cosf_16xf32=%.7f,cephes_cos=%.7f\n",ceph_cosf_16xf32_result,cephes_result);
            printNumber("simd_ceph_cosf_16xf32",ceph_cosf_16xf32_result,0);
            printNumber("_mm512_cephes_cosf_ps_v3",cephes_result,0);
        }
        ceph_cosf_16xf32_result   = ceph_cosf_16xf32_results[i+3ull];
        cephes_result = cephes_cos_results[i+3ull];
        fcomp_result    = almostEqual2sComplement(ceph_cosf_16xf32_result,cephes_result,7);
        if(!fcomp_result) 
        {
            printf(ANSI_COLOR_RED "[PERF-TEST]: ***FAILED***, ceph_cosf_16xf32=%.7f,cephes_cos=%.7f\n",ceph_cosf_16xf32_result,cephes_result);
            printNumber("simd_ceph_cosf_16xf32",ceph_cosf_16xf32_result,0);
            printNumber("_mm512_cephes_cosf_ps_v3",cephes_result,0);
        }
        else
        {
            printf(ANSI_COLOR_CYAN "[PERF-TEST]: ***CORRECT***, ceph_cosf_16xf32=%.7f,cephes_cos=%.7f\n",ceph_cosf_16xf32_result,cephes_result);
            printNumber("simd_ceph_cosf_16xf32",ceph_cosf_16xf32_result,0);
            printNumber("_mm512_cephes_cosf_ps_v3",cephes_result,0);
        }
     }
#endif 
     const bool sz_gt_0 = sz > 0ull;
     if(cos_data_in!=NULL && sz_gt_0)                {gms_mm_free(cos_data_in);}
     if(ceph_cosf_16xf32_results!=NULL && sz_gt_0)   {gms_mm_free(ceph_cosf_16xf32_results);}
     if(cephes_cos_results!=NULL && sz_gt_0)         {gms_mm_free(cephes_cos_results);}
}

/* Currently this test is not working correctly!!*/
__attribute__((hot))
void test_runner_omp_sections_2st_seq();

void test_runner_omp_sections_2st_seq()
{
    using namespace gms::common;
    constexpr int32_t n_runs{10};
    constexpr int32_t n_samples{50};
    constexpr int32_t n_total{n_runs*n_samples};
    constexpr int32_t sz{n_total*16};
    constexpr std::size_t sz_bytes{static_cast<std::size_t>(sz)*sizeof(float)};
    unsigned __int64 ceph_cosf_16xf32_s[n_total]   = {UINT64_MAX};
    unsigned __int64 ceph_cosf_16xf32_e[n_total]   = {UINT64_MAX}; 
    unsigned __int64 ceph_cosf_16xf32_d[n_total]   = {UINT64_MAX};
    unsigned __int64 cephes_cos_ps_s[n_total]      = {UINT64_MAX};
    unsigned __int64 cephes_cos_ps_e[n_total]      = {UINT64_MAX};
    unsigned __int64 cephes_cos_ps_d[n_total]      = {UINT64_MAX};
    std::uint32_t    ceph_cosf_16xf32_tsc_aux_s[n_total]   = {UINT32_MAX};
    std::uint32_t    ceph_cosf_16xf32_tsc_aux_e[n_total]   = {UINT32_MAX};
    std::uint32_t    cephes_cos_tsc_aux_s[n_total] = {UINT32_MAX};
    std::uint32_t    cephes_cos_tsc_aux_e[n_total] = {UINT32_MAX};
    float       * __restrict__ cos_data_in         = NULL;
    float       * __restrict__ cephes_cos_results  = NULL;
    float       * __restrict__ ceph_cosf_16xf32_results   = NULL;
    cos_data_in        = reinterpret_cast<float * __restrict__>(gms_mm_malloc(sz_bytes,64ull));
    cephes_cos_results = reinterpret_cast<float * __restrict__>(gms_mm_malloc(sz_bytes,64ull));
    ceph_cosf_16xf32_results   = reinterpret_cast<float * __restrict__>(gms_mm_malloc(sz_bytes,64ull));
    thread_local std::uniform_real_distribution<float> cos_rand_in;
    thread_local std::uniform_int_distribution<std::int32_t> rand_rad_range;
    thread_local std::uint64_t seed{};
    thread_local std::uint64_t seed_range{};
    float lo;
    float hi;
    seed_range = __rdtsc();
    rand_rad_range = std::uniform_int_distribution<std::int32_t>(0,1);
    auto rand_rad_range_gen{std::mt19937(seed_range)};
    const std::int32_t which_range{rand_rad_range.operator()(rand_rad_range_gen)};
    if(which_range==1)
    {
        lo = -3.14159265358979324f;
        hi = +3.14159265358979324f;
        std::printf("[UNIT_TEST]: Input range: lo=%.7f(rad),hi=%.7f(rad)\n",lo,hi);
    }
    else if(which_range==0)
    {
        lo = -3.14159265358979324f*8192.0f;
        hi = +3.14159265358979324f*8192.0f;
        std::printf("[UNIT_TEST]: Input range: lo=%.7f(rad),hi=%.7f(rad)\n",lo,hi);
    }
    seed = __rdtsc();
    cos_rand_in = std::uniform_real_distribution<float>(lo,hi);
    auto rand_uni_gen{std::mt19937(seed)};
    for(std::size_t i = 0ull; i != ROUND_TO_EIGHT(sz,7); i += 8ull)
    {
         const float in_cos_arg_1{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+0ull] = in_cos_arg_1;
         const float in_cos_arg_2{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+1ull] = in_cos_arg_2;
         const float in_cos_arg_3{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+2ull] = in_cos_arg_3;
         const float in_cos_arg_4{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+3ull] = in_cos_arg_4;
         const float in_cos_arg_5{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+4ull] = in_cos_arg_5;
         const float in_cos_arg_6{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+5ull] = in_cos_arg_6;
         const float in_cos_arg_7{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+6ull] = in_cos_arg_7;
         const float in_cos_arg_8{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+7ull] = in_cos_arg_8;
    }
    
    uint32_t ceph_cosf_16xf32_tid;
    uint32_t cephes_cos_tid;
    int32_t setenv_ret;
    setenv_ret = setenv("OMP_PROC_BIND","true",1);
    if(setenv_ret==-1)
    {
        printf("[**ERROR**]: -- setenv reported an error=%d\n",setenv_ret);
    }
    setenv_ret = setenv("OMP_PROC_BIND","spread",1);
    if(setenv_ret==-1)
    {
        printf("[**ERROR**]: -- setenv reported an error=%d\n",setenv_ret);
    }
#pragma omp parallel sections 
{
    #pragma omp section 
    {
        std::uniform_int_distribution<std::int32_t> rand_kernel_picker;
        std::uint64_t seed{};
        std::int32_t which_kernel_pick{};
        seed = __rdtsc();
        rand_kernel_picker = std::uniform_int_distribution<std::int32_t>(0,1);
        std::mt19937 rand_gen(seed);
        which_kernel_pick = rand_kernel_picker.operator()(rand_gen);
        if(which_kernel_pick==1)
        {
           perf_test_simd_ceph_cosf_16xf32_32kiB(&cos_data_in[0],
                                        &ceph_cosf_16xf32_results[0],
                                        &ceph_cosf_16xf32_s[0],
                                        &ceph_cosf_16xf32_e[0],
                                        &ceph_cosf_16xf32_d[0],
                                        &ceph_cosf_16xf32_tsc_aux_s[0],
                                        &ceph_cosf_16xf32_tsc_aux_e[0],
                                        n_runs,
                                        n_samples,
                                        ceph_cosf_16xf32_tid);
        }
        else if(which_kernel_pick==0)
        {
            perf_test_mm512_cephes_cos_ps_v3_32kiB(&cos_data_in[0],
                                          &cephes_cos_results[0],
                                          &cephes_cos_ps_s[0],
                                          &cephes_cos_ps_e[0],
                                          &cephes_cos_ps_d[0],
                                          &cephes_cos_tsc_aux_s[0],
                                          &cephes_cos_tsc_aux_e[0],
                                          n_runs,
                                          n_samples,
                                          cephes_cos_tid);
        }
        print_thread_affinity();
    }

    #pragma omp section 
    {
        std::uniform_int_distribution<std::int32_t> rand_kernel_picker;
        std::uint64_t seed{};
        std::int32_t which_kernel_pick{};
        seed = __rdtsc();
        rand_kernel_picker = std::uniform_int_distribution<std::int32_t>(0,1);
        std::mt19937 rand_gen(seed);
        which_kernel_pick = rand_kernel_picker.operator()(rand_gen);
        if(which_kernel_pick==1)
        {
            perf_test_mm512_cephes_cos_ps_v3_32kiB(&cos_data_in[0],
                                          &cephes_cos_results[0],
                                          &cephes_cos_ps_s[0],
                                          &cephes_cos_ps_e[0],
                                          &cephes_cos_ps_d[0],
                                          &cephes_cos_tsc_aux_s[0],
                                          &cephes_cos_tsc_aux_e[0],
                                          n_runs,
                                          n_samples,
                                          cephes_cos_tid);
        }
        else if(which_kernel_pick==0)
        {
            perf_test_simd_ceph_cosf_16xf32_32kiB(&cos_data_in[0],
                                        &ceph_cosf_16xf32_results[0],
                                        &ceph_cosf_16xf32_s[0],
                                        &ceph_cosf_16xf32_e[0],
                                        &ceph_cosf_16xf32_d[0],
                                        &ceph_cosf_16xf32_tsc_aux_s[0],
                                        &ceph_cosf_16xf32_tsc_aux_e[0],
                                        n_runs,
                                        n_samples,
                                        ceph_cosf_16xf32_tid);
        }
        print_thread_affinity();
    }
}
     
     printf(ANSI_COLOR_GREEN "[PERF-TEST]:  simd_ceph_cosf_16xf32: Started!!\n");
     printf(ANSI_COLOR_GREEN "[PERF-TEST] -- Executed by Core=%d \n",ceph_cosf_16xf32_tid);
     for(int32_t __i{0}; __i != n_runs; ++__i)
     {
         for(int32_t __j{0}; __j != n_samples; ++__j) 
         {  
            unsigned __int64 s{ceph_cosf_16xf32_s[__i*n_samples+__j]};
            unsigned __int64 e{ceph_cosf_16xf32_e[__i*n_samples+__j]};
            unsigned __int64 d{ceph_cosf_16xf32_d[__i*n_samples+__j]};
            std::uint32_t tas{ceph_cosf_16xf32_tsc_aux_s[__i*n_samples+__j]}; // tsc_aux_s
            std::uint32_t tac{ceph_cosf_16xf32_tsc_aux_e[__i*n_samples+__j]}; // tsc_aux_e
            printf(ANSI_COLOR_GREEN "[PMC: RDTSCP] -- Run=%d, start=%llu,end=%llu,delta=%llu,tas=%d,tac=%d\n",__i,s,e,d,tas,tac);
         }
     }
     printf(ANSI_COLOR_GREEN "[PERF-TEST]: simd_ceph_cosf_16xf32: Finished -- dumping-results" ANSI_RESET_ALL"\n\n");

     printf(ANSI_COLOR_WHITE "[PERF-TEST]:  _mm512_cephes_cosf_ps_v3: Started!!\n");
     printf(ANSI_COLOR_WHITE "[PERF-TEST]: -- Executed by Core=%d \n",cephes_cos_tid);
     for(int32_t __i{0}; __i != n_runs; ++__i)
     {
         for(int32_t __j{0}; __j != n_samples; ++__j) 
         {  
            unsigned __int64 s{cephes_cos_ps_s[__i*n_samples+__j]};
            unsigned __int64 e{cephes_cos_ps_e[__i*n_samples+__j]};
            unsigned __int64 d{cephes_cos_ps_d[__i*n_samples+__j]};
            std::uint32_t tas{cephes_cos_tsc_aux_s[__i*n_samples+__j]}; // tsc_aux_s
            std::uint32_t tac{cephes_cos_tsc_aux_e[__i*n_samples+__j]}; // tsc_aux_e
            printf(ANSI_COLOR_WHITE "[PMC: RDTSCP] -- Run=%d, start=%llu,end=%llu,delta=%llu,tas=%d,tac=%d\n",__i,s,e,d,tas,tac);
         }
     }
     printf(ANSI_COLOR_BLUE "[PERF-TEST]: _mm512_cephes_cosf_ps_v3: Finished -- dumping-results" ANSI_RESET_ALL"\n\n");
#if 1
     printf(ANSI_COLOR_MAGENTA "[PERF-TEST]: Result Comparison (using safe floating-point compare as an integer)\n");
     float ceph_cosf_16xf32_result{};
     float cephes_result{};
     bool fcomp_result{};
     for(std::size_t i = 0;i != ROUND_TO_FOUR(sz,3ull); i += 4ull) 
     {
        ceph_cosf_16xf32_result   = ceph_cosf_16xf32_results[i+0ull];
        cephes_result = cephes_cos_results[i+0ull];
        fcomp_result    = almostEqual2sComplement(ceph_cosf_16xf32_result,cephes_result,7);
        if(!fcomp_result) 
        {
            printf(ANSI_COLOR_RED "[PERF-TEST]: ***FAILED***, ceph_cosf_16xf32=%.7f,cephes_cos=%.7f\n",ceph_cosf_16xf32_result,cephes_result);
            printNumber("simd_ceph_cosf_16xf32",ceph_cosf_16xf32_result,0);
            printNumber("_mm512_cephes_cosf_ps_v3",cephes_result,0);
        }
        else
        {
            printf(ANSI_COLOR_CYAN "[PERF-TEST]: ***CORRECT***, ceph_cosf_16xf32=%.7f,cephes_cos=%.7f\n",ceph_cosf_16xf32_result,cephes_result);
            printNumber("simd_ceph_cosf_16xf32",ceph_cosf_16xf32_result,0);
            printNumber("_mm512_cephes_cosf_ps_v3",cephes_result,0);
        }
        ceph_cosf_16xf32_result   = ceph_cosf_16xf32_results[i+1ull];
        cephes_result = cephes_cos_results[i+1ull];
        fcomp_result    = almostEqual2sComplement(ceph_cosf_16xf32_result,cephes_result,7);
        if(!fcomp_result) 
        {
            printf(ANSI_COLOR_RED "[PERF-TEST]: ***FAILED***, ceph_cosf_16xf32=%.7f,cephes_cos=%.7f\n",ceph_cosf_16xf32_result,cephes_result);
            printNumber("simd_ceph_cosf_16xf32",ceph_cosf_16xf32_result,0);
            printNumber("_mm512_cephes_cosf_ps_v3",cephes_result,0);
        }
        else
        {
            printf(ANSI_COLOR_CYAN "[PERF-TEST]: ***CORRECT***, ceph_cosf_16xf32=%.7f,cephes_cos=%.7f\n",ceph_cosf_16xf32_result,cephes_result);
            printNumber("simd_ceph_cosf_16xf32",ceph_cosf_16xf32_result,0);
            printNumber("_mm512_cephes_cosf_ps_v3",cephes_result,0);
        }
        ceph_cosf_16xf32_result   = ceph_cosf_16xf32_results[i+2ull];
        cephes_result = cephes_cos_results[i+2ull];
        fcomp_result    = almostEqual2sComplement(ceph_cosf_16xf32_result,cephes_result,7);
        if(!fcomp_result) 
        {
            printf(ANSI_COLOR_RED "[PERF-TEST]: ***FAILED***, ceph_cosf_16xf32=%.7f,cephes_cos=%.7f\n",ceph_cosf_16xf32_result,cephes_result);
            printNumber("simd_ceph_cosf_16xf32",ceph_cosf_16xf32_result,0);
            printNumber("_mm512_cephes_cosf_ps_v3",cephes_result,0);
        }
        else
        {
            printf(ANSI_COLOR_CYAN "[PERF-TEST]: ***CORRECT***, ceph_cosf_16xf32=%.7f,cephes_cos=%.7f\n",ceph_cosf_16xf32_result,cephes_result);
            printNumber("simd_ceph_cosf_16xf32",ceph_cosf_16xf32_result,0);
            printNumber("_mm512_cephes_cosf_ps_v3",cephes_result,0);
        }
        ceph_cosf_16xf32_result   = ceph_cosf_16xf32_results[i+3ull];
        cephes_result = cephes_cos_results[i+3ull];
        fcomp_result    = almostEqual2sComplement(ceph_cosf_16xf32_result,cephes_result,7);
        if(!fcomp_result) 
        {
            printf(ANSI_COLOR_RED "[PERF-TEST]: ***FAILED***, ceph_cosf_16xf32=%.7f,cephes_cos=%.7f\n",ceph_cosf_16xf32_result,cephes_result);
            printNumber("simd_ceph_cosf_16xf32",ceph_cosf_16xf32_result,0);
            printNumber("_mm512_cephes_cosf_ps_v3",cephes_result,0);
        }
        else
        {
            printf(ANSI_COLOR_CYAN "[PERF-TEST]: ***CORRECT***, ceph_cosf_16xf32=%.7f,cephes_cos=%.7f\n",ceph_cosf_16xf32_result,cephes_result);
            printNumber("simd_ceph_cosf_16xf32",ceph_cosf_16xf32_result,0);
            printNumber("_mm512_cephes_cosf_ps_v3",cephes_result,0);
        }
     }
#endif 
     const bool sz_gt_0 = sz > 0ull;
     if(cos_data_in!=NULL && sz_gt_0)        {gms_mm_free(cos_data_in);}
     if(ceph_cosf_16xf32_results!=NULL && sz_gt_0)   {gms_mm_free(ceph_cosf_16xf32_results);}
     if(cephes_cos_results!=NULL && sz_gt_0) {gms_mm_free(cephes_cos_results);}
}

__attribute__((hot))
void test_runner_simd_ceph_cosf_16xf32_single_thread(const int32_t, const int32_t);

void test_runner_simd_ceph_cosf_16xf32_single_thread(const int32_t cpu_core, const int32_t priority)
{
    using namespace gms::common;
    constexpr int32_t n_runs{10};
    constexpr int32_t n_samples{50};
    constexpr int32_t n_total{n_runs*n_samples};
    constexpr int32_t sz{n_total*16};
    constexpr std::size_t sz_bytes{static_cast<std::size_t>(sz)*sizeof(float)};
    unsigned __int64 ceph_cosf_16xf32_s[n_total]        = {UINT64_MAX};
    unsigned __int64 ceph_cosf_16xf32_e[n_total]        = {UINT64_MAX}; 
    unsigned __int64 ceph_cosf_16xf32_d[n_total]        = {UINT64_MAX};
    std::uint32_t    ceph_cosf_16xf32_tsc_aux_s[n_total]   = {UINT32_MAX};
    std::uint32_t    ceph_cosf_16xf32_tsc_aux_e[n_total]   = {UINT32_MAX};
    float       * __restrict__ cos_data_in        = NULL;
    float       * __restrict__ ceph_cosf_16xf32_results   = NULL;
    cos_data_in        = reinterpret_cast<float * __restrict__>(gms_mm_malloc(sz_bytes,64ull));
    ceph_cosf_16xf32_results   = reinterpret_cast<float * __restrict__>(gms_mm_malloc(sz_bytes,64ull));
    thread_local std::uniform_real_distribution<float> cos_rand_in;
    thread_local std::uniform_int_distribution<std::int32_t> rand_rad_range;
    thread_local std::uint64_t seed{};
    thread_local std::uint64_t seed_range{};
    float lo;
    float hi;
    seed_range = __rdtsc();
    rand_rad_range = std::uniform_int_distribution<std::int32_t>(0,1);
    auto rand_rad_range_gen{std::mt19937(seed_range)};
    const std::int32_t which_range{rand_rad_range.operator()(rand_rad_range_gen)};
    if(which_range==1)
    {
        lo = -3.14159265358979324f;
        hi = +3.14159265358979324f;
        std::printf("[UNIT_TEST]: Input range: lo=%.7f(rad),hi=%.7f(rad)\n",lo,hi);
    }
    else if(which_range==0)
    {
        lo = -3.14159265358979324f*8192.0f;
        hi = +3.14159265358979324f*8192.0f;
        std::printf("[UNIT_TEST]: Input range: lo=%.7f(rad),hi=%.7f(rad)\n",lo,hi);
    }
    seed = __rdtsc();
    cos_rand_in = std::uniform_real_distribution<float>(lo,hi);
    auto rand_uni_gen{std::mt19937(seed)};
    for(std::size_t i = 0ull; i != ROUND_TO_EIGHT(sz,7); i += 8ull)
    {
         const float in_cos_arg_1{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+0ull] = in_cos_arg_1;
         const float in_cos_arg_2{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+1ull] = in_cos_arg_2;
         const float in_cos_arg_3{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+2ull] = in_cos_arg_3;
         const float in_cos_arg_4{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+3ull] = in_cos_arg_4;
         const float in_cos_arg_5{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+4ull] = in_cos_arg_5;
         const float in_cos_arg_6{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+5ull] = in_cos_arg_6;
         const float in_cos_arg_7{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+6ull] = in_cos_arg_7;
         const float in_cos_arg_8{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+7ull] = in_cos_arg_8;
    }
    
    uint32_t ceph_cosf_16xf32_tid;
    int32_t setenv_ret;
    int32_t status;
    setenv_ret = setenv("OMP_PROC_BIND","true",1);
    if(setenv_ret==-1)
    {
        printf("[**ERROR**]: -- setenv reported an error=%d\n",setenv_ret);
    }
    setenv_ret = setenv("OMP_PROC_BIND","spread",1);
    if(setenv_ret==-1)
    {
        printf("[**ERROR**]: -- setenv reported an error=%d\n",setenv_ret);
    }
    status = set_affinity_and_priority(cpu_core,priority);
    if(status > 0)
    { 
       int32_t which = PRIO_PROCESS;
       id_t pid      = getpid();
       errno         = 0;
       int32_t nice  = getpriority(which,pid);
       printf("[***ERROR***]: set_affinity_and_priority status:%d,errno=%d,pid=%d,nice=%d\n",status,errno,pid,nice);
    }

    perf_test_simd_ceph_cosf_16xf32_32kiB(&cos_data_in[0],
                                    &ceph_cosf_16xf32_results[0],
                                    &ceph_cosf_16xf32_s[0],
                                    &ceph_cosf_16xf32_e[0],
                                    &ceph_cosf_16xf32_d[0],
                                    &ceph_cosf_16xf32_tsc_aux_s[0],
                                    &ceph_cosf_16xf32_tsc_aux_e[0],
                                    n_runs,
                                    n_samples,
                                    ceph_cosf_16xf32_tid);
    print_thread_affinity();
     
    printf(ANSI_COLOR_GREEN "[PERF-TEST]:  simd_ceph_cosf_16xf32: Started!!\n");
    printf(ANSI_COLOR_GREEN "[PERF-TEST] -- Executed by Core=%d \n",ceph_cosf_16xf32_tid);
    for(int32_t __i{0}; __i != n_runs; ++__i)
    {
        for(int32_t __j{0}; __j != n_samples; ++__j) 
        {  
            unsigned __int64 s{ceph_cosf_16xf32_s[__i*n_samples+__j]};
            unsigned __int64 e{ceph_cosf_16xf32_e[__i*n_samples+__j]};
            unsigned __int64 d{ceph_cosf_16xf32_d[__i*n_samples+__j]};
            std::uint32_t tas{ceph_cosf_16xf32_tsc_aux_s[__i*n_samples+__j]}; // tsc_aux_s
            std::uint32_t tac{ceph_cosf_16xf32_tsc_aux_e[__i*n_samples+__j]}; // tsc_aux_e
            printf(ANSI_COLOR_GREEN "[PMC: RDTSCP] -- Run=%d, start=%llu,end=%llu,delta=%llu,tas=%d,tac=%d\n",__i,s,e,d,tas,tac);
        }
    }
     printf(ANSI_COLOR_GREEN "[PERF-TEST]: simd_ceph_cosf_16xf32: Finished -- dumping-results" ANSI_RESET_ALL"\n\n");

#if 1
     float ceph_cosf_16xf32_result{};
     for(std::size_t i = 0;i != ROUND_TO_FOUR(sz,3ull); i += 4ull) 
     {
        ceph_cosf_16xf32_result   = ceph_cosf_16xf32_results[i+0ull];
        printNumber("simd_ceph_cosf_16xf32",ceph_cosf_16xf32_result,0);
        ceph_cosf_16xf32_result   = ceph_cosf_16xf32_results[i+1ull];       
        printNumber("simd_ceph_cosf_16xf32",ceph_cosf_16xf32_result,0);           
        ceph_cosf_16xf32_result   = ceph_cosf_16xf32_results[i+2ull];     
        printNumber("simd_ceph_cosf_16xf32",ceph_cosf_16xf32_result,0);      
        ceph_cosf_16xf32_result   = ceph_cosf_16xf32_results[i+3ull];       
        printNumber("simd_ceph_cosf_16xf32",ceph_cosf_16xf32_result,0);
                  
     }
#endif 
     const bool sz_gt_0 = sz > 0ull;
     if(cos_data_in!=NULL && sz_gt_0)        {gms_mm_free(cos_data_in);}
     if(ceph_cosf_16xf32_results!=NULL && sz_gt_0)   {gms_mm_free(ceph_cosf_16xf32_results);}
}

__attribute__((hot))
void test_runner_mm512_cephes_cosf_ps_v3_single_thread(const int32_t,const int32_t);

void test_runner_mm512_cephes_cosf_ps_v3_single_thread(const int32_t cpu_core,const int32_t priority)
{
    using namespace gms::common;
    constexpr int32_t n_runs{10};
    constexpr int32_t n_samples{50};
    constexpr int32_t n_total{n_runs*n_samples};
    constexpr int32_t sz{n_total*16};
    constexpr std::size_t sz_bytes{static_cast<std::size_t>(sz)*sizeof(float)};
    unsigned __int64 cephes_cos_ps_s[n_total]      = {UINT64_MAX};
    unsigned __int64 cephes_cos_ps_e[n_total]      = {UINT64_MAX};
    unsigned __int64 cephes_cos_ps_d[n_total]      = {UINT64_MAX};
    std::uint32_t    cephes_cos_tsc_aux_s[n_total] = {UINT32_MAX};
    std::uint32_t    cephes_cos_tsc_aux_e[n_total] = {UINT32_MAX};
    float       * __restrict__ cos_data_in        = NULL;
    float       * __restrict__ cephes_cos_results = NULL;
    cos_data_in        = reinterpret_cast<float * __restrict__>(gms_mm_malloc(sz_bytes,64ull));
    cephes_cos_results = reinterpret_cast<float * __restrict__>(gms_mm_malloc(sz_bytes,64ull));
    thread_local std::uniform_real_distribution<float> cos_rand_in;
    thread_local std::uniform_int_distribution<std::int32_t> rand_rad_range;
    thread_local std::uint64_t seed{};
    thread_local std::uint64_t seed_range{};
    float lo;
    float hi;
    seed_range = __rdtsc();
    rand_rad_range = std::uniform_int_distribution<std::int32_t>(0,1);
    auto rand_rad_range_gen{std::mt19937(seed_range)};
    const std::int32_t which_range{rand_rad_range.operator()(rand_rad_range_gen)};
    if(which_range==1)
    {
        lo = -3.14159265358979324f;
        hi = +3.14159265358979324f;
        std::printf("[UNIT_TEST]: Input range: lo=%.7f(rad),hi=%.7f(rad)\n",lo,hi);
    }
    else if(which_range==0)
    {
        lo = -3.14159265358979324f*8192.0f;
        hi = +3.14159265358979324f*8192.0f;
        std::printf("[UNIT_TEST]: Input range: lo=%.7f(rad),hi=%.7f(rad)\n",lo,hi);
    }
    seed = __rdtsc();
    cos_rand_in = std::uniform_real_distribution<float>(lo,hi);
    auto rand_uni_gen{std::mt19937(seed)};
    for(std::size_t i = 0ull; i != ROUND_TO_EIGHT(sz,7); i += 8ull)
    {
         const float in_cos_arg_1{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+0ull] = in_cos_arg_1;
         const float in_cos_arg_2{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+1ull] = in_cos_arg_2;
         const float in_cos_arg_3{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+2ull] = in_cos_arg_3;
         const float in_cos_arg_4{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+3ull] = in_cos_arg_4;
         const float in_cos_arg_5{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+4ull] = in_cos_arg_5;
         const float in_cos_arg_6{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+5ull] = in_cos_arg_6;
         const float in_cos_arg_7{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+6ull] = in_cos_arg_7;
         const float in_cos_arg_8{cos_rand_in.operator()(rand_uni_gen)};
         cos_data_in[i+7ull] = in_cos_arg_8;
    }
    
    uint32_t cephes_cos_tid;
    int32_t setenv_ret;
    int32_t status{};
    setenv_ret = setenv("OMP_PROC_BIND","true",1);
    if(setenv_ret==-1)
    {
        printf("[**ERROR**]: -- setenv reported an error=%d\n",setenv_ret);
    }
    setenv_ret = setenv("OMP_PROC_BIND","spread",1);
    if(setenv_ret==-1)
    {
        printf("[**ERROR**]: -- setenv reported an error=%d\n",setenv_ret);
    }
    status = set_affinity_and_priority(cpu_core,priority);
    if(status > 0)
    { 
       int32_t which = PRIO_PROCESS;
       id_t pid      = getpid();
       errno         = 0;
       int32_t nice  = getpriority(which,pid);
       printf("[***ERROR***]: set_affinity_and_priority status:%d,errno=%d,pid=%d,nice=%d\n",status,errno,pid,nice);
    }
    perf_test_mm512_cephes_cos_ps_v3_32kiB(&cos_data_in[0],
                                          &cephes_cos_results[0],
                                          &cephes_cos_ps_s[0],
                                          &cephes_cos_ps_e[0],
                                          &cephes_cos_ps_d[0],
                                          &cephes_cos_tsc_aux_s[0],
                                          &cephes_cos_tsc_aux_e[0],
                                          n_runs,
                                          n_samples,
                                          cephes_cos_tid);    
    print_thread_affinity();
    printf(ANSI_COLOR_WHITE "[PERF-TEST]:  _mm512_cephes_cosf_ps_v3: Started!!\n");
    printf(ANSI_COLOR_WHITE "[PERF-TEST]: -- Executed by Core=%d \n",cephes_cos_tid);
    for(int32_t __i{0}; __i != n_runs; ++__i)
    {
         for(int32_t __j{0}; __j != n_samples; ++__j) 
         {  
            unsigned __int64 s{cephes_cos_ps_s[__i*n_samples+__j]};
            unsigned __int64 e{cephes_cos_ps_e[__i*n_samples+__j]};
            unsigned __int64 d{cephes_cos_ps_d[__i*n_samples+__j]};
            std::uint32_t tas{cephes_cos_tsc_aux_s[__i*n_samples+__j]}; // tsc_aux_s
            std::uint32_t tac{cephes_cos_tsc_aux_e[__i*n_samples+__j]}; // tsc_aux_e
            printf(ANSI_COLOR_WHITE "[PMC: RDTSCP] -- Run=%d, start=%llu,end=%llu,delta=%llu,tas=%d,tac=%d\n",__i,s,e,d,tas,tac);
         }
     }
     printf(ANSI_COLOR_BLUE "[PERF-TEST]: _mm512_cephes_cosf_ps_v3: Finished -- dumping-results" ANSI_RESET_ALL"\n\n");
     float cephes_result{};
     for(std::size_t i = 0;i != ROUND_TO_FOUR(sz,3ull); i += 4ull) 
     {
      
        cephes_result = cephes_cos_results[i+0ull];       
        printNumber("_mm512_cephes_cosf_ps_v3",cephes_result,0);       
        cephes_result = cephes_cos_results[i+1ull];         
        printNumber("_mm512_cephes_cosf_ps_v3",cephes_result,0);       
        cephes_result = cephes_cos_results[i+2ull];       
        printNumber("_mm512_cephes_cosf_ps_v3",cephes_result,0);       
        cephes_result = cephes_cos_results[i+3ull];
        printNumber("_mm512_cephes_cosf_ps_v3",cephes_result,0);
       
     }
     const bool sz_gt_0 = sz > 0ull;
     if(cos_data_in!=NULL && sz_gt_0)        {gms_mm_free(cos_data_in);}
     if(cephes_cos_results!=NULL && sz_gt_0) {gms_mm_free(cephes_cos_results);}
}

__attribute__((hot))
void test_runner_mm512_cephes_cosf_ps_v3_parallel();

void test_runner_mm512_cephes_cosf_ps_v3_parallel()
{
    thread_local std::uniform_int_distribution<std::int32_t> rand_thread_prio;
    thread_local std::uint64_t prio_seed{};
    thread_local std::int32_t priority{};
    thread_local std::mt19937 rand_prio_gen;
    int32_t tid{};
#pragma omp parallel private(tid) num_threads(omp_get_max_threads())
{
    rand_thread_prio = std::uniform_int_distribution<std::int32_t>(0,99);
    prio_seed        = __rdtsc();
    rand_prio_gen    = std::mt19937(prio_seed);
    priority         = rand_thread_prio.operator()(rand_prio_gen);
    tid              = omp_get_thread_num();
    test_runner_mm512_cephes_cosf_ps_v3_single_thread(tid,priority);
}

}

__attribute__((hot))
void test_runner_simd_ceph_cosf_16xf32_parallel();

void test_runner_simd_ceph_cosf_16xf32_parallel()
{
    thread_local std::uniform_int_distribution<std::int32_t> rand_thread_prio;
    thread_local std::uint64_t prio_seed{};
    thread_local std::int32_t priority{};
    thread_local std::mt19937 rand_prio_gen;
    int32_t tid{};
#pragma omp parallel private(tid) num_threads(omp_get_max_threads())
{
    rand_thread_prio = std::uniform_int_distribution<std::int32_t>(0,99);
    prio_seed        = __rdtsc();
    rand_prio_gen    = std::mt19937(prio_seed);
    priority         = rand_thread_prio.operator()(rand_prio_gen);
    tid              = omp_get_thread_num();
    test_runner_simd_ceph_cosf_16xf32_single_thread(tid,priority);
}

}

int main()
{
    test_runner_omp_sections_1st_seq();
    //test_runner_omp_sections_2st_seq();
    test_runner_simd_ceph_cosf_16xf32_single_thread(4,99);
    test_runner_mm512_cephes_cosf_ps_v3_single_thread(5,99);
    test_runner_mm512_cephes_cosf_ps_v3_parallel();
    test_runner_simd_ceph_cosf_16xf32_parallel();
    return 0;
}
