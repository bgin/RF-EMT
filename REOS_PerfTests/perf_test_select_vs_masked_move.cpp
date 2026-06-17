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
#include "GMS_ansi_colors.h"
#include "GMS_simd_cephes_cos.h"
#include "GMS_vec_cephes_cos_ps.h"


/*
   icpc -o perf_test_select_vs_masked_move_64B -O3 -fp-model fast=2 -fno-exceptions -fopenmp  -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_ansi_colors.h GMS_vec_cephes_cos_ps GMS_simd_cephes_cos.h perf_test_select_vs_masked_move.cpp
   ASM: 
   icpc -S -O3 -fverbose-asm -masm=intel  -fno-exceptions -fopenmp  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_ansi_colors.h GMS_vec_cephes_cos_ps GMS_simd_cephes_cos.h perf_test_select_vs_masked_move.cpp

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

/*Masked-move based*/
__attribute__((hot))
__attribute__((noinline))
void perf_test_mm512_cephes_cos_ps_v3_64B(__m512,
                                       __m512 * __restrict__,
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
void perf_test_mm512_cephes_cos_ps_v3_64B(__m512 cos_args,
                                        __m512 * __restrict__ cos_values,
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
    [[maybe_unused]]
    volatile unsigned __int64 warmup_start;
    [[maybe_unused]]
    volatile unsigned __int64 warmup_end;
    [[maybe_unused]]
    volatile __m512 warmup_cos_values;
    register __m512 result_cos_values;
    uint32_t mem_start{0};
    uint32_t mem_end{0}; 
    // warmup of RDTSCP
    warmup_start = __rdtscp(&mem_start);
    warmup_end   = __rdtscp(&mem_end);
    // warmup call of _mm512_cephes_cosf_ps_v2
    __m512 a_cos_args = cos_args;
    warmup_cos_values = gms::math::_mm512_ceph_cosf_ps_v3<true>(a_cos_args);

    // Main testing loop begins.
    for(int32_t __i{0}; __i != n_runs; ++__i)
    {
        for(int32_t __j{0}; __j != n_samples; ++__j) 
        {   
            __asm__ __volatile__ ("lfence");
            unsigned __int64 start_curr{__rdtscp(&mem_start)};       
            result_cos_values = gms::math::_mm512_ceph_cosf_ps_v3<false>(a_cos_args);
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
        }
    }
    tid = omp_get_thread_num();
    *cos_values = result_cos_values;
}


__attribute__((hot))
__attribute__((noinline))
void perf_test_simd_ceph_cosf_16xf32_64B(__m512,
                                    __m512 * __restrict__,
                                    unsigned __int64 * __restrict__ ,
                                    unsigned __int64 * __restrict__ ,
                                    unsigned __int64 * __restrict__ ,
                                    std::uint32_t    * __restrict__ ,
                                    std::uint32_t    * __restrict__ ,
                                    const int32_t ,
                                    const int32_t,
                                    uint32_t &);
/* IF-Conversion select-based (merge)*/
__attribute__((hot))
__attribute__((noinline))
void perf_test_simd_ceph_cosf_16xf32_64B(__m512 cos_args,
                                    __m512 * __restrict__ cos_values,
                                    unsigned __int64 * __restrict__ ceph_v2_cos_s,
                                    unsigned __int64 * __restrict__ ceph_v2_cos_e,
                                    unsigned __int64 * __restrict__ ceph_v2_cos_d,
                                    std::uint32_t    * __restrict__ ceph_v2_cos_tsc_aux_s,
                                    std::uint32_t    * __restrict__ ceph_v2_cos_tsc_aux_e,
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
    [[maybe_unused]]
    volatile unsigned __int64 warmup_start;
    [[maybe_unused]]
    volatile unsigned __int64 warmup_end;
    [[maybe_unused]]
    volatile __m512 warmup_cos_values;
    register __m512 result_cos_values;
    uint32_t mem_start{0};
    uint32_t mem_end{0}; 
    // warmup of RDTSCP
    warmup_start = __rdtscp(&mem_start);
    warmup_end   = __rdtscp(&mem_end);
    // warmup call of _mm512_cos_ps
    const __m512 a_cos_args = cos_args;
    warmup_cos_values = gms::math::simd_ceph_cosf_16xf32(a_cos_args);

    // Main testing loop begins.
    for(int32_t __i{0}; __i != n_runs; ++__i)
    {
        for(int32_t __j{0}; __j != n_samples; ++__j) 
        {   
            __asm__ __volatile__ ("lfence");
            unsigned __int64 start_curr{__rdtscp(&mem_start)};
            result_cos_values = gms::math::simd_ceph_cosf_16xf32(a_cos_args);
            unsigned __int64 end_curr{__rdtscp(&mem_end)};
            __asm__ __volatile__ ("lfence");
            //remove latency
            unsigned __int64 start_corrected{start_curr-RDTSCP_LAT};
            ceph_v2_cos_s[__i*n_samples+__j] = start_corrected;
            unsigned __int64 end_corrected{end_curr-RDTSCP_LAT};
            ceph_v2_cos_e[__i*n_samples+__j]   = end_corrected;
            ceph_v2_cos_d[__i*n_samples+__j] = end_corrected-start_corrected;
            ceph_v2_cos_tsc_aux_s[__i*n_samples+__j] = mem_start;
            ceph_v2_cos_tsc_aux_e[__i*n_samples+__j] = mem_end;
        }
    }
    tid = omp_get_thread_num();
    *cos_values = result_cos_values;
}


__attribute__((hot))
void test_runner_omp_sections_1st_seq();

void test_runner_omp_sections_1st_seq()
{
    
    constexpr int32_t n_runs{10};
    constexpr int32_t n_samples{50};
    constexpr int32_t n_total{n_runs*n_samples};
    unsigned __int64 ceph_v2_cos_ps_s[n_total]     = {UINT64_MAX};
    unsigned __int64 ceph_v2_cos_ps_e[n_total]     = {UINT64_MAX}; 
    unsigned __int64 ceph_v2_cos_ps_d[n_total]     = {UINT64_MAX};
    unsigned __int64 cephes_cos_ps_s[n_total]      = {UINT64_MAX};
    unsigned __int64 cephes_cos_ps_e[n_total]      = {UINT64_MAX};
    unsigned __int64 cephes_cos_ps_d[n_total]      = {UINT64_MAX};
    std::uint32_t    ceph_v2_cos_tsc_aux_s[n_total]= {UINT32_MAX};
    std::uint32_t    ceph_v2_cos_tsc_aux_e[n_total]= {UINT32_MAX};
    std::uint32_t    cephes_cos_tsc_aux_s[n_total] = {UINT32_MAX};
    std::uint32_t    cephes_cos_tsc_aux_e[n_total] = {UINT32_MAX};
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
    const float in_cos_arg_1{cos_rand_in.operator()(rand_uni_gen)};
    const float in_cos_arg_2{cos_rand_in.operator()(rand_uni_gen)};
    const float in_cos_arg_3{cos_rand_in.operator()(rand_uni_gen)};
    const float in_cos_arg_4{cos_rand_in.operator()(rand_uni_gen)};
    const float in_cos_arg_5{cos_rand_in.operator()(rand_uni_gen)};
    const float in_cos_arg_6{cos_rand_in.operator()(rand_uni_gen)};
    const float in_cos_arg_7{cos_rand_in.operator()(rand_uni_gen)};
    const float in_cos_arg_8{cos_rand_in.operator()(rand_uni_gen)};
    const float in_cos_arg_9{cos_rand_in.operator()(rand_uni_gen)};
    const float in_cos_arg_10{cos_rand_in.operator()(rand_uni_gen)};
    const float in_cos_arg_11{cos_rand_in.operator()(rand_uni_gen)};
    const float in_cos_arg_12{cos_rand_in.operator()(rand_uni_gen)};
    const float in_cos_arg_13{cos_rand_in.operator()(rand_uni_gen)};
    const float in_cos_arg_14{cos_rand_in.operator()(rand_uni_gen)};
    const float in_cos_arg_15{cos_rand_in.operator()(rand_uni_gen)};
    const float in_cos_arg_16{cos_rand_in.operator()(rand_uni_gen)};
    __m512  in_vec_cos_arg{_mm512_setr_ps( in_cos_arg_1,in_cos_arg_2,in_cos_arg_3,in_cos_arg_4,
                                           in_cos_arg_5,in_cos_arg_6,in_cos_arg_7,in_cos_arg_8,
                                           in_cos_arg_9,in_cos_arg_10,in_cos_arg_11,in_cos_arg_12,
                                           in_cos_arg_13,in_cos_arg_14,in_cos_arg_15,in_cos_arg_16)};
    __m512  ceph_v2_cos_results;
    __m512  cephes_cos_results;
    uint32_t ceph_v2_cos_tid;
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
        perf_test_simd_ceph_cosf_16xf32_64B(in_vec_cos_arg,
                                        &ceph_v2_cos_results,
                                        &ceph_v2_cos_ps_s[0],
                                        &ceph_v2_cos_ps_e[0],
                                        &ceph_v2_cos_ps_d[0],
                                        &ceph_v2_cos_tsc_aux_s[0],
                                        &ceph_v2_cos_tsc_aux_e[0],
                                        n_runs,
                                        n_samples,
                                        ceph_v2_cos_tid);
        print_thread_affinity();
    }

    #pragma omp section 
    {
        perf_test_mm512_cephes_cos_ps_v3_64B(in_vec_cos_arg,
                                          &cephes_cos_results,
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
     const float * __restrict__ ptr_in_vec_cos_arg = reinterpret_cast<const float * __restrict__>(&in_vec_cos_arg);
     printf(ANSI_COLOR_GREEN "[PERF-TEST]:  simd_ceph_cosf_16xf32: Started!!\n");
     printf(ANSI_COLOR_GREEN "[PERF-TEST] -- Executed by Core=%d \n",ceph_v2_cos_tid);
     for(std::int32_t i = 0;i != 15; ++i) { std::printf("[PERF-TEST]: Index=%d: input cosine arguments=%.7f(rad)\n",i,ptr_in_vec_cos_arg[i]);}
     for(int32_t __i{0}; __i != n_runs; ++__i)
     {
         for(int32_t __j{0}; __j != n_samples; ++__j) 
         {  
            unsigned __int64 s{ceph_v2_cos_ps_s[__i*n_samples+__j]};
            unsigned __int64 e{ceph_v2_cos_ps_e[__i*n_samples+__j]};
            unsigned __int64 d{ceph_v2_cos_ps_d[__i*n_samples+__j]};
            std::uint32_t tas{ceph_v2_cos_tsc_aux_s[__i*n_samples+__j]}; // tsc_aux_s
            std::uint32_t tac{ceph_v2_cos_tsc_aux_e[__i*n_samples+__j]}; // tsc_aux_e
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
     printf(ANSI_COLOR_MAGENTA "[PERF-TEST]: Result Comparison (using safe floating-point compare as an integer)\n");
     const float * __restrict__ p_ceph_v2_cos_results   = reinterpret_cast<const float * __restrict__>(&ceph_v2_cos_results);
     const float * __restrict__ p_cephes_cos_results    = reinterpret_cast<const float * __restrict__>(&cephes_cos_results);
     for(std::uint32_t i = 0;i != 15; ++i) 
     {
        const float ceph_v2_result   = p_ceph_v2_cos_results[i];
        const float cephes_result    = p_cephes_cos_results[i];
        const bool fcomp_result    = almostEqual2sComplement(ceph_v2_result,cephes_result,7);
        if(!fcomp_result) 
        {
            printf(ANSI_COLOR_RED "[PERF-TEST]: ***FAILED***, ceph_v2_cos=%.7f,cephes_cos=%.7f\n",ceph_v2_result,cephes_result);
            printNumber("simd_ceph_cosf_16xf32",ceph_v2_result,0);
            printNumber("_mm512_cephes_cosf_ps_v3",cephes_result,0);
        }
        else
        {
            printf(ANSI_COLOR_CYAN "[PERF-TEST]: ***CORRECT***, ceph_v2_cos=%.7f,cephes_cos=%.7f\n",ceph_v2_result,cephes_result);
            printNumber("simd_ceph_cosf_16xf32",ceph_v2_result,0);
            printNumber("_mm512_cephes_cosf_ps_v3",cephes_result,0);
        }
     }
}

int main()
{
    test_runner_omp_sections_1st_seq();
    return 0;
}