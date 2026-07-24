
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
#include "GMS_machine_utils.h"
#include "GMS_dsp_simd_fix_float_kernels.h"

/*
   icpc -o perf_test_dsp_simd_fix_float_kernels -fopenmp -std=c++17 -fp-model fast=2 -ftz -ggdb -ipo -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_machine_utils.h GMS_dsp_simd_fix_float_kernels.h perf_test_dsp_simd_fix_float_kernels.cpp 
   ASM: 
   icpc -S -fverbose-asm -masm=intel -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_machine_utils.h GMS_dsp_simd_fix_float_kernels.h unit_test_dsp_simd_fix_float_kernels.cpp
*/

#define BUFFER_STORE_SIZE 80
#define FORMAT_STORE_SIZE 80 

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

static inline 
void 
printNumber(const char * val_name,float f, const std::int32_t offset) 
{
  (*((int *)&f)) += offset;
  std::printf("value: %s,%+1.11g,0x%08X,%d\n",val_name,f, *(int *)&f, *(int *)&f);
}

__attribute__((hot))
__attribute__((aligned(32)))
void perf_test_unpack_zmm16r4_16x8_optimized(const float * __restrict__,
                                             float * __restrict__,
                                             std::uint64_t * __restrict__,
                                             std::uint64_t * __restrict__,
                                             std::uint64_t * __restrict__,
                                             const std::int32_t,
                                             const std::int32_t,
                                             std::uint32_t &);

__attribute__((hot))
__attribute__((aligned(32)))
void perf_test_unpack_zmm16r4_16x8_optimized(const float * __restrict__ data_in,
                                             float * __restrict__ data_out,
                                             std::uint64_t * __restrict__ optim_tsc_s,
                                             std::uint64_t * __restrict__ optim_tsc_e,
                                             std::uint64_t * __restrict__ optim_tsc_d,
                                             const std::int32_t n_runs,
                                             const std::int32_t n_samples,
                                             std::uint32_t & tid)
{
    constexpr std::size_t approx_meter_lat{31ull};
    const float * __restrict__ ptr_data_in     = data_in;
    float * __restrict__ ptr_data_out          = data_out;
    std::uint64_t * __restrict__ p_optim_tsc_s = optim_tsc_s;
    std::uint64_t * __restrict__ p_optim_tsc_e = optim_tsc_e;
    std::uint64_t * __restrict__ p_optim_tsc_d = optim_tsc_d;
    [[maybe_unused]] volatile std::uint64_t rdtsc_warmup_start{};
    [[maybe_unused]] volatile std::uint64_t rdtsc_warmup_end{};
    rdtsc_warmup_start = gms::common::rdtsc_serialized_start();
    rdtsc_warmup_end   = gms::common::rdtsc_serialized_stop();
    //warmup call to prefetch the rodata stored data.
    // during the subsequent calls the prefetching branch is disabled.
    gms::l1_phy::unpack_zmm16r4_16x8<true>(ptr_data_in,ptr_data_out);
    for(std::int32_t i{0}; i < n_runs; ++i) 
    {
        const std::int32_t outer_idx{i*n_samples};
        for(std::int32_t j{0}; j < n_samples; ++j)   
        {
            const std::int32_t inner_idx{outer_idx+j};
            register std::uint64_t start{gms::common::rdtsc_serialized_start()};
            gms::l1_phy::unpack_zmm16r4_16x8<false>(ptr_data_in,ptr_data_out);
            register std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            // remove approx latency
            std::uint64_t start_lat_corr{start-approx_meter_lat};
            p_optim_tsc_s[inner_idx] = start_lat_corr;
            std::uint64_t end_lat_corr{end-approx_meter_lat};
            p_optim_tsc_e[inner_idx] = end_lat_corr;
            p_optim_tsc_d[inner_idx] = end_lat_corr-start_lat_corr;
        }
    }
}

__attribute__((hot))
__attribute__((aligned(32)))
void perf_test_unpack_zmm16r4_16x8_unoptimized(const float * __restrict__,
                                             float * __restrict__,
                                             std::uint64_t * __restrict__,
                                             std::uint64_t * __restrict__,
                                             std::uint64_t * __restrict__,
                                             const std::int32_t,
                                             const std::int32_t,
                                             std::uint32_t &);

__attribute__((hot))
__attribute__((aligned(32)))
void perf_test_unpack_zmm16r4_16x8_unoptimized(const float * __restrict__ data_in,
                                             float * __restrict__ data_out,
                                             std::uint64_t * __restrict__ unoptim_tsc_s,
                                             std::uint64_t * __restrict__ unoptim_tsc_e,
                                             std::uint64_t * __restrict__ unoptim_tsc_d,
                                             const std::int32_t n_runs,
                                             const std::int32_t n_samples,
                                             std::uint32_t & tid)
{
    constexpr std::size_t approx_meter_lat{31ull};
    const float * __restrict__ ptr_data_in     = data_in;
    float * __restrict__ ptr_data_out          = data_out;
    std::uint64_t * __restrict__ p_unoptim_tsc_s = unoptim_tsc_s;
    std::uint64_t * __restrict__ p_unoptim_tsc_e = unoptim_tsc_e;
    std::uint64_t * __restrict__ p_unoptim_tsc_d = unoptim_tsc_d;
    [[maybe_unused]] volatile std::uint64_t rdtsc_warmup_start{};
    [[maybe_unused]] volatile std::uint64_t rdtsc_warmup_end{};
    rdtsc_warmup_start = gms::common::rdtsc_serialized_start();
    rdtsc_warmup_end   = gms::common::rdtsc_serialized_stop();
    // Warmup call.
    gms::l1_phy::unpack_zmm16r4_16x8<false>(ptr_data_in,ptr_data_out);
    for(std::int32_t i{0}; i < n_runs; ++i) 
    {
        const std::int32_t outer_idx{i*n_samples};
        for(std::int32_t j{0}; j < n_samples; ++j)   
        {
            const std::int32_t inner_idx{outer_idx+j};
            register std::uint64_t start{gms::common::rdtsc_serialized_start()};
            gms::l1_phy::unpack_zmm16r4_16x8<false>(ptr_data_in,ptr_data_out);
            register std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            // remove approx latency
            std::uint64_t start_lat_corr{start-approx_meter_lat};
            p_unoptim_tsc_s[inner_idx] = start_lat_corr;
            std::uint64_t end_lat_corr{end-approx_meter_lat};
            p_unoptim_tsc_e[inner_idx] = end_lat_corr;
            p_unoptim_tsc_d[inner_idx] = end_lat_corr-start_lat_corr;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////

__attribute__((hot))
__attribute__((aligned(32)))
void perf_test_transpose_zmm16r4_16x16_optimized(const float * __restrict__,
                                             float * __restrict__,
                                             std::uint64_t * __restrict__,
                                             std::uint64_t * __restrict__,
                                             std::uint64_t * __restrict__,
                                             const std::int32_t,
                                             const std::int32_t,
                                             std::uint32_t &);

__attribute__((hot))
__attribute__((aligned(32)))
void perf_test_transpose_zmm16r4_16x16_optimized(const float * __restrict__ data_in,
                                             float * __restrict__ data_out,
                                             std::uint64_t * __restrict__ optim_tsc_s,
                                             std::uint64_t * __restrict__ optim_tsc_e,
                                             std::uint64_t * __restrict__ optim_tsc_d,
                                             const std::int32_t n_runs,
                                             const std::int32_t n_samples,
                                             std::uint32_t & tid)
{
    constexpr std::size_t approx_meter_lat{31ull};
    const float * __restrict__ ptr_data_in     = data_in;
    float * __restrict__ ptr_data_out          = data_out;
    std::uint64_t * __restrict__ p_optim_tsc_s = optim_tsc_s;
    std::uint64_t * __restrict__ p_optim_tsc_e = optim_tsc_e;
    std::uint64_t * __restrict__ p_optim_tsc_d = optim_tsc_d;
    [[maybe_unused]] volatile std::uint64_t rdtsc_warmup_start{};
    [[maybe_unused]] volatile std::uint64_t rdtsc_warmup_end{};
    rdtsc_warmup_start = gms::common::rdtsc_serialized_start();
    rdtsc_warmup_end   = gms::common::rdtsc_serialized_stop();
    //warmup call to prefetch the rodata stored data.
    // during the subsequent calls the prefetching branch is disabled.
    gms::l1_phy::transpose_zmm16r4_16x16<true>(ptr_data_in,ptr_data_out);
    for(std::int32_t i{0}; i < n_runs; ++i) 
    {
        const std::int32_t outer_idx{i*n_samples};
        for(std::int32_t j{0}; j < n_samples; ++j)   
        {
            const std::int32_t inner_idx{outer_idx+j};
            register std::uint64_t start{gms::common::rdtsc_serialized_start()};
            gms::l1_phy::transpose_zmm16r4_16x16<false>(ptr_data_in,ptr_data_out);
            register std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            // remove approx latency
            std::uint64_t start_lat_corr{start-approx_meter_lat};
            p_optim_tsc_s[inner_idx] = start_lat_corr;
            std::uint64_t end_lat_corr{end-approx_meter_lat};
            p_optim_tsc_e[inner_idx] = end_lat_corr;
            p_optim_tsc_d[inner_idx] = end_lat_corr-start_lat_corr;
        }
    }
}

__attribute__((hot))
__attribute__((aligned(32)))
void perf_test_transpose_zmm16r4_16x16_unoptimized(const float * __restrict__,
                                             float * __restrict__,
                                             std::uint64_t * __restrict__,
                                             std::uint64_t * __restrict__,
                                             std::uint64_t * __restrict__,
                                             const std::int32_t,
                                             const std::int32_t,
                                             std::uint32_t &);

__attribute__((hot))
__attribute__((aligned(32)))
void perf_test_transpose_zmm16r4_16x16_unoptimized(const float * __restrict__ data_in,
                                             float * __restrict__ data_out,
                                             std::uint64_t * __restrict__ unoptim_tsc_s,
                                             std::uint64_t * __restrict__ unoptim_tsc_e,
                                             std::uint64_t * __restrict__ unoptim_tsc_d,
                                             const std::int32_t n_runs,
                                             const std::int32_t n_samples,
                                             std::uint32_t & tid)
{
    constexpr std::size_t approx_meter_lat{31ull};
    const float * __restrict__ ptr_data_in     = data_in;
    float * __restrict__ ptr_data_out          = data_out;
    std::uint64_t * __restrict__ p_unoptim_tsc_s = unoptim_tsc_s;
    std::uint64_t * __restrict__ p_unoptim_tsc_e = unoptim_tsc_e;
    std::uint64_t * __restrict__ p_unoptim_tsc_d = unoptim_tsc_d;
    [[maybe_unused]] volatile std::uint64_t rdtsc_warmup_start{};
    [[maybe_unused]] volatile std::uint64_t rdtsc_warmup_end{};
    rdtsc_warmup_start = gms::common::rdtsc_serialized_start();
    rdtsc_warmup_end   = gms::common::rdtsc_serialized_stop();
    // Warmup call.
    gms::l1_phy::transpose_zmm16r4_16x16<false>(ptr_data_in,ptr_data_out);
    for(std::int32_t i{0}; i < n_runs; ++i) 
    {
        const std::int32_t outer_idx{i*n_samples};
        for(std::int32_t j{0}; j < n_samples; ++j)   
        {
            const std::int32_t inner_idx{outer_idx+j};
            register std::uint64_t start{gms::common::rdtsc_serialized_start()};
            gms::l1_phy::transpose_zmm16r4_16x16<false>(ptr_data_in,ptr_data_out);
            register std::uint64_t end{gms::common::rdtsc_serialized_stop()};
            // remove approx latency
            std::uint64_t start_lat_corr{start-approx_meter_lat};
            p_unoptim_tsc_s[inner_idx] = start_lat_corr;
            std::uint64_t end_lat_corr{end-approx_meter_lat};
            p_unoptim_tsc_e[inner_idx] = end_lat_corr;
            p_unoptim_tsc_d[inner_idx] = end_lat_corr-start_lat_corr;
        }
    }
}

__attribute__((hot))
__attribute__((aligned(32)))
void test_runner_omp_sections_1st_seq();

__attribute__((hot))
__attribute__((aligned(32)))
void test_runner_omp_sections_1st_seq()
{
    constexpr std::int32_t n_runs{20};
    constexpr std::int32_t n_samples{100};
    constexpr std::int32_t nfloats{128};
    // prefetching is on
    std::uint64_t optim_res_s[n_runs*n_samples] = {UINT64_MAX};
    std::uint64_t optim_res_e[n_runs*n_samples] = {UINT64_MAX};
    std::uint64_t optim_res_d[n_runs*n_samples] = {UINT64_MAX};
    // prefetching is off 
    std::uint64_t unoptim_res_s[n_runs*n_samples] = {UINT64_MAX};
    std::uint64_t unoptim_res_e[n_runs*n_samples] = {UINT64_MAX};
    std::uint64_t unoptim_res_d[n_runs*n_samples] = {UINT64_MAX};
    __ATTR_ALIGN__(64) float buf_data_in[nfloats];
    __ATTR_ALIGN__(64) float optim_data_out[nfloats];
    __ATTR_ALIGN__(64) float unoptim_data_out[nfloats];
    unsigned long long seed{0ull};
    [[maybe_unused]] std::int32_t printf_ret{};
    seed = __rdtsc();
    auto rand_re{std::bind(std::uniform_real_distribution<float>(0.001f,5.0f),std::mt19937(seed))};
    for(std::int32_t i{0}; i < nfloats; ++i) 
    {
        buf_data_in[i] = rand_re();
    }
    const float * __restrict__ ptr_data_in{&buf_data_in[0]};
    float * __restrict__ ptr_optim_data_out{&optim_data_out[0]};
    float * __restrict__ ptr_unoptim_data_out{&unoptim_data_out[0]};
    uint32_t optim_tid;
    uint32_t unoptim_tid;
    int32_t setenv_ret;
    setenv_ret = setenv("OMP_PROC_BIND","true",1);
    if(setenv_ret==-1)
    {
        printf_ret = printf("[**ERROR**]: -- setenv reported an error=%d\n",setenv_ret);
    }
    setenv_ret = setenv("OMP_PROC_BIND","spread",1);
    if(setenv_ret==-1)
    {
        printf_ret = printf("[**ERROR**]: -- setenv reported an error=%d\n",setenv_ret);
    }

#pragma omp parallel sections 
{
    #pragma omp section 
    {
        perf_test_unpack_zmm16r4_16x8_optimized(ptr_data_in,ptr_optim_data_out,
                                                &optim_res_s[0],&optim_res_e[0],&optim_res_d[0],
                                                n_runs,n_samples,optim_tid);
         print_thread_affinity();
    }

    #pragma omp section 
    {
        perf_test_unpack_zmm16r4_16x8_unoptimized(ptr_data_in,ptr_unoptim_data_out,
                                                 &unoptim_res_s[0],&unoptim_res_e[0],&unoptim_res_d[0],
                                                 n_runs,n_samples,unoptim_tid);
        print_thread_affinity();
    }
}

    printf_ret = printf(ANSI_COLOR_GREEN "[PERF-TEST]:  unpack_zmm16r4_16x8<true>: RUN-COMPLETED!!\n");
    printf_ret = printf(ANSI_COLOR_GREEN "[PERF-TEST] -- Executed by Core=%d \n",optim_tid);
    for(int32_t i{0}; i < n_runs; ++i)
    {
         register std::int32_t outer_idx = i*n_samples;
         for(int32_t j{0}; j < n_samples; ++j) 
         {  
            register std::int32_t inner_idx = outer_idx+j;
            std::uint64_t s{optim_res_s[inner_idx]};
            std::uint64_t e{optim_res_e[inner_idx]};
            std::uint64_t d{optim_res_d[inner_idx]};
            printf(ANSI_COLOR_GREEN "[PMC: RDTSCP] -- Run=%d, start=%llu,end=%llu,delta=%llu\n",i,s,e,d);
         }
    }
    printf_ret = printf(ANSI_COLOR_GREEN "[PERF-TEST]: unpack_zmm16r4_16x8<true>: Finished -- dumping-results" ANSI_RESET_ALL"\n\n");
    
    printf_ret = printf(ANSI_COLOR_WHITE "[PERF-TEST]:  unpack_zmm16r4_16x8<false>: RUN-COMPLETED!!\n");
    printf_ret = printf(ANSI_COLOR_WHITE "[PERF-TEST] -- Executed by Core=%d \n",unoptim_tid);
    for(int32_t i{0}; i < n_runs; ++i)
    {
        register std::int32_t outer_idx = i*n_samples;
        for(int32_t j{0}; j < n_samples; ++j) 
        {  
            register std::int32_t inner_idx = outer_idx+j;
            std::uint64_t s{unoptim_res_s[inner_idx]};
            std::uint64_t e{unoptim_res_e[inner_idx]};
            std::uint64_t d{unoptim_res_d[inner_idx]};
            printf(ANSI_COLOR_WHITE "[PMC: RDTSCP] -- Run=%d, start=%llu,end=%llu,delta=%llu\n",i,s,e,d);
        }
    }
    printf_ret = printf(ANSI_COLOR_GREEN "[PERF-TEST]: unpack_zmm16r4_16x8<false>: Finished -- dumping-results" ANSI_RESET_ALL"\n\n");
}

__attribute__((hot))
__attribute__((aligned(32)))
void test_runner_omp_sections_1st_seq_transpose_16x16_krnl();

__attribute__((hot))
__attribute__((aligned(32)))
void test_runner_omp_sections_1st_seq_transpose_16x16_krnl()
{
    constexpr std::int32_t n_runs{20};
    constexpr std::int32_t n_samples{100};
    constexpr std::int32_t nfloats{256};
    // prefetching is on
    std::uint64_t optim_res_s[n_runs*n_samples] = {UINT64_MAX};
    std::uint64_t optim_res_e[n_runs*n_samples] = {UINT64_MAX};
    std::uint64_t optim_res_d[n_runs*n_samples] = {UINT64_MAX};
    // prefetching is off 
    std::uint64_t unoptim_res_s[n_runs*n_samples] = {UINT64_MAX};
    std::uint64_t unoptim_res_e[n_runs*n_samples] = {UINT64_MAX};
    std::uint64_t unoptim_res_d[n_runs*n_samples] = {UINT64_MAX};
    __ATTR_ALIGN__(64) float buf_data_in[nfloats];
    __ATTR_ALIGN__(64) float optim_data_out[nfloats];
    __ATTR_ALIGN__(64) float unoptim_data_out[nfloats];
    unsigned long long seed{0ull};
    [[maybe_unused]] std::int32_t printf_ret{};
    seed = __rdtsc();
    auto rand_re{std::bind(std::uniform_real_distribution<float>(0.001f,5.0f),std::mt19937(seed))};
    for(std::int32_t i{0}; i < nfloats; ++i) 
    {
        buf_data_in[i] = rand_re();
    }
    const float * __restrict__ ptr_data_in{&buf_data_in[0]};
    float * __restrict__ ptr_optim_data_out{&optim_data_out[0]};
    float * __restrict__ ptr_unoptim_data_out{&unoptim_data_out[0]};
    uint32_t optim_tid;
    uint32_t unoptim_tid;
    int32_t setenv_ret;
    setenv_ret = setenv("OMP_PROC_BIND","true",1);
    if(setenv_ret==-1)
    {
        printf_ret = printf("[**ERROR**]: -- setenv reported an error=%d\n",setenv_ret);
    }
    setenv_ret = setenv("OMP_PROC_BIND","spread",1);
    if(setenv_ret==-1)
    {
        printf_ret = printf("[**ERROR**]: -- setenv reported an error=%d\n",setenv_ret);
    }

#pragma omp parallel sections 
{
    #pragma omp section 
    {
        perf_test_transpose_zmm16r4_16x16_optimized(ptr_data_in,ptr_optim_data_out,
                                                &optim_res_s[0],&optim_res_e[0],&optim_res_d[0],
                                                n_runs,n_samples,optim_tid);
         print_thread_affinity();
    }

    #pragma omp section 
    {
        perf_test_transpose_zmm16r4_16x16_unoptimized(ptr_data_in,ptr_unoptim_data_out,
                                                 &unoptim_res_s[0],&unoptim_res_e[0],&unoptim_res_d[0],
                                                 n_runs,n_samples,unoptim_tid);
        print_thread_affinity();
    }
}

    printf_ret = printf(ANSI_COLOR_GREEN "[PERF-TEST]:  transpose_zmm16r4_16x16<true>: RUN-COMPLETED!!\n");
    printf_ret = printf(ANSI_COLOR_GREEN "[PERF-TEST] -- Executed by Core=%d \n",optim_tid);
    for(int32_t i{0}; i < n_runs; ++i)
    {
         register std::int32_t outer_idx = i*n_samples;
         for(int32_t j{0}; j < n_samples; ++j) 
         {  
            register std::int32_t inner_idx = outer_idx+j;
            std::uint64_t s{optim_res_s[inner_idx]};
            std::uint64_t e{optim_res_e[inner_idx]};
            std::uint64_t d{optim_res_d[inner_idx]};
            printf(ANSI_COLOR_GREEN "[PMC: RDTSCP] -- Run=%d, start=%llu,end=%llu,delta=%llu\n",i,s,e,d);
         }
    }
    printf_ret = printf(ANSI_COLOR_GREEN "[PERF-TEST]: transpose_zmm16r4_16x16<true>: Finished -- dumping-results" ANSI_RESET_ALL"\n\n");
    
    printf_ret = printf(ANSI_COLOR_WHITE "[PERF-TEST]:  transpose_zmm16r4_16x16<false>: RUN-COMPLETED!!\n");
    printf_ret = printf(ANSI_COLOR_WHITE "[PERF-TEST] -- Executed by Core=%d \n",unoptim_tid);
    for(int32_t i{0}; i < n_runs; ++i)
    {
        register std::int32_t outer_idx = i*n_samples;
        for(int32_t j{0}; j < n_samples; ++j) 
        {  
            register std::int32_t inner_idx = outer_idx+j;
            std::uint64_t s{unoptim_res_s[inner_idx]};
            std::uint64_t e{unoptim_res_e[inner_idx]};
            std::uint64_t d{unoptim_res_d[inner_idx]};
            printf(ANSI_COLOR_WHITE "[PMC: RDTSCP] -- Run=%d, start=%llu,end=%llu,delta=%llu\n",i,s,e,d);
        }
    }
    printf_ret = printf(ANSI_COLOR_GREEN "[PERF-TEST]: transpose_zmm16r4_16x16<false>: Finished -- dumping-results" ANSI_RESET_ALL"\n\n");
}

int main()
{
    (void)gms::common::check_prefetch_data_size();
    (void)test_runner_omp_sections_1st_seq();
    (void)test_runner_omp_sections_1st_seq_transpose_16x16_krnl();
    return 0;
}



