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
#include "GMS_ansi_colors.h"
#include "GMS_machine_utils.h"

/*
   icpc -o perf_test_rdtsc_overhead -std=c++17 -fp-model fast=2 -fopenmp -ftz -ggdb -ipo -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_machine_utils.h perf_test_rdtsc_overhead.cpp 
   ASM: 
   icpc -S -fverbose-asm -masm=intel -fopenmp -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_machine_utils.h perf_test_rdtsc_overhead.cpp 

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

typedef void (*PERF_TEST_RDTSC_FPTR) (std::uint64_t * __restrict__,std::uint64_t * __restrict__,std::uint64_t * __restrict__,const std::int32_t,const std::int32_t,std::uint32_t&);

__attribute__((hot))
__attribute__((noinline))
void perf_test_rdtsc_intrinsic(std::uint64_t * __restrict__ ,
                               std::uint64_t * __restrict__ ,
                               std::uint64_t * __restrict__ ,
                               const std::int32_t ,
                               const std::int32_t ,
                               std::uint32_t &);

__attribute__((hot))
__attribute__((noinline))
void perf_test_rdtsc_intrinsic(std::uint64_t * __restrict__ samples_s,
                               std::uint64_t * __restrict__ samples_e,
                               std::uint64_t * __restrict__ samples_d,
                               const std::int32_t n_runs,
                               const std::int32_t n_samples,
                               std::uint32_t &tid)
{
    tid = 9999;
    std::uint64_t * __restrict__ p_samples_s = samples_s;
    std::uint64_t * __restrict__ p_samples_e = samples_e;
    std::uint64_t * __restrict__ p_samples_d = samples_d;
    [[maybe_unused]] volatile unsigned __int64 warmup1_rdtsc_intrin{},
                                               warmup2_rdtsc_intrin{};
    warmup1_rdtsc_intrin    = __rdtsc();
    warmup2_rdtsc_intrin    = __rdtsc();
    for(std::int32_t i = 0; i < n_runs; ++i) 
    {   
        for(std::int32_t j = 0;j < n_samples; ++j)  
        {
            __asm__ __volatile__ ("lfence");
            register unsigned __int64 start = __rdtsc();
            register unsigned __int64 end   = __rdtsc();
            __asm__ __volatile__ ("mfence");
            p_samples_s[i*n_samples+j] = start;
            p_samples_e[i*n_samples+j] = end;
            p_samples_d[i*n_samples+j] = end-start;
        }
    }
    tid = omp_get_thread_num();
}

__attribute__((hot))
__attribute__((noinline))
void perf_test_rdtscp_intrinsic(std::uint64_t * __restrict__ ,
                               std::uint64_t * __restrict__ ,
                               std::uint64_t * __restrict__ ,
                               const std::int32_t ,
                               const std::int32_t ,
                               std::uint32_t &);

__attribute__((hot))
__attribute__((noinline))
void perf_test_rdtscp_intrinsic(std::uint64_t * __restrict__ samples_s,
                               std::uint64_t * __restrict__ samples_e,
                               std::uint64_t * __restrict__ samples_d,
                               const std::int32_t n_runs,
                               const std::int32_t n_samples,
                               std::uint32_t &tid)
{
    tid = 9999;
    std::uint64_t * __restrict__ p_samples_s = samples_s;
    std::uint64_t * __restrict__ p_samples_e = samples_e;
    std::uint64_t * __restrict__ p_samples_d = samples_d;
    [[maybe_unused]] volatile unsigned __int64 warmup1_rdtsc_intrin{},
                                               warmup2_rdtsc_intrin{};
    [[maybe_unused]] std::uint32_t dummy_mem0{},dummy_mem1{};
    warmup1_rdtsc_intrin    = __rdtscp(&dummy_mem0);
    warmup2_rdtsc_intrin    = __rdtscp(&dummy_mem1);
    for(std::int32_t i = 0; i < n_runs; ++i) 
    {   
        for(std::int32_t j = 0;j < n_samples; ++j)  
        {
            __asm__ __volatile__ ("lfence");
            register unsigned __int64 start = __rdtscp(&dummy_mem0);
            register unsigned __int64 end   = __rdtscp(&dummy_mem1);
            __asm__ __volatile__ ("mfence");
            gms::common::do_not_optimize(dummy_mem0);
            gms::common::do_not_optimize(dummy_mem1);
            p_samples_s[i*n_samples+j] = start;
            p_samples_e[i*n_samples+j] = end;
            p_samples_d[i*n_samples+j] = end-start;
        }
    }
    tid = omp_get_thread_num();
}

__attribute__((hot))
__attribute__((noinline))
void perf_test_rdtscp_fenced(std::uint64_t * __restrict__ ,
                               std::uint64_t * __restrict__ ,
                               std::uint64_t * __restrict__ ,
                               const std::int32_t ,
                               const std::int32_t ,
                               std::uint32_t &);

__attribute__((hot))
__attribute__((noinline))
void perf_test_rdtscp_fenced(std::uint64_t * __restrict__ samples_s,
                               std::uint64_t * __restrict__ samples_e,
                               std::uint64_t * __restrict__ samples_d,
                               const std::int32_t n_runs,
                               const std::int32_t n_samples,
                               std::uint32_t &tid)
{
    tid = 9999;
    std::uint64_t * __restrict__ p_samples_s = samples_s;
    std::uint64_t * __restrict__ p_samples_e = samples_e;
    std::uint64_t * __restrict__ p_samples_d = samples_d;
    [[maybe_unused]] volatile std::uint64_t warmup1_rdtscp_fenced{},
                                            warmup2_rdtscp_fenced{};
    warmup1_rdtscp_fenced    = gms::common::rdtscp_fenced();
    warmup2_rdtscp_fenced    = gms::common::rdtscp_fenced();
    for(std::int32_t i = 0; i < n_runs; ++i) 
    {   
        register std::int32_t outer_idx = i*n_samples;
        for(std::int32_t j = 0;j < n_samples; ++j)  
        {   
            register std::uint64_t start = gms::common::rdtscp_fenced();
            register std::uint64_t end   = gms::common::rdtscp_fenced();
            register std::int32_t inner_idx = outer_idx+j;
            p_samples_s[inner_idx] = start;
            p_samples_e[inner_idx] = end;
            p_samples_d[inner_idx] = end-start;
        }
    }
    tid = omp_get_thread_num();
}

__attribute__((hot))
__attribute__((noinline))
void perf_test_rdtsc_rdtscp_macro(std::uint64_t * __restrict__ ,
                                  std::uint64_t * __restrict__ ,
                                  std::uint64_t * __restrict__ ,
                                  const std::int32_t ,
                                  const std::int32_t ,
                                  std::uint32_t &);

__attribute__((hot))
__attribute__((noinline))
void perf_test_rdtsc_rdtscp_macro(std::uint64_t * __restrict__ samples_s,
                                  std::uint64_t * __restrict__ samples_e,
                                  std::uint64_t * __restrict__ samples_d,
                                  const std::int32_t n_runs,
                                  const std::int32_t n_samples,
                                  std::uint32_t &tid)
{
    tid = 9999;
    std::uint64_t * __restrict__ p_samples_s = samples_s;
    std::uint64_t * __restrict__ p_samples_e = samples_e;
    std::uint64_t * __restrict__ p_samples_d = samples_d;
    [[maybe_unused]] volatile std::uint64_t warmup_cycles_start{};
    [[maybe_unused]] volatile std::uint64_t warmup_cycles_end{};
    register std::uint64_t cycles_start{UINT64_MAX};
    register std::uint64_t cycles_end{UINT64_MAX};
    RDTSC_START_v1(warmup_cycles_start);
    RDTSC_STOP_v1(warmup_cycles_end);
    for(std::int32_t i = 0; i < n_runs; ++i) 
    {   
        register std::int32_t outer_idx = i*n_samples;
        for(std::int32_t j = 0;j < n_samples; ++j)  
        {   
            RDTSC_START_v1(cycles_start);
            RDTSC_STOP_v1(cycles_end);
            register std::int32_t inner_idx = outer_idx+j;
            p_samples_s[inner_idx] = cycles_start;
            p_samples_e[inner_idx] = cycles_end;
            p_samples_d[inner_idx] = cycles_end-cycles_start;
        }
    }
    tid = omp_get_thread_num();
}


__attribute__((hot))
__attribute__((noinline))
void perf_test_rdtsc_rdtscp_func(std::uint64_t * __restrict__ ,
                                  std::uint64_t * __restrict__ ,
                                  std::uint64_t * __restrict__ ,
                                  const std::int32_t ,
                                  const std::int32_t ,
                                  std::uint32_t &);

__attribute__((hot))
__attribute__((noinline))
void perf_test_rdtsc_rdtscp_func(std::uint64_t * __restrict__ samples_s,
                                  std::uint64_t * __restrict__ samples_e,
                                  std::uint64_t * __restrict__ samples_d,
                                  const std::int32_t n_runs,
                                  const std::int32_t n_samples,
                                  std::uint32_t &tid)
{
    tid = 9999;
    std::uint64_t * __restrict__ p_samples_s = samples_s;
    std::uint64_t * __restrict__ p_samples_e = samples_e;
    std::uint64_t * __restrict__ p_samples_d = samples_d;
    [[maybe_unused]] volatile std::uint64_t warmup_cycles_start{};
    [[maybe_unused]] volatile std::uint64_t warmup_cycles_end{};
    register std::uint64_t cycles_start{UINT64_MAX};
    register std::uint64_t cycles_end{UINT64_MAX};
    warmup_cycles_start = gms::common::rdtsc_serialized_start();
    warmup_cycles_end   = gms::common::rdtsc_serialized_stop();
    for(std::int32_t i = 0; i < n_runs; ++i) 
    {   
        register std::int32_t outer_idx = i*n_samples;
        for(std::int32_t j = 0;j < n_samples; ++j)  
        {   
            cycles_start = gms::common::rdtsc_serialized_start();
            cycles_end   = gms::common::rdtsc_serialized_stop();
            register std::int32_t inner_idx = outer_idx+j;
            p_samples_s[inner_idx] = cycles_start;
            p_samples_e[inner_idx] = cycles_end;
            p_samples_d[inner_idx] = cycles_end-cycles_start;
        }
    }
    tid = omp_get_thread_num();
}


__attribute__((hot))
void 
test_runner_single_thread(PERF_TEST_RDTSC_FPTR,
                          std::string &&,
                          const std::int32_t,
                          const std::int32_t);

void 
test_runner_single_thread(PERF_TEST_RDTSC_FPTR fptr,
                          std::string &&intrin_name,
                          const std::int32_t cpu_core,
                          const std::int32_t priority)
{
    constexpr std::int32_t n_runs{20};
    constexpr std::int32_t n_samples{100};
    constexpr std::int32_t tot_samples{n_runs*n_samples};
    __ATTR_ALIGN__(64) std::uint64_t samples_s_core_x[tot_samples] = {UINT64_MAX};
    __ATTR_ALIGN__(64) std::uint64_t samples_e_core_x[tot_samples] = {UINT64_MAX};
    __ATTR_ALIGN__(64) std::uint64_t samples_d_core_x[tot_samples] = {UINT64_MAX};
    //PERF_TEST_RDTSC_FPTR fptr = nullptr;
    [[maybe_unused]] std::int32_t printf_ret{};
    std::uint32_t core_x_tid{9999};
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
    std::int32_t prio_stat{};
    prio_stat = set_affinity_and_priority(cpu_core,priority);
    if(prio_stat>0)
    {
       std::int32_t which = PRIO_PROCESS;
       id_t pid      = getpid();
       errno         = 0;
       std::int32_t nice  = getpriority(which,pid);
       printf_ret = printf("[***ERROR***]: set_affinity_and_priority status:%d,errno=%d,pid=%d,nice=%d\n",prio_stat,errno,pid,nice);
    }
    //fptr = perf_test_rdtsc_intrinsic;
    printf_ret = printf(ANSI_COLOR_WHITE "[PERF-TEST]: START of %s -- overhead measurement.\n",intrin_name.c_str());
    fptr(&samples_s_core_x[0],&samples_e_core_x[0],&samples_d_core_x[0],
             n_runs,n_samples,core_x_tid);
    print_thread_affinity();
    printf_ret = printf(ANSI_COLOR_WHITE "[PERF-TEST]: The results of '__rdtsc' overhead measurement.\n");
    for(std::int32_t i{0}; i != n_runs; ++i) 
    {
        for(std::int32_t j{0}; j != n_samples; ++j)   
        {
            const std::int32_t loop_idx{i*n_samples+j};
            std::uint64_t core_x_s{samples_s_core_x[loop_idx]};
            std::uint64_t core_x_e{samples_e_core_x[loop_idx]};
            std::uint64_t core_x_d{samples_d_core_x[loop_idx]};
            printf_ret = printf(ANSI_COLOR_GREEN "[PERF-TEST]: X_TID=%d,s=%llu,e=%llu,d=%llu\n",core_x_tid,core_x_s,core_x_e,core_x_d);
                                                                                                            
        }
    }
    printf_ret = printf(ANSI_COLOR_WHITE "[PERF-TEST]: END of %s -- overhead measurement.\n",intrin_name.c_str());
}

__attribute__((hot))
void 
test_runner_omp_2_sections(PERF_TEST_RDTSC_FPTR,
                           std::string &&);

void 
test_runner_omp_2_sections(PERF_TEST_RDTSC_FPTR fptr,
                           std::string &&intrin_name)
{
    constexpr std::int32_t n_runs{20};
    constexpr std::int32_t n_samples{100};
    constexpr std::int32_t tot_samples{n_runs*n_samples};
    __ATTR_ALIGN__(64) std::uint64_t samples_s_core_x[tot_samples] = {UINT64_MAX};
    __ATTR_ALIGN__(64) std::uint64_t samples_e_core_x[tot_samples] = {UINT64_MAX};
    __ATTR_ALIGN__(64) std::uint64_t samples_d_core_x[tot_samples] = {UINT64_MAX};
    __ATTR_ALIGN__(64) std::uint64_t samples_s_core_y[tot_samples] = {UINT64_MAX};
    __ATTR_ALIGN__(64) std::uint64_t samples_e_core_y[tot_samples] = {UINT64_MAX};
    __ATTR_ALIGN__(64) std::uint64_t samples_d_core_y[tot_samples] = {UINT64_MAX};
    //PERF_TEST_RDTSC_FPTR fptr = nullptr;
    [[maybe_unused]] std::int32_t printf_ret{};
    std::uint32_t core_x_tid{9999};
    std::uint32_t core_y_tid{9999};
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
    //fptr = perf_test_rdtsc_intrinsic;
    printf_ret = printf(ANSI_COLOR_WHITE "[PERF-TEST]: START of %s -- overhead measurement.\n",intrin_name.c_str());
#pragma omp parallel sections 
{
    #pragma omp section 
    {
        fptr(&samples_s_core_x[0],&samples_e_core_x[0],&samples_d_core_x[0],
             n_runs,n_samples,core_x_tid);
        print_thread_affinity();
    }

    #pragma omp section 
    {
        fptr(&samples_s_core_y[0],&samples_e_core_y[0],&samples_d_core_y[0],
             n_runs,n_samples,core_y_tid);
        print_thread_affinity();
    }
}

    printf_ret = printf(ANSI_COLOR_WHITE "[PERF-TEST]: The results of '__rdtsc' overhead measurement.\n");
    for(std::int32_t i{0}; i != n_runs; ++i) 
    {
        for(std::int32_t j{0}; j != n_samples; ++j)   
        {
            const std::int32_t loop_idx{i*n_samples+j};
            std::uint64_t core_x_s{samples_s_core_x[loop_idx]};
            std::uint64_t core_y_s{samples_s_core_y[loop_idx]};
            std::uint64_t core_x_e{samples_e_core_x[loop_idx]};
            std::uint64_t core_y_e{samples_e_core_y[loop_idx]};
            std::uint64_t core_x_d{samples_d_core_x[loop_idx]};
            std::uint64_t core_y_d{samples_d_core_y[loop_idx]};
            printf_ret = printf(ANSI_COLOR_GREEN "[PERF-TEST]: X_TID=%d,s=%llu,e=%llu,d=%llu,Y_TID=%d,s=%llu,e=%llu,d=%llu\n",core_x_tid,core_x_s,core_x_e,core_x_d,
                                                                                                             core_y_tid,core_y_s,core_y_e,core_y_d);
        }
    }
    printf_ret = printf(ANSI_COLOR_WHITE "[PERF-TEST]: END of %s -- overhead measurement.\n",intrin_name.c_str());
}

int main()
{
    test_runner_omp_2_sections(&perf_test_rdtsc_intrinsic,std::string("__rdtsc"));
    test_runner_omp_2_sections(&perf_test_rdtscp_intrinsic,std::string("__rdtscp"));
    test_runner_omp_2_sections(&perf_test_rdtscp_fenced,std::string("rdtscp_fenced"));
    test_runner_omp_2_sections(&perf_test_rdtsc_rdtscp_macro,std::string("__rdtsc(p)_serialized-macro"));
    test_runner_omp_2_sections(&perf_test_rdtsc_rdtscp_func,std::string("__rdtsc(p)_serialized-function"));
    return 0;
}

