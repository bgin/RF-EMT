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
#include <signal.h>
#include <unistd.h>
#include <immintrin.h>
#include <omp.h> 
#include <algorithm>
#include <random>
#include <ctime>
#include <functional>
#include "GMS_malloc.h"
#include "GMS_sse_memset.h"
#include "GMS_sse_memcpy.h"

/*
    icpc -o perf_test_sse_memcpy_ps_64kiB -fp-model fast=2 -fno-exceptions -fopenmp -qopenmp -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5 \
    GMS_config.h  GMS_malloc.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_sse_memcpy.h GMS_sse_memcpy.cpp memmove-sse2-unaligned-erms.S perf_test_sse_memcpy_ps_64kiB.cpp  
    
   


   ASM: 
   icpc -S -fverbose-asm -masm=intel  -march=skylake-avx512 -mavx512f -fno-exceptions -falign-functions=32 -fopenmp -qopenmp\ 
   -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5\  
    GMS_config.h  GMS_malloc.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_sse_memcpy.h GMS_sse_memcpy.cpp memmove-sse2-unaligned-erms.S perf_test_sse_memcpy_ps_64kiB.cpp 
    

*/

#define BUFFER_STORE_SIZE 80
#define FORMAT_STORE_SIZE 80 

extern "C" void * __memcpy_sse2_unaligned(void * __restrict destination, const void * __restrict source, size_t size);

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
     printf("tid=%d affinity=%s\n",omp_get_thread_num(),buffer);
}




__attribute__((hot))
__attribute__((noinline))
void perf_test_memcpy_sse2_unaligned_64kiB(float * __restrict__ ,
                                          const float * __restrict__ ,
                                          unsigned __int64 * __restrict__ ,
                                          unsigned __int64 * __restrict__ ,
                                          unsigned __int64 * __restrict__ ,
                                          std::size_t,
                                          const int32_t ,
                                          const int32_t,
                                          uint32_t &);



void perf_test_memcpy_sse2_unaligned_64kiB(float * __restrict__ dst,
                                          const float * __restrict__ src,
                                          unsigned __int64 * __restrict__ glibc_memcpy_s,
                                          unsigned __int64 * __restrict__ glibc_memcpy_e,
                                          unsigned __int64 * __restrict__ glibc_memcpy_d,
                                          std::size_t sz,
                                          const int32_t n_runs,
                                          const int32_t n_samples,
                                          uint32_t & tid )
{
        constexpr std::size_t RDTSCP_LAT{42ull};
        const std::size_t data_size{reinterpret_cast<std::size_t>(sizeof(float)*sz)};
        [[maybe_unused]]
        volatile unsigned __int64 warmup_start;
        [[maybe_unused]]
        volatile unsigned __int64 warmup_end;
        uint32_t mem_start{0};
        uint32_t mem_end{0}; 

         // warmup of RDTSCP
        warmup_start = __rdtscp(&mem_start);
        warmup_end   = __rdtscp(&mem_end);
        // Warmup of __memcpy_sse2_unaligned
        __memcpy_sse2_unaligned((float*)&dst[0],(float*)&src[0],data_size);
     
     
        // Begin  looped timing execution 
        for(int32_t __i{0}; __i != n_runs; ++__i)
        {
            for(int32_t __j{0}; __j != n_samples; ++__j) 
            {   
                 __asm__ __volatile__ ("lfence");
                 unsigned __int64 start_curr{__rdtscp(&mem_start)};
                 __memcpy_sse2_unaligned((float*)&dst[0],(float*)&src[0],data_size);
                 unsigned __int64 end_curr{__rdtscp(&mem_end)};
                 __asm__ __volatile__ ("lfence");
            //remove latency
                 unsigned __int64 start_corrected{start_curr-RDTSCP_LAT};
                 glibc_memcpy_s[__i*n_samples+__j] = start_corrected;
                 unsigned __int64 end_corrected{end_curr-RDTSCP_LAT};
                 glibc_memcpy_e[__i*n_samples+__j]   = end_corrected;
                 glibc_memcpy_d[__i*n_samples+__j] = end_corrected-start_corrected;
            }
       }
       tid = omp_get_thread_num();
}

__attribute__((hot))
__attribute__((noinline))
void perf_test_sse_memcpy_unroll8x_ps_64kiB(float * __restrict__ ,
                                            float * __restrict__ ,
                                            unsigned __int64 * __restrict__ ,
                                            unsigned __int64 * __restrict__ ,
                                            unsigned __int64 * __restrict__ ,
                                            std::size_t,
                                            const int32_t ,
                                            const int32_t,
                                            uint32_t &);


void perf_test_sse_memcpy_unroll8x_ps_64kiB(float * __restrict__ dst,
                                            float * __restrict__ src,
                                            unsigned __int64 * __restrict__ my_memcpy_s,
                                            unsigned __int64 * __restrict__ my_memcpy_e,
                                            unsigned __int64 * __restrict__ my_memcpy_d,
                                            std::size_t sz,
                                            const int32_t n_runs,
                                            const int32_t n_samples,
                                            uint32_t & tid)
{
       using namespace gms::common;
       constexpr std::size_t RDTSCP_LAT{42ull};
       [[maybe_unused]]
       volatile unsigned __int64 warmup_start;
       [[maybe_unused]]
       volatile unsigned __int64 warmup_end;
       uint32_t mem_start{0};
       uint32_t mem_end{0};
          
       // warmup of RDTSCP
       warmup_start = __rdtscp(&mem_start);
       warmup_end   = __rdtscp(&mem_end);
       // Warmup of __memcpy_sse2_unaligned
       sse_memcpy_unroll8x_ps(&dst[0],&src[0],sz);

       // Begin  looped timing execution 
       for(int32_t __i{0}; __i != n_runs; ++__i)
       {
           for(int32_t __j{0}; __j != n_samples; ++__j) 
           {   
               __asm__ __volatile__ ("lfence");
               unsigned __int64 start_curr{__rdtscp(&mem_start)};
               sse_memcpy_unroll8x_ps(&dst[0],&src[0],sz);
               unsigned __int64 end_curr{__rdtscp(&mem_end)};
               __asm__ __volatile__ ("lfence");
               //remove latency
               unsigned __int64 start_corrected{start_curr-RDTSCP_LAT};
               my_memcpy_s[__i*n_samples+__j] = start_corrected;
               unsigned __int64 end_corrected{end_curr-RDTSCP_LAT};
               my_memcpy_e[__i*n_samples+__j]   = end_corrected;
               my_memcpy_d[__i*n_samples+__j] = end_corrected-start_corrected;
           }
       }
       tid = omp_get_thread_num();
    
}

// TID #0 handles glibc memcpy
__attribute__((hot))
void test_runner_omp_sections_1st_seq();

void test_runner_omp_sections_1st_seq()
{
    using namespace gms::common;
    constexpr std::size_t sz{16000ull};
    constexpr std::size_t nbytes{(std::size_t)sizeof(float)*sz};
    constexpr std::size_t alignment{16ull};
    constexpr float filler{3.14159265358979323846264338328F};
    constexpr int32_t n_runs{10};
    constexpr int32_t n_samples{50};
    constexpr int32_t n_total{n_runs*n_samples};
    unsigned __int64 glibc_memcpy_s[n_total] = {UINT64_MAX};
    unsigned __int64 glibc_memcpy_e[n_total] = {UINT64_MAX}; 
    unsigned __int64 glibc_memcpy_d[n_total] = {UINT64_MAX};
    unsigned __int64 my_memcpy_s[n_total]    = {UINT64_MAX};
    unsigned __int64 my_memcpy_e[n_total]    = {UINT64_MAX};
    unsigned __int64 my_memcpy_d[n_total]    = {UINT64_MAX};
    float * __restrict__ dst_my_memcpy{NULL};
    float * __restrict__ src_my_memcpy{NULL};
    float * __restrict__ dst_glibc_memcpy{NULL};
    float * __restrict__ src_glibc_memcpy{NULL};
    // Allocation 
    dst_my_memcpy    = reinterpret_cast<float*>(gms_mm_malloc(nbytes,alignment));
    src_my_memcpy    = reinterpret_cast<float*>(gms_mm_malloc(nbytes,alignment));
    dst_glibc_memcpy = reinterpret_cast<float*>(gms_mm_malloc(nbytes,alignment));
    src_glibc_memcpy = reinterpret_cast<float*>(gms_mm_malloc(nbytes,alignment));
    sse_memset_unroll8x_ps(&src_my_memcpy[0],filler,sz);
    sse_memset_unroll8x_ps(&src_glibc_memcpy[0],filler,sz);
    uint32_t glibc_memcpy_tid;
    uint32_t my_memcpy_tid;
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
        perf_test_memcpy_sse2_unaligned_64kiB(dst_glibc_memcpy,
                                              src_glibc_memcpy,
                                              glibc_memcpy_s,
                                              glibc_memcpy_e,
                                              glibc_memcpy_d,
                                              sz,
                                              n_runs,
                                              n_samples,
                                              glibc_memcpy_tid);
        print_thread_affinity();
    }

    #pragma omp section 
    {
        perf_test_sse_memcpy_unroll8x_ps_64kiB(dst_my_memcpy,
                                               src_my_memcpy,
                                               my_memcpy_s,
                                               my_memcpy_e,
                                               my_memcpy_d,
                                               sz,
                                               n_runs,
                                               n_samples,
                                               my_memcpy_tid);
        print_thread_affinity();
    }
}

     printf(ANSI_COLOR_GREEN "[PERF-TEST: -- size=%llu]: __memcpy_sse2_unaligned: Started -- dumping-results\n",nbytes );
     printf(ANSI_COLOR_GREEN "[PERF-TEST] -- Executed by Core=%d \n",glibc_memcpy_tid);
     for(int32_t __i{0}; __i != n_runs; ++__i)
     {
         for(int32_t __j{0}; __j != n_samples; ++__j) 
         {  
            unsigned __int64 s{glibc_memcpy_s[__i*n_samples+__j]};
            unsigned __int64 e{glibc_memcpy_e[__i*n_samples+__j]};
            unsigned __int64 d{glibc_memcpy_d[__i*n_samples+__j]};
            printf(ANSI_COLOR_GREEN "[PMC: RDTSCP] -- Run=%d, start=%llu,end=%llu,delta=%llu\n",__i,s,e,d);
         }
     }
     printf(ANSI_COLOR_GREEN "[PERF-TEST: -- size=%llu]: __memcpy_sse2_unaligned: Finished -- dumping-results" ANSI_RESET_ALL"\n\n",nbytes);

     printf(ANSI_COLOR_RED "[PERF-TEST: -- size=%llu]: sse_memcpy_unroll8x_ps: Started -- dumping-results\n",nbytes );
     printf(ANSI_COLOR_RED "[PERF-TEST] -- Executed by Core=%d \n",my_memcpy_tid);
     for(int32_t __i{0}; __i != n_runs; ++__i)
     {
         for(int32_t __j{0}; __j != n_samples; ++__j) 
         {  
            unsigned __int64 s{my_memcpy_s[__i*n_samples+__j]};
            unsigned __int64 e{my_memcpy_e[__i*n_samples+__j]};
            unsigned __int64 d{my_memcpy_d[__i*n_samples+__j]};
            printf(ANSI_COLOR_RED "[PMC: RDTSCP] -- Run=%d, start=%llu,end=%llu,delta=%llu\n",__i,s,e,d);
         }
     }
     printf(ANSI_COLOR_RED "[PERF-TEST: -- size=%llu]: sse_memcpy_unroll8x_ps: Finished -- dumping-results" ANSI_RESET_ALL"\n",nbytes);
     if(src_glibc_memcpy != NULL && nbytes > 0ull) gms_mm_free(src_glibc_memcpy);
     if(dst_glibc_memcpy != NULL && nbytes > 0ull) gms_mm_free(dst_glibc_memcpy);
     if(src_my_memcpy    != NULL && nbytes > 0ull) gms_mm_free(src_my_memcpy);
     if(dst_my_memcpy    != NULL && nbytes > 0ull) gms_mm_free(dst_my_memcpy);

}

__attribute__((hot))
void test_runner_single_thread(const int32_t,const int32_t);

void test_runner_single_thread(const int32_t cpu,const int32_t priority)
{
    using namespace gms::common;
    constexpr std::size_t sz{64000ull};
    constexpr std::size_t nbytes{(std::size_t)sizeof(float)*sz};
    constexpr std::size_t alignment{16ull};
    constexpr float filler{3.14159265358979323846264338328F};
    constexpr int32_t n_runs{10};
    constexpr int32_t n_samples{50};
    constexpr int32_t n_total{n_runs*n_samples};
    unsigned __int64 glibc_memcpy_s[n_total] = {UINT64_MAX};
    unsigned __int64 glibc_memcpy_e[n_total] = {UINT64_MAX}; 
    unsigned __int64 glibc_memcpy_d[n_total] = {UINT64_MAX};
    unsigned __int64 my_memcpy_s[n_total]    = {UINT64_MAX};
    unsigned __int64 my_memcpy_e[n_total]    = {UINT64_MAX};
    unsigned __int64 my_memcpy_d[n_total]    = {UINT64_MAX};
    float * __restrict__ dst_my_memcpy{NULL};
    float * __restrict__ src_my_memcpy{NULL};
    float * __restrict__ dst_glibc_memcpy{NULL};
    float * __restrict__ src_glibc_memcpy{NULL};
    // Allocation 
    dst_my_memcpy    = reinterpret_cast<float*>(gms_mm_malloc(nbytes,alignment));
    src_my_memcpy    = reinterpret_cast<float*>(gms_mm_malloc(nbytes,alignment));
    dst_glibc_memcpy = reinterpret_cast<float*>(gms_mm_malloc(nbytes,alignment));
    src_glibc_memcpy = reinterpret_cast<float*>(gms_mm_malloc(nbytes,alignment));
    sse_memset_unroll8x_ps(&src_my_memcpy[0],filler,sz);
    sse_memset_unroll8x_ps(&src_glibc_memcpy[0],filler,sz);
    uint32_t glibc_memcpy_tid;
    uint32_t my_memcpy_tid;
    int32_t status;
    status = set_affinity_and_priority(cpu,priority);
    if(status > 0)
    {
       printf("[***ERROR***]: set_affinity_and_priority status:%d, errno=%d\n",status,errno);
    }
    perf_test_memcpy_sse2_unaligned_64kiB(dst_glibc_memcpy,
                                              src_glibc_memcpy,
                                              glibc_memcpy_s,
                                              glibc_memcpy_e,
                                              glibc_memcpy_d,
                                              sz,
                                              n_runs,
                                              n_samples,
                                              glibc_memcpy_tid);
    print_thread_affinity();
    

    
    perf_test_sse_memcpy_unroll8x_ps_64kiB(dst_my_memcpy,
                                               src_my_memcpy,
                                               my_memcpy_s,
                                               my_memcpy_e,
                                               my_memcpy_d,
                                               sz,
                                               n_runs,
                                               n_samples,
                                               my_memcpy_tid);
    print_thread_affinity();
    


     printf(ANSI_COLOR_GREEN "[PERF-TEST: -- size=%llu]: __memcpy_sse2_unaligned: Started -- dumping-results\n",nbytes );
     printf(ANSI_COLOR_GREEN "[PERF-TEST] -- Executed by Core=%d \n",glibc_memcpy_tid);
     for(int32_t __i{0}; __i != n_runs; ++__i)
     {
         for(int32_t __j{0}; __j != n_samples; ++__j) 
         {  
            unsigned __int64 s{glibc_memcpy_s[__i*n_samples+__j]};
            unsigned __int64 e{glibc_memcpy_e[__i*n_samples+__j]};
            unsigned __int64 d{glibc_memcpy_d[__i*n_samples+__j]};
            printf(ANSI_COLOR_GREEN "[PMC: RDTSCP] -- Run=%d, start=%llu,end=%llu,delta=%llu\n",__i,s,e,d);
         }
     }
     printf(ANSI_COLOR_GREEN "[PERF-TEST: -- size=%llu]: __memcpy_sse2_unaligned: Finished -- dumping-results" ANSI_RESET_ALL"\n\n",nbytes);

     printf(ANSI_COLOR_RED "[PERF-TEST: -- size=%llu]: sse_memcpy_unroll8x_ps: Started -- dumping-results\n",nbytes );
     printf(ANSI_COLOR_RED "[PERF-TEST] -- Executed by Core=%d \n",my_memcpy_tid);
     for(int32_t __i{0}; __i != n_runs; ++__i)
     {
         for(int32_t __j{0}; __j != n_samples; ++__j) 
         {  
            unsigned __int64 s{my_memcpy_s[__i*n_samples+__j]};
            unsigned __int64 e{my_memcpy_e[__i*n_samples+__j]};
            unsigned __int64 d{my_memcpy_d[__i*n_samples+__j]};
            printf(ANSI_COLOR_RED "[PMC: RDTSCP] -- Run=%d, start=%llu,end=%llu,delta=%llu\n",__i,s,e,d);
         }
     }
     printf(ANSI_COLOR_RED "[PERF-TEST: -- size=%llu]: sse_memcpy_unroll8x_ps: Finished -- dumping-results" ANSI_RESET_ALL"\n",nbytes);
     if(src_glibc_memcpy != NULL && nbytes > 0ull) gms_mm_free(src_glibc_memcpy);
     if(dst_glibc_memcpy != NULL && nbytes > 0ull) gms_mm_free(dst_glibc_memcpy);
     if(src_my_memcpy    != NULL && nbytes > 0ull) gms_mm_free(src_my_memcpy);
     if(dst_my_memcpy    != NULL && nbytes > 0ull) gms_mm_free(dst_my_memcpy);

}

__attribute__((hot))
__attribute__((noinline))
void test_runner_omp_parallel();

void test_runner_omp_parallel()
{
const int32_t priority{99};
int32_t tid;
#pragma omp parallel shared(priority) private(tid) num_threads(omp_get_max_threads())
{
      tid = omp_get_thread_num();
      test_runner_single_thread(tid,priority);
}
}



int main()
{
    std::clock_t seed{};
    int32_t      which{-1};
    seed = std::clock();
    auto rand{std::bind(std::uniform_int_distribution<int32_t>{0,2},
                          std::mt19937(seed))};
    which = rand();
    switch (which)
    {
        case 0 : 
        test_runner_omp_sections_1st_seq();
        break;
        case 1 :
        test_runner_single_thread(5,99);
        break;
        case 2 : 
        test_runner_omp_parallel();
        break;
        default:
        printf("[**ERROR**] -- Invalid switch case=%d\n",which);
        std::terminate();
    }
    
    
    
    return 0;
}

