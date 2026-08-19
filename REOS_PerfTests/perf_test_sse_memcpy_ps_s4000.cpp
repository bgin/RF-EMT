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
#include <immintrin.h>
#include <omp.h> 
#include <algorithm>
#include <random>
#include <ctime>
#include "GMS_sse_memset.h"
#include "GMS_sse_memcpy.h"

/*
    icpc -o perf_test_sse_memcpy_ps_size4000 -fp-model fast=2 -fno-exceptions -fopenmp -qopenmp\ 
   -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5\  
    GMS_config.h  GMS_sse_memset.h GMS_sse_memset.cpp GMS_sse_memcpy.h GMS_sse_memcpy.cpp memmove-sse2-unaligned-erms.S perf_test_sse_memcpy_ps_size4000.cpp  
    
   


   ASM: 
   icpc -S -fverbose-asm -masm=intel  -march=skylake-avx512 -mavx512f -fno-exceptions -falign-functions=32 -fopenmp -qopenmp\ 
   -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5\  
    GMS_config.h  GMS_sse_memset.h GMS_sse_memset.cpp GMS_sse_memcpy.h GMS_sse_memcpy.cpp memmove-sse2-unaligned-erms.S perf_test_sse_memcpy_ps_size4000.cpp 
    

*/

extern "C" void * __memcpy_sse2_unaligned(void * __restrict destination, const void * __restrict source, size_t size);

// No correctness checking implemented.
__attribute__((noinline))
void perf_test_memcpy_sse2_unaligned_size4000(unsigned __int64 * __restrict__ ,
                                              unsigned __int64 * __restrict__ ,
                                              unsigned __int64 * __restrict__ ,
                                              const int32_t ,
                                              const int32_t,
                                              uint32_t & tid);

void perf_test_memcpy_sse2_unaligned_size4000(unsigned __int64 * __restrict__ glibc_memcpy_s,
                                              unsigned __int64 * __restrict__ glibc_memcpy_e,
                                              unsigned __int64 * __restrict__ glibc_memcpy_d,
                                              const int32_t n_runs,
                                              const int32_t n_samples,
                                              uint32_t & tid)
{
    
     constexpr std::size_t RDTSCP_LAT{42ull};
     constexpr std::size_t size_4000{4000ull};
     [[maybe_unused]]
     volatile unsigned __int64 warmup_start;
     [[maybe_unused]]
     volatile unsigned __int64 warmup_end;
     uint32_t mem_start{0};
     uint32_t mem_end{0};
     float data4000_src[size_4000];
     float data4000_dst[size_4000]  = {}; 
     const float filler{3.14159265358979323846264338328F};
     std::fill((float*)&data4000_src[0],(float*)&data4000_src[0]+size_4000,filler);
     // warmup of RDTSCP
     warmup_start = __rdtscp(&mem_start);
     warmup_end   = __rdtscp(&mem_end);
     // Warmup of __memcpy_sse2_unaligned
     __memcpy_sse2_unaligned((float*)&data4000_dst[0],(float*)&data4000_src[0],(std::size_t)(sizeof(float)*size_4000));
     
     
     // Begin  looped timing execution 
     for(int32_t __i{0}; __i != n_runs; ++__i)
     {
         for(int32_t __j{0}; __j != n_samples; ++__j) 
         {   
            __asm__ __volatile__ ("lfence");
             unsigned __int64 start_curr{__rdtscp(&mem_start)};
           // __asm__ __volatile__ ("lfence");
            __memcpy_sse2_unaligned((float*)&data4000_dst[0],(float*)&data4000_src[0],(std::size_t)(sizeof(float)*size_4000));
           // __asm__ __volatile__ ("lfence");
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


// No correctness checking implemented.
__attribute__((noinline))
void perf_test_sse_memcpy_unroll8x_ps_size4000(unsigned __int64 * __restrict__ ,
                                                unsigned __int64 * __restrict__ ,
                                                unsigned __int64 * __restrict__ ,
                                                const int32_t ,
                                                const int32_t,
                                                uint32_t & tid);

void perf_test_sse_memcpy_unroll8x_ps_size4000(unsigned __int64 * __restrict__ my_memcpy_s,
                                                unsigned __int64 * __restrict__ my_memcpy_e,
                                                unsigned __int64 * __restrict__ my_memcpy_d,
                                                const int32_t n_runs,
                                                const int32_t n_samples,
                                                uint32_t & tid)
{
     using namespace gms::common;
     constexpr std::size_t RDTSCP_LAT{42ull};
     constexpr std::size_t size_4000{4000ULL};
     [[maybe_unused]]
     volatile unsigned __int64 warmup_start;
     [[maybe_unused]]
     volatile unsigned __int64 warmup_end;
     uint32_t mem_start{0};
     uint32_t mem_end{0};
     float data4000_src[size_4000];
     float data4000_dst[size_4000]  = {}; 
     const float filler{3.14159265358979323846264338328F};
     sse_memset_unroll8x_ps(&data4000_src[0],filler,size_4000);
     // warmup of RDTSCP
     warmup_start = __rdtscp(&mem_start);
     warmup_end   = __rdtscp(&mem_end);
     // Warmup of __memcpy_sse2_unaligned
     sse_memcpy_unroll8x_ps(&data4000_dst[0],&data4000_src[0],size_4000);

     // Begin  looped timing execution 
     for(int32_t __i{0}; __i != n_runs; ++__i)
     {
         for(int32_t __j{0}; __j != n_samples; ++__j) 
         {   
            __asm__ __volatile__ ("lfence");
             unsigned __int64 start_curr{__rdtscp(&mem_start)};
            //__asm__ __volatile__ ("lfence");
            sse_memcpy_unroll8x_ps(&data4000_dst[0],&data4000_src[0],size_4000);
           // __asm__ __volatile__ ("lfence");
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
void test_runner_omp_sections_1st_seq();

void test_runner_omp_sections_1st_seq()
{
    constexpr int32_t n_runs{10};
    constexpr int32_t n_samples{50};
    constexpr int32_t n_total{n_runs*n_samples};
    unsigned __int64 glibc_memcpy_s[n_total] = {UINT64_MAX};
    unsigned __int64 glibc_memcpy_e[n_total] = {UINT64_MAX}; 
    unsigned __int64 glibc_memcpy_d[n_total] = {UINT64_MAX};
    unsigned __int64 my_memcpy_s[n_total]    = {UINT64_MAX};
    unsigned __int64 my_memcpy_e[n_total]    = {UINT64_MAX};
    unsigned __int64 my_memcpy_d[n_total]    = {UINT64_MAX};
    uint32_t glibc_memcpy_tid;
    uint32_t my_memcpy_tid;
    int32_t setenv_ret;
    setenv_ret = setenv("OMP_PROC_BIND","true",1);
    if(setenv_ret==-1)
    {
        printf("[**ERROR**]: -- setenv reported an error=%d\n",setenv_ret);
    }
    
#pragma omp parallel sections 
{
   #pragma omp section 
   {
       perf_test_memcpy_sse2_unaligned_size4000(glibc_memcpy_s,
                                                 glibc_memcpy_e,
                                                 glibc_memcpy_d,
                                                 n_runs,
                                                 n_samples,
                                                 glibc_memcpy_tid );
   }

   #pragma omp section 
   {
       perf_test_sse_memcpy_unroll8x_ps_size4000(my_memcpy_s,
                                                  my_memcpy_e,
                                                  my_memcpy_d,
                                                  n_runs,
                                                  n_samples,
                                                  my_memcpy_tid);
   }
}

     printf(ANSI_COLOR_GREEN "[PERF-TEST: -- size: 4000 floats]: __memcpy_sse2_unaligned: Started -- dumping-results\n" );
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
     printf(ANSI_COLOR_GREEN "[PERF-TEST: -- size: 4000 floats]: __memcpy_sse2_unaligned: Finished -- dumping-results" ANSI_RESET_ALL"\n\n");

     printf(ANSI_COLOR_RED "[PERF-TEST: -- size: 4000 floats]: sse_memcpy_unroll8x_ps: Started -- dumping-results\n" );
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
     printf(ANSI_COLOR_RED "[PERF-TEST: -- size: 4000 floats]: sse_memcpy_unroll8x_ps: Finished -- dumping-results" ANSI_RESET_ALL"\n");

}

// TID #0 handles sse_memcpy_unroll8x_ps
void test_runner_omp_sections_2nd_seq();

void test_runner_omp_sections_2nd_seq()
{
    constexpr int32_t n_runs{10};
    constexpr int32_t n_samples{50};
    constexpr int32_t n_total{n_runs*n_samples};
    unsigned __int64 glibc_memcpy_s[n_total] = {UINT64_MAX};
    unsigned __int64 glibc_memcpy_e[n_total] = {UINT64_MAX}; 
    unsigned __int64 glibc_memcpy_d[n_total] = {UINT64_MAX};
    unsigned __int64 my_memcpy_s[n_total]    = {UINT64_MAX};
    unsigned __int64 my_memcpy_e[n_total]    = {UINT64_MAX};
    unsigned __int64 my_memcpy_d[n_total]    = {UINT64_MAX};
    uint32_t glibc_memcpy_tid;
    uint32_t my_memcpy_tid;
    int32_t setenv_ret;
    setenv_ret = setenv("OMP_PROC_BIND","true",1);
    if(setenv_ret==-1)
    {
        printf("[**ERROR**]: -- setenv reported an error=%d\n",setenv_ret);
    }
    
#pragma omp parallel sections 
{
   #pragma omp section 
   {
       

       perf_test_sse_memcpy_unroll8x_ps_size4000(my_memcpy_s,
                                                  my_memcpy_e,
                                                  my_memcpy_d,
                                                  n_runs,
                                                  n_samples,
                                                  my_memcpy_tid);
   }

   #pragma omp section 
   {
       perf_test_memcpy_sse2_unaligned_size4000(glibc_memcpy_s,
                                                 glibc_memcpy_e,
                                                 glibc_memcpy_d,
                                                 n_runs,
                                                 n_samples,
                                                 glibc_memcpy_tid );
   }
}

     printf(ANSI_COLOR_GREEN "[PERF-TEST: -- size: 4000 floats]: __memcpy_sse2_unaligned: Started -- dumping-results\n" );
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
     printf(ANSI_COLOR_GREEN "[PERF-TEST: -- size: 4000 floats]: __memcpy_sse2_unaligned: Finished -- dumping-results" ANSI_RESET_ALL"\n\n");

     printf(ANSI_COLOR_RED "[PERF-TEST: -- size: 4000 floats]: sse_memcpy_unroll8x_ps: Started -- dumping-results\n" );
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
     printf(ANSI_COLOR_RED "[PERF-TEST: -- size: 4000 floats]: sse_memcpy_unroll8x_ps: Finished -- dumping-results" ANSI_RESET_ALL"\n");

}


int main()
{   
    std::clock_t seed{};
    int32_t      which{-1};
    seed = std::clock();
    auto rand{std::bind(std::uniform_int_distribution<int32_t>{0,1},
                          std::mt19937(seed))};
    which = rand();
    switch (which)
    {
        case 0 : 
        test_runner_omp_sections_1st_seq();
        break;
        case 1 : 
        test_runner_omp_sections_2nd_seq();
        break;
        default : 
        printf("[ERROR:] -- switch variable=%d\n",which);
        std::exit(EXIT_FAILURE);
    }
    return 0;
}
