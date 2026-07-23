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
#include "GMS_matrix_inv_cholesky.h"
#include "unit_test_mat_inv_cholesky_defines.h"
#include "GMS_machine_utils.h"

/*
   icpc -o perf_test_mat_inv_cholesky_16x16_instr -std=c++17 -fp-model fast=2 -fopenmp -ftz -ggdb -ipo -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_machine_utils.h unit_test_mat_inv_cholesky_defines.h GMS_matrix_inv_cholesky.h GMS_matrix_inv_cholesky.cpp perf_test_mat_inv_cholesky_16x16_instr.cpp 
   ASM: 
   icpc -S -fverbose-asm -masm=intel -fopenmp -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_machine_utils.h GMS_matrix_inv_cholesky.h GMS_matrix_inv_cholesky.cpp unit_test_mat_inv_cholesky_defines.h perf_test_mat_inv_cholesky_16x16_instr.cpp

*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

static inline
bool detect_subnormal(const float x)
{
    constexpr std::int32_t subnormal_code = 3;
    const std::int32_t subnormal_value = std::fpclassify(x);
    bool b_ret{};
    (subnormal_code==subnormal_value) ? b_ret = true : b_ret = false;
    return (b_ret);
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

#define VEC_PS_LEN 16
#define PRINT_COMPUTATION_RESULTS 0

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
void 
perf_test_matrix_inv_cholesky_16x16_rf_emt_instr(const std::int32_t,
                                                 const std::int32_t);

void 
perf_test_matrix_inv_cholesky_16x16_rf_emt_instr(const std::int32_t cpu_core,
                                                 const std::int32_t priority)
{
    using namespace gms::common;
    using namespace gms::math;
    constexpr int32_t nelems{4096};
    thread_local mat_inv_chol_tsc_instr::tsc_instrumentation_block_t<MAT_SQR_SIZE_16> tsc_instr_block;
    __ATTR_ALIGN__(64) float buf_init_re[nelems];
    __ATTR_ALIGN__(64) float buf_init_im[nelems];
    __ATTR_ALIGN__(64) __m512 matBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16];
    __ATTR_ALIGN__(64) __m512 matBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16];
    __ATTR_ALIGN__(64) __m512 matInvBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16];//flexran
    __ATTR_ALIGN__(64) __m512 matInvBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16];//flexran    
    const float * __restrict ptr_invBRe  = nullptr;
    const float * __restrict ptr_invBIm  = nullptr;
    unsigned long long seed_re{0ull};
    unsigned long long seed_im{0ull};
    [[maybe_unused]] std::int32_t printf_ret{};
    seed_re = __rdtsc();
    auto rand_re{std::bind(std::uniform_real_distribution<float>(0.001f,5.0f),std::mt19937(seed_re))};
    seed_im = __rdtsc();
    auto rand_im{std::bind(std::uniform_real_distribution<float>(0.001f,5.0f),std::mt19937(seed_im))};
    char * distro_name{NULL};
    int32_t status{9999};
    distro_name = demangle(typeid(rand_re).name(),status);
    if(distro_name != NULL && status == 0)
    {
        printf_ret = printf("[UNIT-TEST]: Instantiated distribution of type: %s\n", distro_name);
        gms::common::gms_mm_free(distro_name);
    }
    else
    {
        printf_ret = printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n", typeid(rand_re).name());
        if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
    }
    distro_name = demangle(typeid(rand_im).name(),status);
    if(distro_name != NULL && status == 0)
    {
        printf_ret = printf("[UNIT-TEST]: Instantiated distribution of type: %s\n", distro_name);
        gms::common::gms_mm_free(distro_name);
    }
    else
    {
        printf_ret = printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n", typeid(rand_im).name());
        if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
    }
    for(int32_t __i{0}; __i != nelems; ++__i)
    {
        float rd_re = rand_re();
        buf_init_re[__i] = rd_re;
        float rd_im = rand_im();
        buf_init_im[__i] = rd_im;
    }
    for(std::int32_t i = 0;i != MAT_SQR_SIZE_16; ++i)
    {   
        const std::int32_t outer_idx = i*MAT_SQR_SIZE_16;
        for(std::int32_t j = 0;j != MAT_SQR_SIZE_16; ++j)  
        {
            const std::int32_t inner_idx = (outer_idx+j)*VEC_PS_LEN;
            register __m512 init_vre = _mm512_load_ps((const float*)&buf_init_re[inner_idx]);                    
            matBRe[i][j]             = init_vre;
            //matInvBReRef[i][j]       = init_vre;
            register __m512 init_vim = _mm512_load_ps((const float*)&buf_init_im[inner_idx]);
            matBIm[i][j]             = init_vim;
            //matInvBImRef[i][j]       = init_vim;
        }
    }
    std::uint32_t rf_emt_tid{};
    std::int32_t setenv_ret{};
    std::int32_t status_prio{};
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
    status_prio = set_affinity_and_priority(cpu_core,priority);
    if(status_prio > 0)
    { 
       int32_t which = PRIO_PROCESS;
       id_t pid      = getpid();
       errno         = 0;
       int32_t nice  = getpriority(which,pid);
       printf_ret = printf("[***ERROR***]: set_affinity_and_priority status:%d,errno=%d,pid=%d,nice=%d\n",status_prio,errno,pid,nice);
    }
    
    print_thread_affinity();
    printf_ret = printf(ANSI_COLOR_GREEN "[PERF-TEST]:  matrix_inv_cholesky_16x16_tsc_instr<true,true> -- (TSC-INSTRUMENTED): STARTED!!\n");
    mat_inv_cholesky_16x16_16xf32_tsc_instr<true,true,3>(matBRe,matBIm,
                                                       matInvBRe,matInvBIm,
                                                       tsc_instr_block);
    printf_ret = printf(ANSI_COLOR_GREEN "[PERF-TEST]:  matrix_inv_cholesky_16x16_tsc_instr<true,true> -- (TSC-INSTRUMENTED): COMPLETED!!\n");
    tsc_instr_block.m_region_prefetch_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_prefetch_s,
                                                                                tsc_instr_block.m_region_prefetch_e,
                                                                                true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_prefetch"),
                                          tsc_instr_block.m_region_prefetch_s,
                                          tsc_instr_block.m_region_prefetch_e,
                                          tsc_instr_block.m_region_prefetch_d);
    tsc_instr_block.m_region_g00_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_g00_s,
                                                                           tsc_instr_block.m_region_g00_e,
                                                                           true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_g00"),
                                          tsc_instr_block.m_region_g00_s,
                                          tsc_instr_block.m_region_g00_e,
                                          tsc_instr_block.m_region_g00_d);
    tsc_instr_block.m_region_g11_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_g11_s,
                                                                           tsc_instr_block.m_region_g11_e,
                                                                           true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_g11"),
                                          tsc_instr_block.m_region_g11_s,
                                          tsc_instr_block.m_region_g11_e,
                                          tsc_instr_block.m_region_g11_d);    
    tsc_instr_block.m_region_g22_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_g22_s,
                                                                           tsc_instr_block.m_region_g22_e,
                                                                           true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_g22"),
                                          tsc_instr_block.m_region_g22_s,
                                          tsc_instr_block.m_region_g22_e,
                                          tsc_instr_block.m_region_g22_d);
    tsc_instr_block.m_region_g33_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_g33_s,
                                                                           tsc_instr_block.m_region_g33_e,
                                                                           true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_g33"),
                                          tsc_instr_block.m_region_g33_s,
                                          tsc_instr_block.m_region_g33_e,
                                          tsc_instr_block.m_region_g33_d);
    tsc_instr_block.m_region_g44_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_g44_s,
                                                                           tsc_instr_block.m_region_g44_e,
                                                                           true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_g44"),
                                          tsc_instr_block.m_region_g44_s,
                                          tsc_instr_block.m_region_g44_e,
                                          tsc_instr_block.m_region_g44_d);
    tsc_instr_block.m_region_g55_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_g55_s,
                                                                           tsc_instr_block.m_region_g55_e,
                                                                           true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_g55"),
                                          tsc_instr_block.m_region_g55_s,
                                          tsc_instr_block.m_region_g55_e,
                                          tsc_instr_block.m_region_g55_d);
    tsc_instr_block.m_region_g66_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_g66_s,
                                                                           tsc_instr_block.m_region_g66_e,
                                                                           true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_g66"),
                                          tsc_instr_block.m_region_g66_s,
                                          tsc_instr_block.m_region_g66_e,
                                          tsc_instr_block.m_region_g66_d);
    tsc_instr_block.m_region_g77_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_g77_s,
                                                                           tsc_instr_block.m_region_g77_e,
                                                                           true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_g77"),
                                          tsc_instr_block.m_region_g77_s,
                                          tsc_instr_block.m_region_g77_e,
                                          tsc_instr_block.m_region_g77_d);
    tsc_instr_block.m_region_gii_even_col8_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_gii_even_col8_s,
                                                                                     tsc_instr_block.m_region_gii_even_col8_e,
                                                                                     true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_gii_even_col8"),
                                          tsc_instr_block.m_region_gii_even_col8_s,
                                          tsc_instr_block.m_region_gii_even_col8_e,
                                          tsc_instr_block.m_region_gii_even_col8_d);
    tsc_instr_block.m_region_gii_odd_col9_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_gii_odd_col9_s,
                                                                                    tsc_instr_block.m_region_gii_odd_col9_e,
                                                                                    true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_gii_odd_col9"),
                                          tsc_instr_block.m_region_gii_odd_col9_s,
                                          tsc_instr_block.m_region_gii_odd_col9_e,
                                          tsc_instr_block.m_region_gii_odd_col9_d);
    tsc_instr_block.m_region_gii_even_col10_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_gii_even_col10_s,
                                                                                      tsc_instr_block.m_region_gii_even_col10_e,
                                                                                      true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_gii_even_col10"),
                                          tsc_instr_block.m_region_gii_even_col10_s,
                                          tsc_instr_block.m_region_gii_even_col10_e,
                                          tsc_instr_block.m_region_gii_even_col10_d);
    tsc_instr_block.m_region_gii_odd_col11_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_gii_odd_col11_s,
                                                                                     tsc_instr_block.m_region_gii_odd_col11_e,
                                                                                     true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_gii_odd_col11"),
                                          tsc_instr_block.m_region_gii_odd_col11_s,
                                          tsc_instr_block.m_region_gii_odd_col11_e,
                                          tsc_instr_block.m_region_gii_odd_col11_d);     
    tsc_instr_block.m_region_gii_even_col12_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_gii_even_col12_s,
                                                                                      tsc_instr_block.m_region_gii_even_col12_e,
                                                                                      true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_gii_even_col12"),
                                          tsc_instr_block.m_region_gii_even_col12_s,
                                          tsc_instr_block.m_region_gii_even_col12_e,
                                          tsc_instr_block.m_region_gii_even_col12_d);
    tsc_instr_block.m_region_gii_odd_col13_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_gii_odd_col13_s,
                                                                                     tsc_instr_block.m_region_gii_odd_col13_e,
                                                                                     true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_gii_odd_col13"),
                                          tsc_instr_block.m_region_gii_odd_col13_s,
                                          tsc_instr_block.m_region_gii_odd_col13_e,
                                          tsc_instr_block.m_region_gii_odd_col13_d);
    tsc_instr_block.m_region_gii_even_col14_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_gii_even_col14_s,
                                                                                     tsc_instr_block.m_region_gii_even_col14_e,
                                                                                     true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_gii_even_col14"),
                                          tsc_instr_block.m_region_gii_even_col14_s,
                                          tsc_instr_block.m_region_gii_even_col14_e,
                                          tsc_instr_block.m_region_gii_even_col14_d);
    tsc_instr_block.m_region_gii_odd_col15_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_gii_odd_col15_s,
                                                                                     tsc_instr_block.m_region_gii_odd_col15_e,
                                                                                     true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_gii_odd_col15"),
                                          tsc_instr_block.m_region_gii_odd_col15_s,
                                          tsc_instr_block.m_region_gii_odd_col15_e,
                                          tsc_instr_block.m_region_gii_odd_col15_d);
    tsc_instr_block.m_region_lii_col0_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_lii_col0_s,
                                                                                tsc_instr_block.m_region_lii_col0_e,
                                                                                true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_lii_col0"),
                                          tsc_instr_block.m_region_lii_col0_s,
                                          tsc_instr_block.m_region_lii_col0_e,
                                          tsc_instr_block.m_region_lii_col0_d);
    tsc_instr_block.m_region_lii_col1_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_lii_col1_s,
                                                                                tsc_instr_block.m_region_lii_col1_e,
                                                                                true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_lii_col1"),
                                          tsc_instr_block.m_region_lii_col1_s,
                                          tsc_instr_block.m_region_lii_col1_e,
                                          tsc_instr_block.m_region_lii_col1_d);
    tsc_instr_block.m_region_lii_col2_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_lii_col2_s,
                                                                                tsc_instr_block.m_region_lii_col2_e,
                                                                                true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_lii_col2"),
                                          tsc_instr_block.m_region_lii_col2_s,
                                          tsc_instr_block.m_region_lii_col2_e,
                                          tsc_instr_block.m_region_lii_col2_d);
    tsc_instr_block.m_region_lii_col3_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_lii_col3_s,
                                                                                tsc_instr_block.m_region_lii_col3_e,
                                                                                true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_lii_col3"),
                                          tsc_instr_block.m_region_lii_col3_s,
                                          tsc_instr_block.m_region_lii_col3_e,
                                          tsc_instr_block.m_region_lii_col3_d);
    tsc_instr_block.m_region_lii_col4_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_lii_col4_s,
                                                                                tsc_instr_block.m_region_lii_col4_e,
                                                                                true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_lii_col4"),
                                          tsc_instr_block.m_region_lii_col4_s,
                                          tsc_instr_block.m_region_lii_col4_e,
                                          tsc_instr_block.m_region_lii_col4_d);
    tsc_instr_block.m_region_lii_col5_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_lii_col5_s,
                                                                                tsc_instr_block.m_region_lii_col5_e,
                                                                                true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_lii_col5"),
                                          tsc_instr_block.m_region_lii_col5_s,
                                          tsc_instr_block.m_region_lii_col5_e,
                                          tsc_instr_block.m_region_lii_col5_d);
    tsc_instr_block.m_region_lii_col6_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_lii_col6_s,
                                                                                tsc_instr_block.m_region_lii_col6_e,
                                                                                true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_lii_col6"),
                                          tsc_instr_block.m_region_lii_col6_s,
                                          tsc_instr_block.m_region_lii_col6_e,
                                          tsc_instr_block.m_region_lii_col6_d);
    tsc_instr_block.m_region_lii_col7_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_lii_col7_s,
                                                                                tsc_instr_block.m_region_lii_col7_e,
                                                                                true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_lii_col7"),
                                          tsc_instr_block.m_region_lii_col7_s,
                                          tsc_instr_block.m_region_lii_col7_e,
                                          tsc_instr_block.m_region_lii_col7_d);
    tsc_instr_block.m_region_lii_col8_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_lii_col8_s,
                                                                                tsc_instr_block.m_region_lii_col8_e,
                                                                                true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_lii_col8"),
                                          tsc_instr_block.m_region_lii_col8_s,
                                          tsc_instr_block.m_region_lii_col8_e,
                                          tsc_instr_block.m_region_lii_col8_d);
    tsc_instr_block.m_region_lii_col9_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_lii_col9_s,
                                                                                tsc_instr_block.m_region_lii_col9_e,
                                                                                true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_lii_col9"),
                                          tsc_instr_block.m_region_lii_col9_s,
                                          tsc_instr_block.m_region_lii_col9_e,
                                          tsc_instr_block.m_region_lii_col9_d);
    tsc_instr_block.m_region_lii_col10_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_lii_col10_s,
                                                                                tsc_instr_block.m_region_lii_col10_e,
                                                                                true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_lii_col10"),
                                          tsc_instr_block.m_region_lii_col10_s,
                                          tsc_instr_block.m_region_lii_col10_e,
                                          tsc_instr_block.m_region_lii_col10_d);
    tsc_instr_block.m_region_lii_col11_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_lii_col11_s,
                                                                                tsc_instr_block.m_region_lii_col11_e,
                                                                                true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_lii_col11"),
                                          tsc_instr_block.m_region_lii_col11_s,
                                          tsc_instr_block.m_region_lii_col11_e,
                                          tsc_instr_block.m_region_lii_col11_d);
    tsc_instr_block.m_region_lii_col12_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_lii_col12_s,
                                                                                tsc_instr_block.m_region_lii_col12_e,
                                                                                true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_lii_col12"),
                                          tsc_instr_block.m_region_lii_col12_s,
                                          tsc_instr_block.m_region_lii_col12_e,
                                          tsc_instr_block.m_region_lii_col12_d);
    tsc_instr_block.m_region_lii_col13_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_lii_col13_s,
                                                                                tsc_instr_block.m_region_lii_col13_e,
                                                                                true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_lii_col13"),
                                          tsc_instr_block.m_region_lii_col13_s,
                                          tsc_instr_block.m_region_lii_col13_e,
                                          tsc_instr_block.m_region_lii_col13_d);
    tsc_instr_block.m_region_lii_col14_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_lii_col14_s,
                                                                                tsc_instr_block.m_region_lii_col14_e,
                                                                                true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_lii_col14"),
                                          tsc_instr_block.m_region_lii_col14_s,
                                          tsc_instr_block.m_region_lii_col14_e,
                                          tsc_instr_block.m_region_lii_col14_d);
    tsc_instr_block.m_region_lii_col15_d = mat_inv_chol_tsc_instr::compute_delta(tsc_instr_block.m_region_lii_col15_s,
                                                                                tsc_instr_block.m_region_lii_col15_e,
                                                                                true);
    mat_inv_chol_tsc_instr::print_results(std::string("region_lii_col15"),
                                          tsc_instr_block.m_region_lii_col15_s,
                                          tsc_instr_block.m_region_lii_col15_e,
                                          tsc_instr_block.m_region_lii_col15_d);
    for(std::int32_t i = 0;i != MAT_SQR_SIZE_16; ++i)
    {   
        const std::uint64_t loop2D_s = tsc_instr_block.m_region_loop2D_s[i];
        const std::uint64_t loop2D_e = tsc_instr_block.m_region_loop2D_e[i];
        tsc_instr_block.m_region_loop2D_d[i] = mat_inv_chol_tsc_instr::compute_delta(loop2D_s,loop2D_e,true);
        const std::uint64_t loop3D_s = tsc_instr_block.m_region_loop3D_s[i];
        const std::uint64_t loop3D_e = tsc_instr_block.m_region_loop3D_e[i];
        tsc_instr_block.m_region_loop3D_d[i] = mat_inv_chol_tsc_instr::compute_delta(loop3D_s,loop3D_e,true);
        const std::uint64_t loop2D_d = tsc_instr_block.m_region_loop2D_d[i];
        mat_inv_chol_tsc_instr::print_results(std::string("region_loop2D"),loop2D_s,loop2D_e,loop2D_d);
        const std::uint64_t loop3D_d = tsc_instr_block.m_region_loop3D_d[i];
        mat_inv_chol_tsc_instr::print_results(std::string("region_loop3D"),loop3D_s,loop3D_e,loop3D_d);
    }                                
                                                                                                                         
    printf_ret = printf(ANSI_COLOR_GREEN "[PERF-TEST]: matrix_inv_cholesky_16x16<true,true>(RF-EMT): Finished -- dumping-results" ANSI_RESET_ALL"\n\n");
#if (PRINT_COMPUTATION_RESULTS) == 1
    ptr_invBRe = reinterpret_cast<const float * __restrict>(&matInvBRe);
    ptr_invBIm = reinterpret_cast<const float * __restrict>(&matInvBIm);
    for(int32_t __i{0}; __i != nelems; ++__i)
    {  
        const float bre     = ptr_invBRe[__i];
        const float bim     = ptr_invBIm[__i];
        printf_ret =  printf(ANSI_COLOR_RED "[PERF-TEST]: (RF-EMT)bre=%.7f\n",bre);
        printNumber("(RF-EMT)bre",bre,0);
        printf_ret =  printf(ANSI_COLOR_RED "[PERF-TEST]: (RF-EMT)bim=%.7f\n",bim);
        printNumber("(RF-EMT)bim",bim,0);
        if(__builtin_expect(detect_subnormal(bre),0)) { printf_ret = printf(ANSI_COLOR_RED "[PERF-TEST:] -- (RF-EMT) ***SUBNORMAL*** in I-ch -- %.17f\n",bre);}
        if(__builtin_expect(detect_subnormal(bim),0)) { printf_ret = printf(ANSI_COLOR_RED "[PERF-TEST:] -- (RF-EMT) ***SUBNORMAL*** in Q-ch -- %.17f\n",bim);}       
    } 
#endif 
}

int main()
{
    perf_test_matrix_inv_cholesky_16x16_rf_emt_instr(5,25);
    return 0;
}
