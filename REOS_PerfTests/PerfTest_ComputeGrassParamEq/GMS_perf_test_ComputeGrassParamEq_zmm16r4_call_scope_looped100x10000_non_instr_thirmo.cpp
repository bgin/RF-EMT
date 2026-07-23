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
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>
#include <sched.h>
#include <errno.h>
#include "GMS_grass_scatterer_AVX512.h"
#include "GMS_fast_pmc_access.h"
#include "GMS_timsac_iface"
#include "GMS_avx512_warmup_loops.h"

#if !defined(MALLOC_FAILED)
#define MALLOC_FAILED                                                                                          	\
    do {									                                \
        printf(" %s -- _mm_malloc failed to allocate memory!! at line: %d\n", __PRETTY_FUNCTION__, __LINE__);   \
        exit(EXIT_FAILURE);                                                                                     \
   } while(0)
#endif

#if !defined(TIMSAC_THIRMO_STATIC_ARRAYS_DEFINE_BLOCK)
    #define TIMSAC_THIRMO_STATIC_ARRAYS_DEFINE_BLOCK(a1,a2)                     	\
       do {                                                                            	\
            __attribute__((aligned(64))) double a1[lag1+7];                          	\
            __attribute__((aligned(64))) double a2[lag1+7];                          	\
     } while(0)
#endif

#if !defined(TIMSAC_THIRMO_STATIC_ARRAYS_ZERO_INIT_LOOP)
    #define TIMSAC_THIRMO_STATIC_ARRAYS_ZERO_INIT_LOOP(a1,a2)                           	\
       do {                                                                                    	\       
            for(int32_t i = 0; i != (lag1+7); i += 16) {                                       	\
                _mm512_stream_pd(&a1[i+0], vzf64);                                           	\
		_mm512_stream_pd(&a2[i+0], vzf64);                                           	\
		_mm512_stream_pd(&a1[i+8], vzf64);                                           	\
		_mm512_stream_pd(&a2[i+8], vzf64);                                           	\
	    }                                                                                   \
       } while(0)
#endif

#if !defined(TIMSAC_THIRMO_DYNAMIC_ARRAYS_ZERO_INIT_LOOP)
    #define TIMSAC_THIRMO_DYNAMIC_ARRAYS_ZERO_INIT_LOOP(a)                                    	\
       do {                                                                                     \
            for(int32_t i = 0; i != (lag1*lag1+7); i += 16) {                                  	\
                _mm512_stream_pd(&a[i+0], vzf64);                                            	\
		_mm512_stream_pd(&a[i+8], vzf64);                                            	\
	    }                                                                                   \
       } while(0)
#endif

        // Timsac subroutine THIRMOF called to compute Time-Series 3rd moments. 

void perf_test_ComputeGrassParamEq_zmm16r4_call_scope_looped100x10000_non_instr_thirmo() {

     // By Intel Intrinsic Guide (SKX CPU)
     constexpr int32_t rdtscp_latency = 42; //stored in register or precomputed at runtime.
     constexpr std::size_t alignment = 64ULL;
     constexpr int32_t n_runs = 100;
     constexpr int32_t n_samples = 10000;
     constexpr std::size_t tot_samples = reinterpret_cast<std::size_t>(n_runs*n_samples);
     const int32_t data_len = n_runs*n_samples; // Used by F77 wrapper
     char msr_file_name[64];
     //uint64_t core_counters_start[4] = {};
     //uint64_t core_counters_end[4]   = {};
     uint64_t data[4] = {};
     uint32_t regs[4] = {0x186,0x187,0x188,0x189};
     uint64_t count;
     // Measurement arrays deltas
       // Allocation of arrays containing measurements results (deltas)
     uint64_t * __restrict PMC0  = NULL;
     uint64_t * __restrict PMC1  = NULL;
     uint64_t * __restrict PMC2  = NULL;
     uint64_t * __restrict PMC3  = NULL;
     uint64_t * __restrict FCRef = NULL; // fixed counter -- reference cycles
     uint64_t * __restrict FCAct = NULL; // fixed counter -- actual cycles
     uint64_t * __restrict FCIns = NULL; // fixed counter -- retired instructions
     uint64_t * __restrict TSC   = NULL;
    // Arrays holding converted to float measurements deltas
     double * __restrict PMC0_f64  = NULL;
     double * __restrict PMC1_f64  = NULL;
     double * __restrict PMC2_f64  = NULL;
     double * __restrict PMC3_f64  = NULL;
     double * __restrict FCRef_f64 = NULL;
     double * __restrict FCAct_f64 = NULL;
     double * __restrict FCIns_f64 = NULL;
     double * __restrict TSC_f64   = NULL;
     double   * __restrict core_util   = NULL;  // core utilization
     double   * __restrict core_avg_f  = NULL; // core average frequency
     //uint64_t fix_ins_start,fix_ins_end;
     //uint64_t gen_cyc_start,gen_cyc_end;
     //uint64_t gen_ref_start,gen_ref_end;
     //uint64_t tsc_start,tsc_end;
     //uint64_t lvl2_cyc_start,lvl2_cyc_end;
     //uint64_t core_halt_delta;
     //uint64_t volatile dummy0;
     //uint64_t volatile dummy1;
     //uint64_t volatile dummy2;
     //double utilization,avg_ghz;
     uint64_t volatile PMC0_start,PMC0_end;
     uint64_t volatile PMC1_start,PMC1_end;
     uint64_t volatile PMC2_start,PMC2_end;
     uint64_t volatile PMC3_start,PMC3_end;
     uint64_t volatile FCRef_start,FCRef_end;
     uint64_t volatile FCAct_start,FCAct_end;
     uint64_t volatile FCIns_start,FCIns_end;
     uint64_t volatile TSC_start,TSC_end;
     uint64_t ref_warmup_start,ref_warmup_end;
     uint64_t lat_cpuid_start,lat_cpuid_end,cpuid_bias;
     uint64_t lat_rdpmc_start,lat_rdpmc_end,rdpmc_bias;
     uint64_t a,d,c; // dummy args for rdpmc bias.
     uint64_t dumm0,dumm1,dumm2,dumm3; // dummies for uop-caching rdpmc.
     float nominal_ghz;
     int32_t fd;
     int32_t cpu;
     int32_t chip,core;
     int32_t counter;
     int32_t core_ctr_width;
     int32_t fix_ctr_width;
     uint32_t cyc1_lo,cyc1_hi;
     uint32_t cyc2_lo,cyc2_hi;
     uint32_t cyc3_lo,cyc3_hi;
     uint32_t cyc4_lo,cyc4_hi;
     uint32_t eax,ebx,ecx,edx; // cpuid
     // Begin an allocation of measurement arrays and measurement deltas converted to double precision.
     // Not using  GMS allocation wrappers here.
     PMC0 = reinterpret_cast<uint64_t*>(_mm_malloc(tot_samples*sizeof(uint64_t),alignment));
     if(NULL==PMC0) { MALLOC_FAILED}
     PMC1 = reinterpret_cast<uint64_t*>(_mm_malloc(tot_samples*sizeof(uint64_t),alignment));
     if(NULL==PMC1) { MALLOC_FAILED}
     PMC2 = reinterpret_cast<uint64_t*>(_mm_malloc(tot_samples*sizeof(uint64_t),alignment));
     if(NULL==PMC2) { MALLOC_FAILED}
     PMC3 = reinterpret_cast<uint64_t*>(_mm_malloc(tot_samples*sizeof(uint64_t),alignment));
     if(NULL==PMC3) { MALLOC_FAILED}
     FCRef = reinterpret_cast<uint64_t*>(_mm_malloc(tot_samples*sizeof(uint64_t),alignment));
     if(NULL==FCRef) { MALLOC_FAILED}
     FCAct = reinterpret_cast<uint64_t*>(_mm_malloc(tot_samples*sizeof(uint64_t),alignment));
     if(NULL==FCAct) { MALLOC_FAILED}
     FCIns = reinterpret_cast<uint64_t*>(_mm_malloc(tot_samples*sizeof(uint64_t),alignment));
     if(NULL==FCIns) { MALLOC_FAILED}
     TSC = reinterpret_cast<uint64_t*>(_mm_malloc(tot_samples*sizeof(uint64_t),alignment));
     if(NULL==TSC)  { MALLOC_FAILED}
     // Converted to double precision deltas
     PMC0_f64 = reinterpret_cast<double*>(_mm_malloc(tot_samples*sizeof(double),alignment));
     if(NULL==PMC0_f64) { MALLOC_FAILED}
     PMC1_f64 = reinterpret_cast<double*>(_mm_malloc(tot_samples*sizeof(double),alignment));
     if(NULL==PMC1_f64) { MALLOC_FAILED}
     PMC2_f64 = reinterpret_cast<double*>(_mm_malloc(tot_samples*sizeof(double),alignment));
     if(NULL==PMC2_f64) { MALLOC_FAILED}
     PMC3_f64 = reinterpret_cast<double*>(_mm_malloc(tot_samples*sizeof(double),alignment));
     if(NULL==PMC3_f64) { MALLOC_FAILED}
     FCRef_f64 = reinterpret_cast<double*>(_mm_malloc(tot_samples*sizeof(double),alignment));
     if(NULL==FCRef_f64) { MALLOC_FAILED}
     FCAct_f64 = reinterpret_cast<double*>(_mm_malloc(tot_samples*sizeof(double),alignment));
     if(NULL==FCAct_f64) { MALLOC_FAILED}
     FCIns_f64 = reinterpret_cast<double*>(_mm_malloc(tot_samples*sizeof(double),alignment));
     if(NULL==FCIns_f64) { MALLOC_FAILED}
     TSC_f64 = reinterpret_cast<double*>(_mm_malloc(tot_samples*sizeof(double),alignment));
     if(NULL==TSC_f64)  { MALLOC_FAILED}
     core_util = reinterpret_cast<double*>(_mm_malloc(tot_samples*sizeof(double),alignment));
     if(NULL==core_util) { MALLOC_FAILED}
     core_avg_f = reinterpret_cast<double*>(_mm_malloc(tot_samples*sizeof(double),alignment));
     if(NULL==core_avg_f) { MALLOC_FAILED}
     //Initialize the results arrays
     const __m512i vzi64 = _mm512_setzero_si512();
     const __m512d vzf64 = _mm512_setzero_pd();
     _mm512_stream_si512(&PMC0[0],vzi64);
     _mm512_stream_si512(&PMC1[0],vzi64);
     _mm512_stream_si512(&PMC2[0],vzi64);
     _mm512_stream_si512(&PMC3[0],vzi64);
     _mm512_stream_si512(&FCRef[0],vzi64);
     _mm512_stream_si512(&FCAct[0],vzi64);
     _mm512_stream_si512(&FCIns[0],vzi64);
     _mm512_store_epi64(&TSC[0],vzi64);
     _mm512_stream_pd(&core_util[0],vzf64);
     _mm512_stream_pd(&core_avg_f[0],vzf64);
     if(int32_t i = 7; i != tot_samples-3; i += 32) {
     
#include "call_scope_looped_tests_zero_init_arrays_loop_body.c"

     }
     c = (1UL<<30)+1; // Fixed counter.
     // Measurement Bias calculation (ahead of any test execution phase)
       // First a dry run.
     __asm__ __volatile__ (
                 "CPUID\n\t"
                 "RDTSC\n\t"
                 "mov %%edx, %0\n\t"
                 "mov %%eax, %1\n\t" : "=r" (cyc1_hi),"=r" (cyc1_lo) :: "%rax" "%rbx" "%rcx" "%rdx" 
                 );
     __asm__ __volatile__ ("cpuid" : \
	              "=a" (eax), "=b" (ebx), "=c" (ecx), "=d" (edx) );
     __asm__ __volatile__ ("rdpmc" : "=a" (a), "=d" (d) : "c" (c) );
     __asm__ __volatile__ (
                 "CPUID\n\t"
                 "RDTSC\n\t"
                 "mov %%edx, %0\n\t"
                 "mov %%eax, %1\n\t" : "=r" (cyc2_hi),"=r" (cyc2_lo) :: "%rax" "%rbx" "%rcx" "%rdx" 
                 );
    // Second run is less prone to random noise (IF and ID stages)
    // THis is still prone to random variation unfortunately
    // more accurate analysis may introduce unwanted cache polution
    // and prolonged execution serialization.
    
     __asm__ __volatile__ (
                 "CPUID\n\t"
                 "RDTSC\n\t"
                 "mov %%edx, %0\n\t"
                 "mov %%eax, %1\n\t" : "=r" (cyc1_hi),"=r" (cyc1_lo) :: "%rax" "%rbx" "%rcx" "%rdx" 
                 );
     __asm__ __volatile__ ("cpuid" : \
	              "=a" (eax), "=b" (ebx), "=c" (ecx), "=d" (edx) );
     __asm__ __volatile__ (
                 "CPUID\n\t"
                 "RDTSC\n\t"
                 "mov %%edx, %0\n\t"
                 "mov %%eax, %1\n\t" : "=r" (cyc2_hi),"=r" (cyc2_lo) :: "%rax" "%rbx" "%rcx" "%rdx" 
                 );
      // First a "dry run"
     __asm__ __volatile__ (
                 "CPUID\n\t"
                 "RDTSC\n\t"
                 "mov %%edx, %0\n\t"
                 "mov %%eax, %1\n\t" : "=r" (cyc3_hi),"=r" (cyc3_lo) :: "%rax" "%rbx" "%rcx" "%rdx" 
                 );
     __asm__ __volatile__ ("rdpmc" : "=a" (a), "=d" (d) : "c" (c) );
     __asm__ __volatile__ (
                 "CPUID\n\t"
                 "RDTSC\n\t"
                 "mov %%edx, %0\n\t"
                 "mov %%eax, %1\n\t" : "=r" (cyc4_hi),"=r" (cyc4_lo) :: "%rax" "%rbx" "%rcx" "%rdx" 
                 );
      // Second run is less prone to random noise (IF and ID stages)
    __asm__ __volatile__ (
                 "CPUID\n\t"
                 "RDTSC\n\t"
                 "mov %%edx, %0\n\t"
                 "mov %%eax, %1\n\t" : "=r" (cyc3_hi),"=r" (cyc3_lo) :: "%rax" "%rbx" "%rcx" "%rdx" 
                 );
     __asm__ __volatile__ ("rdpmc" : "=a" (a), "=d" (d) : "c" (c) );
     __asm__ __volatile__ (
                 "CPUID\n\t"
                 "RDTSC\n\t"
                 "mov %%edx, %0\n\t"
                 "mov %%eax, %1\n\t" : "=r" (cyc4_hi),"=r" (cyc4_lo) :: "%rax" "%rbx" "%rcx" "%rdx" 
                 );
      lat_cpuid_start =   (((uint64_t)cyc1_high << 32) | cyc1_low );
      lat_cpuid_end   =   (((uint64_t)cyc2_high << 32) | cyc2_low );
      cpuid_bias = lat_cpuid_end-lat_cpuid_start;
      lat_rdpmc_start =   (((uint64_t)cyc3_high << 32) | cyc3_low );
      lat_rdpmc_end   =   (((uint64_t)cyc4_high << 32) | cyc4_low );
      rdpmc_bias = lat_rdpmc_end-lat_rdpmc_start;
      nominal_ghz = get_TSC_frequency() / 1.0e9;
      printf("Nominal Frequency (GHz): %f\n",nominal_ghz);
      core_ctr_width = get_core_counter_width();
      fixed_ctr_width = get_core_counter_width();
      if (core_ctr_width != fixed_ctr_width) {
		printf("Warning -- programmable counter width does not match fixed-function counter width -- this code may have errors!\n");
      }
      cpu = 3;
      cpu_set_t cpu_set;
      CPU_ZERO(&cpu_set);
      CPU_SET(cpu,&cpu_set);
      if(sched_setaffinity(0,sizeof(cpu_set), &cpu_set) < 0) {
         printf("***FATAL*** -- Failed to set CPU affinity\n");
         std::exit(EXIT_FAILURE);
      }
      printf("Affinity set to cpu: %d\n",cpu);
      struct sched_param sp;
      memset(&sp,0,sizeof(sp));
      sp.sched_priority = 99;
      if((sched_setscheduler(0,SCHED_FIFO,&sp)) == -1) {
          printf("Failed to set thread priority to 99!!, errno=%d\n",errno);
      }
      count = full_rdtscp(&chip,&core);
     
     // Object instantiation parameters.
     const int32_t nsteps  = 128;
     const int32_t ordinal = 1;
     const int32_t param_npts = 360;
     const float lat = 43.040833; // degree N
     const float lon = 77.241833; // degree W
     const float elev = 124.4; // meter above sea level
     const std::complex<float> ceps{0.1f,0.2f};
     const bool lockm = true;
     // Check nominal frequency.
     //nominal_ghz = get_TSC_frequency()/1.0e9;
     //printf("Nominal frequency (Ghz): %f\n",nominal_ghz);
     //core_ctr_width = get_core_counter_width();
     //fix_ctr_width  = get_fix_counter_width();
     //printf("Programmable counters width: %d bits\n",core_ctr_width);
     //printf("Fixed counters width:        %d bits\n",fix_ctr_width);
     //if(core_ctr_width != fix_ctr_width) {
     //   printf("Warning -- programmable counters width does not match the fixed counters width!!\n");
     //}
       // Set core affinity
     //cpu = 3;
     //cpu_set_t cpu_set;
     //CPU_ZERO(&cpu_set);
     //CPU_SET(cpu,&cpu_set);
     //if(sched_setaffinity(0,sizeof(cpu_set),&cpu_set) < 0) {
     //   printf("***FATAL*** -- Failed to set core affinity: errno=%d\n",errno);
     //	exit(EXIT_FAILURE);
     //}
     //count = full_rdtscp(&chip,&core);
     
     GrassScattererAVX512 ScattererZMM16r4_1 = GrassScattererAVX512(nsteps,
                                                             ordinal,
							     param_npts,
							     lat,
							     lon,
							     elev,
							     ceps,
							     mlock);
     // Forcing core to fetch and decode ahead of first usage a CPUID and RDPMC
     // instructions.
     // Warmup-loop ahead
     fix_ins_start = rdpmc_instructions();
     tsc_start = rdtscp();
     ref_warmup_start = rdpmc_reference_cycles();
     avx512_warmup_loop3_ps(); // An attempt to trigger L2 slow stage ramp-up
     tsc_end = rdtscp();
     fix_ins_end = rdpmc_instructions();
     ref_warmup_end = rdpmc_reference_cycles();
      // Enable hot L1I$, IF and ID stages ahead of first usage.
     dumm0 = rdpmc(0);
     dumm1 = rdpmc(1);
     dumm2 = rdpmc(2);
     dumm3 = rdpmc(3);
     for(int32_t i = 0; i != n_runs; ++i) {
         for(int32_t j = 0; j != n_samples; ++j) {
               __asm__ __volatile__ ("lfence"); // lfence overhead
               PMC0_start = rdpmc(0);
               PMC1_start = rdpmc(1);
               PMC2_start = rdpmc(2);
               PMC3_start = rdpmc(3);
               FCRef_start = rdpmc_reference_cycles();
               FCAct_start = rdpmc_actual_cycles();
               FCIns_start = rdpmc_instructions();
               TSC_start   = rdtscp();
	      // Begin measurement double loop.
              // No error checking is inserted in this loop, it was
              // done in order to minimize the pollution of unwanted instruction stream.
              // Call ComputeGrassParamEq_zmm16r4
              // Ideally the overhead of call instruction and whole chain
              // of L1I$ miss and ITLB miss cost should be computed
              // and subtracted.
              // THe results here are plagued by noise.
              // For real world measurement it is a needed scenario
              // For idealized measurement the overhead should be measured and subtracted.
	       __asm__ __volatile__("lfence");
               ScattererZMM16r4_1.ComputeGrassParamEq_zmm16r4();
	       __asm__ __volatile__("lfence");
	       TSC_end  = rdtscp();
               FCIns_end = rdpmc_instructions();
               FCAct_end = rdpmc_actual_cycles();
               FCRef_end = rdpmc_reference_cycles();
               PMC3_end  = rdpmc(3);
               PMC2_end  = rdpmc(2);
               PMC1_end  = rdpmc(1);
               PMC0_end  = rdpmc(0);
               __asm__ __volatile__ ("lfence");
               TSC[i*n_samples+j]     = corrected_pmc_delta(TSC_end,TSC_start,core_ctr_width)-cpuid_bias;
               FCIns[i*n_samples+j]   = corrected_pmc_delta(FCIns_end,FCIns_start,core_ctr_width);
               FCAct[i*n_samples+j]   = corrected_pmc_delta(FCAct_end,FCAct_start,core_ctr_width);
               FCRef[i*n_samples+j]   = corrected_pmc_delta(FCRef_end,FCRef_start,core_ctr_width);
               PMC3[i*n_samples+j]    = corrected_pmc_delta(PMC0_end,PMC0_start,core_ctr_width);
               PMC2[i*n_samples+j]    = corrected_pmc_delta(PMC1_end,PMC1_start,core_ctr_width);
               PMC1[i*n_samples+j]    = corrected_pmc_delta(PMC2_end,PMC2_start,core_ctr_width);
               PMC0[i*n_samples+j]    = corrected_pmc_delta(PMC3_end,PMC3_start,core_ctr_width);
	 }
     }

     // Compute core utilization and core average frequency over test period span
     double res_copy  = 0.0;
     for(int32_t i = 0; i != tot_samples; ++i) {
         res_copy = (double)(FCRef[i]/TSC[i]);
	 core_util[i] = res_copy;
	 double tmp = res_copy*nominal_ghz;
	 core_avg_f[i] = tmp;
     }

     for(int32_t i = 0; i != tot_samples-3; i += 32) {
     
#include "call_scope_looped_tests_zero_init_deltas_loop_body.c"

     }

     //Convert results delta arrays to double precision.
     _mm512_store_pd(&PMC0_f64[0], _mm512_castsi512_pd(_mm512_load_epi64(&PMC0[0])));
      _mm512_store_pd(&PMC1_f64[0], _mm512_castsi512_pd(&PMC1[0]));
      _mm512_store_pd(&PMC2_f64[0], _mm512_castsi512_pd(&PMC2[0]));
      _mm512_store_pd(&PMC3_f64[0], _mm512_castsi512_pd(&PMC3[0]));
      _mm512_store_pd(&FCRef_f64[0],_mm512_castsi512_pd(&FCRef[0]));
      _mm512_store_pd(&FCAct_f64[0],_mm512_castsi512_pd(&FCAct[0]));
      _mm512_store_pd(&FCIns_f64[0],_mm512_castsi512_pd(&FCIns[0]));
      _mm512_store_pd(&TSC_f64[0],  _mm512_castsi512_pd(&TSC[0]));
      for(int32_t i = 7; i != tot_samples-3; i += 32) {
      
#include "call_scope_looped_tests_cast_to_f64_loop_body.c"

      }

       // Deallocating deltas array (uint64_t) type
     _mm_free(TSC);  _mm_free(FCIns); _mm_free(FCAct); _mm_free(FCRef);
     _mm_free(PMC3); _mm_free(PMC2); _mm_free(PMC1); _mm_free(PMC0);

     // Begin series of allocation->computation->deallocation stages
     // Prepare to call Timsac "autocor" 11 times for each parameter
     // time series.
     // 3rd order moments of PMC0 data
     const char * fname1 = "PMC0_3rd-order-moments.csv";
     FILE * fp1;
     TIMSAC_THIRMO_STATIC_ARRAYS_DEFINE_BLOCK((PMC0_acor),(PMC0_acov))
     double * __restrict PMC0_mnt  = NULL; 
     double PMC0_xmean = 0.0;
     constexpr int32_t lag1 = 2001;
     PMC0_mnt = reinterpret_cast<double*>(_mm_malloc((std::size_t)(lag1*lag1+7)*sizeof(double),alignment));
     if(NULL==PMC0_mnt) { MALLOC_FAILED};
     // Initialize the result arrays
     TIMSAC_THIRMO_STATIC_ARRAYS_ZERO_INIT_LOOP((PMC0_acor),(PMC0_acov))

     TIMSAC_THIRMO_DYNAMIC_ARRAYS_ZERO_INIT_LOOP((PMC0_mnt))
     // Call F77 C-interface to THIRMOF subroutine
     THIRMOF(&data_len,&lag1,&PMC0_f64[0],&PMC0_xmean,&PMC0_acov[0],&PMC0_acor[0],&PMC0_mnt[0]);
     if(fopen(&fp1,fname1,"a+") != 0) {
        printf("File open error: %s\n",fname1);
	std::exit(EXIT_FAILURE);
     }
     fprintf(fp1,"PMC1_xmean=%.16f\n",PMC0_xmean);
     for(int32_t i = 0; i != data_len; ++i) { fprintf(fp1,"%.16f\n",PMC0_f64[i]);}
     for(int32_t i = 0; i != lag1; ++i) { fprintf(fp1,"%.16f %.16f\n",PMC0_acov[i],PMC0_acor[i]);}
     for(int32_t i = 0; i != (lag1*lag1); ++i) {fprintf(fp1, "%.16f\n",PMC0_mnt[i] );}
     fclose(fp1);
     const char * fname2 = "PMC1_3rd-order-moments.csv";
     FILE * fp2;
      TIMSAC_THIRMO_STATIC_ARRAYS_DEFINE_BLOCK((PMC1_acor),(PMC1_acov))
     double * __restrict PMC1_mnt  = NULL;
     double PMC1_xmean = 0.0;
     PMC1_mnt = reinterpret_cast<double*>(_mm_malloc((std::size_t)(lag1*lag1+7)*sizeof(double),alignment));
     if(NULL==PMC1_mnt) { MALLOC_FAILED};
     // Initialize the result arrays
      TIMSAC_THIRMO_STATIC_ARRAYS_ZERO_INIT_LOOP((PMC1_acor),(PMC1_acov))

     TIMSAC_THIRMO_DYNAMIC_ARRAYS_ZERO_INIT_LOOP((PMC1_mnt))
     // Call F77 C-interface to THIRMOF subroutine
     THIRMOF(&data_len,&lag1,&PMC1_f64[0],&PMC1_xmean,&PMC1_acov[0],&PMC1_acor[0],&PMC1_mnt[0]);
     if(fopen(&fp2,fname2,"a+") != 0) {
        printf("File open error: %s\n",fname2);
	std::exit(EXIT_FAILURE);
     }
     fprintf(fp2,"PMC1_xmean=%.16f\n",PMC1_xmean);
     for(int32_t i = 0; i != data_len; ++i) { fprintf(fp2,"%.16f\n",PMC1_f64[i]);}
     for(int32_t i = 0; i != lag1; ++i) { fprintf(fp2,"%.16f %.16f\n",PMC1_acov[i],PMC1_acor[i]);}
     for(int32_t i = 0; i != (lag1*lag1); ++i) {fprintf(fp2, "%.16f\n",PMC1_mnt[i] );}
     fclose(fp2);
     const char * fname3 = "PMC2_3rd_order-moments.csv";
     FILE * fp3;
      TIMSAC_THIRMO_STATIC_ARRAYS_DEFINE_BLOCK((PMC2_acor),(PMC2_acov))
     double * __restrict PMC2_mnt  = NULL; // 3rd moment
     double PMC2_xmean = 0.0;
     PMC2_mnt = reinterpret_cast<double*>(_mm_malloc((std::size_t)(lag1*lag1+7)*sizeof(double),alignment));
     if(NULL==PMC2_mnt) { MALLOC_FAILED};
     // Initialize the result arrays

      TIMSAC_THIRMO_STATIC_ARRAYS_ZERO_INIT_LOOP((PMC2_acor),(PMC2_acov))

     TIMSAC_THIRMO_DYNAMIC_ARRAYS_ZERO_INIT_LOOP((PMC2_mnt))
     
     // Call F77 C-interface to THIRMOF subroutine
     THIRMOF(&data_len,&lag1,&PMC2_f64[0],&PMC2_xmean,&PMC2_acov[0],&PMC2_acor[0],&PMC2_mnt[0]);
     if(fopen(&fp3,fname3,"a+") != 0) {
        printf("File open error: %s\n",fname3);
	std::exit(EXIT_FAILURE);
     }
     fprintf(fp3,"PMC2_xmean=%.16f\n",PMC2_xmean);
     for(int32_t i = 0; i != data_len; ++i) { fprintf(fp3,"%.16f\n",PMC2_f64[i]);}
     for(int32_t i = 0; i != lag1; ++i) { fprintf(fp3,"%.16f %.16f\n",PMC2_acov[i],PMC2_acor[i]);}
     for(int32_t i = 0; i != (lag1*lag1); ++i) {fprintf(fp3, "%.16f\n",PMC3_mnt[i] );}
     fclose(fp3);
     const char * fname4 = "PMC3_3rd-order-moments.csv";
     FILE * fp4;
     TIMSAC_THIRMO_STATIC_ARRAYS_DEFINE_BLOCK((PMC3_acor),(PMC3_acov))
     double * __restrict PMC3_mnt  = NULL;
     double PMC3_xmean = 0.0;
     PMC3_mnt = reinterpret_cast<double*>(_mm_malloc((std::size_t)(lag1*lag1+7)*sizeof(double),alignment));
     if(NULL==PMC3_mnt) { MALLOC_FAILED};

     TIMSAC_THIRMO_STATIC_ARRAYS_ZERO_INIT_LOOP((PMC3_acor),(PMC3_acov))

     TIMSAC_THIRMO_DYNAMIC_ARRAYS_ZERO_INIT_LOOP((PMC3_mnt))
     
     // Call F77 C-interface to THIRMOF subroutine
     THIRMOF(&data_len,&lag1,&PMC3_f64[0],&PMC3_xmean,&PMC3_acov[0],&PMC3_acor[0],&PMC3_mnt[0]);
     if(fopen(&fp4,fname4,"a+") != 0) {
        printf("File open error: %s\n",fname4);
	std::exit(EXIT_FAILURE);
     }
     fprintf(fp4,"PMC3_xmean=%.16f\n",PMC3_xmean);
     for(int32_t i = 0; i != data_len; ++i) { fprintf(fp4,"%.16f\n",PMC3_f64[i]);}
     for(int32_t i = 0; i != lag1; ++i) { fprintf(fp4,"%.16f %.16f\n",PMC3_acov[i],PMC3_acor[i]);}
     for(int32_t i = 0; i != (lag1*lag1); ++i) {fprintf(fp4, "%.16f\n",PMC4_mnt[i] );}
     fclose(fp4);
     const char * fname5 = "FCRef_3rd-order-moments.csv";
     FILE * fp5;
     TIMSAC_THIRMO_STATIC_ARRAYS_DEFINE_BLOCK((FCRef_acor),(FCRef_acov))
     double * __restrict FCRef_mnt  = NULL; // 
     double FCRef_xmean = 0.0;
     FCRef_mnt = reinterpret_cast<double*>(_mm_malloc((std::size_t)(lag1*lag1+7)*sizeof(double),alignment));
     if(NULL==FCRef_mnt) { MALLOC_FAILED};
     // Initialize the result arrays

      TIMSAC_THIRMO_STATIC_ARRAYS_ZERO_INIT_LOOP((FCRef_acor),(FCRef_acov))

     TIMSAC_THIRMO_DYNAMIC_ARRAYS_ZERO_INIT_LOOP((FCRef_mnt))
     
     // Call F77 C-interface to THIRMOF subroutine
     THIRMOF(&data_len,&lag1,&FCRef_f64[0],&FCRef_xmean,&FCRef_acov[0],&FCRef_acor[0],&FCRef_mnt[0]);
     if(fopen(&fp5,fname5,"a+") != 0) {
        printf("File open error: %s\n",fname5);
	std::exit(EXIT_FAILURE);
     }
     fprintf(fp5,"FCRef_xmean=%.16f\n",FCRef_xmean);
     for(int32_t i = 0; i != data_len; ++i) { fprintf(fp5,"%.16f\n",FCRef_f64[i]);}
     for(int32_t i = 0; i != lag1; ++i) { fprintf(fp5,"%.16f %.16f\n",FCRef_acov[i],FCRef_acor[i]);}
     for(int32_t i = 0; i != (lag1*lag1); ++i) {fprintf(fp5, "%.16f\n",FCRef_mnt[i] );}
     fclose(fp5);
     const char * fname6 = "FCAct_3rd-order-moments.csv";
     FILE * fp6;
      TIMSAC_THIRMO_STATIC_ARRAYS_DEFINE_BLOCK((FCAct_acor),(FCAct_acov))
     double * __restrict FCAct_mnt  = NULL; //
     double FCAct_xmean = 0.0;
     FCAct_mnt = reinterpret_cast<double*>(_mm_malloc((std::size_t)(lag1*lag1+7)*sizeof(double),alignment));
     if(NULL==FCAct_mnt) { MALLOC_FAILED};
     // Initialize the result arrays

     TIMSAC_THIRMO_STATIC_ARRAYS_ZERO_INIT_LOOP((FCAct_acor),(FCAct_acov))

     TIMSAC_THIRMO_DYNAMIC_ARRAYS_ZERO_INIT_LOOP((FCAct_mnt))
     
     // Call F77 C-interface to THIRMOF subroutine
     THIRMOF(&data_len,&lag1,&FCAct_f64[0],&FCAct_xmean,&FCAct_acov[0],&FCAct_acor[0],&FCAct_mnt[0]);
     if(fopen(&fp6,fname6,"a+") != 0) {
        printf("File open error: %s\n",fname6);
	std::exit(EXIT_FAILURE);
     }
     fprintf(fp6,"FCAct_xmean=%.16f\n",FCAct_xmean);
     for(int32_t i = 0; i != data_len; ++i) { fprintf(fp6,"%.16f\n",FCAct_f64[i]);}
     for(int32_t i = 0; i != lag1; ++i) { fprintf(fp6,"%.16f %.16f\n",FCAct_acov[i],FCAct_acor[i]);}
     for(int32_t i = 0; i != (lag1*lag1); ++i) {fprintf(fp6, "%.16f\n",FCAct_mnt[i] );}
     fclose(fp6);
     const char * fname7 = "FCIns_3rd-order-moments.csv";
     FILE * fp7;
      TIMSAC_THIRMO_STATIC_ARRAYS_DEFINE_BLOCK((FCIns_acor),(FCIns_acov))
     double * __restrict FCIns_mnt  = NULL; // 
     double FCIns_xmean = 0.0;
     FCIns_mnt = reinterpret_cast<double*>(_mm_malloc((std::size_t)(lag1*lag1+7)*sizeof(double),alignment));
     if(NULL==FCIns_mnt) { MALLOC_FAILED};
     // Initialize the result arrays

       TIMSAC_THIRMO_STATIC_ARRAYS_ZERO_INIT_LOOP((FCIns_acor),(FCIns_acov))

     TIMSAC_THIRMO_DYNAMIC_ARRAYS_ZERO_INIT_LOOP((FCIns_mnt))
     
     // Call F77 C-interface to THIRMOF subroutine
     THIRMOF(&data_len,&lag1,&FCIns_f64[0],&FCIns_xmean,&FCIns_acov[0],&FCIns_acor[0],&FCIns_mnt[0]);
     if(fopen(&fp7,fname7,"a+") != 0) {
        printf("File open error: %s\n",fname7);
	std::exit(EXIT_FAILURE);
     }
     fprintf(fp7,"FCIns_xmean=%.16f\n",FCIns_xmean);
     for(int32_t i = 0; i != data_len; ++i) { fprintf(fp7,"%.16f\n",FCIns_f64[i]);}
     for(int32_t i = 0; i != lag1; ++i) { fprintf(fp7,"%.16f %.16f\n",FCIns_acov[i],FCIns_acor[i]);}
     for(int32_t i = 0; i != (lag1*lag1); ++i) {fprintf(fp7, "%.16f\n",FCIns_mnt[i] );}
     fclose(fp7);
     const char * fname8 = "TSC_3rd-order-moments.csv";
     FILE * fp8;
     TIMSAC_THIRMO_STATIC_ARRAYS_DEFINE_BLOCK((TSC_acor),(TSC_acov))
     double * __restrict TSC_mnt  = NULL; //
     double TSC_xmean = 0.0;
     TSC_mnt = reinterpret_cast<double*>(_mm_malloc((std::size_t)(lag1*lag1+7)*sizeof(double),alignment));
     if(NULL==TSC_mnt) { MALLOC_FAILED};
     // Initialize the result arrays

       TIMSAC_THIRMO_STATIC_ARRAYS_ZERO_INIT_LOOP((TSC_acor),(TSC_acov))

     TIMSAC_THIRMO_DYNAMIC_ARRAYS_ZERO_INIT_LOOP((TSC_mnt))
     
     // Call F77 C-interface to THIRMOF subroutine
     THIRMOF(&data_len,&lag1,&TSC_f64[0],&TSC_xmean,&TSC_acov[0],&TSC_acor[0],&TSC_mnt[0]);
     if(fopen(&fp8,fname8,"a+") != 0) {
        printf("File open error: %s\n",fname8);
	std::exit(EXIT_FAILURE);
     }
     fprintf(fp8,"TSC_xmean=%.16f\n",TSC_xmean);
     for(int32_t i = 0; i != data_len; ++i) { fprintf(fp8,"%.16f\n",TSC_f64[i]);}
     for(int32_t i = 0; i != lag1; ++i) { fprintf(fp8,"%.16f %.16f\n",TSC_acov[i],TSC_acor[i]);}
     for(int32_t i = 0; i != (lag1*lag1); ++i) {fprintf(fp8, "%.16f\n",TSC_mnt[i] );}
     fclose(fp8);
     const char * fname9 = "Core_util_3rd-order-moments.csv";
     FILE * fp9;
      TIMSAC_THIRMO_STATIC_ARRAYS_DEFINE_BLOCK((Core_util_acor),(Core_util_acov))
     double * __restrict Core_util_mnt  = NULL; // 
     double Core_util_xmean = 0.0;
     Core_util_mnt = reinterpret_cast<double*>(_mm_malloc((std::size_t)(lag1*lag1+7)*sizeof(double),alignment));
     if(NULL==Core_util_mnt) { MALLOC_FAILED};
     // Initialize the result arrays

      TIMSAC_THIRMO_STATIC_ARRAYS_ZERO_INIT_LOOP((Core_util_acor),(Core_util_acov))

     TIMSAC_THIRMO_DYNAMIC_ARRAYS_ZERO_INIT_LOOP((Core_util_mnt))
     
     // Call F77 C-interface to THIRMOF subroutine
     THIRMOF(&data_len,&lag1,&core_util[0],&Core_util_xmean,&Core_util_acov[0],&Core_util_acor[0],&Core_util_mnt[0]);
     if(fopen(&fp9,fname9,"a+") != 0) {
        printf("File open error: %s\n",fname9);
	std::exit(EXIT_FAILURE);
     }
     fprintf(fp9,"Core_util_xmean=%.16f\n",Core_util_xmean);
     for(int32_t i = 0; i != data_len; ++i) { fprintf(fp9,"%.16f\n",Core_util_f64[i]);}
     for(int32_t i = 0; i != lag1; ++i) { fprintf(fp9,"%.16f %.16f\n",Core_util_acov[i],Core_util_acor[i]);}
     for(int32_t i = 0; i != (lag1*lag1); ++i) {fprintf(fp9, "%.16f\n",Core_util_mnt[i] );}
     fclose(fp9);
     const char * fname10 = "Core_avg_f_3rd-order-moments.csv";
     FILE * fp10;
      TIMSAC_THIRMO_STATIC_ARRAYS_DEFINE_BLOCK((Core_avg_f_acor),(Core_avg_f_acov))
     double * __restrict Core_avg_f_mnt  = NULL; // 
     double Core_avg_f_xmean = 0.0;
     Core_avg_f_mnt = reinterpret_cast<double*>(_mm_malloc((std::size_t)(lag1*lag1+7)*sizeof(double),alignment));
     if(NULL==Core_avg_f_mnt) { MALLOC_FAILED};
     // Initialize the result arrays

       TIMSAC_THIRMO_STATIC_ARRAYS_ZERO_INIT_LOOP((Core_avg_f_acor),(Core_avg_f_acov))

     TIMSAC_THIRMO_DYNAMIC_ARRAYS_ZERO_INIT_LOOP((Core_avg_f_mnt))
     
     // Call F77 C-interface to THIRMOF subroutine
     THIRMOF(&data_len,&lag1,&core_avg_f[0],&Core_avg_f_xmean,&Core_avg_f_acov[0],&Core_avg_f_acor[0],&Core_avg_f_mnt[0]);
     if(fopen(&fp10,fname10,"a+") != 0) {
        printf("File open error: %s\n",fname10);
	std::exit(EXIT_FAILURE);
     }
     fprintf(fp10,"Core_avg_f_xmean=%.16f\n",Core_avg_f_xmean);
     for(int32_t i = 0; i != data_len; ++i) { fprintf(fp10,"%.16f\n",Core_avg_f_f64[i]);}
     for(int32_t i = 0; i != lag1; ++i) { fprintf(fp10,"%.16f %.16f\n",Core_avg_f_acov[i],Core_avg_f_acor[i]);}
     for(int32_t i = 0; i != (lag1*lag1); ++i) {fprintf(fp10, "%.16f\n",Core_avg_f_mnt[i] );}
     fclose(fp10);
     // Core utilization
     //utilization = (double)(corrected_pmc_delta(gen_ref_end,gen_ref_start,core_ctr_width)) /
     //                      (tsc_end-tsc_start-rdtscp_latency);
     //avg_ghz     = (double)(corrected_pmc_delta(gen_ref_end,gen_ref_start,core_ctr_width)) /
     //                      (tsc_end-tsc_start-rdtscp_latency)*nominal_ghz;
     sprintf(msr_file,"/dev/cpu/%d/msr",cpu);
     fd = open(msr_file,O_RDONLY);
     if(fd >= 0) {
        for(int32_t i = 0; i != 4; ++i) {
            if(pread(fd,&data[i],sizeof(data[i]),regs[i]) != sizeof(data)) {
               printf("***WARNING*** -- pread: CPU %d cannot read MSR: 0x%08\n",cpu,regs[i]);
	    }
	}
     }
    
     // Print test results and some statistics
     printf("The core count:      %d\n", get_core_number());
     printf("The socket number:   %d\n", get_socket_number());
     printf("CPU socket executing now: %d,  CPU core executing now: %d\n",chip,core);
     //printf("Core utilization: %f,  Core average frequency: %f\n",utilization,avg_ghz);
     //printf("Unhalt-to-Halt discrepancy: %lu\n",core_halt_delta);
     //printf("L2 cycles: %lu\n", lvl2_cyc_end-lvl_cyc_start);
     printf("Warmup-loop statistics: TSC=%lu,Ins=%lu,Ref_cyc=%lu\n",tsc_end-tsc_start,
	    fix_ins_end-fix_ins_start,ref_warmup_end-ref_warmup_start);
     if(fd >= 0) {
        printf(" 0x186: 0x%08, 0x187: 0x%08, 0x188: 0x%08, 0x189: 0x%08\n",data[0],data[1],data[2],data[3]);
     }
     //for(int32_t i = 0; i != 4; ++i) {
     //    printf(" %lu ", corrected_pmc_delta(core_counters_end[i],core_counters_start[i],core_ctr_width));
     //}
     // Begin deallocation of arrays
     _mm_free(core_avg_f);        _mm_free(core_util);
     _mm_free(TSC_f64);           _mm_free(FCIns_f64);
     _mm_free(FCAct_f64);         _mm_free(FCRef_f64);
     _mm_free(PMC3_f64);          _mm_free(PMC2_f64);
     _mm_free(PMC1_f64);          _mm_free(PMC0_f64);
      _mm_free(Core_avg_f_mnt);
      _mm_free(Core_util_mnt);
      _mm_free(TSC_mnt);
      _mm_free(FCIns_mnt);
      _mm_free(FCAct_mnt);
      _mm_free(FCRef_mnt);
      _mm_free(PMC3_mnt);
      _mm_free(PMC2_mnt);
      _mm_free(PMC1_mnt);
      _mm_free(PMC0_mnt);
}



int main(int argc, char * argv[]) {

    perf_test_ComputeGrassParamEq_zmm16r4_call_scope_looped100x10000_non_instr_thirmo();
    return (0);
}
