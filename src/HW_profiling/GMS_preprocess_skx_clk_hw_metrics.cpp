
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
#include <omp.h>
#include "GMS_skx_clk_hw_events_metrics.hpp"



/*
   CPU operating frequency (in GHz)
*/
void
skx_cpu_operating_freq_samples(const double * __restrict __attribute__((aligned(64))) a,
			       const double * __restrict __attribute__((aligned(64))) b,
			       const double * __restrict __attribute__((aligned(64))) c,
			       double * __restrict __attribute__((aligned(64))) d,
			       const int32_t data_len) {

       
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,d,data_len) private(i) \
        aligned(a:64,b,c,d) linear(i:1) unroll partial(2) if(data_len>=100000)
           for(int32_t i = 0; i != data_len; ++i) {
               d[i] = skx_cpu_operating_freq(a[i],
	                                     b[i],
					     c[i]);
	   }       
}


/*
   CPU utilization (percentage) all cores.
*/
void
skx_cpu_utilization_samples(const double * __restrict __attribute__((aligned(64))) a,
                            const double * __restrict __attribute__((aligned(64))) b,
			    double * __restrict __attribute__((aligned(64))) c,
			    const int32_t data_len) {

      
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_cpu_utilization(a[i],
	                                  b[i]);
					
          }
}


/*
    CPU utilization (percentage) in kernel mode (all cores).
*/
void
skx_cpu_utlization_kernel_samples(const double * __restrict __attribute__((aligned(64))) a,
                                  const double * __restrict __attribute__((aligned(64))) b,
			          double * __restrict c,
			          const int32_t) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif        
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_cpu_utilization_kernel(a[i],
	                                         b[i]);
					
          }
}


/*
    Cycles Per Instruction (CPI).
*/
void
skx_cycles_per_instr_samples(const double * __restrict __attribute__((aligned(64))) a,
                             const double * __restrict __attribute__((aligned(64))) b,
			     double * __restrict __attribute__((aligned(64))) c,
			     const int32_t data_len) {

 #if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif      
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_cycles_per_instr(a[i],
	                                   b[i]);
					
          }	 
}


/*
    Cycles Per Instruction (CPI) kernel mode.
*/
void
skx_cycles_per_instr_kernel_samples(const double * __restrict __attribute__((aligned(64))) a,
                                    const double * __restrict __attribute__((aligned(64))) b,
			            double * __restrict __attribute__((aligned(64))) c,
			            const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif       
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_cycles_per_instr_kernel(a[i],
	                                          b[i]);
					
          }	 
}


/*
    EMON event multiplexing reliability (>95% --- means satisfying ratio).
*/
void
skx_mux_reliability_samples(const double * __restrict  __attribute__((aligned(64))) a,
                            const double * __restrict  __attribute__((aligned(64))) b,
			    double * __restrict  __attribute__((aligned(64))) c,
			    const int32_t data_len) {
			    
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_emon_mux_reliability(a[i],
	                                       b[i]);
					
          }	       
}


/*
      Branch mispredict ratio.
*/
void
skx_branch_mispred_ratio_samples(const double * __restrict __attribute__((aligned(64))) a,
                                 const double * __restrict __attribute__((aligned(64))) b,
			         double * __restrict __attribute__((aligned(64))) c,
			         const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_branch_mispredict_ratio(a[i],
	                                          b[i]);
					
          }	
}


/*
    Loads per instruction.
*/
void
skx_loads_per_instr_samples(const double * __restrict  __attribute__((aligned(64))) a,
                            const double * __restrict  __attribute__((aligned(64))) b,
			    double * __restrict  __attribute__((aligned(64))) c,
			    const int32_t data_len data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_loads_per_instr(a[i],
	                                  b[i]);
					
          }
}


/*
    Stores per instruction.
*/
void
skx_stores_per_instr_samples(const double * __restrict  __attribute__((aligned(64))) a,
                             const double * __restrict  __attribute__((aligned(64))) b,
			     double * __restrict  __attribute__((aligned(64))) c,
			     const int32_t data_len data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_stores_per_instr(a[i],
	                                  b[i]);
					
          }
}


/*
   Memory operations per instruction.
*/
void
skx_mem_ops_instr_samples(const double * __restrict  __attribute__((aligned(64))) a,
                          const double * __restrict  __attribute__((aligned(64))) b,
		          const double * __restrict  __attribute__((aligned(64))) c,
		          double * __restrict  __attribute__((aligned(64))) d,
		          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,d,data_len) private(i) \
        aligned(a:64,b,c,d) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               d[i] = skx_mem_ops_per_instr(a[i],
	                                    b[i],
					    c[i]);
					
          }
}


/*
    Locks retired per instruction.
*/
void
skx_locks_per_instr_samples(const double * __restrict  __attribute__((aligned(64))) a,
                           const double * __restrict  __attribute__((aligned(64))) b,
			   double * __restrict  __attribute__((aligned(64))) c,
			   const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_locks_per_instr(a[i],
	                                  b[i]);
					
          }
}


/*
   Uncacheable reads per instruction.
*/
void
skx_uncacheable_reads_instr_samples(const double * __restrict __attribute__((aligned(64))) a,
                                    const double * __restrict __attribute__((aligned(64))) b,
			            double * __restrict __attribute__((aligned(64))) c,
			            const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_uncacheable_reads_per_instr(a[i],
	                                              b[i]);
					
          }
}


/*
    Streaming-stores (full line) per instruction.
*/
void
skx_stream_stores_fl_instr_samples(const double * __restrict __attribute__((aligned(64))) a,
                                   const double * __restrict __attribute__((aligned(64))) b,
			           double * __restrict __attribute__((aligned(64))) c,
			           const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_streaming_stores_fl_per_instr(a[i],
	                                              b[i]);
					
          }
}


/*
     Streaming-stores (partial line) per instruction.
*/
void
skx_stream_stores_pl_instr_samples(const double * __restrict __attribute__((aligned(64))) a,
                                   const double * __restrict __attribute__((aligned(64))) b,
			           double * __restrict __attribute__((aligned(64))) c,
			           const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_streaming_stores_pl_per_instr(a[i],
	                                              b[i]);
					
          }
}


/*
     L1D$ Misses per instruction including:
     -- data
     -- reads for ownership (rfo) with prefetches
*/
void
skx_L1D_misses_instr_samples(const double * __restrict __attribute__((aligned(64))) a,
                             const double * __restrict __attribute__((aligned(64))) b,
			     double * __restrict __attribute__((aligned(64))) c,
			     const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L1D_misses_per_instr(a[i],
	                                       b[i]);
					
          }
}


/*
     L1D$ Demand data read hit per instruction.
*/
void
skx_L1D_hits_instr_samples(  const double * __restrict __attribute__((aligned(64))) a,
                             const double * __restrict __attribute__((aligned(64))) b,
			     double * __restrict __attribute__((aligned(64))) c,
			     const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L1D_hits_per_instr(a[i],
	                                     b[i]);
					
          }
}


/*
    L1$ code read misses (including prefetches) per instruction.
*/
void
skx_L1I_read_misses_instr_samples(const double * __restrict __attribute__((aligned(64))) a,
                                  const double * __restrict __attribute__((aligned(64))) b,
			          double * __restrict __attribute__((aligned(64))) c,
			          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L1I_read_misses_per_instr(a[i],
	                                            b[i]);
					
          }
}


/*
     L2 Demand data read hits per instruction.
*/
void
skx_L2_data_read_hits_instr_samples(const double * __restrict __attribute__((aligned(64))) a,
                                      const double * __restrict __attribute__((aligned(64))) b,
			              double * __restrict __attribute__((aligned(64))) c,
			              const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L2_data_read_hits_per_instr(a[i],
	                                            b[i]);
					
          }
}


/*
      L2 Misses per instruction including:
     -- code
     -- data
     -- read for ownership and prefetches
*/
void
skx_L2_all_misses_instr_samples(const double * __restrict __attribute__((aligned(64))) a,
                                const double * __restrict __attribute__((aligned(64))) b,
			        double * __restrict __attribute__((aligned(64))) c,
			        const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L2_all_misses_per_instr(a[i],
	                                          b[i]);
					
          }
}


/*
       L2 demand data read misses per instruction.
*/
void
skx_L2_demand_data_read_mpi_samples(const double * __restrict __attribute__((aligned(64))) a,
                                    const double * __restrict __attribute__((aligned(64))) b,
			            double * __restrict __attribute__((aligned(64))) c,
			            const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L2_demand_data_read_mpi(a[i],
	                                          b[i]);
					
          }
}


/*
     L2 demand data code misses per instruction.
*/
void
skx_L2_demand_code_mpi_samples(const double * __restrict __attribute__((aligned(64))) a,
                               const double * __restrict __attribute__((aligned(64))) b,
			       double * __restrict __attribute__((aligned(64))) c,
			       const int32_t) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L2_demand_code_mpi(a[i],
	                                     b[i]);
					
          }
}


/*
    L2 Any local request that HITM in a sibling core (per instruction).
*/
void
skx_L2_req_hitm_sibling_core_samples(const double * __restrict  __attribute__((aligned(64))) a,
                                     const double * __restrict  __attribute__((aligned(64))) b,
			             double * __restrict  __attribute__((aligned(64))) c,
			             const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L2_request_hitm_sibling_core(a[i],
	                                               b[i]);
					
          }
}


/*
     L2 (percentage) of all lines evicted that are unused prefetches.
*/
void
skx_L2_lines_evict_unused_prefetch_samples(const double * __restrict  __attribute__((aligned(64))) a,
                                           const double * __restrict  __attribute__((aligned(64))) b,
					   const double * __restrict  __attribute__((aligned(64))) c,
					   const double * __restrict  __attribute__((aligned(64))) d,
					   double * __restrict  __attribute__((aligned(64))) e,
					   const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,,d,e,data_len) private(i) \ 
        aligned(a:64,b,c,d,e) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               e[i] = skx_L2_lines_evict_unused_prefetch(a[i],
	                                                 b[i],
							 c[i],
							 d[i]);
					
          }
}


/*
    L2 percentage of L2 evictions that are allocated into L3$.
*/
void
skx_L2_evict_L3_alloc_samples(const double * __restrict  __attribute__((aligned(64))) a,
                              const double * __restrict  __attribute__((aligned(64))) b,
			      double * __restrict  __attribute__((aligned(64))) c,
			      const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L2_evict_L3_allocated(a[i],
	                                        b[i]);
					
          }
}


/*
    L2 percentage of L2 evictions that are allocated into L3$.
*/
void
skx_L2_evict_L3_no_alloc_samples(const double * __restrict  __attribute__((aligned(64))) a,
                                 const double * __restrict  __attribute__((aligned(64))) b,
			         double * __restrict  __attribute__((aligned(64))) c,
			         const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) shared(a,b,c,data_len) private(i) \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L2_evict_L3_not_allocated(a[i],
	                                            b[i]);
					
          }
}


/*
     LLC code references per instruction (L3 prefetch excluded).
*/
void
skx_L3_code_references_instr_samples(const double * __restrict  __attribute__((aligned(64))) a,
                                  const double * __restrict  __attribute__((aligned(64))) b,
			          double * __restrict  __attribute__((aligned(64))) c,
			          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L3_code_references_per_instr(a[i],
	                                               b[i]);
					
          }
}


/*
     LLC data read references per instruction (L3 prefetch excluded).
*/
void
skx_L3_data_references_instr_samples(const double * __restrict  __attribute__((aligned(64))) a,
                                  const double * __restrict  __attribute__((aligned(64))) b,
			          double * __restrict  __attribute__((aligned(64))) c,
			          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L3_data_references_per_instr(a[i],
	                                               b[i]);
					
          }
}


/*
      LLC RFO references per instrctions (L3 prefetch excluded).
*/
void
skx_L3_rfo_references_instr_samples( const double * __restrict  __attribute__((aligned(64))) a,
                                  const double * __restrict  __attribute__((aligned(64))) b,
			          double * __restrict  __attribute__((aligned(64))) c,
			          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L3_rfo_references_per_instr(a[i],
	                                              b[i]);
					
          }
}


/*
     LLC Misses Per Instruction (includes code and data and rfo with prefetches).
*/
void
skx_L3_all_mpi_samples(const double * __restrict __attribute__((aligned(64))) a,
               const double * __restrict __attribute__((aligned(64))) b,
	       const double * __restrict __attribute__((aligned(64))) c,
               const double * __restrict __attribute__((aligned(64))) d,
	       double * __restrict __attribute__((aligned(64))) e,
	       const int32_t data_len)  {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,,d,e,data_len) private(i) \ 
        aligned(a:64,b,c,d,e) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               e[i] = skx_L3_all_mpi(a[i],
	                             b[i],
				     c[i],
				     d[i]);
					
          }
}


/*
     LLC data read Misses Per Instruction (demand and prefetch).
*/
void
skx_L3_data_read_mpi_samples(const double * __restrict __attribute__((aligned(64))) a,
                          const double * __restrict __attribute__((aligned(64))) b,
			  double * __restrict  __attribute__((aligned(64))) c,
			  const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L3_data_read_mpi(a[i],
	                                   b[i]);
					
          }
}


/*
      LLC RFO read Misses Per Instruction (demand and prefetch).  
*/
void
skx_L3_rfo_data_read_mpi_samples(const double * __restrict __attribute__((aligned(64))) a,
                              const double * __restrict __attribute__((aligned(64))) b,
			      double * __restrict  __attribute__((aligned(64))) c,
			      const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L3_rfo_read_mpi(a[i],
	                                  b[i]);
					
          }
}


/*
     LLC code read MPI (demand and prefetch).
*/
void
skx_L3_total_hitm_instr_samples(const double * __restrict  __attribute__((aligned(64))) a,
                                const double * __restrict  __attribute__((aligned(64))) b,
			        double * __restrict   __attribute__((aligned(64))) c,
			        const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L3_total_hitm_per_instr(a[i],
	                                          b[i]);
					
          }
}


/*
     LLC total HIT clean line forwards (per instr) (excludes LLC prefetches)
*/
void
skx_L3_total_hitm_clean_lines_instr_samples(const double * __restrict __attribute__((aligned(64))) a,
                                            const double * __restrict __attribute__((aligned(64))) b,
			                    double * __restrict __attribute__((aligned(64))) c,
			                    const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L3_total_hitm_clean_lines_per_instr(a[i],
	                                                      b[i]);
					
          }
}


/*
    Average LLC data read (demand+prefetch) miss latency (in ns).  
*/
void
skx_L3_avg_data_read_ns_samples(const double * __restrict __attribute__((aligned(64))) a,
                                const double * __restrict __attribute__((aligned(64))) b,
			        const double * __restrict __attribute__((aligned(64))) c,
                                const double * __restrict  __attribute__((aligned(64))) d,
			        double * __restrict __attribute__((aligned(64))) e,
			        const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,,d,e,data_len) private(i) \ 
        aligned(a:64,b,c,d,e) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               e[i] = skx_L3_avg_data_read_ns(a[i],
	                                      b[i],
				              c[i],
				              d[i]);
					
          }
}


/*
      Average LLC data read (demand and prefetch) miss latency (in UNCORE clk).
*/
void
skx_L3_avg_data_read_unc_clk_samples(const double * __restrict __attribute__((aligned(64))) a,
                                     const double * __restrict __attribute__((aligned(64))) b,
			             double * __restrict __attribute__((aligned(64))) c,
			             const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L3_avg_data_read_unc_clk(a[i],
	                                           b[i]);
					
          }
}


/*
      Average LLC data read (demand+prefetch) miss latency for LOCAL requests (in ns)
*/
void
skx_L3_avg_data_read_loc_ns_samples(const double * __restrict  __attribute__((aligned(64))) a,
                                    const double * __restrict  __attribute__((aligned(64))) b,
			            const double * __restrict  __attribute__((aligned(64))) c,
                                    const double * __restrict  __attribute__((aligned(64))) d,
			            double * __restrict  __attribute__((aligned(64))) e,
			            const int32_t) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,,d,e,data_len) private(i) \ 
        aligned(a:64,b,c,d,e) linear(i:1) unroll partial(10) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               e[i] = skx_L3_avg_data_read_loc_req_ns(a[i],
	                                              b[i],
				                      c[i],
				                      d[i]);
					
          }
}


/*
    Average LLC data read (demand+prefetch) miss latency for REMOTE requests (in ns).
*/
void
skx_L3_avg_data_read_rem_ns_samples(const double * __restrict  __attribute__((aligned(64))) a,
                                    const double * __restrict  __attribute__((aligned(64))) b,
			            const double * __restrict  __attribute__((aligned(64))) c,
                                    const double * __restrict  __attribute__((aligned(64))) d,
			            double * __restrict  __attribute__((aligned(64))) e,
			            const int32_t) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,,d,e,data_len) private(i) \ 
        aligned(a:64,b,c,d,e) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               e[i] = skx_L3_avg_data_read_rem_req_ns(a[i],
	                                              b[i],
				                      c[i],
				                      d[i]);
					
          }
}


/*
      Average LLC data read (demand+prefetch) miss latency  for LOCAL requests (in UNCORE clk)
*/
void
skx_L3_avg_data_read_loc_unc_clk_samples(const double * __restrict __attribute__((aligned(64))) a,
                                         const double * __restrict __attribute__((aligned(64))) b,
			                 double * __restrict __attribute__((aligned(64))) c,
			                 const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L3_avg_data_read_loc_unc_clk(a[i],
	                                               b[i]);
					
          }
}


/*
     Average LLC data read (demand+prefetch) miss latency  for REMOTE requests (in UNCORE clk)
*/
void
skx_L3_avg_data_read_rem_unc_clk_samples(const double * __restrict __attribute__((aligned(64))) a,
                                         const double * __restrict __attribute__((aligned(64))) b,
			                 double * __restrict __attribute__((aligned(64))) c,
			                 const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_L3_avg_data_read_rem_unc_clk(a[i],
	                                               b[i]);
					
          }
}


/*
    ITLB MPI
*/
void
skx_ITLB_mpi_data_samples(const double * __restrict __attribute__((aligned(64))) a,
                          const double * __restrict __attribute__((aligned(64))) b,
		          double * __restrict __attribute__((aligned(64))) c,
		          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_ITLB_mpi(a[i],
	                           b[i]);
					
          }
}


/*
   ITLB large page MPI
*/
void
skx_ITLB_2M4M_mpi_samples(const double * __restrict __attribute__((aligned(64))) a,
                          const double * __restrict __attribute__((aligned(64))) b,
		          double * __restrict __attribute__((aligned(64))) c,
		          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_ITLB_2M_4M_mpi(a[i],
	                                 b[i]);
					
          }
}


/*
     DTLB load MPI.
*/
void
skx_DTLB_load_mpi_samples(const double * __restrict __attribute__((aligned(64))) a,
                          const double * __restrict __attribute__((aligned(64))) b,
		          double * __restrict __attribute__((aligned(64))) c,
		          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_DTLB_load_mpi(a[i],
	                                b[i]);
					
          }
}


/*
    DTLB large page load MPI.
*/
void
skx_DTLB_2M4M_mpi_samples(const double * __restrict __attribute__((aligned(64))) a,
                          const double * __restrict __attribute__((aligned(64))) b,
		          double * __restrict __attribute__((aligned(64))) c,
		          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_DTLB_2M_4M_mpi(a[i],
	                                b[i]);
					
          }
}


/*
   DTLB store MPI.
*/
void
skx_DTLB_store_mpi_samples(const double * __restrict __attribute__((aligned(64))) a,
                          const double * __restrict __attribute__((aligned(64))) b,
		          double * __restrict __attribute__((aligned(64))) c,
		          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_DTLB_store_mpi(a[i],
	                                b[i]);
					
          }
}


/*
    DTLB load miss latency (in core clks)
*/
void
skx_DTLB_load_miss_clks_samples(const double * __restrict  __attribute__((aligned(64))) a,
                                const double * __restrict  __attribute__((aligned(64))) b,
		                double * __restrict  __attribute__((aligned(64))) c,
		                const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_DTLB_load_miss_clks(a[i],
	                                      b[i]);
					
          }
}


/*
     DTLB store miss latency (in core clks).
*/
void
skx_DTLB_store_miss_clks_samples(const double * __restrict  __attribute__((aligned(64))) a,
                                const double * __restrict  __attribute__((aligned(64))) b,
		                double * __restrict  __attribute__((aligned(64))) c,
		                const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_DTLB_store_miss_clks(a[i],
	                                      b[i]);
					
          }
}


/*
      ITLB miss latency (in core clks).
*/
void
skx_ITLB_miss_latency_clks_samples(const double * __restrict  __attribute__((aligned(64))) a,
                                   const double * __restrict   __attribute__((aligned(64))) b,
		                   double * __restrict  __attribute__((aligned(64))) c,
		                   const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_ITLB_miss_latency_clks(a[i],
	                                      b[i]);
					
          }
}


/*
    NUMA percentage of Reads addressed to local DRAM.
*/
void
skx_numa_reads_local_dram_samples(const double * __restrict   __attribute__((aligned(64))) a,
                                  const double * __restrict   __attribute__((aligned(64))) b,
		                  double * __restrict  __attribute__((aligned(64))) c,
		                  const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_nume_reads_local_dram(a[i],
	                                        b[i]);
					
          }
}


/*
    NUMA percentage of Reads addressed to remote  DRAM.
*/
void
skx_numa_reads_remote_dram_samples(const double * __restrict   __attribute__((aligned(64))) a,
                                  const double * __restrict   __attribute__((aligned(64))) b,
		                  double * __restrict  __attribute__((aligned(64))) c,
		                  const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_nume_reads_remote_dram(a[i],
	                                        b[i]);
					
          }
}


/*
     Uncore Frequency Ghz.
*/
void
skx_uncore_frequency_samples(const double * __restrict  __attribute__((aligned(64))) a,
                             const double * __restrict  __attribute__((aligned(64))) b,
			     const double * __restrict  __attribute__((aligned(64))) c,
		             double * __restrict d,
		             const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(i)					\
        aligned(a:64,b,c,d) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               d[i] = skx_uncore_frequency_ghz(a[i],
	                                       b[i],
					       c[i]);
					
          }
}


/*
    UPI speed - GT/s.
*/
void
skx_UPI_speed_samples(const double * __restrict __attribute__((aligned(64))) a,
                      const double * __restrict __attribute__((aligned(64))) b,
	              double * __restrict __attribute__((aligned(64))) c,
	              const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_UPI_speed(a[i],
	                            b[i]);
					
          }
}


/*
    UPI Data transmit BW (MB/sec) (only data)
*/
void
skx_UPI_data_bw_samples(const double * __restrict __attribute__((aligned(64))) a,
                        const double * __restrict __attribute__((aligned(64))) b,
	                double * __restrict __attribute__((aligned(64))) c,
	                const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i) 	 \
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_UPI_data_bw(a[i],
	                            b[i]);
					
          }
}


/*
     UPI Total transmit BW (MB/sec) (includes control)
*/
void
skx_UPI_tx_total_bw_samples(const double * __restrict __attribute__((aligned(64))) a,
                            const double * __restrict __attribute__((aligned(64))) b,
			    const double * __restrict __attribute__((aligned(64))) c,
	                    double * __restrict __attribute__((aligned(64))) d,
	                    const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(i)					\
        aligned(a:64,b,c,d) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               d[i] = skx_UPI_transmit_total_bw(a[i],
	                                       b[i],
					       c[i]);
					
          }
}


/*
     UPI Transmit utilization percentage (includes control).
   Percentage of time the processor is communicating with sibling processors.
*/
void
skx_UPI_tx_utlization_samples(const double * __restrict __attribute__((aligned(64))) a,
                              const double * __restrict __attribute__((aligned(64))) b,
			      const double * __restrict __attribute__((aligned(64))) c,
			      const double * __restrict __attribute__((aligned(64))) d,
			      const double * __restrict __attribute__((aligned(64))) e,
			      const double * __restrict __attribute__((aligned(64))) f,
	                      double * __restrict __attribute__((aligned(64))) g,
	                      const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,,d,e,f,g,data_len) private(i) \ 
        aligned(a:64,b,c,d,e,f,g) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               g[i] = skx_UPI_transmit_utilization(a[i],
	                                           b[i],
				                   c[i],
				                   d[i],
						   e[i],
						   f[i]);
					
          }
}


/*
   UPI percentage of  cycles transmit link is half-width (L0p) 
*/
void
skx_UPI_half_width_link_tx_cycles_samples( const double * __restrict __attribute__((aligned(64))) a,
			                   const double * __restrict __attribute__((aligned(64))) b,
			                   const double * __restrict __attribute__((aligned(64))) c,
	                                   double * __restrict __attribute__((aligned(64))) d,
	                                   const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(i)					\
        aligned(a:64,b,c,d) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               d[i] = skx_half_width_link_tx_cycles(a[i],
	                                            b[i],
					            c[i]);
					
          }
}


/*
    UPI percentage of  cycles receive link is half-width (L0p)
*/
void
skx_UPI_half_width_link_rx_cycles_samples( const double * __restrict __attribute__((aligned(64))) a,
			                   const double * __restrict __attribute__((aligned(64))) b,
			                   const double * __restrict __attribute__((aligned(64))) c,
	                                   double * __restrict __attribute__((aligned(64))) d,
	                                   const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(i)					\
        aligned(a:64,b,c,d) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               d[i] = skx_half_width_link_rx_cycles(a[i],
	                                            b[i],
					            c[i]);
					
          }
}


/*
    HA - Reads vs. all requests
*/
void
skx_HA_reads_vs_all_reqs_samples(const double * __restrict __attribute__((aligned(64))) a,
                                 const double * __restrict __attribute__((aligned(64))) b,
	                         double * __restrict __attribute__((aligned(64))) c,
	                         const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
  aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_HA_reads_vs_all_reqs(a[i],
	                                       b[i]);
					       
					
          }
}


/*
    HA - Writes vs. all requests
*/
void
skx_HA_writes_vs_all_reqs_samples(const double * __restrict __attribute__((aligned(64))) a,
                                 const double * __restrict __attribute__((aligned(64))) b,
	                         double * __restrict __attribute__((aligned(64))) c,
	                         const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_HA_writes_vs_all_reqs(a[i],
	                                       b[i]);
					       
					
          }
}


/*
      HA percentage of all reads that are local.
*/
void
skx_HA_all_reads_local_samples(const double * __restrict  __attribute__((aligned(64))) a,
                               const double * __restrict  __attribute__((aligned(64))) b,
	                       double * __restrict  __attribute__((aligned(64))) c,
	                       const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_HA_all_reads_local(a[i],
	                                       b[i]);
					       
					
          }
}


/*
     HA percentage of all writes that are local.
*/
void
skx_HA_all_writes_local_samples(const double * __restrict  __attribute__((aligned(64))) a,
                               const double * __restrict  __attribute__((aligned(64))) b,
	                       double * __restrict  __attribute__((aligned(64))) c,
	                       const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_HA_all_writes_local(a[i],
	                                       b[i]);
					       
					
          }
}


/*
     HA conflict responses per instruction.
*/
void
skx_HA_conflict_resp_samples(const double * __restrict __attribute__((aligned(64))) a,
                             const double * __restrict __attribute__((aligned(64))) b,
	                     double * __restrict __attribute__((aligned(64))) c,
	                     const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_HA_conflict_resp(a[i],
	                                   b[i])
					       
					
          }
}


/*
    HA directory lookups that spawned a snoop (per instruction)
*/
void
skx_HA_dir_lookup_snoop_samples(const double * __restrict __attribute__((aligned(64))) a,
                               const double * __restrict __attribute__((aligned(64))) b,
	                       double * __restrict __attribute__((aligned(64))) c,
	                       const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_HA_dir_lookup_snoop(a[i],
	                                     b[i]);
					       
					
          }
}


/*
    HA directory lookups that did not spawn a snoop (per instruction).
*/
void
skx_HA_dir_lookup_no_snoop_samples(const double * __restrict __attribute__((aligned(64))) a,
                                   const double * __restrict __attribute__((aligned(64))) b,
	                           double * __restrict __attribute__((aligned(64))) c,
	                           const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_HA_dir_lookup_no_snoop(a[i],
	                                     b[i]);
					       
					
          }
}


/*
    M2M directory updates (per instruction). 
*/
void
skx_M2M_dir_updates_samples(const double * __restrict __attribute__((aligned(64))) a,
                            const double * __restrict __attribute__((aligned(64))) b,
			    const double * __restrict __attribute__((aligned(64))) c,
			    const double * __restrict __attribute__((aligned(64))) d,
	                    double * __restrict __attribute__((aligned(64))) e,
	                    const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,data_len) private(i)					\
        aligned(a:64,b,c,d,e) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               e[i] = skx_M2M_dir_update(a[i],
	                                 b[i],
					 c[i],
					 d[i]);
					
          }
}


/*
    M2M extra reads from XPT-UPI prefetches (per instruction).
*/
void
skx_M2M_reads_XPT_UPI_prefetch_samples(const double * __restrict __attribute__((aligned(64))) a,
                                       const double * __restrict __attribute__((aligned(64))) b,
			               const double * __restrict __attribute__((aligned(64))) c,
			               double * __restrict __attribute__((aligned(64))) d,
	                               const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(i)					\
        aligned(a:64,b,c,d) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               d[i] = skx_M2M_reads_XPT_UPI_prefetch(a[i],
	                                             b[i],
					             c[i]);
				
					
          }
}


/*
     DDR data rate (MT/sec).
*/
void
skx_DDR_date_rate_samples(const double * __restrict  __attribute__((aligned(64))) a,
                          const double * __restrict  __attribute__((aligned(64))) b,
	                  double * __restrict  __attribute__((aligned(64))) c,
	                  const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_DDR_data_rate(a[i],
	                                b[i]);
					       
					
          }
}


/*
    Memory bandwidth read (MB/sec).
*/
void
skx_memory_read_bw_samples(const double * __restrict  __attribute__((aligned(64))) a,
                           const double * __restrict  __attribute__((aligned(64))) b,
	                   double * __restrict  __attribute__((aligned(64))) c,
	                  const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_memory_read_bw(a[i],
	                                b[i]);
					       
					
          }
}


/*
    Load instructions per memory bandwidth.
*/
void
skx_load_mem_inst_mem_bw_samples(const double * __restrict __attribute__((aligned(64))) a,
                                 const double * __restrict  __attribute__((aligned(64))) b,
			         const double * __restrict  __attribute__((aligned(64))) c,
			         double * __restrict  __attribute__((aligned(64))) d,
	                         const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(i)					\
        aligned(a:64,b,c,d) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               d[i] = skx_load_mem_instr_mem_bw(a[i],
	                                        b[i],
					        c[i]);
				
					
          }
}


/*
    Memory bandwidth write  (MB/sec).
*/
void
skx_memory_write_bw_samples(const double * __restrict  __attribute__((aligned(64))) a,
                            const double * __restrict __attribute__((aligned(64))) b,
	                    double * __restrict __attribute__((aligned(64))) c,
	                    const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_memory_write_bw(a[i],
	                                b[i]);
					       
					
          }
}


/*
    Store instructions per memory bandwidth.
*/
void
skx_load_mem_inst_mem_bw_samples(const double * __restrict __attribute__((aligned(64))) a,
                                 const double * __restrict  __attribute__((aligned(64))) b,
			         const double * __restrict  __attribute__((aligned(64))) c,
			         double * __restrict  __attribute__((aligned(64))) d,
	                         const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(i)					\
        aligned(a:64,b,c,d) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               d[i] = skx_load_mem_instr_mem_bw(a[i],
	                                        b[i],
					        c[i]);
				
					
          }
}


/*
    Memory bandwidth total (MB/sec).
*/
void
skx_mem_bw_total_samples(const double * __restrict  __attribute__((aligned(64))) a,
                         const double * __restrict  __attribute__((aligned(64))) b,
			 const double * __restrict  __attribute__((aligned(64))) c,
			 double * __restrict  __attribute__((aligned(64))) d,
	                 const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(i)					\
        aligned(a:64,b,c,d) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               d[i] = skx_memory_bw_total(a[i],
	                                  b[i],
					  c[i]);
				
					
          }
}


/*
     Load and store instructions per total memory bandwidth.
*/
void
skx_total_mem_inst_mem_bw_samples(const double * __restrict __attribute__((aligned(64))) a,
                                  const double * __restrict __attribute__((aligned(64))) b,
			          const double * __restrict __attribute__((aligned(64))) c,
			          const double * __restrict __attribute__((aligned(64))) d,
			          const double * __restrict __attribute__((aligned(64))) e,
			          double * __restrict __attribute__((aligned(64))) f,
	                          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,f,data_len) private(i)					\
        aligned(a:64,b,c,d,e,f) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               f[i] = skx_total_mem_instr_mem_bw(a[i],
	                                 b[i],
					 c[i],
					 d[i],
					 e[i]);
					
          }
}


/*
    Memory extra read b/w due to XPT prefetches (MB/sec).
*/
void
skx_XPT_mem_bw_prefetch_samples(const double * __restrict __attribute__((aligned(64))) a,
                                const double * __restrict __attribute__((aligned(64))) b,
			        const double * __restrict __attribute__((aligned(64))) c,
			        double * __restrict __attribute__((aligned(64))) d,
	                        const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(i)					\
        aligned(a:64,b,c,d) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               d[i] = skx_mem_bw_xpt_prefetch(a[i],
	                                  b[i],
					  c[i]);
				
					
          }
}


/*
    Memory extra write b/w due to directory updates (MB/sec).
*/
void
skx_mem_bw_dir_update_samples(const double * __restrict   __attribute__((aligned(64))) a,
                              const double * __restrict   __attribute__((aligned(64))) b,
			      const double * __restrict   __attribute__((aligned(64))) c,
			      const double * __restrict   __attribute__((aligned(64))) d,
			      double * __restrict  __attribute__((aligned(64))) e,
	                      const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,data_len) private(i)					\
        aligned(a:64,b,c,d,e) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               e[i] = skx_mem_bw_dir_update(a[i],
	                                  b[i],
					  c[i],
					  d[i]);
				
					
          }
}


/*
      DRAM RPQ read latency (ns).
*/
void
skx_DRAM_RPQ_read_latency_samples(const double * __restrict  __attribute__((aligned(64))) a,
                                  const double * __restrict  __attribute__((aligned(64))) b,
			          const double * __restrict  __attribute__((aligned(64))) c,
			          const double * __restrict  __attribute__((aligned(64))) d,
			          double * __restrict  __attribute__((aligned(64))) e,
	                          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,data_len) private(i)					\
        aligned(a:64,b,c,d,e) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               e[i] = skx_DRAM_rpq_read_latency_ns(a[i],
	                                           b[i],
					           c[i],
					           d[i]);
				
					
          }
}


/*
     DRAM RPQ write latency (ns).
*/
void
skx_DRAM_RPQ_write_latency_samples(const double * __restrict  __attribute__((aligned(64))) a,
                                  const double * __restrict  __attribute__((aligned(64))) b,
			          const double * __restrict  __attribute__((aligned(64))) c,
			          const double * __restrict  __attribute__((aligned(64))) d,
			          double * __restrict  __attribute__((aligned(64))) e,
	                          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,data_len) private(i)					\
        aligned(a:64,b,c,d,e) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               e[i] = skx_DRAM_rpq_write_latency_ns(a[i],
	                                           b[i],
					           c[i],
					           d[i]);
				
					
          }
}


/*
    Memory average number of entries in each read Q (RPQ)
*/
void
skx_RPQ_mem_avg_writes_samples(const double * __restrict  __attribute__((aligned(64))) a,
                               const double * __restrict  __attribute__((aligned(64))) b,
	                       double * __restrict  __attribute__((aligned(64))) c,
	                       const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_mem_avg_rpq_read(a[i],
	                                b[i]);
					       
					
          }
}


/*
    memory average number of entries in each write Q (WPQ).
*/
void
skx_RPQ_mem_avg_reads_samples(const double * __restrict  __attribute__((aligned(64))) a,
                               const double * __restrict  __attribute__((aligned(64))) b,
	                       double * __restrict  __attribute__((aligned(64))) c,
	                       const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_mem_avg_rpq_write(a[i],
	                                b[i]);
					       
					
          }
}


/*
     I/O bandwidth disk or network writes (MB/sec).
*/
void
skx_IO_disk_or_net_writes_bw_samples(const double * __restrict __attribute__((aligned(64))) a,
                                     const double * __restrict __attribute__((aligned(64))) b,
			             const double * __restrict __attribute__((aligned(64))) c,
			             const double * __restrict __attribute__((aligned(64))) d,
			             const double * __restrict __attribute__((aligned(64))) e,
			             double * __restrict __attribute__((aligned(64))) f,
	                             const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,f,data_len) private(i)					\
        aligned(a:64,b,c,d,e,f) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               f[i] = skx_IO_disk_or_net_bw_writes(a[i],
	                                           b[i],
					           c[i],
					           d[i],
					           e[i]);
					
          }
}


/*
     I/O bandwidth disk or network reads (MB/sec).
*/
void
skx_IO_disk_or_net_reads_bw_samples( const double * __restrict __attribute__((aligned(64))) a,
                                     const double * __restrict __attribute__((aligned(64))) b,
			             const double * __restrict __attribute__((aligned(64))) c,
			             const double * __restrict __attribute__((aligned(64))) d,
			             const double * __restrict __attribute__((aligned(64))) e,
			             double * __restrict __attribute__((aligned(64))) f,
	                             const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,f,data_len) private(i)					\
        aligned(a:64,b,c,d,e,f) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               f[i] = skx_IO_disk_or_net_bw_reads(a[i],
	                                           b[i],
					           c[i],
					           d[i],
					           e[i]);
					
          }
}


/*
   I/O bandwidth disk or network (MB/sec)
*/
void
skx_IO_total_bw_samples(const double * __restrict  __attribute__((aligned(64))) a,
                        const double * __restrict  __attribute__((aligned(64))) b,
			const double * __restrict  __attribute__((aligned(64))) c,
			const double * __restrict  __attribute__((aligned(64))) d,
			const double * __restrict  __attribute__((aligned(64))) e,
			const double * __restrict  __attribute__((aligned(64))) f,
			const double * __restrict  __attribute__((aligned(64))) g,
			const double * __restrict  __attribute__((aligned(64))) h,
			const double * __restrict  __attribute__((aligned(64))) i,
			double * __restrict  __attribute__((aligned(64))) j,
	                const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,f,g,h,i,j,data_len) private(idx)					\
        aligned(a:64,b,c,d,e,f,g,h,i,j) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               j[idx] = skx_IO_disk__net_total(a[idx],
	                                     b[idx],
					     c[idx],
					     d[idx],
					     e[idx],
					     g[idx],
					     h[idx],
					     i[idx]);
					
          }
}


/*
   I/O number of partial PCI writes per second.
*/
void
skx_PCI_part_writes_sec_samples(const double * __restrict __attribute__((aligned(64))) a,
			        const double * __restrict __attribute__((aligned(64))) b,
			        const double * __restrict __attribute__((aligned(64))) c,
			        double * __restrict __attribute__((aligned(64))) d,
	                        const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(i)					\
        aligned(a:64,b,c,d) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               d[i] = skx_PCI_part_writes_sec(a[i],
	                                      b[i],
					      c[i]);
					       
				
					
          }
}


/*
    I/O write cache miss(disk/network reads) bandwidth (MB/sec)
*/
void
skx_IO_writes_cache_miss_samples(const double * __restrict __attribute__((aligned(64))) a,
			         const double * __restrict __attribute__((aligned(64))) b,
			         const double * __restrict __attribute__((aligned(64))) c,
			         double * __restrict __attribute__((aligned(64))) d,
	                         const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(i)					\
        aligned(a:64,b,c,d) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               d[i] = skx_IO_write_cache_miss_bw(a[i],
	                                      b[i],
					      c[i]);
					       
				
					
          }				 
}


/*
    I/O write cache miss(disk/network writes) bandwidth (MB/sec)
*/
void
skx_IO_reads_cache_miss_samples(const double * __restrict __attribute__((aligned(64))) a,
			         const double * __restrict __attribute__((aligned(64))) b,
			         const double * __restrict __attribute__((aligned(64))) c,
			         double * __restrict __attribute__((aligned(64))) d,
	                         const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(i)					\
        aligned(a:64,b,c,d) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               d[i] = skx_IO_read_cache_miss_bw(a[i],
	                                      b[i],
					      c[i]);
					       
				
					
          }				 
}


/*
     IO cache miss(disk/network) bandwidth (MB/sec)
*/
void
skx_IO_cache_miss_total_bw_samples(const double * __restrict __attribute__((aligned(64))) a,
			           const double * __restrict __attribute__((aligned(64))) b,
				   const double * __restrict __attribute__((aligned(64))) c,
				   const double * __restrict __attribute__((aligned(64))) d,
			       	   double * __restrict __attribute__((aligned(64))) e,
	                           const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,data_len) private(i)					\
        aligned(a:64,b,c,d,e) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               e[i] = skx_IO_cache_miss_total_bw(a[i],
	                                      b[i],
					      c[i],
					      d[i]);
					       
				
					
          }
}


/*
    MMIO reads per second.
*/
void
skx_MMIO_reads_sec_samples(const double * __restrict  __attribute__((aligned(64))) a,
			   const double * __restrict  __attribute__((aligned(64))) b,
			   double * __restrict  __attribute__((aligned(64))) c,
	                   const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_MMIO_reads_sec(a[i],
	                                 b[i]);
					    
					       
				
					
          }
}


/*
    MMIO writes per second.
*/
void
skx_MMIO_writes_sec_samples(const double * __restrict  __attribute__((aligned(64))) a,
			   const double * __restrict  __attribute__((aligned(64))) b,
			   double * __restrict  __attribute__((aligned(64))) c,
	                   const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_MMIO_writes_sec(a[i],
	                                 b[i]);
					    
					       
				
					
          }
}


/*
    Memory Page Empty vs. all requests
*/
void
skx_mem_page_empty_all_reqs_samples(const double * __restrict __attribute__((aligned(64))) a,
			            const double * __restrict __attribute__((aligned(64))) b,
				    const double * __restrict __attribute__((aligned(64))) c,
				    const double * __restrict __attribute__((aligned(64))) d,
			       	    double * __restrict __attribute__((aligned(64))) e,
	                            const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,data_len) private(i)					\
        aligned(a:64,b,c,d,e) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               e[i] = skx_mem_page_empty_all_reqs(a[i],
	                                      b[i],
					      c[i],
					      d[i]);
					       
				
					
          }
}


/*
    Memory Page Misses vs. all requests
*/
void
skx_mem_page_misses_all_req_samples( const double * __restrict __attribute__((aligned(64))) a,
				     const double * __restrict __attribute__((aligned(64))) b,
				     const double * __restrict __attribute__((aligned(64))) c,
			       	     double * __restrict __attribute__((aligned(64))) d,
	                             const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(i)					\
        aligned(a:64,b,c,d) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               d[i] = skx_mem_page_misses_all_reqs(a[i],
	                                      b[i],
					      c[i]);
					    
	}				       
				
}


/*
     Memory Page Hits vs. all requests
*/
void
skx_mem_page_hits_all_req_samples(   const double * __restrict __attribute__((aligned(64))) a,
				     const double * __restrict __attribute__((aligned(64))) b,
				     const double * __restrict __attribute__((aligned(64))) c,
			       	     double * __restrict __attribute__((aligned(64))) d,
	                             const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(i)					\
        aligned(a:64,b,c,d) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               d[i] = skx_mem_page_hits_all_reqs(a[i],
	                                      b[i],
					      c[i]);
					    
	}
}


/*
    Memory percentage  Cycles where all DRAM ranks are in PPD mode.
*/
void
skx_PPD_DRAM_cycles_samples( const double * __restrict  __attribute__((aligned(64))) a,
			     const double * __restrict  __attribute__((aligned(64))) b,
			     double * __restrict  __attribute__((aligned(64))) c,
	                     const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_DRAM_PPD_mode_cycles(a[i],
	                                       b[i]);
					     
					    
	}
}


/*
     Memory percentage Cycles all ranks in critical thermal throttle.
*/
void
skx_mem_cycles_thermal_throttle_samples(const double * __restrict __attribute__((aligned(64))) a,
			                const double * __restrict __attribute__((aligned(64))) b,
			                double * __restrict __attribute__((aligned(64))) c,
	                                const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_mem_cycles_thermal_throttled(a[i],
	                                               b[i]);
					     
					    
	}
}


/*
    Memory  Cycles Memory is in self refresh power mode
*/
void
skx_mem_cycles_self_refresh_samples(const double * __restrict  __attribute__((aligned(64))) a,
			            const double * __restrict  __attribute__((aligned(64))) b,
			            double * __restrict  __attribute__((aligned(64))) c,
	                            const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
  aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_mem_cycles_self_refresh(a[i],
	                                          b[i]);
					     
					    
	}
}


/*
   Uops delivered from decoded Icache (DSB).
*/
void
skx_DSB_uops_samples(const double * __restrict __attribute__((aligned(64))) a,
		     const double * __restrict __attribute__((aligned(64))) b,
		     const double * __restrict __attribute__((aligned(64))) c,
		     const double * __restrict __attribute__((aligned(64))) d,
		     double * __restrict __attribute__((aligned(64))) e,
	             const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,data_len) private(i)					\
        aligned(a:64,b,c,d,e) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               e[i] = skx_DSB_uops_delivered(a[i],
	                                     b[i],
					     c[i],
					     d[i]);
					    
	}
}


/*
   Uops delivered from MITE.
*/
void
skx_MITE_uops_samples(const double * __restrict __attribute__((aligned(64))) a,
		     const double * __restrict __attribute__((aligned(64))) b,
		     const double * __restrict __attribute__((aligned(64))) c,
		     const double * __restrict __attribute__((aligned(64))) d,
		     double * __restrict __attribute__((aligned(64))) e,
	             const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,data_len) private(i)					\
        aligned(a:64,b,c,d,e) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               e[i] = skx_MITE_uops_delivered(a[i],
	                                     b[i],
					     c[i],
					     d[i]);
					    
	}
}


/*
   Uops delivered from MS.
*/
void
skx_MS_uops_samples(const double * __restrict __attribute__((aligned(64))) a,
		     const double * __restrict __attribute__((aligned(64))) b,
		     const double * __restrict __attribute__((aligned(64))) c,
		     const double * __restrict __attribute__((aligned(64))) d,
		     double * __restrict __attribute__((aligned(64))) e,
	             const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,data_len) private(i)					\
        aligned(a:64,b,c,d,e) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               e[i] = skx_Ms_uops_delivered(a[i],
	                                     b[i],
					     c[i],
					     d[i]);
					    
	}
}


/*
   Uops delivered from LSD.
*/
void
skx_LSD_uops_samples(const double * __restrict __attribute__((aligned(64))) a,
		     const double * __restrict __attribute__((aligned(64))) b,
		     const double * __restrict __attribute__((aligned(64))) c,
		     const double * __restrict __attribute__((aligned(64))) d,
		     double * __restrict __attribute__((aligned(64))) e,
	             const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,data_len) private(i)					\
        aligned(a:64,b,c,d,e) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               e[i] = skx_LSD_uops_delivered(a[i],
	                                     b[i],
					     c[i],
					     d[i]);
					    
	}
}


/*
     FP scalar single-precision FP instructions retired per instruction.
*/
void
skx_fp32_scalar_retired_samples(const double * __restrict __attribute__((aligned(64))) a,
		                const double * __restrict __attribute__((aligned(64))) b,
		                double * __restrict __attribute__((aligned(64))) c,
	                        const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_fp32_scalar_retired(a[i],
	                                      b[i]);
					     
					    
	}
}


/*
     FP scalar double-precision FP instructions retired per instruction.
*/
void
skx_fp64_scalar_retired_samples(const double * __restrict __attribute__((aligned(64))) a,
		                const double * __restrict __attribute__((aligned(64))) b,
		                double * __restrict __attribute__((aligned(64))) c,
	                        const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_fp64_scalar_retired(a[i],
	                                      b[i]);
					     
					    
	}
}


/*
   FP 128-bit packed single-precision FP instructions retired per instruction.
*/
void
skx_fp32_vec128b_retired_samples(const double * __restrict __attribute__((aligned(64))) a,
		                const double * __restrict __attribute__((aligned(64))) b,
		                double * __restrict __attribute__((aligned(64))) c,
	                        const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_fp32_vec128b_retired(a[i],
	                                      b[i]);
					     
					    
	}
}


/*
   FP 128-bit packed double-precision FP instructions retired per instruction.
*/
void
skx_fp64_vec128b_retired_samples(const double * __restrict __attribute__((aligned(64))) a,
		                const double * __restrict __attribute__((aligned(64))) b,
		                double * __restrict __attribute__((aligned(64))) c,
	                        const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_fp64_vec128b_retired(a[i],
	                                      b[i]);
					     
					    
	}
}


/*
     FP 256-bit packed single-precision FP instructions retired per instruction.
*/
void
skx_fp32_vec256b_retired_samples(const double * __restrict __attribute__((aligned(64))) a,
		                const double * __restrict __attribute__((aligned(64))) b,
		                double * __restrict __attribute__((aligned(64))) c,
	                        const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_fp32_vec256b_retired(a[i],
	                                      b[i]);
					     
					    
	}
}


/*
     FP 256-bit packed single-precision FP instructions retired per instruction.
*/
void
skx_fp64_vec256b_retired_samples(const double * __restrict __attribute__((aligned(64))) a,
		                const double * __restrict __attribute__((aligned(64))) b,
		                double * __restrict __attribute__((aligned(64))) c,
	                        const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_fp64_vec256b_retired(a[i],
	                                      b[i]);
					     
					    
	}
}


/*
   FP 512-bit packed single-precision FP instructions retired per instruction.
*/
void
skx_fp32_vec512b_retired_samples(const double * __restrict __attribute__((aligned(64))) a,
		                const double * __restrict __attribute__((aligned(64))) b,
		                double * __restrict __attribute__((aligned(64))) c,
	                        const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_fp32_vec512b_retired(a[i],
	                                      b[i]);
					     
					    
	}
}


/*
   FP 512-bit packed double-precision FP instructions retired per instruction.
*/
void
skx_fp64_vec512b_retired_samples(const double * __restrict __attribute__((aligned(64))) a,
		                const double * __restrict __attribute__((aligned(64))) b,
		                double * __restrict __attribute__((aligned(64))) c,
	                        const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_fp64_vec512b_retired(a[i],
	                                      b[i]);
					     
					    
	}
}


/*
    FP instruction density (percentage).
*/
void
skx_fp_instr_density_samples(const double * __restrict  __attribute__((aligned(64))) a,
		             const double * __restrict  __attribute__((aligned(64))) b,
			     const double * __restrict  __attribute__((aligned(64))) c,
			     const double * __restrict  __attribute__((aligned(64))) d,
			     const double * __restrict  __attribute__((aligned(64))) e,
			     const double * __restrict  __attribute__((aligned(64))) f,
			     const double * __restrict  __attribute__((aligned(64))) g,
			     const double * __restrict  __attribute__((aligned(64))) h,
			     const double * __restrict  __attribute__((aligned(64))) i,
		             double * __restrict  __attribute__((aligned(64))) j,
	                     const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,f,g,h,i,j,data_len) private(idx)					\
        aligned(a:64,b,c,d,e,f,g,h,i,j) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               j[idx] = skx_fp_instructions_density(a[idx],
	                                     b[idx],
					     c[idx],
					     d[idx],
					     e[idx],
					     g[idx],
					     h[idx],
					     i[idx]);
					
          }
}


/*
   Branch instructions density.
*/
void
skx_branch_instr_ratio_samples(const double * __restrict __attribute__((aligned(64))) a,
		               const double * __restrict __attribute__((aligned(64))) b,
		               double * __restrict __attribute__((aligned(64))) c,
	                        const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_branch_instr_density(a[i],
	                                      b[i]);
					     
					    
	}
}


/*
    DRAM power (watts).
*/
void
skx_DRAM_power_samples(const double * __restrict __attribute__((aligned(64))) a,
		       const double * __restrict __attribute__((aligned(64))) b,
		       double * __restrict __attribute__((aligned(64))) c,
	               const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_DRAM_power(a[i],
	                             b[i]);
					     
					    
	}
}


/*
   Package power (watts)
*/
void
skx_Package_power_samples(const double * __restrict __attribute__((aligned(64))) a,
		       const double * __restrict __attribute__((aligned(64))) b,
		       double * __restrict __attribute__((aligned(64))) c,
	               const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_Package_power(a[i],
	                             b[i]);
					     
					    
	}
}


/*
   Core c3 residency
*/
void
skx_Core_C3_residency_samples(const double * __restrict __attribute__((aligned(64))) a,
		       const double * __restrict __attribute__((aligned(64))) b,
		       double * __restrict __attribute__((aligned(64))) c,
	               const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_Core_C3_residency(a[i],
	                             b[i]);
					     
					    
	}
}


/*
   Core c6 residency
*/
void
skx_Core_C6_residency_samples(const double * __restrict __attribute__((aligned(64))) a,
		       const double * __restrict __attribute__((aligned(64))) b,
		       double * __restrict __attribute__((aligned(64))) c,
	               const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_Core_C6_residency(a[i],
	                             b[i]);
					     
					    
	}
}


/*
    Package C2 residency
*/
void
skx_Package_C2_residency_samples(const double * __restrict __attribute__((aligned(64))) a,
		       const double * __restrict __attribute__((aligned(64))) b,
		       double * __restrict __attribute__((aligned(64))) c,
	               const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_Package_C2_residency(a[i],
	                             b[i]);
					     
					    
	}
}


/*
    Package C3 residency
*/
void
skx_Package_C3_residency_samples(const double * __restrict __attribute__((aligned(64))) a,
		       const double * __restrict __attribute__((aligned(64))) b,
		       double * __restrict __attribute__((aligned(64))) c,
	               const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_Package_C3_residency(a[i],
	                             b[i]);
					     
					    
	}
}


/*
    Package C6 residency
*/
void
skx_Package_C6_residency_samples(const double * __restrict __attribute__((aligned(64))) a,
		       const double * __restrict __attribute__((aligned(64))) b,
		       double * __restrict __attribute__((aligned(64))) c,
	               const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_Package_C6_residency(a[i],
	                             b[i]);
					     
					    
	}
}


/*
     core SW prefetch NTA per instruction.
*/
void
skx_NTA_sw_prefetches_samples(const double * __restrict __attribute__((aligned(64))) a,
		              const double * __restrict __attribute__((aligned(64))) b,
		              double * __restrict __attribute__((aligned(64))) c,
	                      const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_sw_nta_prefetches(a[i],
	                                    b[i]);
					     
					    
	}
}


/*
    Core cycles power throttled
*/
void
skx_Core_power_throttled_samples(const double * __restrict  __attribute__((aligned(64))) a,
				 const double * __restrict   __attribute__((aligned(64))) b,
		                 double * __restrict  __attribute__((aligned(64))) c,
	                         const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_Core_power_throttled(a[i],
	                                       b[i]);
					     
					    
	}
}


/*
   Core IPC
*/
void
skx_Core_IPC_samples(const double * __restrict  __attribute__((aligned(64))) a,
		     const double * __restrict  __attribute__((aligned(64))) b,
		     const double * __restrict  __attribute__((aligned(64))) c,
		     double * __restrict  __attribute__((aligned(64))) d,
	             const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(i)					\
        aligned(a:64,b,c,d) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               d[i] = skx_TMAM_core_ipc(a[i],
	                                b[i],
					c[i]);
					     
					    
	}
}


/*
    Info Memory Level Parallelism
*/
void
skx_mem_lvl_parallelism_samples(const double * __restrict  __attribute__((aligned(64))) a,
		                const double * __restrict  __attribute__((aligned(64))) b,
		                const double * __restrict  __attribute__((aligned(64))) c,
		                double * __restrict  __attribute__((aligned(64))) d,
	                        const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(i)					\
        aligned(a:64,b,c,d) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               d[i] = skx_TMAM_mem_level_parallelism(a[i],
	                                             b[i],
					             c[i]);
					     
					    
	}
}


/*
    Info cycles both threads active
*/
void
skx_SMT_activity_samples(const double * __restrict __attribute__((aligned(64))) a,
		         const double * __restrict __attribute__((aligned(64))) b,
		         const double * __restrict __attribute__((aligned(64))) c,
		         double * __restrict __attribute__((aligned(64))) d,
	                 const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(i)					\
        aligned(a:64,b,c,d) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               d[i] = skx_TMAM_SMT_activity(a[i],
	                                    b[i],
					    c[i]);
					     
					    
	}
}


/*
   Frontend bound
*/
void
skx_FrontEnd_bound_samples(const double * __restrict __attribute__((aligned(64))) a,
		           const double * __restrict __attribute__((aligned(64))) b,
		           double * __restrict __attribute__((aligned(64))) c,
	                   const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_frontend_bound(a[i],
	                                 b[i]);
					   
					     
					    
	}
}


/*
    Frontend latency
*/
void
skx_FrontEnd_latency_samples(const double * __restrict __attribute__((aligned(64))) a,
		             const double * __restrict __attribute__((aligned(64))) b,
		             double * __restrict __attribute__((aligned(64))) c,
	                     const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_frontend_bound(a[i],
	                                 b[i]);
					   
					     
					    
	}
}


/*
   ICache misses
*/
void
skx_ICache_misses_samples(const double * __restrict __attribute__((aligned(64))) a,
		          const double * __restrict __attribute__((aligned(64))) b,
			  const double * __restrict __attribute__((aligned(64))) c,
		          double * __restrict __attribute__((aligned(64))) d,
	                  const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(i)					\
        aligned(a:64,b,c,d) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               d[i] = skx_ICache_misses(a[i],
	                                b[i],
					c[i]);
					     
					    
	}
}


/*
   ITLB misses
*/
void
skx_ITLB_misses_samples(const double * __restrict __attribute__((aligned(64))) a,
		          const double * __restrict __attribute__((aligned(64))) b,
			  const double * __restrict __attribute__((aligned(64))) c,
		          double * __restrict __attribute__((aligned(64))) d,
	                  const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(i)					\
        aligned(a:64,b,c,d) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               d[i] = skx_ITLB_misses(a[i],
	                                b[i],
					c[i]);
					     
					    
	}
}


/*
   Branch resteers
*/
void
skx_branch_resteers_samples(const double * __restrict __attribute__((aligned(64))) a,
		            const double * __restrict __attribute__((aligned(64))) b,
			    const double * __restrict __attribute__((aligned(64))) c,
			    const double * __restrict __attribute__((aligned(64))) d,
			    const double * __restrict __attribute__((aligned(64))) e,
			    const double * __restrict __attribute__((aligned(64))) f,
			    const double * __restrict __attribute__((aligned(64))) g,
		            double * __restrict __attribute__((aligned(64))) h,
	                    const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,f,g,h,data_len) private(idx)					\
        aligned(a:64,b,c,d,e,f,g,h) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               h[idx] = skx_branch_resteers(a[idx],
	                                    b[idx],
					    c[idx],
					    d[idx],
					    e[idx],
					    f[idx],
					    g[idx]);
					         
					
          }
}


/*
   DSB switches
*/
void
skx_DSB_switches_samples( const double * __restrict  __attribute__((aligned(64))) a,
			  const double * __restrict  __attribute__((aligned(64))) b,
		          double * __restrict  __attribute__((aligned(64))) c,
	                  const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_DSB_switches(a[i],
	                               b[i]);
				
					     
					    
	}
}


/*
    MS switches
*/
void
skx_MS_switches_samples( const double * __restrict  __attribute__((aligned(64))) a,
			  const double * __restrict  __attribute__((aligned(64))) b,
		          double * __restrict  __attribute__((aligned(64))) c,
	                  const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_MS_switches(a[i],
	                               b[i]);
				
					     
					    
	}
}


/*
   Frontend bandwidth (SMT enabled)
*/
void
skx_FrontEnd_smt_bw_samples(const double * __restrict  __attribute__((aligned(64))) a,
			    const double * __restrict  __attribute__((aligned(64))) b,
			    const double * __restrict  __attribute__((aligned(64))) c,
			    const double * __restrict  __attribute__((aligned(64))) d,
		            double * __restrict  __attribute__((aligned(64))) e,
	                    const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,data_len) private(i)					\
        aligned(a:64,b,c,d,e) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               e[i] = skx_smt_frontend_bw(a[i],
	                                    b[i],
					    c[i],
					    d[i]);
					     
					    
	}
}


/*
    Frontend bandwidth (SMT disabled).
*/
void
skx_FrontEnd_n_smt_bw_samples(const double * __restrict  __attribute__((aligned(64))) a,
			      const double * __restrict  __attribute__((aligned(64))) b,
			      const double * __restrict  __attribute__((aligned(64))) c,
		              double * __restrict  __attribute__((aligned(64))) d,
	                      const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,data_len) private(i)					\
        aligned(a:64,b,c,d,e) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               d[i] = skx_smt_frontend_bw(a[i],
	                                    b[i],
					    c[i]);
					   
					     
					    
	}
}


/*
    Bad Speculation (SMT enabled).
*/
void
skx_smt_bad_speculate_samples(  const double * __restrict __attribute__((aligned(64))) a,
			        const double * __restrict __attribute__((aligned(64))) b,
			        const double * __restrict __attribute__((aligned(64))) c,
			        const double * __restrict __attribute__((aligned(64))) d,
				const double * __restrict __attribute__((aligned(64))) e,
		                double * __restrict __attribute__((aligned(64))) f,
	                        const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,f,data_len) private(i)					\
        aligned(a:64,b,c,d,e,f) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               f[i] = skx_smt_bad_speculation(a[i],
	                                        b[i],
					        c[i],
					        d[i],
					        e[i]);
					  
					         
					
          }
}


/*
    Bad Speculation (SMT disabled).
*/
void
skx_no_smt_bad_speculate_samples(  const double * __restrict __attribute__((aligned(64))) a,
			           const double * __restrict __attribute__((aligned(64))) b,
			           const double * __restrict __attribute__((aligned(64))) c,
			           const double * __restrict __attribute__((aligned(64))) d,
				   double * __restrict __attribute__((aligned(64))) e,
	                           const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,data_len) private(i)					\
        aligned(a:64,b,c,d,e) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               e[i] = skx_no_smt_bad_speculation(a[i],
	                                        b[i],
					        c[i],
					        d[i]);
					    
					  
					         
					
          }
}


/*
    Branch Mispredicts (STM enabled).
*/
void
skx_smt_branch_mispredicts_sample(const double * __restrict __attribute__((aligned(64))) a,
			          const double * __restrict __attribute__((aligned(64))) b,
			          const double * __restrict __attribute__((aligned(64))) c,
			          const double * __restrict __attribute__((aligned(64))) d,
				  const double * __restrict __attribute__((aligned(64))) e,
				  const double * __restrict __attribute__((aligned(64))) f,
				  const double * __restrict __attribute__((aligned(64))) g,
			          double * __restrict __attribute__((aligned(64))) h,
	                          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,f,g,h,data_len) private(idx)					\
        aligned(a:64,b,c,d,e,f,g,h) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               h[idx] = skx_smt_branch_mispredicts(a[idx],
	                                           b[idx],
					           c[idx],
					           d[idx],
					           e[idx],
					           f[idx],
					           g[idx]);
					         
					
          }
}


/*
     Branch Misprediction (SMT disabled).
*/
void
skx_no_smt_branch_mispredicts_sample(const double * __restrict  __attribute__((aligned(64))) a,
			          const double * __restrict __attribute__((aligned(64))) b,
			          const double * __restrict __attribute__((aligned(64))) c,
			          const double * __restrict __attribute__((aligned(64))) d,
				  const double * __restrict __attribute__((aligned(64))) e,
				  const double * __restrict __attribute__((aligned(64))) f,
				  double * __restrict __attribute__((aligned(64))) g,
	                          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,f,g,data_len) private(idx)					\
        aligned(a:64,b,c,d,e,f,g) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               g[idx] = skx_no_smt_branch_mispredicts(a[idx],
	                                           b[idx],
					           c[idx],
					           d[idx],
					           e[idx],
					           f[idx]);
					         
					         
					
          }
}


/*
    Machine Clears (SMT enabled).
*/
void
skx_smt_machine_clears_samples(   const double * __restrict __attribute__((aligned(64))) a,
			          const double * __restrict __attribute__((aligned(64))) b,
			          const double * __restrict __attribute__((aligned(64))) c,
			          const double * __restrict __attribute__((aligned(64))) d,
				  const double * __restrict __attribute__((aligned(64))) e,
				  const double * __restrict __attribute__((aligned(64))) f,
				  const double * __restrict __attribute__((aligned(64))) g,
			          double * __restrict __attribute__((aligned(64))) h,
	                          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,f,g,h,data_len) private(idx)					\
        aligned(a:64,b,c,d,e,f,g,h) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               h[idx] = skx_smt_machine_clears(    a[idx],
	                                           b[idx],
					           c[idx],
					           d[idx],
					           e[idx],
					           f[idx],
					           g[idx]);
					         
					
          }
}


/*
    Machine clears (SMT disabled).
*/
void
skx_smt_machine_clears_samples(   const double * __restrict __attribute__((aligned(64))) a,
			          const double * __restrict __attribute__((aligned(64))) b,
			          const double * __restrict __attribute__((aligned(64))) c,
			          const double * __restrict __attribute__((aligned(64))) d,
				  const double * __restrict __attribute__((aligned(64))) e,
				  const double * __restrict __attribute__((aligned(64))) f,
				  const double * __restrict __attribute__((aligned(64))) g,
			          double * __restrict __attribute__((aligned(64))) h,
	                          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,f,g,h,data_len) private(idx)					\
        aligned(a:64,b,c,d,e,f,g,h) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               h[idx] = skx_smt_machine_clears(    a[idx],
	                                           b[idx],
					           c[idx],
					           d[idx],
					           e[idx],
					           f[idx],
					           g[idx]);
					         
					
          }
}


/*
    Machine clears (SMT disabled).
*/
void
skx_smt_machine_clears_samples(   const double * __restrict __attribute__((aligned(64))) a,
			          const double * __restrict __attribute__((aligned(64))) b,
			          const double * __restrict __attribute__((aligned(64))) c,
			          const double * __restrict __attribute__((aligned(64))) d,
				  const double * __restrict __attribute__((aligned(64))) e,
				  const double * __restrict __attribute__((aligned(64))) f,
				  const double * __restrict __attribute__((aligned(64))) g,
			          double * __restrict __attribute__((aligned(64))) h,
	                          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,f,g,h,data_len) private(idx)					\
        aligned(a:64,b,c,d,e,f,g,h) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               h[idx] = skx_smt_machine_clears(    a[idx],
	                                           b[idx],
					           c[idx],
					           d[idx],
					           e[idx],
					           f[idx],
					           g[idx]);
					         
					
          }
}


/*
    Machine clears (SMT disabled).
*/
void
skx_no_smt_machine_clears_samples(   const double * __restrict __attribute__((aligned(64))) a,
			          const double * __restrict __attribute__((aligned(64))) b,
			          const double * __restrict __attribute__((aligned(64))) c,
			          const double * __restrict __attribute__((aligned(64))) d,
				  const double * __restrict __attribute__((aligned(64))) e,
				  const double * __restrict __attribute__((aligned(64))) f,
				  double * __restrict __attribute__((aligned(64))) g,
	                          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,f,g,data_len) private(idx)					\
        aligned(a:64,b,c,d,e,f,g) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               g[idx] = skx_no_smt_machine_clears(    a[idx],
	                                           b[idx],
					           c[idx],
					           d[idx],
					           e[idx],
					           f[idx]);
					         
					         
					
          }
}


/*
    Backend Bound (SMT enabled)
*/
void
skx_smt_BackEnd_bound_samples(    const double * __restrict  __attribute__((aligned(64))) a,
			          const double * __restrict  __attribute__((aligned(64))) b,
			          const double * __restrict  __attribute__((aligned(64))) c,
			          const double * __restrict  __attribute__((aligned(64))) d,
				  const double * __restrict  __attribute__((aligned(64))) e,
				  const double * __restrict  __attribute__((aligned(64))) f,
				  double * __restrict  __attribute__((aligned(64))) g,
	                          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,f,g,data_len) private(idx)					\
        aligned(a:64,b,c,d,e,f,g) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               g[idx] = skx_smt_backend_bound(     a[idx],
	                                           b[idx],
					           c[idx],
					           d[idx],
					           e[idx],
					           f[idx]);
					         
					         
					
          }
}


/*
     Backend Boud (SMT disabled).
*/
void
skx_smt_BackEnd_bound_samples(    const double * __restrict  __attribute__((aligned(64))) a,
			          const double * __restrict  __attribute__((aligned(64))) b,
			          const double * __restrict  __attribute__((aligned(64))) c,
			          const double * __restrict  __attribute__((aligned(64))) d,
				  const double * __restrict  __attribute__((aligned(64))) e,
				  double * __restrict  __attribute__((aligned(64))) f,
	                          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,f,data_len) private(idx)					\
        aligned(a:64,b,c,d,e,f) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               f[idx] = skx_smt_backend_bound(     a[idx],
	                                           b[idx],
					           c[idx],
					           d[idx],
					           e[idx]);
					        
					         
					         
					
          }
}


/*
   L1D bound
*/
void
skx_L1D_bound_samples(const double * __restrict  __attribute__((aligned(64))) a,
		      const double * __restrict  __attribute__((aligned(64))) b,
		      const double * __restrict  __attribute__((aligned(64))) c,
		      double * __restrict  __attribute__((aligned(64))) d,
	              const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(i)					\
        aligned(a:64,b,c,d) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               d[i] = skx_L1D_bound(a[i],
	                            b[i],
				    c[i]);
					   
					     
					    
	}
}


/*
   DTLB Load
*/
void
skx_DTLB_load_samples(const double * __restrict __attribute__((aligned(64))) a,
		      const double * __restrict __attribute__((aligned(64))) b,
		      const double * __restrict __attribute__((aligned(64))) c,
		      double * __restrict __attribute__((aligned(64))) d,
	              const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(i)					\
        aligned(a:64,b,c,d) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               d[i] = skx_DTLB_load(a[i],
	                            b[i],
				    c[i]);
					   
					     
					    
	}
}


/*
   Stores forward blocked
*/
void
skx_stores_fwd_block_samples( const double * __restrict __attribute__((aligned(64))) a,
		              const double * __restrict __attribute__((aligned(64))) b,
		              double * __restrict __attribute__((aligned(64))) c,
	                      const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(i)					\
        aligned(a:64,b,c) linear(i:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t i = 0; i != data_len; ++i) {
               c[i] = skx_stores_fwd_blocked(a[i],
	                                     b[i]);
				 
	 }				   
					     

}


/*
   Lock latency
*/
void
skx_lock_latency_samples(const double * __restrict __attribute__((aligned(64))) a,
		         const double * __restrict __attribute__((aligned(64))) b,
		         const double * __restrict __attribute__((aligned(64))) c,
			 const double * __restrict __attribute__((aligned(64))) d,
		         double * __restrict __attribute__((aligned(64))) e,
	                 const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,data_len) private(idx)					\
        aligned(a:64,b,c,d,e) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               e[idx] = skx_lock_latency(          a[idx],
	                                           b[idx],
					           c[idx],
					           d[idx]);
					        
					        
	}				         
					         
}


/*
   L2 bound
*/
void
skx_L2_bound_samples(    const double * __restrict __attribute__((aligned(64))) a,
		         const double * __restrict __attribute__((aligned(64))) b,
			 const double * __restrict __attribute__((aligned(64))) c,
		         double * __restrict __attribute__((aligned(64))) d,
	                 const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(idx)					\
        aligned(a:64,b,c,d) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               d[idx] = skx_L2_bound(  a[idx],
	                               b[idx],
				       c[idx]);
					       
					        
					        
	}	
}


/*
    L3 bound.
*/
void
skx_L3_bound_samples(    const double * __restrict __attribute__((aligned(64))) a,
		         const double * __restrict __attribute__((aligned(64))) b,
			 const double * __restrict __attribute__((aligned(64))) c,
		         double * __restrict __attribute__((aligned(64))) d,
	                 const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(idx)					\
        aligned(a:64,b,c,d) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               d[idx] = skx_L3_bound(  a[idx],
	                               b[idx],
				       c[idx]);
					       
					        
					        
	}
}


/*
    Contested accessses
*/
void
skx_contested_accesses_samples( const double * __restrict __attribute__((aligned(64))) a,
		                const double * __restrict __attribute__((aligned(64))) b,
			        const double * __restrict __attribute__((aligned(64))) c,
		                double * __restrict __attribute__((aligned(64))) d,
	                        const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(idx)					\
        aligned(a:64,b,c,d) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               d[idx] = skx_contested_accesses(  a[idx],
	                                         b[idx],
				                 c[idx]);
					       
					        
					        
	}
}


/*
    Data sharing
*/
void
skx_data_sharing_samples( const double * __restrict  __attribute__((aligned(64))) a,
			  const double * __restrict  __attribute__((aligned(64))) b,
		          double * __restrict  __attribute__((aligned(64))) c,
	                  const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(idx)					\
        aligned(a:64,b,c) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               c[idx] = skx_data_sharing(  a[idx],
	                                   b[idx]);
				               
	}				       
					        
}


/*
    L3 latency
*/
void
skx_L3_latency_samples(const double * __restrict __attribute__((aligned(64))) a,
		       const double * __restrict __attribute__((aligned(64))) b,
		       const double * __restrict __attribute__((aligned(64))) c,
		       const double * __restrict __attribute__((aligned(64))) d,
		       const double * __restrict __attribute__((aligned(64))) e,
		       double * __restrict __attribute__((aligned(64))) f,
	               const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,f,data_len) private(idx)					\
        aligned(a:64,b,c,d,e,f) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               f[idx] = skx_L3_latency(           a[idx],
	                                           b[idx],
					           c[idx],
					           d[idx],
						   e[idx]);
					        
					        
	}
}


/*
  L3 bandwidth
*/
void
skx_L3_bw_samples(     const double * __restrict  __attribute__((aligned(64))) a,
		       const double * __restrict  __attribute__((aligned(64))) b,
		       const double * __restrict  __attribute__((aligned(64))) c,
		       double * __restrict  __attribute__((aligned(64))) d,
	               const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(idx)					\
        aligned(a:64,b,c,d) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               d[idx] = skx_L3_bandwidth(  a[idx],
	                                   b[idx],
				           c[idx]);
					       
					        
					        
	}
}


/*
   SuperQueue (SMT enabled)
*/
void
skx_smt_SQ_full_samples(const double * __restrict __attribute__((aligned(64))) a,
		        const double * __restrict __attribute__((aligned(64))) b,
		        const double * __restrict __attribute__((aligned(64))) c,
		        double * __restrict __attribute__((aligned(64))) d,
	                const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(idx)					\
        aligned(a:64,b,c,d) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               d[idx] = skx_smt_SQ_full(  a[idx],
	                                   b[idx],
				           c[idx]);
					       
					        
					        
	}
}



/*
   SuperQueue (SMT disabled)
*/
void
skx_no_smt_SQ_full_samples(const double * __restrict __attribute__((aligned(64))) a,
		           const double * __restrict __attribute__((aligned(64))) b,
		           double * __restrict __attribute__((aligned(64))) c,
		      	   const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(idx)					\
        aligned(a:64,b,c) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               c[idx] = skx_smt_SQ_full(  a[idx],
	                                   b[idx]);
				      
					       
					        
					        
	}
}


/*
   Memory bound.
*/
void
skx_memory_bound_samples(  const double * __restrict __attribute__((aligned(64))) a,
		           const double * __restrict __attribute__((aligned(64))) b,
		       	   double * __restrict __attribute__((aligned(64))) c,
	                   const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(idx)					\
        aligned(a:64,b,c) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               c[idx] = skx_memory_bound(  a[idx],
	                                   b[idx]);
				      
					       
					        
					        
	}
}


/*
   Memory BW.
*/
void
skx_memory_bw_samples(  const double * __restrict __attribute__((aligned(64))) a,
		           const double * __restrict __attribute__((aligned(64))) b,
		       	   double * __restrict __attribute__((aligned(64))) c,
	                   const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(idx)					\
        aligned(a:64,b,c) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               c[idx] = skx_memory_bandwidth(  a[idx],
	                                   b[idx]);
				      
					       
					        
					        
	}
}


/*
   Memory latency
*/
void
skx_memory_latency_samples(const double * __restrict  __attribute__((aligned(64))) a,
		           const double * __restrict  __attribute__((aligned(64))) b,
			   const double * __restrict  __attribute__((aligned(64))) c,
		       	   double * __restrict  __attribute__((aligned(64))) d,
	                   const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(idx)					\
        aligned(a:64,b,c,d) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               d[idx] = skx_memory_latency(  a[idx],
	                                   b[idx],
				           c[idx]);
					       
					        
					        
	}
}


/*
   Stores bound
*/
void
skx_stores_bound_samples(  const double * __restrict  __attribute__((aligned(64))) a,
			   const double * __restrict   __attribute__((aligned(64))) b,
		       	   double * __restrict  __attribute__((aligned(64))) c,
	                   const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(idx)					\
        aligned(a:64,b,c) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               c[idx] = skx_stores_bound(  a[idx],
	                                   b[idx]);
				      
					       
					        
					        
	}
}


/*
   DTLB stores (SMT enabled)
*/
void
skx_smt_DTLB_stores_samples(const double * __restrict __attribute__((aligned(64))) a,
		            const double * __restrict __attribute__((aligned(64))) b,
			    const double * __restrict __attribute__((aligned(64))) c,
			    const double * __restrict __attribute__((aligned(64))) d,
		       	    double * __restrict __attribute__((aligned(64))) e,
	                    const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,data_len) private(idx)					\
        aligned(a:64,b,c,d,e,) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               e[idx] = skx_smt_DTLB_stores(       a[idx],
	                                           b[idx],
					           c[idx],
					           d[idx]);
						 
					        
					        
	}
}


/*
   DTLB stores (SMT enabled)
*/
void
skx_no_smt_DTLB_stores_samples(const double * __restrict __attribute__((aligned(64))) a,
		            const double * __restrict __attribute__((aligned(64))) b,
			    const double * __restrict __attribute__((aligned(64))) c,
			   double * __restrict __attribute__((aligned(64))) d,
	                    const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(idx)					\
        aligned(a:64,b,c,d,) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               d[idx] = skx_no_smt_DTLB_stores(       a[idx],
	                                           b[idx],
					           c[idx]);
					        
						 
					        
					        
	}
}


/*
   Divider
*/
void
skx_divider_samples(const double * __restrict  __attribute__((aligned(64))) a,
		    const double * __restrict  __attribute__((aligned(64))) b,
		    double * __restrict  __attribute__((aligned(64))) b,
	            const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(idx)					\
        aligned(a:64,b,c) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               c[idx] = skx_Divider(  a[idx],
	                                   b[idx]);
				      
					       
					        
					        
	}
}


/*
   Ports Utilization
*/
void
skx_Ports_utilization_samples(const double * __restrict  __attribute__((aligned(64))) a,
		              const double * __restrict  __attribute__((aligned(64))) b,
			      const double * __restrict  __attribute__((aligned(64))) c,
			      const double * __restrict  __attribute__((aligned(64))) d,
			      const double * __restrict  __attribute__((aligned(64))) e,
			      const double * __restrict  __attribute__((aligned(64))) f,
			      const double * __restrict  __attribute__((aligned(64))) g,
			      double * __restrict  __attribute__((aligned(64))) h,
	                      const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,f,g,h,data_len) private(idx)					\
        aligned(a:64,b,c,d,e,f,g,h) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               h[idx] = skx_Ports_utlization(      a[idx],
	                                           b[idx],
					           c[idx],
					           d[idx],
						   e[idx],
						   f[idx],
						   g[idx]);
					        
					        
	}
}


/*
    Port_0 Utilized.
*/
void
skx_Port_0_utlized_samples(   const double * __restrict __attribute__((aligned(64))) a,
			      const double * __restrict __attribute__((aligned(64))) b,
			      const double * __restrict __attribute__((aligned(64))) c,
			      const double * __restrict __attribute__((aligned(64))) d,
			      double * __restrict __attribute__((aligned(64))) e,
	                      const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,data_len) private(idx)					\
        aligned(a:64,b,c,d,e) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               e[idx] = skx_Port_0_utilized(       a[idx],
	                                           b[idx],
					           c[idx],
					           d[idx]);
						 
					        
	     
	}
}


/*
   Port_1 utlized
*/
void
skx_Port_1_utlized_samples(   const double * __restrict __attribute__((aligned(64))) a,
			      const double * __restrict __attribute__((aligned(64))) b,
			      const double * __restrict __attribute__((aligned(64))) c,
			      const double * __restrict __attribute__((aligned(64))) d,
			      const double * __restrict __attribute__((aligned(64))) e,
			      double * __restrict __attribute__((aligned(64))) f,
	                      const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,f,data_len) private(idx)					\
        aligned(a:64,b,c,d,e,f) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               f[idx] = skx_Port_1_utilized(       a[idx],
	                                           b[idx],
					           c[idx],
					           d[idx],
						   e[idx]);
						 
					        
	     
	}
}


/*
    Port_2 utilized
*/
void
skx_Port_2_utlized_samples(   const double * __restrict __attribute__((aligned(64))) a,
			      const double * __restrict __attribute__((aligned(64))) b,
			      const double * __restrict __attribute__((aligned(64))) c,
			      const double * __restrict __attribute__((aligned(64))) d,
			      const double * __restrict __attribute__((aligned(64))) e,
			      double * __restrict __attribute__((aligned(64))) f,
	                      const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,f,data_len) private(idx)					\
        aligned(a:64,b,c,d,e,f) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               f[idx] = skx_Port_2_utilized(       a[idx],
	                                           b[idx],
					           c[idx],
					           d[idx],
						   e[idx]);
						 
					        
	     
	}
}


/*
    Port_3m utilized
*/
void
skx_Port_3m_utlized_samples(  const double * __restrict __attribute__((aligned(64))) a,
			      const double * __restrict __attribute__((aligned(64))) b,
			      const double * __restrict __attribute__((aligned(64))) c,
			      double * __restrict __attribute__((aligned(64))) d,
	                      const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(idx)					\
        aligned(a:64,b,c,d) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               d[idx] = skx_Port_3m_utilized(       a[idx],
	                                           b[idx],
					           c[idx]);
					        
						 
					        
					        
	}
}


/*
   Retiring (SMT enabled)
*/
void
skx_smt_retiring_samples(     const double * __restrict __attribute__((aligned(64))) a,
			      const double * __restrict __attribute__((aligned(64))) b,
			      const double * __restrict __attribute__((aligned(64))) c,
			      double * __restrict __attribute__((aligned(64))) d,
	                      const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(idx)					\
        aligned(a:64,b,c,d) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               d[idx] = skx_smt_retiring(          a[idx],
	                                           b[idx],
					           c[idx]);
					        
						 
					        
					        
	}
}


/*
   Retiring (SMT disabled)
*/
void
skx_no_smt_retiring_samples(  const double * __restrict __attribute__((aligned(64))) a,
			      const double * __restrict __attribute__((aligned(64))) b,
			      double * __restrict __attribute__((aligned(64))) c,
			      const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(idx)					\
        aligned(a:64,b,c) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               c[idx] = skx_no_smt_retiring(          a[idx],
	                                           b[idx]);
					          
					        
						 
					        
					        
	}
}


/*
    Basic activity (SMT enabled).
*/
void
skx_smt_basic_activity_samples( const double * __restrict  __attribute__((aligned(64))) a,
			        const double * __restrict  __attribute__((aligned(64))) b,
			        const double * __restrict  __attribute__((aligned(64))) c,
				const double * __restrict  __attribute__((aligned(64))) d,
				const double * __restrict  __attribute__((aligned(64))) f,
			        double * __restrict  __attribute__((aligned(64))) g,
	                        const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,f,g,data_len) private(idx)					\
        aligned(a:64,b,c,d,e,f,g) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               g[idx] = skx_smt_base_activity(     a[idx],
	                                           b[idx],
					           c[idx],
					           d[idx],
						   e[idx],
						   f[idx]);
						 
					        
	     
	}
}


/*
    Basic activity (SMT disabled).
*/
void
skx_smt_basic_activity_samples( const double * __restrict  __attribute__((aligned(64))) a,
			        const double * __restrict  __attribute__((aligned(64))) b,
			        const double * __restrict  __attribute__((aligned(64))) c,
				const double * __restrict  __attribute__((aligned(64))) d,
				const double * __restrict  __attribute__((aligned(64))) e,
				double * __restrict __attribute__((aligned(64))) f,
			        const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,f,data_len) private(idx)					\
        aligned(a:64,b,c,d,e,f) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               f[idx] = skx_smt_base_activity(     a[idx],
	                                           b[idx],
					           c[idx],
					           d[idx],
						   e[idx]);
						 
						 
					        
	     
	}
}


/*
      FP scalar retiring fraction.
*/
void
skx_fp_scalar_fract_samples(           const double * __restrict  __attribute__((aligned(64))) a,
			               const double * __restrict  __attribute__((aligned(64))) b,
				       const double * __restrict  __attribute__((aligned(64))) c,
			               double * __restrict __attribute__((aligned(64))) d,
	                               const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,data_len) private(idx)					\
  aligned(a:64,b,c,d) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               d[idx] = skx_fp_scalar_fraction(    a[idx],
	                                           b[idx],
					           c[idx]);
					        
						 
					        
					        
	}
}


/*
     FP vector retiring fraction.
*/
void
skx_fp_vector_fract_samples(    const double * __restrict __attribute__((aligned(64))) a,
			        const double * __restrict __attribute__((aligned(64))) b,
			        const double * __restrict __attribute__((aligned(64))) c,
				const double * __restrict __attribute__((aligned(64))) d,
				const double * __restrict __attribute__((aligned(64))) e,
				const double * __restrict __attribute__((aligned(64))) f,
				const double * __restrict __attribute__((aligned(64))) g,
			        double * __restrict __attribute__((aligned(64))) h,
	                        const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,f,g,h,data_len) private(idx)					\
        aligned(a:64,b,c,d,e,f,g,h) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               h[idx] = skx_fp_vector_fraction(     a[idx],
	                                           b[idx],
					           c[idx],
					           d[idx],
						   e[idx],
						   f[idx],
						   g[idx]);
						 
						 
					        
	     
	}
}


/*
     Microcode Sequencer (SMT enabled).
*/
void
skx_smt_MS_samples(const double * __restrict  __attribute__((aligned(64))) a,
		   const double * __restrict  __attribute__((aligned(64))) b,
		   const double * __restrict  __attribute__((aligned(64))) c,
		   const double * __restrict  __attribute__((aligned(64))) d,
		   const double * __restrict  __attribute__((aligned(64))) e,
		   double * __restrict  __attribute__((aligned(64))) f,
	           const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,f,data_len) private(idx)					\
        aligned(a:64,b,c,d,e,f) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               f[idx] = skx_smt_MS(                a[idx],
	                                           b[idx],
					           c[idx],
					           d[idx],
						   e[idx]);
						 
						 
					        
	     
	}
}


/*
     Microcode Sequencer (SMT disabled).
*/
void
skx_smt_MS_samples(const double * __restrict  __attribute__((aligned(64))) a,
		   const double * __restrict  __attribute__((aligned(64))) b,
		   const double * __restrict  __attribute__((aligned(64))) c,
		   const double * __restrict  __attribute__((aligned(64))) d,
		   double * __restrict  __attribute__((aligned(64))) e,
	           const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,d,e,data_len) private(idx)					\
        aligned(a:64,b,c,d,e) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               e[idx] = skx_smt_MS(                a[idx],
	                                           b[idx],
					           c[idx],
					           d[idx]);
					
						 
						 
					        
	     
	}
}


/*
     LLC Local code/data reads hitting in S state in snoop filter per instruction.
*/
void
skx_L3_code_data_read_S_hit_samples(  const double * __restrict __attribute__((aligned(64))) a,
		                      const double * __restrict __attribute__((aligned(64))) b,
		                      double * __restrict __attribute__((aligned(64))) c,
	                              const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(idx)					\
        aligned(a:64,b,c) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               c[idx] = skx_L3_code_data_read_S_hit( a[idx],
	                                             b[idx]);
	}				          
					        
}


/*
     LLC Local code/data reads hitting in E state in snoop filter per instruction.
*/
void
skx_L3_code_data_read_E_hit_samples(  const double * __restrict __attribute__((aligned(64))) a,
		                      const double * __restrict __attribute__((aligned(64))) b,
		                      double * __restrict __attribute__((aligned(64))) c,
	                              const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(idx)					\
        aligned(a:64,b,c) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               c[idx] = skx_L3_code_data_read_E_hit( a[idx],
	                                             b[idx]);
	}				          
					        
}


/*
     LLC Local code/data reads hitting in M/E/F state in snoop filter per instruction.
*/
void
skx_L3_code_data_read_MEF_hit_samples(const double * __restrict __attribute__((aligned(64))) a,
		                      const double * __restrict __attribute__((aligned(64))) b,
		                      double * __restrict __attribute__((aligned(64))) c,
	                              const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(idx)					\
        aligned(a:64,b,c) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               c[idx] = skx_L3_code_data_read_MEF_hit( a[idx],
	                                             b[idx]);
	}				          
					        
}


/*
     INT_MISC.RECOVERY_CYCLES ratio
*/
void
skx_int_misc_recovery_cycles_samples( const double * __restrict  __attribute__((aligned(64))) a,
		                      const double * __restrict  __attribute__((aligned(64))) b,
		                      double * __restrict  __attribute__((aligned(64))) c,
	                              const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(idx)					\
        aligned(a:64,b,c) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               c[idx] = skx_int_misc_recovery_cycles_ratio( a[idx],
	                                             b[idx]);
	}
}


/*
    INT_MISC.CLEAR_RESTEER_CYCLES ratio
*/
void
skx_int_misc_clear_resteer_cycles_samples(const double * __restrict __attribute__((aligned(64))) a,
		                          const double * __restrict __attribute__((aligned(64))) b,
		                          double * __restrict __attribute__((aligned(64))) c,
	                                  const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(idx)					\
        aligned(a:64,b,c) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               c[idx] = skx_int_misc_clear_resteers_cycles_ratio( a[idx],
	                                                          b[idx]);
	}
}


/*
   RS_EVENTS.EMPTY_CYCLES ratio
*/
void
skx_rs_events_empty_cycles_samples(const double * __restrict __attribute__((aligned(64))) a,
		                   const double * __restrict __attribute__((aligned(64))) b,
		                   double * __restrict __attribute__((aligned(64))) c,
	                           const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(idx)					\
        aligned(a:64,b,c) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               c[idx] = skx_rs_events_empty_cycles_ratio( a[idx],
	                                                  b[idx]);
	}
}


/*
    BR_INST_RETIRED.ALL_BRANCHES ratio
*/
void
skx_br_inst_retired_all_branches_samples(const double * __restrict __attribute__((aligned(64))) a,
		                         const double * __restrict __attribute__((aligned(64))) b,
		                         double * __restrict __attribute__((aligned(64))) c,
	                                 const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(idx)					\
        aligned(a:64,b,c) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               c[idx] = skx_br_inst_retired_all_branches_ratio( a[idx],
	                                                        b[idx]);
	}
}


/*
     BR_INST_RETIRED.CONDITIONAL
*/
void
skx_br_inst_retired_cond_samples( const double * __restrict __attribute__((aligned(64))) a,
		                  const double * __restrict __attribute__((aligned(64))) b,
		                  double * __restrict __attribute__((aligned(64))) c,
	                          const int32_t data_len) {

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
#pragma omp parallel for simd schedule(static,8) default(none) \
        shared(a,b,c,data_len) private(idx)					\
        aligned(a:64,b,c) linear(idx:1) unroll partial(2) if(data_len>=100000)
	  for(int32_t idx = 0; idx != data_len; ++idx) {
               c[idx] = skx_br_inst_retired_conditional_ratio( a[idx],
	                                                        b[idx]);
	}
}


					        
              
	
















	

          










