

#ifndef __GMS_SKX_CLK_HW_EVENTS_METRICS_HPP__
#define __GMS_SKX_CLK_HW_EVENTS_METRICS_HPP__


#include <omp.h>

/*
   CPU operating frequency (in GHz)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_cpu_operating_freq(const double CPU_CLK_UNHALTED_THREAD,
                              const double CPU_CLK_UNHALTED_REF_TSC,
                              const double TSC_FREQUENCY) {

  return ((CPU_CLK_UNHALTED_THREAD/
                           CPU_CLK_UNHALTED_REF_TSC*TSC_FREQUENCY)/
                                                      1000000000.0);
}

/*
   CPU utilization (percentage) all cores.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_cpu_utilization(const double CPU_CLK_UNHALTED_REF_TSC,
                           const double TIME_STAMP_CYCLES) {

        return (100.0*CPU_CLK_UNHALTED_REF_TSC/
                                      TIME_STAMP_CYCLES);
}


/*
    CPU utilization (percentage) in kernel mode (all cores).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_cpu_utilization_kernel(const double CPU_CLK_UNHALTED_REF_TSC_SUP,
                                  const double TIME_STAMP_CYCLES) {

        return (100.0*CPU_CLK_UNHALTED_REF_TSC_SUP/
                                        TIME_STAMP_CYCLES);
}


/*
    Cycles Per Instruction (CPI).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_cycles_per_instr(const double CPU_CLK_UNHALTED_THREAD,
                            const double INST_RETIRED_ANY) {

        return (CPU_CLK_UNHALTED_THREAD/
                                    INST_RETIRED_ANY);
}


/*
     Cycles Per Instruction (CPI) kernel mode.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_cycles_per_instr_kernel(const double CPU_CLK_UNHALTED_THREAD_SUP,
                                   const double INST_RETIRED_ANY_SUP) {

        return (CPU_CLK_UNHALTED_THREAD_SUP/
                                      INST_RETIRED_ANY_SUP);
}


/*
    EMON event multiplexing reliability (>95% --- means satisfying ratio).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_emon_mux_reliability(const double CPU_CLK_UNHALTED_THREAD_P,
                                const double CPU_CLK_UNHALTED_THREAD) {

        return (CPU_CLK_UNHALTED_THREAD_P-
                            CPU_CLK_UNHALTED_THREAD < 0.0 ? 
                                               CPU_CLK_UNHALTED_THREAD_P/
                                                             CPU_CLK_UNHALTED_THREAD : 
                                                                   CPU_CLK_UNHALTED_THREAD/
                                                                                CPU_CLK_UNHALTED_THREAD_P);
}


/*
    Branch mispredict ratio.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_branch_mispredict_ratio(const double BR_MISP_RETIRED_ALL_BRANCHES,
                                   const double BR_INST_RETIRED_ALL_BRANCHES) {

        return (BR_MISP_RETIRED_ALL_BRANCHES/
                                  BR_INST_RETIRED_ALL_BRANCHES);
}


/*
    Loads per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_loads_per_instr(const double MEM_INST_RETIRED_ALL_LOADS,
                           const double INST_RETIRED_ANY) {

        return (MEM_INST_RETIRED_ALL_LOADS/
                                  INST_RETIRED_ANY);
}


/*
    Stores per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_stores_per_instr(const double MEM_INST_RETIRED_ALL_STORES,
                            const double INST_RETIRED_ANY) {

        return (MEM_INST_RETIRED_ALL_STORES/
                                   INST_RETIRED_ANY);
}


/*
    Memory operations per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_mem_ops_per_instr(const double MEM_INST_RETIRED_ALL_LOADS,
                             const double MEM_INST_RETIRED_ALL_STORES,
                             const double INST_RETIRED_ANY) {

        const double mem_ops = MEM_INST_RETIRED_ALL_LOADS+
                                 MEM_INST_RETIRED_ALL_STORES;
        return (mem_ops/INST_RETIRED_ANY);
}


/*
    Locks retired per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_locks_per_instr(const double MEM_INST_RETIRED_LOCK_LOADS,
                           const double INST_RETIRED_ANY) {

        return (MEM_INST_RETIRED_LOCK_LOADS/
                                       INST_RETIRED_ANY);
}


/*
    Uncacheable reads per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_uncacheable_reads_per_instr(const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40e33,
                                       const double INST_RETIRED_ANY) {
        
        return (UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40e33/
                                                     INST_RETIRED_ANY);
}


/*
    Streaming-stores (full line) per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_streaming_stores_fl_per_instr(const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x41833,
                                         const double INST_RETIRED_ANY) {

        return (UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x41833/
                                                      INST_RETIRED_ANY);
}


/*
     Streaming-stores (partial line) per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_streaming_stores_pl_per_instr(const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x41a33,
                                         const double INST_RETIRED_ANY) {

        return (UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x41a33/
                                                        INST_RETIRED_ANY);
}


/*
     L1D$ Misses per instruction including:
     -- data
     -- reads for ownership (rfo) with prefetches
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L1D_misses_per_instr(const double L1D_REPLACEMENT,
                                const double INST_RETIRED_ANY) {

        return (L1D_REPLACEMENT/INST_RETIRED_ANY);
}


/*
    L1D$ Demand data read hit per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L1D_hits_per_instr(const double MEM_LOAD_RETIRED_L1_HIT,
                              const double INST_RETIRED_ANY) {

        return (MEM_LOAD_RETIRED_L1_HIT/INST_RETIRED_ANY);
}


/*
     L1$ code read misses (including prefetches) per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L1I_read_misses_per_instr(const double L2_RQSTS_ALL_CODE_RD,
                                     const double INST_RETIRED_ANY) {

        return (L2_RQSTS_ALL_CODE_RD/INST_RETIRED_ANY);
}


/*
     L2D Demand data read hits per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L2_data_read_hits_per_instr(const double MEM_LOAD_RETIRED_L2_HIT,
                                       const double INST_RETIRED_ANY) {

        return (MEM_LOAD_RETIRED_L2_HIT/INST_RETIRED_ANY);
}


/*
     L2 Misses per instruction including:
     -- code
     -- data
     -- read for ownership and prefetches
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L2_all_misses_per_instr(const double L2_LINES_IN_ALL,
                                   const double INST_RETIRED_ANY) {

        return (L2_LINES_IN_ALL/INST_RETIRED_ANY);
}


/*
    L2 demand data read misses per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L2_demand_data_read_mpi(const double MEM_LOAD_RETIRED_L2_MISS,
                                   const double INST_RETIRED_ANY) {

        return (MEM_LOAD_RETIRED_L2_MISS/INST_RETIRED_ANY);
}


/*
    L2 demand data code misses per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L2_demand_code_mpi(const double L2_RQSTS_CODE_RD_MISS,
                              const double INST_RETIRED_ANY) {

        return (L2_RQSTS_CODE_RD_MISS/INST_RETIRED_ANY);
}


/*
    L2 Any local request that HITM in a sibling core (per instruction).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L2_request_hitm_sibling_core(const double OFFCORE_RESPONSE_request_ALL_READS_response_L3_HIT_HITM_OTHER_CORE,
                                        const double INST_RETIRED_ANY) {

        return (OFFCORE_RESPONSE_request_ALL_READS_response_L3_HIT_HITM_OTHER_CORE/
                                                                    INST_RETIRED_ANY);
}


/*
    L2 (percentage) of all lines evicted that are unused prefetches.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L2_lines_evict_unused_prefetch(const double L2_LINES_OUT_USELESS_HWPF,
                                          const double LINES_OUT_NON_SILENT,
                                          const double LINES_OUT_SILENT,
                                          const double HW_Thread_Count) {

        const double term1 = LINES_OUT_SILENT/HW_Thread_Count;
        return (100.0*L2_LINES_OUT_USELESS_HWPF/
                                          (LINES_OUT_NON_SILENT+term1));
}


/*
    L2 percentage of L2 evictions that are allocated into L3$.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double  skx_L2_evict_L3_allocated(const double L2_LINES_OUT_NON_SILENT,
                                  const double IDI_MISC_WB_DOWNGRADE) {

        const double term1 = L2_LINES_OUT_NON_SILENT-
                                                 IDI_MISC_WB_DOWNGRADE;
        return (100.0*term1/L2_LINES_OUT_NON_SILENT);
}


/*
    L2 percentage of L2 evictions that are NOT allocated into L3$.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L2_evict_L3_not_allocated(const double L2_LINES_OUT_NON_SILENT,
                                     const double IDI_MISC_WB_DOWNGRADE) {

        return (100.0*IDI_MISC_WB_DOWNGRADE/
                                        L2_LINES_OUT_NON_SILENT);
}


/*
    LLC code references per instruction (L3 prefetch excluded).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double  skx_L3_code_references_per_instr(const double UNC_CHA_TOR_INSERTS_IA_filter1_0x40433,
                                         const double INST_RETIRED_ANY) {

        return (UNC_CHA_TOR_INSERTS_IA_filter1_0x40433/
                                              INST_RETIRED_ANY);
}


/*
    LLC data read references per instruction (L3 prefetch excluded).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_data_references_per_instr(const double UNC_CHA_TOR_INSERTS_IA_filter1_0x40433,
                                        const double INST_RETIRED_ANY) {

        return (UNC_CHA_TOR_INSERTS_IA_filter1_0x40433/
                                                  INST_RETIRED_ANY);
}


/*
    LLC RFO references per instrctions (L3 prefetch excluded).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_rfo_references_per_instr(const double UNC_CHA_TOR_INSERTS_IA_filter1_0x40033,
                                      const double INST_RETIRED_ANY) {

        return (UNC_CHA_TOR_INSERTS_IA_filter1_0x40033/
                                                     INST_RETIRED_ANY);
}


/*
    LLC Misses Per Instruction (includes code and data and rfo with prefetches).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_all_mpi(const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x12D40433,
                   const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x12CC0233,
                   const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x12C40033,
                   const double INST_RETIRED_ANY) {

        const double sum = UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x12D40433+
                           UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x12CC0233+
                           UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x12C40033;
        return (sum/INST_RETIRED_ANY);
}


/*
    LLC data read Misses Per Instruction (demand and prefetch).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_data_read_mpi(const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x12D40433,
                            const double INST_RETIRED_ANY) {

        return (UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x12D40433/
                                                      INST_RETIRED_ANY);
}


/*
    LLC RFO read Misses Per Instruction (demand and prefetch).  
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_rfo_read_mpi(const uint64_t UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x12C40033,
                           const uint64_t INST_RETIRED_ANY) {

        return ((double)(UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x12C40033/
                                                        INST_RETIRED_ANY));
}


/*
     LLC code read Misses Per Instruction (demand and prefetch).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_code_read_mpi(const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x12CC0233,
                            const double INST_RETIRED_ANY) {

        return (UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x12CC0233/
                                                           INST_RETIRED_ANY);
}


/*
      LLC code read MPI (demand and prefetch).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_total_hitm_per_instr(const double OFFCORE_RESPONSE_request_ALL_READS_response_L3_MISS_REMOTE_HITM,
                                const double INST_RETIRED_ANY) {

       return (OFFCORE_RESPONSE_request_ALL_READS_response_L3_MISS_REMOTE_HITM/
                                                                   INST_RETIRED_ANY);
}


/*
        LLC total HIT clean line forwards (per instr) (excludes LLC prefetches)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_total_hitm_clean_lines_per_instr(const double OFFCORE_RESPONSE_request_ALL_READS_response_L3_MISS_REMOTE_HIT_FORWARD,
                                            const double INST_RETIRED_ANY) {

       return (OFFCORE_RESPONSE_request_ALL_READS_response_L3_MISS_REMOTE_HIT_FORWARD/
                                                                          INST_RETIRED_ANY);
}


/*
        Average LLC data read (demand+prefetch) miss latency (in ns). 
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_avg_data_read_ns(const double UNC_CHA_TOR_OCCUPANCY_IA_MISS_filter1_0x40433,
                               const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40433,
			       const double UNC_CHA_CLOCKTICKS,
			       const double TIME_INTERVAL_SEC) {

        const double term1 = UNC_CHA_TOR_OCCUPANCY_IA_MISS_filter1_0x40433/
	                     UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40433;
	return (1000000000.0*TIME_INTERVAL_SEC*term1/UNC_CHA_CLOCKTICKS);
}


/*
      Average LLC data read (demand and prefetch) miss latency (in UNCORE clk).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_avg_data_read_unc_clk(const double UNC_CHA_TOR_OCCUPANCY_IA_MISS_filter1_0x40433,
                                    const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40433) {

       return (UNC_CHA_TOR_OCCUPANCY_IA_MISS_filter1_0x40433/
               UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40433);
}


/*
    Average LLC data read (demand+prefetch) miss latency for LOCAL requests (in ns)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_avg_data_read_loc_req_ns(const double UNC_CHA_TOR_OCCUPANCY_IA_MISS_filter1_0x40432,
                                       const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40432,
				       const double UNC_CHA_CLOCKTICKS,
			               const double TIME_INTERVAL_SEC) {
    
        const double term1 = UNC_CHA_TOR_OCCUPANCY_IA_MISS_filter1_0x40432/
	                     UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40432;
	return (1000000000.0*TIME_INTERVAL_SEC*term1/UNC_CHA_CLOCKTICKS);
}


/*
     Average LLC data read (demand+prefetch) miss latency for REMOTE requests (in ns).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_avg_data_read_rem_ns(const double UNC_CHA_TOR_OCCUPANCY_IA_MISS_filter1_0x40431,
                                       const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40431,
				       const double UNC_CHA_CLOCKTICKS,
			               const double TIME_INTERVAL_SEC) {

        const double term1 = UNC_CHA_TOR_OCCUPANCY_IA_MISS_filter1_0x40431/
	                     UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40431;
	return (1000000000.0*TIME_INTERVAL_SEC*term1/UNC_CHA_CLOCKTICKS);
}


/*
      Average LLC data read (demand+prefetch) miss latency  for LOCAL requests (in UNCORE clk)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_avg_data_read_loc_unc_clk(const double UNC_CHA_TOR_OCCUPANCY_IA_MISS_filter1_0x40432,
                                            const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40432) {

       return ( UNC_CHA_TOR_OCCUPANCY_IA_MISS_filter1_0x40432/
                UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40432);
}


/*
      Average LLC data read (demand+prefetch) miss latency  for REMOTE requests (in UNCORE clk)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_avg_data_read_rem_unc_clk(const double UNC_CHA_TOR_OCCUPANCY_IA_MISS_filter1_0x40432,
                                        const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40432) {

       return (UNC_CHA_TOR_OCCUPANCY_IA_MISS_filter1_0x40432/
               UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40432);
}


/*
   ITLB MPI. 
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_ITLB_mpi(const double ITLB_MISSES_WALK_COMPLETED,
                    const double INST_RETIRED_ANY) {

       return (ITLB_MISSES_WALK_COMPLETED/INST_RETIRED_ANY);
}


/*
    ITLB large page MPI
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_ITLB_2M_4M_mpi(const double ITLB_MISSES_WALK_COMPLETED_2M_4M,
                          const double INST_RETIRED_ANY) {

       return (ITLB_MISSES_WALK_COMPLETED_2M_4M/INST_RETIRED_ANY);
}


/*
    DTLB load MPI.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_DTLB_load_mpi(const double DTLB_LOAD_MISSES_WALK_COMPLETED,
                         const double INST_RETIRED_ANY) {

       return (DTLB_LOAD_MISSES_WALK_COMPLETED/INST_RETIRED_ANY);
}


/*
    DTLB large page load MPI.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_DTLB_2M_4M_mpi(const double DTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M,
                          const double INST_RETIRED_ANY) {

       return (DTLB_LOAD_MISSES_WALK_COMPLETED_2M_4M/
                                        INST_RETIRED_ANY);
}


/*
    DTLB store MPI.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_DTLB_store_mpi(const double DTLB_STORE_MISSES_WALK_COMPLETED,
                          const double INST_RETIRED_ANY) {

       return (DTLB_STORE_MISSES_WALK_COMPLETED/
                                INST_RETIRED_ANY);
}


/*
    DTLB load miss latency (in core clks)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_DTLB_load_miss_clks(const double DTLB_LOAD_MISSES_WALK_ACTIVE,
                               const double DTLB_LOAD_MISSES_WALK_COMPLETED) {

       return (DTLB_LOAD_MISSES_WALK_ACTIVE/
                DTLB_LOAD_MISSES_WALK_COMPLETED);
}


/*
    DTLB store miss latency (in core clks).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_DTLB_store_miss_clks(const double DTLB_STORE_MISSES_WALK_ACTIVE,
                                const double DTLB_STORE_MISSES_WALK_COMPLETED) {

       return ( DTLB_STORE_MISSES_WALK_ACTIVE/
                DTLB_STORE_MISSES_WALK_COMPLETED);
}          


/*
    ITLB miss latency (in core clks).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_ITLB_miss_latency_clks(const double ITLB_MISSES_WALK_ACTIVE,
                                  const double ITLB_MISSES_WALK_COMPLETED) {

       return (ITLB_MISSES_WALK_ACTIVE/
                 ITLB_MISSES_WALK_COMPLETED);
}


/*
     NUMA percentage of Reads addressed to local DRAM.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_numa_reads_local_dram(const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40432,
                                 const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40431) {

       const double sum = UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40432+
                          UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40431;
       return (100.0*UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40432/sum);
}


/*
     NUMA percentage of Reads addressed to remote  DRAM.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_numa_reads_remote_dram(const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40432,
                                  const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40431){

       const double sum = UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40432+
                          UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40431;
       return (100.0*UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40431/sum);
}


/*
      Uncore Frequency Ghz.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_uncore_frequency_ghz(const double UNC_CHA_CLOCKTICKS,
                                const double HW_CORE_COUNT,
				const double TIME_INTERVAL_SECONDS) {
                              
       const double term1 = UNC_CHA_CLOCKTICKS/HW_CORE_COUNT;
       const double term2 = 1000000000.0/TIME_INTERVAL_SECONDS;
       return (term1/term2);
}


/*
    UPI speed - GT/s.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_UPI_speed(const double UNC_UPI_CLOCKTICKS,
                     const double TIME_INTERVAL_SECONDS) {

       const double term1 = 8.0*UNC_UPI_CLOCKTICKS;
       return (term1/1000000000.0/TIME_INTERVAL_SECONDS);
}


/*
    UPI Data transmit BW (MB/sec) (only data)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_UPI_data_bw(const double UNC_UPI_TxL_FLITS_ALL_DATA,
                       const double TIME_INTERVAL_SECONDS) {

       const double term1 = UNC_UPI_TxL_FLITS_ALL_DATA*
                            7.1111111111111111111111;
       const double term2 = 1000000/TIME_INTERVAL_SECONDS;
       return (term1/term2);
}


/*
    UPI Total transmit BW (MB/sec) (includes control)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_UPI_transmit_total_bw(const double UNC_UPI_TxL_FLITS_ALL_DATA,
                                 const double UNC_UPI_TxL_FLITS_NON_DATA,
				 const double TIME_INTERVAL_SECONDS) {
       
       const double term1 =  UNC_UPI_TxL_FLITS_ALL_DATA+
                             UNC_UPI_TxL_FLITS_NON_DATA*
			     7.1111111111111111111111;
       const double term2 = 1000000/TIME_INTERVAL_SECONDS;
       return (term1/term2);
}


/*
   UPI Transmit utilization percentage (includes control).
   Percentage of time the processor is communicating with sibling processors.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_UPI_transmit_utilization(const double UNC_UPI_CLOCKTICKS,
                                    const double UNC_UPI_L1_POWER_CYCLES,
				    const double FREERUN_PKG_C6_RESIDENCY,
				    const double UNC_UPI_TxL_FLITS_ALL_DATA,
				    const double UNC_UPI_TxL_FLITS_NON_DATA,
				    const double TIME_INTERVAL_SECONDS) {

        const double term1 = UNC_UPI_TxL_FLITS_ALL_DATA+
                             UNC_UPI_TxL_FLITS_NON_DATA*
			     0.333333333333333333333333;
	const double term2 = TIME_INTERVAL_SECONDS/(TIME_INTERVAL_SECONDS-
	                                            FREERUN_PKG_C6_RESIDENCY);
	const double term3 = UNC_UPI_CLOCKTICKS-UNC_UPI_L1_POWER_CYCLES;
	return (100.0*term1/term2*term3*0.8333333333333333333333);
}


/*
   UPI percentage of  cycles transmit link is half-width (L0p) 
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_half_width_link_tx_cycles(const double UNC_UPI_TxL0P_POWER_CYCLES,
                                     const double UNC_UPI_CLOCKTICKS,
				     const double UNC_UPI_L1_POWER_CYCLES) {

        const double term1 = UNC_UPI_TxL0P_POWER_CYCLES/
	                     (UNC_UPI_CLOCKTICKS-UNC_UPI_L1_POWER_CYCLES);
	return (100.0*term1);
}


/*
    UPI percentage of  cycles receive link is half-width (L0p)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_half_width_link_rx_cycles(const double UNC_UPI_RxL0P_POWER_CYCLES,
                                     const double UNC_UPI_CLOCKTICKS,
				     const double UNC_UPI_L1_POWER_CYCLES) {

        const double term1 = UNC_UPI_RxL0P_POWER_CYCLES/
	                     (UNC_UPI_CLOCKTICKS-UNC_UPI_L1_POWER_CYCLES);
	return (100.0*term1);
}


/*
     HA - Reads vs. all requests
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_HA_reads_vs_all_req(const double UNC_CHA_REQUESTS_READS,
                               const double UNC_CHA_REQUESTS_WRITES) {

       return (UNC_CHA_REQUESTS_READS/
              (UNC_CHA_REQUESTS_READS+UNC_CHA_REQUESTS_WRITES));
}


/*
     HA - Writes vs. all requests
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_HA_writes_vs_all_req(const double UNC_CHA_REQUESTS_READS,
                                const double UNC_CHA_REQUESTS_WRITES) {

       return (UNC_CHA_REQUESTS_WRITES/
              (UNC_CHA_REQUESTS_READS+UNC_CHA_REQUESTS_WRITES));
}


/*
     HA percentage of all reads that are local.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_HA_all_reads_local(const double UNC_CHA_REQUESTS_READS_LOCAL,
			      const double UNC_CHA_REQUESTS_READS) {

       return (100.0*UNC_CHA_REQUESTS_READS_LOCAL/
                       UNC_CHA_REQUESTS_READS);
}


/*
     HA percentage of all writes that are local.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_HA_all_writes_local(const double UNC_CHA_REQUESTS_WRITES_LOCAL,
			      const double UNC_CHA_REQUESTS_WRITES) {

       return (100.0*UNC_CHA_REQUESTS_WRITES_LOCAL/
                       UNC_CHA_REQUESTS_WRITES);
}


/*
    HA conflict responses per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_HA_conflict_resp(const double UNC_CHA_SNOOP_RESP_RSPCNFLCTS,
                            const double INST_RETIRED_ANY) {

       return (UNC_CHA_SNOOP_RESP_RSPCNFLCTS/
                                  INST_RETIRED_ANY);
}


/*
    HA directory lookups that spawned a snoop (per instruction)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_HA_dir_lookup_snoop(const double UNC_CHA_DIR_LOOKUP_SNP,
                               const double INST_RETIRED_ANY) {

       return (UNC_CHA_DIR_LOOKUP_SNP/
                          INST_RETIRED_ANY);
}


/*
    HA directory lookups that did not spawn a snoop (per instruction).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_HA_dir_lookup_no_snoop(const double UNC_CHA_DIR_LOOKUP_NO_SNP,
                               const double INST_RETIRED_ANY) {

       return (UNC_CHA_DIR_LOOKUP_NO_SNP/
                          INST_RETIRED_ANY);
}


/*
   M2M directory updates (per instruction). 
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_M2M_dir_updates(const double UNC_CHA_DIR_UPDATE_HA,
                           const double UNC_CHA_DIR_UPDATE_TOR,
			   const double UNC_M2M_DIRECTORY_UPDATE_ANY,
			   const double INST_RETIRED_ANY) {

        const double term1 =  UNC_CHA_DIR_UPDATE_HA+
	                      UNC_CHA_DIR_UPDATE_TOR+
			      UNC_M2M_DIRECTORY_UPDATE_ANY;
        return (term1/INST_RETIRED_ANY);
}


/*
    M2M extra reads from XPT-UPI prefetches (per instruction).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_M2M_reads_XPT_UPI_prefetch(const double UNC_M2M_PREFCAM_INSERTS,
                                      const double UNC_M2M_PREFCAM_DEMAND_PROMOTIONS,
				      const double INST_RETIRED_ANY) {

       return (UNC_M2M_PREFCAM_INSERTS-
                 UNC_M2M_PREFCAM_DEMAND_PROMOTIONS/
		                       INST_RETIRED_ANY);
}


/*
    DDR data rate (MT/sec).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_DDR_data_rate(const double UNC_M_CLOCKTICKS,
                         const double TIME_INTERVAL_SECONDS) {

       return(2.0*UNC_M_CLOCKTICKS/1000000.0/
                           TIME_INTERVAL_SECONDS);
}


/*
    Memory bandwidth read (MB/sec).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_memory_read_bw(const double UNC_M_CAS_COUNT_RD,
                          const double TIME_INTERVAL_SECONDS) {

       return (64.0*UNC_M_CAS_COUNT_RD/1000000.0/
                            TIME_INTERVAL_SECONDS);
}


/*
   Load instructions per memory bandwidth.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_load_mem_instr_mem_bw(const double UNC_M_CAS_COUNT_RD,
				 const double MEM_INST_RETIRED_ALL_LOADS,
                                 const double TIME_INTERVAL_SECONDS) {

       const double mem_read_bw = skx_memory_read_bw(UNC_M_CAS_COUNT_RD,
                                                     TIME_INTERVAL_SECONDS);
       return (MEM_INST_RETIRED_ALL_LOADS/mem_read_bw);
}


/*
    Memory bandwidth write  (MB/sec).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_memory_write_bw(const double UNC_M_CAS_COUNT_WR,
                          const double TIME_INTERVAL_SECONDS) {

       return (64.0*UNC_M_CAS_COUNT_WR/1000000.0/
                            TIME_INTERVAL_SECONDS);
}


/*
    Store instructions per memory bandwidth.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_store_mem_instr_mem_bw(const double UNC_M_CAS_COUNT_WR,
				 const double MEM_INST_RETIRED_ALL_STORES,
                                 const double TIME_INTERVAL_SECONDS) {

       const double mem_write_bw = skx_memory_write_bw(UNC_M_CAS_COUNT_WR,
                                                       TIME_INTERVAL_SECONDS);
       return (MEM_INST_RETIRED_ALL_STORES/mem_write_bw);
}


/*
    Memory bandwidth total (MB/sec).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_memory_bw_total(const double UNC_M_CAS_COUNT_RD,
                           const double UNC_M_CAS_COUNT_WR,
			   const double TIME_INTERVAL_SECONDS) {

        const double term1 = UNC_M_CAS_COUNT_RD+
	                     UNC_M_CAS_COUNT_WR;
	return (term1*64.0/1000000.0/TIME_INTERVAL_SECONDS);
}


/*
    Load and store instructions per total memory bandwidth.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_total_mem_instr_mem_bw(const double UNC_M_CAS_COUNT_WR,
				  const double MEM_INST_RETIRED_ALL_STORES,
				  const double UNC_M_CAS_COUNT_RD,
				  const double MEM_INST_RETIRED_ALL_LOADS,
                                  const double TIME_INTERVAL_SECONDS) {

        const double mem_write_bw = skx_memory_write_bw(UNC_M_CAS_COUNT_WR,
                                                       TIME_INTERVAL_SECONDS);
	const double mem_read_bw =  skx_memory_read_bw(UNC_M_CAS_COUNT_RD,
                                                     TIME_INTERVAL_SECONDS);
	return ((MEM_INST_RETIRED_ALL_STORES/mem_write_bw)+
	        (MEM_INST_RETIRED_ALL_LOADS/mem_read_bw));
}


/*
   Memory extra read b/w due to XPT prefetches (MB/sec).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_mem_bw_xpt_prefetch(const double UNC_M2M_PREFCAM_INSERTS,
                               const double UNC_M2M_PREFCAM_DEMAND_PROMOTIONS,
			       const double TIME_INTERVAL_SECONDS) {

       const double term1 = UNC_M2M_PREFCAM_INSERTS-
                            UNC_M2M_PREFCAM_DEMAND_PROMOTIONS;
       return (term1*64.0/1000000.0/TIME_INTERVAL_SECONDS);
}


/*
     Memory extra write b/w due to directory updates (MB/sec).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_mem_bw_dir_update(const double UNC_CHA_DIR_UPDATE_HA,
                             const double UNC_CHA_DIR_UPDATE_TOR,
			     const double UNC_M2M_DIRECTORY_UPDATE_ANY,
			     const double TIME_INTERVAL_SECONDS) {

        const double term1 = UNC_CHA_DIR_UPDATE_HA+
	                     UNC_CHA_DIR_UPDATE_TOR+
			     UNC_M2M_DIRECTORY_UPDATE_ANY;
	return (term1*64.0/1000000.0/TIME_INTERVAL_SECONDS);
}


/*
   DRAM RPQ read latency (ns).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_DRAM_rpq_read_latency_ns(const double UNC_M_CLOCKTICKS,
                                    const double UNC_M_RPQ_INSERTS,
				    const double UNC_M_RPQ_OCCUPANCY,
				    const double TIME_INTERVAL_SECONDS) {

       const double term1 = (UNC_M_RPQ_OCCUPANCY/UNC_M_RPQ_INSERTS)/
                                                   UNC_M_CLOCKTICKS;
       return (term1*TIME_INTERVAL_SECONDS*1000000000.0);			 
}


/*
     DRAM WPQ write latency (ns)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_DRAM_rpq_write_latency_ns(const double UNC_M_CLOCKTICKS,
                                    const double UNC_M_WPQ_INSERTS,
				    const double UNC_M_WPQ_OCCUPANCY,
				    const double TIME_INTERVAL_SECONDS) {

       const double term1 = (UNC_M_WPQ_OCCUPANCY/UNC_M_WPQ_INSERTS)/
                                                   UNC_M_CLOCKTICKS;
       return (term1*TIME_INTERVAL_SECONDS*1000000000.0);			 
}


/*
    Memory average number of entries in each read Q (RPQ)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_mem_avg_rpq_read(const double UNC_M_CLOCKTICKS,
                            const double UNC_M_RPQ_OCCUPANCY) {

       return (UNC_M_RPQ_OCCUPANCY/UNC_M_CLOCKTICKS);
}


/*
   memory average number of entries in each write Q (WPQ).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_mem_avg_rpq_write(const double UNC_M_CLOCKTICKS,
                            const double UNC_M_WPQ_OCCUPANCY) {

       return (UNC_M_WPQ_OCCUPANCY/UNC_M_CLOCKTICKS);
}


/*
    I/O bandwidth disk or network writes (MB/sec).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_IO_disk_or_net_bw_writes(const double UNC_IIO_DATA_REQ_OF_CPU_MEM_READ_PART0,
				    const double UNC_IIO_DATA_REQ_OF_CPU_MEM_READ_PART1,
				    const double UNC_IIO_DATA_REQ_OF_CPU_MEM_READ_PART2,
				    const double UNC_IIO_DATA_REQ_OF_CPU_MEM_READ_PART3,
				    const double TIME_INTERVAL_SECONDS) {

       const double term1 = UNC_IIO_DATA_REQ_OF_CPU_MEM_READ_PART0+
                            UNC_IIO_DATA_REQ_OF_CPU_MEM_READ_PART1+
			    UNC_IIO_DATA_REQ_OF_CPU_MEM_READ_PART2+
			    UNC_IIO_DATA_REQ_OF_CPU_MEM_READ_PART3;
       return (term1*4.0/1000000.0/TIME_INTERVAL_SECONDS);
}


/*
     I/O bandwidth disk or network reads (MB/sec).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_IO_disk_or_net_bw_reads(const double UNC_IIO_DATA_REQ_OF_CPU_MEM_WRITE_PART0,
				    const double UNC_IIO_DATA_REQ_OF_CPU_MEM_WRITE_PART1,
				    const double UNC_IIO_DATA_REQ_OF_CPU_MEM_WRITE_PART2,
				    const double UNC_IIO_DATA_REQ_OF_CPU_MEM_WRITE_PART3,
				    const double TIME_INTERVAL_SECONDS) {

       const double term1 = UNC_IIO_DATA_REQ_OF_CPU_MEM_WRITE_PART0+
                            UNC_IIO_DATA_REQ_OF_CPU_MEM_WRITE_PART1+
			    UNC_IIO_DATA_REQ_OF_CPU_MEM_WRITE_PART2+
			    UNC_IIO_DATA_REQ_OF_CPU_MEM_WRITE_PART3;
       return (term1*4.0/1000000.0/TIME_INTERVAL_SECONDS);
}


/*
    I/O bandwidth disk or network (MB/sec)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_IO_bw_disk_net_total(const double UNC_IIO_DATA_REQ_OF_CPU_MEM_READ_PART0,
				const double UNC_IIO_DATA_REQ_OF_CPU_MEM_READ_PART1,
				const double UNC_IIO_DATA_REQ_OF_CPU_MEM_READ_PART2,
				const double UNC_IIO_DATA_REQ_OF_CPU_MEM_READ_PART3,
				const double UNC_IIO_DATA_REQ_OF_CPU_MEM_WRITE_PART0,
				const double UNC_IIO_DATA_REQ_OF_CPU_MEM_WRITE_PART1,
				const double UNC_IIO_DATA_REQ_OF_CPU_MEM_WRITE_PART2,
				const double UNC_IIO_DATA_REQ_OF_CPU_MEM_WRITE_PART3,
				const double TIME_INTERVAL_SECONDS) {

        const double term1 = UNC_IIO_DATA_REQ_OF_CPU_MEM_WRITE_PART0+
                             UNC_IIO_DATA_REQ_OF_CPU_MEM_WRITE_PART1+
			     UNC_IIO_DATA_REQ_OF_CPU_MEM_WRITE_PART2+
			     UNC_IIO_DATA_REQ_OF_CPU_MEM_WRITE_PART3;
	const double term2 = UNC_IIO_DATA_REQ_OF_CPU_MEM_READ_PART0+
                             UNC_IIO_DATA_REQ_OF_CPU_MEM_READ_PART1+
			     UNC_IIO_DATA_REQ_OF_CPU_MEM_READ_PART2+
			     UNC_IIO_DATA_REQ_OF_CPU_MEM_READ_PART3;
	return ((term1+term2)*4.0/1000000.0/TIME_INTERVAL_SECONDS);
}


/*
     I/O number of partial PCI writes per second.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_PCI_part_writes_sec(const double UNC_CHA_TOR_INSERTS_IO_HIT_filter1_0x40033,
                               const double UNC_CHA_TOR_INSERTS_IO_MISS_filter1_0x40033,
			       const double TIME_INTERVAL_SECONDS) {

        return ((UNC_CHA_TOR_INSERTS_IO_HIT_filter1_0x40033+
	        UNC_CHA_TOR_INSERTS_IO_MISS_filter1_0x40033)/
		TIME_INTERVAL_SECONDS);
}


/*
    I/O write cache miss(disk/network reads) bandwidth (MB/sec)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_IO_write_cache_miss_bw(const double UNC_CHA_TOR_INSERTS_IO_MISS_filter1_0x49033,
                                  const double UNC_CHA_TOR_INSERTS_IO_MISS_filter1_0x40033,
				  const double TIME_INTERVAL_SECONDS) {

       const double term1 = UNC_CHA_TOR_INSERTS_IO_MISS_filter1_0x49033+
                            UNC_CHA_TOR_INSERTS_IO_MISS_filter1_0x40033;
       return (term1*64.0/1000000.0/TIME_INTERVAL_SECONDS);
}


/*
   I/O read cache miss(disk/network writes) bandwidth (MB/sec).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_IO_read_cache_miss_bw(const double UNC_CHA_TOR_INSERTS_IO_MISS_filter1_0x43c33,
                                 const double TIME_INTERVAL_SECONDS) {

       return (UNC_CHA_TOR_INSERTS_IO_MISS_filter1_0x43c33*64.0/
               1000000.0/TIME_INTERVAL_SECONDS);
}


/*
    IO cache miss(disk/network) bandwidth (MB/sec)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_IO_cache_miss_total_bw(const double UNC_CHA_TOR_INSERTS_IO_MISS_filter1_0x49033,
                                  const double UNC_CHA_TOR_INSERTS_IO_MISS_filter1_0x40033,
				  const double UNC_CHA_TOR_INSERTS_IO_MISS_filter1_0x43c33,
				  const double TIME_INTERVAL_SECONDS) {

       const double term1 = UNC_CHA_TOR_INSERTS_IO_MISS_filter1_0x49033+
                            UNC_CHA_TOR_INSERTS_IO_MISS_filter1_0x40033+
			    UNC_CHA_TOR_INSERTS_IO_MISS_filter1_0x43c33;
       return (term1*64.0/1000000.0/TIME_INTERVAL_SECONDS);
}


/*
     MMIO reads per second.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_MMIO_reads_sec(const double UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40040e33,
                          const double TIME_INTERVAL_SECONDS) {

       return (UNC_CHA_TOR_INSERTS_IA_MISS_filter1_0x40040e33/
                TIME_INTERVAL_SECONDS);
}


/*
    Memory Page Empty vs. all requests
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_mem_page_empty_all_reqs(const double UNC_M_PRE_COUNT_RD_u0xc,
                                   const double UNC_M_PRE_COUNT_PAGE_MISS,
				   const double UNC_M_CAS_COUNT_RD,
				   const double UNC_M_CAS_COUNT_WR) {

       const double term1 = UNC_M_PRE_COUNT_RD_u0xc-
                             UNC_M_PRE_COUNT_PAGE_MISS;
       const double term2 = UNC_M_CAS_COUNT_RD+
                            UNC_M_CAS_COUNT_WR;
       return (term1/term2);
}


/*
    Memory Page Misses vs. all requests
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_mem_page_misses_all_req( const double UNC_M_PRE_COUNT_PAGE_MISS,
				    const double UNC_M_CAS_COUNT_RD,
				    const double UNC_M_CAS_COUNT_WR) {

       return ( UNC_M_PRE_COUNT_PAGE_MISS/(UNC_M_CAS_COUNT_RD+
                                           UNC_M_CAS_COUNT_WR));
}


/*
   Memory Page Hits vs. all requests
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_mem_page_hits_all_req(const double UNC_M_PRE_COUNT_RD_u0xc,
                                 const double UNC_M_CAS_COUNT_RD,
				 const double UNC_M_CAS_COUNT_WR) {

       const double term1 = UNC_M_PRE_COUNT_RD_u0xc/
                            ( UNC_M_CAS_COUNT_RD+UNC_M_CAS_COUNT_WR);
       return (1.0-term1);
}


/*
   Memory percentage  Cycles where all DRAM ranks are in PPD mode.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_DRAM_PPD_mode_cycles(const double UNC_M_POWER_CHANNEL_PPD,
                                const double UNC_M_CLOCKTICKS) {

       return (100.0*UNC_M_POWER_CHANNEL_PPD/
                                UNC_M_CLOCKTICKS);
}


/*
   Memory percentage Cycles all ranks in critical thermal throttle.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_mem_cycles_thermal_throttled(const double UNC_M_POWER_CRITICAL_THROTTLE_CYCLES,
                                        const double UNC_M_CLOCKTICKS) {

       return (100.0*UNC_M_POWER_CRITICAL_THROTTLE_CYCLES/
                                           UNC_M_CLOCKTICKS);
}


/*
   Memory  Cycles Memory is in self refresh power mode
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_mem_cycles_self_refresh(const double UNC_M_POWER_SELF_REFRESH,
                                   const double UNC_M_CLOCKTICKS) {

       return (100.0*UNC_M_POWER_SELF_REFRESH/
                                  UNC_M_CLOCKTICKS);
}


/*
   Uops delivered from decoded Icache (DSB).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_DSB_uops_delivered(const double IDQ_DSB_UOPS,
                              const double IDQ_MITE_UOPS,
			      const double IDQ_MS_UOPS,
			      const double LSD_UOPS) {

        const double term1 = IDQ_DSB_UOPS+
	                     IDQ_MITE_UOPS+
			     IDQ_MS_UOPS+
			     LSD_UOPS;
	return (IDQ_DSB_UOPS/term1);
}


/*
   Uops delivered from legacy decode pipeline (MITE).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_MITE_uops_delivered(const double IDQ_DSB_UOPS,
                              const double IDQ_MITE_UOPS,
			      const double IDQ_MS_UOPS,
			      const double LSD_UOPS) {

        const double term1 = IDQ_DSB_UOPS+
	                     IDQ_MITE_UOPS+
			     IDQ_MS_UOPS+
			     LSD_UOPS;
        return (IDQ_MITE_UOPS/term1);
}


/*
    Uops delivered from microcode sequencer (MS).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_MS_uops_delivered(const double IDQ_DSB_UOPS,
                              const double IDQ_MITE_UOPS,
			      const double IDQ_MS_UOPS,
			      const double LSD_UOPS) {

        const double term1 = IDQ_DSB_UOPS+
	                     IDQ_MITE_UOPS+
			     IDQ_MS_UOPS+
			     LSD_UOPS;
        return (IDQ_MS_UOPS/term1);
}


/*
   Uops delivered from loop stream detector (LSD). 
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_LSD_uops_delivered(const double IDQ_DSB_UOPS,
                              const double IDQ_MITE_UOPS,
			      const double IDQ_MS_UOPS,
			      const double LSD_UOPS) {

        const double term1 = IDQ_DSB_UOPS+
	                     IDQ_MITE_UOPS+
			     IDQ_MS_UOPS+
			     LSD_UOPS;
        return (LSD_UOPS/term1);
}


/*
    FP scalar single-precision FP instructions retired per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_fp32_scalar_retired(const double FP_ARITH_INST_RETIRED_SCALAR_SINGLE,
                                const double INST_RETIRED_ANY) {

       return (FP_ARITH_INST_RETIRED_SCALAR_SINGLE/
                                        INST_RETIRED_ANY);
}


/*
   FP scalar double-precision FP instructions retired per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_fp64_scalar_retired(const double FP_ARITH_INST_RETIRED_SCALAR_DOUBLE,
                                const double INST_RETIRED_ANY) {

       return (FP_ARITH_INST_RETIRED_SCALAR_DOUBLE/
                                        INST_RETIRED_ANY);
}


/*
   FP 128-bit packed single-precision FP instructions retired per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_fp32_vec128b_retired(const double FP_ARITH_INST_RETIRED_128B_PACKED_SINGLE,
                                const double INST_RETIRED_ANY) {

       return (FP_ARITH_INST_RETIRED_128B_PACKED_SINGLE/
                                    INST_RETIRED_ANY);
}


/*
   FP 128-bit packed double-precision FP instructions retired per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_fp64_vec128b_retired(const double FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE,
                                const double INST_RETIRED_ANY) {

       return (FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE/
                                    INST_RETIRED_ANY);
}


/*
   FP 256-bit packed single-precision FP instructions retired per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_fp32_vec256b_retired(const double FP_ARITH_INST_RETIRED_256B_PACKED_SINGLE,
                                const double INST_RETIRED_ANY) {

       return (FP_ARITH_INST_RETIRED_256B_PACKED_SINGLE/
                                    INST_RETIRED_ANY);
}


/*
    FP 256-bit packed double-precision FP instructions retired per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_fp64_vec256b_retired(const double FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE,
                                const double INST_RETIRED_ANY) {

       return (FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE/
                                    INST_RETIRED_ANY);
}


/*
    FP 512-bit packed single-precision FP instructions retired per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_fp32_vec512b_retired(const double FP_ARITH_INST_RETIRED_512B_PACKED_SINGLE,
                                const double INST_RETIRED_ANY) {

       return (FP_ARITH_INST_RETIRED_512B_PACKED_SINGLE/
                                    INST_RETIRED_ANY);
}


/*
    FP 512-bit packed double-precision FP instructions retired per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_fp64_vec512b_retired(const double FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE,
                                const double INST_RETIRED_ANY) {

       return (FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE/
                                    INST_RETIRED_ANY);
}


/*
    FP instruction density (percentage).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_fp_instructions_density(const double FP_ARITH_INST_RETIRED_SCALAR_SINGLE,
                                   const double FP_ARITH_INST_RETIRED_SCALAR_DOUBLE,
				   const double FP_ARITH_INST_RETIRED_128B_PACKED_SINGLE,
				   const double FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE,
				   const double FP_ARITH_INST_RETIRED_256B_PACKED_SINGLE,
				   const double FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE,
				   const double FP_ARITH_INST_RETIRED_512B_PACKED_SINGLE,
				   const double FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE,
				   const double INST_RETIRED_ANY) {

       const double term1 = FP_ARITH_INST_RETIRED_SCALAR_SINGLE+
                            FP_ARITH_INST_RETIRED_SCALAR_DOUBLE+
			    FP_ARITH_INST_RETIRED_128B_PACKED_SINGLE+
			    FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE;
       const double term2 = FP_ARITH_INST_RETIRED_256B_PACKED_SINGLE+
                            FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE+
			    FP_ARITH_INST_RETIRED_512B_PACKED_SINGLE+
			    FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE;
       return ((term1+term2)/INST_RETIRED_ANY);
}


/*
    Branch instruction density
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_branch_instr_density(const double BR_INST_RETIRED_ALL_BRANCHES,
                                const double INST_RETIRED_ANY) {

       return (100.0*BR_INST_RETIRED_ALL_BRANCHES/
                                      INST_RETIRED_ANY);
}


/*
    DRAM power (watts).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_DRAM_power(const double FREERUN_DRAM_ENERGY_STATUS,
                      const double TIME_INTERVAL_SECONDS) {

       const double term1 =  FREERUN_DRAM_ENERGY_STATUS*15.3;
       return (term1/1000000.0/TIME_INTERVAL_SECONDS);
}


/*
    Package power (watts)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_Package_power(const double FREERUN_PKG_ENERGY_STATUS,
                         const double TIME_INTERVAL_SECONDS) {

       const double term1 = FREERUN_PKG_ENERGY_STATUS*61.0;
       return (term1/1000000.0/TIME_INTERVAL_SECONDS);
}


/*
    Core c3 residency
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_Core_C3_residency(const double FREERUN_CORE_C3_RESIDENCY,
                        const double TSC) {

       return (100.0*FREERUN_CORE_C3_RESIDENCY/TSC);
}


/*
   Core c6 residency
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_Core_C6_residency(const double FREERUN_CORE_C6_RESIDENCY,
                             const double TSC) {

       return (100.0*FREERUN_CORE_C6_RESIDENCY/TSC);
}


/*
   Package c2 residency
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_Package_C2_residency(const double FREERUN_PKG_C2_RESIDENCY,
                                const double TSC) {

       return (100.0*FREERUN_PKG_C2_RESIDENCY/TSC);
}


/*
   Package c3 residency
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_Package_C3_residency(const double FREERUN_PKG_C3_RESIDENCY,
                                const double TSC) {

       return (100.0*FREERUN_PKG_C3_RESIDENCY/TSC);
}


/*
   Package c6 residency
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_Package_C6_residency(const double FREERUN_PKG_C6_RESIDENCY,
                                const double TSC) {

       return (100.0*FREERUN_PKG_C6_RESIDENCY/TSC);
}


/*
    core SW prefetch NTA per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_sw_nta_prefetches(const double SW_PREFETCH_ACCESS_NTA,
                             const double INST_RETIRED_ANY) {

       return(SW_PREFETCH_ACCESS_NTA/INST_RETIRED_ANY);
}


/*
   core  cycles core power throttled.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_Core_power_throttled(const double CORE_POWER_THROTTLED,
                                const double CPU_CLK_UNHALTED_THREAD) {

       return (100.0*CORE_POWER_THROTTLED/
                     CPU_CLK_UNHALTED_THREAD);
}


/*
     Core IPC.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_TMAM_core_ipc(const double INST_RETIRED_ANY,
                         const double CPU_CLK_UNHALTED_THREAD_ANY,
			 const double HW_THREAD_COUNT) {

       const double term1 = CPU_CLK_UNHALTED_THREAD_ANY/
                            HW_THREAD_COUNT;
       return (INST_RETIRED_ANY/term1);
}


/*
    TMAM Info Memory Level Parallelism
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_TMAM_mem_lvl_parallelism(const double L1D_PEND_MISS_PENDING,
                                    const double L1D_PEND_MISS_PENDING_CYCLES_ANY,
				    const double HW_THREAD_COUNT) {

       const double term1 = L1D_PEND_MISS_PENDING_CYCLES_ANY/
                            HW_THREAD_COUNT;
       return (L1D_PEND_MISS_PENDING/term1);
}


/*
   TMAM Info cycles both threads active
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_TMAM_SMT_activity(const double CPU_CLK_THREAD_UNHALTED_ONE_THREAD_ACTIVE,
                             const double CPU_CLK_THREAD_UNHALTED_REF_XCLK_ANY,
			     const double HW_THREAD_COUNT) {

       const double term1 = 1.0-CPU_CLK_THREAD_UNHALTED_ONE_THREAD_ACTIVE/
                            (CPU_CLK_THREAD_UNHALTED_REF_XCLK_ANY*0.5);
       return ((HW_THREAD_COUNT < 2.0) ? 0.0 : 100.0*term1);
}


/*
    TMAM Frontend Bound
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_frontend_bound(const double CPU_CLK_UNHALTED_THREAD,
                          const double IDQ_UOPS_NOT_DELIVERED_CORE){
			
        const double term1 = 4.0*CPU_CLK_UNHALTED_THREAD;
	return (100.0*IDQ_UOPS_NOT_DELIVERED_CORE/term1);
}


/*
   TMAM Frontend_Latency
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_frontend_latency(const double CPU_CLK_UNHALTED_THREAD,
                            const double IDQ_UOPS_NOT_DELIVERED_CYCLES_0_UOPS_DELIV_CORE){
			  
       return (100.0*IDQ_UOPS_NOT_DELIVERED_CYCLES_0_UOPS_DELIV_CORE/
                 CPU_CLK_UNHALTED_THREAD); 
}


/*
    TMAM ICache Misses.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_ICache_misses(const double ICACHE_16B_IFDATA_STALL,
                         const double ICACHE_16B_IFDATA_STALL_c1_e1,
			 const double CPU_CLK_UNHALTED_THREAD) {

       const double term1 = ICACHE_16B_IFDATA_STALL+2.0*
                            ICACHE_16B_IFDATA_STALL_c1_e1;
       return (100.0*term1/CPU_CLK_UNHALTED_THREAD);
}


/*
   ITLB Misses 
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_ITLB_misses(const double ITLB_MISSES_STLB_HIT,
                       const double ITLB_MISSES_WALK_PENDING_c1,
		       const double CPU_CLK_UNHALTED_THREAD) {

       const double term1 = 9.0*ITLB_MISSES_STLB_HIT+
                            ITLB_MISSES_WALK_PENDING_c1;
       return (100.0*term1/CPU_CLK_UNHALTED_THREAD);
}


/*
   Branch Resteers
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_branch_resteers(const double INT_MISC_CLEAR_RESTEER_CYCLES,
                           const double BACLEARS_ANY,
			   const double RS_EVENTS_EMPTY_CYCLES,
			   const double RS_EVENTS_EMPTY_END,
			   const double  ICACHE_16B_IFDATA_STALL,
			   const double ICACHE_64B_IFTAG_STALL,
			   const double CPU_CLK_UNHALTED_THREAD) {
   
       const double term1 = (RS_EVENTS_EMPTY_CYCLES-ICACHE_16B_IFDATA_STALL-
                            ICACHE_64B_IFTAG_STALL)/RS_EVENTS_EMPTY_END;
       return (100.0*(INT_MISC_CLEAR_RESTEER_CYCLES+term1*BACLEARS_ANY)/
                      CPU_CLK_UNHALTED_THREAD);
}


/*
    DSB switches.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_DSB_switches(const double DSB2MITE_SWITCHES_PENALTY_CYCLES,
                        const double CPU_CLK_UNHALTED_THREAD) {

       return (100.0*DSB2MITE_SWITCHES_PENALTY_CYCLES/
                   CPU_CLK_UNHALTED_THREAD);
}


/*
    MS switches.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_MS_switches(const double IDQ_MS_SWITCHES,
                       const double CPU_CLK_UNHALTED_THREAD) {

       return (200.0*IDQ_MS_SWITCHES/CPU_CLK_UNHALTED_THREAD);
}


/*
   Frontend Bandwidth (SMT is enabled)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_smt_frontend_bw(const double CPU_CLK_UNHALTED_THREAD_ANY,
                       const double IDQ_UOPS_NOT_DELIVERED_CORE,
		       const double IDQ_UOPS_NOT_DELIVERED_CYCLES_0_UOPS_DELIV_CORE,
		       const double HW_THREAD_COUNT){
	

        const double term1 = IDQ_UOPS_NOT_DELIVERED_CORE-4.0*
	                     IDQ_UOPS_NOT_DELIVERED_CYCLES_0_UOPS_DELIV_CORE;
        const double term2 = 4.0*(CPU_CLK_UNHALTED_THREAD_ANY/HW_THREAD_COUNT);
        return (100.0*term1/term2);
}


/*
    Frontend Bandwidth (SMT is disabled)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_no_smt_frontend_bw(const double CPU_CLK_UNHALTED_THREAD_ANY,
                              const double IDQ_UOPS_NOT_DELIVERED_CORE,
		              const double IDQ_UOPS_NOT_DELIVERED_CYCLES_0_UOPS_DELIV_CORE) {

       const double term1 = IDQ_UOPS_NOT_DELIVERED_CORE-4.0*
	                    IDQ_UOPS_NOT_DELIVERED_CYCLES_0_UOPS_DELIV_CORE;
       const double term2 = 4.0*CPU_CLK_UNHALTED_THREAD_ANY;
        return (100.0*term1/term2);    
}


/*
    Bad Speculation (SMT enabled).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_smt_bad_speculation(const double UOPS_ISSUED_ANY,
			   const double UOPS_RETIRED_RETIRED_SLOTS,
			   const double INT_MISC_RECOVERY_CYCLES_ANY,
			   const double CPU_CLK_UNHALTED_THREAD_ANY,
			   const double HW_THREAD_COUNT) {

       const double term1 = (UOPS_ISSUED_ANY-UOPS_RETIRED_RETIRED_SLOTS)+
                            4.0*INT_MISC_RECOVERY_CYCLES_ANY/HW_THREAD_COUNT;
       const double term2 = 4.0*CPU_CLK_UNHALTED_THREAD_ANY/HW_THREAD_COUNT;
       return (100.0*term1/term2);
}


/*
    Bad Speculation (SMT disabled).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_no_smt_bad_speculation(const double UOPS_ISSUED_ANY,
			          const double UOPS_RETIRED_RETIRED_SLOTS,
			          const double INT_MISC_RECOVERY_CYCLES_ANY,
			          const double CPU_CLK_UNHALTED_THREAD) {

        const double term1 = (UOPS_ISSUED_ANY-UOPS_RETIRED_RETIRED_SLOTS)+
                             4.0*INT_MISC_RECOVERY_CYCLES_ANY;
	const double term2 = 4.0*CPU_CLK_UNHALTED_THREAD;
	return (100.0*term1/term2);
}


/*
    Branch Mispredicts (STM enabled).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_smt_branch_mispredicts(const double UOPS_ISSUED_ANY,
			          const double UOPS_RETIRED_RETIRED_SLOTS,
			          const double INT_MISC_RECOVERY_CYCLES_ANY,
			          const double CPU_CLK_UNHALTED_THREAD_ANY,
				  const double BR_MISP_RETIRED_ALL_BRANCHES,
				  const double MACHINE_CLEARS_COUNT,
				  const double HW_THREAD_COUNT) {

        const double term1 =  BR_MISP_RETIRED_ALL_BRANCHES/(BR_MISP_RETIRED_ALL_BRANCHES+
	                                                    MACHINE_CLEARS_COUNT);
	const double term2 =  (UOPS_ISSUED_ANY-UOPS_RETIRED_RETIRED_SLOTS)+
                               4.0*INT_MISC_RECOVERY_CYCLES_ANY/HW_THREAD_COUNT;
	const double term3 =  4.0*CPU_CLK_UNHALTED_THREAD_ANY/HW_THREAD_COUNT;
	return (term1*100.0*term2/term3);
}


/*
   Branch Misprediction (SMT disabled).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_no_smt_branch_misprediction(const double UOPS_ISSUED_ANY,
			          const double UOPS_RETIRED_RETIRED_SLOTS,
			          const double INT_MISC_RECOVERY_CYCLES_ANY,
			          const double CPU_CLK_UNHALTED_THREAD,
				  const double BR_MISP_RETIRED_ALL_BRANCHES,
				  const double MACHINE_CLEARS_COUNT) {

        const double term1 =  BR_MISP_RETIRED_ALL_BRANCHES/(BR_MISP_RETIRED_ALL_BRANCHES+
	                                                    MACHINE_CLEARS_COUNT);
	const double term2 =  (UOPS_ISSUED_ANY-UOPS_RETIRED_RETIRED_SLOTS)+
                               4.0*INT_MISC_RECOVERY_CYCLES_ANY;
	const double term3 =  4.0*CPU_CLK_UNHALTED_THREAD;
	return (term1*100.0*term2/term3);  
}


/*
    Machine Clears (SMT enabled).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_smt_machine_clears(const double UOPS_ISSUED_ANY,
			      const double UOPS_RETIRED_RETIRED_SLOTS,
			      const double INT_MISC_RECOVERY_CYCLES_ANY,
			      const double CPU_CLK_UNHALTED_THREAD_ANY,
			      const double BR_MISP_RETIRED_ALL_BRANCHES,
			      const double MACHINE_CLEARS_COUNT,
			      const double HW_THREAD_COUNT) {

        const double term1 =  MACHINE_CLEARS_COUNT/(BR_MISP_RETIRED_ALL_BRANCHES+
	                                                    MACHINE_CLEARS_COUNT);
	const double term2 =  (UOPS_ISSUED_ANY-UOPS_RETIRED_RETIRED_SLOTS)+
                               4.0*INT_MISC_RECOVERY_CYCLES_ANY/HW_THREAD_COUNT;
	const double term3 =  4.0*CPU_CLK_UNHALTED_THREAD_ANY/HW_THREAD_COUNT;
	return (term1*100.0*term2/term3);
}


/*
    Machine Clears (SMT disabled).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_no_smt_machine_clears(const double UOPS_ISSUED_ANY,
			          const double UOPS_RETIRED_RETIRED_SLOTS,
			          const double INT_MISC_RECOVERY_CYCLES_ANY,
			          const double CPU_CLK_UNHALTED_THREAD,
				  const double BR_MISP_RETIRED_ALL_BRANCHES,
				  const double MACHINE_CLEARS_COUNT) {

        const double term1 =  MACHINE_CLEARS_COUNT/(BR_MISP_RETIRED_ALL_BRANCHES+
	                                                    MACHINE_CLEARS_COUNT);
	const double term2 =  (UOPS_ISSUED_ANY-UOPS_RETIRED_RETIRED_SLOTS)+
                               4.0*INT_MISC_RECOVERY_CYCLES_ANY;
	const double term3 =  4.0*CPU_CLK_UNHALTED_THREAD;
	return (term1*100.0*term2/term3);
}


/*
   Backend Bound (SMT enabled)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_smt_backend_bound(const double IDQ_UOPS_NOT_DELIVERED_CORE,
                             const double UOPS_ISSUED_ANY,
			     const double INT_MISC_RECOVERY_CYCLES_ANY,
			     const double CPU_CLK_UNHALTED_THREAD_ANY,
			     const double UOPS_RETIRED_RETIRE_SLOTS,
			     const double HW_THREAD_COUNT) {

        const double term1 = UOPS_ISSUED_ANY-UOPS_RETIRED_RETIRE_SLOTS+4.0*
	                     (INT_MISC_RECOVERY_CYCLES_ANY/HW_THREAD_COUNT);
	const double term2 = IDQ_UOPS_NOT_DELIVERED_CORE+UOPS_RETIRED_RETIRE_SLOTS;
	const double term3 = 4.0*CPU_CLK_UNHALTED_THREAD_ANY/HW_THREAD_COUNT;
	return (100.0-(100.0*(term1+term2)/term3));
}


/*
    Backend Boud (SMT disabled).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_no_smt_backend_bound(const double IDQ_UOPS_NOT_DELIVERED_CORE,
                                const double UOPS_ISSUED_ANY,
			        const double INT_MISC_RECOVERY_CYCLES_ANY,
			        const double CPU_CLK_UNHALTED_THREAD,
			        const double UOPS_RETIRED_RETIRE_SLOTS) {

        const double term1 = UOPS_ISSUED_ANY-UOPS_RETIRED_RETIRE_SLOTS+4.0*
	                     INT_MISC_RECOVERY_CYCLES_ANY;
	const double term2 = IDQ_UOPS_NOT_DELIVERED_CORE+UOPS_RETIRED_RETIRE_SLOTS;
	const double term3 = 4.0*CPU_CLK_UNHALTED_THREAD;
	return (100.0-(100.0*(term1+term2)/term3));  
}


/*
    L1 Bound
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L1D_bound(const double CYCLE_ACTIVITY_STALLS_MEM_ANY,
                    const double CYCLE_ACTIVITY_STALLS_L1D_MISS,
		    const double CPU_CLK_UNHALTED_THREAD) {

        return (100.0*(CYCLE_ACTIVITY_STALLS_MEM_ANY-
	               CYCLE_ACTIVITY_STALLS_L1D_MISS)/
		       CPU_CLK_UNHALTED_THREAD);
}


/*
   DTLB Load
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_DTLB_load(const double DTLB_LOAD_MISSES_STLB_HIT,
                     const double DTLB_LOAD_MISSES_WALK_ACTIVE,
		     const double CPU_CLK_UNHALTED_THREAD) {

       const double term1 = 7.0*DTLB_LOAD_MISSES_STLB_HIT+
                            DTLB_LOAD_MISSES_WALK_ACTIVE;
       return (100.0*term1/CPU_CLK_UNHALTED_THREAD);
}


/*
    Stores forward blocked.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_stores_fwd_blocked(const double LD_BLOCKS_STORE_FORWARD,
                              const double CPU_CLK_UNHALTED_THREAD) {

       return (100.0*(13.0*LD_BLOCKS_STORE_FORWARD)/
                      CPU_CLK_UNHALTED_THREAD);
}

#include <algorithm>


/*
    Lock latency.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_lock_latency(const double MEM_INST_RETIRED_LOCK_LOADS,
                        const double CPU_CLK_UNHALTED_THREAD,
			const double MEM_INST_RETIRED_ALL_STORES,
			const double OFFCORE_REQUESTS_OUTSTANDING_CYCLES_WITH_DEMAND_RFO) {

        const double term1 = MEM_INST_RETIRED_LOCK_LOADS/
	                     MEM_INST_RETIRED_ALL_STORES;
	const double term2 = std::min(CPU_CLK_UNHALTED_THREAD,
	                              OFFCORE_REQUESTS_OUTSTANDING_CYCLES_WITH_DEMAND_RFO);
	return (100.0*(term1*term2)/CPU_CLK_UNHALTED_THREAD);
}


/*
   L2 Bounds
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L2_bounds(const double CYCLE_ACTIVITY_STALLS_L1D_MISS,
                     const double CYCLE_ACTIVITY_STALLS_L2_MISS,
		     const double CPU_CLK_UNHALTED_THREAD) {

       return (100.0*(CYCLE_ACTIVITY_STALLS_L1D_MISS-
                      CYCLE_ACTIVITY_STALLS_L2_MISS)/
		      CPU_CLK_UNHALTED_THREAD);
}


/*
   L3 Bounds.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_bounds(const double CYCLE_ACTIVITY_STALLS_L2_MISS,
                     const double CYCLE_ACTIVITY_STALLS_L3_MISS,
		     const double CPU_CLK_UNHALTED_THREAD) {

       return (100.0*(CYCLE_ACTIVITY_STALLS_L2_MISS-
                      CYCLE_ACTIVITY_STALLS_L3_MISS)/
		      CPU_CLK_UNHALTED_THREAD);
}


/*
    Contested Accesses
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_contested_accesses(const double MEM_LOAD_L3_HIT_RETIRED_XSNP_HITM,
                              const double CPU_CLK_UNHALTED_THREAD,
			      const double MEM_LOAD_L3_HIT_RETIRED_XSNP_MISS) {

       const double term1 =  MEM_LOAD_L3_HIT_RETIRED_XSNP_HITM+
                              MEM_LOAD_L3_HIT_RETIRED_XSNP_MISS;
       return (100.0*(60.0*term1/CPU_CLK_UNHALTED_THREAD));
}


/*
    Data Sharing
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_data_sharing(const double MEM_LOAD_L3_HIT_RETIRED_XSNP_HIT,
                        const double CPU_CLK_UNHALTED_THREAD) {

       const double term1 = 43.0*MEM_LOAD_L3_HIT_RETIRED_XSNP_HIT/
                            CPU_CLK_UNHALTED_THREAD;
       return (100.0*term1);
}


/*
   L3 latency.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_latency(const double OFFCORE_REQUESTS_OUTSTANDING_DEMAND_DATA_RD_GE_6,
                      const double OFFCORE_REQUESTS_OUTSTANDING_L3_MISS_DEMAND_DATA_RD_GE_6,
		      const double CPU_CLK_UNHALTED_THREAD,
		      const double OFFCORE_REQUESTS_OUTSTANDING_CYCLES_WITH_DEMAND_DATA_RD,
		      const double OFFCORE_REQUESTS_OUTSTANDING_CYCLES_WITH_L3_MISS_DEMAND_DATA_RD) {

        const double term1 = (std::min(OFFCORE_REQUESTS_OUTSTANDING_CYCLES_WITH_DEMAND_DATA_RD,
	                              CPU_CLK_UNHALTED_THREAD)-std::min(OFFCORE_REQUESTS_OUTSTANDING_CYCLES_WITH_L3_MISS_DEMAND_DATA_RD,
				      CPU_CLK_UNHALTED_THREAD))/CPU_CLK_UNHALTED_THREAD;
	const double term2 = (std::min(OFFCORE_REQUESTS_OUTSTANDING_DEMAND_DATA_RD_GE_6,
	                               CPU_CLK_UNHALTED_THREAD)-std::min(OFFCORE_REQUESTS_OUTSTANDING_L3_MISS_DEMAND_DATA_RD_GE_6,
				       CPU_CLK_UNHALTED_THREAD))/CPU_CLK_UNHALTED_THREAD;
	return (100*(term1-term2));
}


/*
   L3 Bandwidth
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_bandwidth(const double OFFCORE_REQUESTS_OUTSTANDING_DEMAND_DATA_RD_GE_6,
                        const double OFFCORE_REQUESTS_OUTSTANDING_L3_MISS_DEMAND_DATA_RD_GE_6,
			const double CPU_CLK_UNHALTED_THREAD) {

       const double term1 = (std::min(OFFCORE_REQUESTS_OUTSTANDING_DEMAND_DATA_RD_GE_6,
                                      CPU_CLK_UNHALTED_THREAD)-std::min(OFFCORE_REQUESTS_OUTSTANDING_L3_MISS_DEMAND_DATA_RD_GE_6,
				      CPU_CLK_UNHALTED_THREAD))/CPU_CLK_UNHALTED_THREAD;
       return (100.0*term1);
}


/*
    SuperQueue Full (SMT enabled)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_smt_SQ_full(const double CPU_CLK_UNHALTED_THREAD_ANY,
                       const double OFFCORE_REQUESTS_BUFFER_SQ_FULL,
		       const double HW_THREAD_COUNT) {

       const double term1 = OFFCORE_REQUESTS_BUFFER_SQ_FULL/
                            HW_THREAD_COUNT;
       const double term2 = CPU_CLK_UNHALTED_THREAD_ANY/
                            HW_THREAD_COUNT;
       return (100.0*term1/term2);
}


/*
    SuperQueue Full (SMT disabled)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_no_smt_SQ_full(const double CPU_CLK_UNHALTED_THREAD,
                       const double OFFCORE_REQUESTS_BUFFER_SQ_FULL){
	 
       return (100.0*OFFCORE_REQUESTS_BUFFER_SQ_FULL/
                           CPU_CLK_UNHALTED_THREAD);
}


/*
     Memory bound.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_memory_bound(const double CYCLE_ACTIVITY_STALLS_L3_MISS,
                        const double CPU_CLK_UNHALTED_THREAD) {

       return (100.0*CYCLE_ACTIVITY_STALLS_L3_MISS/
                      CPU_CLK_UNHALTED_THREAD);
}


/*
    Memory bandwidth.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_memory_bandwidth(const double OFFCORE_REQUESTS_OUTSTANDING_L3_MISS_DEMAND_DATA_RD_GE_6,
                            const double CPU_CLK_UNHALTED_THREAD) {

       return (100.0*std::min(OFFCORE_REQUESTS_OUTSTANDING_L3_MISS_DEMAND_DATA_RD_GE_6,
                              CPU_CLK_UNHALTED_THREAD)/CPU_CLK_UNHALTED_THREAD);
}


/*
   Memory latency
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_memory_latency(const double OFFCORE_REQUESTS_OUTSTANDING_L3_MISS_DEMAND_DATA_RD_GE_6,
                          const double OFFCORE_REQUESTS_OUTSTANDING_CYCLES_WITH_L3_MISS_DEMAND_DATA_RD,
			  const double CPU_CLK_UNHALTED_THREAD) {

       const double term1 = (std::min(OFFCORE_REQUESTS_OUTSTANDING_CYCLES_WITH_L3_MISS_DEMAND_DATA_RD,
                                      CPU_CLK_UNHALTED_THREAD)-std::min(OFFCORE_REQUESTS_OUTSTANDING_L3_MISS_DEMAND_DATA_RD_GE_6,
				      CPU_CLK_UNHALTED_THREAD))/CPU_CLK_UNHALTED_THREAD;
       return (100.0*term1);
}


/*
   Stores bound.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_stores_bound(const double EXE_ACTIVITY_BOUND_ON_STORES,
                        const double CPU_CLK_UNHALTED_THREAD) {

       return (100.0*( EXE_ACTIVITY_BOUND_ON_STORES/
                       CPU_CLK_UNHALTED_THREAD));
}


/*
   DTLB stores (SMT enabled).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_smt_DTLB_stores(const double DTLB_STORE_MISSES_STLB_HIT,
                       const double DTLB_STORE_MISSES_WALK_ACTIVE,
		       const double CPU_CLK_UNHALTED_THREAD_ANY,
		       const double HW_THREAD_COUNT) {

       const double term1 = 7.0*DTLB_STORE_MISSES_STLB_HIT+
                             DTLB_STORE_MISSES_WALK_ACTIVE;
       const double term2 = CPU_CLK_UNHALTED_THREAD_ANY/
                            HW_THREAD_COUNT;
       return (100.0*term1/term2);
}


/*
    DTLB stores (SMT disabled).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_no_smt_DTLB_stores(const double DTLB_STORE_MISSES_STLB_HIT,
                       const double DTLB_STORE_MISSES_WALK_ACTIVE,
		       const double CPU_CLK_UNHALTED_THREAD){
		    

       const double term1 = 7.0*DTLB_STORE_MISSES_STLB_HIT+
                             DTLB_STORE_MISSES_WALK_ACTIVE;
       return (100.0*term1/CPU_CLK_UNHALTED_THREAD);
}


/*
    Divider.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_Divider(const double ARITH_DIVIDER_ACTIVE,
                   const double CPU_CLK_UNHALTED_THREAD) {

       return (100.0*ARITH_DIVIDER_ACTIVE/
                     CPU_CLK_UNHALTED_THREAD);
}


/*
   Ports Utilization
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_Ports_utilization(const double EXE_ACTIVITY_EXE_BOUND_0_PORTS,
                             const double EXE_ACTIVITY_1_PORTS_UTIL,
			     const double CYCLE_ACTIVITY_STALLS_MEM_ANY,
			     const double EXE_ACTIVITY_BOUND_ON_STORES,
			     const double EXE_ACTIVITY_2_PORTS_UTIL,
			     const double CPU_CLK_UNHALTED_THREAD,
			     const double INST_RETIRED_ANY) {

       double term1 = 0.0;
       if(INST_RETIRED_ANY/
          CPU_CLK_UNHALTED_THREAD > 1.8) {
          term1 = EXE_ACTIVITY_2_PORTS_UTIL;
       }
       const double term2 = EXE_ACTIVITY_EXE_BOUND_0_PORTS+
                            EXE_ACTIVITY_1_PORTS_UTIL+
			    term1+CYCLE_ACTIVITY_STALLS_MEM_ANY+
			    EXE_ACTIVITY_BOUND_ON_STORES;
       const double term3 = CYCLE_ACTIVITY_STALLS_MEM_ANY-
                            EXE_ACTIVITY_BOUND_ON_STORES;
       return (100.0*(term2-term3/CPU_CLK_UNHALTED_THREAD));
}


/*
     Port_0 Utilized.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_Port_0_utilized(const double UOPS_EXECUTED_CORE_CYCLES_NONE,
                               const double CPU_CLK_UNHALTED_THREAD_ANY,
			       const double EXE_ACTIVITY_EXE_BOUND_0_PORTS,
			       const double HW_THREAD_COUNT) {

        double term1 = 0.0;
	HW_THREAD_COUNT > 1.0 ? term1 = UOPS_EXECUTED_CORE_CYCLES_NONE/2.0 :
	                                term1 = EXE_ACTIVITY_EXE_BOUND_0_PORTS;
        return (100.0*(term1/(CPU_CLK_UNHALTED_THREAD_ANY/HW_THREAD_COUNT)));
}


/*
   Port_1 Utilized
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_Port_1_utilized(const double UOPS_EXECUTED_CORE_CYCLES_GE_1,
                           const double CPU_CLK_UNHALTED_THREAD_ANY,
			   const double UOPS_EXECUTED_CORE_CYCLES_GE_2,
			   const double EXE_ACTIVITY_1_PORTS_UTIL,
			   const double HW_THREAD_COUNT) {

       double term1 = 0.0;
       HW_THREAD_COUNT > 1.0 ? term1 = (UOPS_EXECUTED_CORE_CYCLES_GE_1-
                                        UOPS_EXECUTED_CORE_CYCLES_GE_2)/2.0 :
					              term1 =  EXE_ACTIVITY_1_PORTS_UTIL;
       return (100.0*(term1/(CPU_CLK_UNHALTED_THREAD_ANY/HW_THREAD_COUNT)));
}


/*
   Port_2 Utilized.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_Port_2_utilized(const double UOPS_EXECUTED_CORE_CYCLES_GE_2,
                           const double CPU_CLK_UNHALTED_THREAD_ANY,
			   const double UOPS_EXECUTED_CORE_CYCLES_GE_3,
			   const double EXE_ACTIVITY_2_PORTS_UTIL,
			   const double HW_THREAD_COUNT) {

       double term1 = 0.0;
       HW_THREAD_COUNT > 1.0 ? term1 = (UOPS_EXECUTED_CORE_CYCLES_GE_2-
                                        UOPS_EXECUTED_CORE_CYCLES_GE_3)/2.0 :
					              term1 =  EXE_ACTIVITY_2_PORTS_UTIL;
       return (100.0*(term1/(CPU_CLK_UNHALTED_THREAD_ANY/HW_THREAD_COUNT)));
}


/*
   Port_3m Utilized.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_Port_3m_utilized(const double UOPS_EXECUTED_CORE_CYCLES_GE_3,
                            const double CPU_CLK_UNHALTED_THREAD_ANY,
			    const double HW_THREAD_COUNT) {

       const double term1 = UOPS_EXECUTED_CORE_CYCLES_GE_3/
                            HW_THREAD_COUNT;
       const double term2 = CPU_CLK_UNHALTED_THREAD_ANY/
                            HW_THREAD_COUNT;
       return (100.0*(term1/term2));
}


/*
   Retiring (SMT enabled)
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_smt_retiring(const double CPU_CLK_UNHALTED_THREAD_ANY,
                        const double UOPS_RETIRED_RETIRE_SLOTS,
			const double HW_THREAD_COUNT) {

       const double term1 = 4.0*(CPU_CLK_UNHALTED_THREAD_ANY/
                                            HW_THREAD_COUNT);
       return (100.0*(UOPS_RETIRED_RETIRE_SLOTS/term1));
}


/*
   Retiring (SMT disabled).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_no_smt_retiring(const double CPU_CLK_UNHALTED_THREAD,
                           const double UOPS_RETIRED_RETIRE_SLOTS) {

       const double term1 = 4.0*CPU_CLK_UNHALTED_THREAD;
       return (100.0*(UOPS_RETIRED_RETIRE_SLOTS/term1));
}


/*
    Basic activit (SMT enabled).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_smt_base_activity(const double CPU_CLK_UNHALTED_THREAD_ANY,
                             const double UOPS_RETIRED_RETIRE_SLOTS,
			     const double IDQ_MS_UOPS,
			     const double UOPS_ISSUED_ANY,
			     const double HW_THREAD_COUNT) {

        const double term1 = 4.0*(CPU_CLK_UNHALTED_THREAD_ANY/
	                                   HW_THREAD_COUNT);
	const double term2 = UOPS_RETIRED_RETIRE_SLOTS/UOPS_ISSUED_ANY;
	return (100.0*(UOPS_RETIRED_RETIRE_SLOTS/term1)-(term2*IDQ_MS_UOPS/term1));
}


/*
    Basic activit (SMT disabled).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_no_smt_base_activity(const double CPU_CLK_UNHALTED_THREAD,
                             const double UOPS_RETIRED_RETIRE_SLOTS,
			     const double IDQ_MS_UOPS,
			     const double UOPS_ISSUED_ANY) {
			    

        const double term1 = 4.0*CPU_CLK_UNHALTED_THREAD;
	const double term2 = UOPS_RETIRED_RETIRE_SLOTS/UOPS_ISSUED_ANY;
	return (100.0*(UOPS_RETIRED_RETIRE_SLOTS/term1)-(term2*IDQ_MS_UOPS/term1));
}


/*
    FP scalar retiring fraction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_fp_scalar_fraction(const double UOPS_RETIRED_RETIRE_SLOTS,
                              const double FP_ARITH_INST_RETIRED_SCALAR_SINGLE,
			      const double FP_ARITH_INST_RETIRED_SCALAR_DOUBLE) {

       return (100.0*(FP_ARITH_INST_RETIRED_SCALAR_SINGLE+
                      FP_ARITH_INST_RETIRED_SCALAR_DOUBLE)/
		       UOPS_RETIRED_RETIRE_SLOTS);
}


/*
    FP vector retiring fraction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_fp_vector_fraction(const double UOPS_RETIRED_RETIRE_SLOTS,
                              const double FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE,
			      const double FP_ARITH_INST_RETIRED_128B_PACKED_SINGLE,
			      const double FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE,
			      const double FP_ARITH_INST_RETIRED_256B_PACKED_SINGLE,
			      const double FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE,
			      const double FP_ARITH_INST_RETIRED_512B_PACKED_SINGLE) {

        const double term1 = FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE+
	                     FP_ARITH_INST_RETIRED_128B_PACKED_SINGLE+
			     FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE;
	const double term2 = FP_ARITH_INST_RETIRED_256B_PACKED_SINGLE+
	                     FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE+
			     FP_ARITH_INST_RETIRED_512B_PACKED_SINGLE;
	return (100.0*(term1+term2)/UOPS_RETIRED_RETIRE_SLOTS);
}


/*
     Microcode Sequencer (SMT enabled).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_smt_MS(const double CPU_CLK_UNHALTED_THREAD_ANY,
                  const double UOPS_RETIRED_RETIRE_SLOTS,
		  const double IDQ_MS_UOPS,
		  const double UOPS_ISSUED_ANY,
		  const double HW_THREAD_COUNT) {

       const double term1 = IDQ_MS_UOPS/(4.0*(CPU_CLK_UNHALTED_THREAD_ANY/
                                              HW_THREAD_COUNT));
       return (100.0*(UOPS_RETIRED_RETIRE_SLOTS/UOPS_ISSUED_ANY)*term1);
}


/*
     Microcode Sequencer (SMT disabled).
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_no_smt_MS(const double CPU_CLK_UNHALTED_THREAD,
                  const double UOPS_RETIRED_RETIRE_SLOTS,
		  const double IDQ_MS_UOPS,
		  const double UOPS_ISSUED_ANY) {
	

       const double term1 = IDQ_MS_UOPS/(4.0*CPU_CLK_UNHALTED_THREAD);
                                             
       return (100.0*(UOPS_RETIRED_RETIRE_SLOTS/UOPS_ISSUED_ANY)*term1);
}


/*
    LLC Local code/data reads hitting in S state in snoop filter per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_code_data_read_S_hit(const double UNC_CHA_LLC_LOOKUP_DATA_READ_state_0x2,
                                   const double INST_RETIRED_ANY) {

       return (UNC_CHA_LLC_LOOKUP_DATA_READ_state_0x2/
                                    INST_RETIRED_ANY);
}


/*
     LLC Local code/data reads hitting in E state in snoop filter per instruction.  
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_code_data_read_E_hit(const double UNC_CHA_LLC_LOOKUP_DATA_READ_state_0x4,
                                   const double INST_RETIRED_ANY) {

        return (UNC_CHA_LLC_LOOKUP_DATA_READ_state_0x4/
	                           INST_RETIRED_ANY);
}


/*
     LLC Local code/data reads hitting in M/E/F states in LLC per instruction.
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_L3_code_data_read_MEF_hit(const double UNC_CHA_LLC_LOOKUP_DATA_READ_state_0xE0,
                                     const double INST_RETIRED_ANY) {

        return (UNC_CHA_LLC_LOOKUP_DATA_READ_state_0xE0/
	                                INST_RETIRED_ANY);
}


/*
     INT_MISC.RECOVERY_CYCLES ratio
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_int_misc_recovery_cycles_ratio(const double INT_MISC_RECOVERY_CYCLES,
                                          const double CPU_CLK_UNHALTED_THREAD) {

       return (100.0*( INT_MISC_RECOVERY_CYCLES/
                           CPU_CLK_UNHALTED_THREAD));
}


/*
   INT_MISC.CLEAR_RESTEER_CYCLES ratio
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_int_misc_clear_resteers_cycles_ratio(const double  INT_MISC_CLEAR_RESTEER_CYCLES,
                                                const double CPU_CLK_UNHALTED_THREAD) {

       return (100.0*( INT_MISC_CLEAR_RESTEER_CYCLES/
                           CPU_CLK_UNHALTED_THREAD));
}


/*
   RS_EVENTS.EMPTY_CYCLES ratio
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_rs_events_empty_cycles_ratio(const double RS_EVENTS_EMPTY_CYCLES,
                                        const double CPU_CLK_UNHALTED_THREAD) {

       return (100.0*(RS_EVENTS_EMPTY_CYCLES/
                      CPU_CLK_UNHALTED_THREAD));
}


/*
    BR_INST_RETIRED.ALL_BRANCHES ratio
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_br_inst_retired_all_branches_ratio(const double  BR_INST_RETIRED_ALL_BRANCHES,
                                              const double  INST_RETIRED_ANY) {

       return (100.0*(BR_INST_RETIRED_ALL_BRANCHES/
                      INST_RETIRED_ANY));
}


/*
    BR_INST_RETIRED.CONDITIONAL
*/
__attribute__((always_inline))
__attribute__((pure))
__attribute__((aligned(32)))
#pragma omp declare simd
static inline
double skx_br_inst_retired_conditional_ratio(const double BR_INST_RETIRED_CONDITIONAL,
                                             const double BR_INST_RETIRED_ALL_BRANCHES) {

       return (100.0*(BR_INST_RETIRED_CONDITIONAL/
                      BR_INST_RETIRED_ALL_BRANCHES));
}


/*
*/

#endif /*__GMS_SKX_CLK_HW_EVENTS_METRICS_HPP__*/
