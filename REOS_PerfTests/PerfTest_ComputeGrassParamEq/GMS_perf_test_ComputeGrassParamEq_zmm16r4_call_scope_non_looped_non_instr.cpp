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



void perf_test_ComputeGrassParamEq_zmm16r4_call_scope_non_looped_non_instr() {

     // By Intel Intrinsic Guide (SKX CPU)
     constexpr int32_t rdtscp_latency = 42;
     char msr_file_name[64];
     uint64_t core_counters_start[4] = {};
     uint64_t core_counters_end[4]   = {};
     uint64_t data[4] = {};
     uint32_t regs[4] = {0x186,0x187,0x188,0x189};
     uint64_t count;
     uint64_t fix_ins_start,fix_ins_end;
     uint64_t gen_cyc_start,gen_cyc_end;
     uint64_t gen_ref_start,gen_ref_end;
     uint64_t tsc_start,tsc_end;
     uint64_t lvl2_cyc_start,lvl2_cyc_end;
     uint64_t core_halt_delta;
     uint64_t volatile dummy0;
     uint64_t volatile dummy1;
     uint64_t volatile dummy2;
     double utilization,avg_ghz;
     float nominal_ghz;
     int32_t fd;
     int32_t cpu;
     int32_t chip,core;
     int32_t counter;
     int32_t core_ctr_width;
     int32_t fix_ctr_width;
     
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
     nominal_ghz = get_TSC_frequency()/1.0e9;
     //printf("Nominal frequency (Ghz): %f\n",nominal_ghz);
     core_ctr_width = get_core_counter_width();
     fix_ctr_width  = get_fix_counter_width();
     //printf("Programmable counters width: %d bits\n",core_ctr_width);
     //printf("Fixed counters width:        %d bits\n",fix_ctr_width);
     if(core_ctr_width != fix_ctr_width) {
        printf("Warning -- programmable counters width does not match the fixed counters width!!\n");
     }
       // Set core affinity
     cpu = 3;
     cpu_set_t cpu_set;
     CPU_ZERO(&cpu_set);
     CPU_SET(cpu,&cpu_set);
     if(sched_setaffinity(0,sizeof(cpu_set),&cpu_set) < 0) {
        printf("***FATAL*** -- Failed to set core affinity: errno=%d\n",errno);
	exit(EXIT_FAILURE);
     }
     count = full_rdtscp(&chip,&core);
     // Call ctor as close as possible to tested function call.
     GrassScattererAVX512 ScattererZMM16r4_1 = GrassScattererAVX512(nsteps,
                                                             ordinal,
							     param_npts,
							     lat,
							     lon,
							     elev,
							     ceps,
							     mlock);
     dummy0 = rdpmc(0);
     // Begin of test block
     //dummy1 = get_core_counter_width();
     // Forcing core to fetch and decode ahead of first usage a lfence
     // instructions.
     __asm__ __volatile__ ("lfence");
     for(int32_t i = 0; i != 4; ++i) {
         core_counters_start[i] = rdpmc(i);
     }
     lvl2_cyc_start = rdpmc(4);
     fix_ins_start = rdpmc_instructions();
     gen_cyc_start = rdpmc_actual_cycles();
     gen_ref_start = rdpmc_reference_cycles();
     tsc_start = rdtscp();
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
     //dummy2 = get_core_counter_width();
     for(int32_t i = 0; i != 4; ++i) {
         core_counters_end[i] = rdpmc(i);
     }
     lvl2_cyc_end = rdpmc(4);
     fix_ins_end  = rdpmc_instructions();
     gen_cyc_end  = rdpmc_actual_cycles();
     gen_ref_end  = rdpmc_reference_cycles();
     tsc_end      = rdtscp();
     __asm__ __volatile__("lfence");
     // Measurements end.
     core_halt_delta = (tsc_end-tsc_start-rdtscp_latency)-(gen_ref_end-gen_ref_start);
    
     // Core utilization
     utilization = (double)(corrected_pmc_delta(gen_ref_end,gen_ref_start,core_ctr_width)) /
                           (tsc_end-tsc_start-rdtscp_latency);
     avg_ghz     = (double)(corrected_pmc_delta(gen_ref_end,gen_ref_start,core_ctr_width)) /
                           (tsc_end-tsc_start-rdtscp_latency)*nominal_ghz;
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
     printf("Core utilization: %f,  Core average frequency: %f\n",utilization,avg_ghz);
     printf("Unhalt-to-Halt discrepancy: %lu\n",core_halt_delta);
     printf("L2 cycles: %lu\n", lvl2_cyc_end-lvl_cyc_start);
     if(fd >= 0) {
        printf(" 0x186: 0x%08, 0x187: 0x%08, 0x188: 0x%08, 0x189: 0x%08\n",data[0],data[1],data[2],data[3]);
     }
     for(int32_t i = 0; i != 4; ++i) {
         printf(" %lu ", corrected_pmc_delta(core_counters_end[i],core_counters_start[i],core_ctr_width));
     }
}



int main(int argc, char * argv[]) {

    perf_test_ComputeGrassParamEq_zmm16r4_call_scope_non_looped_non_instr()
    return (0);
}
