
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <cpuid.h>
#include <iostream>
#include <iomanip>
#include <string>
#include "GMS_skx_cpu_operating_freq_metric.h"

/*
   icpc -o unit_test_skx_cpu_operating_freq -fp-model fast=2 -ftz -std=c++17 -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_skx_clk_hw_events_metrics.hpp GMS_skx_cpu_operating_freq_metric.h GMS_skx_cpu_operating_freq_metric.cpp unit_test_skx_cpu_operating_freq.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 -fopenmp GMS_config.h GMS_skx_clk_hw_events_metrics.hpp GMS_skx_cpu_operating_freq_metric.h GMS_skx_cpu_operating_freq_metric.cpp unit_test_skx_cpu_operating_freq.cpp

*/

// Cross-platform compiler intrinsic configuration for CPUID

static inline
void call_cpuid(uint32_t leaf, uint32_t registers[4]) 
{
    // GCC/Clang intrinsic populates EAX, EBX, ECX, EDX
    __get_cpuid(leaf, &registers[0], &registers[1], &registers[2], &registers[3]);
}

__ATTR_HOT__
__ATTR_ALIGN__(32)
void show_intel_cpu_info(std::uint32_t cpu_uarch_basic[7])
{
    uint32_t regs[4] = {0};

    // Leaf 1 fetches processor signature (EAX) and features
    call_cpuid(1, regs);
    uint32_t eax = regs[0];

    // Decode fields using standard Intel bitmasks
    uint32_t base_family     = (eax >> 8) & 0xF;
    uint32_t base_model      = (eax >> 4) & 0xF;
    uint32_t stepping        = eax & 0xF;
    uint32_t extended_model  = (eax >> 16) & 0xF;
    uint32_t extended_family = (eax >> 20) & 0xFF;

    // Calculate final Display Family
    // If base_family is 0xF, Display Family = base_family + extended_family
    uint32_t display_family = (base_family == 0xF) ? (base_family + extended_family) : base_family;

    // Calculate final Display Model
    // For Intel processors, if display_family is 6 or 15, 
    // Display Model = (extended_model << 4) | base_model
    uint32_t display_model = base_model;
    if (display_family == 6 || display_family == 15) {
        display_model = (extended_model << 4) | base_model;
    }

    // Output raw bit allocations and translated numbers
    std::cout << "====================================================\n";
    std::cout << "      INTEL CPUID DECODER (MICROARCHITECTURE)       \n";
    std::cout << "====================================================\n";
    std::cout << "Raw EAX Signature: 0x" << std::hex << std::uppercase << std::setw(8) << std::setfill('0') << eax << std::dec << "\n\n";
    
    std::cout << "Parsed Raw Fields:\n";
    std::cout << "  - Base Family:     " << base_family << "\n";
    std::cout << "  - Base Model:      " << base_model << " (0x" << std::hex << base_model << std::dec << ")\n";
    std::cout << "  - Extended Model:  " << extended_model << " (0x" << std::hex << extended_model << std::dec << ")\n";
    std::cout << "  - Extended Family: " << extended_family << "\n";
    std::cout << "  - Stepping ID:     " << stepping << "\n\n";

    std::cout << "Calculated Values for Signature Matching:\n";
    std::cout << "  -> Display Family: " << display_family << "\n";
    std::cout << "  -> Display Model:  " << display_model << " (0x" << std::hex << display_model << std::dec << ")\n\n";

    // Explicit validation check against Skylake-X Server attributes
    std::cout << "Microarchitecture Identification Target:\n";
    if (display_family == 6 && display_model == 85) { // 85 == 0x55
        std::cout << "  [MATCH] Detected: Intel Skylake-X Server / Skylake-SP Microarchitecture.\n";
    } else {
        std::cout << "  [MISMATCH] Host processor is not a Skylake-X Server chip.\n";
    }
    std::cout << "====================================================\n";

    cpu_uarch_basic[0] = base_family;
    cpu_uarch_basic[1] = base_model;
    cpu_uarch_basic[2] = stepping;
    cpu_uarch_basic[3] = extended_model;
    cpu_uarch_basic[4] = extended_family;
    cpu_uarch_basic[5] = display_family;
    cpu_uarch_basic[6] = display_model;
}




__ATTR_HOT__
__ATTR_ALIGN__(32)
void unit_test_constructors();

void unit_test_constructors()
{
    using namespace gms::hw_profiling;
    [[nodiscard]] std::int32_t printf_retv;
    printf_retv = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED at day=%s,time=%s!!.\n",__func__,__DATE__,__TIME__);
    printf_retv = std::printf("===============================================================\n");
    printf_retv = std::printf("   SKX_cpu_operating_frequency_t class constructors test\n");
    printf_retv = std::printf("===============================================================\n");

    printf_retv = std::printf("[1] ============== Test Default Constructor =============\n");
    {
        SKX_cpu_operating_frequency_t default_ctor;
        printf_retv = std::printf("m_nruns=%llu\n",default_ctor.m_nruns);
        printf_retv = std::printf("m_nsamples  =%llu\n",default_ctor.m_nsamples);
        printf_retv = std::printf("m_totsamples=%llu\n",default_ctor.m_totsamples);
        printf_retv = std::printf("m_a.size=%llu\n",default_ctor.m_a.size());
        printf_retv = std::printf("m_b.size=%llu\n",default_ctor.m_b.size());
        printf_retv = std::printf("m_c.size=%llu\n",default_ctor.m_c.size());
        printf_retv = std::printf("m_d.size=%llu\n",default_ctor.m_d.size());
        printf_retv = std::printf("m_cpu_arch=%s\n",default_ctor.m_cpu_arch.c_str());
        printf_retv = std::printf("m_date=%s\n",default_ctor.m_date.c_str());
        printf_retv = std::printf("m_time=%s\n",default_ctor.m_time.c_str());
        printf_retv = std::printf("m_mean_a=%.17f\n",default_ctor.m_a_robust_mean);
        printf_retv = std::printf("m_mean_b=%.17f\n",default_ctor.m_b_robust_mean);
        printf_retv = std::printf("m_mean_c=%.17f\n",default_ctor.m_c_robust_mean);
        printf_retv = std::printf("m_mean_d=%.17f\n",default_ctor.m_d_robust_mean);
        printf_retv = std::printf("m_a_normality=%.17f\n",default_ctor.m_a_normality);
        printf_retv = std::printf("m_b_normality=%.17f\n",default_ctor.m_b_normality);
        printf_retv = std::printf("m_c_normality=%.17f\n",default_ctor.m_c_normality);
        printf_retv = std::printf("m_d_normality=%.17f\n",default_ctor.m_d_normality);
        printf_retv = std::printf("m_brief_descr=%s\n",default_ctor.metadata_t.m_brief_descr.c_str());
        printf_retv = std::printf("m_metric_expr=%s\n",default_ctor.metadata_t.m_metric_expr.c_str());
        printf_retv = std::printf("m_metric_group=%s\n",default_ctor.metadata_t.m_metric_group.c_str());
        printf_retv = std::printf("m_metric_name=%s\n",default_ctor.metadata_t.m_metric_name.c_str());
        printf_retv = std::printf("m_scale_unit=%s\n",default_ctor.metadata_t.m_scale_unit.c_str());
    }

     printf_retv = std::printf("[2] ============== Test Multi-argument Constructor =============\n");
    {
        constexpr std::size_t nruns = 50ull;
        constexpr std::size_t nsamples = 1000ull;
        std::uint32_t cpu_basic_info[7] = {0};
        show_intel_cpu_info(cpu_basic_info);
        const std::string cpu_uarch = std::to_string(cpu_basic_info[6]);
        const std::string date      = __DATE__;
        const std::string time      = __TIME__;
        SKX_cpu_operating_frequency_t multiarg_ctor(nruns,nsamples,cpu_uarch,date,time);
        printf_retv = std::printf("m_nruns=%llu\n",multiarg_ctor.m_nruns);
        printf_retv = std::printf("m_nsamples  =%llu\n",multiarg_ctor.m_nsamples);
        printf_retv = std::printf("m_totsamples=%llu\n",multiarg_ctor.m_totsamples);
        printf_retv = std::printf("m_a.size=%llu\n",multiarg_ctor.m_a.size());
        printf_retv = std::printf("m_b.size=%llu\n",multiarg_ctor.m_b.size());
        printf_retv = std::printf("m_c.size=%llu\n",multiarg_ctor.m_c.size());
        printf_retv = std::printf("m_d.size=%llu\n",multiarg_ctor.m_d.size());
        printf_retv = std::printf("m_cpu_arch=%s\n",multiarg_ctor.m_cpu_arch.c_str());
        printf_retv = std::printf("m_date=%s\n",multiarg_ctor.m_date.c_str());
        printf_retv = std::printf("m_time=%s\n",multiarg_ctor.m_time.c_str());
        printf_retv = std::printf("m_mean_a=%.17f\n",multiarg_ctor.m_a_robust_mean);
        printf_retv = std::printf("m_mean_b=%.17f\n",multiarg_ctor.m_b_robust_mean);
        printf_retv = std::printf("m_mean_c=%.17f\n",multiarg_ctor.m_c_robust_mean);
        printf_retv = std::printf("m_mean_d=%.17f\n",multiarg_ctor.m_d_robust_mean);
        printf_retv = std::printf("m_a_normality=%.17f\n",multiarg_ctor.m_a_normality);
        printf_retv = std::printf("m_b_normality=%.17f\n",multiarg_ctor.m_b_normality);
        printf_retv = std::printf("m_c_normality=%.17f\n",multiarg_ctor.m_c_normality);
        printf_retv = std::printf("m_d_normality=%.17f\n",multiarg_ctor.m_d_normality);
        printf_retv = std::printf("m_brief_descr=%s\n",multiarg_ctor.metadata_t.m_brief_descr.c_str());
        printf_retv = std::printf("m_metric_expr=%s\n",multiarg_ctor.metadata_t.m_metric_expr.c_str());
        printf_retv = std::printf("m_metric_group=%s\n",multiarg_ctor.metadata_t.m_metric_group.c_str());
        printf_retv = std::printf("m_metric_name=%s\n",multiarg_ctor.metadata_t.m_metric_name.c_str());
        printf_retv = std::printf("m_scale_unit=%s\n",multiarg_ctor.metadata_t.m_scale_unit.c_str());
    }
    printf_retv = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED at day=%s,time=%s!!.\n",__func__,__DATE__,__TIME__);
}


int main()
{
    unit_test_constructors();
    return 0;
}