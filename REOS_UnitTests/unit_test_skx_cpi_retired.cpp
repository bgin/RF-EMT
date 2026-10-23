
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <string>
#include "GMS_skx_cpi_retired_metric.h"

/*
   icpc -o unit_test_skx_cpi_retired -fp-model fast=2 -ftz -std=c++17 -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_skx_clk_hw_events_metrics.hpp GMS_skx_cpi_retired_metric.h GMS_skx_cpi_retired_metric.cpp unit_test_skx_cpi_retired.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 -fopenmp GMS_config.h GMS_skx_clk_hw_events_metrics.hpp GMS_skx_cpi_retired_metric.h GMS_skx_cpi_retired_metric.cpp unit_test_skx_cpi_retired.cpp

*/


__ATTR_HOT__
__ATTR_ALIGN__(32)
void unit_test_constructors();

void unit_test_constructors()
{
    using namespace gms::hw_profiling;
    [[nodiscard]] std::int32_t printf_retv;
    printf_retv = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED at day=%s,time=%s!!.\n",__func__,__DATE__,__TIME__);
    printf_retv = std::printf("===============================================================\n");
    printf_retv = std::printf("    SKX_cpu_cpi_retired_t class constructors test\n");
    printf_retv = std::printf("===============================================================\n");

    printf_retv = std::printf("[1] ============== Test Default Constructor =============\n");
    {
        SKX_cpu_cpi_retired_t default_ctor;
        printf_retv = std::printf("m_nruns=%llu\n",default_ctor.m_nruns);
        printf_retv = std::printf("m_nsamples  =%llu\n",default_ctor.m_nsamples);
        printf_retv = std::printf("m_totsamples=%llu\n",default_ctor.m_totsamples);
        printf_retv = std::printf("m_a.size=%llu\n",default_ctor.m_a.size());
        printf_retv = std::printf("m_b.size=%llu\n",default_ctor.m_b.size());
        printf_retv = std::printf("m_c.size=%llu\n",default_ctor.m_c.size());
        printf_retv = std::printf("m_date=%s\n",default_ctor.m_date.c_str());
        printf_retv = std::printf("m_time=%s\n",default_ctor.m_time.c_str());
        printf_retv = std::printf("m_mean_a=%.17f\n",default_ctor.m_a_robust_mean);
        printf_retv = std::printf("m_mean_b=%.17f\n",default_ctor.m_b_robust_mean);
        printf_retv = std::printf("m_mean_c=%.17f\n",default_ctor.m_c_robust_mean);
        printf_retv = std::printf("m_a_normality=%.17f\n",default_ctor.m_a_normality);
        printf_retv = std::printf("m_b_normality=%.17f\n",default_ctor.m_b_normality);
        printf_retv = std::printf("m_c_normality=%.17f\n",default_ctor.m_c_normality);
        printf_retv = std::printf("m_brief_descr=%s\n",default_ctor.metadata_t.m_brief_descr.c_str());
        printf_retv = std::printf("m_metric_expr=%s\n",default_ctor.metadata_t.m_metric_expr.c_str());
        printf_retv = std::printf("m_metric_group=%s\n",default_ctor.metadata_t.m_metric_group.c_str());
        printf_retv = std::printf("m_metric_name=%s\n",default_ctor.metadata_t.m_metric_name.c_str());
        printf_retv = std::printf("m_scale_unit=%s\n",default_ctor.metadata_t.m_scale_unit.c_str());
        for(std::int32_t i = 0; i<8; ++i) { printf_retv = std::printf("m_cpuid_basic[%d]=%d\n",i,default_ctor.m_cpuid_basic[i]);}
    }

     printf_retv = std::printf("[2] ============== Test Multi-argument Constructor =============\n");
    {
        constexpr std::size_t nruns = 50ull;
        constexpr std::size_t nsamples = 1000ull;
        const std::string date      = __DATE__;
        const std::string time      = __TIME__;
        SKX_cpu_cpi_retired_t multiarg_ctor(nruns,nsamples,date,time,nullptr);
        printf_retv = std::printf("m_nruns=%llu\n",multiarg_ctor.m_nruns);
        printf_retv = std::printf("m_nsamples  =%llu\n",multiarg_ctor.m_nsamples);
        printf_retv = std::printf("m_totsamples=%llu\n",multiarg_ctor.m_totsamples);
        printf_retv = std::printf("m_a.size=%llu\n",multiarg_ctor.m_a.size());
        printf_retv = std::printf("m_b.size=%llu\n",multiarg_ctor.m_b.size());
        printf_retv = std::printf("m_c.size=%llu\n",multiarg_ctor.m_c.size());
        printf_retv = std::printf("m_date=%s\n",multiarg_ctor.m_date.c_str());
        printf_retv = std::printf("m_time=%s\n",multiarg_ctor.m_time.c_str());
        printf_retv = std::printf("m_mean_a=%.17f\n",multiarg_ctor.m_a_robust_mean);
        printf_retv = std::printf("m_mean_b=%.17f\n",multiarg_ctor.m_b_robust_mean);
        printf_retv = std::printf("m_mean_c=%.17f\n",multiarg_ctor.m_c_robust_mean);
        printf_retv = std::printf("m_a_normality=%.17f\n",multiarg_ctor.m_a_normality);
        printf_retv = std::printf("m_b_normality=%.17f\n",multiarg_ctor.m_b_normality);
        printf_retv = std::printf("m_c_normality=%.17f\n",multiarg_ctor.m_c_normality);
        printf_retv = std::printf("m_brief_descr=%s\n",multiarg_ctor.metadata_t.m_brief_descr.c_str());
        printf_retv = std::printf("m_metric_expr=%s\n",multiarg_ctor.metadata_t.m_metric_expr.c_str());
        printf_retv = std::printf("m_metric_group=%s\n",multiarg_ctor.metadata_t.m_metric_group.c_str());
        printf_retv = std::printf("m_metric_name=%s\n",multiarg_ctor.metadata_t.m_metric_name.c_str());
        printf_retv = std::printf("m_scale_unit=%s\n",multiarg_ctor.metadata_t.m_scale_unit.c_str());
        for(std::int32_t i = 0; i<8; ++i) { printf_retv = std::printf("m_cpuid_basic[%d]=%d\n",i,multiarg_ctor.m_cpuid_basic[i]);}
    }
    printf_retv = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED at day=%s,time=%s!!.\n",__func__,__DATE__,__TIME__);
}


int main()
{
    unit_test_constructors();
    return 0;
}