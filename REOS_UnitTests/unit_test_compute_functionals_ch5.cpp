#include <cstdio>
#include <cstdlib>
#include "GMS_compute_functionals_ch5.h"

/*
   icpc -o unit_test_compute_functionals_ch5 -fopenmp -O3 -fp-model fast=2 -fno-exceptions -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_cephes_double.h GMS_integrands_func_ch5.h GMS_integrands_func_ch5.cpp GMS_cquadpack.h GMS_cquadpack.c GMS_tabulated_quadrature.h GMS_tabulated_quadrature.cpp GMS_machine_utils.h GMS_compute_functionals_ch5.h GMS_compute_functionals_ch5.cpp unit_test_compute_functionals_ch5.cpp
   ASM: 
   icpc -S -O3 -fopenmp -fverbose-asm -masm=intel -fno-exceptions -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_cephes_double.h GMS_integrands_func_ch5.h GMS_integrands_func_ch5.cpp GMS_machine_utils.h GMS_compute_functionals_ch5.h GMS_compute_functionals_ch5.cpp unit_test_compute_functionals_ch5.cpp

*/

namespace 
{
static inline 
std::int32_t 
print_double(const char * val_name,double f, const std::int32_t offset) 
{
  (*((std::uint64_t *)&f)) += offset;
  std::int32_t printf_ret = std::printf("value: %s,%1.17g,%017llx,%llu\n",val_name,f, *(std::uint64_t *)&f, *(std::uint64_t *)&f);
  return (printf_ret);
}
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_compute_functional_Rayleigh_chan_5_6();

void unit_test_compute_functional_Rayleigh_chan_5_6()
{
    //using namespace gms::fading_channel;
    constexpr std::int32_t n_func_args{100};
    [[maybe_unused]] __ATTR_ALIGN__(64) gms::fading_channel::func_args_ch5_payload_t func_args[n_func_args]; // unused
    __ATTR_ALIGN__(16) double Rayleigh_Channel_functional[n_func_args];
    __ATTR_ALIGN__(16) double abser[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work1[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work2[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_results[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_start[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_end[n_func_args];
    __ATTR_ALIGN__(16) std::int32_t neval[n_func_args];
    __ATTR_ALIGN__(16) std::int32_t ier[n_func_args];
    __ATTR_ALIGN__(16) std::int32_t last[n_func_args];
    double epsabs[1];
    double epsrel[1];
    std::int32_t inf[1];
    std::int32_t irule[1];
    [[maybe_unused]] std::int32_t print_retv{};
    print_retv = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!.\n",__func__);
    epsabs[0] = 0.0;
    epsrel[0] = 0.001;
    inf[0]    = 1;
    irule[0]  = 4;
    gms::fading_channel::quadpack_integrator_payload_ch5_t integrator_payload;
    integrator_payload.integrand = &gms::fading_channel::integrand_Rayleigh_chan_5_6_iface;
    integrator_payload.func_args_payload = &func_args[0];
    integrator_payload.tmp_work1 = &tmp_work1[0];
    integrator_payload.tmp_work2 = &tmp_work2[0];
    integrator_payload.inf   = &inf[0];
    integrator_payload.irule = &irule[0];
    integrator_payload.epsabs = &epsabs[0];
    integrator_payload.epsrel = &epsrel[0];
    integrator_payload.abser  = &abser[0];
    integrator_payload.functional = &Rayleigh_Channel_functional[0];
    integrator_payload.neval = &neval[0];
    integrator_payload.ier = &ier[0];
    integrator_payload.last = &last[0];
    integrator_payload.crude_tsc_start = &crude_tsc_start[0];
    integrator_payload.crude_tsc_end   = &crude_tsc_end[0];
    integrator_payload.crude_tsc_measurement = &crude_tsc_results[0];
    integrator_payload.rand_lo1 = +0.1;
    integrator_payload.rand_hi1 = +1.0;
    integrator_payload.rand_lo2 = +5.0;
    integrator_payload.rand_hi2 = +15.0;
    integrator_payload.n_func_vals = n_func_args;
    integrator_payload.randomly_generate_inputs = true;
    //integrator_payload.which_integrator = 3;
    gms::fading_channel::quadpack_integrator_payload_ch5_t * __restrict__ p_payload = &integrator_payload;
    for(std::int32_t ii{1}; ii<5; ++ii) 
    {
        integrator_payload.which_integrator = ii;
        std::int32_t integrator_ret = gms::fading_channel::compute_functional_Rayleigh_chan_5_6(p_payload);
        if(p_payload->which_integrator==1)
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",gms::fading_channel::integrators_names_ch5[0].c_str());
        else if(p_payload->which_integrator==2)
           print_retv = std::printf("[UNIT-TEST:] -- 'DQAGI' Not-Supported!!\n");
        else if(p_payload->which_integrator==3)
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",gms::fading_channel::integrators_names_ch5[2].c_str());
        else if(p_payload->which_integrator==4)
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",gms::fading_channel::integrators_names_ch5[3].c_str());
        if(p_payload->which_integrator!=2)
        {
           for(std::int32_t i{0}; i < p_payload->n_func_vals; ++i) 
           {
               print_retv = print_double("Rayleigh-Channel (f:5.6)",p_payload->functional[i],0);
           }
           for(std::int32_t i{0}; i < p_payload->n_func_vals; ++i) 
           {
               const std::uint64_t start = crude_tsc_start[i];
               const std::uint64_t end   = crude_tsc_end[i];
               print_retv = std::printf("Iter=%d,abser=%.17f,neval=%d,ier=%d,last=%d,tsc_start=%llu,tsc_end=%llu,total_tsc=%llu,ret=%d\n",i,abser[i],neval[i],ier[i],last[i],start,end,crude_tsc_results[i],integrator_ret);
           }
        }
    }
    print_retv = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_compute_functional_Hoyt_chan_5_6();

void unit_test_compute_functional_Hoyt_chan_5_6()
{
    //using namespace gms::fading_channel;
    constexpr std::int32_t n_func_args{100};
    [[maybe_unused]] __ATTR_ALIGN__(64) gms::fading_channel::func_args_ch5_payload_t func_args[n_func_args]; // unused
    __ATTR_ALIGN__(16) double Rayleigh_Channel_functional[n_func_args];
    __ATTR_ALIGN__(16) double abser[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work1[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work2[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work3[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_results[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_start[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_end[n_func_args];
    __ATTR_ALIGN__(16) std::int32_t neval[n_func_args];
    __ATTR_ALIGN__(16) std::int32_t ier[n_func_args];
    __ATTR_ALIGN__(16) std::int32_t last[n_func_args];
    double epsabs[1];
    double epsrel[1];
    std::int32_t inf[1];
    std::int32_t irule[1];
    [[maybe_unused]] std::int32_t print_retv{};
    print_retv = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!.\n",__func__);
    epsabs[0] = 0.0;
    epsrel[0] = 0.001;
    inf[0]    = 1;
    irule[0]  = 4;
    gms::fading_channel::quadpack_integrator_payload_ch5_t integrator_payload;
    integrator_payload.integrand = &gms::fading_channel::integrand_Hoyt_chan_5_9_iface;
    integrator_payload.func_args_payload = &func_args[0];
    integrator_payload.tmp_work1 = &tmp_work1[0];
    integrator_payload.tmp_work2 = &tmp_work2[0];
    integrator_payload.tmp_work3 = &tmp_work3[0];
    integrator_payload.inf   = &inf[0];
    integrator_payload.irule = &irule[0];
    integrator_payload.epsabs = &epsabs[0];
    integrator_payload.epsrel = &epsrel[0];
    integrator_payload.abser  = &abser[0];
    integrator_payload.functional = &Rayleigh_Channel_functional[0];
    integrator_payload.neval = &neval[0];
    integrator_payload.ier = &ier[0];
    integrator_payload.last = &last[0];
    integrator_payload.crude_tsc_start = &crude_tsc_start[0];
    integrator_payload.crude_tsc_end   = &crude_tsc_end[0];
    integrator_payload.crude_tsc_measurement = &crude_tsc_results[0];
    integrator_payload.rand_lo1 = +0.1;
    integrator_payload.rand_hi1 = +1.0;
    integrator_payload.rand_lo2 = +5.0;
    integrator_payload.rand_hi2 = +15.0;
    integrator_payload.rand_lo3 = +0.1;
    integrator_payload.rand_hi3 = +1.0;
    integrator_payload.n_func_vals = n_func_args;
    integrator_payload.randomly_generate_inputs = true;
    //integrator_payload.which_integrator = 3;
    gms::fading_channel::quadpack_integrator_payload_ch5_t * __restrict__ p_payload = &integrator_payload;
    for(std::int32_t ii{1}; ii<5; ++ii) 
    {
        integrator_payload.which_integrator = ii;
        std::int32_t integrator_ret = gms::fading_channel::compute_functional_Hoyt_chan_5_9(p_payload);
        if(p_payload->which_integrator==1)
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",gms::fading_channel::integrators_names_ch5[0].c_str());
        else if(p_payload->which_integrator==2)
           print_retv = std::printf("[UNIT-TEST:] -- 'DQAGI' Not-Supported!!\n");
        else if(p_payload->which_integrator==3)
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",gms::fading_channel::integrators_names_ch5[2].c_str());
        else if(p_payload->which_integrator==4)
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",gms::fading_channel::integrators_names_ch5[3].c_str());
        if(p_payload->which_integrator!=2)
        {
           for(std::int32_t i{0}; i < p_payload->n_func_vals; ++i) 
           {
               print_retv = print_double("Hoyt-Channel (f:5.9)",p_payload->functional[i],0);
           }
           for(std::int32_t i{0}; i < p_payload->n_func_vals; ++i) 
           {
               const std::uint64_t start = crude_tsc_start[i];
               const std::uint64_t end   = crude_tsc_end[i];
               print_retv = std::printf("Iter=%d,abser=%.17f,neval=%d,ier=%d,last=%d,tsc_start=%llu,tsc_end=%llu,total_tsc=%llu,ret=%d\n",i,abser[i],neval[i],ier[i],last[i],start,end,crude_tsc_results[i],integrator_ret);
           }
        }
    }
    print_retv = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_compute_functional_Rice_chan_5_12();

void unit_test_compute_functional_Rice_chan_5_12()
{
    //using namespace gms::fading_channel;
    constexpr std::int32_t n_func_args{100};
    [[maybe_unused]] __ATTR_ALIGN__(64) gms::fading_channel::func_args_ch5_payload_t func_args[n_func_args]; // unused
    __ATTR_ALIGN__(16) double Rice_Channel_functional[n_func_args];
    __ATTR_ALIGN__(16) double abser[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work1[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work2[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work3[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_results[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_start[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_end[n_func_args];
    __ATTR_ALIGN__(16) std::int32_t neval[n_func_args];
    __ATTR_ALIGN__(16) std::int32_t ier[n_func_args];
    __ATTR_ALIGN__(16) std::int32_t last[n_func_args];
    double epsabs[1];
    double epsrel[1];
    std::int32_t inf[1];
    std::int32_t irule[1];
    [[maybe_unused]] std::int32_t print_retv{};
    print_retv = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!.\n",__func__);
    epsabs[0] = 0.0;
    epsrel[0] = 0.001;
    inf[0]    = 1;
    irule[0]  = 4;
    gms::fading_channel::quadpack_integrator_payload_ch5_t integrator_payload;
    integrator_payload.integrand = &gms::fading_channel::integrand_Rice_chan_5_12_iface;
    integrator_payload.func_args_payload = &func_args[0];
    integrator_payload.tmp_work1 = &tmp_work1[0];
    integrator_payload.tmp_work2 = &tmp_work2[0];
    integrator_payload.tmp_work3 = &tmp_work3[0];
    integrator_payload.inf   = &inf[0];
    integrator_payload.irule = &irule[0];
    integrator_payload.epsabs = &epsabs[0];
    integrator_payload.epsrel = &epsrel[0];
    integrator_payload.abser  = &abser[0];
    integrator_payload.functional = &Rice_Channel_functional[0];
    integrator_payload.neval = &neval[0];
    integrator_payload.ier = &ier[0];
    integrator_payload.last = &last[0];
    integrator_payload.crude_tsc_start = &crude_tsc_start[0];
    integrator_payload.crude_tsc_end   = &crude_tsc_end[0];
    integrator_payload.crude_tsc_measurement = &crude_tsc_results[0];
    integrator_payload.rand_lo1 = +0.1;
    integrator_payload.rand_hi1 = +1.0;
    integrator_payload.rand_lo2 = +5.0;
    integrator_payload.rand_hi2 = +15.0;
    integrator_payload.rand_lo3 = +0.1;
    integrator_payload.rand_hi3 = +1.0;
    integrator_payload.n_func_vals = n_func_args;
    integrator_payload.randomly_generate_inputs = true;
    //integrator_payload.which_integrator = 3;
    gms::fading_channel::quadpack_integrator_payload_ch5_t * __restrict__ p_payload = &integrator_payload;
    for(std::int32_t ii{1}; ii<5; ++ii) 
    {
        integrator_payload.which_integrator = ii;
        std::int32_t integrator_ret = gms::fading_channel::compute_functional_Rice_chan_5_12(p_payload);
        if(p_payload->which_integrator==1)
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",gms::fading_channel::integrators_names_ch5[0].c_str());
        else if(p_payload->which_integrator==2)
           print_retv = std::printf("[UNIT-TEST:] -- 'DQAGI' Not-Supported!!\n");
        else if(p_payload->which_integrator==3)
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",gms::fading_channel::integrators_names_ch5[2].c_str());
        else if(p_payload->which_integrator==4)
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",gms::fading_channel::integrators_names_ch5[3].c_str());
        if(p_payload->which_integrator!=2)
        {
           for(std::int32_t i{0}; i < p_payload->n_func_vals; ++i) 
           {
               print_retv = print_double("Rice-Channel (f:5.12)",p_payload->functional[i],0);
           }
           for(std::int32_t i{0}; i < p_payload->n_func_vals; ++i) 
           {
               const std::uint64_t start = crude_tsc_start[i];
               const std::uint64_t end   = crude_tsc_end[i];
               print_retv = std::printf("Iter=%d,abser=%.17f,neval=%d,ier=%d,last=%d,tsc_start=%llu,tsc_end=%llu,total_tsc=%llu,ret=%d\n",i,abser[i],neval[i],ier[i],last[i],start,end,crude_tsc_results[i],integrator_ret);
           }
        }
    }
    print_retv = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_compute_functional_Nakagami_m_chan_5_16();

void unit_test_compute_functional_Nakagami_m_chan_5_16()
{
    //using namespace gms::fading_channel;
    constexpr std::int32_t n_func_args{100};
    [[maybe_unused]] __ATTR_ALIGN__(64) gms::fading_channel::func_args_ch5_payload_t func_args[n_func_args]; // unused
    __ATTR_ALIGN__(16) double Nakagami_m_Channel_functional[n_func_args];
    __ATTR_ALIGN__(16) double abser[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work1[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work2[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work3[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_results[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_start[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_end[n_func_args];
    __ATTR_ALIGN__(16) std::int32_t neval[n_func_args];
    __ATTR_ALIGN__(16) std::int32_t ier[n_func_args];
    __ATTR_ALIGN__(16) std::int32_t last[n_func_args];
    double epsabs[1];
    double epsrel[1];
    std::int32_t inf[1];
    std::int32_t irule[1];
    [[maybe_unused]] std::int32_t print_retv{};
    print_retv = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!.\n",__func__);
    epsabs[0] = 0.0;
    epsrel[0] = 0.001;
    inf[0]    = 1;
    irule[0]  = 4;
    gms::fading_channel::quadpack_integrator_payload_ch5_t integrator_payload;
    integrator_payload.integrand = &gms::fading_channel::integrand_Nakagami_m_chan_5_16_iface;
    integrator_payload.func_args_payload = &func_args[0];
    integrator_payload.tmp_work1 = &tmp_work1[0];
    integrator_payload.tmp_work2 = &tmp_work2[0];
    integrator_payload.tmp_work3 = &tmp_work3[0];
    integrator_payload.inf   = &inf[0];
    integrator_payload.irule = &irule[0];
    integrator_payload.epsabs = &epsabs[0];
    integrator_payload.epsrel = &epsrel[0];
    integrator_payload.abser  = &abser[0];
    integrator_payload.functional = &Nakagami_m_Channel_functional[0];
    integrator_payload.neval = &neval[0];
    integrator_payload.ier = &ier[0];
    integrator_payload.last = &last[0];
    integrator_payload.crude_tsc_start = &crude_tsc_start[0];
    integrator_payload.crude_tsc_end   = &crude_tsc_end[0];
    integrator_payload.crude_tsc_measurement = &crude_tsc_results[0];
    integrator_payload.rand_lo1 = +0.1;
    integrator_payload.rand_hi1 = +1.0;
    integrator_payload.rand_lo2 = +5.0;
    integrator_payload.rand_hi2 = +15.0;
    integrator_payload.rand_lo3 = +0.5;
    integrator_payload.rand_hi3 = +3.0;
    integrator_payload.n_func_vals = n_func_args;
    integrator_payload.randomly_generate_inputs = true;
    //integrator_payload.which_integrator = 3;
    gms::fading_channel::quadpack_integrator_payload_ch5_t * __restrict__ p_payload = &integrator_payload;
    for(std::int32_t ii{1}; ii<5; ++ii) 
    {
        integrator_payload.which_integrator = ii;
        std::int32_t integrator_ret = gms::fading_channel::compute_functional_Nakagami_m_chan_5_16(p_payload);
        if(p_payload->which_integrator==1)
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",gms::fading_channel::integrators_names_ch5[0].c_str());
        else if(p_payload->which_integrator==2)
           print_retv = std::printf("[UNIT-TEST:] -- 'DQAGI' Not-Supported!!\n");
        else if(p_payload->which_integrator==3)
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",gms::fading_channel::integrators_names_ch5[2].c_str());
        else if(p_payload->which_integrator==4)
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",gms::fading_channel::integrators_names_ch5[3].c_str());
        if(p_payload->which_integrator!=2)
        {
           for(std::int32_t i{0}; i < p_payload->n_func_vals; ++i) 
           {
               print_retv = print_double("Rice-Channel (f:5.16)",p_payload->functional[i],0);
           }
           for(std::int32_t i{0}; i < p_payload->n_func_vals; ++i) 
           {
               const std::uint64_t start = crude_tsc_start[i];
               const std::uint64_t end   = crude_tsc_end[i];
               print_retv = std::printf("Iter=%d,abser=%.17f,neval=%d,ier=%d,last=%d,tsc_start=%llu,tsc_end=%llu,total_tsc=%llu,ret=%d\n",i,abser[i],neval[i],ier[i],last[i],start,end,crude_tsc_results[i],integrator_ret);
           }
        }
    }
    print_retv = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_compute_functional_LogNormShadow_chan_5_20();

void unit_test_compute_functional_LogNormShadow_chan_5_20()
{
    //using namespace gms::fading_channel;
    constexpr std::int32_t n_func_args{100};
    [[maybe_unused]] __ATTR_ALIGN__(64) gms::fading_channel::func_args_ch5_payload_t func_args[n_func_args]; // unused
    __ATTR_ALIGN__(16) double LogNormShadow_Channel_functional[n_func_args];
    __ATTR_ALIGN__(16) double abser[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work1[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work2[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work3[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work4[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_results[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_start[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_end[n_func_args];
    __ATTR_ALIGN__(16) std::int32_t neval[n_func_args];
    __ATTR_ALIGN__(16) std::int32_t ier[n_func_args];
    __ATTR_ALIGN__(16) std::int32_t last[n_func_args];
    double epsabs[1];
    double epsrel[1];
    std::int32_t inf[1];
    std::int32_t irule[1];
    [[maybe_unused]] std::int32_t print_retv{};
    print_retv = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!.\n",__func__);
    epsabs[0] = 0.0;
    epsrel[0] = 0.001;
    inf[0]    = 1;
    irule[0]  = 4;
    gms::fading_channel::quadpack_integrator_payload_ch5_t integrator_payload;
    integrator_payload.integrand = &gms::fading_channel::integrand_Log_Norm_Shadow_chan_5_20_iface;
    integrator_payload.func_args_payload = &func_args[0];
    integrator_payload.tmp_work1 = &tmp_work1[0];
    integrator_payload.tmp_work2 = &tmp_work2[0];
    integrator_payload.tmp_work3 = &tmp_work3[0];
    integrator_payload.tmp_work4 = &tmp_work4[0]; 
    integrator_payload.inf   = &inf[0];
    integrator_payload.irule = &irule[0];
    integrator_payload.epsabs = &epsabs[0];
    integrator_payload.epsrel = &epsrel[0];
    integrator_payload.abser  = &abser[0];
    integrator_payload.functional = &LogNormShadow_Channel_functional[0];
    integrator_payload.neval = &neval[0];
    integrator_payload.ier = &ier[0];
    integrator_payload.last = &last[0];
    integrator_payload.crude_tsc_start = &crude_tsc_start[0];
    integrator_payload.crude_tsc_end   = &crude_tsc_end[0];
    integrator_payload.crude_tsc_measurement = &crude_tsc_results[0];
    integrator_payload.rand_lo1 = +0.0;
    integrator_payload.rand_hi1 = +1.570796326794896619231321692;
    integrator_payload.rand_lo2 = +0.1;
    integrator_payload.rand_hi2 = +1.1;
    integrator_payload.rand_lo3 = +10.0*std::log10(5.0);
    integrator_payload.rand_hi3 = +10.0*std::log10(15.0);
    integrator_payload.rand_lo4 = +4.5;
    integrator_payload.rand_hi4 = +7.58;
    integrator_payload.rand_lo5 = +0.0;
    integrator_payload.rand_hi5 = +7.5;
    integrator_payload.n_func_vals = n_func_args;
    integrator_payload.randomly_generate_inputs = true;
    //integrator_payload.which_integrator = 3;
    gms::fading_channel::quadpack_integrator_payload_ch5_t * __restrict__ p_payload = &integrator_payload;
    for(std::int32_t ii{1}; ii<5; ++ii) 
    {
        integrator_payload.which_integrator = ii;
        std::int32_t integrator_ret = gms::fading_channel::compute_functional_LogNormShadow_chan_5_20(p_payload);
        if(p_payload->which_integrator==1)
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",gms::fading_channel::integrators_names_ch5[0].c_str());
        else if(p_payload->which_integrator==2)
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",gms::fading_channel::integrators_names_ch5[1].c_str());
        else if(p_payload->which_integrator==3)
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",gms::fading_channel::integrators_names_ch5[2].c_str());
        else if(p_payload->which_integrator==4)
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",gms::fading_channel::integrators_names_ch5[3].c_str());
         for(std::int32_t i{0}; i < p_payload->n_func_vals; ++i) 
         {
               print_retv = print_double("Log-Normal Shadowing-Channel (f:5.20)",p_payload->functional[i],0);
           }
           for(std::int32_t i{0}; i < p_payload->n_func_vals; ++i) 
           {
               const std::uint64_t start = crude_tsc_start[i];
               const std::uint64_t end   = crude_tsc_end[i];
               print_retv = std::printf("Iter=%d,abser=%.17f,neval=%d,ier=%d,last=%d,tsc_start=%llu,tsc_end=%llu,total_tsc=%llu,ret=%d\n",i,abser[i],neval[i],ier[i],last[i],start,end,crude_tsc_results[i],integrator_ret);
         }
    }
    print_retv = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_compute_outer_functional_LogNormShadow_chan_5_20();

void unit_test_compute_outer_functional_LogNormShadow_chan_5_20()
{
    //using namespace gms::fading_channel;
    constexpr std::int32_t n_func_args{100};
    constexpr std::int32_t n_outer_func_args{50};
    [[maybe_unused]] __ATTR_ALIGN__(64) gms::fading_channel::func_args_ch5_payload_t func_args[n_func_args]; // unused
    __ATTR_ALIGN__(16) double LogNormShadow_Channel_functional[n_func_args];
    __ATTR_ALIGN__(16) double outer_LogNormShadow_Channel_functional[n_outer_func_args];
    __ATTR_ALIGN__(16) double abser[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work1[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work2[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work3[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work4[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work5[2*n_func_args-1]; // hiordq workspace.
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_results[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_start[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_end[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_meas_outer[n_outer_func_args];
    __ATTR_ALIGN__(16) std::int32_t neval[n_func_args];
    __ATTR_ALIGN__(16) std::int32_t ier[n_func_args];
    __ATTR_ALIGN__(16) std::int32_t last[n_func_args];
    double epsabs[1];
    double epsrel[1];
    std::int32_t inf[1];
    std::int32_t irule[1];
    [[maybe_unused]] std::int32_t print_retv{};
    print_retv = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!.\n",__func__);
    epsabs[0] = 0.0;
    epsrel[0] = 0.001;
    inf[0]    = 1;
    irule[0]  = 4;
    gms::fading_channel::quadpack_integrator_payload_ch5_t integrator_payload;
    integrator_payload.integrand = &gms::fading_channel::integrand_Log_Norm_Shadow_chan_5_20_iface;
    integrator_payload.func_args_payload = &func_args[0];
    integrator_payload.tmp_work1 = &tmp_work1[0];
    integrator_payload.tmp_work2 = &tmp_work2[0];
    integrator_payload.tmp_work3 = &tmp_work3[0];
    integrator_payload.tmp_work4 = &tmp_work4[0]; 
    integrator_payload.tmp_work5 = &tmp_work5[0];
    integrator_payload.inf   = &inf[0];
    integrator_payload.irule = &irule[0];
    integrator_payload.epsabs = &epsabs[0];
    integrator_payload.epsrel = &epsrel[0];
    integrator_payload.abser  = &abser[0];
    integrator_payload.functional = &LogNormShadow_Channel_functional[0];
    integrator_payload.neval = &neval[0];
    integrator_payload.ier = &ier[0];
    integrator_payload.last = &last[0];
    integrator_payload.crude_tsc_start = &crude_tsc_start[0];
    integrator_payload.crude_tsc_end   = &crude_tsc_end[0];
    integrator_payload.crude_tsc_measurement = &crude_tsc_results[0];
    integrator_payload.crude_tsc_meas_outer  = &crude_tsc_meas_outer[0];
    integrator_payload.rand_lo1 = +0.0;
    integrator_payload.rand_hi1 = +1.570796326794896619231321692;
    integrator_payload.rand_lo2 = +0.1;
    integrator_payload.rand_hi2 = +1.1;
    integrator_payload.rand_lo3 = +10.0*std::log10(5.0);
    integrator_payload.rand_hi3 = +10.0*std::log10(15.0);
    integrator_payload.rand_lo4 = +4.5;
    integrator_payload.rand_hi4 = +7.58;
    integrator_payload.rand_lo5 = +0.0;
    integrator_payload.rand_hi5 = +7.5;
    integrator_payload.n_func_vals = n_func_args;
    integrator_payload.randomly_generate_inputs = true;
    integrator_payload.which_integrator = 1;
    gms::fading_channel::quadpack_integrator_payload_ch5_t * __restrict__ p_payload = &integrator_payload;
    integrator_payload.which_tabulated_integrator = 1;
    for(std::int32_t ii{0}; ii<n_outer_func_args; ++ii) 
    {
            std::int32_t integrator_ret = gms::fading_channel::compute_outer_functional_LogNormShadow_chan_5_20(p_payload);
            const double outer_func_res = p_payload->outer_func_tmp_res;
            outer_LogNormShadow_Channel_functional[ii] = outer_func_res;
            print_retv = print_double("Outer LogNormShadow-Channel Functional (f:5.20)",outer_func_res,0);
            print_retv = std::printf("[UNIT_TEST]: tsc-delta=%llu\n",p_payload->crude_tsc_meas_outer[ii]);
    }
    print_retv = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_compute_functional_compositeLogNormShadow_chan_5_25();

void unit_test_compute_functional_compositeLogNormShadow_chan_5_25()
{
    //using namespace gms::fading_channel;
    constexpr std::int32_t n_func_args{100};
    [[maybe_unused]] __ATTR_ALIGN__(64) gms::fading_channel::func_args_ch5_payload_t func_args[n_func_args]; // unused
    __ATTR_ALIGN__(16) double compositeLogNormShadow_Channel_functional[n_func_args];
    __ATTR_ALIGN__(16) double abser[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work1[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work2[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work3[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work4[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work5[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_results[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_start[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_end[n_func_args];
    __ATTR_ALIGN__(16) std::int32_t neval[n_func_args];
    __ATTR_ALIGN__(16) std::int32_t ier[n_func_args];
    __ATTR_ALIGN__(16) std::int32_t last[n_func_args];
    double epsabs[1];
    double epsrel[1];
    std::int32_t inf[1];
    std::int32_t irule[1];
    [[maybe_unused]] std::int32_t print_retv{};
    print_retv = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!.\n",__func__);
    epsabs[0] = 0.0;
    epsrel[0] = 0.001;
    inf[0]    = 1;
    irule[0]  = 4;
    gms::fading_channel::quadpack_integrator_payload_ch5_t integrator_payload;
    integrator_payload.integrand = &gms::fading_channel::integrand_LNSh_Nakagami_m_chan_5_25_iface;
    integrator_payload.func_args_payload = &func_args[0];
    integrator_payload.tmp_work1 = &tmp_work1[0];
    integrator_payload.tmp_work2 = &tmp_work2[0];
    integrator_payload.tmp_work3 = &tmp_work3[0];
    integrator_payload.tmp_work4 = &tmp_work4[0]; 
    integrator_payload.tmp_work5 = &tmp_work5[0];
    integrator_payload.inf   = &inf[0];
    integrator_payload.irule = &irule[0];
    integrator_payload.epsabs = &epsabs[0];
    integrator_payload.epsrel = &epsrel[0];
    integrator_payload.abser  = &abser[0];
    integrator_payload.functional = &compositeLogNormShadow_Channel_functional[0];
    integrator_payload.neval = &neval[0];
    integrator_payload.ier = &ier[0];
    integrator_payload.last = &last[0];
    integrator_payload.crude_tsc_start = &crude_tsc_start[0];
    integrator_payload.crude_tsc_end   = &crude_tsc_end[0];
    integrator_payload.crude_tsc_measurement = &crude_tsc_results[0];
    integrator_payload.rand_lo1 = +0.0;
    integrator_payload.rand_hi1 = +1.570796326794896619231321692;
    integrator_payload.rand_lo2 = +0.1;
    integrator_payload.rand_hi2 = +1.1;
    integrator_payload.rand_lo3 = +10.0*std::log10(5.0);
    integrator_payload.rand_hi3 = +10.0*std::log10(15.0);
    integrator_payload.rand_lo4 = +4.5;
    integrator_payload.rand_hi4 = +7.58;
    integrator_payload.rand_lo5 = +0.0;
    integrator_payload.rand_hi5 = +7.5;
    integrator_payload.rand_lo6 = +0.51;
    integrator_payload.rand_hi6 = +2.78956423555547;
    integrator_payload.n_func_vals = n_func_args;
    integrator_payload.randomly_generate_inputs = true;
    //integrator_payload.which_integrator = 3;
    gms::fading_channel::quadpack_integrator_payload_ch5_t * __restrict__ p_payload = &integrator_payload;
    for(std::int32_t ii{1}; ii<5; ++ii) 
    {
        integrator_payload.which_integrator = ii;
        std::int32_t integrator_ret = gms::fading_channel::compute_functional_LNSh_Nakagami_m_chan_5_25(p_payload);
        if(p_payload->which_integrator==1)
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",gms::fading_channel::integrators_names_ch5[0].c_str());
        else if(p_payload->which_integrator==2)
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",gms::fading_channel::integrators_names_ch5[1].c_str());
        else if(p_payload->which_integrator==3)
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",gms::fading_channel::integrators_names_ch5[2].c_str());
        else if(p_payload->which_integrator==4)
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",gms::fading_channel::integrators_names_ch5[3].c_str());
         for(std::int32_t i{0}; i < p_payload->n_func_vals; ++i) 
         {
               print_retv = print_double("Composite Log-Normal Shadowing-Channel (f:5.25)",p_payload->functional[i],0);
           }
           for(std::int32_t i{0}; i < p_payload->n_func_vals; ++i) 
           {
               const std::uint64_t start = crude_tsc_start[i];
               const std::uint64_t end   = crude_tsc_end[i];
               print_retv = std::printf("Iter=%d,abser=%.17f,neval=%d,ier=%d,last=%d,tsc_start=%llu,tsc_end=%llu,total_tsc=%llu,ret=%d\n",i,abser[i],neval[i],ier[i],last[i],start,end,crude_tsc_results[i],integrator_ret);
         }
    }
    print_retv = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_compute_outer_functional_compositeLogNormShadow_chan_5_25();

void unit_test_compute_outer_functional_compositeLogNormShadow_chan_5_25()
{
    //using namespace gms::fading_channel;
    constexpr std::int32_t n_func_args{100};
    constexpr std::int32_t n_outer_func_args{50};
    [[maybe_unused]] __ATTR_ALIGN__(64) gms::fading_channel::func_args_ch5_payload_t func_args[n_func_args]; // unused
    __ATTR_ALIGN__(16) double LogNormShadow_Channel_functional[n_func_args];
    __ATTR_ALIGN__(16) double outer_LogNormShadow_Channel_functional[n_outer_func_args];
    __ATTR_ALIGN__(16) double abser[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work1[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work2[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work3[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work4[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work5[2*n_func_args-1]; // hiordq workspace.
    __ATTR_ALIGN__(16) double tmp_work6[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_results[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_start[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_end[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_meas_outer[n_outer_func_args];
    __ATTR_ALIGN__(16) std::int32_t neval[n_func_args];
    __ATTR_ALIGN__(16) std::int32_t ier[n_func_args];
    __ATTR_ALIGN__(16) std::int32_t last[n_func_args];
    double epsabs[1];
    double epsrel[1];
    std::int32_t inf[1];
    std::int32_t irule[1];
    [[maybe_unused]] std::int32_t print_retv{};
    print_retv = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!.\n",__func__);
    epsabs[0] = 0.0;
    epsrel[0] = 0.001;
    inf[0]    = 1;
    irule[0]  = 4;
    gms::fading_channel::quadpack_integrator_payload_ch5_t integrator_payload;
    integrator_payload.integrand = &gms::fading_channel::integrand_LNSh_Nakagami_m_chan_5_25_iface;
    integrator_payload.func_args_payload = &func_args[0];
    integrator_payload.tmp_work1 = &tmp_work1[0];
    integrator_payload.tmp_work2 = &tmp_work2[0];
    integrator_payload.tmp_work3 = &tmp_work3[0];
    integrator_payload.tmp_work4 = &tmp_work4[0]; 
    integrator_payload.tmp_work5 = &tmp_work5[0];
    integrator_payload.tmp_work6 = &tmp_work6[0];
    integrator_payload.inf   = &inf[0];
    integrator_payload.irule = &irule[0];
    integrator_payload.epsabs = &epsabs[0];
    integrator_payload.epsrel = &epsrel[0];
    integrator_payload.abser  = &abser[0];
    integrator_payload.functional = &LogNormShadow_Channel_functional[0];
    integrator_payload.neval = &neval[0];
    integrator_payload.ier = &ier[0];
    integrator_payload.last = &last[0];
    integrator_payload.crude_tsc_start = &crude_tsc_start[0];
    integrator_payload.crude_tsc_end   = &crude_tsc_end[0];
    integrator_payload.crude_tsc_measurement = &crude_tsc_results[0];
    integrator_payload.crude_tsc_meas_outer  = &crude_tsc_meas_outer[0];
    integrator_payload.rand_lo1 = +0.0;
    integrator_payload.rand_hi1 = +1.570796326794896619231321692;
    integrator_payload.rand_lo2 = +0.1;
    integrator_payload.rand_hi2 = +1.1;
    integrator_payload.rand_lo3 = +10.0*std::log10(5.0);
    integrator_payload.rand_hi3 = +10.0*std::log10(15.0);
    integrator_payload.rand_lo4 = +4.5;
    integrator_payload.rand_hi4 = +7.58;
    integrator_payload.rand_lo5 = +0.0;
    integrator_payload.rand_hi5 = +7.5;
    integrator_payload.rand_lo6 = +0.51;
    integrator_payload.rand_hi6 = +2.78956423555547;
    integrator_payload.n_func_vals = n_func_args;
    integrator_payload.randomly_generate_inputs = true;
    integrator_payload.which_integrator = 1;
    gms::fading_channel::quadpack_integrator_payload_ch5_t * __restrict__ p_payload = &integrator_payload;
    integrator_payload.which_tabulated_integrator = 1;
    for(std::int32_t ii{0}; ii<n_outer_func_args; ++ii) 
    {
            std::int32_t integrator_ret = gms::fading_channel::compute_outer_functional_LNSh_Nakagami_m_chan_5_25(p_payload);
            const double outer_func_res = p_payload->outer_func_tmp_res;
            outer_LogNormShadow_Channel_functional[ii] = outer_func_res;
            print_retv = print_double("Outer Composite LogNormShadow-Channel Functional (f:5.25)",outer_func_res,0);
            print_retv = std::printf("[UNIT_TEST]: tsc-delta=%llu\n",p_payload->crude_tsc_meas_outer[ii]);
    }
    print_retv = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_compute_functional_Rayleigh_LaplaceT_chan_5_39();

void unit_test_compute_functional_Rayleigh_LaplaceT_chan_5_39()
{
    //using namespace gms::fading_channel;
    constexpr std::int32_t n_func_args{100};
    [[maybe_unused]] __ATTR_ALIGN__(64) gms::fading_channel::func_args_ch5_payload_t func_args[n_func_args]; // unused
    __ATTR_ALIGN__(16) double Rayleigh_LaplaceT_Channel_functional[n_func_args];
    __ATTR_ALIGN__(16) double abser[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work1[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work2[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work3[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work4[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_results[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_start[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_end[n_func_args];
    __ATTR_ALIGN__(16) std::int32_t neval[n_func_args];
    __ATTR_ALIGN__(16) std::int32_t ier[n_func_args];
    __ATTR_ALIGN__(16) std::int32_t last[n_func_args];
    double epsabs[1];
    double epsrel[1];
    std::int32_t inf[1];
    std::int32_t irule[1];
    [[maybe_unused]] std::int32_t print_retv{};
    print_retv = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!.\n",__func__);
    epsabs[0] = 0.0;
    epsrel[0] = 0.001;
    inf[0]    = 1;
    irule[0]  = 4;
    gms::fading_channel::quadpack_integrator_payload_ch5_t integrator_payload;
    integrator_payload.integrand = &gms::fading_channel::integrand_Rayleigh_lt_chan_5_39_iface;
    integrator_payload.func_args_payload = &func_args[0];
    integrator_payload.tmp_work1 = &tmp_work1[0];
    integrator_payload.tmp_work2 = &tmp_work2[0];
    integrator_payload.tmp_work3 = &tmp_work3[0];
    integrator_payload.tmp_work4 = &tmp_work4[0];
    integrator_payload.inf   = &inf[0];
    integrator_payload.irule = &irule[0];
    integrator_payload.epsabs = &epsabs[0];
    integrator_payload.epsrel = &epsrel[0];
    integrator_payload.abser  = &abser[0];
    integrator_payload.functional = &Rayleigh_LaplaceT_Channel_functional[0];
    integrator_payload.neval = &neval[0];
    integrator_payload.ier = &ier[0];
    integrator_payload.last = &last[0];
    integrator_payload.crude_tsc_start = &crude_tsc_start[0];
    integrator_payload.crude_tsc_end   = &crude_tsc_end[0];
    integrator_payload.crude_tsc_measurement = &crude_tsc_results[0];
    integrator_payload.rand_lo1 = +0.1;
    integrator_payload.rand_hi1 = +1.0;
    integrator_payload.rand_lo2 = +0.1;
    integrator_payload.rand_hi2 = +1.0;
    integrator_payload.rand_lo3 = +5.5;
    integrator_payload.rand_hi3 = +15.0;
    integrator_payload.rand_lo4 = +2.0;
    integrator_payload.rand_lo4 = 10.0;
    integrator_payload.n_func_vals = n_func_args;
    integrator_payload.randomly_generate_inputs = true;
    //integrator_payload.which_integrator = 3;
    gms::fading_channel::quadpack_integrator_payload_ch5_t * __restrict__ p_payload = &integrator_payload;
    std::int32_t integrators_cnt = 0;
    for(std::int32_t ii{1}; ii<5; ++ii) 
    {   
        ++integrators_cnt;
        integrator_payload.which_integrator = ii;
        std::int32_t integrator_ret = gms::fading_channel::compute_functional_Rayleigh_LaplaceT_chan_5_39(p_payload);
        if(p_payload->which_integrator==1)
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",gms::fading_channel::integrators_names_ch5[0].c_str());
        else if(p_payload->which_integrator==2)
           print_retv = std::printf("[UNIT-TEST:] -- 'DQAGI' Not-Supported!!\n");
        else if(p_payload->which_integrator==3)
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",gms::fading_channel::integrators_names_ch5[2].c_str());
        else if(p_payload->which_integrator==4)
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",gms::fading_channel::integrators_names_ch5[3].c_str());
        if(p_payload->which_integrator!=2)
        {
           const std::string append_integrator_name{gms::fading_channel::integrators_names_ch5[integrators_cnt].c_str()};
           gms::fading_channel::create_functional_ch5_plot(p_payload->n_func_vals,nullptr,&p_payload->functional[0],
                                                           "unit_test_compute_functional_Rayleigh_LaplaceT_chan_5_39_"+append_integrator_name,
                                                           "Functional of Rayleigh Fading Channel computed by:"+append_integrator_name,false);
                                                           
           for(std::int32_t i{0}; i < p_payload->n_func_vals; ++i) 
           {
               print_retv = print_double("Rayleigh-Fading Channel (f:5.39)",p_payload->functional[i],0);
           }
           for(std::int32_t i{0}; i < p_payload->n_func_vals; ++i) 
           {
               const std::uint64_t start = crude_tsc_start[i];
               const std::uint64_t end   = crude_tsc_end[i];
               print_retv = std::printf("Iter=%d,abser=%.17f,neval=%d,ier=%d,last=%d,tsc_start=%llu,tsc_end=%llu,total_tsc=%llu,ret=%d\n",i,abser[i],neval[i],ier[i],last[i],start,end,crude_tsc_results[i],integrator_ret);
           }
        }
    }
    print_retv = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_compute_functional_Hoyt_LaplaceT_chan_5_40();

void unit_test_compute_functional_Hoyt_LaplaceT_chan_5_40()
{
    //using namespace gms::fading_channel;
    constexpr std::int32_t n_func_args{100};
    [[maybe_unused]] __ATTR_ALIGN__(64) gms::fading_channel::func_args_ch5_payload_t func_args[n_func_args]; // unused
    __ATTR_ALIGN__(16) double Hoyt_LaplaceT_Channel_functional[n_func_args];
    __ATTR_ALIGN__(16) double abser[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work1[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work2[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work3[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work4[n_func_args];
    __ATTR_ALIGN__(16) double tmp_work5[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_results[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_start[n_func_args];
    __ATTR_ALIGN__(16) std::uint64_t crude_tsc_end[n_func_args];
    __ATTR_ALIGN__(16) std::int32_t neval[n_func_args];
    __ATTR_ALIGN__(16) std::int32_t ier[n_func_args];
    __ATTR_ALIGN__(16) std::int32_t last[n_func_args];
    double epsabs[1];
    double epsrel[1];
    std::int32_t inf[1];
    std::int32_t irule[1];
    [[maybe_unused]] std::int32_t print_retv{};
    print_retv = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!.\n",__func__);
    epsabs[0] = 0.0;
    epsrel[0] = 0.001;
    inf[0]    = 1;
    irule[0]  = 4;
    gms::fading_channel::quadpack_integrator_payload_ch5_t integrator_payload;
    integrator_payload.integrand = &gms::fading_channel::integrand_Hoyt_lt_chan_5_40_iface;
    integrator_payload.func_args_payload = &func_args[0];
    integrator_payload.tmp_work1 = &tmp_work1[0];
    integrator_payload.tmp_work2 = &tmp_work2[0];
    integrator_payload.tmp_work3 = &tmp_work3[0];
    integrator_payload.tmp_work4 = &tmp_work4[0];
    integrator_payload.tmp_work5 = &tmp_work5[0];
    integrator_payload.inf   = &inf[0];
    integrator_payload.irule = &irule[0];
    integrator_payload.epsabs = &epsabs[0];
    integrator_payload.epsrel = &epsrel[0];
    integrator_payload.abser  = &abser[0];
    integrator_payload.functional = &Hoyt_LaplaceT_Channel_functional[0];
    integrator_payload.neval = &neval[0];
    integrator_payload.ier = &ier[0];
    integrator_payload.last = &last[0];
    integrator_payload.crude_tsc_start = &crude_tsc_start[0];
    integrator_payload.crude_tsc_end   = &crude_tsc_end[0];
    integrator_payload.crude_tsc_measurement = &crude_tsc_results[0];
    integrator_payload.rand_lo1 = +0.1;
    integrator_payload.rand_hi1 = +1.0;
    integrator_payload.rand_lo2 = +0.1;
    integrator_payload.rand_hi2 = +1.0;
    integrator_payload.rand_lo3 = +5.5;
    integrator_payload.rand_hi3 = +15.0;
    integrator_payload.rand_lo4 = +2.0;
    integrator_payload.rand_lo4 = 10.0;
    integrator_payload.rand_lo5 = +0.1;
    integrator_payload.rand_hi5 = +1.0;
    integrator_payload.n_func_vals = n_func_args;
    integrator_payload.randomly_generate_inputs = true;
    //integrator_payload.which_integrator = 3;
    gms::fading_channel::quadpack_integrator_payload_ch5_t * __restrict__ p_payload = &integrator_payload;
    std::int32_t integrators_cnt = 0;
    for(std::int32_t ii{1}; ii<5; ++ii) 
    {   
        ++integrators_cnt;
        integrator_payload.which_integrator = ii;
        std::int32_t integrator_ret = gms::fading_channel::compute_functional_Hoyt_LaplaceT_chan_5_40(p_payload);
        if(p_payload->which_integrator==1)
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",gms::fading_channel::integrators_names_ch5[0].c_str());
        else if(p_payload->which_integrator==2)
           print_retv = std::printf("[UNIT-TEST:] -- 'DQAGI' Not-Supported!!\n");
        else if(p_payload->which_integrator==3)
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",gms::fading_channel::integrators_names_ch5[2].c_str());
        else if(p_payload->which_integrator==4)
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",gms::fading_channel::integrators_names_ch5[3].c_str());
        if(p_payload->which_integrator!=2)
        {
           const std::string append_integrator_name{gms::fading_channel::integrators_names_ch5[integrators_cnt].c_str()};
           gms::fading_channel::create_functional_ch5_plot(p_payload->n_func_vals,nullptr,&p_payload->functional[0],
                                                           "unit_test_compute_functional_Hoyt_LaplaceT_chan_5_40_"+append_integrator_name,
                                                           "Functional of Hoyt Fading Channel computed by:"+append_integrator_name,false);
                                                           
           for(std::int32_t i{0}; i < p_payload->n_func_vals; ++i) 
           {
               print_retv = print_double("Hoyt-Fading Channel (f:5.40)",p_payload->functional[i],0);
           }
           for(std::int32_t i{0}; i < p_payload->n_func_vals; ++i) 
           {
               const std::uint64_t start = crude_tsc_start[i];
               const std::uint64_t end   = crude_tsc_end[i];
               print_retv = std::printf("Iter=%d,abser=%.17f,neval=%d,ier=%d,last=%d,tsc_start=%llu,tsc_end=%llu,total_tsc=%llu,ret=%d\n",i,abser[i],neval[i],ier[i],last[i],start,end,crude_tsc_results[i],integrator_ret);
           }
        }
    }
    print_retv = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}



int main()
{
   //(void)unit_test_compute_functional_Rayleigh_chan_5_6();
   //(void)unit_test_compute_functional_Hoyt_chan_5_6();
   //(void)unit_test_compute_functional_Rice_chan_5_12();
   //(void)unit_test_compute_functional_Nakagami_m_chan_5_16();
   //(void)unit_test_compute_functional_LogNormShadow_chan_5_20();
   //(void)unit_test_compute_outer_functional_LogNormShadow_chan_5_20();
   //(void)unit_test_compute_functional_compositeLogNormShadow_chan_5_25();
   //(void)unit_test_compute_outer_functional_compositeLogNormShadow_chan_5_25();
   //(void)unit_test_compute_functional_Rayleigh_LaplaceT_chan_5_39();
   (void)unit_test_compute_functional_Hoyt_LaplaceT_chan_5_40();
    return 0;
}