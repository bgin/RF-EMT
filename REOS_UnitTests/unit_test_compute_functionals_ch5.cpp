#include <cstdio>
#include <cstdlib>
#include "GMS_compute_functionals_ch5.h"

/*
   icpc -o unit_test_compute_functionals_ch5 -fopenmp -O3 -fp-model fast=2 -fno-exceptions -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_cephes_double.h GMS_integrands_func_ch5.h GMS_integrands_func_ch5.cpp GMS_cquadpack.h GMS_cquadpack.c GMS_machine_utils.h GMS_compute_functionals_ch5.h GMS_compute_functionals_ch5.cpp unit_test_compute_functionals_ch5.cpp
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
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",p_payload->integrators_names[0].c_str());
        else if(p_payload->which_integrator==2)
           print_retv = std::printf("[UNIT-TEST:] -- 'DQAGI' Not-Supported!!\n");
        else if(p_payload->which_integrator==3)
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",p_payload->integrators_names[2].c_str());
        else if(p_payload->which_integrator==4)
           print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",p_payload->integrators_names[3].c_str());
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


int main()
{
    (void)unit_test_compute_functional_Rayleigh_chan_5_6();
    return 0;
}