
#include <cstdio>
#include <cstdlib>
#include "GMS_compute_functionals_ch4.h"

/*
   icpc -o unit_test_compute_functionals_ch4 -fopenmp -O3 -fp-model fast=2 -fno-exceptions -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_cephes_double.h GMS_integrands_func_ch4.h GMS_integrands_func_ch4.cpp GMS_cquadpack.h GMS_cquadpack.c GMS_compute_functionals_ch4.h GMS_compute_functionals_ch4.cpp unit_test_compute_functionals_ch4.cpp
   ASM: 
   icpc -S -O3 -fopenmp -fverbose-asm -masm=intel -fno-exceptions -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_cephes_double.h GMS_integrands_func_ch4.h GMS_integrands_func_ch4.cpp GMS_compute_functionals_ch4.h GMS_compute_functionals_ch4.cpp unit_test_compute_functionals_ch4.cpp

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


/*

    double (*integrand)(double,void * __restrict__);
    func_args_payload_t * __restrict__ func_args_payload{nullptr};
    double              * __restrict__ lo{nullptr}; // lower limit of integration
    double              * __restrict__ hi{nullptr}; // upper limit of integration
    double              * __restrict__ bound{nullptr}; //optional finite bound on integral.
    std::int32_t        * __restrict__ inf{nullptr}; // range of integration
    std::int32_t        * __restrict__ irule{nullptr};
    double              * __restrict__ epsabs{nullptr};
    double              * __restrict__ epsrel{nullptr};
    double              * __restrict__ abser{nullptr}; // integrator abosulute error
    double              * __restrict__ functional{nullptr}; // computed functional values
    std::int32_t        * __restrict__ neval{nullptr};      // number of evaluation
    std::int32_t        * __restrict__ ier{nullptr};        // integrator error indicator
    std::int32_t        * __restrict__ last{nullptr};
    std::int32_t                       n_func_vals;
    std::int32_t                       which_integrator; //currently 1=dqage,2=dqagi
 
*/

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_compute_functional_4_1_gauss_Q();

void unit_test_compute_functional_4_1_gauss_Q()
{
    //using namespace gms::fading_channel;
    constexpr std::int32_t n_func_args{100};
    [[maybe_unused]] __ATTR_ALIGN__(64) gms::fading_channel::func_args_payload_t func_args[1]; // unused
    __ATTR_ALIGN__(16) double gauss_Q_functional[n_func_args];
    __ATTR_ALIGN__(16) double lo_lim[n_func_args];
    __ATTR_ALIGN__(16) double hi_lim[n_func_args];
    __ATTR_ALIGN__(16) double bound[n_func_args];
    __ATTR_ALIGN__(16) double abser[n_func_args];
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
    gms::fading_channel::quadpack_integrator_payload_t integrator_payload;
    integrator_payload.integrand = &gms::fading_channel::integrand_4_1_gauss_Q_quad_iface;
    integrator_payload.func_args_payload = nullptr;
    integrator_payload.lo = &lo_lim[0];
    integrator_payload.hi = &hi_lim[0];
    integrator_payload.bound = &bound[0];
    integrator_payload.inf   = &inf[0];
    integrator_payload.irule = &irule[0];
    integrator_payload.epsabs = &epsabs[0];
    integrator_payload.epsrel = &epsrel[0];
    integrator_payload.abser  = &abser[0];
    integrator_payload.functional = &gauss_Q_functional[0];
    integrator_payload.neval = &neval[0];
    integrator_payload.ier = &ier[0];
    integrator_payload.last = &last[0];
    integrator_payload.rand_lo1 = -4.0;
    integrator_payload.rand_hi1 = +4.0;
    integrator_payload.rand_lo2 = +5.0;
    integrator_payload.rand_hi2 = +100.0;
    integrator_payload.rand_lo3 = -4.0;
    integrator_payload.rand_hi3 = +10.0;
    integrator_payload.n_func_vals = n_func_args;
    integrator_payload.which_integrator = 3;
    gms::fading_channel::quadpack_integrator_payload_t * __restrict__ p_payload = &integrator_payload;
    std::int32_t integrator_ret = gms::fading_channel::compute_functional_gauss_Q_4_1(p_payload);
    if(p_payload->which_integrator==1)
       print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",p_payload->integrators_names[0].c_str());
    else if(p_payload->which_integrator==2)
       print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",p_payload->integrators_names[1].c_str());
    else if(p_payload->which_integrator==3)
       print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",p_payload->integrators_names[2].c_str());
    for(std::int32_t i{0}; i < p_payload->n_func_vals; ++i) 
    {
        print_retv = print_double("Gauss-Q (f:4.1)",p_payload->functional[i],0);
    }
    for(std::int32_t i{0}; i < p_payload->n_func_vals; ++i) 
    {
        print_retv = std::printf("Iter=%d,abser=%.17f,neval=%d,ier=%d,last=%d\n",i,abser[i],neval[i],ier[i],last[i]);
    }
    print_retv = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_compute_functional_4_2_gauss_Q();

void unit_test_compute_functional_4_2_gauss_Q()
{
    //using namespace gms::fading_channel;
    constexpr std::int32_t n_func_args{100};
    [[maybe_unused]] __ATTR_ALIGN__(64) gms::fading_channel::func_args_payload_t func_args[n_func_args]; // unused
    __ATTR_ALIGN__(16) double gauss_Q_functional[n_func_args];
    __ATTR_ALIGN__(16) double abser[n_func_args];
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
    gms::fading_channel::quadpack_integrator_payload_t integrator_payload;
    integrator_payload.integrand = &gms::fading_channel::integrand_4_2_gauss_Q_quad_iface;
    integrator_payload.func_args_payload = &func_args[0];
    integrator_payload.inf   = &inf[0];
     integrator_payload.irule = &irule[0];
    integrator_payload.epsabs = &epsabs[0];
    integrator_payload.epsrel = &epsrel[0];
    integrator_payload.abser  = &abser[0];
    integrator_payload.functional = &gauss_Q_functional[0];
    integrator_payload.neval = &neval[0];
    integrator_payload.ier = &ier[0];
    integrator_payload.last = &last[0];
    integrator_payload.rand_lo1 = +0.1;
    integrator_payload.rand_hi1 = +1.570796326794896619231321692;
    integrator_payload.rand_lo2 = +0.1;
    integrator_payload.rand_hi2 = +10.0;
    integrator_payload.n_func_vals = n_func_args;
    integrator_payload.which_integrator = 1;
    gms::fading_channel::quadpack_integrator_payload_t * __restrict__ p_payload = &integrator_payload;
    std::int32_t integrator_ret = gms::fading_channel::compute_functional_gauss_Q_4_2(p_payload);
    if(p_payload->which_integrator==1)
       print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",p_payload->integrators_names[0].c_str());
    else if(p_payload->which_integrator==2)
       print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",p_payload->integrators_names[1].c_str());
    else if(p_payload->which_integrator==3)
       print_retv = std::printf("[UNIT-TEST:] -- currently executing: %s integrator\n",p_payload->integrators_names[2].c_str());
    for(std::int32_t i{0}; i < p_payload->n_func_vals; ++i) 
    {
        print_retv = print_double("Gauss-Q (f:4.1)",p_payload->functional[i],0);
    }
    for(std::int32_t i{0}; i < p_payload->n_func_vals; ++i) 
    {
        print_retv = std::printf("Iter=%d,abser=%.17f,neval=%d,ier=%d,last=%d\n",i,abser[i],neval[i],ier[i],last[i]);
    }
    print_retv = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

int main()
{
    (void)unit_test_compute_functional_4_1_gauss_Q();
    return 0;
}