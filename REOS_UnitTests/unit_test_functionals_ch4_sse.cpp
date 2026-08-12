
#include <cstdio>
#include <cstdlib>
#include "GMS_functionals_ch4_sse.h"

/*
   icpc -o unit_test_functionals_ch4_sse -O3 -fp-model fast=2 -fno-exceptions -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_cephes_double.h GMS_integrands_ch4.h GMS_integrands_ch4.cpp GMS_bessel_i0_sse.h GMS_integrands_func_ch4_sse.h GMS_simpne_quad.h GMS_simpne_quad.cpp GMS_functionals_ch4_sse.h GMS_functionals_ch4_sse.cpp unit_test_functionals_ch4_sse.cpp
   ASM: 
   icpc -S -O3 -fverbose-asm -masm=intel -fno-exceptions -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_cephes_double.h GMS_integrands_ch4.h GMS_integrands_ch4.cpp GMS_bessel_i0_sse.h GMS_integrands_func_ch4_sse.h GMS_simpne_quad.h GMS_simpne_quad.cpp unit_test_functionals_ch4_sse.cpp

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
void unit_test_functional_4_1_gauss_Q_sse();

void unit_test_functional_4_1_gauss_Q_sse()
{
    //using namespace gms::fading_channel;
    constexpr std::int32_t n_func_args{10};
    constexpr std::int32_t n_gauss_q_vals{50};
    constexpr std::int32_t tot_elems{n_func_args*n_gauss_q_vals};
    constexpr double lo1_x{0.1};
    constexpr double hi1_x{1.0};
    constexpr double lo2_x{lo1_x};
    constexpr double hi2_x{hi1_x};
    __ATTR_ALIGN__(16) double in_buf[tot_elems];
    __ATTR_ALIGN__(16) double out_buf[tot_elems];
    __ATTR_ALIGN__(16) double gauss_Q_functional[n_func_args];
    [[maybe_unused]] std::int32_t print_retv{};
    print_retv = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!.\n",__func__);
    gms::fading_channel::integrator_payload_sse_t gauss_Q_payload;
    gauss_Q_payload.in_buf  = &in_buf[0];
    gauss_Q_payload.out_buf = &out_buf[0];
    gauss_Q_payload.functional = &gauss_Q_functional[0];
    gauss_Q_payload.lo1 = lo1_x;
    gauss_Q_payload.hi1 = hi1_x;
    gauss_Q_payload.lo2 = lo2_x;
    gauss_Q_payload.hi2 = hi2_x;
    gauss_Q_payload.n_func_args = n_func_args;
    gauss_Q_payload.n_integrand_vals = n_gauss_q_vals;
    gms::fading_channel::integrator_payload_sse_t * __restrict__ p_payload = &gauss_Q_payload;
    std::int32_t ret_val = gms::fading_channel::integrate_4_1_gauss_Q_func_sse(p_payload);
#if 1
    for(std::int32_t i{0}; i < p_payload->n_func_args; ++i) 
    {
        [[maybe_unused]] std::int32_t printf_ret = print_double("Gauss-Q (f:4.1)",p_payload->functional[i],0);
    }
#endif 
    print_retv = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}


int main()
{
    unit_test_functional_4_1_gauss_Q_sse();
    return 0;
}