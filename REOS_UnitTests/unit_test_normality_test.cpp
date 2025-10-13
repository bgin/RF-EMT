#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <algorithm>
#include "GMS_dyn_array.h"
#include "GMS_normality_test.h"

/*
   icpc -o unit_test_normality -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h  GMS_normality_test.h GMS_normality_test.cpp unit_test_normality_test.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h  GMS_normality_test.h GMS_normality_test.cpp unit_test_normality_test.cpp

*/

__attribute__((hot))
__attribute__((noinline))
void unit_test_normality();

void unit_test_normality()
{
    //using namespace gms::math;
    using namespace gms;

    constexpr std::size_t n{4000ull};
    constexpr std::size_t n1{n};
    constexpr std::size_t n2{n/2ull};
    float w;
    float pw;
    int32_t ier;
    bool init;
    darray_r4_t x;
    darray_r4_t a;
    std::clock_t seed;
    x = darray_r4_t(n);
    a = darray_r4_t(n2);
    seed = std::clock();
    auto rand{std::bind(std::normal_distribution<float>(5.0f,10.0f),std::mt19937(seed))};
    for(std::size_t __i{0}; __i != n; ++__i) 
    {
         x.m_data[__i] = rand();
    }
    std::sort(x.m_data,x.m_data+n, [](float x,float y){ return x < y;});
    init  = false;
    gms::math::shapiro_wilk(init,x.m_data,static_cast<std::int32_t>(n),static_cast<std::int32_t>(n1),static_cast<std::int32_t>(n2),a.m_data,w,pw,ier);
    if(ier==0)
    {
        std::printf("ier=%d -- shapiro_wilk normality test executed correctly.\n",ier);
        std::printf("Shapiro-Wilk (W)=%.7f, P-value=%.7f\n",w,pw);
    }
    else 
    {
        std::printf("ier=%d -- shapiro_wilk normality test erred.\n",ier);
    }

}

int main()
{
    unit_test_normality();
    return 0;
}