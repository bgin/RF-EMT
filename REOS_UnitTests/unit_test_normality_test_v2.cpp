#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_dyn_array.h"
#include "GMS_normality_test_v2.h"

/*
   icpc -o unit_test_normality_v2 -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h  GMS_normality_test_v2.h GMS_normality_test_v2.cpp unit_test_normality_test_v2.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h  GMS_normality_test_v2.h GMS_normality_test_v2.cpp unit_test_normality_test_v2.cpp

*/

__attribute__((hot))
__attribute__((noinline))
void unit_test_normality_v2_1();

void unit_test_normality_v2_1()
{
    //using namespace gms::math;
    using namespace gms;

    constexpr std::size_t n{4000ull};
    double w{};
    darray_r8_t x;
    std::clock_t seed;
    x = darray_r8_t(n);
    
    seed = std::clock();
    auto rand1{std::bind(std::normal_distribution<double>(0.0,1.0),std::mt19937(seed))};
    for(std::size_t __i{0}; __i != n; ++__i) 
    {
         x.m_data[__i] = rand1();
    }
    
    w = gms::math::calculate_W_statistics(x);
    std::printf("[Distribution: Normal]: -- shapiro_wilk normality test (v2) executed correctly.\n");
    std::printf("Shapiro-Wilk (W)=%.7f\n",w);
    seed = std::clock();
    auto rand2{std::bind(std::weibull_distribution<double>(1.0f,1.0f),std::mt19937(seed))};
    for(std::size_t __i{0}; __i != n; ++__i) 
    {
         x.m_data[__i] = rand2();
    }
    w = 0.0;
    w = gms::math::calculate_W_statistics(x);    
    std::printf("[Distribution: Weibull]: shapiro_wilk normality test (v2) executed correctly.\n");
    std::printf("Shapiro-Wilk (W)=%.7f\n",w); 

}


int main()
{
    unit_test_normality_v2_1();
    return 0;
}