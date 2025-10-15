#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_dyn_array.h"
#include "GMS_normality_test_v2.h"

/*
   icpc -o unit_test_normality_v2 -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_half.h GMS_dyn_array.h  GMS_normality_test_v2.h GMS_normality_test_v2.cpp unit_test_normality_test_v2.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_half.h GMS_dyn_array.h  GMS_normality_test_v2.h GMS_normality_test_v2.cpp unit_test_normality_test_v2.cpp

*/

__attribute__((hot))
__attribute__((noinline))
void unit_test_normality_v2_1();

void unit_test_normality_v2_1()
{
    //using namespace gms::math;
    using namespace gms;

    constexpr std::size_t n{200ull};
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

__attribute__((hot))
__attribute__((noinline))
void unit_test_normality_v2_iterated();

void unit_test_normality_v2_iterated()
{
     using namespace gms;
     
     std::size_t n{};
     std::clock_t seed_data{};
     std::clock_t seed_size{};
     darray_r8_t  x;
     double w{0.0};
     double a{1.0};
     double b{1.0};
     double a_inc{0.25};
     double b_inc{0.1};
     std::int32_t n_iter{0};
     do 
     {
          n_iter += 1;
          seed_data = std::clock();
          auto rand_data{std::bind(std::weibull_distribution<double>(a,b),std::mt19937(seed_data))};
          seed_size = std::clock();
          auto rand_size{std::bind(std::uniform_int_distribution<std::size_t>(100ull,4000ull),std::mt19937(seed_size))};
          n = rand_size.operator()();
          if(__builtin_expect(n<100ull,0)) n = 100ull;
          x = darray_r8_t(n);
          for(std::size_t __i{0ull}; __i != x.mnx; ++__i) 
          {
            x.m_data[__i] = rand_data.operator()();
          }
          w = 0.0;
          w = gms::math::calculate_W_statistics(x);
          a += a_inc;
          b += b_inc;
          std::printf("[Distribution: Weibull] -- Data size=%llu,Iteration=%d,W=%.7f\n",n,n_iter,w);
     }while(w<0.0);
}


int main()
{
    //unit_test_normality_v2_1();
    unit_test_normality_v2_iterated();
    return 0;
}