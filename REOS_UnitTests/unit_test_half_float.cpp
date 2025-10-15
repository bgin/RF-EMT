#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <algorithm>
#include <valarray>
#include <iostream>
#include "GMS_half.h"

/*
   icpc -o unit_test_half_float -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_half.h unit_test_half_float.cpp 
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_half.h unit_test_half_float.cpp 

*/

/*
__attribute__((hot))
__attribute__((noinline))
void unit_test_half_float();

void unit_test_half_float()
{
    using half_float::half;
    using namespace half_float::literal
    constexpr std::uint32_t n_runs{10};
    std::uniform_real_distribution<double> distros[10]{};
    std::uint64_t data_seeds[10];
    std::clock_t buf_len_seed;
    half lo{0.0_h};
    half hi{10.0_h};
    buf_len_seed = std::clock();
    auto rand{std::bind(std::uniform_int_distribution<std::size_t>(1000ull,10000ull),std::mt19937(buf_len_seed))};
    for(std::uint32_t __i{0u}; __i != n_runs; ++__i) 
    {
        //const half h__i{static_cast<half>(__i)};
        //lo *= h__i;
        //hi *= h__i;
        
        std::size_t vec_len{rand()};
        if(__builtin_expect(vec_len==0ull,0)) vec_len = 1000ull;
        std::valarray<half> vec(vec_len);
        std::valarray<half> res(vec_len);
        data_seeds[__i] = __rdtsc();
        distros[__i]    = std::uniform_real_distribution<double>(0.0,1.0);
        auto rdev{std::mt19937(data_seeds[__i])};
        for(std::uint64_t __j{0ull}; __j != vec_len; ++__j) 
        {
               
                res.operator[](__j).operator=(std::sin<half>(vec.operator[](__j).operator=(half_cast<half>(distros.operator()(rdev())))));
        }
    }
}
*/

void unit_test_half_float_simple();

void unit_test_half_float_simple()
{
      using half_float::half;
      using namespace half_float::literal;
      {
        half a(3.4), b(5);
        half c = a * b;
        c += 3.0_h;
        std::cout << sizeof(c) << std::endl;
      }
      {
        half a(-3.14159);
        half s = sin(abs(a));
        long l = lround(s);
      }
      
}



int main()
{
    //unit_test_half_float();
    unit_test_half_float_simple();
    return 0;
}
