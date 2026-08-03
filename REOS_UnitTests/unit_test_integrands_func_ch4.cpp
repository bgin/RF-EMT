#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <cmath>
#include "GMS_integrands_func_ch4.h"

/*
   icpc -o unit_test_integrands_func_ch4 -O3 -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_cephes_double.h GMS_integrands_func_ch4.h GMS_integrands_func_ch4.cpp unit_test_integrands_func_ch4.cpp
   ASM: 
   icpc -S -O3 -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_cephes_double.h GMS_integrands_func_ch4.h GMS_integrands_func_ch4.cpp unit_test_integrands_func_ch4.cpp

*/

#define UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT 0

static inline
bool 
almostEqual2sComplement(float A,float B, std::int32_t maxUlps) 
{
  // Make sure maxUlps is non-negative and small enough that the
  // default NAN won't compare as equal to anything.
  // This check disabled for now so that my tests run properly.
  // This assert should be enabled for normal use.
  // assert(maxUlps > 0 && maxUlps < 4 * 1024 * 1024);
  int aInt =  *reinterpret_cast<std::int32_t*>(&A);   //*(int *)&A;
  // Make aInt lexicographically ordered as a twos-complement int
  if (aInt < 0)
    aInt = 0x80000000 - aInt;
  // Make bInt lexicographically ordered as a twos-complement int
  int bInt =   *reinterpret_cast<std::int32_t*>(&B);  //*(int *)&B;
  if (bInt < 0)
    bInt = 0x80000000 - bInt;
  int intDiff = std::abs(aInt - bInt);
  if (intDiff <= maxUlps)
    return true;
  return false;
}

static inline 
std::int32_t  
printNumber(const char * val_name,float f, const std::int32_t offset) 
{
  (*((int *)&f)) += offset;
  std::int32_t printf_ret = std::printf("value: %s,%+1.11g,0x%08X,%d\n",val_name,f, *(int *)&f, *(int *)&f);
  return (printf_ret);
}

static inline 
std::int32_t 
print_double(const char * val_name,double f, const std::int32_t offset) 
{
  (*((std::uint64_t *)&f)) += offset;
  std::int32_t printf_ret = std::printf("value: %s,%1.17g,%017llx,%llu\n",val_name,f, *(std::uint64_t *)&f, *(std::uint64_t *)&f);
  return (printf_ret);
}

static inline
bool detect_subnormal(const float x)
{
    constexpr std::int32_t subnormal_code = 3;
    const std::int32_t subnormal_value = std::fpclassify(x);
    bool b_ret{};
    (subnormal_code==subnormal_value) ? b_ret = true : b_ret = false;
    return (b_ret);
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_integrand_4_1_gauss_Q_func();

void unit_test_integrand_4_1_gauss_Q_func()
{
    constexpr double lo_x{0.1};
    constexpr double hi_x{10.0};
    thread_local std::uniform_real_distribution<double> rv_func_arg;
    thread_local std::mt19937 rv_func_arg_gen;
    thread_local std::uint64_t seed_func_arg{};
    thread_local std::uniform_real_distribution<double> rv_gauss_q;
    thread_local std::mt19937 rv_gauss_q_gen;
    thread_local std::uint64_t seed_gauss_q{};
    double y=lo_x;
    constexpr std::int32_t n_func_args{10};
    constexpr std::int32_t n_gauss_q_vals{50};
    [[maybe_unused]] std::int32_t printf_ret{};
    
    rv_func_arg = std::uniform_real_distribution<double>(lo_x,hi_x);
    seed_func_arg = __rdtsc();
    rv_func_arg_gen = std::mt19937(seed_func_arg);
    rv_gauss_q = std::uniform_real_distribution<double>(y,5.0);
    seed_gauss_q = __rdtsc();
    rv_gauss_q_gen = std::mt19937(seed_gauss_q);
#if (INTEGRANDS_FUNC_CH4_USE_CEPHES_DOUBLE) == 1
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!! (IMPL: CEPHES-DOUBLE)\n",__func__);
    for(std::int32_t i {0}; i < n_func_args; ++i) 
    {
        const double x = rv_func_arg.operator()(rv_func_arg_gen);
        for(std::int32_t j{0}; j < n_gauss_q_vals; ++j)   
        {
            const double y = rv_gauss_q.operator()(rv_gauss_q_gen);
            const double q_gauss_res = gms::fading_channel::integrand_4_1_gauss_Q_func(x+y);
            printf_ret = print_double("q_gauss_result",q_gauss_res,0);
        }
    }
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!! (IMPL: CEPHES-DOUBLE)\n",__func__);
#else 
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!! (IMPL: STD::LIB -- CMATH)\n",__func__);
    for(std::int32_t i {0}; i < n_func_args; ++i) 
    {
        const double x = rv_func_arg.operator()(rv_func_arg_gen);
        for(std::int32_t j{0}; j < n_gauss_q_vals; ++j)   
        {
            const double y = rv_gauss_q.operator()(rv_gauss_q_gen);
            const double q_gauss_res = gms::fading_channel::integrand_4_1_gauss_Q_func(x+y);
            printf_ret = print_double("q_gauss_result",q_gauss_res,0);
        }
    }
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!! (IMPL: STD::LIB -- CMATH)\n",__func__);
#endif
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_integrand_4_2_gauss_Q_func();

void unit_test_integrand_4_2_gauss_Q_func()
{
    constexpr double lo_x{std::numeric_limits<double>::epsilon()};
    constexpr double hi_x{1.0};
    thread_local std::uniform_real_distribution<double> rv_func_arg;
    thread_local std::mt19937 rv_func_arg_gen;
    thread_local std::uint64_t seed_func_arg{};
    thread_local std::uniform_real_distribution<double> rv_gauss_q;
    thread_local std::mt19937 rv_gauss_q_gen;
    thread_local std::uint64_t seed_gauss_q{};
    constexpr std::int32_t n_func_args{10};
    constexpr std::int32_t n_gauss_q_vals{50};
    [[maybe_unused]] std::int32_t printf_ret{};
    
    rv_func_arg = std::uniform_real_distribution<double>(lo_x,hi_x);
    seed_func_arg = __rdtsc();
    rv_func_arg_gen = std::mt19937(seed_func_arg);
    rv_gauss_q = std::uniform_real_distribution<double>(0.0,1.57079632679489662);
    seed_gauss_q = __rdtsc();
    rv_gauss_q_gen = std::mt19937(seed_gauss_q);
#if (INTEGRANDS_FUNC_CH4_USE_CEPHES_DOUBLE) == 1
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!! (IMPL: CEPHES-DOUBLE)\n",__func__);
    for(std::int32_t i {0}; i < n_func_args; ++i) 
    {
        const double x = rv_func_arg.operator()(rv_func_arg_gen);
        for(std::int32_t j{0}; j < n_gauss_q_vals; ++j)   
        {
            const double y = rv_gauss_q.operator()(rv_gauss_q_gen);
            const double q_gauss_res = gms::fading_channel::integrand_4_2_gauss_Q_func(x,y);
            printf_ret = print_double("q_gauss_result",q_gauss_res,0);
        }
    }
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!! (IMPL: CEPHES-DOUBLE)\n",__func__);
#else 
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!! (IMPL: STD::LIB -- CMATH)\n",__func__);
    for(std::int32_t i {0}; i < n_func_args; ++i) 
    {
        const double x = rv_func_arg.operator()(rv_func_arg_gen);
        for(std::int32_t j{0}; j < n_gauss_q_vals; ++j)   
        {
            const double y = rv_gauss_q.operator()(rv_gauss_q_gen);
            const double q_gauss_res = gms::fading_channel::integrand_4_2_gauss_Q_func(x,y);
            printf_ret = print_double("q_gauss_result",q_gauss_res,0);
        }
    }
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!! (IMPL: STD::LIB -- CMATH)\n",__func__);
#endif
}

#include "GMS_cephes_double.h"

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_integrand_4_6_gauss_Q_func();

void unit_test_integrand_4_6_gauss_Q_func()
{
    constexpr double lo_x1{0.001};
    constexpr double hi_x1{1.0};
    constexpr double lo_y1{lo_x1};
    constexpr double hi_y1{hi_x1};
    constexpr double lo_rho{0.0};
    constexpr double hi_rho{1.0};
    thread_local std::uniform_real_distribution<double> rv_func_arg_x1;
    thread_local std::mt19937 rv_func_arg_x1_gen;
    thread_local std::uint64_t seed_func_x1_arg{};
    thread_local std::uniform_real_distribution<double> rv_func_arg_y1;
    thread_local std::mt19937 rv_func_arg_y1_gen;
    thread_local std::uint64_t seed_func_y1_arg{};
    thread_local std::uniform_real_distribution<double> rv_func_arg_rho;
    thread_local std::mt19937 rv_func_arg_rho_gen;
    thread_local std::uint64_t seed_func_rho_arg{};
    thread_local std::uniform_real_distribution<double> rv_gauss_q;
    thread_local std::mt19937 rv_gauss_q_gen;
    thread_local std::uint64_t seed_gauss_q{};
    constexpr std::int32_t n_func_args{10};
    constexpr std::int32_t n_gauss_q_vals{50};
    [[maybe_unused]] std::int32_t printf_ret{};
    
    rv_func_arg_x1 = std::uniform_real_distribution<double>(lo_x1,hi_x1);
    seed_func_x1_arg = __rdtsc();
    rv_func_arg_x1_gen = std::mt19937(seed_func_x1_arg);
    rv_func_arg_y1 = std::uniform_real_distribution<double>(lo_y1,hi_y1);
    seed_func_y1_arg = __rdtsc();
    rv_func_arg_y1_gen = std::mt19937(seed_func_y1_arg);
    rv_func_arg_rho = std::uniform_real_distribution<double>(lo_rho,hi_rho);
    seed_func_rho_arg = __rdtsc();
    rv_func_arg_rho_gen = std::mt19937(seed_func_rho_arg);
#if (INTEGRANDS_FUNC_CH4_USE_CEPHES_DOUBLE) == 1
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!! (IMPL: CEPHES-DOUBLE)\n",__func__);
    for(std::int32_t i {0}; i < n_func_args; ++i) 
    {
        const double x1 = rv_func_arg_x1.operator()(rv_func_arg_x1_gen);
        const double y1 = rv_func_arg_y1.operator()(rv_func_arg_y1_gen);
        const double rho= rv_func_arg_rho.operator()(rv_func_arg_rho_gen);
        const double phi_arg = y1/x1;
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 1
        printf_ret = print_double("phi_arg",phi_arg,0);
#endif 
        const double phi_s = gms::math::cephes_d::atan(phi_arg);
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 1
        printf_ret = print_double("phi_s",phi_s,0);
#endif 
        rv_gauss_q = std::uniform_real_distribution<double>(0.0,1.57079632679489662-phi_s);
        seed_gauss_q = __rdtsc();
        rv_gauss_q_gen = std::mt19937(seed_gauss_q);
        for(std::int32_t j{0}; j < n_gauss_q_vals; ++j)   
        {
            const double theta = rv_gauss_q.operator()(rv_gauss_q_gen);
            const double q_gauss_res = gms::fading_channel::integrand_4_6_gauss_Q_func(x1,y1,rho,theta,2);
            printf_ret = print_double("q_gauss_2D_result",q_gauss_res,0);
        }
    }
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!! (IMPL: CEPHES-DOUBLE)\n",__func__);
#else 
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!! (IMPL: STD::LIB -- CMATH)\n",__func__);
    for(std::int32_t i {0}; i < n_func_args; ++i) 
    {
        const double x1 = rv_func_arg_x1.operator()(rv_func_arg_x1_gen);
        const double y1 = rv_func_arg_y1.operator()(rv_func_arg_y1_gen);
        const double rho= rv_func_arg_rho.operator()(rv_func_arg_rho_gen);
        const double phi_arg = y1/x1;
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 1
        printf_ret = print_double("phi_arg",phi_arg,0);
#endif 
        const double phi_s = std::atan(phi_arg);
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 1
        printf_ret = print_double("phi_arg",phi_arg,0);
#endif 
        rv_gauss_q = std::uniform_real_distribution<double>(0.0,1.57079632679489662-phi_s);
        seed_gauss_q = __rdtsc();
        rv_gauss_q_gen = std::mt19937(seed_gauss_q);
        for(std::int32_t j{0}; j < n_gauss_q_vals; ++j)   
        {
            const double theta = rv_gauss_q.operator()(rv_gauss_q_gen);
            const double q_gauss_res = gms::fading_channel::integrand_4_6_gauss_Q_func(x1,y1,rho,theta,1);
            printf_ret = print_double("q_gauss_2D_result",q_gauss_res,0)
        }
    }
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!! (IMPL: STD::LIB -- CMATH)\n",__func__);
#endif
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_integrand_4_7_gauss_Q_func();

void unit_test_integrand_4_7_gauss_Q_func()
{
    constexpr double lo_x1{0.001};
    constexpr double hi_x1{1.0};
    constexpr double lo_y1{lo_x1};
    constexpr double hi_y1{hi_x1};
    constexpr double lo_rho{0.0};
    constexpr double hi_rho{1.0};
    thread_local std::uniform_real_distribution<double> rv_func_arg_x1;
    thread_local std::mt19937 rv_func_arg_x1_gen;
    thread_local std::uint64_t seed_func_x1_arg{};
    thread_local std::uniform_real_distribution<double> rv_func_arg_y1;
    thread_local std::mt19937 rv_func_arg_y1_gen;
    thread_local std::uint64_t seed_func_y1_arg{};
    thread_local std::uniform_real_distribution<double> rv_func_arg_rho;
    thread_local std::mt19937 rv_func_arg_rho_gen;
    thread_local std::uint64_t seed_func_rho_arg{};
    thread_local std::uniform_real_distribution<double> rv_gauss_q;
    thread_local std::mt19937 rv_gauss_q_gen;
    thread_local std::uint64_t seed_gauss_q{};
    constexpr std::int32_t n_func_args{10};
    constexpr std::int32_t n_gauss_q_vals{50};
    [[maybe_unused]] std::int32_t printf_ret{};
    
    rv_func_arg_x1 = std::uniform_real_distribution<double>(lo_x1,hi_x1);
    seed_func_x1_arg = __rdtsc();
    rv_func_arg_x1_gen = std::mt19937(seed_func_x1_arg);
    rv_func_arg_y1 = std::uniform_real_distribution<double>(lo_y1,hi_y1);
    seed_func_y1_arg = __rdtsc();
    rv_func_arg_y1_gen = std::mt19937(seed_func_y1_arg);
    rv_func_arg_rho = std::uniform_real_distribution<double>(lo_rho,hi_rho);
    seed_func_rho_arg = __rdtsc();
    rv_func_arg_rho_gen = std::mt19937(seed_func_rho_arg);
#if (INTEGRANDS_FUNC_CH4_USE_CEPHES_DOUBLE) == 1
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!,(IMPL: CEPHES-DOUBLE),arg=x1\n",__func__);
    for(std::int32_t i {0}; i < n_func_args; ++i) 
    {
        const double x1 = rv_func_arg_x1.operator()(rv_func_arg_x1_gen);
        const double y1 = rv_func_arg_y1.operator()(rv_func_arg_y1_gen);
        const double rho= rv_func_arg_rho.operator()(rv_func_arg_rho_gen);
        const double phi_arg = y1/x1; 
        const double phi_s = gms::math::cephes_d::atan(phi_arg);
        rv_gauss_q = std::uniform_real_distribution<double>(0.0,1.57079632679489662-phi_s);
        seed_gauss_q = __rdtsc();
        rv_gauss_q_gen = std::mt19937(seed_gauss_q);
        for(std::int32_t j{0}; j < n_gauss_q_vals; ++j)   
        {
            const double theta = rv_gauss_q.operator()(rv_gauss_q_gen);
            const double q_gauss_res = gms::fading_channel::integrand_4_7_gauss_Q_func(x1,y1,rho,theta,1);
            printf_ret = print_double("q_gauss_2D_result",q_gauss_res,0);
        }
    }
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!! (IMPL: CEPHES-DOUBLE),arg=x1\n",__func__);
#else 
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!! (IMPL: STD::LIB-CMATH),arg=x1\n",__func__);
    for(std::int32_t i {0}; i < n_func_args; ++i) 
    {
        const double x1 = rv_func_arg_x1.operator()(rv_func_arg_x1_gen);
        const double y1 = rv_func_arg_y1.operator()(rv_func_arg_y1_gen);
        const double rho= rv_func_arg_rho.operator()(rv_func_arg_rho_gen);
        const double phi_arg = y1/x1;
        const double phi_s = std::atan(phi_arg); 
        rv_gauss_q = std::uniform_real_distribution<double>(0.0,1.57079632679489662-phi_s);
        seed_gauss_q = __rdtsc();
        rv_gauss_q_gen = std::mt19937(seed_gauss_q);
        for(std::int32_t j{0}; j < n_gauss_q_vals; ++j)   
        {
            const double theta = rv_gauss_q.operator()(rv_gauss_q_gen);
            const double q_gauss_res = gms::fading_channel::integrand_4_7_gauss_Q_func(x1,y1,rho,theta,1);
            printf_ret = print_double("q_gauss_2D_result",q_gauss_res,0)
        }
    }
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!! (IMPL: STD::LIB-CMATH),arg=x1\n",__func__);
#endif
#if (INTEGRANDS_FUNC_CH4_USE_CEPHES_DOUBLE) == 1
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!,(IMPL: CEPHES-DOUBLE),arg=y1\n",__func__);
    for(std::int32_t i {0}; i < n_func_args; ++i) 
    {
        const double x1 = rv_func_arg_x1.operator()(rv_func_arg_x1_gen);
        const double y1 = rv_func_arg_y1.operator()(rv_func_arg_y1_gen);
        const double rho= rv_func_arg_rho.operator()(rv_func_arg_rho_gen);
        const double phi_arg = y1/x1; 
        const double phi_s = gms::math::cephes_d::atan(phi_arg);
        rv_gauss_q = std::uniform_real_distribution<double>(0.0,1.57079632679489662-phi_s);
        seed_gauss_q = __rdtsc();
        rv_gauss_q_gen = std::mt19937(seed_gauss_q);
        for(std::int32_t j{0}; j < n_gauss_q_vals; ++j)   
        {
            const double theta = rv_gauss_q.operator()(rv_gauss_q_gen);
            const double q_gauss_res = gms::fading_channel::integrand_4_7_gauss_Q_func(x1,y1,rho,theta,2);
            printf_ret = print_double("q_gauss_2D_result",q_gauss_res,0);
        }
    }
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!! (IMPL: CEPHES-DOUBLE),arg=y1\n",__func__);
#else 
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!! (IMPL: STD::LIB-CMATH),arg=y1\n",__func__);
    for(std::int32_t i {0}; i < n_func_args; ++i) 
    {
        const double x1 = rv_func_arg_x1.operator()(rv_func_arg_x1_gen);
        const double y1 = rv_func_arg_y1.operator()(rv_func_arg_y1_gen);
        const double rho= rv_func_arg_rho.operator()(rv_func_arg_rho_gen);
        const double phi_arg = y1/x1;
        const double phi_s = std::atan(phi_arg); 
        rv_gauss_q = std::uniform_real_distribution<double>(0.0,1.57079632679489662-phi_s);
        seed_gauss_q = __rdtsc();
        rv_gauss_q_gen = std::mt19937(seed_gauss_q);
        for(std::int32_t j{0}; j < n_gauss_q_vals; ++j)   
        {
            const double theta = rv_gauss_q.operator()(rv_gauss_q_gen);
            const double q_gauss_res = gms::fading_channel::integrand_4_7_gauss_Q_func(x1,y1,rho,theta,2);
            printf_ret = print_double("q_gauss_2D_result",q_gauss_res,0)
        }
    }
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!! (IMPL: STD::LIB-CMATH),arg=y1\n",__func__);
#endif
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_integrand_4_8_gauss_Q_func();

void unit_test_integrand_4_8_gauss_Q_func()
{
    constexpr double lo_x1{0.001};
    constexpr double hi_x1{1.0};
    constexpr double lo_y1{lo_x1};
    constexpr double hi_y1{hi_x1};
    thread_local std::uniform_real_distribution<double> rv_func_arg_x1;
    thread_local std::mt19937 rv_func_arg_x1_gen;
    thread_local std::uint64_t seed_func_x1_arg{};
    thread_local std::uniform_real_distribution<double> rv_func_arg_y1;
    thread_local std::mt19937 rv_func_arg_y1_gen;
    thread_local std::uint64_t seed_func_y1_arg{};
    thread_local std::uniform_real_distribution<double> rv_gauss_q;
    thread_local std::mt19937 rv_gauss_q_gen;
    thread_local std::uint64_t seed_gauss_q{};
    constexpr std::int32_t n_func_args{10};
    constexpr std::int32_t n_gauss_q_vals{50};
    [[maybe_unused]] std::int32_t printf_ret{};
    rv_func_arg_x1 = std::uniform_real_distribution<double>(lo_x1,hi_x1);
    seed_func_x1_arg = __rdtsc();
    rv_func_arg_x1_gen = std::mt19937(seed_func_x1_arg);
    rv_func_arg_y1 = std::uniform_real_distribution<double>(lo_y1,hi_y1);
    seed_func_y1_arg = __rdtsc();
    rv_func_arg_y1_gen = std::mt19937(seed_func_y1_arg);
#if (INTEGRANDS_FUNC_CH4_USE_CEPHES_DOUBLE) == 1
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!,(IMPL: CEPHES-DOUBLE),arg=x1\n",__func__);
    for(std::int32_t i {0}; i < n_func_args; ++i) 
    {
        const double x1 = rv_func_arg_x1.operator()(rv_func_arg_x1_gen);
        const double y1 = rv_func_arg_y1.operator()(rv_func_arg_y1_gen);
        const double phi_arg = y1/x1; 
        const double phi_s = gms::math::cephes_d::atan(phi_arg);
        rv_gauss_q = std::uniform_real_distribution<double>(0.0,1.57079632679489662-phi_s);
        seed_gauss_q = __rdtsc();
        rv_gauss_q_gen = std::mt19937(seed_gauss_q);
        for(std::int32_t j{0}; j < n_gauss_q_vals; ++j)   
        {
            const double theta = rv_gauss_q.operator()(rv_gauss_q_gen);
            const double q_gauss_res = gms::fading_channel::integrand_4_8_gauss_Q_func(x1,y1,theta,1);
            printf_ret = print_double("q_gauss_2D_result",q_gauss_res,0);
        }
    }
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!! (IMPL: CEPHES-DOUBLE),arg=x1\n",__func__);
#else 
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!! (IMPL: STD::LIB-CMATH),arg=x1\n",__func__);
    for(std::int32_t i {0}; i < n_func_args; ++i) 
    {
        const double x1 = rv_func_arg_x1.operator()(rv_func_arg_x1_gen);
        const double y1 = rv_func_arg_y1.operator()(rv_func_arg_y1_gen);
        const double phi_arg = y1/x1;
        const double phi_s = std::atan(phi_arg); 
        rv_gauss_q = std::uniform_real_distribution<double>(0.0,1.57079632679489662-phi_s);
        seed_gauss_q = __rdtsc();
        rv_gauss_q_gen = std::mt19937(seed_gauss_q);
        for(std::int32_t j{0}; j < n_gauss_q_vals; ++j)   
        {
            const double theta = rv_gauss_q.operator()(rv_gauss_q_gen);
            const double q_gauss_res = gms::fading_channel::integrand_4_8_gauss_Q_func(x1,y1,theta,1);
            printf_ret = print_double("q_gauss_2D_result",q_gauss_res,0)
        }
    }
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!! (IMPL: STD::LIB-CMATH),arg=x1\n",__func__);
#endif
#if (INTEGRANDS_FUNC_CH4_USE_CEPHES_DOUBLE) == 1
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!,(IMPL: CEPHES-DOUBLE),arg=y1\n",__func__);
    for(std::int32_t i {0}; i < n_func_args; ++i) 
    {
        const double x1 = rv_func_arg_x1.operator()(rv_func_arg_x1_gen);
        const double y1 = rv_func_arg_y1.operator()(rv_func_arg_y1_gen);
        const double phi_arg = y1/x1; 
        const double phi_s = gms::math::cephes_d::atan(phi_arg);
        rv_gauss_q = std::uniform_real_distribution<double>(0.0,1.57079632679489662-phi_s);
        seed_gauss_q = __rdtsc();
        rv_gauss_q_gen = std::mt19937(seed_gauss_q);
        for(std::int32_t j{0}; j < n_gauss_q_vals; ++j)   
        {
            const double theta = rv_gauss_q.operator()(rv_gauss_q_gen);
            const double q_gauss_res = gms::fading_channel::integrand_4_8_gauss_Q_func(x1,y1,theta,2);
            printf_ret = print_double("q_gauss_2D_result",q_gauss_res,0);
        }
    }
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!! (IMPL: CEPHES-DOUBLE),arg=y1\n",__func__);
#else 
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!! (IMPL: STD::LIB-CMATH),arg=y1\n",__func__);
    for(std::int32_t i {0}; i < n_func_args; ++i) 
    {
        const double x1 = rv_func_arg_x1.operator()(rv_func_arg_x1_gen);
        const double y1 = rv_func_arg_y1.operator()(rv_func_arg_y1_gen);
        const double phi_arg = y1/x1;
        const double phi_s = std::atan(phi_arg); 
        rv_gauss_q = std::uniform_real_distribution<double>(0.0,1.57079632679489662-phi_s);
        seed_gauss_q = __rdtsc();
        rv_gauss_q_gen = std::mt19937(seed_gauss_q);
        for(std::int32_t j{0}; j < n_gauss_q_vals; ++j)   
        {
            const double theta = rv_gauss_q.operator()(rv_gauss_q_gen);
            const double q_gauss_res = gms::fading_channel::integrand_4_8_gauss_Q_func(x1,y1,theta,2);
            printf_ret = print_double("q_gauss_2D_result",q_gauss_res,0)
        }
    }
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!! (IMPL: STD::LIB-CMATH),arg=y1\n",__func__);
#endif
}

int main()
{   
    (void)unit_test_integrand_4_1_gauss_Q_func();
    (void)unit_test_integrand_4_2_gauss_Q_func();
    (void)unit_test_integrand_4_6_gauss_Q_func();
    (void)unit_test_integrand_4_7_gauss_Q_func();
    (void)unit_test_integrand_4_8_gauss_Q_func();
    return 0;
}