#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <cmath>
#include <type_traits>
#include "GMS_analytic_bep_sep_ch8.h"

/*
   icpc -o unit_test_analytic_bep_sep_ch8 -O3 -fasm-blocks -fp-model fast=2 -fno-exceptions -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_cephes_double.h GMS_fast_simd_funcs_approx.h GMS_gaussian_Q_approx.h GMS_marcum_Q_approx.h GMS_analytic_bep_sep_ch8.h unit_test_analytic_bep_sep_ch8.cpp
   ASM: 
   icpc -S -O3 -fverbose-asm -masm=intel -fno-exceptions -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_cephes_double.h GMS_fast_sim_func_approx.h GMS_gaussian_Q_approx.h GMS_marcum_Q_approx.h GMS_analytic_bep_sep_ch8.h unit_test_analytic_bep_sep_ch8.cpp

*/

static inline 
std::int32_t  
print_float(const char * val_name,float f, const std::int32_t offset) 
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

template<typename T>
bool is_number_subnormal(const T x)
{
    constexpr std::int32_t subnormal_code = 3;
    const std::int32_t subnormal_value = std::fpclassify(x);
    bool b_ret{};
    (subnormal_code==subnormal_value) ? b_ret = true : b_ret = false;
    return (b_ret);
}

/*
__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_BPSK_param_a_8_61();

void unit_test_BPSK_param_a_8_61()
{
    constexpr double Bl = 100000.0;
    constexpr double Tb = 0.000001;
    constexpr double lo_Ac{0.1};
    constexpr double hi_Ac{1.0};
    constexpr double lo_N0{0.001};
    constexpr double hi_N0{1.5};
    constexpr double Ts=Tb;
    constexpr double M = 2.0;
    thread_local std::uniform_real_distribution<double> rv_Ac_arg;
    thread_local std::mt19937 rv_Ac_arg_gen;
    thread_local std::uint64_t seed_Ac_arg{};
    thread_local std::uniform_real_distribution<double> rv_N0_arg;
    thread_local std::mt19937 rv_N0_arg_gen;
    thread_local std::uint64_t seed_N0_arg{};
    [[maybe_unused]] std::int32_t printf_ret{};
    rv_Ac_arg = std::uniform_real_distribution<double>(lo_Ac,hi_Ac);
    seed_Ac_arg = __rdtsc();
    rv_Ac_arg_gen = std::mt19937(seed_Ac_arg);
    rv_N0_arg = std::uniform_real_distribution<double>(lo_N0,hi_N0);
    seed_N0_arg = __rdtsc();
    rv_N0_arg_gen = std::mt19937(seed_N0_arg);
    const double rv_Ac = rv_Ac_arg.operator()(rv_Ac_arg_gen);
    const double rv_N0 = rv_N0_arg.operator()(rv_N0_arg_gen);
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!\n",__func__);
    const double bpsk_param_a = gms::fading_channel::BPSK_param_a_8_61<double>(rv_Ac,Ts,M,rv_N0,Bl,Tb);
    if(is_number_subnormal<double>(bpsk_param_a))
    {
        std::printf("[UNIT-TEST]: detected subnormal in: %s, of value=%.17f\n",__func__,bpsk_param_a);
    }
    printf_ret = print_double("bpsk_param_a",bpsk_param_a,0);
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}
*/



template<typename T> 
__attribute__((aligned(32)))
void unit_test_BPSK_param_a_8_61(const T lo1,const T hi1,
                                 const T lo2,const T hi2)
{
    constexpr T Bl = 100000.0;
    constexpr T Tb = static_cast<T>(0.000001);
    const T lo_Ac{lo1}; // 0.1
    const T hi_Ac{hi1}; // 1.0
    const T lo_N0{lo2}; // 0.001
    const T hi_N0{hi2}; // 1.5
    constexpr T Ts=Tb;
    constexpr T M = 2.0;
    [[maybe_unused]] std::int32_t printf_ret{};
    thread_local std::uniform_real_distribution<T> rv_Ac_arg;
    thread_local std::mt19937 rv_Ac_arg_gen;
    thread_local std::uint64_t seed_Ac_arg{};
    thread_local std::uniform_real_distribution<T> rv_N0_arg;
    thread_local std::mt19937 rv_N0_arg_gen;
    thread_local std::uint64_t seed_N0_arg{};
    rv_Ac_arg = std::uniform_real_distribution<T>(lo_Ac,hi_Ac);
    seed_Ac_arg = __rdtsc();
    rv_Ac_arg_gen = std::mt19937(seed_Ac_arg);
    rv_N0_arg = std::uniform_real_distribution<T>(lo_N0,hi_N0);
    seed_N0_arg = __rdtsc();
    rv_N0_arg_gen = std::mt19937(seed_N0_arg);
    const T rv_Ac = rv_Ac_arg.operator()(rv_Ac_arg_gen);
    const T rv_N0 = rv_N0_arg.operator()(rv_N0_arg_gen);
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!\n",__func__);
    const T bpsk_param_a = gms::fading_channel::BPSK_param_a_8_61<T>(rv_Ac,Ts,M,rv_N0,Bl,Tb);
    if(is_number_subnormal<T>(bpsk_param_a))
    {
        std::printf("[UNIT-TEST]: detected subnormal in: %s, of value=%.17f\n",__func__,bpsk_param_a);
    }
    if constexpr(std::is_same_v<T,float>)
    {
        printf_ret = print_float("bpsk_param_a",bpsk_param_a,0);
    }
    else 
    {
    printf_ret = print_double("bpsk_param_a",bpsk_param_a,0);
    }
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

template<typename T> 
__attribute__((aligned(32)))
void unit_test_BPSK_param_b_8_61(const T lo1,const T hi1,
                                 const T lo2,const T hi2)
{
    constexpr T Bl = 100000.0;
    constexpr T Tb = static_cast<T>(0.000001);
    const T lo_Ac{lo1}; // 0.1
    const T hi_Ac{hi1}; // 1.0
    const T lo_N0{lo2}; // 0.001
    const T hi_N0{hi2}; // 1.5
    constexpr T Ts=Tb;
    constexpr T M = 2.0;
    [[maybe_unused]] std::int32_t printf_ret{};
    thread_local std::uniform_real_distribution<T> rv_Ac_arg;
    thread_local std::mt19937 rv_Ac_arg_gen;
    thread_local std::uint64_t seed_Ac_arg{};
    thread_local std::uniform_real_distribution<T> rv_N0_arg;
    thread_local std::mt19937 rv_N0_arg_gen;
    thread_local std::uint64_t seed_N0_arg{};
    rv_Ac_arg = std::uniform_real_distribution<T>(lo_Ac,hi_Ac);
    seed_Ac_arg = __rdtsc();
    rv_Ac_arg_gen = std::mt19937(seed_Ac_arg);
    rv_N0_arg = std::uniform_real_distribution<T>(lo_N0,hi_N0);
    seed_N0_arg = __rdtsc();
    rv_N0_arg_gen = std::mt19937(seed_N0_arg);
    const T rv_Ac = rv_Ac_arg.operator()(rv_Ac_arg_gen);
    const T rv_N0 = rv_N0_arg.operator()(rv_N0_arg_gen);
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!\n",__func__);
    const T bpsk_param_b = gms::fading_channel::BPSK_param_b_8_61<T>(rv_Ac,Ts,M,rv_N0,Bl,Tb);
    if(is_number_subnormal<T>(bpsk_param_b))
    {
        std::printf("[UNIT-TEST]: detected subnormal in: %s, of value=%.17f\n",__func__,bpsk_param_b);
    }
    if constexpr(std::is_same_v<T,float>)
    {
        printf_ret = print_float("bpsk_param_b",bpsk_param_b,0);
    }
    else 
    {
    printf_ret = print_double("bpsk_param_a",bpsk_param_b,0);
    }
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

template<typename T> 
__attribute__((aligned(32)))
void unit_test_OQPSK_param_a2_8_64(const T lo1,const T hi1,
                                 const T lo2,const T hi2)
{
    constexpr T Bl = 100000.0;
    constexpr T Tb = static_cast<T>(0.000001);
    const T lo_Ac{lo1}; // 0.1
    const T hi_Ac{hi1}; // 1.0
    const T lo_N0{lo2}; // 0.001
    const T hi_N0{hi2}; // 1.5
    constexpr T Ts=Tb;
    constexpr T M = 2.0;
    [[maybe_unused]] std::int32_t printf_ret{};
    thread_local std::uniform_real_distribution<T> rv_Ac_arg;
    thread_local std::mt19937 rv_Ac_arg_gen;
    thread_local std::uint64_t seed_Ac_arg{};
    thread_local std::uniform_real_distribution<T> rv_N0_arg;
    thread_local std::mt19937 rv_N0_arg_gen;
    thread_local std::uint64_t seed_N0_arg{};
    rv_Ac_arg = std::uniform_real_distribution<T>(lo_Ac,hi_Ac);
    seed_Ac_arg = __rdtsc();
    rv_Ac_arg_gen = std::mt19937(seed_Ac_arg);
    rv_N0_arg = std::uniform_real_distribution<T>(lo_N0,hi_N0);
    seed_N0_arg = __rdtsc();
    rv_N0_arg_gen = std::mt19937(seed_N0_arg);
    const T rv_Ac = rv_Ac_arg.operator()(rv_Ac_arg_gen);
    const T rv_N0 = rv_N0_arg.operator()(rv_N0_arg_gen);
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!\n",__func__);
    const T oqpsk_param_a2 = gms::fading_channel::OQPSK_param_a2_8_64<T>(rv_Ac,Ts,M,rv_N0,Bl,Tb);
    if(is_number_subnormal<T>(oqpsk_param_a2))
    {
        std::printf("[UNIT-TEST]: detected subnormal in: %s, of value=%.17f\n",__func__,oqpsk_param_a2);
    }
    if constexpr(std::is_same_v<T,float>)
    {
        printf_ret = print_float("oqpsk_param_a2",oqpsk_param_a2,0);
    }
    else 
    {
    printf_ret = print_double("oqpsk_param_a2",oqpsk_param_a2,0);
    }
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

template<typename T> 
__attribute__((aligned(32)))
void unit_test_OQPSK_param_b2_8_64(const T lo1,const T hi1,
                                 const T lo2,const T hi2)
{
    constexpr T Bl = 100000.0;
    constexpr T Tb = static_cast<T>(0.000001);
    const T lo_Ac{lo1}; // 0.1
    const T hi_Ac{hi1}; // 1.0
    const T lo_N0{lo2}; // 0.001
    const T hi_N0{hi2}; // 1.5
    constexpr T Ts=Tb;
    constexpr T M = 2.0;
    [[maybe_unused]] std::int32_t printf_ret{};
    thread_local std::uniform_real_distribution<T> rv_Ac_arg;
    thread_local std::mt19937 rv_Ac_arg_gen;
    thread_local std::uint64_t seed_Ac_arg{};
    thread_local std::uniform_real_distribution<T> rv_N0_arg;
    thread_local std::mt19937 rv_N0_arg_gen;
    thread_local std::uint64_t seed_N0_arg{};
    rv_Ac_arg = std::uniform_real_distribution<T>(lo_Ac,hi_Ac);
    seed_Ac_arg = __rdtsc();
    rv_Ac_arg_gen = std::mt19937(seed_Ac_arg);
    rv_N0_arg = std::uniform_real_distribution<T>(lo_N0,hi_N0);
    seed_N0_arg = __rdtsc();
    rv_N0_arg_gen = std::mt19937(seed_N0_arg);
    const T rv_Ac = rv_Ac_arg.operator()(rv_Ac_arg_gen);
    const T rv_N0 = rv_N0_arg.operator()(rv_N0_arg_gen);
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!\n",__func__);
    const T oqpsk_param_b2 = gms::fading_channel::OQPSK_param_b2_8_64<T>(rv_Ac,Ts,M,rv_N0,Bl,Tb);
    if(is_number_subnormal<T>(oqpsk_param_b2))
    {
        std::printf("[UNIT-TEST]: detected subnormal in: %s, of value=%.17f\n",__func__,oqpsk_param_b2);
    }
    if constexpr(std::is_same_v<T,float>)
    {
        printf_ret = print_float("oqpsk_param_b2",oqpsk_param_b2,0);
    }
    else 
    {
    printf_ret = print_double("oqpsk_param_b2",oqpsk_param_b2,0);
    }
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

template<typename T> 
__attribute__((aligned(32)))
void unit_test_MSK_param_a2_8_65(const T lo1,const T hi1,
                                 const T lo2,const T hi2)
{
    constexpr T Bl = 100000.0;
    constexpr T Tb = static_cast<T>(0.000001);
    const T lo_Ac{lo1}; // 0.1
    const T hi_Ac{hi1}; // 1.0
    const T lo_N0{lo2}; // 0.001
    const T hi_N0{hi2}; // 1.5
    constexpr T Ts=Tb;
    constexpr T M = 2.0;
    [[maybe_unused]] std::int32_t printf_ret{};
    thread_local std::uniform_real_distribution<T> rv_Ac_arg;
    thread_local std::mt19937 rv_Ac_arg_gen;
    thread_local std::uint64_t seed_Ac_arg{};
    thread_local std::uniform_real_distribution<T> rv_N0_arg;
    thread_local std::mt19937 rv_N0_arg_gen;
    thread_local std::uint64_t seed_N0_arg{};
    rv_Ac_arg = std::uniform_real_distribution<T>(lo_Ac,hi_Ac);
    seed_Ac_arg = __rdtsc();
    rv_Ac_arg_gen = std::mt19937(seed_Ac_arg);
    rv_N0_arg = std::uniform_real_distribution<T>(lo_N0,hi_N0);
    seed_N0_arg = __rdtsc();
    rv_N0_arg_gen = std::mt19937(seed_N0_arg);
    const T rv_Ac = rv_Ac_arg.operator()(rv_Ac_arg_gen);
    const T rv_N0 = rv_N0_arg.operator()(rv_N0_arg_gen);
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!\n",__func__);
    const T msk_param_a2 = gms::fading_channel::MSK_param_a2_8_65<T>(rv_Ac,Ts,M,rv_N0,Bl,Tb);
    if(is_number_subnormal<T>(msk_param_a2))
    {
        std::printf("[UNIT-TEST]: detected subnormal in: %s, of value=%.17f\n",__func__,msk_param_a2);
    }
    if constexpr(std::is_same_v<T,float>)
    {
        printf_ret = print_float("msk_param_a2",msk_param_a2,0);
    }
    else 
    {
    printf_ret = print_double("msk_param_a2",msk_param_a2,0);
    }
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

template<typename T> 
__attribute__((aligned(32)))
void unit_test_MSK_param_b2_8_65(const T lo1,const T hi1,
                                 const T lo2,const T hi2)
{
    constexpr T Bl = 100000.0;
    constexpr T Tb = static_cast<T>(0.000001);
    const T lo_Ac{lo1}; // 0.1
    const T hi_Ac{hi1}; // 1.0
    const T lo_N0{lo2}; // 0.001
    const T hi_N0{hi2}; // 1.5
    constexpr T Ts=Tb;
    constexpr T M = 2.0;
    [[maybe_unused]] std::int32_t printf_ret{};
    thread_local std::uniform_real_distribution<T> rv_Ac_arg;
    thread_local std::mt19937 rv_Ac_arg_gen;
    thread_local std::uint64_t seed_Ac_arg{};
    thread_local std::uniform_real_distribution<T> rv_N0_arg;
    thread_local std::mt19937 rv_N0_arg_gen;
    thread_local std::uint64_t seed_N0_arg{};
    rv_Ac_arg = std::uniform_real_distribution<T>(lo_Ac,hi_Ac);
    seed_Ac_arg = __rdtsc();
    rv_Ac_arg_gen = std::mt19937(seed_Ac_arg);
    rv_N0_arg = std::uniform_real_distribution<T>(lo_N0,hi_N0);
    seed_N0_arg = __rdtsc();
    rv_N0_arg_gen = std::mt19937(seed_N0_arg);
    const T rv_Ac = rv_Ac_arg.operator()(rv_Ac_arg_gen);
    const T rv_N0 = rv_N0_arg.operator()(rv_N0_arg_gen);
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!\n",__func__);
    const T msk_param_b2 = gms::fading_channel::MSK_param_b2_8_65<T>(rv_Ac,Ts,M,rv_N0,Bl,Tb);
    if(is_number_subnormal<T>(msk_param_b2))
    {
        std::printf("[UNIT-TEST]: detected subnormal in: %s, of value=%.17f\n",__func__,msk_param_b2);
    }
    if constexpr(std::is_same_v<T,float>)
    {
        printf_ret = print_float("msk_param_b2",msk_param_b2,0);
    }
    else 
    {
    printf_ret = print_double("msk_param_b2",msk_param_b2,0);
    }
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

 static thread_local std::uniform_real_distribution<double> g_rv_Ac_arg;
 static thread_local std::mt19937 g_rv_Ac_arg_gen;
 static thread_local std::uint64_t g_seed_Ac_arg{};
 static thread_local std::uniform_real_distribution<double> g_rv_N0_arg;
 static thread_local std::mt19937 g_rv_N0_arg_gen;
 static thread_local std::uint64_t g_seed_N0_arg{};

 static thread_local std::uniform_real_distribution<float> g_rv_Ac_arg_f;
 static thread_local std::mt19937 g_rv_Ac_arg_gen_f;
 static thread_local std::uint64_t g_seed_Ac_arg_f{};
 static thread_local std::uniform_real_distribution<float> g_rv_N0_arg_f;
 static thread_local std::mt19937 g_rv_N0_arg_gen_f;
 static thread_local std::uint64_t g_seed_N0_arg_f{};


 __attribute__((constructor(545)))
 bool init_g_rand_generators();

 bool init_g_rand_generators()
 {
    constexpr double Ac_lo{0.1};
    constexpr double Ac_hi{5.0};
    constexpr double N0_lo{0.01};
    constexpr double N0_hi{4.568};
    thread_local std::uniform_real_distribution<double> rv_Ac_range;
    thread_local std::mt19937 rv_Ac_range_gen;
    thread_local std::uint64_t seed_Ac_range{};
    thread_local std::uniform_real_distribution<double> rv_N0_range;
    thread_local std::mt19937 rv_N0_range_gen;
    thread_local std::uint64_t seed_N0_range{};
    bool ret_state = false;
    rv_Ac_range = std::uniform_real_distribution<double>(Ac_lo,Ac_hi);
    seed_Ac_range = __rdtsc();
    rv_Ac_range_gen = std::mt19937(seed_Ac_range);
    rv_N0_range = std::uniform_real_distribution<double>(N0_lo,N0_hi);
    seed_N0_range = __rdtsc();
    rv_N0_range_gen = std::mt19937(seed_N0_range);
    double rv_Ac_val_lo;
    double rv_Ac_val_hi;
    double rv_N0_val_lo;
    double rv_N0_val_hi;
    bool   is_Ac_lo_greater{};
    bool   is_N0_lo_greater{};
    do 
    {
        rv_Ac_val_lo = rv_Ac_range.operator()(rv_Ac_range_gen);
        rv_Ac_val_hi = rv_Ac_range.operator()(rv_Ac_range_gen);
        rv_N0_val_lo = rv_N0_range.operator()(rv_N0_range_gen);
        rv_N0_val_hi = rv_N0_range.operator()(rv_N0_range_gen);
        is_Ac_lo_greater = rv_Ac_val_lo>rv_Ac_val_hi;
        is_N0_lo_greater = rv_N0_val_lo>rv_N0_val_hi;
    }while(is_Ac_lo_greater && is_N0_lo_greater);

    g_rv_Ac_arg = std::uniform_real_distribution<double>(rv_Ac_val_lo,rv_Ac_val_hi);
    g_seed_Ac_arg = __rdtsc();
    g_rv_Ac_arg_gen = std::mt19937(g_seed_Ac_arg);
    g_rv_N0_arg = std::uniform_real_distribution<double>(rv_N0_val_lo,rv_N0_val_hi);
    g_seed_N0_arg = __rdtsc();
    g_rv_N0_arg_gen = std::mt19937(g_seed_N0_arg);
    ret_state = true;
    return (ret_state);
 }

 __attribute__((constructor(546)))
 bool init_g_rand_generators_f();

 bool init_g_rand_generators_f()
 {
    constexpr float Ac_lo{0.1f};
    constexpr float Ac_hi{5.0f};
    constexpr float N0_lo{0.01f};
    constexpr float N0_hi{4.568f};
    thread_local std::uniform_real_distribution<float> rv_Ac_range;
    thread_local std::mt19937 rv_Ac_range_gen;
    thread_local std::uint64_t seed_Ac_range{};
    thread_local std::uniform_real_distribution<float> rv_N0_range;
    thread_local std::mt19937 rv_N0_range_gen;
    thread_local std::uint64_t seed_N0_range{};
    bool ret_state = false;
    rv_Ac_range = std::uniform_real_distribution<float>(Ac_lo,Ac_hi);
    seed_Ac_range = __rdtsc();
    rv_Ac_range_gen = std::mt19937(seed_Ac_range);
    rv_N0_range = std::uniform_real_distribution<float>(N0_lo,N0_hi);
    seed_N0_range = __rdtsc();
    rv_N0_range_gen = std::mt19937(seed_N0_range);
    float rv_Ac_val_lo;
    float rv_Ac_val_hi;
    float rv_N0_val_lo;
    float rv_N0_val_hi;
    bool   is_Ac_lo_greater{};
    bool   is_N0_lo_greater{};
    do 
    {
        rv_Ac_val_lo = rv_Ac_range.operator()(rv_Ac_range_gen);
        rv_Ac_val_hi = rv_Ac_range.operator()(rv_Ac_range_gen);
        rv_N0_val_lo = rv_N0_range.operator()(rv_N0_range_gen);
        rv_N0_val_hi = rv_N0_range.operator()(rv_N0_range_gen);
        is_Ac_lo_greater = rv_Ac_val_lo>rv_Ac_val_hi;
        is_N0_lo_greater = rv_N0_val_lo>rv_N0_val_hi;
    }while(is_Ac_lo_greater && is_N0_lo_greater);

    g_rv_Ac_arg_f = std::uniform_real_distribution<float>(rv_Ac_val_lo,rv_Ac_val_hi);
    g_seed_Ac_arg_f = __rdtsc();
    g_rv_Ac_arg_gen_f = std::mt19937(g_seed_Ac_arg);
    g_rv_N0_arg_f = std::uniform_real_distribution<float>(rv_N0_val_lo,rv_N0_val_hi);
    g_seed_N0_arg_f = __rdtsc();
    g_rv_N0_arg_gen_f = std::mt19937(g_seed_N0_arg);
    ret_state = true;
    return (ret_state);
 }

#define DOUBLE_FUNC_BODY_BLOCK(name1,name2,name3,name4,name5,name6,name7,name8)\
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!\n",__func__);\
    const double sep_mam_1 = (name1);\
    if(is_number_subnormal<double>(sep_mam_1))\
    {\
        printf_ret = std::printf("[UNIT-TEST]: detected subnormal in: %s, of value=%.17f\n",__func__,sep_mam_1);\
    }\
    printf_ret = print_double("sep_mam_1",sep_mam_1,0);\
    printf_ret = std::printf("===========================================================\n");\
    const double sep_mam_2 = (name2);\
    if(is_number_subnormal<double>(sep_mam_2))\
    {\
        printf_ret = std::printf("[UNIT-TEST]: detected subnormal in: %s, of value=%.17f\n",__func__,sep_mam_2);\
    }\
    printf_ret = print_double("sep_mam_2",sep_mam_2,0);\
    printf_ret = std::printf("===========================================================\n");\
    const double sep_mam_3 = (name3);\
    if(is_number_subnormal<double>(sep_mam_3))\
    {\
        printf_ret = std::printf("[UNIT-TEST]: detected subnormal in: %s, of value=%.17f\n",__func__,sep_mam_3);\
    }\
    printf_ret = print_double("sep_mam_3",sep_mam_3,0);\
    printf_ret = std::printf("===========================================================\n");\
    const double sep_mam_4 = (name4);\
    if(is_number_subnormal<double>(sep_mam_4))\
    {\
        printf_ret = std::printf("[UNIT-TEST]: detected subnormal in: %s, of value=%.17f\n",__func__,sep_mam_4);\
    }\
    printf_ret = print_double("sep_mam_4",sep_mam_4,0);\
    printf_ret = std::printf("===========================================================\n");\
    const double sep_mam_5 = (name5);\
    if(is_number_subnormal<double>(sep_mam_5))\
    {\
        printf_ret = std::printf("[UNIT-TEST]: detected subnormal in: %s, of value=%.17f\n",__func__,sep_mam_5);\
    }\
    printf_ret = print_double("sep_mam_5",sep_mam_5,0);\
    printf_ret = std::printf("===========================================================\n");\
    const double sep_mam_6 = (name6);\
    if(is_number_subnormal<double>(sep_mam_6))\
    {\
        printf_ret = std::printf("[UNIT-TEST]: detected subnormal in: %s, of value=%.17f\n",__func__,sep_mam_6);\
    }\
    printf_ret = print_double("sep_mam_6",sep_mam_6,0);\
    printf_ret = std::printf("===========================================================\n");\
    const double sep_mam_7 = (name7);\
    if(is_number_subnormal<double>(sep_mam_7))\
    {\
        printf_ret = std::printf("[UNIT-TEST]: detected subnormal in: %s, of value=%.17f\n",__func__,sep_mam_7);\
    }\
    printf_ret = print_double("sep_mam_7",sep_mam_7,0);\
    printf_ret = std::printf("===========================================================\n");\
    const double sep_mam_8 = (name8);\
    if(is_number_subnormal<double>(sep_mam_8))\
    {\
        printf_ret = std::printf("[UNIT-TEST]: detected subnormal in: %s, of value=%.17f\n",__func__,sep_mam_8);\
    }\
    printf_ret = print_double("sep_mam_8",sep_mam_8,0);\
    printf_ret = std::printf("===========================================================\n");\
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n\n",__func__);

#define FLOAT_FUNC_BODY_BLOCK(name1,name2,name3,name4,name5,name6,name7,name8)\
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!\n",__func__);\
    const float sep_mam_1 = (name1);\
    if(is_number_subnormal<float>(sep_mam_1))\
    {\
        printf_ret = std::printf("[UNIT-TEST]: detected subnormal in: %s, of value=%.7f\n",__func__,sep_mam_1);\
    }\
    printf_ret = print_float("sep_mam_1",sep_mam_1,0);\
    printf_ret = std::printf("===========================================================\n");\
    const float sep_mam_2 = (name2);\
    if(is_number_subnormal<float>(sep_mam_2))\
    {\
        printf_ret = std::printf("[UNIT-TEST]: detected subnormal in: %s, of value=%.7f\n",__func__,sep_mam_2);\
    }\
    printf_ret = print_float("sep_mam_2",sep_mam_2,0);\
    printf_ret = std::printf("===========================================================\n");\
    const float sep_mam_3 = (name3);\
    if(is_number_subnormal<float>(sep_mam_3))\
    {\
        printf_ret = std::printf("[UNIT-TEST]: detected subnormal in: %s, of value=%.7f\n",__func__,sep_mam_3);\
    }\
    printf_ret = print_float("sep_mam_3",sep_mam_3,0);\
    printf_ret = std::printf("===========================================================\n");\
    const float sep_mam_4 = (name4);\
    if(is_number_subnormal<float>(sep_mam_4))\
    {\
        printf_ret = std::printf("[UNIT-TEST]: detected subnormal in: %s, of value=%.7f\n",__func__,sep_mam_4);\
    }\
    printf_ret = print_float("sep_mam_4",sep_mam_4,0);\
    printf_ret = std::printf("===========================================================\n");\
    const float sep_mam_5 = (name5);\
    if(is_number_subnormal<float>(sep_mam_5))\
    {\
        printf_ret = std::printf("[UNIT-TEST]: detected subnormal in: %s, of value=%.7f\n",__func__,sep_mam_5);\
    }\
    printf_ret = print_float("sep_mam_5",sep_mam_5,0);\
    printf_ret = std::printf("===========================================================\n");\
    const float sep_mam_6 = (name6);\
    if(is_number_subnormal<float>(sep_mam_6))\
    {\
        printf_ret = std::printf("[UNIT-TEST]: detected subnormal in: %s, of value=%.7f\n",__func__,sep_mam_6);\
    }\
    printf_ret = print_float("sep_mam_6",sep_mam_6,0);\
    printf_ret = std::printf("===========================================================\n");\
    const float sep_mam_7 = (name7);\
    if(is_number_subnormal<float>(sep_mam_7))\
    {\
        printf_ret = std::printf("[UNIT-TEST]: detected subnormal in: %s, of value=%.7f\n",__func__,sep_mam_7);\
    }\
    printf_ret = print_float("sep_mam_7",sep_mam_7,0);\
    printf_ret = std::printf("===========================================================\n");\
    const float sep_mam_8 = (name8);\
    if(is_number_subnormal<float>(sep_mam_8))\
    {\
        printf_ret = std::printf("[UNIT-TEST]: detected subnormal in: %s, of value=%.7f\n",__func__,sep_mam_8);\
    }\
    printf_ret = print_float("sep_mam_8",sep_mam_8,0);\
    printf_ret = std::printf("===========================================================\n");\
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n\n",__func__);

__attribute__((aligned(32)))
void unit_test_analytic_SEP_MAM_8_1(const double,const std::int32_t);
 
void unit_test_analytic_SEP_MAM_8_1(const double M,const std::int32_t n)      
{
    constexpr double Ts = 0.000001;
    [[maybe_unused]] std::int32_t printf_ret{};
    const double rv_Ac = g_rv_Ac_arg.operator()(g_rv_Ac_arg_gen);
    const double rv_N0 = g_rv_N0_arg.operator()(g_rv_N0_arg_gen);
    DOUBLE_FUNC_BODY_BLOCK(gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani>(rv_Ac,Ts,M,rv_N0,n),
                           gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T>(rv_Ac,Ts,M,rv_N0,n),
                           gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T>(rv_Ac,Ts,M,rv_N0,n),
                           gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T>(rv_Ac,Ts,M,rv_N0,n),
                           gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T>(rv_Ac,Ts,M,rv_N0,n),
                           gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T>(rv_Ac,Ts,M,rv_N0,n),
                           gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson>(rv_Ac,Ts,M,rv_N0,n),
                           gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed>(rv_Ac,Ts,M,rv_N0,n))
                           
}

__attribute__((aligned(32)))
void unit_test_analytic_SEP_MAM_8_1(const float,const std::int32_t);
 
void unit_test_analytic_SEP_MAM_8_1(const float M,const std::int32_t n)      
{
    constexpr double Ts = 0.000001f;
    [[maybe_unused]] std::int32_t printf_ret{};
    const float rv_Ac = g_rv_Ac_arg.operator()(g_rv_Ac_arg_gen);
    const float rv_N0 = g_rv_N0_arg.operator()(g_rv_N0_arg_gen);
    FLOAT_FUNC_BODY_BLOCK(gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani>(rv_Ac,Ts,M,rv_N0,n),
                           gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T>(rv_Ac,Ts,M,rv_N0,n),
                           gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T>(rv_Ac,Ts,M,rv_N0,n),
                           gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T>(rv_Ac,Ts,M,rv_N0,n),
                           gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T>(rv_Ac,Ts,M,rv_N0,n),
                           gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T>(rv_Ac,Ts,M,rv_N0,n),
                           gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson>(rv_Ac,Ts,M,rv_N0,n),
                           gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed>(rv_Ac,Ts,M,rv_N0,n))

}

__attribute__((aligned(32)))
void unit_test_analytic_SEP_MAM_8_3(const double,const std::int32_t);
 
void unit_test_analytic_SEP_MAM_8_3(const double M,const std::int32_t n)      
{
    constexpr double Ts = 0.000001;
    [[maybe_unused]] std::int32_t printf_ret{};
    const double rv_Ac = g_rv_Ac_arg.operator()(g_rv_Ac_arg_gen);
    const double rv_N0 = g_rv_N0_arg.operator()(g_rv_N0_arg_gen);
    DOUBLE_FUNC_BODY_BLOCK(gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani>(rv_Ac,Ts,M,rv_N0,n),
                           gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T>(rv_Ac,Ts,M,rv_N0,n),
                           gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T>(rv_Ac,Ts,M,rv_N0,n),
                           gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T>(rv_Ac,Ts,M,rv_N0,n),
                           gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T>(rv_Ac,Ts,M,rv_N0,n),
                           gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T>(rv_Ac,Ts,M,rv_N0,n),
                           gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson>(rv_Ac,Ts,M,rv_N0,n),
                           gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed>(rv_Ac,Ts,M,rv_N0,n))
                           
}

__attribute__((aligned(32)))
void unit_test_analytic_SEP_MAM_8_3(const float,const std::int32_t);
 
void unit_test_analytic_SEP_MAM_8_3(const float M,const std::int32_t n)      
{
    constexpr double Ts = 0.000001f;
    [[maybe_unused]] std::int32_t printf_ret{};
    const float rv_Ac = g_rv_Ac_arg.operator()(g_rv_Ac_arg_gen);
    const float rv_N0 = g_rv_N0_arg.operator()(g_rv_N0_arg_gen);
    FLOAT_FUNC_BODY_BLOCK(gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_chiani>(rv_Ac,Ts,M,rv_N0,n),
                           gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_2T>(rv_Ac,Ts,M,rv_N0,n),
                           gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_loskot_3T>(rv_Ac,Ts,M,rv_N0,n),
                           gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_1T>(rv_Ac,Ts,M,rv_N0,n),
                           gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_2T>(rv_Ac,Ts,M,rv_N0,n),
                           gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_4T>(rv_Ac,Ts,M,rv_N0,n),
                           gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_borjesson>(rv_Ac,Ts,M,rv_N0,n),
                           gms::fading_channel::analytic_SEP_MAM_8_1<gms::fading_channel::Gaussian_Q_approxmations_t::Gaussian_Q_approx_sadhwani_summed>(rv_Ac,Ts,M,rv_N0,n))

}

int main()
{
/*
   (void)unit_test_BPSK_param_a_8_61<double>(1.0,5.0,0.1,1.0);
   (void)unit_test_BPSK_param_a_8_61<float>(1.0f,5.0f,0.1f,1.0f);
   (void)unit_test_BPSK_param_b_8_61<double>(1.0,4.5,0.01,0.9);
   (void)unit_test_BPSK_param_b_8_61<float>(1.0f,4.5f,0.01f,1.0f);
   (void)unit_test_OQPSK_param_a2_8_64<double>(1.0,5.0,0.1,1.0);
   (void)unit_test_OQPSK_param_a2_8_64<float>(1.0f,4.5f,0.01f,1.0f);
   (void)unit_test_OQPSK_param_b2_8_64<double>(1.0,5.0,0.1,1.0);
   (void)unit_test_OQPSK_param_b2_8_64<float>(1.0f,4.5f,0.01f,1.0f);
   (void)unit_test_MSK_param_a2_8_65<double>(1.0,5.0,0.1,1.0);
   (void)unit_test_MSK_param_a2_8_65<float>(1.0f,4.5f,0.01f,1.0f);
   (void)unit_test_MSK_param_b2_8_65<double>(1.0,5.0,0.1,1.0);
   (void)unit_test_MSK_param_b2_8_65<float>(1.0f,4.5f,0.01f,1.0f);
*/

   bool ret_g_init = init_g_rand_generators();
   if(!ret_g_init) 
   { 
     std::printf("[***FATAL***] in: File=%s,function=%s,line=%d\n",__FILE__,__func__,__LINE__);
     std::exit(EXIT_FAILURE);
   }
   bool ret_g_init_f = init_g_rand_generators_f();
   if(!ret_g_init_f)
   {
     std::printf("[***FATAL***] in: File=%s,function=%s,line=%d\n",__FILE__,__func__,__LINE__);
     std::exit(EXIT_FAILURE);
   }
   (void)unit_test_analytic_SEP_MAM_8_1(32.0,128);
   (void)unit_test_analytic_SEP_MAM_8_1(32.0f,128);
   (void)unit_test_analytic_SEP_MAM_8_3(32.0,128);
   (void)unit_test_analytic_SEP_MAM_8_3(32.0f,128);
   return (0);
}