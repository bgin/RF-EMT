#include <cstdio>
#include <limits>
#include <random>
#include "GMS_integrands_helpers.h"

/*
   icpc -o unit_test_integrands_helpers -O3 -fp-model fast=2 -fno-exceptions -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_integrands_helpers.h GMS_integrands_helpers.cpp unit_test_integrands_helpers.cpp
   ASM: 
   icpc -S -O3 -fverbose-asm -masm=intel -fno-exceptions -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_integrands_helpers.h GMS_integrands_helpers.cpp unit_test_integrands_helpers.cpp

*/

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_integrand_4_1_gauss_Q_arg_checker();

void unit_test_integrand_4_1_gauss_Q_arg_checker()
{
    constexpr const char * fut_name = "integrand_4_1_gauss_Q_arg_checker";
    thread_local std::uniform_real_distribution<double> rv_in_params;
    thread_local std::mt19937 rv_in_params_gen;
    thread_local std::uint64_t seed_rv_in_params{};
    rv_in_params = std::uniform_real_distribution<double>(0.1,1.5);
    seed_rv_in_params = __rdtsc();
    rv_in_params_gen = std::mt19937(seed_rv_in_params);
    const double lo1{rv_in_params.operator()(rv_in_params_gen)};
    const double hi1{rv_in_params.operator()(rv_in_params_gen)};
    const double hi2{1.0};
    constexpr std::int32_t n_func_args{10};
    constexpr std::int32_t n_gauss_q_vals{50};
    [[maybe_unused]] std::int32_t printf_ret{};
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!\n",__func__);
    std::int32_t ret_val = gms::fading_channel::integrand_4_1_gauss_Q_arg_checker(nullptr,n_func_args,n_gauss_q_vals,lo1,hi1,hi2);
    printf_ret = std::printf("[UNIT_TEST:] -- tested function=%s, returned value=%d\n",fut_name,ret_val);
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s  -- ENDED!!\n",__func__);
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_integrand_4_2_gauss_Q_arg_checker();

void unit_test_integrand_4_2_gauss_Q_arg_checker()
{
    constexpr const char * fut_name = "integrand_4_2_gauss_Q_arg_checker";
    thread_local std::uniform_real_distribution<double> rv_in_params1;
    thread_local std::mt19937 rv_in_params1_gen;
    thread_local std::uint64_t seed_rv_in_params1{};
    rv_in_params1 = std::uniform_real_distribution<double>(0.1,0.3);
    seed_rv_in_params1 = __rdtsc();
    rv_in_params1_gen = std::mt19937(seed_rv_in_params1);
    thread_local std::uniform_real_distribution<double> rv_in_params2;
    thread_local std::mt19937 rv_in_params2_gen;
    thread_local std::uint64_t seed_rv_in_params2{};
    rv_in_params2 = std::uniform_real_distribution<double>(+0.1,+1.5707963267948966192313216916398);
    seed_rv_in_params2 = __rdtsc();
    rv_in_params2_gen = std::mt19937(seed_rv_in_params2);
    const double lo1{rv_in_params1.operator()(rv_in_params1_gen)};
    const double hi1{rv_in_params1.operator()(rv_in_params1_gen)};
    const double lo2{rv_in_params1.operator()(rv_in_params2_gen)};
    const double hi2{rv_in_params1.operator()(rv_in_params2_gen)};
    constexpr std::int32_t n_func_args{10};
    constexpr std::int32_t n_gauss_q_vals{50};
    [[maybe_unused]] std::int32_t printf_ret{};
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!\n",__func__);
    std::int32_t ret_val = gms::fading_channel::integrand_4_2_gauss_Q_arg_checker(nullptr,n_func_args,n_gauss_q_vals,lo1,hi1,lo2,hi2);
    printf_ret = std::printf("[UNIT_TEST:] -- tested function=%s, returned value=%d\n",fut_name,ret_val);
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s  -- ENDED!!\n",__func__);
}

#define RANDOM_GEN_PARAMS_BLOCK\
    thread_local std::uniform_real_distribution<double> rv_in_params1;\
    thread_local std::mt19937 rv_in_params1_gen;\
    thread_local std::uint64_t seed_rv_in_params1{};\
    rv_in_params1 = std::uniform_real_distribution<double>(0.1,1.0);\
    seed_rv_in_params1 = __rdtsc();\
    rv_in_params1_gen = std::mt19937(seed_rv_in_params1);\
    thread_local std::uniform_real_distribution<double> rv_in_params2;\
    thread_local std::mt19937 rv_in_params2_gen;\
    thread_local std::uint64_t seed_rv_in_params2{};\
    rv_in_params2 = std::uniform_real_distribution<double>(0.2,+1.2);\
    seed_rv_in_params2 = __rdtsc();\
    rv_in_params2_gen = std::mt19937(seed_rv_in_params2);\
    thread_local std::uniform_real_distribution<double> rv_in_params3;\
    thread_local std::mt19937 rv_in_params3_gen;\
    thread_local std::uint64_t seed_rv_in_params3{};\
    rv_in_params3 = std::uniform_real_distribution<double>(0.1,0.5);\
    seed_rv_in_params3 = __rdtsc();\
    rv_in_params3_gen = std::mt19937(seed_rv_in_params3);\
    thread_local std::uniform_real_distribution<double> rv_in_params4;\
    thread_local std::mt19937 rv_in_params4_gen;\
    thread_local std::uint64_t seed_rv_in_params4{};\
    rv_in_params4 = std::uniform_real_distribution<double>(+0.001,+1.5707963267948966192313216916398);\
    seed_rv_in_params4 = __rdtsc();\
    rv_in_params4_gen = std::mt19937(seed_rv_in_params4);\
    const double lo1{rv_in_params1.operator()(rv_in_params1_gen)};\
    const double hi1{rv_in_params1.operator()(rv_in_params1_gen)};\
    const double lo2{rv_in_params2.operator()(rv_in_params2_gen)};\
    const double hi2{rv_in_params2.operator()(rv_in_params2_gen)};\
    const double lo3{rv_in_params3.operator()(rv_in_params3_gen)};\
    const double hi3{rv_in_params3.operator()(rv_in_params3_gen)};\
    const double lo4{rv_in_params4.operator()(rv_in_params4_gen)};\
    const double hi4{rv_in_params4.operator()(rv_in_params4_gen)};

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_integrand_4_6_sin_gauss_Q_arg_checker();

void unit_test_integrand_4_6_sin_gauss_Q_arg_checker()
{
    constexpr const char * fut_name = "integrand_4_6_sin_gauss_Q_arg_checker";
    RANDOM_GEN_PARAMS_BLOCK;
    constexpr std::int32_t n_func_args{10};
    constexpr std::int32_t n_gauss_q_vals{50};
    [[maybe_unused]] std::int32_t printf_ret{};
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!\n",__func__);
    std::int32_t ret_val = gms::fading_channel::integrand_4_6_gauss_Q_sin_arg_checker(nullptr,n_func_args,n_gauss_q_vals,lo1,hi1,
                                                                                      lo2,hi2,lo3,hi3,lo4,hi4);
    printf_ret = std::printf("[UNIT_TEST:] -- tested function=%s, returned value=%d\n",fut_name,ret_val);
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s  -- ENDED!!\n",__func__);
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_integrand_4_6_cos_gauss_Q_arg_checker();

void unit_test_integrand_4_6_cos_gauss_Q_arg_checker()
{
    constexpr const char * fut_name = "integrand_4_6_cos_gauss_Q_arg_checker";
    RANDOM_GEN_PARAMS_BLOCK;
    constexpr std::int32_t n_func_args{10};
    constexpr std::int32_t n_gauss_q_vals{50};
    [[maybe_unused]] std::int32_t printf_ret{};
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!\n",__func__);
    std::int32_t ret_val = gms::fading_channel::integrand_4_6_gauss_Q_cos_arg_checker(nullptr,n_func_args,n_gauss_q_vals,lo1,hi1,
                                                                                      lo2,hi2,lo3,hi3,lo4,hi4);
    printf_ret = std::printf("[UNIT_TEST:] -- tested function=%s, returned value=%d\n",fut_name,ret_val);
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s  -- ENDED!!\n",__func__);
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_integrand_4_10_marcum_Q_arg_checker();

void unit_test_integrand_4_10_marcum_Q_arg_checker()
{
    constexpr const char * fut_name = "integrand_4_10_marcum_Q_arg_checker";
    thread_local std::uniform_real_distribution<double> rv_in_params1;
    thread_local std::mt19937 rv_in_params1_gen;
    thread_local std::uint64_t seed_rv_in_params1{};
    rv_in_params1 = std::uniform_real_distribution<double>(0.1,10.3);
    seed_rv_in_params1 = __rdtsc();
    rv_in_params1_gen = std::mt19937(seed_rv_in_params1);
    thread_local std::uniform_real_distribution<double> rv_in_params2;
    thread_local std::mt19937 rv_in_params2_gen;
    thread_local std::uint64_t seed_rv_in_params2{};
    rv_in_params2 = std::uniform_real_distribution<double>(0.0,5.0);
    seed_rv_in_params2 = __rdtsc();
    rv_in_params2_gen = std::mt19937(seed_rv_in_params2);
    const double lo1{rv_in_params1.operator()(rv_in_params1_gen)};
    const double hi1{rv_in_params1.operator()(rv_in_params1_gen)};
    const double lo2{rv_in_params2.operator()(rv_in_params2_gen)};
    const double hi2{rv_in_params2.operator()(rv_in_params2_gen)};
    constexpr std::int32_t n_func_args{10};
    constexpr std::int32_t n_gauss_q_vals{50};
    [[maybe_unused]] std::int32_t printf_ret{};
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!\n",__func__);
    std::int32_t ret_val = gms::fading_channel::integrand_4_10_marcum_Q_arg_checker(nullptr,n_func_args,n_gauss_q_vals,lo1,hi1,lo2,hi2);
    printf_ret = std::printf("[UNIT_TEST:] -- tested function=%s, returned value=%d\n",fut_name,ret_val);
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s  -- ENDED!!\n",__func__);
}

int main()
{
    (void)unit_test_integrand_4_1_gauss_Q_arg_checker();
    (void)unit_test_integrand_4_2_gauss_Q_arg_checker();
    (void)unit_test_integrand_4_6_sin_gauss_Q_arg_checker();
    (void)unit_test_integrand_4_6_cos_gauss_Q_arg_checker();
    (void)unit_test_integrand_4_10_marcum_Q_arg_checker();
    return 0;
}