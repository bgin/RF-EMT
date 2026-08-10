#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <cmath>
#include <algorithm>
#include "GMS_integrands_func_ch4_sse.h"

/*
   icpc -o unit_test_integrands_func_ch4_sse -O3 -fp-model fast=2 -fno-exceptions -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_bessel_i0_sse.h GMS_integrands_func_ch4_sse.h unit_test_integrands_func_ch4_sse.cpp
   ASM: 
   icpc -S -O3 -fverbose-asm -masm=intel -fno-exceptions -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_simd_utils.h GMS_bessel_i0_sse.h GMS_integrands_func_ch4_sse.h unit_test_integrands_func_ch4_sse.cpp

*/

#define UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT 1

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
void unit_test_integrand_4_1_gauss_Q_func_sse_pd();

void unit_test_integrand_4_1_gauss_Q_func_sse_pd()
{
    constexpr std::int32_t n_func_args{10};
    constexpr std::int32_t n_gauss_q_vals{50};
    constexpr std::int32_t tot_elems{n_func_args*n_gauss_q_vals};
    constexpr double lo_x{0.1};
    constexpr double hi_x{10.0};
    __ATTR_ALIGN__(16) double in_buf[tot_elems];
    __ATTR_ALIGN__(16) double out_buf[tot_elems];
    thread_local std::uniform_real_distribution<double> rv_func_arg;
    thread_local std::mt19937 rv_func_arg_gen;
    thread_local std::uint64_t seed_func_arg{};
    thread_local std::uniform_real_distribution<double> rv_gauss_q;
    thread_local std::mt19937 rv_gauss_q_gen;
    thread_local std::uint64_t seed_gauss_q{};
    double y=lo_x;
    [[maybe_unused]] std::int32_t printf_ret{};
    rv_func_arg = std::uniform_real_distribution<double>(lo_x,hi_x);
    seed_func_arg = __rdtsc();
    rv_func_arg_gen = std::mt19937(seed_func_arg);
    rv_gauss_q = std::uniform_real_distribution<double>(y,5.0);
    seed_gauss_q = __rdtsc();
    rv_gauss_q_gen = std::mt19937(seed_gauss_q);

    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!.\n",__func__);
    for(std::int32_t i {0}; i < n_func_args; ++i) 
    {
        const std::int32_t outer_idx = i*n_gauss_q_vals;
        const double x = rv_func_arg.operator()(rv_func_arg_gen);
        for(std::int32_t j{0}; j < n_gauss_q_vals; ++j)   
        {
            const std::int32_t inner_idx = outer_idx+j;
            const double y = rv_gauss_q.operator()(rv_gauss_q_gen);
            in_buf[inner_idx] = x+y;
        }
    }
    std::sort(&in_buf[0],&in_buf[tot_elems-1]);
    for(std::int32_t i = 0; i != ROUND_DOWN(tot_elems,1); i += 2)
    {
        const __m128d vy{_mm_load_pd(&in_buf[i])};
        const __m128d gauss_q_res{gms::fading_channel::integrand_4_1_gauss_Q_func_sse_pd(vy)};
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 1
        const double * __restrict__ p_res{reinterpret_cast<const double * __restrict__>(&gauss_q_res)};
        printf_ret = print_double("p_res[0]",p_res[0],0);
        printf_ret = print_double("p_res[1]",p_res[1],0);
#endif 
        _mm_store_pd(&out_buf[i],gauss_q_res);
    }
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_integrand_4_2_gauss_Q_func_sse_pd();

void unit_test_integrand_4_2_gauss_Q_func_sse_pd()
{
    constexpr std::int32_t n_func_args{10};
    constexpr std::int32_t n_gauss_q_vals{50};
    constexpr std::int32_t tot_elems{n_func_args*n_gauss_q_vals};
    constexpr double lo_x{0.01};
    constexpr double hi_x{1.0};
    __ATTR_ALIGN__(16) double in_buf[tot_elems];
    __ATTR_ALIGN__(16) double out_buf[tot_elems];
    __ATTR_ALIGN__(16) double in_rv_func_buf[n_func_args];
    thread_local std::uniform_real_distribution<double> rv_func_arg;
    thread_local std::mt19937 rv_func_arg_gen;
    thread_local std::uint64_t seed_func_arg{};
    thread_local std::uniform_real_distribution<double> rv_gauss_q;
    thread_local std::mt19937 rv_gauss_q_gen;
    thread_local std::uint64_t seed_gauss_q{};
    [[maybe_unused]] std::int32_t printf_ret{};
    rv_func_arg = std::uniform_real_distribution<double>(lo_x,hi_x);
    seed_func_arg = __rdtsc();
    rv_func_arg_gen = std::mt19937(seed_func_arg);
    rv_gauss_q = std::uniform_real_distribution<double>(+0.0,+1.5707963267948966192313216916398);
    seed_gauss_q = __rdtsc();
    rv_gauss_q_gen = std::mt19937(seed_gauss_q);
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!.\n",__func__);
    for(std::int32_t i {0}; i < n_func_args; ++i) 
    {
        const std::int32_t outer_idx = i*n_gauss_q_vals;
        const double x = rv_func_arg.operator()(rv_func_arg_gen);
        in_rv_func_buf[i] = x;
        for(std::int32_t j{0}; j < n_gauss_q_vals; ++j)   
        {
            const std::int32_t inner_idx = outer_idx+j;
            const double y = rv_gauss_q.operator()(rv_gauss_q_gen);
            in_buf[inner_idx] = y;
        }
    }
    std::sort(&in_buf[0],&in_buf[tot_elems-1]);
    std::sort(&in_rv_func_buf[0],&in_rv_func_buf[n_func_args-1]);
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 1
     for(std::int32_t i {0}; i < tot_elems; ++i) 
     {  
        if(i>=0 && i<=(n_func_args-1)) {print_double("sorted in_rv_func_buf",in_rv_func_buf[i],0);}
        print_double("sorted: in_buf=",in_buf[i],0);
     }
#endif 
    for(std::int32_t i{0}; i < n_func_args; ++i) 
    {  
      const double x{in_rv_func_buf[i]};
      const std::int32_t outer_idx = i*n_gauss_q_vals;
      const __m128d vx{_mm_set1_pd(x)};
      for(std::int32_t j{0}; j < ROUND_DOWN(n_gauss_q_vals,1); j += 2) 
      {
          const std::int32_t inner_idx = outer_idx+j;
          const __m128d theta{_mm_load_pd(&in_buf[inner_idx])};
          const __m128d gauss_Q_res{gms::fading_channel::integrand_4_2_gauss_Q_func_sse_pd(vx,theta)};
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 1
          const double * __restrict__ p_res{reinterpret_cast<const double * __restrict__>(&gauss_Q_res)};
          printf_ret = print_double("p_res[0]",p_res[0],0);
          printf_ret = print_double("p_res[1]",p_res[1],0);
#endif 
        _mm_store_pd(&out_buf[inner_idx],gauss_Q_res);
      }
    } 
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_integrand_4_6_sin_gauss_Q_func_sse_pd();

void unit_test_integrand_4_6_sin_gauss_Q_func_sse_pd()
{
    constexpr std::int32_t n_func_args{10};
    constexpr std::int32_t n_gauss_q_vals{50};
    constexpr std::int32_t tot_elems{n_func_args*n_gauss_q_vals};
    constexpr double lo_x1{0.001};
    constexpr double hi_x1{1.0};
    constexpr double lo_y1{lo_x1};
    constexpr double hi_y1{hi_x1};
    constexpr double lo_rho{0.0};
    constexpr double hi_rho{1.0};
    __ATTR_ALIGN__(16) double in_buf[tot_elems];
    __ATTR_ALIGN__(16) double out_buf[tot_elems];
    __ATTR_ALIGN__(16) double rv_func_x1_buf[n_func_args];
    __ATTR_ALIGN__(16) double rv_func_y1_buf[n_func_args];
    __ATTR_ALIGN__(16) double rv_func_rho_buf[n_func_args];
    __ATTR_ALIGN__(16) double rv_phi_arg_buf[n_func_args];
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
    rv_func_arg_x1 = std::uniform_real_distribution<double>(lo_x1,hi_x1);
    seed_func_x1_arg = __rdtsc();
    rv_func_arg_x1_gen = std::mt19937(seed_func_x1_arg);
    rv_func_arg_y1 = std::uniform_real_distribution<double>(lo_y1,hi_y1);
    seed_func_y1_arg = __rdtsc();
    rv_func_arg_y1_gen = std::mt19937(seed_func_y1_arg);
    rv_func_arg_rho = std::uniform_real_distribution<double>(lo_rho,hi_rho);
    seed_func_rho_arg = __rdtsc();
    rv_func_arg_rho_gen = std::mt19937(seed_func_rho_arg);
    [[maybe_unused]] std::int32_t printf_ret{};
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!.\n",__func__);

    for(std::int32_t i {0}; i < n_func_args; ++i) 
    {
        const std::int32_t outer_idx = i*n_gauss_q_vals;
        const double x1 = rv_func_arg_x1.operator()(rv_func_arg_x1_gen);
        rv_func_x1_buf[i] = x1;
        const double y1 = rv_func_arg_y1.operator()(rv_func_arg_y1_gen);
        rv_func_y1_buf[i] = y1;
        const double rho= rv_func_arg_rho.operator()(rv_func_arg_rho_gen);
        rv_func_rho_buf[i] = rho;
        const double phi_arg = y1/x1;
        rv_phi_arg_buf[i] = phi_arg;
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 1
        printf_ret = print_double("phi_arg",phi_arg,0);
#endif 
        const double phi_s = std::atan(phi_arg);
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 1
        printf_ret = print_double("phi_s",phi_s,0);
#endif 
        rv_gauss_q = std::uniform_real_distribution<double>(+0.0,+1.57079632679489662-phi_s);
        seed_gauss_q = __rdtsc();
        rv_gauss_q_gen = std::mt19937(seed_gauss_q);
        for(std::int32_t j{0}; j < n_gauss_q_vals; ++j)   
        {
            const std::int32_t inner_idx = outer_idx+j;
            const double theta = rv_gauss_q.operator()(rv_gauss_q_gen);
            in_buf[inner_idx] = theta;
        }
    }
    std::sort(&in_buf[0],&in_buf[tot_elems-1]);
    std::sort(&rv_func_x1_buf[0],&rv_func_x1_buf[n_func_args-1]);
    std::sort(&rv_func_y1_buf[0],&rv_func_y1_buf[n_func_args-1]);
    std::sort(&rv_func_rho_buf[0],&rv_func_rho_buf[n_func_args-1]);
    std::sort(&rv_phi_arg_buf[0],&rv_phi_arg_buf[n_func_args-1]);

#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 0
     for(std::int32_t i {0}; i < tot_elems; ++i) 
     {  
        if(i>=0 && i<=(n_func_args-1)) 
        {
          print_double("sort-check: rv_func_x1_buf",rv_func_x1_buf[i],0);
          print_double("sort-check: rv_func_y1_buf",rv_func_y1_buf[i],0);
          print_double("sort-check: rv_func_rho_buf",rv_func_rho_buf[i],0);
          print_double("sort-check: rv_phi_arg_buf",rv_phi_arg_buf[i],0);
        }
        print_double("sort-check: in_buf=",in_buf[i],0);
     }
#endif 

    for(std::int32_t i{0}; i < n_func_args; ++i) 
    {  
      const std::int32_t outer_idx = i*n_gauss_q_vals;
      const double x1{rv_func_x1_buf[i]};
      const __m128d vx1{_mm_set1_pd(x1)};
      const double y1{rv_func_y1_buf[i]};
      const __m128d vy1{_mm_set1_pd(y1)};
      const double rho{rv_func_rho_buf[i]};
      const __m128d vrho{_mm_set1_pd(rho)};
      const double phiarg{rv_phi_arg_buf[i]};
      const __m128d vphis{_mm_set1_pd(phiarg)};
      for(std::int32_t j{0}; j < ROUND_DOWN(n_gauss_q_vals,1); j += 2) 
      {
          const std::int32_t inner_idx = outer_idx+j;
          const __m128d theta{_mm_load_pd(&in_buf[inner_idx])};
          const __m128d gauss_Q_res{gms::fading_channel::integrand_4_6_sin_gauss_Q_func_sse_pd(vx1,vy1,vrho,theta,vphis)};
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 1
          const double * __restrict__ p_res{reinterpret_cast<const double * __restrict__>(&gauss_Q_res)};
          printf_ret = print_double("p_res[0]",p_res[0],0);
          printf_ret = print_double("p_res[1]",p_res[1],0);
#endif 
        _mm_store_pd(&out_buf[inner_idx],gauss_Q_res);
      }
    } 
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_integrand_4_6_cos_gauss_Q_func_sse_pd();

void unit_test_integrand_4_6_cos_gauss_Q_func_sse_pd()
{
    constexpr std::int32_t n_func_args{10};
    constexpr std::int32_t n_gauss_q_vals{50};
    constexpr std::int32_t tot_elems{n_func_args*n_gauss_q_vals};
    constexpr double lo_x1{0.001};
    constexpr double hi_x1{1.0};
    constexpr double lo_y1{lo_x1};
    constexpr double hi_y1{hi_x1};
    constexpr double lo_rho{0.0};
    constexpr double hi_rho{1.0};
    __ATTR_ALIGN__(16) double in_buf[tot_elems];
    __ATTR_ALIGN__(16) double out_buf[tot_elems];
    __ATTR_ALIGN__(16) double rv_func_x1_buf[n_func_args];
    __ATTR_ALIGN__(16) double rv_func_y1_buf[n_func_args];
    __ATTR_ALIGN__(16) double rv_func_rho_buf[n_func_args];
    __ATTR_ALIGN__(16) double rv_phi_arg_buf[n_func_args];
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
    rv_func_arg_x1 = std::uniform_real_distribution<double>(lo_x1,hi_x1);
    seed_func_x1_arg = __rdtsc();
    rv_func_arg_x1_gen = std::mt19937(seed_func_x1_arg);
    rv_func_arg_y1 = std::uniform_real_distribution<double>(lo_y1,hi_y1);
    seed_func_y1_arg = __rdtsc();
    rv_func_arg_y1_gen = std::mt19937(seed_func_y1_arg);
    rv_func_arg_rho = std::uniform_real_distribution<double>(lo_rho,hi_rho);
    seed_func_rho_arg = __rdtsc();
    rv_func_arg_rho_gen = std::mt19937(seed_func_rho_arg);
    [[maybe_unused]] std::int32_t printf_ret{};
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!.\n",__func__);

    for(std::int32_t i {0}; i < n_func_args; ++i) 
    {
        const std::int32_t outer_idx = i*n_gauss_q_vals;
        const double x1 = rv_func_arg_x1.operator()(rv_func_arg_x1_gen);
        rv_func_x1_buf[i] = x1;
        const double y1 = rv_func_arg_y1.operator()(rv_func_arg_y1_gen);
        rv_func_y1_buf[i] = y1;
        const double rho= rv_func_arg_rho.operator()(rv_func_arg_rho_gen);
        rv_func_rho_buf[i] = rho;
        const double phi_arg = y1/x1;
        rv_phi_arg_buf[i] = phi_arg;
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 1
        printf_ret = print_double("phi_arg",phi_arg,0);
#endif 
        const double phi_s = std::atan(phi_arg);
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 1
        printf_ret = print_double("phi_s",phi_s,0);
#endif 
        rv_gauss_q = std::uniform_real_distribution<double>(+0.0,+1.57079632679489662-phi_s);
        seed_gauss_q = __rdtsc();
        rv_gauss_q_gen = std::mt19937(seed_gauss_q);
        for(std::int32_t j{0}; j < n_gauss_q_vals; ++j)   
        {
            const std::int32_t inner_idx = outer_idx+j;
            const double theta = rv_gauss_q.operator()(rv_gauss_q_gen);
            in_buf[inner_idx] = theta;
        }
    }
    std::sort(&in_buf[0],&in_buf[tot_elems-1]);
    std::sort(&rv_func_x1_buf[0],&rv_func_x1_buf[n_func_args-1]);
    std::sort(&rv_func_y1_buf[0],&rv_func_y1_buf[n_func_args-1]);
    std::sort(&rv_func_rho_buf[0],&rv_func_rho_buf[n_func_args-1]);
    std::sort(&rv_phi_arg_buf[0],&rv_phi_arg_buf[n_func_args-1]);

#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 0
     for(std::int32_t i {0}; i < tot_elems; ++i) 
     {  
        if(i>=0 && i<=(n_func_args-1)) 
        {
          print_double("sort-check: rv_func_x1_buf",rv_func_x1_buf[i],0);
          print_double("sort-check: rv_func_y1_buf",rv_func_y1_buf[i],0);
          print_double("sort-check: rv_func_rho_buf",rv_func_rho_buf[i],0);
          print_double("sort-check: rv_phi_arg_buf",rv_phi_arg_buf[i],0);
        }
        print_double("sort-check: in_buf=",in_buf[i],0);
     }
#endif 

    for(std::int32_t i{0}; i < n_func_args; ++i) 
    {  
      const std::int32_t outer_idx = i*n_gauss_q_vals;
      const double x1{rv_func_x1_buf[i]};
      const __m128d vx1{_mm_set1_pd(x1)};
      const double y1{rv_func_y1_buf[i]};
      const __m128d vy1{_mm_set1_pd(y1)};
      const double rho{rv_func_rho_buf[i]};
      const __m128d vrho{_mm_set1_pd(rho)};
      const double phiarg{rv_phi_arg_buf[i]};
      const __m128d vphis{_mm_set1_pd(phiarg)};
      for(std::int32_t j{0}; j < ROUND_DOWN(n_gauss_q_vals,1); j += 2) 
      {
          const std::int32_t inner_idx = outer_idx+j;
          const __m128d theta{_mm_load_pd(&in_buf[inner_idx])};
          const __m128d gauss_Q_res{gms::fading_channel::integrand_4_6_cos_gauss_Q_func_sse_pd(vx1,vy1,vrho,theta,vphis)};
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 1
          const double * __restrict__ p_res{reinterpret_cast<const double * __restrict__>(&gauss_Q_res)};
          printf_ret = print_double("p_res[0]",p_res[0],0);
          printf_ret = print_double("p_res[1]",p_res[1],0);
#endif 
        _mm_store_pd(&out_buf[inner_idx],gauss_Q_res);
      }
    } 
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_integrand_4_7_x1_gauss_Q_func_sse_pd();

void unit_test_integrand_4_7_x1_gauss_Q_func_sse_pd()
{
    constexpr std::int32_t n_func_args{10};
    constexpr std::int32_t n_gauss_q_vals{50};
    constexpr std::int32_t tot_elems{n_func_args*n_gauss_q_vals};
    constexpr double lo_x1{0.001};
    constexpr double hi_x1{1.0};
    constexpr double lo_y1{lo_x1};
    constexpr double hi_y1{hi_x1};
    constexpr double lo_rho{0.0};
    constexpr double hi_rho{1.0};
    __ATTR_ALIGN__(16) double in_buf[tot_elems];
    __ATTR_ALIGN__(16) double out_buf[tot_elems];
    __ATTR_ALIGN__(16) double rv_func_x1_buf[n_func_args];
    __ATTR_ALIGN__(16) double rv_func_y1_buf[n_func_args];
    __ATTR_ALIGN__(16) double rv_func_rho_buf[n_func_args];
    __ATTR_ALIGN__(16) double rv_phi_arg_buf[n_func_args];
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
    rv_func_arg_x1 = std::uniform_real_distribution<double>(lo_x1,hi_x1);
    seed_func_x1_arg = __rdtsc();
    rv_func_arg_x1_gen = std::mt19937(seed_func_x1_arg);
    rv_func_arg_y1 = std::uniform_real_distribution<double>(lo_y1,hi_y1);
    seed_func_y1_arg = __rdtsc();
    rv_func_arg_y1_gen = std::mt19937(seed_func_y1_arg);
    rv_func_arg_rho = std::uniform_real_distribution<double>(lo_rho,hi_rho);
    seed_func_rho_arg = __rdtsc();
    rv_func_arg_rho_gen = std::mt19937(seed_func_rho_arg);
    [[maybe_unused]] std::int32_t printf_ret{};
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!.\n",__func__);

    for(std::int32_t i {0}; i < n_func_args; ++i) 
    {
        const std::int32_t outer_idx = i*n_gauss_q_vals;
        const double x1 = rv_func_arg_x1.operator()(rv_func_arg_x1_gen);
        rv_func_x1_buf[i] = x1;
        const double y1 = rv_func_arg_y1.operator()(rv_func_arg_y1_gen);
        rv_func_y1_buf[i] = y1;
        const double rho= rv_func_arg_rho.operator()(rv_func_arg_rho_gen);
        rv_func_rho_buf[i] = rho;
        const double phi_arg = y1/x1;
        rv_phi_arg_buf[i] = phi_arg;
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 1
        printf_ret = print_double("phi_arg",phi_arg,0);
#endif 
        const double phi_s = std::atan(phi_arg);
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 1
        printf_ret = print_double("phi_s",phi_s,0);
#endif 
        rv_gauss_q = std::uniform_real_distribution<double>(+0.0,+1.57079632679489662-phi_s);
        seed_gauss_q = __rdtsc();
        rv_gauss_q_gen = std::mt19937(seed_gauss_q);
        for(std::int32_t j{0}; j < n_gauss_q_vals; ++j)   
        {
            const std::int32_t inner_idx = outer_idx+j;
            const double theta = rv_gauss_q.operator()(rv_gauss_q_gen);
            in_buf[inner_idx] = theta;
        }
    }
    std::sort(&in_buf[0],&in_buf[tot_elems-1]);
    std::sort(&rv_func_x1_buf[0],&rv_func_x1_buf[n_func_args-1]);
    // Not used and left for verification.
    std::sort(&rv_func_y1_buf[0],&rv_func_y1_buf[n_func_args-1]);
    std::sort(&rv_func_rho_buf[0],&rv_func_rho_buf[n_func_args-1]);
    // Not used and left for verification.
    std::sort(&rv_phi_arg_buf[0],&rv_phi_arg_buf[n_func_args-1]);

#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 0
     for(std::int32_t i {0}; i < tot_elems; ++i) 
     {  
        if(i>=0 && i<=(n_func_args-1)) 
        {
          print_double("sort-check: rv_func_x1_buf",rv_func_x1_buf[i],0);
          print_double("sort-check: rv_func_y1_buf",rv_func_y1_buf[i],0);
          print_double("sort-check: rv_func_rho_buf",rv_func_rho_buf[i],0);
          print_double("sort-check: rv_phi_arg_buf",rv_phi_arg_buf[i],0);
        }
        print_double("sort-check: in_buf=",in_buf[i],0);
     }
#endif 

    for(std::int32_t i{0}; i < n_func_args; ++i) 
    {  
      const std::int32_t outer_idx = i*n_gauss_q_vals;
      const double x1{rv_func_x1_buf[i]};
      const __m128d vx1{_mm_set1_pd(x1)};
      const double rho{rv_func_rho_buf[i]};
      const __m128d vrho{_mm_set1_pd(rho)};
      for(std::int32_t j{0}; j < ROUND_DOWN(n_gauss_q_vals,1); j += 2) 
      {
          const std::int32_t inner_idx = outer_idx+j;
          const __m128d theta{_mm_load_pd(&in_buf[inner_idx])};
          const __m128d gauss_Q_res{gms::fading_channel::integrand_4_7_x1_gauss_Q_func_sse_pd(vx1,vrho,theta)};
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 1
          const double * __restrict__ p_res{reinterpret_cast<const double * __restrict__>(&gauss_Q_res)};
          printf_ret = print_double("p_res[0]",p_res[0],0);
          printf_ret = print_double("p_res[1]",p_res[1],0);
#endif 
        _mm_store_pd(&out_buf[inner_idx],gauss_Q_res);
      }
    } 
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_integrand_4_7_y1_gauss_Q_func_sse_pd();

void unit_test_integrand_4_7_y1_gauss_Q_func_sse_pd()
{
    constexpr std::int32_t n_func_args{10};
    constexpr std::int32_t n_gauss_q_vals{50};
    constexpr std::int32_t tot_elems{n_func_args*n_gauss_q_vals};
    constexpr double lo_x1{0.001};
    constexpr double hi_x1{1.0};
    constexpr double lo_y1{lo_x1};
    constexpr double hi_y1{hi_x1};
    constexpr double lo_rho{0.0};
    constexpr double hi_rho{1.0};
    __ATTR_ALIGN__(16) double in_buf[tot_elems];
    __ATTR_ALIGN__(16) double out_buf[tot_elems];
    __ATTR_ALIGN__(16) double rv_func_x1_buf[n_func_args];
    __ATTR_ALIGN__(16) double rv_func_y1_buf[n_func_args];
    __ATTR_ALIGN__(16) double rv_func_rho_buf[n_func_args];
    __ATTR_ALIGN__(16) double rv_phi_arg_buf[n_func_args];
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
    rv_func_arg_x1 = std::uniform_real_distribution<double>(lo_x1,hi_x1);
    seed_func_x1_arg = __rdtsc();
    rv_func_arg_x1_gen = std::mt19937(seed_func_x1_arg);
    rv_func_arg_y1 = std::uniform_real_distribution<double>(lo_y1,hi_y1);
    seed_func_y1_arg = __rdtsc();
    rv_func_arg_y1_gen = std::mt19937(seed_func_y1_arg);
    rv_func_arg_rho = std::uniform_real_distribution<double>(lo_rho,hi_rho);
    seed_func_rho_arg = __rdtsc();
    rv_func_arg_rho_gen = std::mt19937(seed_func_rho_arg);
    [[maybe_unused]] std::int32_t printf_ret{};
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!.\n",__func__);

    for(std::int32_t i {0}; i < n_func_args; ++i) 
    {
        const std::int32_t outer_idx = i*n_gauss_q_vals;
        const double x1 = rv_func_arg_x1.operator()(rv_func_arg_x1_gen);
        rv_func_x1_buf[i] = x1;
        const double y1 = rv_func_arg_y1.operator()(rv_func_arg_y1_gen);
        rv_func_y1_buf[i] = y1;
        const double rho= rv_func_arg_rho.operator()(rv_func_arg_rho_gen);
        rv_func_rho_buf[i] = rho;
        const double phi_arg = y1/x1;
        rv_phi_arg_buf[i] = phi_arg;
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 1
        printf_ret = print_double("phi_arg",phi_arg,0);
#endif 
        const double phi_s = std::atan(phi_arg);
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 1
        printf_ret = print_double("phi_s",phi_s,0);
#endif 
        rv_gauss_q = std::uniform_real_distribution<double>(+0.0,+1.57079632679489662-phi_s);
        seed_gauss_q = __rdtsc();
        rv_gauss_q_gen = std::mt19937(seed_gauss_q);
        for(std::int32_t j{0}; j < n_gauss_q_vals; ++j)   
        {
            const std::int32_t inner_idx = outer_idx+j;
            const double theta = rv_gauss_q.operator()(rv_gauss_q_gen);
            in_buf[inner_idx] = theta;
        }
    }
    std::sort(&in_buf[0],&in_buf[tot_elems-1]);
    std::sort(&rv_func_x1_buf[0],&rv_func_x1_buf[n_func_args-1]);
    // Not used and left for verification.
    std::sort(&rv_func_y1_buf[0],&rv_func_y1_buf[n_func_args-1]);
    std::sort(&rv_func_rho_buf[0],&rv_func_rho_buf[n_func_args-1]);
    // Not used and left for verification.
    std::sort(&rv_phi_arg_buf[0],&rv_phi_arg_buf[n_func_args-1]);

#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 0
     for(std::int32_t i {0}; i < tot_elems; ++i) 
     {  
        if(i>=0 && i<=(n_func_args-1)) 
        {
          print_double("sort-check: rv_func_x1_buf",rv_func_x1_buf[i],0);
          print_double("sort-check: rv_func_y1_buf",rv_func_y1_buf[i],0);
          print_double("sort-check: rv_func_rho_buf",rv_func_rho_buf[i],0);
          print_double("sort-check: rv_phi_arg_buf",rv_phi_arg_buf[i],0);
        }
        print_double("sort-check: in_buf=",in_buf[i],0);
     }
#endif 

    for(std::int32_t i{0}; i < n_func_args; ++i) 
    {  
      const std::int32_t outer_idx = i*n_gauss_q_vals;
      const double y1{rv_func_y1_buf[i]};
      const __m128d vy1{_mm_set1_pd(y1)};
      const double rho{rv_func_rho_buf[i]};
      const __m128d vrho{_mm_set1_pd(rho)};
      for(std::int32_t j{0}; j < ROUND_DOWN(n_gauss_q_vals,1); j += 2) 
      {
          const std::int32_t inner_idx = outer_idx+j;
          const __m128d theta{_mm_load_pd(&in_buf[inner_idx])};
          const __m128d gauss_Q_res{gms::fading_channel::integrand_4_7_y1_gauss_Q_func_sse_pd(vy1,vrho,theta)};
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 1
          const double * __restrict__ p_res{reinterpret_cast<const double * __restrict__>(&gauss_Q_res)};
          printf_ret = print_double("p_res[0]",p_res[0],0);
          printf_ret = print_double("p_res[1]",p_res[1],0);
#endif 
        _mm_store_pd(&out_buf[inner_idx],gauss_Q_res);
      }
    } 
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_integrand_4_8_x1_gauss_Q_func_sse_pd();

void unit_test_integrand_4_8_x1_gauss_Q_func_sse_pd()
{
    constexpr std::int32_t n_func_args{10};
    constexpr std::int32_t n_gauss_q_vals{50};
    constexpr std::int32_t tot_elems{n_func_args*n_gauss_q_vals};
    constexpr double lo_x1{0.1};
    constexpr double hi_x1{1.0};
    constexpr double lo_y1{lo_x1};
    constexpr double hi_y1{hi_x1};
    __ATTR_ALIGN__(16) double in_buf[tot_elems];
    __ATTR_ALIGN__(16) double out_buf[tot_elems];
    __ATTR_ALIGN__(16) double rv_func_x1_buf[n_func_args];
    __ATTR_ALIGN__(16) double rv_func_y1_buf[n_func_args];
    __ATTR_ALIGN__(16) double rv_phi_arg_buf[n_func_args];
    thread_local std::uniform_real_distribution<double> rv_func_arg_x1;
    thread_local std::mt19937 rv_func_arg_x1_gen;
    thread_local std::uint64_t seed_func_x1_arg{};
    thread_local std::uniform_real_distribution<double> rv_func_arg_y1;
    thread_local std::mt19937 rv_func_arg_y1_gen;
    thread_local std::uint64_t seed_func_y1_arg{};
    thread_local std::uniform_real_distribution<double> rv_gauss_q;
    thread_local std::mt19937 rv_gauss_q_gen;
    thread_local std::uint64_t seed_gauss_q{};
    rv_func_arg_x1 = std::uniform_real_distribution<double>(lo_x1,hi_x1);
    seed_func_x1_arg = __rdtsc();
    rv_func_arg_x1_gen = std::mt19937(seed_func_x1_arg);
    rv_func_arg_y1 = std::uniform_real_distribution<double>(lo_y1,hi_y1);
    seed_func_y1_arg = __rdtsc();
    rv_func_arg_y1_gen = std::mt19937(seed_func_y1_arg);
    [[maybe_unused]] std::int32_t printf_ret{};
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!.\n",__func__);

    for(std::int32_t i {0}; i < n_func_args; ++i) 
    {
        const std::int32_t outer_idx = i*n_gauss_q_vals;
        const double x1 = rv_func_arg_x1.operator()(rv_func_arg_x1_gen);
        rv_func_x1_buf[i] = x1;
        const double y1 = rv_func_arg_y1.operator()(rv_func_arg_y1_gen);
        rv_func_y1_buf[i] = y1;
        const double phi_arg = y1/x1;
        rv_phi_arg_buf[i] = phi_arg;
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 0
        printf_ret = print_double("phi_arg",phi_arg,0);
#endif 
        const double phi_s = std::atan(phi_arg);
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 0
        printf_ret = print_double("phi_s",phi_s,0);
#endif 
        rv_gauss_q = std::uniform_real_distribution<double>(+0.0,+1.57079632679489662-phi_s);
        seed_gauss_q = __rdtsc();
        rv_gauss_q_gen = std::mt19937(seed_gauss_q);
        for(std::int32_t j{0}; j < n_gauss_q_vals; ++j)   
        {
            const std::int32_t inner_idx = outer_idx+j;
            const double theta = rv_gauss_q.operator()(rv_gauss_q_gen);
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 0
            print_double("theta:",theta,0);
#endif
            in_buf[inner_idx] = theta;
        }
    }
    std::sort(&in_buf[0],&in_buf[tot_elems-1]);
    std::sort(&rv_func_x1_buf[0],&rv_func_x1_buf[n_func_args-1]);
    // Not used and left for verification.
    std::sort(&rv_func_y1_buf[0],&rv_func_y1_buf[n_func_args-1]);
    // Not used and left for verification.
    std::sort(&rv_phi_arg_buf[0],&rv_phi_arg_buf[n_func_args-1]);

#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 0
     for(std::int32_t i {0}; i < tot_elems; ++i) 
     {  
        if(i>=0 && i<=(n_func_args-1)) 
        {
          print_double("sort-check: rv_func_x1_buf",rv_func_x1_buf[i],0);
          print_double("sort-check: rv_func_y1_buf",rv_func_y1_buf[i],0);
          print_double("sort-check: rv_phi_arg_buf",rv_phi_arg_buf[i],0);
        }
        print_double("sort-check: in_buf=",in_buf[i],0);
     }
#endif 

    for(std::int32_t i{0}; i < n_func_args; ++i) 
    {  
      const std::int32_t outer_idx = i*n_gauss_q_vals;
      const double x1{rv_func_x1_buf[i]}; 
      const __m128d vx1{_mm_set1_pd(x1)};
      for(std::int32_t j{0}; j < ROUND_DOWN(n_gauss_q_vals,1); j += 2) 
      {
          const std::int32_t inner_idx = outer_idx+j;
          const __m128d theta{_mm_load_pd(&in_buf[inner_idx])};
          const __m128d gauss_Q_res{gms::fading_channel::integrand_4_8_x1_gauss_Q_func_sse_pd(vx1,theta)};
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 1
          const double * __restrict__ p_res{reinterpret_cast<const double * __restrict__>(&gauss_Q_res)};
          printf_ret = print_double("p_res[0]",p_res[0],0);
          printf_ret = print_double("p_res[1]",p_res[1],0);
#endif 
        _mm_store_pd(&out_buf[inner_idx],gauss_Q_res);
      }
    } 
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_integrand_4_8_y1_gauss_Q_func_sse_pd();

void unit_test_integrand_4_8_y1_gauss_Q_func_sse_pd()
{
    constexpr std::int32_t n_func_args{10};
    constexpr std::int32_t n_gauss_q_vals{50};
    constexpr std::int32_t tot_elems{n_func_args*n_gauss_q_vals};
    constexpr double lo_x1{0.1};
    constexpr double hi_x1{1.0};
    constexpr double lo_y1{lo_x1};
    constexpr double hi_y1{hi_x1};
    __ATTR_ALIGN__(16) double in_buf[tot_elems];
    __ATTR_ALIGN__(16) double out_buf[tot_elems];
    __ATTR_ALIGN__(16) double rv_func_x1_buf[n_func_args];
    __ATTR_ALIGN__(16) double rv_func_y1_buf[n_func_args];
    __ATTR_ALIGN__(16) double rv_phi_arg_buf[n_func_args];
    thread_local std::uniform_real_distribution<double> rv_func_arg_x1;
    thread_local std::mt19937 rv_func_arg_x1_gen;
    thread_local std::uint64_t seed_func_x1_arg{};
    thread_local std::uniform_real_distribution<double> rv_func_arg_y1;
    thread_local std::mt19937 rv_func_arg_y1_gen;
    thread_local std::uint64_t seed_func_y1_arg{};
    thread_local std::uniform_real_distribution<double> rv_gauss_q;
    thread_local std::mt19937 rv_gauss_q_gen;
    thread_local std::uint64_t seed_gauss_q{};
    rv_func_arg_x1 = std::uniform_real_distribution<double>(lo_x1,hi_x1);
    seed_func_x1_arg = __rdtsc();
    rv_func_arg_x1_gen = std::mt19937(seed_func_x1_arg);
    rv_func_arg_y1 = std::uniform_real_distribution<double>(lo_y1,hi_y1);
    seed_func_y1_arg = __rdtsc();
    rv_func_arg_y1_gen = std::mt19937(seed_func_y1_arg);
    [[maybe_unused]] std::int32_t printf_ret{};
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!.\n",__func__);

    for(std::int32_t i {0}; i < n_func_args; ++i) 
    {
        const std::int32_t outer_idx = i*n_gauss_q_vals;
        const double x1 = rv_func_arg_x1.operator()(rv_func_arg_x1_gen);
        rv_func_x1_buf[i] = x1;
        const double y1 = rv_func_arg_y1.operator()(rv_func_arg_y1_gen);
        rv_func_y1_buf[i] = y1;
        const double phi_arg = y1/x1;
        rv_phi_arg_buf[i] = phi_arg;
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 0
        printf_ret = print_double("phi_arg",phi_arg,0);
#endif 
        const double phi_s = std::atan(phi_arg);
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 0
        printf_ret = print_double("phi_s",phi_s,0);
#endif 
        rv_gauss_q = std::uniform_real_distribution<double>(+0.0,+1.57079632679489662-phi_s);
        seed_gauss_q = __rdtsc();
        rv_gauss_q_gen = std::mt19937(seed_gauss_q);
        for(std::int32_t j{0}; j < n_gauss_q_vals; ++j)   
        {
            const std::int32_t inner_idx = outer_idx+j;
            const double theta = rv_gauss_q.operator()(rv_gauss_q_gen);
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 0
            print_double("theta:",theta,0);
#endif
            in_buf[inner_idx] = theta;
        }
    }
    std::sort(&in_buf[0],&in_buf[tot_elems-1]);
    std::sort(&rv_func_x1_buf[0],&rv_func_x1_buf[n_func_args-1]);
    // Not used and left for verification.
    std::sort(&rv_func_y1_buf[0],&rv_func_y1_buf[n_func_args-1]);
    // Not used and left for verification.
    std::sort(&rv_phi_arg_buf[0],&rv_phi_arg_buf[n_func_args-1]);

#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 0
     for(std::int32_t i {0}; i < tot_elems; ++i) 
     {  
        if(i>=0 && i<=(n_func_args-1)) 
        {
          print_double("sort-check: rv_func_x1_buf",rv_func_x1_buf[i],0);
          print_double("sort-check: rv_func_y1_buf",rv_func_y1_buf[i],0);
          print_double("sort-check: rv_phi_arg_buf",rv_phi_arg_buf[i],0);
        }
        print_double("sort-check: in_buf=",in_buf[i],0);
     }
#endif 

    for(std::int32_t i{0}; i < n_func_args; ++i) 
    {  
      const std::int32_t outer_idx = i*n_gauss_q_vals;
      const double y1{rv_func_y1_buf[i]}; 
      const __m128d vy1{_mm_set1_pd(y1)};
      for(std::int32_t j{0}; j < ROUND_DOWN(n_gauss_q_vals,1); j += 2) 
      {
          const std::int32_t inner_idx = outer_idx+j;
          const __m128d theta{_mm_load_pd(&in_buf[inner_idx])};
          const __m128d gauss_Q_res{gms::fading_channel::integrand_4_8_y1_gauss_Q_func_sse_pd(vy1,theta)};
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 1
          const double * __restrict__ p_res{reinterpret_cast<const double * __restrict__>(&gauss_Q_res)};
          printf_ret = print_double("p_res[0]",p_res[0],0);
          printf_ret = print_double("p_res[1]",p_res[1],0);
#endif 
        _mm_store_pd(&out_buf[inner_idx],gauss_Q_res);
      }
    } 
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_integrand_4_10_marcum_Q_func_sse_pd();

void unit_test_integrand_4_10_marcum_Q_func_sse_pd()
{
    constexpr std::int32_t n_func_args{10};
    constexpr std::int32_t n_marcum_q_vals{50};
    constexpr std::int32_t tot_elems{n_func_args*n_marcum_q_vals};
    constexpr double lo_s{0.1};
    constexpr double hi_s{10.0};
    __ATTR_ALIGN__(16) double in_buf[tot_elems];
    __ATTR_ALIGN__(16) double out_buf[tot_elems];
    __ATTR_ALIGN__(16) double rv_func_s_buf[n_func_args];
    thread_local std::uniform_real_distribution<double> rv_func_arg_s;
    thread_local std::mt19937 rv_func_arg_s_gen;
    thread_local std::uint64_t seed_func_s_arg{};
    thread_local std::uniform_real_distribution<double> rv_marcum_q;
    thread_local std::mt19937 rv_marcum_q_gen;
    thread_local std::uint64_t seed_marcum_q{};
    rv_func_arg_s = std::uniform_real_distribution<double>(lo_s,hi_s);
    seed_func_s_arg = __rdtsc();
    rv_func_arg_s_gen = std::mt19937(seed_func_s_arg);
    [[maybe_unused]] std::int32_t printf_ret{};
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!.\n",__func__);

    for(std::int32_t i {0}; i < n_func_args; ++i) 
    {
        const std::int32_t outer_idx = i*n_marcum_q_vals;
        const double s = rv_func_arg_s.operator()(rv_func_arg_s_gen);
        rv_func_s_buf[i] = s;
        rv_marcum_q = std::uniform_real_distribution<double>(+0.0,+5.0);
        seed_marcum_q = __rdtsc();
        rv_marcum_q_gen = std::mt19937(seed_marcum_q);
        for(std::int32_t j{0}; j < n_marcum_q_vals; ++j)   
        {
            const std::int32_t inner_idx = outer_idx+j;
            const double x = rv_marcum_q.operator()(rv_marcum_q_gen);
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 0
            print_double("x:",x,0);
#endif
            in_buf[inner_idx] = x;
        }
    }
    std::sort(&in_buf[0],&in_buf[tot_elems-1]);
    std::sort(&rv_func_s_buf[0],&rv_func_s_buf[n_func_args-1]);
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 0
     for(std::int32_t i {0}; i < tot_elems; ++i) 
     {  
        if(i>=0 && i<=(n_func_args-1)) 
        {
          print_double("sort-check: rv_func_s_buf",rv_func_s_buf[i],0);
        }
        print_double("sort-check: in_buf=",in_buf[i],0);
     }
#endif 

    for(std::int32_t i{0}; i < n_func_args; ++i) 
    {  
      const std::int32_t outer_idx = i*n_marcum_q_vals;
      const double s{rv_func_s_buf[i]}; 
      const __m128d vs{_mm_set1_pd(s)};
      for(std::int32_t j{0}; j < ROUND_DOWN(n_marcum_q_vals,1); j += 2) 
      {
          const std::int32_t inner_idx = outer_idx+j;
          const __m128d vx{_mm_load_pd(&in_buf[inner_idx])};
          const __m128d marcum_Q_res{gms::fading_channel::integrand_4_10_marcum_Q_func_sse_pd(vx,vs)};
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 1
          const double * __restrict__ p_res{reinterpret_cast<const double * __restrict__>(&marcum_Q_res)};
          printf_ret = print_double("p_res[0]",p_res[0],0);
          printf_ret = print_double("p_res[1]",p_res[1],0);
#endif 
        _mm_store_pd(&out_buf[inner_idx],marcum_Q_res);
      }
    } 
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_integrand_4_16_marcum_Q_func_sse_pd();

void unit_test_integrand_4_16_marcum_Q_func_sse_pd()
{
    constexpr std::int32_t n_func_args{10};
    constexpr std::int32_t n_marcum_q_vals{50};
    constexpr std::int32_t tot_elems{n_func_args*n_marcum_q_vals};
    constexpr double lo_theta{-3.141592653589793238462643383};
    constexpr double hi_theta{+3.141592653589793238462643383};
    __ATTR_ALIGN__(16) double in_buf[tot_elems];
    __ATTR_ALIGN__(16) double out_buf[tot_elems];
    __ATTR_ALIGN__(16) double rv_func_psi_buf[n_func_args];
    __ATTR_ALIGN__(16) double rv_func_beta_buf[n_func_args];
    thread_local std::uniform_real_distribution<double> rv_func_arg_theta;
    thread_local std::mt19937 rv_func_arg_theta_gen;
    thread_local std::uint64_t seed_func_theta_arg{};
    thread_local std::uniform_real_distribution<double> rv_func_arg_psi;
    thread_local std::mt19937 rv_func_arg_psi_gen;
    thread_local std::uint64_t seed_func_psi_arg{};
    thread_local std::uniform_real_distribution<double> rv_func_arg_beta;
    thread_local std::mt19937 rv_func_arg_beta_gen;
    thread_local std::uint64_t seed_func_beta_arg{};
    rv_func_arg_theta = std::uniform_real_distribution<double>(lo_theta,hi_theta);
    seed_func_theta_arg = __rdtsc();
    rv_func_arg_theta_gen = std::mt19937(seed_func_theta_arg);
    rv_func_arg_psi = std::uniform_real_distribution<double>(0.0,1.0);
    seed_func_psi_arg = __rdtsc();
    rv_func_arg_psi_gen = std::mt19937(seed_func_psi_arg);
    rv_func_arg_beta    = std::uniform_real_distribution<double>(0.1,1.0);
    seed_func_beta_arg  = __rdtsc();
    rv_func_arg_beta_gen     = std::mt19937(seed_func_beta_arg);
    [[maybe_unused]] std::int32_t printf_ret{};
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!.\n",__func__);

    for(std::int32_t i {0}; i < n_func_args; ++i) 
    {
        const std::int32_t outer_idx = i*n_marcum_q_vals;
        const double psi = rv_func_arg_psi.operator()(rv_func_arg_psi_gen);
        rv_func_psi_buf[i] = psi;
        const double beta = rv_func_arg_beta.operator()(rv_func_arg_beta_gen);
        rv_func_beta_buf[i] = beta;
        for(std::int32_t j{0}; j < n_marcum_q_vals; ++j)   
        {
            const std::int32_t inner_idx = outer_idx+j;
            const double theta = rv_func_arg_theta.operator()(rv_func_arg_theta_gen);
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 0
            print_double("theta:",theta,0);
#endif
            in_buf[inner_idx] = theta;
        }
    }
    std::sort(&in_buf[0],&in_buf[tot_elems-1]);
    std::sort(&rv_func_psi_buf[0],&rv_func_psi_buf[n_func_args-1]);
    std::sort(&rv_func_beta_buf[0],&rv_func_beta_buf[n_func_args-1]);
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 0
     for(std::int32_t i {0}; i < tot_elems; ++i) 
     {  
        if(i>=0 && i<=(n_func_args-1)) 
        {
          print_double("sort-check: rv_func_beta_buf",rv_func_beta_buf[i],0);
          print_double("sort-check: rv_func_psi_buf",rv_func_psi_buf[i],0);
         
        }
        print_double("sort-check: in_buf=",in_buf[i],0);
     }
#endif 

    for(std::int32_t i{0}; i < n_func_args; ++i) 
    {  
      const std::int32_t outer_idx = i*n_marcum_q_vals;
      const double beta{rv_func_beta_buf[i]}; 
      const __m128d vbeta{_mm_set1_pd(beta)};
      const double psi{rv_func_psi_buf[i]};
      const __m128d vpsi{_mm_set1_pd(psi)};
      for(std::int32_t j{0}; j < ROUND_DOWN(n_marcum_q_vals,1); j += 2) 
      {
          const std::int32_t inner_idx = outer_idx+j;
          const __m128d theta{_mm_load_pd(&in_buf[inner_idx])};
          const __m128d marcum_Q_res{gms::fading_channel::integrand_4_16_marcum_Q_func_sse_pd(vbeta,vpsi,theta)};
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 1
          const double * __restrict__ p_res{reinterpret_cast<const double * __restrict__>(&marcum_Q_res)};
          printf_ret = print_double("p_res[0]",p_res[0],0);
          printf_ret = print_double("p_res[1]",p_res[1],0);
#endif 
        _mm_store_pd(&out_buf[inner_idx],marcum_Q_res);
      }
    } 
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_integrand_4_20_marcum_Q_func_lo_sse_pd();

void unit_test_integrand_4_20_marcum_Q_func_lo_sse_pd()
{
    constexpr std::int32_t n_func_args{10};
    constexpr std::int32_t n_marcum_q_vals{50};
    constexpr std::int32_t tot_elems{n_func_args*n_marcum_q_vals};
    constexpr double lo_theta{-3.141592653589793238462643383};
    constexpr double hi_theta{+0.0};
    __ATTR_ALIGN__(16) double in_buf[tot_elems];
    __ATTR_ALIGN__(16) double out_buf[tot_elems];
    __ATTR_ALIGN__(16) double rv_func_psi_buf[n_func_args];
    __ATTR_ALIGN__(16) double rv_func_beta_buf[n_func_args];
    thread_local std::uniform_real_distribution<double> rv_func_arg_theta;
    thread_local std::mt19937 rv_func_arg_theta_gen;
    thread_local std::uint64_t seed_func_theta_arg{};
    thread_local std::uniform_real_distribution<double> rv_func_arg_psi;
    thread_local std::mt19937 rv_func_arg_psi_gen;
    thread_local std::uint64_t seed_func_psi_arg{};
    thread_local std::uniform_real_distribution<double> rv_func_arg_beta;
    thread_local std::mt19937 rv_func_arg_beta_gen;
    thread_local std::uint64_t seed_func_beta_arg{};
    rv_func_arg_theta = std::uniform_real_distribution<double>(lo_theta,hi_theta);
    seed_func_theta_arg = __rdtsc();
    rv_func_arg_theta_gen = std::mt19937(seed_func_theta_arg);
    rv_func_arg_psi = std::uniform_real_distribution<double>(0.0,1.0);
    seed_func_psi_arg = __rdtsc();
    rv_func_arg_psi_gen = std::mt19937(seed_func_psi_arg);
    rv_func_arg_beta    = std::uniform_real_distribution<double>(0.1,1.0);
    seed_func_beta_arg  = __rdtsc();
    rv_func_arg_beta_gen     = std::mt19937(seed_func_beta_arg);
    [[maybe_unused]] std::int32_t printf_ret{};
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!.\n",__func__);

    for(std::int32_t i {0}; i < n_func_args; ++i) 
    {
        const std::int32_t outer_idx = i*n_marcum_q_vals;
        const double psi = rv_func_arg_psi.operator()(rv_func_arg_psi_gen);
        rv_func_psi_buf[i] = psi;
        const double beta = rv_func_arg_beta.operator()(rv_func_arg_beta_gen);
        rv_func_beta_buf[i] = beta;
        for(std::int32_t j{0}; j < n_marcum_q_vals; ++j)   
        {
            const std::int32_t inner_idx = outer_idx+j;
            const double theta = rv_func_arg_theta.operator()(rv_func_arg_theta_gen);
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 0
            print_double("theta:",theta,0);
#endif
            in_buf[inner_idx] = theta;
        }
    }
    std::sort(&in_buf[0],&in_buf[tot_elems-1]);
    std::sort(&rv_func_psi_buf[0],&rv_func_psi_buf[n_func_args-1]);
    std::sort(&rv_func_beta_buf[0],&rv_func_beta_buf[n_func_args-1]);
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 0
     for(std::int32_t i {0}; i < tot_elems; ++i) 
     {  
        if(i>=0 && i<=(n_func_args-1)) 
        {
          print_double("sort-check: rv_func_beta_buf",rv_func_beta_buf[i],0);
          print_double("sort-check: rv_func_psi_buf",rv_func_psi_buf[i],0);
         
        }
        print_double("sort-check: in_buf=",in_buf[i],0);
     }
#endif 

    for(std::int32_t i{0}; i < n_func_args; ++i) 
    {  
      const std::int32_t outer_idx = i*n_marcum_q_vals;
      const double beta{rv_func_beta_buf[i]}; 
      const __m128d vbeta{_mm_set1_pd(beta)};
      const double psi{rv_func_psi_buf[i]};
      const __m128d vpsi{_mm_set1_pd(psi)};
      for(std::int32_t j{0}; j < ROUND_DOWN(n_marcum_q_vals,1); j += 2) 
      {
          const std::int32_t inner_idx = outer_idx+j;
          const __m128d theta{_mm_load_pd(&in_buf[inner_idx])};
          const __m128d marcum_Q_res{gms::fading_channel::integrand_4_20_marcum_Q_func_lo_sse_pd(vbeta,vpsi,theta)};
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 1
          const double * __restrict__ p_res{reinterpret_cast<const double * __restrict__>(&marcum_Q_res)};
          printf_ret = print_double("p_res[0]",p_res[0],0);
          printf_ret = print_double("p_res[1]",p_res[1],0);
#endif 
        _mm_store_pd(&out_buf[inner_idx],marcum_Q_res);
      }
    } 
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_integrand_4_20_marcum_Q_func_hi_sse_pd();

void unit_test_integrand_4_20_marcum_Q_func_hi_sse_pd()
{
    constexpr std::int32_t n_func_args{10};
    constexpr std::int32_t n_marcum_q_vals{50};
    constexpr std::int32_t tot_elems{n_func_args*n_marcum_q_vals};
    constexpr double lo_theta{+0.0};
    constexpr double hi_theta{+3.141592653589793238462643383};
    __ATTR_ALIGN__(16) double in_buf[tot_elems];
    __ATTR_ALIGN__(16) double out_buf[tot_elems];
    __ATTR_ALIGN__(16) double rv_func_psi_buf[n_func_args];
    __ATTR_ALIGN__(16) double rv_func_beta_buf[n_func_args];
    thread_local std::uniform_real_distribution<double> rv_func_arg_theta;
    thread_local std::mt19937 rv_func_arg_theta_gen;
    thread_local std::uint64_t seed_func_theta_arg{};
    thread_local std::uniform_real_distribution<double> rv_func_arg_psi;
    thread_local std::mt19937 rv_func_arg_psi_gen;
    thread_local std::uint64_t seed_func_psi_arg{};
    thread_local std::uniform_real_distribution<double> rv_func_arg_beta;
    thread_local std::mt19937 rv_func_arg_beta_gen;
    thread_local std::uint64_t seed_func_beta_arg{};
    rv_func_arg_theta = std::uniform_real_distribution<double>(lo_theta,hi_theta);
    seed_func_theta_arg = __rdtsc();
    rv_func_arg_theta_gen = std::mt19937(seed_func_theta_arg);
    rv_func_arg_psi = std::uniform_real_distribution<double>(0.0,1.0);
    seed_func_psi_arg = __rdtsc();
    rv_func_arg_psi_gen = std::mt19937(seed_func_psi_arg);
    rv_func_arg_beta    = std::uniform_real_distribution<double>(0.1,1.0);
    seed_func_beta_arg  = __rdtsc();
    rv_func_arg_beta_gen     = std::mt19937(seed_func_beta_arg);
    [[maybe_unused]] std::int32_t printf_ret{};
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!.\n",__func__);

    for(std::int32_t i {0}; i < n_func_args; ++i) 
    {
        const std::int32_t outer_idx = i*n_marcum_q_vals;
        const double psi = rv_func_arg_psi.operator()(rv_func_arg_psi_gen);
        rv_func_psi_buf[i] = psi;
        const double beta = rv_func_arg_beta.operator()(rv_func_arg_beta_gen);
        rv_func_beta_buf[i] = beta;
        for(std::int32_t j{0}; j < n_marcum_q_vals; ++j)   
        {
            const std::int32_t inner_idx = outer_idx+j;
            const double theta = rv_func_arg_theta.operator()(rv_func_arg_theta_gen);
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 0
            print_double("theta:",theta,0);
#endif
            in_buf[inner_idx] = theta;
        }
    }
    std::sort(&in_buf[0],&in_buf[tot_elems-1]);
    std::sort(&rv_func_psi_buf[0],&rv_func_psi_buf[n_func_args-1]);
    std::sort(&rv_func_beta_buf[0],&rv_func_beta_buf[n_func_args-1]);
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 0
     for(std::int32_t i {0}; i < tot_elems; ++i) 
     {  
        if(i>=0 && i<=(n_func_args-1)) 
        {
          print_double("sort-check: rv_func_beta_buf",rv_func_beta_buf[i],0);
          print_double("sort-check: rv_func_psi_buf",rv_func_psi_buf[i],0);
         
        }
        print_double("sort-check: in_buf=",in_buf[i],0);
     }
#endif 

    for(std::int32_t i{0}; i < n_func_args; ++i) 
    {  
      const std::int32_t outer_idx = i*n_marcum_q_vals;
      const double beta{rv_func_beta_buf[i]}; 
      const __m128d vbeta{_mm_set1_pd(beta)};
      const double psi{rv_func_psi_buf[i]};
      const __m128d vpsi{_mm_set1_pd(psi)};
      for(std::int32_t j{0}; j < ROUND_DOWN(n_marcum_q_vals,1); j += 2) 
      {
          const std::int32_t inner_idx = outer_idx+j;
          const __m128d theta{_mm_load_pd(&in_buf[inner_idx])};
          const __m128d marcum_Q_res{gms::fading_channel::integrand_4_20_marcum_Q_func_hi_sse_pd(vbeta,vpsi,theta)};
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 1
          const double * __restrict__ p_res{reinterpret_cast<const double * __restrict__>(&marcum_Q_res)};
          printf_ret = print_double("p_res[0]",p_res[0],0);
          printf_ret = print_double("p_res[1]",p_res[1],0);
#endif 
        _mm_store_pd(&out_buf[inner_idx],marcum_Q_res);
      }
    } 
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}

__attribute__((hot))
__attribute__((aligned(32)))
void unit_test_integrand_4_32_marcum_Q_m_func_sse_pd();

void unit_test_integrand_4_32_marcum_Q_m_func_sse_pd()
{
    constexpr std::int32_t n_func_args{10};
    constexpr std::int32_t n_marcum_q_vals{50};
    constexpr std::int32_t tot_elems{n_func_args*n_marcum_q_vals};
    constexpr double lo_s{+0.1};
    constexpr double hi_s{+10.0};
    __ATTR_ALIGN__(16) double in_buf[tot_elems];
    __ATTR_ALIGN__(16) double out_buf[tot_elems];
    __ATTR_ALIGN__(16) double rv_func_s_buf[n_func_args];
    thread_local std::uniform_real_distribution<double> rv_s_param;
    thread_local std::mt19937 rv_s_param_gen;
    thread_local std::uint64_t seed_s_param{};
    thread_local std::uniform_real_distribution<double> rv_marcum_q;
    thread_local std::mt19937 rv_marcum_q_gen;
    thread_local std::uint64_t seed_marcum_q{};
    rv_s_param = std::uniform_real_distribution<double>(lo_s,hi_s);
    seed_s_param = __rdtsc();
    rv_s_param_gen = std::mt19937(seed_s_param);
    rv_marcum_q = std::uniform_real_distribution<double>(0.0,5.0);
    seed_marcum_q = __rdtsc();
    rv_marcum_q_gen = std::mt19937(seed_marcum_q);
    [[maybe_unused]] std::int32_t printf_ret{};
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED!!.\n",__func__);

    for(std::int32_t i {0}; i < n_func_args; ++i) 
    {
        const std::int32_t outer_idx = i*n_marcum_q_vals;
        const double s = rv_s_param.operator()(rv_s_param_gen);
        rv_func_s_buf[i] = s;
        for(std::int32_t j{0}; j < n_marcum_q_vals; ++j)   
        {
            const std::int32_t inner_idx = outer_idx+j;
            const double x = rv_marcum_q.operator()(rv_marcum_q_gen);
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 0
            print_double("x:",x,0);
#endif
            in_buf[inner_idx] = x;
        }
    }
    std::sort(&in_buf[0],&in_buf[tot_elems-1]);
    std::sort(&rv_func_s_buf[0],&rv_func_s_buf[n_func_args-1]);
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 0
     for(std::int32_t i {0}; i < tot_elems; ++i) 
     {  
        if(i>=0 && i<=(n_func_args-1)) 
        {
          print_double("sort-check: rv_func_s_buf",rv_func_s_buf[i],0);       
        }
        print_double("sort-check: in_buf=",in_buf[i],0);
     }
#endif 
    const __m128d d_m{_mm_set1_pd(2.0)};
    for(std::int32_t i{0}; i < n_func_args; ++i) 
    {  
      const std::int32_t outer_idx = i*n_marcum_q_vals;
      const double s{rv_func_s_buf[i]}; 
      const __m128d vs{_mm_set1_pd(s)};
      for(std::int32_t j{0}; j < ROUND_DOWN(n_marcum_q_vals,1); j += 2) 
      {
          const std::int32_t inner_idx = outer_idx+j;
          const __m128d vx{_mm_load_pd(&in_buf[inner_idx])};
          const __m128d marcum_Q_res{gms::fading_channel::integrand_4_32_marcum_Q_m_func_sse_pd(vx,vs,d_m,1)};
#if (UNIT_TEST_INTEGRANDS_FUNC_CH4_DBG_PRINT) == 1
          const double * __restrict__ p_res{reinterpret_cast<const double * __restrict__>(&marcum_Q_res)};
          printf_ret = print_double("p_res[0]",p_res[0],0);
          printf_ret = print_double("p_res[1]",p_res[1],0);
#endif 
        _mm_store_pd(&out_buf[inner_idx],marcum_Q_res);
      }
    } 
    printf_ret = std::printf("[UNIT-TEST:] -- of function=%s -- ENDED!!\n",__func__);
}



int main()
{

   (void)unit_test_integrand_4_1_gauss_Q_func_sse_pd();
   (void)unit_test_integrand_4_2_gauss_Q_func_sse_pd();
   (void)unit_test_integrand_4_6_sin_gauss_Q_func_sse_pd();
   (void)unit_test_integrand_4_6_cos_gauss_Q_func_sse_pd();
   (void)unit_test_integrand_4_7_x1_gauss_Q_func_sse_pd();
   (void)unit_test_integrand_4_7_y1_gauss_Q_func_sse_pd();
   (void)unit_test_integrand_4_8_x1_gauss_Q_func_sse_pd();
   (void)unit_test_integrand_4_8_y1_gauss_Q_func_sse_pd();
   (void)unit_test_integrand_4_10_marcum_Q_func_sse_pd();
   (void)unit_test_integrand_4_16_marcum_Q_func_sse_pd();
   (void)unit_test_integrand_4_20_marcum_Q_func_lo_sse_pd();
   (void)unit_test_integrand_4_20_marcum_Q_func_hi_sse_pd();
   (void)unit_test_integrand_4_32_marcum_Q_m_func_sse_pd();
    return 0;
}