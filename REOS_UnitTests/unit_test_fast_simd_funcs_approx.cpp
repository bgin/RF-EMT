#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <cmath>
#include "GMS_fast_simd_funcs_approx.h"

/*
   icpc -o unit_test_fast_simd_funcs_approx -O3 -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_fast_simd_funcs_approx.h unit_test_fast_simd_funcs_approx.cpp
   ASM: 
   icpc -S -O3 -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_fast_simd_funcs_approx.h unit_test_fast_simd_funcs_approx.cpp

*/

static inline 
bool approximatelyEqual(const float a,
		                const float b,
				        const float epsilon) 
{
	const float fabsa = std::fabs(a);
	const float fabsb = std::fabs(b);
     return std::fabs(a - b) <=
            ((fabsa < fabsb ? fabsb : fabsa) * epsilon);
}

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
void 
printNumber(const char * val_name,float f, const std::int32_t offset) 
{
  (*((int *)&f)) += offset;
  std::printf("value: %s,%+1.11g,0x%08X,%d\n",val_name,f, *(int *)&f, *(int *)&f);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_simd_fast_exp_approx_4xf32_trivial();

void unit_test_simd_fast_exp_approx_4xf32_trivial()
{
    using namespace gms::math;
    thread_local std::uniform_real_distribution<float> exp_rand_in;
    thread_local std::uint64_t seed{};
    constexpr float lo{-4.999999f};
    constexpr float hi{2.718281828459045235360287471352f};
    seed = __rdtsc();
    exp_rand_in = std::uniform_real_distribution<float>(lo,hi);
    auto rand_uni_gen{std::mt19937(seed)};
    const float in_ref_1{exp_rand_in.operator()(rand_uni_gen)};
    const float in_ref_2{exp_rand_in.operator()(rand_uni_gen)};
    const float in_ref_3{exp_rand_in.operator()(rand_uni_gen)};
    const float in_ref_4{exp_rand_in.operator()(rand_uni_gen)};
    const __m128    in_vref_vec{_mm_setr_ps(in_ref_1,in_ref_2,in_ref_3,in_ref_4)};
    __m128 vec_actual_vref{_mm_setzero_ps()};
    const float * __restrict__ p_vref_in = reinterpret_cast<const float*>(&in_vref_vec);
    float       * __restrict__ p_vref_out;
    float sca_expected_ref_1{};
    float sca_expected_ref_2{};
    float sca_expected_ref_3{};
    float sca_expected_ref_4{};
    printf("[UNIT_TEST]: function=%s -- **START**\n\n", __PRETTY_FUNCTION__);
    printf("[UNIT-TEST]: -- Calling reference implementation: expapprox.\n");
    printf("[UNIT-TEST]: Inputs are: ref_1=%.7f,ref_2=%.7f,ref_3=%.7f,ref_4=%.7f\n",in_ref_1,in_ref_2,in_ref_3,in_ref_4);
    sca_expected_ref_1 = expapprox(in_ref_1);
    sca_expected_ref_2 = expapprox(in_ref_2);
    sca_expected_ref_3 = expapprox(in_ref_3);
    sca_expected_ref_4 = expapprox(in_ref_4);
    printf("[UNIT-TEST]: -- Calling SSE implementation: simd_fast_exp_approx_4xf32.\n");
    printf("[UNIT-TEST]: Inputs are: ref_1=%.7f,ref_2=%.7f,ref_3=%.7f,ref_4=%.7f\n",p_vref_in[0],p_vref_in[1],p_vref_in[2],p_vref_in[3]);
    //__asm__("int3");
    vec_actual_vref = simd_fast_exp_approx_4xf32(in_vref_vec);
    p_vref_out      = reinterpret_cast<float*>(&vec_actual_vref);
    //printf("[UNIT-TEST]: Calculating the results approximately equal by the epsilon:%.17f\n",std::numeric_limits<double>::epsilon());
    const bool ref_1_approx_equality{almostEqual2sComplement(sca_expected_ref_1,p_vref_out[0],7)};
    const bool ref_2_approx_equality{almostEqual2sComplement(sca_expected_ref_2,p_vref_out[1],7)};
    const bool ref_3_approx_equality{almostEqual2sComplement(sca_expected_ref_3,p_vref_out[2],7)};
    const bool ref_4_approx_equality{almostEqual2sComplement(sca_expected_ref_4,p_vref_out[3],7)};
    if(!ref_1_approx_equality) 
    {
       std::printf("[UNIT-TEST #1]-- **FAIL**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_1,p_vref_out[0],std::abs(sca_expected_ref_1-p_vref_out[0]));
       printNumber("expapprox",sca_expected_ref_1,0);
       printNumber("simd_fast_exp_approx_4xf32",p_vref_out[0],0);
    }
    else 
    {
       std::printf("[UNIT-TEST #1]-- **SUCCESS**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_1,p_vref_out[0],std::abs(sca_expected_ref_1-p_vref_out[0]));
       printNumber("expapprox",sca_expected_ref_1,0);
       printNumber("simd_fast_exp_approx_4xf32",p_vref_out[0],0);
    }
    if(!ref_2_approx_equality)
    {
       std::printf("[UNIT-TEST #2]-- **FAIL**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_2,p_vref_out[1],std::abs(sca_expected_ref_2-p_vref_out[1]));
       printNumber("expapprox",sca_expected_ref_2,0);
       printNumber("simd_fast_exp_approx_4xf32",p_vref_out[1],0);
    }
    else 
    {
       std::printf("[UNIT-TEST #2]-- **SUCCESS**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_2,p_vref_out[1],std::abs(sca_expected_ref_2-p_vref_out[1]));
       printNumber("expapprox",sca_expected_ref_2,0);
       printNumber("simd_fast_exp_approx_4xf32",p_vref_out[1],0);
    }
    if(!ref_3_approx_equality) 
    {
       std::printf("[UNIT-TEST #3]-- **FAIL**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_3,p_vref_out[2],std::abs(sca_expected_ref_3-p_vref_out[2]));
       printNumber("expapprox",sca_expected_ref_3,0);
       printNumber("simd_fast_exp_approx_4xf32",p_vref_out[2],0);
    }
    else 
    {
       std::printf("[UNIT-TEST #3]-- **SUCCESS**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_3,p_vref_out[2],std::abs(sca_expected_ref_3-p_vref_out[2]));
       printNumber("expapprox",sca_expected_ref_3,0);
       printNumber("simd_fast_exp_approx_4xf32",p_vref_out[2],0);
    }
    if(!ref_4_approx_equality) 
    {
       std::printf("[UNIT-TEST #4]-- **FAIL**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_4,p_vref_out[3],std::abs(sca_expected_ref_4-p_vref_out[3]));
       printNumber("expapprox",sca_expected_ref_4,0);
       printNumber("simd_fast_exp_approx_4xf32",p_vref_out[3],0);
    }
    else 
    {
       std::printf("[UNIT-TEST #4]-- **SUCCESS**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_4,p_vref_out[3],std::abs(sca_expected_ref_4-p_vref_out[3]));
       printNumber("expapprox",sca_expected_ref_4,0);
       printNumber("simd_fast_exp_approx_4xf32",p_vref_out[3],0);
    }
     printf("[UNIT_TEST]: function=%s -- **END**\n\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_simd_fast_exp_approx_2xf64_trivial();

void unit_test_simd_fast_exp_approx_2xf64_trivial()
{
    using namespace gms::math;
    thread_local std::uniform_real_distribution<double> exp_rand_in;
    thread_local std::uint64_t seed{};
    constexpr double lo{-4.999999999999999};
    constexpr double hi{2.718281828459045235360287471352};
    seed = __rdtsc();
    exp_rand_in = std::uniform_real_distribution<double>(lo,hi);
    auto rand_uni_gen{std::mt19937(seed)};
    const double in_ref_1{exp_rand_in.operator()(rand_uni_gen)};
    const double in_ref_2{exp_rand_in.operator()(rand_uni_gen)};
    const __m128d    in_vref_vec{_mm_setr_pd(in_ref_1,in_ref_2)};
    __m128d vec_actual_vref{_mm_setzero_pd()};
    const double * __restrict__ p_vref_in = reinterpret_cast<const double*>(&in_vref_vec);
    double       * __restrict__ p_vref_out;
    double sca_expected_ref_1{};
    double sca_expected_ref_2{};
    printf("[UNIT_TEST]: function=%s -- **START**\n\n", __PRETTY_FUNCTION__);
    printf("[UNIT-TEST]: -- Calling reference implementation: expapprox.\n");
    printf("[UNIT-TEST]: Inputs are: ref_1=%.17f,ref_2=%.17f\n",in_ref_1,in_ref_2);
    sca_expected_ref_1 = expapprox_d(in_ref_1);
    sca_expected_ref_2 = expapprox_d(in_ref_2);    
    printf("[UNIT-TEST]: -- Calling SSE implementation: simd_fast_exp_approx_42xf64.\n");
    printf("[UNIT-TEST]: Inputs are: ref_1=%.17f,ref_2=%.17f\n",p_vref_in[0],p_vref_in[1]);
    //__asm__("int3");
    vec_actual_vref = simd_fast_exp_approx_2xf64(in_vref_vec);
    p_vref_out      = reinterpret_cast<double*>(&vec_actual_vref);
    //printf("[UNIT-TEST]: Calculating the results approximately equal by the epsilon:%.17f\n",std::numeric_limits<double>::epsilon());
    const bool ref_1_approx_equality{almostEqual2sComplement(sca_expected_ref_1,p_vref_out[0],17)};
    const bool ref_2_approx_equality{almostEqual2sComplement(sca_expected_ref_2,p_vref_out[1],17)};
    if(!ref_1_approx_equality) 
    {
       std::printf("[UNIT-TEST #1]-- **FAIL**, sca=%.17f,vec=%.17f,delta=%.17f\n",sca_expected_ref_1,p_vref_out[0],std::abs(sca_expected_ref_1-p_vref_out[0]));
       printNumber("expapprox",sca_expected_ref_1,0);
       printNumber("simd_fast_exp_approx_2xf64",p_vref_out[0],0);
    }
    else 
    {
       std::printf("[UNIT-TEST #1]-- **SUCCESS**, sca=%.17f,vec=%.17f,delta=%.17f\n",sca_expected_ref_1,p_vref_out[0],std::abs(sca_expected_ref_1-p_vref_out[0]));
       printNumber("expapprox",sca_expected_ref_1,0);
       printNumber("simd_fast_exp_approx_2xf64",p_vref_out[0],0);
    }
    if(!ref_2_approx_equality)
    {
       std::printf("[UNIT-TEST #2]-- **FAIL**, sca=%.17f,vec=%.17f,delta=%.17f\n",sca_expected_ref_2,p_vref_out[1],std::abs(sca_expected_ref_2-p_vref_out[1]));
       printNumber("expapprox",sca_expected_ref_2,0);
       printNumber("simd_fast_exp_approx_2xf64",p_vref_out[1],0);
    }
    else 
    {
       std::printf("[UNIT-TEST #2]-- **SUCCESS**, sca=%.17f,vec=%.17f,delta=%.17f\n",sca_expected_ref_2,p_vref_out[1],std::abs(sca_expected_ref_2-p_vref_out[1]));
       printNumber("expapprox",sca_expected_ref_2,0);
       printNumber("simd_fast_exp_approx_2xf64",p_vref_out[1],0);
    }
     printf("[UNIT_TEST]: function=%s -- **END**\n\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_simd_fast_cos_approx_4xf32_trivial();

void unit_test_simd_fast_cos_approx_4xf32_trivial()
{
    using namespace gms::math;
    thread_local std::uniform_real_distribution<float> cos_rand_in;
    thread_local std::uint64_t seed{};
    constexpr float lo{-3.14159265358979324f};
    constexpr float hi{+3.14159265358979324f};
    seed = __rdtsc();
    cos_rand_in = std::uniform_real_distribution<float>(lo,hi);
    auto rand_uni_gen{std::mt19937(seed)};
    const float in_ref_1{cos_rand_in.operator()(rand_uni_gen)};
    const float in_ref_2{cos_rand_in.operator()(rand_uni_gen)};
    const float in_ref_3{cos_rand_in.operator()(rand_uni_gen)};
    const float in_ref_4{cos_rand_in.operator()(rand_uni_gen)};
    const __m128    in_vref_vec{_mm_setr_ps(in_ref_1,in_ref_2,in_ref_3,in_ref_4)};
    __m128 vec_actual_vref{_mm_setzero_ps()};
    const float * __restrict__ p_vref_in = reinterpret_cast<const float*>(&in_vref_vec);
    float       * __restrict__ p_vref_out;
    float sca_expected_ref_1{};
    float sca_expected_ref_2{};
    float sca_expected_ref_3{};
    float sca_expected_ref_4{};
    printf("[UNIT_TEST]: function=%s -- **START**\n\n", __PRETTY_FUNCTION__);
    printf("[UNIT-TEST]: -- Calling reference implementation: cosapprox.\n");
    printf("[UNIT-TEST]: Inputs are: ref_1=%.7f,ref_2=%.7f,ref_3=%.7f,ref_4=%.7f\n",in_ref_1,in_ref_2,in_ref_3,in_ref_4);
    sca_expected_ref_1 = cosapprox(in_ref_1);
    sca_expected_ref_2 = cosapprox(in_ref_2);
    sca_expected_ref_3 = cosapprox(in_ref_3);
    sca_expected_ref_4 = cosapprox(in_ref_4);
    printf("[UNIT-TEST]: -- Calling SSE implementation: simd_fast_cos_approx_4xf32.\n");
    printf("[UNIT-TEST]: Inputs are: ref_1=%.7f,ref_2=%.7f,ref_3=%.7f,ref_4=%.7f\n",p_vref_in[0],p_vref_in[1],p_vref_in[2],p_vref_in[3]);
    //__asm__("int3");
    vec_actual_vref = simd_fast_cos_approx_4xf32(in_vref_vec);
    p_vref_out      = reinterpret_cast<float*>(&vec_actual_vref);
    //printf("[UNIT-TEST]: Calculating the results approximately equal by the epsilon:%.17f\n",std::numeric_limits<double>::epsilon());
    const bool ref_1_approx_equality{almostEqual2sComplement(sca_expected_ref_1,p_vref_out[0],7)};
    const bool ref_2_approx_equality{almostEqual2sComplement(sca_expected_ref_2,p_vref_out[1],7)};
    const bool ref_3_approx_equality{almostEqual2sComplement(sca_expected_ref_3,p_vref_out[2],7)};
    const bool ref_4_approx_equality{almostEqual2sComplement(sca_expected_ref_4,p_vref_out[3],7)};
    if(!ref_1_approx_equality) 
    {
       std::printf("[UNIT-TEST #1]-- **FAIL**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_1,p_vref_out[0],std::abs(sca_expected_ref_1-p_vref_out[0]));
       printNumber("cosapprox",sca_expected_ref_1,0);
       printNumber("simd_fast_cos_approx_4xf32",p_vref_out[0],0);
    }
    else 
    {
       std::printf("[UNIT-TEST #1]-- **SUCCESS**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_1,p_vref_out[0],std::abs(sca_expected_ref_1-p_vref_out[0]));
       printNumber("cosapprox",sca_expected_ref_1,0);
       printNumber("simd_fast_cos_approx_4xf32",p_vref_out[0],0);
    }
    if(!ref_2_approx_equality)
    {
       std::printf("[UNIT-TEST #2]-- **FAIL**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_2,p_vref_out[1],std::abs(sca_expected_ref_2-p_vref_out[1]));
       printNumber("cosapprox",sca_expected_ref_2,0);
       printNumber("simd_fast_cos_approx_4xf32",p_vref_out[1],0);
    }
    else 
    {
       std::printf("[UNIT-TEST #2]-- **SUCCESS**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_2,p_vref_out[1],std::abs(sca_expected_ref_2-p_vref_out[1]));
       printNumber("cosapprox",sca_expected_ref_2,0);
       printNumber("simd_fast_cos_approx_4xf32",p_vref_out[1],0);
    }
    if(!ref_3_approx_equality) 
    {
       std::printf("[UNIT-TEST #3]-- **FAIL**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_3,p_vref_out[2],std::abs(sca_expected_ref_3-p_vref_out[2]));
       printNumber("cosapprox",sca_expected_ref_3,0);
       printNumber("simd_fast_cos_approx_4xf32",p_vref_out[2],0);
    }
    else 
    {
       std::printf("[UNIT-TEST #3]-- **SUCCESS**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_3,p_vref_out[2],std::abs(sca_expected_ref_3-p_vref_out[2]));
       printNumber("cosapprox",sca_expected_ref_3,0);
       printNumber("simd_fast_cos_approx_4xf32",p_vref_out[2],0);
    }
    if(!ref_4_approx_equality) 
    {
       std::printf("[UNIT-TEST #4]-- **FAIL**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_4,p_vref_out[3],std::abs(sca_expected_ref_4-p_vref_out[3]));
       printNumber("cosapprox",sca_expected_ref_4,0);
       printNumber("simd_fast_cos_approx_4xf32",p_vref_out[3],0);
    }
    else 
    {
       std::printf("[UNIT-TEST #4]-- **SUCCESS**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_4,p_vref_out[3],std::abs(sca_expected_ref_4-p_vref_out[3]));
       printNumber("cosapprox",sca_expected_ref_4,0);
       printNumber("simd_fast_cos_approx_4xf32",p_vref_out[3],0);
    }
     printf("[UNIT_TEST]: function=%s -- **END**\n\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_simd_fast_cos_approx_2xf64_trivial();

void unit_test_simd_fast_cos_approx_2xf64_trivial()
{
    using namespace gms::math;
    thread_local std::uniform_real_distribution<double> cos_rand_in;
    thread_local std::uint64_t seed{};
    constexpr double lo{-3.14159265358979324};
    constexpr double hi{+3.14159265358979324};
    seed = __rdtsc();
    cos_rand_in = std::uniform_real_distribution<double>(lo,hi);
    auto rand_uni_gen{std::mt19937(seed)};
    const double in_ref_1{cos_rand_in.operator()(rand_uni_gen)};
    const double in_ref_2{cos_rand_in.operator()(rand_uni_gen)};
    const __m128d    in_vref_vec{_mm_setr_pd(in_ref_1,in_ref_2)};
    __m128d vec_actual_vref{_mm_setzero_pd()};
    const double * __restrict__ p_vref_in = reinterpret_cast<const double*>(&in_vref_vec);
    double       * __restrict__ p_vref_out;
    double sca_expected_ref_1{};
    double sca_expected_ref_2{};
    printf("[UNIT_TEST]: function=%s -- **START**\n\n", __PRETTY_FUNCTION__);
    printf("[UNIT-TEST]: -- Calling reference implementation: cosapprox_d.\n");
    printf("[UNIT-TEST]: Inputs are: ref_1=%.17f,ref_2=%.17f\n",in_ref_1,in_ref_2);
    sca_expected_ref_1 = cosapprox_d(in_ref_1);
    sca_expected_ref_2 = cosapprox_d(in_ref_2);
    printf("[UNIT-TEST]: -- Calling SSE implementation: simd_fast_cos_approx_2xf64.\n");
    printf("[UNIT-TEST]: Inputs are: ref_1=%.17f,ref_2=%.17f\n",p_vref_in[0],p_vref_in[1]);
    //__asm__("int3");
    vec_actual_vref = simd_fast_cos_approx_2xf64(in_vref_vec);
    p_vref_out      = reinterpret_cast<double*>(&vec_actual_vref);
    //printf("[UNIT-TEST]: Calculating the results approximately equal by the epsilon:%.17f\n",std::numeric_limits<double>::epsilon());
    const bool ref_1_approx_equality{almostEqual2sComplement(sca_expected_ref_1,p_vref_out[0],17)};
    const bool ref_2_approx_equality{almostEqual2sComplement(sca_expected_ref_2,p_vref_out[1],17)};
    if(!ref_1_approx_equality) 
    {
       std::printf("[UNIT-TEST #1]-- **FAIL**, sca=%.17f,vec=%.17f,delta=%.17f\n",sca_expected_ref_1,p_vref_out[0],std::abs(sca_expected_ref_1-p_vref_out[0]));
       printNumber("cosapprox_d",sca_expected_ref_1,0);
       printNumber("simd_fast_cos_approx_2xf64",p_vref_out[0],0);
    }
    else 
    {
       std::printf("[UNIT-TEST #1]-- **SUCCESS**, sca=%.17f,vec=%.17f,delta=%.17f\n",sca_expected_ref_1,p_vref_out[0],std::abs(sca_expected_ref_1-p_vref_out[0]));
       printNumber("cosapprox_d",sca_expected_ref_1,0);
       printNumber("simd_fast_cos_approx_2xf64",p_vref_out[0],0);
    }
    if(!ref_2_approx_equality)
    {
       std::printf("[UNIT-TEST #2]-- **FAIL**, sca=%.17f,vec=%.17f,delta=%.17f\n",sca_expected_ref_2,p_vref_out[1],std::abs(sca_expected_ref_2-p_vref_out[1]));
       printNumber("cosapprox_d",sca_expected_ref_2,0);
       printNumber("simd_fast_cos_approx_2xf64",p_vref_out[1],0);
    }
    else 
    {
       std::printf("[UNIT-TEST #2]-- **SUCCESS**, sca=%.17f,vec=%.17f,delta=%.17f\n",sca_expected_ref_2,p_vref_out[1],std::abs(sca_expected_ref_2-p_vref_out[1]));
       printNumber("cosapprox_d",sca_expected_ref_2,0);
       printNumber("simd_fast_cos_approx_2xf64",p_vref_out[1],0);
    }
     printf("[UNIT_TEST]: function=%s -- **END**\n\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_simd_fast_log_approx_4xf32_trivial();

void unit_test_simd_fast_log_approx_4xf32_trivial()
{
    using namespace gms::math;
    thread_local std::uniform_real_distribution<float> log_rand_in;
    thread_local std::uint64_t seed{};
    constexpr float lo{0.1f};
    constexpr float hi{+2.718281828459045235360287471352f};
    seed = __rdtsc();
    log_rand_in = std::uniform_real_distribution<float>(lo,hi);
    auto rand_uni_gen{std::mt19937(seed)};
    const float in_ref_1{log_rand_in.operator()(rand_uni_gen)};
    const float in_ref_2{log_rand_in.operator()(rand_uni_gen)};
    const float in_ref_3{log_rand_in.operator()(rand_uni_gen)};
    const float in_ref_4{log_rand_in.operator()(rand_uni_gen)};
    const __m128    in_vref_vec{_mm_setr_ps(in_ref_1,in_ref_2,in_ref_3,in_ref_4)};
    __m128 vec_actual_vref{_mm_setzero_ps()};
    const float * __restrict__ p_vref_in = reinterpret_cast<const float*>(&in_vref_vec);
    float       * __restrict__ p_vref_out;
    float sca_expected_ref_1{};
    float sca_expected_ref_2{};
    float sca_expected_ref_3{};
    float sca_expected_ref_4{};
    printf("[UNIT_TEST]: function=%s -- **START**\n\n", __PRETTY_FUNCTION__);
    printf("[UNIT-TEST]: -- Calling reference implementation: logapprox.\n");
    printf("[UNIT-TEST]: Inputs are: ref_1=%.7f,ref_2=%.7f,ref_3=%.7f,ref_4=%.7f\n",in_ref_1,in_ref_2,in_ref_3,in_ref_4);
    sca_expected_ref_1 = logapprox(in_ref_1);
    sca_expected_ref_2 = logapprox(in_ref_2);
    sca_expected_ref_3 = logapprox(in_ref_3);
    sca_expected_ref_4 = logapprox(in_ref_4);
    printf("[UNIT-TEST]: -- Calling SSE implementation: simd_fast_log_approx_4xf32.\n");
    printf("[UNIT-TEST]: Inputs are: ref_1=%.7f,ref_2=%.7f,ref_3=%.7f,ref_4=%.7f\n",p_vref_in[0],p_vref_in[1],p_vref_in[2],p_vref_in[3]);
    //__asm__("int3");
    vec_actual_vref = simd_fast_log_approx_4xf32(in_vref_vec);
    p_vref_out      = reinterpret_cast<float*>(&vec_actual_vref);
    //printf("[UNIT-TEST]: Calculating the results approximately equal by the epsilon:%.17f\n",std::numeric_limits<double>::epsilon());
    const bool ref_1_approx_equality{almostEqual2sComplement(sca_expected_ref_1,p_vref_out[0],7)};
    const bool ref_2_approx_equality{almostEqual2sComplement(sca_expected_ref_2,p_vref_out[1],7)};
    const bool ref_3_approx_equality{almostEqual2sComplement(sca_expected_ref_3,p_vref_out[2],7)};
    const bool ref_4_approx_equality{almostEqual2sComplement(sca_expected_ref_4,p_vref_out[3],7)};
    if(!ref_1_approx_equality) 
    {
       std::printf("[UNIT-TEST #1]-- **FAIL**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_1,p_vref_out[0],std::abs(sca_expected_ref_1-p_vref_out[0]));
       printNumber("logapprox",sca_expected_ref_1,0);
       printNumber("simd_fast_log_approx_4xf32",p_vref_out[0],0);
    }
    else 
    {
       std::printf("[UNIT-TEST #1]-- **SUCCESS**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_1,p_vref_out[0],std::abs(sca_expected_ref_1-p_vref_out[0]));
       printNumber("logapprox",sca_expected_ref_1,0);
       printNumber("simd_fast_log_approx_4xf32",p_vref_out[0],0);
    }
    if(!ref_2_approx_equality)
    {
       std::printf("[UNIT-TEST #2]-- **FAIL**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_2,p_vref_out[1],std::abs(sca_expected_ref_2-p_vref_out[1]));
       printNumber("logapprox",sca_expected_ref_2,0);
       printNumber("simd_fast_log_approx_4xf32",p_vref_out[1],0);
    }
    else 
    {
       std::printf("[UNIT-TEST #2]-- **SUCCESS**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_2,p_vref_out[1],std::abs(sca_expected_ref_2-p_vref_out[1]));
       printNumber("logapprox",sca_expected_ref_2,0);
       printNumber("simd_fast_log_approx_4xf32",p_vref_out[1],0);
    }
    if(!ref_3_approx_equality) 
    {
       std::printf("[UNIT-TEST #3]-- **FAIL**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_3,p_vref_out[2],std::abs(sca_expected_ref_3-p_vref_out[2]));
       printNumber("logapprox",sca_expected_ref_3,0);
       printNumber("simd_fast_log_approx_4xf32",p_vref_out[2],0);
    }
    else 
    {
       std::printf("[UNIT-TEST #3]-- **SUCCESS**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_3,p_vref_out[2],std::abs(sca_expected_ref_3-p_vref_out[2]));
       printNumber("logapprox",sca_expected_ref_3,0);
       printNumber("simd_fast_log_approx_4xf32",p_vref_out[2],0);
    }
    if(!ref_4_approx_equality) 
    {
       std::printf("[UNIT-TEST #4]-- **FAIL**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_4,p_vref_out[3],std::abs(sca_expected_ref_4-p_vref_out[3]));
       printNumber("logapprox",sca_expected_ref_4,0);
       printNumber("simd_fast_log_approx_4xf32",p_vref_out[3],0);
    }
    else 
    {
       std::printf("[UNIT-TEST #4]-- **SUCCESS**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_4,p_vref_out[3],std::abs(sca_expected_ref_4-p_vref_out[3]));
       printNumber("logapprox",sca_expected_ref_4,0);
       printNumber("simd_fast_log_approx_4xf32",p_vref_out[3],0);
    }
     printf("[UNIT_TEST]: function=%s -- **END**\n\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_simd_fast_log_approx_2xf64_trivial();

void unit_test_simd_fast_log_approx_2xf64_trivial()
{
    using namespace gms::math;
    thread_local std::uniform_real_distribution<double> log_rand_in;
    thread_local std::uint64_t seed{};
    constexpr double lo{0.1};
    constexpr double hi{+2.718281828459045235360287471352};
    seed = __rdtsc();
    log_rand_in = std::uniform_real_distribution<double>(lo,hi);
    auto rand_uni_gen{std::mt19937(seed)};
    const double in_ref_1{log_rand_in.operator()(rand_uni_gen)};
    const double in_ref_2{log_rand_in.operator()(rand_uni_gen)};
    const __m128d    in_vref_vec{_mm_setr_pd(in_ref_1,in_ref_2)};
    __m128d vec_actual_vref{_mm_setzero_pd()};
    const double * __restrict__ p_vref_in = reinterpret_cast<const double*>(&in_vref_vec);
    double       * __restrict__ p_vref_out;
    double sca_expected_ref_1{};
    double sca_expected_ref_2{};
    printf("[UNIT_TEST]: function=%s -- **START**\n\n", __PRETTY_FUNCTION__);
    printf("[UNIT-TEST]: -- Calling reference implementation: logapprox_d.\n");
    printf("[UNIT-TEST]: Inputs are: ref_1=%.17f,ref_2=%.17f\n",in_ref_1,in_ref_2);
    sca_expected_ref_1 = logapprox(in_ref_1);
    sca_expected_ref_2 = logapprox(in_ref_2);
    printf("[UNIT-TEST]: -- Calling SSE implementation: simd_fast_log_approx_2xf64.\n");
    printf("[UNIT-TEST]: Inputs are: ref_1=%.17f,ref_2=%.17f\n",p_vref_in[0],p_vref_in[1]);
    //__asm__("int3");
    vec_actual_vref = simd_fast_log_approx_2xf64(in_vref_vec);
    p_vref_out      = reinterpret_cast<double*>(&vec_actual_vref);
    //printf("[UNIT-TEST]: Calculating the results approximately equal by the epsilon:%.17f\n",std::numeric_limits<double>::epsilon());
    const bool ref_1_approx_equality{almostEqual2sComplement(sca_expected_ref_1,p_vref_out[0],17)};
    const bool ref_2_approx_equality{almostEqual2sComplement(sca_expected_ref_2,p_vref_out[1],17)};
    if(!ref_1_approx_equality) 
    {
       std::printf("[UNIT-TEST #1]-- **FAIL**, sca=%.17f,vec=%.17f,delta=%.17f\n",sca_expected_ref_1,p_vref_out[0],std::abs(sca_expected_ref_1-p_vref_out[0]));
       printNumber("logapprox_d",sca_expected_ref_1,0);
       printNumber("simd_fast_log_approx_2xf64",p_vref_out[0],0);
    }
    else 
    {
       std::printf("[UNIT-TEST #1]-- **SUCCESS**, sca=%.17f,vec=%.17f,delta=%.17f\n",sca_expected_ref_1,p_vref_out[0],std::abs(sca_expected_ref_1-p_vref_out[0]));
       printNumber("logapprox_d",sca_expected_ref_1,0);
       printNumber("simd_fast_log_approx_2xf64",p_vref_out[0],0);
    }
    if(!ref_2_approx_equality)
    {
       std::printf("[UNIT-TEST #2]-- **FAIL**, sca=%.17f,vec=%.17f,delta=%.17f\n",sca_expected_ref_2,p_vref_out[1],std::abs(sca_expected_ref_2-p_vref_out[1]));
       printNumber("logapprox_d",sca_expected_ref_2,0);
       printNumber("simd_fast_log_approx_6xf64",p_vref_out[1],0);
    }
    else 
    {
       std::printf("[UNIT-TEST #2]-- **SUCCESS**, sca=%.17f,vec=%.17f,delta=%.17f\n",sca_expected_ref_2,p_vref_out[1],std::abs(sca_expected_ref_2-p_vref_out[1]));
       printNumber("logapprox_d",sca_expected_ref_2,0);
       printNumber("simd_fast_log_approx_2xf64",p_vref_out[1],0);
    }
     printf("[UNIT_TEST]: function=%s -- **END**\n\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_simd_fast_sin_approx_4xf32_trivial();

void unit_test_simd_fast_sin_approx_4xf32_trivial()
{
    using namespace gms::math;
    thread_local std::uniform_real_distribution<float> sin_rand_in;
    thread_local std::uint64_t seed{};
    constexpr float lo{-3.14159265358979324f};
    constexpr float hi{+3.14159265358979324f};
    seed = __rdtsc();
    sin_rand_in = std::uniform_real_distribution<float>(lo,hi);
    auto rand_uni_gen{std::mt19937(seed)};
    const float in_ref_1{sin_rand_in.operator()(rand_uni_gen)};
    const float in_ref_2{sin_rand_in.operator()(rand_uni_gen)};
    const float in_ref_3{sin_rand_in.operator()(rand_uni_gen)};
    const float in_ref_4{sin_rand_in.operator()(rand_uni_gen)};
    const __m128    in_vref_vec{_mm_setr_ps(in_ref_1,in_ref_2,in_ref_3,in_ref_4)};
    __m128 vec_actual_vref{_mm_setzero_ps()};
    const float * __restrict__ p_vref_in = reinterpret_cast<const float*>(&in_vref_vec);
    float       * __restrict__ p_vref_out;
    float sca_expected_ref_1{};
    float sca_expected_ref_2{};
    float sca_expected_ref_3{};
    float sca_expected_ref_4{};
    printf("[UNIT_TEST]: function=%s -- **START**\n\n", __PRETTY_FUNCTION__);
    printf("[UNIT-TEST]: -- Calling reference implementation: sinapprox.\n");
    printf("[UNIT-TEST]: Inputs are: ref_1=%.7f,ref_2=%.7f,ref_3=%.7f,ref_4=%.7f\n",in_ref_1,in_ref_2,in_ref_3,in_ref_4);
    sca_expected_ref_1 = sinapprox(in_ref_1);
    sca_expected_ref_2 = sinapprox(in_ref_2);
    sca_expected_ref_3 = sinapprox(in_ref_3);
    sca_expected_ref_4 = sinapprox(in_ref_4);
    printf("[UNIT-TEST]: -- Calling SSE implementation: simd_fast_sin_approx_4xf32.\n");
    printf("[UNIT-TEST]: Inputs are: ref_1=%.7f,ref_2=%.7f,ref_3=%.7f,ref_4=%.7f\n",p_vref_in[0],p_vref_in[1],p_vref_in[2],p_vref_in[3]);
    //__asm__("int3");
    vec_actual_vref = simd_fast_sin_approx_4xf32(in_vref_vec);
    p_vref_out      = reinterpret_cast<float*>(&vec_actual_vref);
    //printf("[UNIT-TEST]: Calculating the results approximately equal by the epsilon:%.17f\n",std::numeric_limits<double>::epsilon());
    const bool ref_1_approx_equality{almostEqual2sComplement(sca_expected_ref_1,p_vref_out[0],7)};
    const bool ref_2_approx_equality{almostEqual2sComplement(sca_expected_ref_2,p_vref_out[1],7)};
    const bool ref_3_approx_equality{almostEqual2sComplement(sca_expected_ref_3,p_vref_out[2],7)};
    const bool ref_4_approx_equality{almostEqual2sComplement(sca_expected_ref_4,p_vref_out[3],7)};
    if(!ref_1_approx_equality) 
    {
       std::printf("[UNIT-TEST #1]-- **FAIL**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_1,p_vref_out[0],std::abs(sca_expected_ref_1-p_vref_out[0]));
       printNumber("sinapprox",sca_expected_ref_1,0);
       printNumber("simd_fast_sin_approx_4xf32",p_vref_out[0],0);
    }
    else 
    {
       std::printf("[UNIT-TEST #1]-- **SUCCESS**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_1,p_vref_out[0],std::abs(sca_expected_ref_1-p_vref_out[0]));
       printNumber("sinapprox",sca_expected_ref_1,0);
       printNumber("simd_fast_sin_approx_4xf32",p_vref_out[0],0);
    }
    if(!ref_2_approx_equality)
    {
       std::printf("[UNIT-TEST #2]-- **FAIL**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_2,p_vref_out[1],std::abs(sca_expected_ref_2-p_vref_out[1]));
       printNumber("sinapprox",sca_expected_ref_2,0);
       printNumber("simd_fast_sin_approx_4xf32",p_vref_out[1],0);
    }
    else 
    {
       std::printf("[UNIT-TEST #2]-- **SUCCESS**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_2,p_vref_out[1],std::abs(sca_expected_ref_2-p_vref_out[1]));
       printNumber("sinapprox",sca_expected_ref_2,0);
       printNumber("simd_fast_sin_approx_4xf32",p_vref_out[1],0);
    }
    if(!ref_3_approx_equality) 
    {
       std::printf("[UNIT-TEST #3]-- **FAIL**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_3,p_vref_out[2],std::abs(sca_expected_ref_3-p_vref_out[2]));
       printNumber("sinapprox",sca_expected_ref_3,0);
       printNumber("simd_fast_sin_approx_4xf32",p_vref_out[2],0);
    }
    else 
    {
       std::printf("[UNIT-TEST #3]-- **SUCCESS**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_3,p_vref_out[2],std::abs(sca_expected_ref_3-p_vref_out[2]));
       printNumber("sinapprox",sca_expected_ref_3,0);
       printNumber("simd_fast_sin_approx_4xf32",p_vref_out[2],0);
    }
    if(!ref_4_approx_equality) 
    {
       std::printf("[UNIT-TEST #4]-- **FAIL**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_4,p_vref_out[3],std::abs(sca_expected_ref_4-p_vref_out[3]));
       printNumber("sinapprox",sca_expected_ref_4,0);
       printNumber("simd_fast_sin_approx_4xf32",p_vref_out[3],0);
    }
    else 
    {
       std::printf("[UNIT-TEST #4]-- **SUCCESS**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_4,p_vref_out[3],std::abs(sca_expected_ref_4-p_vref_out[3]));
       printNumber("sinapprox",sca_expected_ref_4,0);
       printNumber("simd_fast_sin_approx_4xf32",p_vref_out[3],0);
    }
     printf("[UNIT_TEST]: function=%s -- **END**\n\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_simd_fast_sin_approx_2xf64_trivial();

void unit_test_simd_fast_sin_approx_2xf64_trivial()
{
    using namespace gms::math;
    thread_local std::uniform_real_distribution<double> sin_rand_in;
    thread_local std::uint64_t seed{};
    constexpr double lo{-3.14159265358979324};
    constexpr double hi{+3.14159265358979324};
    seed = __rdtsc();
    sin_rand_in = std::uniform_real_distribution<double>(lo,hi);
    auto rand_uni_gen{std::mt19937(seed)};
    const double in_ref_1{sin_rand_in.operator()(rand_uni_gen)};
    const double in_ref_2{sin_rand_in.operator()(rand_uni_gen)};
    const __m128d    in_vref_vec{_mm_setr_pd(in_ref_1,in_ref_2)};
    __m128d vec_actual_vref{_mm_setzero_pd()};
    const double * __restrict__ p_vref_in = reinterpret_cast<const double*>(&in_vref_vec);
    double       * __restrict__ p_vref_out;
    double sca_expected_ref_1{};
    double sca_expected_ref_2{};
    printf("[UNIT_TEST]: function=%s -- **START**\n\n", __PRETTY_FUNCTION__);
    printf("[UNIT-TEST]: -- Calling reference implementation: sinapprox_d.\n");
    printf("[UNIT-TEST]: Inputs are: ref_1=%.17f,ref_2=%.17f\n",in_ref_1,in_ref_2);
    sca_expected_ref_1 = sinapprox_d(in_ref_1);
    sca_expected_ref_2 = sinapprox_d(in_ref_2);
    printf("[UNIT-TEST]: -- Calling SSE implementation: simd_fast_sin_approx_2xf64.\n");
    printf("[UNIT-TEST]: Inputs are: ref_1=%.17f,ref_2=%.17f\n",p_vref_in[0],p_vref_in[1]);
    //__asm__("int3");
    vec_actual_vref = simd_fast_sin_approx_2xf64(in_vref_vec);
    p_vref_out      = reinterpret_cast<double*>(&vec_actual_vref);
    //printf("[UNIT-TEST]: Calculating the results approximately equal by the epsilon:%.17f\n",std::numeric_limits<double>::epsilon());
    const bool ref_1_approx_equality{almostEqual2sComplement(sca_expected_ref_1,p_vref_out[0],17)};
    const bool ref_2_approx_equality{almostEqual2sComplement(sca_expected_ref_2,p_vref_out[1],17)};
    if(!ref_1_approx_equality) 
    {
       std::printf("[UNIT-TEST #1]-- **FAIL**, sca=%.17f,vec=%.17f,delta=%.17f\n",sca_expected_ref_1,p_vref_out[0],std::abs(sca_expected_ref_1-p_vref_out[0]));
       printNumber("sinapprox_d",sca_expected_ref_1,0);
       printNumber("simd_fast_sin_approx_2xf64",p_vref_out[0],0);
    }
    else 
    {
       std::printf("[UNIT-TEST #1]-- **SUCCESS**, sca=%.17f,vec=%.17f,delta=%.17f\n",sca_expected_ref_1,p_vref_out[0],std::abs(sca_expected_ref_1-p_vref_out[0]));
       printNumber("sinapprox_d",sca_expected_ref_1,0);
       printNumber("simd_fast_sin_approx_2xf64",p_vref_out[0],0);
    }
    if(!ref_2_approx_equality)
    {
       std::printf("[UNIT-TEST #2]-- **FAIL**, sca=%.17f,vec=%.17f,delta=%.17f\n",sca_expected_ref_2,p_vref_out[1],std::abs(sca_expected_ref_2-p_vref_out[1]));
       printNumber("sinapprox_d",sca_expected_ref_2,0);
       printNumber("simd_fast_sin_approx_2xf64",p_vref_out[1],0);
    }
    else 
    {
       std::printf("[UNIT-TEST #2]-- **SUCCESS**, sca=%.17f,vec=%.17f,delta=%.17f\n",sca_expected_ref_2,p_vref_out[1],std::abs(sca_expected_ref_2-p_vref_out[1]));
       printNumber("sinapprox_d",sca_expected_ref_2,0);
       printNumber("simd_fast_sin_approx_2xf64",p_vref_out[1],0);
    }
     printf("[UNIT_TEST]: function=%s -- **END**\n\n", __PRETTY_FUNCTION__);
}

int main()
{
   unit_test_simd_fast_exp_approx_4xf32_trivial();
   unit_test_simd_fast_exp_approx_2xf64_trivial();
   unit_test_simd_fast_cos_approx_4xf32_trivial();
   unit_test_simd_fast_cos_approx_2xf64_trivial();
   unit_test_simd_fast_log_approx_4xf32_trivial();
   unit_test_simd_fast_log_approx_2xf64_trivial();
   unit_test_simd_fast_sin_approx_4xf32_trivial();
   unit_test_simd_fast_sin_approx_2xf64_trivial();
    return 0;
}