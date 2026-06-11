#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <cmath>
#include "GMS_ansi_colors.h"
#include "GMS_cephes_sin_cos.h"
#include "GMS_simd_cephes_cos.h"
/*
   icpc -o unit_test_mm256_cephes_cos_ps -O3 -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_ansi_colors.h GMS_cephes_sin_cos.h GMS_simd_cephes_cos_ps.h unit_test_mm256_cephes_cos_ps.cpp
   ASM: 
   icpc -S --O3 fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_ansi_colors.h GMS_cephes_sin_cos.h GMS_simd_cephes_cos_ps.h unit_test_mm256_cephes_cos_ps.cpp

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

__attribute__((hot))
__attribute__((noinline))
void unit_test_mm256_cephes_cos_ps_trivial();

void unit_test_mm256_cephes_cos_ps_trivial()
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
    const float in_ref_5{cos_rand_in.operator()(rand_uni_gen)};
    const float in_ref_6{cos_rand_in.operator()(rand_uni_gen)};
    const float in_ref_7{cos_rand_in.operator()(rand_uni_gen)};
    const float in_ref_8{cos_rand_in.operator()(rand_uni_gen)};
    const __m256    in_vref_vec{_mm256_setr_ps(in_ref_1,in_ref_2,in_ref_3,in_ref_4,
                                               in_ref_5,in_ref_6,in_ref_7,in_ref_8)};
    __m256 vec_actual_vref{_mm256_setzero_ps()};
    const float * __restrict__ p_vref_in = reinterpret_cast<const float*>(&in_vref_vec);
    float       * __restrict__ p_vref_out;
    float sca_expected_ref_1{};
    float sca_expected_ref_2{};
    float sca_expected_ref_3{};
    float sca_expected_ref_4{};
    float sca_expected_ref_5{};
    float sca_expected_ref_6{};
    float sca_expected_ref_7{};
    float sca_expected_ref_8{};
    printf("[UNIT_TEST]: function=%s -- **START**\n\n", __PRETTY_FUNCTION__);
    printf("[UNIT-TEST]: -- Calling reference implementation: ceph_cosf.\n");
    printf("[UNIT-TEST]: Inputs (1st set) are: ref_1=%.7f,ref_2=%.7f,ref_3=%.7f,ref_4=%.7f\n",in_ref_1,in_ref_2,in_ref_3,in_ref_4);
    sca_expected_ref_1 = ceph_cosf(in_ref_1);
    sca_expected_ref_2 = ceph_cosf(in_ref_2);
    sca_expected_ref_3 = ceph_cosf(in_ref_3);
    sca_expected_ref_4 = ceph_cosf(in_ref_4);
    printf("[UNIT-TEST]: Inputs (2nd set) are: ref_5=%.7f,ref_6=%.7f,ref_7=%.7f,ref_8=%.7f\n",in_ref_5,in_ref_6,in_ref_7,in_ref_8);
    sca_expected_ref_5 = ceph_cosf(in_ref_5);
    sca_expected_ref_6 = ceph_cosf(in_ref_6);
    sca_expected_ref_7 = ceph_cosf(in_ref_7);
    sca_expected_ref_8 = ceph_cosf(in_ref_8);
    printf("[UNIT-TEST]: -- Calling AVX/AVX2 implementation: _mm256_ceph_cos_ps.\n");
    printf("[UNIT-TEST]: Inputs (low -YMM) are: ref_1=%.7f,ref_2=%.7f,ref_3=%.7f,ref_4=%.7f\n",p_vref_in[0],p_vref_in[1],p_vref_in[2],p_vref_in[3]);
    printf("[UNIT-TEST]: Inputs (high-YMM) are: ref_5=%.7f,ref_6=%.7f,ref_7=%.7f,ref_8=%.7f\n",p_vref_in[4],p_vref_in[5],p_vref_in[6],p_vref_in[7]);
    //__asm__("int3");
    vec_actual_vref = _mm256_ceph_cosf_ps(in_vref_vec);
    p_vref_out      = reinterpret_cast<float*>(&vec_actual_vref);
    //printf("[UNIT-TEST]: Calculating the results approximately equal by the epsilon:%.17f\n",std::numeric_limits<double>::epsilon());
    const bool ref_1_approx_equality{almostEqual2sComplement(sca_expected_ref_1,p_vref_out[0],7)};
    const bool ref_2_approx_equality{almostEqual2sComplement(sca_expected_ref_2,p_vref_out[1],7)};
    const bool ref_3_approx_equality{almostEqual2sComplement(sca_expected_ref_3,p_vref_out[2],7)};
    const bool ref_4_approx_equality{almostEqual2sComplement(sca_expected_ref_4,p_vref_out[3],7)};
    const bool ref_5_approx_equality{almostEqual2sComplement(sca_expected_ref_5,p_vref_out[4],7)};
    const bool ref_6_approx_equality{almostEqual2sComplement(sca_expected_ref_6,p_vref_out[5],7)};
    const bool ref_7_approx_equality{almostEqual2sComplement(sca_expected_ref_7,p_vref_out[6],7)};
    const bool ref_8_approx_equality{almostEqual2sComplement(sca_expected_ref_8,p_vref_out[7],7)};
    if(!ref_1_approx_equality) 
       {printf(ANSI_COLOR_RED "[UNIT-TEST #1]-- **FAILED**, sca=%.7f,vec=%.7f,delta=%.7f" ANSI_RESET_ALL "\n",sca_expected_ref_1,p_vref_out[0],std::abs(sca_expected_ref_1-p_vref_out[0]));}
    if(!ref_2_approx_equality) 
       {printf(ANSI_COLOR_RED "[UNIT-TEST #2]-- **FAILED**, sca=%.7f,vec=%.7f,delta=%.7f" ANSI_RESET_ALL "\n",sca_expected_ref_2,p_vref_out[1],std::abs(sca_expected_ref_2-p_vref_out[1]));}
    if(!ref_3_approx_equality) 
       {printf(ANSI_COLOR_RED "[UNIT-TEST #3]-- **FAILED**, sca=%.7f,vec=%.7f,delta=%.7f" ANSI_RESET_ALL "\n",sca_expected_ref_3,p_vref_out[2],std::abs(sca_expected_ref_3-p_vref_out[2]));}
    if(!ref_4_approx_equality) 
       {printf(ANSI_COLOR_RED "[UNIT-TEST #4]-- **FAILED**, sca=%.7f,vec=%.7f,delta=%.7f" ANSI_RESET_ALL "\n",sca_expected_ref_4,p_vref_out[3],std::abs(sca_expected_ref_4-p_vref_out[3]));}
    if(!ref_5_approx_equality) 
       {printf(ANSI_COLOR_RED "[UNIT-TEST #5]-- **FAILED**, sca=%.7f,vec=%.7f,delta=%.7f" ANSI_RESET_ALL "\n",sca_expected_ref_5,p_vref_out[4],std::abs(sca_expected_ref_5-p_vref_out[4]));}
    if(!ref_6_approx_equality) 
       {printf(ANSI_COLOR_RED "[UNIT-TEST #6]-- **FAILED**, sca=%.7f,vec=%.7f,delta=%.7f" ANSI_RESET_ALL "\n",sca_expected_ref_6,p_vref_out[5],std::abs(sca_expected_ref_6-p_vref_out[5]));}
    if(!ref_7_approx_equality) 
       {printf(ANSI_COLOR_RED "[UNIT-TEST #7]-- **FAILED**, sca=%.7f,vec=%.7f,delta=%.7f" ANSI_RESET_ALL "\n",sca_expected_ref_7,p_vref_out[6],std::abs(sca_expected_ref_7-p_vref_out[6]));}
    if(!ref_8_approx_equality) 
       {printf(ANSI_COLOR_RED "[UNIT-TEST #8]-- **FAILED**, sca=%.7f,vec=%.7f,delta=%.7f" ANSI_RESET_ALL "\n",sca_expected_ref_8,p_vref_out[7],std::abs(sca_expected_ref_8-p_vref_out[7]));}  
     printf("[UNIT_TEST]: function=%s -- **END**\n\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_mm256_cephes_cos_ps_v2_trivial();

void unit_test_mm256_cephes_cos_ps_v2_trivial()
{
    using namespace gms::math;
    thread_local std::uniform_real_distribution<float> cos_rand_in;
    thread_local std::uint64_t seed{};
    constexpr float lo{-3.14159265358979324f};
    constexpr float hi{1000.0f*lo};
    seed = __rdtsc();
    cos_rand_in = std::uniform_real_distribution<float>(lo,hi);
    auto rand_uni_gen{std::mt19937(seed)};
    const float in_ref_1{cos_rand_in.operator()(rand_uni_gen)};
    const float in_ref_2{cos_rand_in.operator()(rand_uni_gen)};
    const float in_ref_3{cos_rand_in.operator()(rand_uni_gen)};
    const float in_ref_4{cos_rand_in.operator()(rand_uni_gen)};
    const float in_ref_5{cos_rand_in.operator()(rand_uni_gen)};
    const float in_ref_6{cos_rand_in.operator()(rand_uni_gen)};
    const float in_ref_7{cos_rand_in.operator()(rand_uni_gen)};
    const float in_ref_8{cos_rand_in.operator()(rand_uni_gen)};
    const __m256    in_vref_vec{_mm256_setr_ps(in_ref_1,in_ref_2,in_ref_3,in_ref_4,
                                               in_ref_5,in_ref_6,in_ref_7,in_ref_8)};
    __m256 vec_actual_vref{_mm256_setzero_ps()};
    const float * __restrict__ p_vref_in = reinterpret_cast<const float*>(&in_vref_vec);
    float       * __restrict__ p_vref_out;
    float sca_expected_ref_1{};
    float sca_expected_ref_2{};
    float sca_expected_ref_3{};
    float sca_expected_ref_4{};
    float sca_expected_ref_5{};
    float sca_expected_ref_6{};
    float sca_expected_ref_7{};
    float sca_expected_ref_8{};
    printf("[UNIT_TEST]: function=%s -- **START**\n\n", __PRETTY_FUNCTION__);
    printf("[UNIT-TEST]: -- Calling reference implementation: ceph_cosf.\n");
    printf("[UNIT-TEST]: Inputs (1st set) are: ref_1=%.7f,ref_2=%.7f,ref_3=%.7f,ref_4=%.7f\n",in_ref_1,in_ref_2,in_ref_3,in_ref_4);
    sca_expected_ref_1 = ceph_cosf(in_ref_1);
    sca_expected_ref_2 = ceph_cosf(in_ref_2);
    sca_expected_ref_3 = ceph_cosf(in_ref_3);
    sca_expected_ref_4 = ceph_cosf(in_ref_4);
    printf("[UNIT-TEST]: Inputs (2nd set) are: ref_5=%.7f,ref_6=%.7f,ref_7=%.7f,ref_8=%.7f\n",in_ref_5,in_ref_6,in_ref_7,in_ref_8);
    sca_expected_ref_5 = ceph_cosf(in_ref_5);
    sca_expected_ref_6 = ceph_cosf(in_ref_6);
    sca_expected_ref_7 = ceph_cosf(in_ref_7);
    sca_expected_ref_8 = ceph_cosf(in_ref_8);
    printf("[UNIT-TEST]: -- Calling AVX/AVX2 implementation: _mm256_ceph_cos_ps.\n");
    printf("[UNIT-TEST]: Inputs (low -YMM) are: ref_1=%.7f,ref_2=%.7f,ref_3=%.7f,ref_4=%.7f\n",p_vref_in[0],p_vref_in[1],p_vref_in[2],p_vref_in[3]);
    printf("[UNIT-TEST]: Inputs (high-YMM) are: ref_5=%.7f,ref_6=%.7f,ref_7=%.7f,ref_8=%.7f\n",p_vref_in[4],p_vref_in[5],p_vref_in[6],p_vref_in[7]);
    //__asm__("int3");
    vec_actual_vref = _mm256_ceph_cosf_ps_v2<true,true>(in_vref_vec);
    p_vref_out      = reinterpret_cast<float*>(&vec_actual_vref);
    //printf("[UNIT-TEST]: Calculating the results approximately equal by the epsilon:%.17f\n",std::numeric_limits<double>::epsilon());
    const bool ref_1_approx_equality{almostEqual2sComplement(sca_expected_ref_1,p_vref_out[0],7)};
    const bool ref_2_approx_equality{almostEqual2sComplement(sca_expected_ref_2,p_vref_out[1],7)};
    const bool ref_3_approx_equality{almostEqual2sComplement(sca_expected_ref_3,p_vref_out[2],7)};
    const bool ref_4_approx_equality{almostEqual2sComplement(sca_expected_ref_4,p_vref_out[3],7)};
    const bool ref_5_approx_equality{almostEqual2sComplement(sca_expected_ref_5,p_vref_out[4],7)};
    const bool ref_6_approx_equality{almostEqual2sComplement(sca_expected_ref_6,p_vref_out[5],7)};
    const bool ref_7_approx_equality{almostEqual2sComplement(sca_expected_ref_7,p_vref_out[6],7)};
    const bool ref_8_approx_equality{almostEqual2sComplement(sca_expected_ref_8,p_vref_out[7],7)};
    if(!ref_1_approx_equality) 
       {printf(ANSI_COLOR_RED "[UNIT-TEST #1]-- **FAILED**, sca=%.7f,vec=%.7f,delta=%.7f" ANSI_RESET_ALL "\n",sca_expected_ref_1,p_vref_out[0],std::abs(sca_expected_ref_1-p_vref_out[0]));}
    if(!ref_2_approx_equality) 
       {printf(ANSI_COLOR_RED "[UNIT-TEST #2]-- **FAILED**, sca=%.7f,vec=%.7f,delta=%.7f" ANSI_RESET_ALL "\n",sca_expected_ref_2,p_vref_out[1],std::abs(sca_expected_ref_2-p_vref_out[1]));}
    if(!ref_3_approx_equality) 
       {printf(ANSI_COLOR_RED "[UNIT-TEST #3]-- **FAILED**, sca=%.7f,vec=%.7f,delta=%.7f" ANSI_RESET_ALL "\n",sca_expected_ref_3,p_vref_out[2],std::abs(sca_expected_ref_3-p_vref_out[2]));}
    if(!ref_4_approx_equality) 
       {printf(ANSI_COLOR_RED "[UNIT-TEST #4]-- **FAILED**, sca=%.7f,vec=%.7f,delta=%.7f" ANSI_RESET_ALL "\n",sca_expected_ref_4,p_vref_out[3],std::abs(sca_expected_ref_4-p_vref_out[3]));}
    if(!ref_5_approx_equality) 
       {printf(ANSI_COLOR_RED "[UNIT-TEST #5]-- **FAILED**, sca=%.7f,vec=%.7f,delta=%.7f" ANSI_RESET_ALL "\n",sca_expected_ref_5,p_vref_out[4],std::abs(sca_expected_ref_5-p_vref_out[4]));}
    if(!ref_6_approx_equality) 
       {printf(ANSI_COLOR_RED "[UNIT-TEST #6]-- **FAILED**, sca=%.7f,vec=%.7f,delta=%.7f" ANSI_RESET_ALL "\n",sca_expected_ref_6,p_vref_out[5],std::abs(sca_expected_ref_6-p_vref_out[5]));}
    if(!ref_7_approx_equality) 
       {printf(ANSI_COLOR_RED "[UNIT-TEST #7]-- **FAILED**, sca=%.7f,vec=%.7f,delta=%.7f" ANSI_RESET_ALL "\n",sca_expected_ref_7,p_vref_out[6],std::abs(sca_expected_ref_7-p_vref_out[6]));}
    if(!ref_8_approx_equality) 
       {printf(ANSI_COLOR_RED "[UNIT-TEST #8]-- **FAILED**, sca=%.7f,vec=%.7f,delta=%.7f" ANSI_RESET_ALL "\n",sca_expected_ref_8,p_vref_out[7],std::abs(sca_expected_ref_8-p_vref_out[7]));}  
     printf("[UNIT_TEST]: function=%s -- **END**\n\n", __PRETTY_FUNCTION__);
}

int main()
{
    unit_test_mm256_cephes_cos_ps_trivial();
    unit_test_mm256_cephes_cos_ps_v2_trivial();
    return 0;
}
