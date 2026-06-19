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

float exp_cst1_f = 2139095040.f;
float exp_cst2_f = 0.f;

/* Relative error bounded by 1e-5 for normalized outputs
   Returns invalid outputs for nan inputs
   Continuous error */
static 
inline float expapprox(float val) {
  union { int32_t i; float f; } xu, xu2;
  float val2, val3, val4, b;
  int32_t val4i;
  val2 = 12102203.1615614f*val+1065353216.f;
  val3 = val2 < exp_cst1_f ? val2 : exp_cst1_f;
  val4 = val3 > exp_cst2_f ? val3 : exp_cst2_f;
  val4i = (int32_t) val4;
  xu.i = val4i & 0x7F800000;
  xu2.i = (val4i & 0x7FFFFF) | 0x3F800000;
  b = xu2.f;

  /* Generated in Sollya with:
     > f=remez(1-x*exp(-(x-1)*log(2)),
               [|(x-1)*(x-2), (x-1)*(x-2)*x, (x-1)*(x-2)*x*x|],
               [1.000001,1.999999], exp(-(x-1)*log(2)));
     > plot(exp((x-1)*log(2))/(f+x)-1, [1,2]);
     > f+x;
  */
  return
    xu.f * (0.509871020f + b * (0.312146713f + b * (0.166617139f + b *
              (-2.190619930e-3f + b * 1.3555747234e-2f))));
}

double exp_cst1_d = 9218868437227405312.;
double exp_cst2_d = 0.;

/* Relative error bounded by 3e-9 for normalized outputs
   Returns invalid outputs for nan inputs
   Continuous error

   Vectorizable only with AVX512dq extensions because of the
   double->int64 cast. On GCC, use option -mavx512dq. */
inline double expapprox_d(double val) {
  union { int64_t i; double f; } xu, xu2;
  double val2, val3, val4, b;
  int64_t val4i;
  val2 = 6497320848556798.092*val+4607182418800017408.;
  val3 = val2 < exp_cst1_d ? val2 : exp_cst1_d;
  val4 = val3 > exp_cst2_d ? val3 : exp_cst2_d;
  val4i = (int64_t) val4;
  xu.i = val4i & 0x7FF0000000000000;
  xu2.i = (val4i & 0xFFFFFFFFFFFFF) | 0x3FF0000000000000;
  b = xu2.f;

  /* Generated in Sollya with:
     > f=remez(1-x*exp(-(x-1)*log(2)),
               [|(x-1)*(x-2), (x-1)*(x-2)*x, (x-1)*(x-2)*x*x,
                 (x-1)*(x-2)*x*x*x, (x-1)*(x-2)*x*x*x*x|],
                 [1.000001,1.999999], exp(-(x-1)*log(2)));
     > plot(exp((x-1)*log(2))/(f+x)-1, [1,2]);
     > f+x;
  */
  return
    xu.f * (0.5002494548377929861 + b * (0.3453457447382168695 + b *
            (0.1226618159032020501 + b * (2.4869768911673215212e-2 + b *
             (6.7148791796145116483e-3 + b * (-5.8813825553864185693e-5 + b *
              2.17150255054231565039e-4))))));
}

/* Correct only in [-pi, pi]
   Absolute error bounded by 5e-5
   Continuous error */
inline float cosapprox(float val) {
  float val2 = val*val;
  /* Generated in Sollya using:
     > f = remez(cos(x)-1, [|x*x, x*x*x*x, x*x*x*x*x*x, x*x*x*x*x*x*x*x|],
                            [0.000001, pi], 1, 1e-8);
     > plot(f-cos(x)+1, [0, pi]);
     > f+1
  */
  return
    1.f + val2 * (-0.4998515820f + val2 * (4.1518035216e-2f + val2 *
      (-1.3422947025e-3f + val2 * 1.8929864824e-5f)));
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
       printNumber("expapprox",sca_expected_ref_1,0);
       printNumber("simd_fast_cos_approx_4xf32",p_vref_out[0],0);
    }
    else 
    {
       std::printf("[UNIT-TEST #1]-- **SUCCESS**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_1,p_vref_out[0],std::abs(sca_expected_ref_1-p_vref_out[0]));
       printNumber("expapprox",sca_expected_ref_1,0);
       printNumber("simd_fast_cos_approx_4xf32",p_vref_out[0],0);
    }
    if(!ref_2_approx_equality)
    {
       std::printf("[UNIT-TEST #2]-- **FAIL**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_2,p_vref_out[1],std::abs(sca_expected_ref_2-p_vref_out[1]));
       printNumber("expapprox",sca_expected_ref_2,0);
       printNumber("simd_fast_cos_approx_4xf32",p_vref_out[1],0);
    }
    else 
    {
       std::printf("[UNIT-TEST #2]-- **SUCCESS**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_2,p_vref_out[1],std::abs(sca_expected_ref_2-p_vref_out[1]));
       printNumber("expapprox",sca_expected_ref_2,0);
       printNumber("simd_fast_cos_approx_4xf32",p_vref_out[1],0);
    }
    if(!ref_3_approx_equality) 
    {
       std::printf("[UNIT-TEST #3]-- **FAIL**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_3,p_vref_out[2],std::abs(sca_expected_ref_3-p_vref_out[2]));
       printNumber("expapprox",sca_expected_ref_3,0);
       printNumber("simd_fast_cos_approx_4xf32",p_vref_out[2],0);
    }
    else 
    {
       std::printf("[UNIT-TEST #3]-- **SUCCESS**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_3,p_vref_out[2],std::abs(sca_expected_ref_3-p_vref_out[2]));
       printNumber("expapprox",sca_expected_ref_3,0);
       printNumber("simd_fast_cos_approx_4xf32",p_vref_out[2],0);
    }
    if(!ref_4_approx_equality) 
    {
       std::printf("[UNIT-TEST #4]-- **FAIL**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_4,p_vref_out[3],std::abs(sca_expected_ref_4-p_vref_out[3]));
       printNumber("expapprox",sca_expected_ref_4,0);
       printNumber("simd_fast_cos_approx_4xf32",p_vref_out[3],0);
    }
    else 
    {
       std::printf("[UNIT-TEST #4]-- **SUCCESS**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_4,p_vref_out[3],std::abs(sca_expected_ref_4-p_vref_out[3]));
       printNumber("expapprox",sca_expected_ref_4,0);
       printNumber("simd_fast_cos_approx_4xf32",p_vref_out[3],0);
    }
     printf("[UNIT_TEST]: function=%s -- **END**\n\n", __PRETTY_FUNCTION__);
}


int main()
{
    unit_test_simd_fast_exp_approx_4xf32_trivial();
    unit_test_simd_fast_exp_approx_2xf64_trivial();
    unit_test_simd_fast_cos_approx_4xf32_trivial();
    return 0;
}