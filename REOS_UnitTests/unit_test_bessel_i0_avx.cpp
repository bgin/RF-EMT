#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <cmath>
#include "GMS_bessel_i0_avx.h"

/*
   icpc -o unit_test_bessel_i0_avx -O3 -fp-model fast=2 -fno-exceptions -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_bessel_i0_avx.h unit_test_bessel_i0_avx.cpp
   ASM: 
   icpc -S -O3 fverbose-asm -masm=intel -fno-exceptions -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_bessel_i0_avx.h unit_test_bessel_i0_avx.cpp

*/
static inline 
double bessel_i0(double arg, int jint) {
    // Mathematical constants
    const double one = 1.0;
    const double one5 = 15.0;
    const double exp40 = 2.353852668370199854e17;
    const double forty = 40.0;
    const double rec15 = 6.6666666666666666666e-2;
    const double two25 = 225.0;

    // Machine-dependent constants
    const double xsmall = 5.55e-17;
    const double xinf = std::numeric_limits<double>::infinity();
    const double xmax = 713.986;

    // Coefficients for XSMALL <= ABS(ARG) < 15.0
    const double p[15] = {
        -5.2487866627945699800e-18, -1.5982226675653184646e-14,
        -2.6843448573468483278e-11, -3.0517226450451067446e-08,
        -2.5172644670688975051e-05, -1.5453977791786851041e-02,
        -7.0935347449210549190e+00, -2.4125195876041896775e+03,
        -5.9545626019847898221e+05, -1.0313066708737980747e+08,
        -1.1912746104985237192e+10, -8.4925101247114157499e+11,
        -3.2940087627407749166e+13, -5.5050369673018427753e+14,
        -2.2335582639474375249e+15
    };
    const double q[5] = {
        -3.7277560179962773046e+03,  6.5158506418655165707e+06,
        -6.5626560740833869295e+09,  3.7604188704092954661e+12,
        -9.7087946179594019126e+14
    };

    // Coefficients for 15.0 <= ABS(ARG)
    const double pp[8] = {
        -3.9843750000000000000e-01,  2.9205384596336793945e+00,
        -2.4708469169133954315e+00,  4.7914889422856814203e-01,
        -3.7384991926068969150e-03, -2.6801520353328635310e-03,
         9.9168777670983678974e-05, -2.1877128189032726730e-06
    };
    const double qq[7] = {
        -3.1446690275135491500e+01,  8.5539563258012929600e+01,
        -6.0228002066743340583e+01,  1.3982595353892851542e+01,
        -1.1151759188741312645e+00,  3.2547697594819615062e-02,
        -5.5194330231005480228e-04
    };

    double x = std::abs(arg);
    double result = 0.0;

    if (x < xsmall) {
        result = one;
    } 
    // XSMALL <= ABS(ARG) < 15.0
    else if (x < one5) {
        double xx = x * x;
        double sump = p[0];
        for (int i = 1; i < 15; ++i) {
            sump = sump * xx + p[i];
        }
        xx = xx - two25;

        double sumq = ((((xx + q[0]) * xx + q[1]) * xx + q[2]) * xx + q[3]) * xx + q[4];

        result = sump / sumq;

        if (jint == 2) {
            result = result * std::exp(-x);
        }
    } 
    // 15.0 <= ABS(ARG)
    else if (one5 <= x) {
        if (jint == 1 && xmax < x) {
            result = xinf;
        } else {
            double xx = one / x - rec15;

            double sump = ((((((pp[0] * xx + pp[1]) * xx + pp[2]) * xx + pp[3]) * xx + pp[4]) * xx + pp[5]) * xx + pp[6]) * xx + pp[7];
            double sumq = ((((((xx + qq[0]) * xx + qq[1]) * xx + qq[2]) * xx + qq[3]) * xx + qq[4]) * xx + qq[5]) * xx + qq[6];

            result = sump / sumq;

            if (jint == 2) {
                result = (result - pp[0]) / std::sqrt(x);
            } else {
                // Calculation reformulated to avoid premature overflow
                double a, b;
                if (x <= (xmax - one5)) {
                    a = std::exp(x);
                    b = one;
                } else {
                    a = std::exp(x - forty);
                    b = exp40;
                }

                result = ((result * a - pp[0] * a) / std::sqrt(x)) * b;
            }
        }
    }

    return result;
}


static inline 
bool approximatelyEqual(const double a,
		                const double b,
				        const double epsilon) 
{
	const double fabsa = std::fabs(a);
	const double fabsb = std::fabs(b);
     return std::fabs(a - b) <=
            ((fabsa < fabsb ? fabsb : fabsa) * epsilon);
}

static inline
bool 
almostEqual2sComplement(double A,double B, std::int32_t maxUlps) 
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
void unit_test_bessel_i0_avx_le_15();

void unit_test_bessel_i0_avx_le_15()
{
    using namespace gms::math;
    thread_local std::uniform_real_distribution<double> bessel_i0_rand_in;
    thread_local std::uint64_t seed{};
    constexpr double lo{5.55e-17};
    constexpr double hi{5.0};
    seed = __rdtsc();
    bessel_i0_rand_in = std::uniform_real_distribution<double>(lo,hi);
    auto rand_uni_gen{std::mt19937(seed)};
    const double in_ref_1{bessel_i0_rand_in.operator()(rand_uni_gen)};
    const double in_ref_2{bessel_i0_rand_in.operator()(rand_uni_gen)};
    const double in_ref_3{bessel_i0_rand_in.operator()(rand_uni_gen)};
    const double in_ref_4{bessel_i0_rand_in.operator()(rand_uni_gen)};
    __m256d    in_vref_vec{_mm256_setr_pd(in_ref_1,in_ref_2,in_ref_3,in_ref_4)};
    __m256d vec_actual_vref{_mm256_setzero_pd()};
    const double * __restrict__ p_vref_in = reinterpret_cast<const double*>(&in_vref_vec);
    double       * __restrict__ p_vref_out;
    double sca_expected_ref_1{};
    double sca_expected_ref_2{};
    double sca_expected_ref_3{};
    double sca_expected_ref_4{};
    printf("[UNIT_TEST]: function=%s -- **START**\n\n", __PRETTY_FUNCTION__);
    printf("[UNIT-TEST]: -- Calling reference implementation: bessel_i0(calci0).\n");
    printf("[UNIT-TEST]: Inputs are: ref_1=%.17f,ref_2=%.17f,ref_3=%.17f,ref_4=%.17f\n",in_ref_1,in_ref_2,in_ref_3,in_ref_4);
    sca_expected_ref_1 = bessel_i0(in_ref_1,1);
    sca_expected_ref_2 = bessel_i0(in_ref_2,1);
    sca_expected_ref_3 = bessel_i0(in_ref_3,1);
    sca_expected_ref_4 = bessel_i0(in_ref_4,1);
    printf("[UNIT-TEST]: -- Calling AVX2 implementation: bessel_i0_le15_avx_pd.\n");
    printf("[UNIT-TEST]: Inputs are: ref_1=%.17f,ref_2=%.17f,ref_3=%.17f,ref_4=%.17f\n",p_vref_in[0],p_vref_in[1],p_vref_in[2],p_vref_in[3]);
    //__asm__("int3");
    vec_actual_vref = bessel_i0_le15_avx_pd(in_vref_vec);
    p_vref_out      = reinterpret_cast<double*>(&vec_actual_vref);
    //printf("[UNIT-TEST]: Calculating the results approximately equal by the epsilon:%.17f\n",std::numeric_limits<double>::epsilon());
    const bool ref_1_approx_equality{approximatelyEqual(sca_expected_ref_1,p_vref_out[0],std::numeric_limits<double>::epsilon())};
    const bool ref_2_approx_equality{approximatelyEqual(sca_expected_ref_2,p_vref_out[1],std::numeric_limits<double>::epsilon())};
    const bool ref_3_approx_equality{approximatelyEqual(sca_expected_ref_3,p_vref_out[2],std::numeric_limits<double>::epsilon())};
    const bool ref_4_approx_equality{approximatelyEqual(sca_expected_ref_4,p_vref_out[3],std::numeric_limits<double>::epsilon())};
    if(!ref_1_approx_equality) 
       {std::printf("[UNIT-TEST #1]-- **FAILED**, sca=%.17f,vec=%.17f,delta=%.17f\n",sca_expected_ref_1,p_vref_out[0],sca_expected_ref_1-p_vref_out[0]);}
    if(!ref_2_approx_equality) 
       {std::printf("[UNIT-TEST #2]-- **FAILED**, sca=%.17f,vec=%.17f,delta=%.17f\n",sca_expected_ref_2,p_vref_out[1],sca_expected_ref_2-p_vref_out[1]);}
    if(!ref_3_approx_equality) 
       {std::printf("[UNIT-TEST #1]-- **FAILED**, sca=%.17f,vec=%.17f,delta=%.17f\n",sca_expected_ref_3,p_vref_out[2],sca_expected_ref_3-p_vref_out[2]);}
    if(!ref_4_approx_equality) 
       {std::printf("[UNIT-TEST #2]-- **FAILED**, sca=%.17f,vec=%.17f,delta=%.17f\n",sca_expected_ref_4,p_vref_out[3],sca_expected_ref_4-p_vref_out[3]);}
     printf("[UNIT_TEST]: function=%s -- **END**\n\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_bessel_i0_avx_ge_15();

void unit_test_bessel_i0_avx_ge_15()
{
    using namespace gms::math;
    thread_local std::uniform_real_distribution<double> bessel_i0_rand_in;
    thread_local std::uint64_t seed{};
    constexpr double lo{+15.0};
    constexpr double hi{+100.0};
    seed = __rdtsc();
    bessel_i0_rand_in = std::uniform_real_distribution<double>(lo,hi);
    auto rand_uni_gen{std::mt19937(seed)};
    const double in_ref_1{bessel_i0_rand_in.operator()(rand_uni_gen)};
    const double in_ref_2{bessel_i0_rand_in.operator()(rand_uni_gen)};
    const double in_ref_3{bessel_i0_rand_in.operator()(rand_uni_gen)};
    const double in_ref_4{bessel_i0_rand_in.operator()(rand_uni_gen)};
    __m256d    in_vref_vec{_mm256_setr_pd(in_ref_1,in_ref_2,in_ref_3,in_ref_4)};
    __m256d vec_actual_vref{_mm256_setzero_pd()};
    const double * __restrict__ p_vref_in = reinterpret_cast<const double*>(&in_vref_vec);
    double       * __restrict__ p_vref_out;
    double sca_expected_ref_1{};
    double sca_expected_ref_2{};
    double sca_expected_ref_3{};
    double sca_expected_ref_4{};
    printf("[UNIT_TEST]: function=%s -- **START**\n\n", __PRETTY_FUNCTION__);
    printf("[UNIT-TEST]: -- Calling reference implementation: bessel_i0(calci0).\n");
    printf("[UNIT-TEST]: Inputs are: ref_1=%.17f,ref_2=%.17f,ref_3=%.17f,ref_4=%.17f\n",in_ref_1,in_ref_2,in_ref_3,in_ref_4);
    sca_expected_ref_1 = bessel_i0(in_ref_1,1);
    sca_expected_ref_2 = bessel_i0(in_ref_2,1);
    sca_expected_ref_3 = bessel_i0(in_ref_3,1);
    sca_expected_ref_4 = bessel_i0(in_ref_4,1);
    printf("[UNIT-TEST]: -- Calling AVX2 implementation: bessel_i0_le15_avx_pd.\n");
    printf("[UNIT-TEST]: Inputs are: ref_1=%.17f,ref_2=%.17f,ref_3=%.17f,ref_4=%.17f\n",p_vref_in[0],p_vref_in[1],p_vref_in[2],p_vref_in[3]);
    //__asm__("int3");
    vec_actual_vref = bessel_i0_le15_avx_pd(in_vref_vec);
    p_vref_out      = reinterpret_cast<double*>(&vec_actual_vref);
    //printf("[UNIT-TEST]: Calculating the results approximately equal by the epsilon:%.17f\n",std::numeric_limits<double>::epsilon());
    const bool ref_1_approx_equality{approximatelyEqual(sca_expected_ref_1,p_vref_out[0],std::numeric_limits<double>::epsilon())};
    const bool ref_2_approx_equality{approximatelyEqual(sca_expected_ref_2,p_vref_out[1],std::numeric_limits<double>::epsilon())};
    const bool ref_3_approx_equality{approximatelyEqual(sca_expected_ref_3,p_vref_out[2],std::numeric_limits<double>::epsilon())};
    const bool ref_4_approx_equality{approximatelyEqual(sca_expected_ref_4,p_vref_out[3],std::numeric_limits<double>::epsilon())};
    if(!ref_1_approx_equality) 
       {std::printf("[UNIT-TEST #1]-- **FAILED**, sca=%.17f,vec=%.17f,delta=%.17f\n",sca_expected_ref_1,p_vref_out[0],sca_expected_ref_1-p_vref_out[0]);}
    if(!ref_2_approx_equality) 
       {std::printf("[UNIT-TEST #2]-- **FAILED**, sca=%.17f,vec=%.17f,delta=%.17f\n",sca_expected_ref_2,p_vref_out[1],sca_expected_ref_2-p_vref_out[1]);}
    if(!ref_3_approx_equality) 
       {std::printf("[UNIT-TEST #1]-- **FAILED**, sca=%.17f,vec=%.17f,delta=%.17f\n",sca_expected_ref_3,p_vref_out[2],sca_expected_ref_3-p_vref_out[2]);}
    if(!ref_4_approx_equality) 
       {std::printf("[UNIT-TEST #2]-- **FAILED**, sca=%.17f,vec=%.17f,delta=%.17f\n",sca_expected_ref_4,p_vref_out[3],sca_expected_ref_4-p_vref_out[3]);}
     printf("[UNIT_TEST]: function=%s -- **END**\n\n", __PRETTY_FUNCTION__);
}


int main()
{
    unit_test_bessel_i0_avx_le_15();
    unit_test_bessel_i0_avx_ge_15();
    return 0;
}