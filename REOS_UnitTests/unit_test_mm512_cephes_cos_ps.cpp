#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <cmath>
#include "GMS_ansi_colors.h"
#include "GMS_cephes_sin_cos.h"
#include "GMS_vec_cephes_cos_ps.h"
#include "GMS_loop_indices_LUT.h"
/*
   icpc -o unit_test_mm512_cephes_cos_ps -O3 -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_ansi_colors.h GMS_cephes_sin_cos.h GMS_vec_cephes_cos_ps.h GMS_loop_indices_LUT.h GMS_loop_indices_huge_LUT.cpp unit_test_mm512_cephes_cos_ps.cpp
   ASM: 
   icpc -S --O3 fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_ansi_colors.h GMS_cephes_sin_cos.h GMS_vec_cephes_cos_ps.h GMS_loop_indices_LUT.h GMS_loop_indices_huge_LUT.cpp unit_test_mm512_cephes_cos_ps.cpp

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

__attribute__((hot))
__attribute__((noinline))
void unit_test_mm512_cephes_cos_ps_trivial();

void unit_test_mm512_cephes_cos_ps_trivial()
{
    using namespace gms::math;
    constexpr float in_ref_1{9.01f};
    constexpr float in_ref_2{9.19f};
    constexpr float in_ref_3{9.458f};
    constexpr float in_ref_4{10.57079632679489661923132169164f};
    constexpr float in_ref_5{10.0f};
    constexpr float in_ref_6{10.000000001f};
    constexpr float in_ref_7{10.999999999f};
    constexpr float in_ref_8{9.14159265358979323846264338328f};
    constexpr float in_ref_9{10.000000009f};
    constexpr float in_ref_10{12.0f};
    constexpr float in_ref_11{1124.458f};
    constexpr float in_ref_12{11254.57079632679489661923132169164f};
    constexpr float in_ref_13{11.0f};
    constexpr float in_ref_14{859.000000001f};
    constexpr float in_ref_15{99.999999999f};
    constexpr float in_ref_16{45.14159265358979323846264338328f};
    const __m512    in_vref_vec{_mm512_setr_ps(in_ref_1,in_ref_2,in_ref_3,in_ref_4,
                                               in_ref_5,in_ref_6,in_ref_7,in_ref_8,
                                               in_ref_9,in_ref_10,in_ref_11,in_ref_12,
                                               in_ref_13,in_ref_14,in_ref_15,in_ref_16)};
    __m512 vec_actual_vref{_mm512_setzero_ps()};
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
    float sca_expected_ref_9{};
    float sca_expected_ref_10{};
    float sca_expected_ref_11{};
    float sca_expected_ref_12{};
    float sca_expected_ref_13{};
    float sca_expected_ref_14{};
    float sca_expected_ref_15{};
    float sca_expected_ref_16{};
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
    printf("[UNIT-TEST]: Inputs (3rd set) are: ref_9=%.7f,ref_10=%.7f,ref_11=%.7f,ref_812=%.7f\n",in_ref_9,in_ref_10,in_ref_11,in_ref_12);
    sca_expected_ref_9  = ceph_cosf(in_ref_9);
    sca_expected_ref_10 = ceph_cosf(in_ref_10);
    sca_expected_ref_11 = ceph_cosf(in_ref_11);
    sca_expected_ref_12 = ceph_cosf(in_ref_12);
    printf("[UNIT-TEST]: Inputs (4th set) are: ref_13=%.7f,ref_14=%.7f,ref_15=%.7f,ref_16=%.7f\n",in_ref_13,in_ref_14,in_ref_15,in_ref_16);
    sca_expected_ref_13 = ceph_cosf(in_ref_13);
    sca_expected_ref_14 = ceph_cosf(in_ref_14);
    sca_expected_ref_15 = ceph_cosf(in_ref_15);
    sca_expected_ref_16 = ceph_cosf(in_ref_16);
    printf("[UNIT-TEST]: -- Calling AVX512 implementation: _mm512_ceph_cos_ps.\n");
    printf("[UNIT-TEST]: Inputs (low1-ZMM) are: ref_1=%.7f,ref_2=%.7f,ref_3=%.7f,ref_4=%.7f\n",p_vref_in[0],p_vref_in[1],p_vref_in[2],p_vref_in[3]);
    printf("[UNIT-TEST]: Inputs (low2-ZMM) are: ref_5=%.7f,ref_6=%.7f,ref_7=%.7f,ref_8=%.7f\n",p_vref_in[4],p_vref_in[5],p_vref_in[6],p_vref_in[7]);
    printf("[UNIT-TEST]: Inputs (high1-ZMM)are: ref_9=%.7f,ref_10=%.7f,ref_11=%.7f,ref_12=%.7f\n",p_vref_in[8],p_vref_in[9],p_vref_in[10],p_vref_in[11]);
    printf("[UNIT-TEST]: Inputs (high2-ZMM)are: ref_13=%.7f,ref_14=%.7f,ref_15=%.7f,ref_16=%.7f\n",p_vref_in[12],p_vref_in[13],p_vref_in[14],p_vref_in[15]);
    //__asm__("int3");
    vec_actual_vref = _mm512_ceph_cosf_ps(in_vref_vec);
    p_vref_out      = reinterpret_cast<float*>(&vec_actual_vref);
    printf("[UNIT-TEST]: Calculating the results approximately equal by the epsilon:%.17f\n",std::numeric_limits<double>::epsilon());
    const bool ref_1_approx_equality{approximatelyEqual(sca_expected_ref_1,p_vref_out[0],std::numeric_limits<double>::epsilon())};
    const bool ref_2_approx_equality{approximatelyEqual(sca_expected_ref_2,p_vref_out[1],std::numeric_limits<double>::epsilon())};
    const bool ref_3_approx_equality{approximatelyEqual(sca_expected_ref_3,p_vref_out[2],std::numeric_limits<double>::epsilon())};
    const bool ref_4_approx_equality{approximatelyEqual(sca_expected_ref_4,p_vref_out[3],std::numeric_limits<double>::epsilon())};
    const bool ref_5_approx_equality{approximatelyEqual(sca_expected_ref_5,p_vref_out[4],std::numeric_limits<double>::epsilon())};
    const bool ref_6_approx_equality{approximatelyEqual(sca_expected_ref_6,p_vref_out[5],std::numeric_limits<double>::epsilon())};
    const bool ref_7_approx_equality{approximatelyEqual(sca_expected_ref_7,p_vref_out[6],std::numeric_limits<double>::epsilon())};
    const bool ref_8_approx_equality{approximatelyEqual(sca_expected_ref_8,p_vref_out[7],std::numeric_limits<double>::epsilon())};
    const bool ref_9_approx_equality{approximatelyEqual(sca_expected_ref_9,p_vref_out[8],std::numeric_limits<double>::epsilon())};
    const bool ref_10_approx_equality{approximatelyEqual(sca_expected_ref_10,p_vref_out[9],std::numeric_limits<double>::epsilon())};
    const bool ref_11_approx_equality{approximatelyEqual(sca_expected_ref_11,p_vref_out[10],std::numeric_limits<double>::epsilon())};
    const bool ref_12_approx_equality{approximatelyEqual(sca_expected_ref_12,p_vref_out[11],std::numeric_limits<double>::epsilon())};
    const bool ref_13_approx_equality{approximatelyEqual(sca_expected_ref_13,p_vref_out[12],std::numeric_limits<double>::epsilon())};
    const bool ref_14_approx_equality{approximatelyEqual(sca_expected_ref_14,p_vref_out[13],std::numeric_limits<double>::epsilon())};
    const bool ref_15_approx_equality{approximatelyEqual(sca_expected_ref_15,p_vref_out[14],std::numeric_limits<double>::epsilon())};
    const bool ref_16_approx_equality{approximatelyEqual(sca_expected_ref_16,p_vref_out[15],std::numeric_limits<double>::epsilon())};
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
    if(!ref_9_approx_equality) 
       {printf(ANSI_COLOR_RED "[UNIT-TEST #9]-- **FAILED**, sca=%.7f,vec=%.7f,delta=%.7f" ANSI_RESET_ALL "\n",sca_expected_ref_9,p_vref_out[8],std::abs(sca_expected_ref_9-p_vref_out[8]));}
    if(!ref_10_approx_equality) 
       {printf(ANSI_COLOR_RED "[UNIT-TEST #10]--**FAILED**, sca=%.7f,vec=%.7f,delta=%.7f" ANSI_RESET_ALL "\n",sca_expected_ref_10,p_vref_out[9],std::abs(sca_expected_ref_10-p_vref_out[9]));}
    if(!ref_11_approx_equality) 
       {printf(ANSI_COLOR_RED "[UNIT-TEST #11]-- **FAILED**, sca=%.7f,vec=%.7f,delta=%.7f" ANSI_RESET_ALL "\n",sca_expected_ref_11,p_vref_out[10],std::abs(sca_expected_ref_11-p_vref_out[10]));}
    if(!ref_12_approx_equality) 
       {printf(ANSI_COLOR_RED "[UNIT-TEST #12]-- **FAILED**, sca=%.7f,vec=%.7f,delta=%.7f" ANSI_RESET_ALL "\n",sca_expected_ref_12,p_vref_out[11],std::abs(sca_expected_ref_12-p_vref_out[11]));}
    if(!ref_13_approx_equality) 
       {printf(ANSI_COLOR_RED "[UNIT-TEST #13]-- **FAILED**, sca=%.7f,vec=%.7f,delta=%.7f" ANSI_RESET_ALL "\n",sca_expected_ref_13,p_vref_out[12],std::abs(sca_expected_ref_13-p_vref_out[12]));}
    if(!ref_14_approx_equality) 
       {printf(ANSI_COLOR_RED "[UNIT-TEST #14]-- **FAILED**, sca=%.7f,vec=%.7f,delta=%.7f" ANSI_RESET_ALL "\n",sca_expected_ref_14,p_vref_out[13],std::abs(sca_expected_ref_14-p_vref_out[13]));}
    if(!ref_15_approx_equality) 
       {printf(ANSI_COLOR_RED "[UNIT-TEST #15]-- **FAILED**, sca=%.7f,vec=%.7f,delta=%.7f" ANSI_RESET_ALL "\n",sca_expected_ref_15,p_vref_out[14],std::abs(sca_expected_ref_15-p_vref_out[14]));}
    if(!ref_16_approx_equality) 
       {printf(ANSI_COLOR_RED "[UNIT-TEST #16]-- **FAILED**, sca=%.7f,vec=%.7f,delta=%.7f" ANSI_RESET_ALL "\n",sca_expected_ref_16,p_vref_out[15],std::abs(sca_expected_ref_16-p_vref_out[15]));}     
     printf("[UNIT_TEST]: function=%s -- **END**\n\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_mm512_cephes_cos_ps_full_range();

void unit_test_mm512_cephes_cos_ps_full_range()
{
     using namespace gms::math;
     std::int32_t found_inequalities{0};
     printf("[UNIT_TEST]: function=%s -- **START**\n\n", __PRETTY_FUNCTION__);
     for(std::int32_t i{0}; i < 62500-16; i += 16) 
     {
                      
            std::int32_t j;
            _mm_prefetch((const char*)&gms::math::LUT_loop_indices_1000000_align64[i],_MM_HINT_T1);
            __m512 vec_cos_arg_in{_mm512_load_ps(&gms::math::LUT_loop_indices_1000000_align64[i])};
             __m512 vec_cos_arg_out{_mm512_ceph_cosf_ps(vec_cos_arg_in)};
            float sca_cos_arg_in[16] = { gms::math::LUT_loop_indices_1000000_align64[i],
                                         gms::math::LUT_loop_indices_1000000_align64[i+1],
                                         gms::math::LUT_loop_indices_1000000_align64[i+2],
                                         gms::math::LUT_loop_indices_1000000_align64[i+3],
                                         gms::math::LUT_loop_indices_1000000_align64[i+4],
                                         gms::math::LUT_loop_indices_1000000_align64[i+5],
                                         gms::math::LUT_loop_indices_1000000_align64[i+6],
                                         gms::math::LUT_loop_indices_1000000_align64[i+7],
                                         gms::math::LUT_loop_indices_1000000_align64[i+8],
                                         gms::math::LUT_loop_indices_1000000_align64[i+9],
                                         gms::math::LUT_loop_indices_1000000_align64[i+10],
                                         gms::math::LUT_loop_indices_1000000_align64[i+11],
                                         gms::math::LUT_loop_indices_1000000_align64[i+12],
                                         gms::math::LUT_loop_indices_1000000_align64[i+13],
                                         gms::math::LUT_loop_indices_1000000_align64[i+14],
                                         gms::math::LUT_loop_indices_1000000_align64[i+15]};
            float sca_cos_arg_out[16] = {};
            float * __restrict__ p_vec_cos_arg_out{reinterpret_cast<float*>(&vec_cos_arg_out)};
            for(j = 0; j != 16; ++j)  
            {
                 float cos_arg_in{sca_cos_arg_in[j]};
                 sca_cos_arg_out[j] = ceph_cosf(cos_arg_in);
                 const bool compare_approx_equal{approximatelyEqual(sca_cos_arg_out[j],p_vec_cos_arg_out[j],std::numeric_limits<float>::epsilon())};
                 const float sca_res{sca_cos_arg_out[j]};
                 const float vec_res{p_vec_cos_arg_out[j]};
                 if(!compare_approx_equal) 
                 { 
                   ++found_inequalities;
                   printf(ANSI_COLOR_RED "[UNIT-TEST]-- **FAILED**,iter=%d,input=%.7f,sca=%.7f,vec=%.7f,delta=%.7f" ANSI_RESET_ALL "\n",i,cos_arg_in,sca_res,vec_res,sca_res-vec_res);
                 }
            }            
     }
     std::printf("[UNIT-TEST] -- Found=%d inequalities\n",found_inequalities);
     std::printf("[UNIT_TEST]: function=%s -- **END**\n\n", __PRETTY_FUNCTION__);
}



int main()
{
    //unit_test_mm512_cephes_cos_ps_trivial();
    unit_test_mm512_cephes_cos_ps_full_range();
    return 0;
}