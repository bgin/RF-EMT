#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <cmath>
#include "GMS_cephes_sin_cos.h"

/*
   icpc -o unit_test_mm_cephes_cos_ps -O3 -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h   GMS_cephes_sin_cos.h unit_test_mm_cephes_cos_ps.cpp
   ASM: 
   icpc -S --O3 fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h   GMS_cephes_sin_cos.h GMS_cephes_sin_cos.h unit_test_mm_cephes_cos_ps.cpp

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
void unit_test_mm_cephes_cos_ps_trivial();

void unit_test_mm_cephes_cos_ps_trivial()
{
    using namespace gms::math;
    constexpr float in_ref_1{0.01f};
    constexpr float in_ref_2{0.19f};
    constexpr float in_ref_3{0.458f};
    constexpr float in_ref_4{1.57079632679489661923132169164f};
    const __m128    in_vref_vec{_mm_setr_ps(in_ref_1,in_ref_2,in_ref_3,in_ref_4)};
    __m128 vec_actual_vref{_mm_setzero_ps()};
    const float * __restrict__ p_vref_in = reinterpret_cast<const float*>(&in_vref_vec);
    float       * __restrict__ p_vref_out;
    float sca_expected_ref_1{};
    float sca_expected_ref_2{};
    float sca_expected_ref_3{};
    float sca_expected_ref_4{};
    printf("[UNIT_TEST]: function=%s -- **START**\n\n", __PRETTY_FUNCTION__);
    printf("[UNIT-TEST]: -- Calling reference implementation: ceph_cosf.\n");
    printf("[UNIT-TEST]: Inputs are: ref_1=%.7f,ref_2=%.7f,ref_3=%.7f,ref_4=%.7f\n",in_ref_1,in_ref_2,in_ref_3,in_ref_4);
    sca_expected_ref_1 = ceph_cosf(in_ref_1);
    sca_expected_ref_2 = ceph_cosf(in_ref_2);
    sca_expected_ref_3 = ceph_cosf(in_ref_3);
    sca_expected_ref_4 = ceph_cosf(in_ref_4);
    printf("[UNIT-TEST]: -- Calling SSE implementation: _mm_ceph_cos_ps.\n");
    printf("[UNIT-TEST]: Inputs are: ref_1=%.7f,ref_2=%.7f,ref_3=%.7f,ref_4=%.7f\n",p_vref_in[0],p_vref_in[1],p_vref_in[2],p_vref_in[3]);
    //__asm__("int3");
    vec_actual_vref = _mm_ceph_cosf_ps(in_vref_vec);
    p_vref_out      = reinterpret_cast<float*>(&vec_actual_vref);
    printf("[UNIT-TEST]: Calculating the results approximately equal by the epsilon:%.17f\n",std::numeric_limits<double>::epsilon());
    const bool ref_1_approx_equality{approximatelyEqual(sca_expected_ref_1,p_vref_out[0],std::numeric_limits<double>::epsilon())};
    const bool ref_2_approx_equality{approximatelyEqual(sca_expected_ref_2,p_vref_out[1],std::numeric_limits<double>::epsilon())};
    const bool ref_3_approx_equality{approximatelyEqual(sca_expected_ref_3,p_vref_out[2],std::numeric_limits<double>::epsilon())};
    const bool ref_4_approx_equality{approximatelyEqual(sca_expected_ref_4,p_vref_out[3],std::numeric_limits<double>::epsilon())};
    if(!ref_1_approx_equality) 
       {std::printf("[UNIT-TEST #1]-- **FAILED**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_1,p_vref_out[0],sca_expected_ref_1-p_vref_out[0]);}
    if(!ref_2_approx_equality) 
       {std::printf("[UNIT-TEST #2]-- **FAILED**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_2,p_vref_out[1],sca_expected_ref_2-p_vref_out[1]);}
    if(!ref_3_approx_equality) 
       {std::printf("[UNIT-TEST #3]-- **FAILED**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_3,p_vref_out[2],sca_expected_ref_3-p_vref_out[2]);}
    if(!ref_4_approx_equality) 
       {std::printf("[UNIT-TEST #4]-- **FAILED**, sca=%.7f,vec=%.7f,delta=%.7f\n",sca_expected_ref_4,p_vref_out[3],sca_expected_ref_4-p_vref_out[3]);}
     printf("[UNIT_TEST]: function=%s -- **END**\n\n", __PRETTY_FUNCTION__);
}



int main()
{
    unit_test_mm_cephes_cos_ps_trivial();
    return 0;
}