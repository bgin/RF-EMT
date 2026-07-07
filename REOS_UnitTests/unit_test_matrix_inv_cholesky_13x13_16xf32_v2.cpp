#include <stdio.h>
#include <cstdlib>
#include <cfloat> // subnormal detection
#include <cmath>  // subnormal detection
#include <ctime>
#include <random>
#include <functional>
#include <immintrin.h> // __rdtsc()
#include "GMS_malloc.h" // for demangle 
#include "GMS_matrix_inv_cholesky.h"
#include "GMS_machine_utils.h"
#include "GMS_triangle_LUT.h"
#include "GMS_ref_sine_waves_LUT.h"

/*
   icpc -o unit_test_matrix_inv_cholesky_13x13_16xf32_v2  -std=c++17 -fp-model fast=2 -ftz -ggdb -ipo -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_machine_utils.h GMS_triangle_LUT.h GMS_ref_sine_waves_LUT.h GMS_ref_sine_waves_LUT.cpp GMS_matrix_inv_cholesky.h GMS_matrix_inv_cholesky.cpp unit_test_matrix_inv_cholesky_13x13_16xf32_v2.cpp 
   ASM: 
   icpc -S -fverbose-asm -masm=intel -std=c++17 -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_machine_utils.h GMS_triangle_LUT.h GMS_ref_sine_waves_LUT.h GMS_ref_sine_waves_LUT.cpp GMS_matrix_inv_cholesky.h GMS_matrix_inv_cholesky.cpp unit_test_matrix_inv_cholesky_13x13_16xf32_v2.cpp 
   
*/
 

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

static inline
bool detect_subnormal(const float x)
{
    constexpr std::int32_t subnormal_code = 3;
    const std::int32_t subnormal_value = std::fpclassify(x);
    bool b_ret{};
    (subnormal_code==subnormal_value) ? b_ret = true : b_ret = false;
    return (b_ret);
}

#define VEC_PS_LEN 16

__attribute__((hot))
__attribute__((noinline))
__attribute__((aligned(32)))
void unit_test_matrix_inv_cholesky_13x13_v2();

void unit_test_matrix_inv_cholesky_13x13_v2()
{
    using namespace gms::math;
    using namespace gms::math::mat_inv_chol_detail;
    __ATTR_ALIGN__(64) __m512 matBRe[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13];
    __ATTR_ALIGN__(64) __m512 matBIm[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13];
    __ATTR_ALIGN__(64) __m512 matInvBRe[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13];
    __ATTR_ALIGN__(64) __m512 matInvBIm[MAT_SQR_SIZE_13][MAT_SQR_SIZE_13];
    constexpr int32_t nelems{2704};
    const float * __restrict ptr_invBRe = nullptr;
    const float * __restrict ptr_invBIm = nullptr;
    const float * __restrict ptr_LUT_sin_1_1      = &LUT_sin_1_1_4096[0];
    const float * __restrict ptr_triang_wave_LUT  = &triangular_wave_LUT[0];   
    std::uint64_t tsc_start{},tsc_end{};
    [[maybe_unused]] std::int32_t printf_ret{};
    for(std::int32_t i = 0;i != MAT_SQR_SIZE_13; ++i)
    {   
        const std::int32_t outer_idx = i*MAT_SQR_SIZE_13;
        for(std::int32_t j = 0;j != MAT_SQR_SIZE_13; ++j)  
        {
            const std::int32_t inner_idx = (outer_idx+j)*VEC_PS_LEN;
            matBRe[i][j] = _mm512_load_ps((const float*)&ptr_LUT_sin_1_1[inner_idx]);
            matBIm[i][j] = _mm512_load_ps((const float*)&ptr_triang_wave_LUT[inner_idx]);
        }
    }
    printf_ret = printf("[UNIT-TEST]: -- START: mat_inv_cholesky_13x13_16xf32_v2<true,true,13>.\n");
    tsc_start = gms::common::rdtscp_fenced();
    gms::math::mat_inv_cholesky_13x13_16xf32_v2<true,true,13>(matBRe,matBIm,
                                                               matInvBRe,matInvBIm);
    tsc_end   = gms::common::rdtscp_fenced();                                                
    ptr_invBRe = reinterpret_cast<const float * __restrict>(&matInvBRe);
    ptr_invBIm = reinterpret_cast<const float * __restrict>(&matInvBIm);
    printf("[UNIT-TEST:] -- Dumping results:\n");
    for(int32_t __i{0}; __i != nelems; ++__i)
    {

        printf_ret =  printf("[UNIT-TEST]: matBRe[%d]=%.7f,matInvBRe[%d]=%.7f\n",__i,ptr_LUT_sin_1_1[__i],__i,ptr_invBRe[__i]);
        printNumber("matBRe",ptr_LUT_sin_1_1[__i],0);
        printNumber("matInvBRe",ptr_invBRe[__i],0);
                     
        printf_ret =  printf("[UNIT-TEST]: matBIm[%d]=%.7f,matInvBIm[%d]=%.7f\n",__i,ptr_triang_wave_LUT[__i],__i,ptr_invBIm[__i]);
        printNumber("matBIm",ptr_triang_wave_LUT[__i],0);
        printNumber("matInvBIm",ptr_invBIm[__i],0);
        if(detect_subnormal(ptr_invBRe[__i])) { std::printf("[UNIT-TEST:] ***SUBNORMAL*** in I-ch -- %.17f\n",ptr_invBRe[__i]);}
        if(detect_subnormal(ptr_invBIm[__i])) { std::printf("[UNIT-TEST:] ***SUBNORMAL*** in Q-ch -- %.17f\n",ptr_invBIm[__i]);}
    } 
    printf_ret = printf("[UNIT-TEST]: Crude-Timing!!: start=%llu,end=%llu,delta=%llu\n",tsc_start,tsc_end,(tsc_end-tsc_start));
    printf_ret = printf("[UNIT-TEST]: -- END: mat_inv_cholesky_13x13_16xf32_v2<true,true,13>.\n");
}            
        

int main()
{   
    unit_test_matrix_inv_cholesky_13x13_v2();
    return 0;
}

