#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <immintrin.h> // __rdtsc()
#include "GMS_malloc.h" // for demangle 
#include "GMS_matrix_inv_cholesky.h"
#include "unit_test_mat_inv_cholesky_defines.h"
#include "GMS_machine_utils.h"
/*
   icpc -o unit_test_matrix_inv_cholesky_11x11_validity -std=c++17 -fp-model fast=2 -ftz -ggdb -ipo -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_machine_utils.h GMS_matrix_inv_cholesky.h GMS_matrix_inv_cholesky.cpp unit_test_mat_inv_cholesky_defines.h unit_test_inv_mat_cholesky_11x11_validity.cpp 
   ASM: 
   icpc -S -fverbose-asm -masm=intel -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_config.h GMS_malloc.h GMS_machine_utils.h GMS_matrix_inv_cholesky.h GMS_matrix_inv_cholesky.cpp unit_test_mat_inv_cholesky_defines.h unit_test_inv_mat_cholesky_11x11_validity.cpp 
*/
 
static __m512 constZero = _mm512_setzero_ps();

void matrix_inv_cholesky_11x11(__m512 matBRe[N_11][N_11], __m512 matBIm[N_11][N_11],
    __m512 matInvBRe[N_11][N_11], __m512 matInvBIm[N_11][N_11])
{
    // temp matrix and variables for matrix inversion
    __m512 matGRe[N_11][N_11], matGIm[N_11][N_11];
    __m512 matLRe[N_11][N_11], matLIm[N_11][N_11];
    __m512 matD[N_11], matND[N_11];
    __m512 temp0, temp1, temp2;
    int32_t i, j, k;

    /////////////////////////////////// get G, B = G*G', G is a lower triangular matrix
    // Column 0
    GET_G00_REF(matGRe, matBRe, matD, matND);
    GET_G_COL0_REF(matGRe, matGIm, matBRe, matBIm, matD, 1);
    GET_G_COL0_REF(matGRe, matGIm, matBRe, matBIm, matD, 2);
    GET_G_COL0_REF(matGRe, matGIm, matBRe, matBIm, matD, 3);
    GET_G_COL0_REF(matGRe, matGIm, matBRe, matBIm, matD, 4);
    GET_G_COL0_REF(matGRe, matGIm, matBRe, matBIm, matD, 5);
    GET_G_COL0_REF(matGRe, matGIm, matBRe, matBIm, matD, 6);
    GET_G_COL0_REF(matGRe, matGIm, matBRe, matBIm, matD, 7);
    GET_G_COL0_REF(matGRe, matGIm, matBRe, matBIm, matD, 8);
    GET_G_COL0_REF(matGRe, matGIm, matBRe, matBIm, matD, 9);
    GET_G_COL0_REF(matGRe, matGIm, matBRe, matBIm, matD, 10);

    // Column 1
    GET_G11_REF(matGRe, matGIm, matBRe, matD, matND, temp0);
    GET_G_COL1_REF(matGRe, matGIm, matBRe, matBIm, matD, 2, temp0, temp1);
    GET_G_COL1_REF(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);
    GET_G_COL1_REF(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL1_REF(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL1_REF(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL1_REF(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL1_REF(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);
    GET_G_COL1_REF(matGRe, matGIm, matBRe, matBIm, matD, 9, temp0, temp1);
    GET_G_COL1_REF(matGRe, matGIm, matBRe, matBIm, matD, 10, temp0, temp1);

    // Column 2
    GET_G22_REF(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    GET_G_COL2_REF(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);
    GET_G_COL2_REF(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL2_REF(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL2_REF(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL2_REF(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL2_REF(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);
    GET_G_COL2_REF(matGRe, matGIm, matBRe, matBIm, matD, 9, temp0, temp1);
    GET_G_COL2_REF(matGRe, matGIm, matBRe, matBIm, matD, 10, temp0, temp1);

    // Column 3
    GET_G33_REF(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    GET_G_COL3_REF(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL3_REF(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL3_REF(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL3_REF(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL3_REF(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);
    GET_G_COL3_REF(matGRe, matGIm, matBRe, matBIm, matD, 9, temp0, temp1);
    GET_G_COL3_REF(matGRe, matGIm, matBRe, matBIm, matD, 10, temp0, temp1);

    // Column 4
    GET_G44_REF(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 5, 4, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 6, 4, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 7, 4, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 8, 4, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 9, 4, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 10, 4, temp0, temp1);

    // Column 5
    GET_G55_REF(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 6, 5, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 7, 5, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 8, 5, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 9, 5, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 10, 5, temp0, temp1);

    // Column 6
    GET_G66_REF(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 7, 6, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 8, 6, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 9, 6, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 10, 6, temp0, temp1);

    // Column 7
    GET_G77_REF(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 8, 7, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 9, 7, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 10, 7, temp0, temp1);

    // Column 8
    GET_Gii_EVEN_REF(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 8);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 9, 8, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 10, 8, temp0, temp1);

    // Column 9
    GET_Gii_ODD_REF(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 9);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 10, 9, temp0, temp1);

    // Column 10
    GET_Gii_EVEN_REF(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 10);

    /////////////////////////////////// get L, L=invG
    // Column 0
    SET_Lii_REF(matLRe, matLIm, matD, 0);
    GET_L_i1i_REF(matLRe, matLIm, matGRe, matGIm, matND, 1, 0);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 2, 0, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 3, 0, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 4, 0, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 5, 0, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 6, 0, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 7, 0, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 8, 0, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 9, 0, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 10, 0, temp0, temp1);

    // Column 1
    SET_Lii_REF(matLRe, matLIm, matD, 1);
    GET_L_i1i_REF(matLRe, matLIm, matGRe, matGIm, matND, 2, 1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 3, 1, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 4, 1, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 5, 1, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 6, 1, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 7, 1, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 8, 1, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 9, 1, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 10, 1, temp0, temp1);

    // Column 2
    SET_Lii_REF(matLRe, matLIm, matD, 2);
    GET_L_i1i_REF(matLRe, matLIm, matGRe, matGIm, matND, 3, 2);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 4, 2, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 5, 2, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 6, 2, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 7, 2, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 8, 2, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 9, 2, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 10, 2, temp0, temp1);

    // Column 3
    SET_Lii_REF(matLRe, matLIm, matD, 3);
    GET_L_i1i_REF(matLRe, matLIm, matGRe, matGIm, matND, 4, 3);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 5, 3, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 6, 3, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 7, 3, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 8, 3, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 9, 3, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 10, 3, temp0, temp1);

    // Column 4
    SET_Lii_REF(matLRe, matLIm, matD, 4);
    GET_L_i1i_REF(matLRe, matLIm, matGRe, matGIm, matND, 5, 4);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 6, 4, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 7, 4, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 8, 4, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 9, 4, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 10, 4, temp0, temp1);

    // Column 5
    SET_Lii_REF(matLRe, matLIm, matD, 5);
    GET_L_i1i_REF(matLRe, matLIm, matGRe, matGIm, matND, 6, 5);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 7, 5, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 8, 5, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 9, 5, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 10, 5, temp0, temp1);

    // Column 6
    SET_Lii_REF(matLRe, matLIm, matD, 6);
    GET_L_i1i_REF(matLRe, matLIm, matGRe, matGIm, matND, 7, 6);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 8, 6, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 9, 6, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 10, 6, temp0, temp1);

    // Column 7
    SET_Lii_REF(matLRe, matLIm, matD, 7);
    GET_L_i1i_REF(matLRe, matLIm, matGRe, matGIm, matND, 8, 7);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 9, 7, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 10, 7, temp0, temp1);

    // Column 8
    SET_Lii_REF(matLRe, matLIm, matD, 8);
    GET_L_i1i_REF(matLRe, matLIm, matGRe, matGIm, matND, 9, 8);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 10, 8, temp0, temp1);

    // Column 9
    SET_Lii_REF(matLRe, matLIm, matD, 9);
    GET_L_i1i_REF(matLRe, matLIm, matGRe, matGIm, matND, 10, 9);

    // Column 10
    SET_Lii_REF(matLRe, matLIm, matD, 10);

    /////////////////////////////////// get invB = L'*L
    for(i = 0; i < N_11; i ++)
    {
        matInvBRe[i][i] = _mm512_mul_ps(matLRe[i][i], matLRe[i][i]);
        for (k = (i+1); k < N_11; k++)
        {
            temp1 = _mm512_add_ps(_mm512_mul_ps(matLRe[k][i], matLRe[k][i]), _mm512_mul_ps(matLIm[k][i], matLIm[k][i]));
            matInvBRe[i][i] = _mm512_add_ps(matInvBRe[i][i], temp1);
        }
        matInvBIm[i][i] = _mm512_setzero_ps();
    }

    for(i = 0; i < N_11; i ++)
    {
        for(j = i+1; j < N_11; j ++)
        {
            matInvBRe[i][j] = _mm512_mul_ps(matLRe[j][i], matLRe[j][j]);
            matInvBIm[i][j] = _mm512_sub_ps(constZero, _mm512_mul_ps(matLIm[j][i], matLRe[j][j]));

            for (k = (j+1); k < N_11; k++)
            {
                GET_AxBH_REF(matLRe[k][j], matLIm[k][j], matLRe[k][i], matLIm[k][i], temp1, temp2);
                matInvBRe[i][j] = _mm512_add_ps(matInvBRe[i][j], temp1);
                matInvBIm[i][j] = _mm512_add_ps(matInvBIm[i][j], temp2);
            }

            // Hermite matrix
            matInvBRe[j][i] = matInvBRe[i][j];
            matInvBIm[j][i] = _mm512_sub_ps(constZero, matInvBIm[i][j]);
        }
    }
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

#define PRINT_CRUDE_TIMING_RESULTS\
  do\
  {\
     printf_ret = printf("[UNIT-TEST]: Crude-Timing!!: [RF-EMT]  -- start=%llu,end=%llu,delta=%llu\n",tsc_start_rf_emt,tsc_end_rf_emt,(tsc_end_rf_emt-tsc_start_rf_emt));\
     printf_ret = printf("[UNIT-TEST]: Crude-Timing!!: [FlexRAN] -- start=%llu,end=%llu,delta=%llu\n",tsc_start_flexr,tsc_end_flexr,(tsc_end_flexr-tsc_start_flexr));\
  }while(0);\

__attribute__((hot))
__attribute__((noinline))
__attribute__((aligned(32)))
void unit_test_matrix_inv_cholesky_11x11_validity();

void unit_test_matrix_inv_cholesky_11x11_validity()
{
    using namespace gms::math;
    constexpr int32_t nelems{1936};
    __ATTR_ALIGN__(64) float buf_init_re[nelems];
    __ATTR_ALIGN__(64) float buf_init_im[nelems];
    const float * __restrict ptr_invBRe     = nullptr;
    const float * __restrict ptr_invBIm     = nullptr;
    const float * __restrict ptr_invBReRef  = nullptr;
    const float * __restrict ptr_invBImRef  = nullptr;
    unsigned long long seed_re{0ull};
    unsigned long long seed_im{0ull};
    std::uint64_t tsc_start_rf_emt{}, tsc_end_rf_emt{};
    std::uint64_t tsc_start_flexr{},  tsc_end_flexr{};
    std::clock_t seed_d{0ull};
    int32_t which(-1);
    [[maybe_unused]] std::int32_t printf_ret{-9999};
    printf("[UNIT-TEST]: -- START: (Impl:RF-EMT) mat_inv_cholesky_11x11_16xf32<true,true> vs. (Impl:FlexRAN) matrix_inv_cholesky_11x11 (Validity).\n");
    seed_d = std::clock();
    auto rand_d{std::bind(std::uniform_int_distribution<int32_t>(0,2),std::mt19937(seed_d))};
    which = rand_d();
    switch (which)
    {
        case 0 : 
        {
             std::uniform_real_distribution<float> uniform_distro_re;
             char * distro_name_re{NULL};
             int32_t status{9999};
             distro_name_re = gms::common::demangle(typeid(uniform_distro_re).name(),status);
             if(distro_name_re != NULL && status == 0)
             {
                   printf_ret = printf("[UNIT-TEST]: Instantiated distribution of type: %s\n", distro_name_re);
                   gms::common::gms_mm_free(distro_name_re);
             }
             else
             {
                   printf_ret = printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n", typeid(uniform_distro_re).name());
                   if(distro_name_re != NULL) gms::common::gms_mm_free(distro_name_re);
             }
             std::uniform_real_distribution<float> uniform_distro_im;
             char * distro_name_im{NULL};
             status = 9999;
             distro_name_im = gms::common::demangle(typeid(uniform_distro_im).name(),status);
             if(distro_name_im!= NULL && status == 0)
             {
                   printf_ret = printf("[UNIT-TEST]: Instantiated distribution of type: %s\n", distro_name_im);
                   gms::common::gms_mm_free(distro_name_im);
             }
             else
             {
                   printf_ret = printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n", typeid(uniform_distro_im).name());
                   if(distro_name_re != NULL) gms::common::gms_mm_free(distro_name_im);
             }
            __ATTR_ALIGN__(64) __m512 matBRe[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11];
            __ATTR_ALIGN__(64) __m512 matBIm[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11];
            __ATTR_ALIGN__(64) __m512 matInvBRe[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11];
            __ATTR_ALIGN__(64) __m512 matInvBIm[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11];
            __ATTR_ALIGN__(64) __m512 matInvBReRef[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11];
            __ATTR_ALIGN__(64) __m512 matInvBImRef[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11];
             seed_re = __rdtsc();
             auto rand_re{std::bind(std::uniform_real_distribution<float>(0.001f,5.0f),std::mt19937(seed_re))};
             seed_im = __rdtsc();
             auto rand_im{std::bind(std::uniform_real_distribution<float>(0.001f,5.0f),std::mt19937(seed_im))};
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                 float rd_re = rand_re();
                 buf_init_re[__i] = rd_re;
                 float rd_im = rand_im();
                 buf_init_im[__i] = rd_im;
             }
            for(std::int32_t i = 0;i != MAT_SQR_SIZE_11; ++i)
            {   
                const std::int32_t outer_idx = i*MAT_SQR_SIZE_11;
                for(std::int32_t j = 0;j != MAT_SQR_SIZE_11; ++j)  
                {
                    const std::int32_t inner_idx = (outer_idx+j)*VEC_PS_LEN;
                    register __m512 init_vre = _mm512_load_ps((const float*)&buf_init_re[inner_idx]);                    
                    matBRe[i][j]             = init_vre;
                    matInvBReRef[i][j]       = init_vre;
                    register __m512 init_vim = _mm512_load_ps((const float*)&buf_init_im[inner_idx]);
                    matBIm[i][j]             = init_vim;
                    matInvBImRef[i][j]       = init_vim;
                }
            }
            printf_ret = printf("[UNIT-TEST]: -- START: mat_inv_cholesky_11x11_16xf32<true,true>.\n");
            tsc_start_rf_emt = gms::common::rdtscp_fenced();
            gms::math::mat_inv_cholesky_11x11_16xf32<true,true>(matBRe,matBIm,
                                                               matInvBRe,matInvBIm); 
            tsc_end_rf_emt   = gms::common::rdtscp_fenced();                                                 
            ptr_invBRe = reinterpret_cast<const float * __restrict>(&matInvBRe);
            ptr_invBIm = reinterpret_cast<const float * __restrict>(&matInvBIm);
            printf_ret = printf("[UNIT-TEST]: -- START: matrix_inv_cholesky_11x11(***FlexRAN(Reference)***).\n");
            tsc_start_flexr  = gms::common::rdtscp_fenced();
            matrix_inv_cholesky_11x11(matBRe,matBIm,matInvBReRef,matInvBImRef);
            tsc_end_flexr    = gms::common::rdtscp_fenced();
            ptr_invBReRef = reinterpret_cast<const float * __restrict>(&matInvBReRef);
            ptr_invBImRef = reinterpret_cast<const float * __restrict>(&matInvBImRef);
            printf_ret = printf("[UNIT-TEST:] -- Dumping results:\n");
            volatile std::int32_t bre_err_count = 0;
            volatile std::int32_t bim_err_count = 0;
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                
                const float bre     = ptr_invBRe[__i];
                const float bre_ref = ptr_invBReRef[__i];
                const float bim     = ptr_invBIm[__i];
                const float bim_ref = ptr_invBImRef[__i];
                const bool  bre_cmp = almostEqual2sComplement(bre,bre_ref,7);
                const bool  bim_cmp = almostEqual2sComplement(bim,bim_ref,7);
                if(!bre_cmp)
                {
                    ++bre_err_count;
                    [[maybe_unused]] std::int32_t 
                    printf_res_re =  printf("[UNIT-TEST]: ***FAILED*** -- (RF-EMT)bre=%.7f,(FlexRAN)bre_ref=%.7f,fail-count=%d\n",bre,bre_ref,bre_err_count);
                    printNumber("(RF-EMT)bre",bre,0);
                    printNumber("(FlexRAN)bre_ref",bre_ref,0);
                }   
                if(!bim_cmp)
                { 
                    ++bim_err_count;
                    [[maybe_unused]] std::int32_t 
                    printf_res_im =  printf("[UNIT-TEST]:  ***FAILED*** -- (RF-EMT)bim=%.7f,(FlexRAN)bim_ref=%.7f,fail-count=%d\n",bim,bim_ref,bim_err_count);
                    printNumber("(RF-EMT)bim",bim,0);
                    printNumber("(FlexRAN)bim_ref",bim_ref,0);
                } 
                if(__builtin_expect(detect_subnormal(bre),0))     { std::printf("[UNIT-TEST:] -- (RF-EMT)  ***SUBNORMAL*** in I-ch -- %.17f\n",bre);}
                if(__builtin_expect(detect_subnormal(bim),0))     { std::printf("[UNIT-TEST:] -- (RF-EMT)  ***SUBNORMAL*** in Q-ch -- %.17f\n",bim);}
                if(__builtin_expect(detect_subnormal(bre_ref),0)) { std::printf("[UNIT-TEST:] -- (FlexRAN) ***SUBNORMAL*** in I-ch -- %.17f\n",bre_ref);}
                if(__builtin_expect(detect_subnormal(bim_ref),0)) { std::printf("[UNIT-TEST:] -- (FlexRAN) ***SUBNORMAL*** in Q-ch -- %.17f\n",bim_ref);}
               
            } 
            if(!bre_err_count) {printf_ret = printf("[UNIT_TEST:] -- Validity Test: (I-channel): ***PASSED***\n");}
            if(!bim_err_count) {printf_ret = printf("[UNIT_TEST:] -- Validity Test: (Q-channel): ***PASSED***\n");}
        }          
            //printf_ret = printf("[UNIT-TEST]: Crude-Timing!!: [RF-EMT]  -- start=%llu,end=%llu,delta=%llu\n",tsc_start_rf_emt,tsc_end_rf_emt,(tsc_end_rf_emt-tsc_start_rf_emt));  
            //printf_ret = printf("[UNIT-TEST]: Crude-Timing!!: [FlexRAN] -- start=%llu,end=%llu,delta=%llu\n",tsc_start_flexr,tsc_end_flexr,(tsc_end_flexr-tsc_start_flexr));
            PRINT_CRUDE_TIMING_RESULTS;
            printf_ret = printf("[UNIT-TEST]: -- END: mat_inv_cholesky_11x11_16xf32<true,true> vs matrix_inv_cholesky_9x9 (Validity).\n");
        
        break;
        case 1 : 
        {
            std::normal_distribution<float> normal_distro_re;
            char * distro_name_re{NULL};
            int32_t status{9999};
            distro_name_re = gms::common::demangle(typeid(normal_distro_re).name(),status);
            if(distro_name_re != NULL && status == 0)
            {
                   printf_ret = printf("[UNIT-TEST]: Instantiated distribution of type: %s\n", distro_name_re);
                   gms::common::gms_mm_free(distro_name_re);
            }
            else
            {
                   printf_ret = printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n", typeid(normal_distro_re).name());
                   if(distro_name_re != NULL) gms::common::gms_mm_free(distro_name_re);
            }
            std::normal_distribution<float> normal_distro_im;
            char * distro_name_im{NULL};
            status = 9999;
            distro_name_im = gms::common::demangle(typeid(normal_distro_im).name(),status);
            if(distro_name_im != NULL && status == 0)
            {
                   printf_ret = printf("[UNIT-TEST]: Instantiated distribution of type: %s\n", distro_name_im);
                   gms::common::gms_mm_free(distro_name_im);
            }
            else
            {
                   printf_ret = printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n", typeid(normal_distro_im).name());
                   if(distro_name_re != NULL) gms::common::gms_mm_free(distro_name_im);
            }
            __ATTR_ALIGN__(64) __m512 matBRe[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11];
            __ATTR_ALIGN__(64) __m512 matBIm[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11];
            __ATTR_ALIGN__(64) __m512 matInvBRe[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11];
            __ATTR_ALIGN__(64) __m512 matInvBIm[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11];
            __ATTR_ALIGN__(64) __m512 matInvBReRef[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11];
            __ATTR_ALIGN__(64) __m512 matInvBImRef[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11];
             seed_re = __rdtsc();
             auto rand_re{std::bind(std::normal_distribution<float>(0.0f,1.0f),std::mt19937(seed_re))};
             seed_im = __rdtsc();
             auto rand_im{std::bind(std::normal_distribution<float>(1.0f,2.54f),std::mt19937(seed_im))};
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                 float rd_re = rand_re();
                 if(rd_re<0.0f)  rd_re = -1.0f*rd_re;
                 if(rd_re==0.0f) rd_re = 0.1f;
                 buf_init_re[__i] = rd_re;
                 float rd_im = rand_im();
                 if(rd_im<0.0f)  rd_im = -1.0f*rd_im;
                 if(rd_im==0.0f) rd_im = 0.1f;
                 buf_init_im[__i] = rd_im;
             }
            for(std::int32_t i = 0;i != MAT_SQR_SIZE_11; ++i)
            {   
                const std::int32_t outer_idx = i*MAT_SQR_SIZE_11;
                for(std::int32_t j = 0;j != MAT_SQR_SIZE_11; ++j)  
                {
                    const std::int32_t inner_idx = (outer_idx+j)*VEC_PS_LEN;
                    register __m512 init_vre = _mm512_load_ps((const float*)&buf_init_re[inner_idx]);                    
                    matBRe[i][j]             = init_vre;
                    matInvBReRef[i][j]       = init_vre;
                    register __m512 init_vim = _mm512_load_ps((const float*)&buf_init_im[inner_idx]);
                    matBIm[i][j]             = init_vim;
                    matInvBImRef[i][j]       = init_vim;
                }
            }
            printf_ret = printf("[UNIT-TEST]: -- START: mat_inv_cholesky_11x11_16xf32<true,true> ***RF-EMT*** .\n");
            tsc_start_rf_emt  = gms::common::rdtscp_fenced();
            gms::math::mat_inv_cholesky_11x11_16xf32<true,true>(matBRe,matBIm,
                                                               matInvBRe,matInvBIm);
            tsc_end_rf_emt    = gms::common::rdtscp_fenced();
            ptr_invBRe = reinterpret_cast<const float * __restrict>(&matInvBRe);
            ptr_invBIm = reinterpret_cast<const float * __restrict>(&matInvBIm);
            printf_ret = printf("[UNIT-TEST]: -- START: matrix_inv_cholesky_11x11 ***FlexRAN(Reference)***.\n");
            tsc_start_flexr   = gms::common::rdtscp_fenced();
            matrix_inv_cholesky_11x11(matBRe,matBIm,matInvBReRef,matInvBImRef);
            tsc_end_flexr     = gms::common::rdtscp_fenced();
            ptr_invBReRef = reinterpret_cast<const float * __restrict>(&matInvBReRef);
            ptr_invBImRef = reinterpret_cast<const float * __restrict>(&matInvBImRef);
            printf_ret = printf("[UNIT-TEST:] -- Dumping results:\n");
            volatile std::int32_t bre_err_count = 0;
            volatile std::int32_t bim_err_count = 0;
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                const float bre     = ptr_invBRe[__i];
                const float bre_ref = ptr_invBReRef[__i];
                const float bim     = ptr_invBIm[__i];
                const float bim_ref = ptr_invBImRef[__i];
                const bool  bre_cmp = almostEqual2sComplement(bre,bre_ref,7);
                const bool  bim_cmp = almostEqual2sComplement(bim,bim_ref,7);
                if(!bre_cmp)
                {   
                    ++bre_err_count;
                    [[maybe_unused]] std::int32_t 
                    printf_res_re =  printf("[UNIT-TEST]: ***FAILED*** -- (RF-EMT)bre=%.7f,(FlexRAN)bre_ref=%.7f,fail-count=%d\n",bre,bre_ref,bre_err_count);
                    printNumber("(RF-EMT)bre",bre,0);
                    printNumber("(FlexRAN)bre_ref",bre_ref,0);
                }   
                if(!bim_cmp)
                {   
                    ++bim_err_count;
                    [[maybe_unused]] std::int32_t 
                    printf_res_im =  printf("[UNIT-TEST]:  ***FAILED*** -- (RF-EMT)bim=%.7f,(FlexRAN)bim_ref=%.7f,fail-count=%d\n",bim,bim_ref,bim_err_count);
                    printNumber("(RF-EMT)bim",bim,0);
                    printNumber("(FlexRAN)bim_ref",bim_ref,0);
                } 
                if(__builtin_expect(detect_subnormal(bre),0))     { std::printf("[UNIT-TEST:] -- (RF-EMT)  ***SUBNORMAL*** in I-ch -- %.17f\n",bre);}
                if(__builtin_expect(detect_subnormal(bim),0))     { std::printf("[UNIT-TEST:] -- (RF-EMT)  ***SUBNORMAL*** in Q-ch -- %.17f\n",bim);}
                if(__builtin_expect(detect_subnormal(bre_ref),0)) { std::printf("[UNIT-TEST:] -- (FlexRAN) ***SUBNORMAL*** in I-ch -- %.17f\n",bre_ref);}
                if(__builtin_expect(detect_subnormal(bim_ref),0)) { std::printf("[UNIT-TEST:] -- (FlexRAN) ***SUBNORMAL*** in Q-ch -- %.17f\n",bim_ref);}
            } 
            if(!bre_err_count) {printf("[UNIT_TEST:] -- Validity Test: (I-channel): ***PASSED***\n");}
            if(!bim_err_count) {printf("[UNIT_TEST:] -- Validity Test: (Q-channel): ***PASSED***\n");}
            PRINT_CRUDE_TIMING_RESULTS;
            printf_ret = printf("[UNIT-TEST]: -- END: mat_inv_cholesky_11x11_16xf32<true,true> vs matrix_inv_cholesky_11x11 (Validity).\n");
        }
        break;
        case 2 : 
        {
            __ATTR_ALIGN__(64) __m512 matBRe[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11];
            __ATTR_ALIGN__(64) __m512 matBIm[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11];
            __ATTR_ALIGN__(64) __m512 matInvBRe[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11];
            __ATTR_ALIGN__(64) __m512 matInvBIm[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11];
            __ATTR_ALIGN__(64) __m512 matInvBReRef[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11];
            __ATTR_ALIGN__(64) __m512 matInvBImRef[MAT_SQR_SIZE_11][MAT_SQR_SIZE_11];
            printf_ret = printf("[UNIT-TEST:] -- ***DETERMINISTIC-INPUT-FILL***\n");
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                 float f_i = static_cast<float>(__i);
                 buf_init_re[__i] = buf_init_im[__i] = f_i;
             }
            for(std::int32_t i = 0;i != MAT_SQR_SIZE_11; ++i)
            {   
                const std::int32_t outer_idx = i*MAT_SQR_SIZE_11;
                for(std::int32_t j = 0;j != MAT_SQR_SIZE_11; ++j)  
                {
                    const std::int32_t inner_idx = (outer_idx+j)*VEC_PS_LEN;
                    register __m512 init_vre = _mm512_load_ps((const float*)&buf_init_re[inner_idx]);                    
                    matBRe[i][j]             = init_vre;
                    matInvBReRef[i][j]       = init_vre;
                    register __m512 init_vim = _mm512_load_ps((const float*)&buf_init_im[inner_idx]);
                    matBIm[i][j]             = init_vim;
                    matInvBImRef[i][j]       = init_vim;
                }
            }
            printf_ret = printf("[UNIT-TEST]: -- START: mat_inv_cholesky_11x11_16xf32<true,true> ***RF-EMT*** .\n");
            tsc_start_rf_emt = gms::common::rdtscp_fenced();
            gms::math::mat_inv_cholesky_11x11_16xf32<true,true>(matBRe,matBIm,
                                                               matInvBRe,matInvBIm);
            tsc_end_rf_emt   = gms::common::rdtscp_fenced();
            ptr_invBRe = reinterpret_cast<const float * __restrict>(&matInvBRe);
            ptr_invBIm = reinterpret_cast<const float * __restrict>(&matInvBIm);
            printf_ret = printf("[UNIT-TEST]: -- START: matrix_inv_cholesky_11x11 ***FlexRAN(Reference)***.\n");
            tsc_start_flexr   = gms::common::rdtscp_fenced();
            matrix_inv_cholesky_11x11(matBRe,matBIm,matInvBReRef,matInvBImRef);
            tsc_end_flexr     = gms::common::rdtscp_fenced();
            ptr_invBReRef = reinterpret_cast<const float * __restrict>(&matInvBReRef);
            ptr_invBImRef = reinterpret_cast<const float * __restrict>(&matInvBImRef);
            printf_ret = printf("[UNIT-TEST:] -- Dumping results:\n");
            volatile std::int32_t bre_err_count = 0;
            volatile std::int32_t bim_err_count = 0;
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                const float bre     = ptr_invBRe[__i];
                const float bre_ref = ptr_invBReRef[__i];
                const float bim     = ptr_invBIm[__i];
                const float bim_ref = ptr_invBImRef[__i];
                const bool  bre_cmp = almostEqual2sComplement(bre,bre_ref,7);
                const bool  bim_cmp = almostEqual2sComplement(bim,bim_ref,7);
                if(!bre_cmp)
                {   
                    ++bre_err_count;
                    [[maybe_unused]] std::int32_t 
                    printf_res_re =  printf("[UNIT-TEST]: ***FAILED*** -- (RF-EMT)bre=%.7f,(FlexRAN)bre_ref=%.7f,fail-count=%d\n",bre,bre_ref,bre_err_count);
                    printNumber("(RF-EMT)bre",bre,0);
                    printNumber("(FlexRAN)bre_ref",bre_ref,0);
                }   
                if(!bim_cmp)
                { 
                    ++bim_err_count;
                    [[maybe_unused]] std::int32_t 
                    printf_res_im =  printf("[UNIT-TEST]:  ***FAILED*** -- (RF-EMT)bim=%.7f,(FlexRAN)bim_ref=%.7f,fail-count=%d\n",bim,bim_ref,bim_err_count);
                    printNumber("(RF-EMT)bim",bim,0);
                    printNumber("(FlexRAN)bim_ref",bim_ref,0);
                } 
                if(__builtin_expect(detect_subnormal(bre),0))     { std::printf("[UNIT-TEST:] -- (RF-EMT)  ***SUBNORMAL*** in I-ch -- %.17f\n",bre);}
                if(__builtin_expect(detect_subnormal(bim),0))     { std::printf("[UNIT-TEST:] -- (RF-EMT)  ***SUBNORMAL*** in Q-ch -- %.17f\n",bim);}
                if(__builtin_expect(detect_subnormal(bre_ref),0)) { std::printf("[UNIT-TEST:] -- (FlexRAN) ***SUBNORMAL*** in I-ch -- %.17f\n",bre_ref);}
                if(__builtin_expect(detect_subnormal(bim_ref),0)) { std::printf("[UNIT-TEST:] -- (FlexRAN) ***SUBNORMAL*** in Q-ch -- %.17f\n",bim_ref);}
            } 
            if(!bre_err_count) {printf("[UNIT_TEST:] -- Validity Test: (I-channel): ***PASSED***\n");}
            if(!bim_err_count) {printf("[UNIT_TEST:] -- Validity Test: (Q-channel): ***PASSED***\n");}
            PRINT_CRUDE_TIMING_RESULTS;
            printf_ret = printf("[UNIT-TEST]: -- END: mat_inv_cholesky_11x11_16xf32<true,true> vs matrix_inv_cholesky_11x11 (Validity).\n");
        }
        break;
        default : 
        printf_ret = printf("[UNIT-TEST:] -- Invalid switch variable=%d\n",which);
                 return;
    }   

}

int main()
{   
    unit_test_matrix_inv_cholesky_11x11_validity();
    return 0;
}

