#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <sched.h>
#include <cerrno>
#include <sys/types.h>
#include <sys/resource.h> //getpid 
#include <signal.h>
#include <unistd.h>
#include <immintrin.h>
#include <omp.h> 
#include <algorithm>
#include <random>
#include <ctime>
#include <functional>
#include "GMS_malloc.h"
#include "GMS_ansi_colors.h"
#include "GMS_matrix_inv_cholesky.h"
#include "unit_test_mat_inv_cholesky_defines.h"
#include "GMS_machine_utils.h"

/*
   icpc -o perf_test_mat_inv_cholesky_16x16_corrected -std=c++17 -fp-model fast=2 -fopenmp -ftz -ggdb -ipo -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_machine_utils.h GMS_matrix_inv_cholesky.h GMS_matrix_inv_cholesky.cpp unit_test_mat_inv_cholesky_defines.h perf_test_mat_inv_cholesky_16x16_corrected.cpp 
   ASM: 
   icpc -S -fverbose-asm -masm=intel -fopenmp -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_machine_utils.h GMS_matrix_inv_cholesky.h GMS_matrix_inv_cholesky.cpp unit_test_mat_inv_cholesky_defines.h perf_test_mat_inv_cholesky_16x16_corrected.cpp

*/

static __m512 constZero = _mm512_setzero_ps();
static __m512 vzero     = constZero;
static __m512 vneg_one   = _mm512_set1_ps(-1.0f);


void matrix_inv_cholesky_16x16(__m512 matBRe[N_16][N_16], __m512 matBIm[N_16][N_16],
    __m512 matInvBRe[N_16][N_16], __m512 matInvBIm[N_16][N_16])
{
    // temp matrix and variables for matrix inversion
    __m512 matGRe[N_16][N_16], matGIm[N_16][N_16];
    __m512 matLRe[N_16][N_16], matLIm[N_16][N_16];
    __m512 matD[N_16], matND[N_16];
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
    GET_G_COL0_REF(matGRe, matGIm, matBRe, matBIm, matD, 11);
    GET_G_COL0_REF(matGRe, matGIm, matBRe, matBIm, matD, 12);
    GET_G_COL0_REF(matGRe, matGIm, matBRe, matBIm, matD, 13);
    GET_G_COL0_REF(matGRe, matGIm, matBRe, matBIm, matD, 14);
    GET_G_COL0_REF(matGRe, matGIm, matBRe, matBIm, matD, 15);

    // Column 1
    GET_G11_REF_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0);
    GET_G_COL1_REF(matGRe, matGIm, matBRe, matBIm, matD, 2, temp0, temp1);
    GET_G_COL1_REF(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);
    GET_G_COL1_REF(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL1_REF(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL1_REF(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL1_REF(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL1_REF(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);
    GET_G_COL1_REF(matGRe, matGIm, matBRe, matBIm, matD, 9, temp0, temp1);
    GET_G_COL1_REF(matGRe, matGIm, matBRe, matBIm, matD, 10, temp0, temp1);
    GET_G_COL1_REF(matGRe, matGIm, matBRe, matBIm, matD, 11, temp0, temp1);
    GET_G_COL1_REF(matGRe, matGIm, matBRe, matBIm, matD, 12, temp0, temp1);
    GET_G_COL1_REF(matGRe, matGIm, matBRe, matBIm, matD, 13, temp0, temp1);
    GET_G_COL1_REF(matGRe, matGIm, matBRe, matBIm, matD, 14, temp0, temp1);
    GET_G_COL1_REF(matGRe, matGIm, matBRe, matBIm, matD, 15, temp0, temp1);

    // Column 2
    GET_G22_REF_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    GET_G_COL2_REF(matGRe, matGIm, matBRe, matBIm, matD, 3, temp0, temp1);
    GET_G_COL2_REF(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL2_REF(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL2_REF(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL2_REF(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL2_REF(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);
    GET_G_COL2_REF(matGRe, matGIm, matBRe, matBIm, matD, 9, temp0, temp1);
    GET_G_COL2_REF(matGRe, matGIm, matBRe, matBIm, matD, 10, temp0, temp1);
    GET_G_COL2_REF(matGRe, matGIm, matBRe, matBIm, matD, 11, temp0, temp1);
    GET_G_COL2_REF(matGRe, matGIm, matBRe, matBIm, matD, 12, temp0, temp1);
    GET_G_COL2_REF(matGRe, matGIm, matBRe, matBIm, matD, 13, temp0, temp1);
    GET_G_COL2_REF(matGRe, matGIm, matBRe, matBIm, matD, 14, temp0, temp1);
    GET_G_COL2_REF(matGRe, matGIm, matBRe, matBIm, matD, 15, temp0, temp1);

    // Column 3
    GET_G33_REF_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    GET_G_COL3_REF(matGRe, matGIm, matBRe, matBIm, matD, 4, temp0, temp1);
    GET_G_COL3_REF(matGRe, matGIm, matBRe, matBIm, matD, 5, temp0, temp1);
    GET_G_COL3_REF(matGRe, matGIm, matBRe, matBIm, matD, 6, temp0, temp1);
    GET_G_COL3_REF(matGRe, matGIm, matBRe, matBIm, matD, 7, temp0, temp1);
    GET_G_COL3_REF(matGRe, matGIm, matBRe, matBIm, matD, 8, temp0, temp1);
    GET_G_COL3_REF(matGRe, matGIm, matBRe, matBIm, matD, 9, temp0, temp1);
    GET_G_COL3_REF(matGRe, matGIm, matBRe, matBIm, matD, 10, temp0, temp1);
    GET_G_COL3_REF(matGRe, matGIm, matBRe, matBIm, matD, 11, temp0, temp1);
    GET_G_COL3_REF(matGRe, matGIm, matBRe, matBIm, matD, 12, temp0, temp1);
    GET_G_COL3_REF(matGRe, matGIm, matBRe, matBIm, matD, 13, temp0, temp1);
    GET_G_COL3_REF(matGRe, matGIm, matBRe, matBIm, matD, 14, temp0, temp1);
    GET_G_COL3_REF(matGRe, matGIm, matBRe, matBIm, matD, 15, temp0, temp1);

    // Column 4
    GET_G44_REF_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 5, 4, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 6, 4, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 7, 4, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 8, 4, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 9, 4, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 10, 4, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 11, 4, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 12, 4, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 13, 4, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 14, 4, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 15, 4, temp0, temp1);

    // Column 5
    GET_G55_REF_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 6, 5, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 7, 5, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 8, 5, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 9, 5, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 10, 5, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 11, 5, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 12, 5, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 13, 5, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 14, 5, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 15, 5, temp0, temp1);

    // Column 6
    GET_G66_REF_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 7, 6, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 8, 6, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 9, 6, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 10, 6, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 11, 6, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 12, 6, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 13, 6, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 14, 6, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 15, 6, temp0, temp1);

    // Column 7
    GET_G77_REF_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 8, 7, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 9, 7, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 10, 7, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 11, 7, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 12, 7, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 13, 7, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 14, 7, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 15, 7, temp0, temp1);

    // Column 8
    GET_Gii_EVEN_REF_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 8);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 9, 8, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 10, 8, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 11, 8, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 12, 8, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 13, 8, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 14, 8, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 15, 8, temp0, temp1);

    // Column 9
    GET_Gii_ODD_REF_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 9);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 10, 9, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 11, 9, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 12, 9, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 13, 9, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 14, 9, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 15, 9, temp0, temp1);

    // Column 10
    GET_Gii_EVEN_REF_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 10);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 11, 10, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 12, 10, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 13, 10, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 14, 10, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 15, 10, temp0, temp1);

    // Column 11
    GET_Gii_ODD_REF_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 11);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 12, 11, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 13, 11, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 14, 11, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 15, 11, temp0, temp1);

    // Column 12
    GET_Gii_EVEN_REF_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 12);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 13, 12, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 14, 12, temp0, temp1);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 15, 12, temp0, temp1);

    // Column 13
    GET_Gii_ODD_REF_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 13);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 14, 13, temp0, temp1);
    GET_G_COL_ODD_REF(matGRe, matGIm, matBRe, matBIm, matD, 15, 13, temp0, temp1);

    // Column 14
    GET_Gii_EVEN_REF_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 14);
    GET_G_COL_EVEN_REF(matGRe, matGIm, matBRe, matBIm, matD, 15, 14, temp0, temp1);

    // Column 15
    GET_Gii_ODD_REF_SAFE(matGRe, matGIm, matBRe, matD, matND, temp0, temp1, temp2, 15);


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
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 11, 0, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 12, 0, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 13, 0, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 14, 0, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 15, 0, temp0, temp1);

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
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 11, 1, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 12, 1, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 13, 1, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 14, 1, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 15, 1, temp0, temp1);

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
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 11, 2, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 12, 2, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 13, 2, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 14, 2, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 15, 2, temp0, temp1);

    // Column 3
    SET_Lii_REF(matLRe, matLIm, matD, 3);
    GET_L_i1i_REF(matLRe, matLIm, matGRe, matGIm, matND, 4, 3);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 5, 3, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 6, 3, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 7, 3, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 8, 3, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 9, 3, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 10, 3, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 11, 3, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 12, 3, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 13, 3, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 14, 3, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 15, 3, temp0, temp1);

    // Column 4
    SET_Lii_REF(matLRe, matLIm, matD, 4);
    GET_L_i1i_REF(matLRe, matLIm, matGRe, matGIm, matND, 5, 4);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 6, 4, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 7, 4, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 8, 4, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 9, 4, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 10, 4, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 11, 4, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 12, 4, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 13, 4, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 14, 4, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 15, 4, temp0, temp1);

    // Column 5
    SET_Lii_REF(matLRe, matLIm, matD, 5);
    GET_L_i1i_REF(matLRe, matLIm, matGRe, matGIm, matND, 6, 5);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 7, 5, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 8, 5, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 9, 5, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 10, 5, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 11, 5, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 12, 5, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 13, 5, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 14, 5, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 15, 5, temp0, temp1);

    // Column 6
    SET_Lii_REF(matLRe, matLIm, matD, 6);
    GET_L_i1i_REF(matLRe, matLIm, matGRe, matGIm, matND, 7, 6);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 8, 6, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 9, 6, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 10, 6, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 11, 6, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 12, 6, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 13, 6, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 14, 6, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 15, 6, temp0, temp1);

    // Column 7
    SET_Lii_REF(matLRe, matLIm, matD, 7);
    GET_L_i1i_REF(matLRe, matLIm, matGRe, matGIm, matND, 8, 7);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 9, 7, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 10, 7, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 11, 7, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 12, 7, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 13, 7, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 14, 7, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 15, 7, temp0, temp1);

    // Column 8
    SET_Lii_REF(matLRe, matLIm, matD, 8);
    GET_L_i1i_REF(matLRe, matLIm, matGRe, matGIm, matND, 9, 8);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 10, 8, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 11, 8, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 12, 8, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 13, 8, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 14, 8, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 15, 8, temp0, temp1);

    // Column 9
    SET_Lii_REF(matLRe, matLIm, matD, 9);
    GET_L_i1i_REF(matLRe, matLIm, matGRe, matGIm, matND, 10, 9);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 11, 9, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 12, 9, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 13, 9, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 14, 9, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 15, 9, temp0, temp1);

    // Column 10
    SET_Lii_REF(matLRe, matLIm, matD, 10);
    GET_L_i1i_REF(matLRe, matLIm, matGRe, matGIm, matND, 11, 10);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 12, 10, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 13, 10, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 14, 10, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 15, 10, temp0, temp1);

    // Column 11
    SET_Lii_REF(matLRe, matLIm, matD, 11);
    GET_L_i1i_REF(matLRe, matLIm, matGRe, matGIm, matND, 12, 11);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 13, 11, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 14, 11, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 15, 11, temp0, temp1);

    // Column 12
    SET_Lii_REF(matLRe, matLIm, matD, 12);
    GET_L_i1i_REF(matLRe, matLIm, matGRe, matGIm, matND, 13, 12);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 14, 12, temp0, temp1);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 15, 12, temp0, temp1);

    // Column 13
    SET_Lii_REF(matLRe, matLIm, matD, 13);
    GET_L_i1i_REF(matLRe, matLIm, matGRe, matGIm, matND, 14, 13);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 15, 13, temp0, temp1);

    // Column 14
    SET_Lii_REF(matLRe, matLIm, matD, 14);
    GET_L_i1i_REF(matLRe, matLIm, matGRe, matGIm, matND, 15, 14);

    // Column 15
    SET_Lii_REF(matLRe, matLIm, matD, 15);


    /////////////////////////////////// get invB = L'*L
    for(i = 0; i < N_16; i ++)
    {
        matInvBRe[i][i] = _mm512_mul_ps(matLRe[i][i], matLRe[i][i]);
        for (k = (i+1); k < N_16; k++)
        {
            temp1 = _mm512_add_ps(_mm512_mul_ps(matLRe[k][i], matLRe[k][i]), _mm512_mul_ps(matLIm[k][i], matLIm[k][i]));
            matInvBRe[i][i] = _mm512_add_ps(matInvBRe[i][i], temp1);
        }
        matInvBIm[i][i] = _mm512_setzero_ps();
    }

    for(i = 0; i < N_16; i ++)
    {
        for(j = i+1; j < N_16; j ++)
        {
            matInvBRe[i][j] = _mm512_mul_ps(matLRe[j][i], matLRe[j][j]);
            matInvBIm[i][j] = _mm512_sub_ps(constZero, _mm512_mul_ps(matLIm[j][i], matLRe[j][j]));

            for (k = (j+1); k < N_16; k++)
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define BUFFER_STORE_SIZE 80
#define FORMAT_STORE_SIZE 80 

// More accurate for the floating-point comparison.
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

__attribute__((hot))
int32_t set_affinity_and_priority(const int32_t,const int32_t);
// 0 -- success, 1 and  2 failure.
int32_t set_affinity_and_priority(const int32_t cpu, const int32_t priority)
{
    cpu_set_t cpu_set;
    sched_param sp;
    int32_t status{-1};
    CPU_ZERO(&cpu_set);
    CPU_SET(cpu,&cpu_set);
    if(sched_setaffinity(0,sizeof(cpu_set), &cpu_set) < 0) 
    {
         status = 1;
         return status;
    }
    printf("Affinity set to cpu: %d\n",cpu);
    __builtin_memset(&sp,0,sizeof(sp));
    sp.sched_priority = priority; //99
    if((sched_setscheduler(0,SCHED_FIFO,&sp)) == -1) 
    {
        status = 2;
        return status;
    }
    status = 0;
    return status;
}

__attribute__((hot))
void print_thread_affinity();

void print_thread_affinity()
{
     char default_format[FORMAT_STORE_SIZE];
     char format_specifier[] = "host=%20H tid=%0.4n binds_to=%A";
     char buffer[BUFFER_STORE_SIZE];
     std::size_t nchars{};
     std::size_t diff;
     nchars = omp_get_affinity_format(default_format,(std::size_t)FORMAT_STORE_SIZE);
     diff   = nchars-(std::size_t)FORMAT_STORE_SIZE;
     if(diff>0ull)
        nchars += diff;
     omp_set_affinity_format(format_specifier);
     nchars = omp_capture_affinity(&buffer[0],(std::size_t)BUFFER_STORE_SIZE,NULL);
     printf("tid=%d affinity:%s\n",omp_get_thread_num(),buffer);
}

__attribute__((hot))
__attribute__((noinline))
void perf_test_mat_inv_cholesky_16x16_flexran(__m512 matBRe[N_16][N_16],    __m512 matBIm[N_16][N_16],
                                              __m512 matInvBReRef[N_16][N_16], __m512 matInvBImRef[N_16][N_16],
                                              unsigned __int64 * __restrict__ ,
                                              unsigned __int64 * __restrict__ ,
                                              unsigned __int64 * __restrict__ ,
                                              std::uint32_t    * __restrict__ ,
                                              std::uint32_t    * __restrict__ ,
                                              const int32_t ,
                                              const int32_t,
                                              uint32_t &);

__attribute__((hot))
__attribute__((noinline))
void perf_test_mat_inv_cholesky_16x16_flexran(__m512 matBRe[N_16][N_16],    __m512 matBIm[N_16][N_16],
                                              __m512 matInvBReRef[N_16][N_16], __m512 matInvBImRef[N_16][N_16],
                                              unsigned __int64 * __restrict__ matinv_flexran_s,
                                              unsigned __int64 * __restrict__ matinv_flexran_e,
                                              unsigned __int64 * __restrict__ matinv_flexran_d,
                                              std::uint32_t    * __restrict__ matinv_flexran_tsc_aux_s,
                                              std::uint32_t    * __restrict__ matinv_flexran_tsc_aux_e,
                                              const int32_t n_runs,
                                              const int32_t n_samples,
                                              uint32_t & tid)
{
    tid = 9999;
    /*
    Latency and Throughput
    Architecture	Latency	Throughput (CPI)
    Skylake	          42	     5.5
*/
    constexpr std::size_t RDTSCP_LAT{42ull};
    constexpr std::int32_t zmm_size{16};
    unsigned __int64 * __restrict__ p_matinv_flexran_s = matinv_flexran_s;
    unsigned __int64 * __restrict__ p_matinv_flexran_e = matinv_flexran_e;
    unsigned __int64 * __restrict__ p_matinv_flexran_d = matinv_flexran_d;
    std::uint32_t    * __restrict__ p_matinv_flexran_tsc_aux_s{nullptr};
    std::uint32_t    * __restrict__ p_matinv_flexran_tsc_aux_e{nullptr};
    const bool cond_p_tsc_aux_s_neq_nullptr = p_matinv_flexran_tsc_aux_s!=nullptr;
    const bool cond_p_tsc_aux_e_neq_nullptr = p_matinv_flexran_tsc_aux_e!=nullptr;
    if(cond_p_tsc_aux_s_neq_nullptr && cond_p_tsc_aux_e_neq_nullptr)
    {
        p_matinv_flexran_tsc_aux_s = matinv_flexran_tsc_aux_s;
        p_matinv_flexran_tsc_aux_e = matinv_flexran_tsc_aux_e;
        [[maybe_unused]]
        volatile unsigned __int64 warmup_rdtscp_s;
        [[maybe_unused]]
        volatile unsigned __int64 warmup_rdtscp_e;
        [[maybe_unused]]
        std::uint32_t mem_start{};
        std::uint32_t mem_end{};
        warmup_rdtscp_s = __rdtscp(&mem_start);
        warmup_rdtscp_e = __rdtscp(&mem_end);
        for(int32_t __i{0}; __i != n_runs; ++__i)
        {
            for(int32_t __j{0}; __j != n_samples; ++__j) 
            {   
                __asm__ __volatile__ ("lfence");
                unsigned __int64 start_curr{__rdtscp(&mem_start)};       
                matrix_inv_cholesky_16x16(matBRe,matBIm,matInvBReRef,matInvBImRef);
                unsigned __int64 end_curr{__rdtscp(&mem_end)};
                __asm__ __volatile__ ("lfence");
                //remove latency
                unsigned __int64 start_corrected{start_curr-RDTSCP_LAT};
                p_matinv_flexran_s[__i*n_samples+__j] = start_corrected;
                unsigned __int64 end_corrected{end_curr-RDTSCP_LAT};
                p_matinv_flexran_e[__i*n_samples+__j]   = end_corrected;
                p_matinv_flexran_d[__i*n_samples+__j] = end_corrected-start_corrected;
                p_matinv_flexran_tsc_aux_s[__i*n_samples+__j] = mem_start;
                p_matinv_flexran_tsc_aux_e[__i*n_samples+__j] = mem_end;
           }
        }
    }
    else 
    {
        [[maybe_unused]]
        volatile unsigned __int64 warmup_rdtscp_fenced;
        warmup_rdtscp_fenced = gms::common::rdtscp_fenced();
        for(int32_t __i{0}; __i != n_runs; ++__i)
        {
            for(int32_t __j{0}; __j != n_samples; ++__j) 
            {   
                unsigned __int64 start_curr{gms::common::rdtscp_fenced()};       
                matrix_inv_cholesky_16x16(matBRe,matBIm,matInvBReRef,matInvBImRef);
                unsigned __int64 end_curr{gms::common::rdtscp_fenced()};
                //remove latency (approximately equal = 42)
                unsigned __int64 start_corrected{start_curr-RDTSCP_LAT};
                p_matinv_flexran_s[__i*n_samples+__j] = start_corrected;
                unsigned __int64 end_corrected{end_curr-RDTSCP_LAT};
                p_matinv_flexran_e[__i*n_samples+__j] = end_corrected;
                p_matinv_flexran_d[__i*n_samples+__j] = end_corrected-start_corrected;
           }
        }
    }
    tid = omp_get_thread_num();
}

using namespace gms::math::mat_inv_chol_detail;

__attribute__((hot))
__attribute__((noinline))
void perf_test_mat_inv_cholesky_16x16_rf_emt(__m512 matBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],__m512 matBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                              __m512 matInvBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16], __m512 matInvBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                              unsigned __int64 * __restrict__ ,
                                              unsigned __int64 * __restrict__ ,
                                              unsigned __int64 * __restrict__ ,
                                              std::uint32_t    * __restrict__ ,
                                              std::uint32_t    * __restrict__ ,
                                              const int32_t ,
                                              const int32_t,
                                              uint32_t &);

__attribute__((hot))
__attribute__((noinline))
void perf_test_mat_inv_cholesky_16x16_rf_emt(__m512 matBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],__m512 matBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                              __m512 matInvBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16], __m512 matInvBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16],
                                              unsigned __int64 * __restrict__ matinv_rf_emt_s,
                                              unsigned __int64 * __restrict__ matinv_rf_emt_e,
                                              unsigned __int64 * __restrict__ matinv_rf_emt_d,
                                              std::uint32_t    * __restrict__ matinv_rf_emt_tsc_aux_s,
                                              std::uint32_t    * __restrict__ matinv_rf_emt_tsc_aux_e,
                                              const int32_t n_runs,
                                              const int32_t n_samples,
                                              uint32_t & tid)
{
    tid = 9999;
    /*
    Latency and Throughput
    Architecture	Latency	Throughput (CPI)
    Skylake	          42	     5.5
*/
    constexpr std::size_t RDTSCP_LAT{42ull};
    constexpr std::int32_t zmm_size{16};
    unsigned __int64 * __restrict__ p_matinv_rf_emt_s = matinv_rf_emt_s;
    unsigned __int64 * __restrict__ p_matinv_rf_emt_e = matinv_rf_emt_e;
    unsigned __int64 * __restrict__ p_matinv_rf_emt_d = matinv_rf_emt_d;
    std::uint32_t    * __restrict__ p_matinv_rf_emt_tsc_aux_s{nullptr};
    std::uint32_t    * __restrict__ p_matinv_rf_emt_tsc_aux_e{nullptr};
    const bool cond_p_tsc_aux_s_neq_nullptr = p_matinv_rf_emt_tsc_aux_s!=nullptr;
    const bool cond_p_tsc_aux_e_neq_nullptr = p_matinv_rf_emt_tsc_aux_e!=nullptr;
    if(cond_p_tsc_aux_s_neq_nullptr && cond_p_tsc_aux_e_neq_nullptr)
    {
        p_matinv_rf_emt_tsc_aux_s = matinv_rf_emt_tsc_aux_s;
        p_matinv_rf_emt_tsc_aux_e = matinv_rf_emt_tsc_aux_e;
        [[maybe_unused]]
        volatile unsigned __int64 warmup_rdtscp_s;
        [[maybe_unused]]
        volatile unsigned __int64 warmup_rdtscp_e;
        [[maybe_unused]]
        std::uint32_t mem_start{};
        std::uint32_t mem_end{};
        warmup_rdtscp_s = __rdtscp(&mem_start);
        warmup_rdtscp_e = __rdtscp(&mem_end);
        for(int32_t __i{0}; __i != n_runs; ++__i)
        {
            for(int32_t __j{0}; __j != n_samples; ++__j) 
            {   
                __asm__ __volatile__ ("lfence");
                unsigned __int64 start_curr{__rdtscp(&mem_start)};       
                gms::math::mat_inv_cholesky_16x16_16xf32<true,true>(matBRe,matBIm,
                                                               matInvBRe,matInvBIm);
                unsigned __int64 end_curr{__rdtscp(&mem_end)};
                __asm__ __volatile__ ("lfence");
                //remove latency
                unsigned __int64 start_corrected{start_curr-RDTSCP_LAT};
                p_matinv_rf_emt_s[__i*n_samples+__j] = start_corrected;
                unsigned __int64 end_corrected{end_curr-RDTSCP_LAT};
                p_matinv_rf_emt_e[__i*n_samples+__j]   = end_corrected;
                p_matinv_rf_emt_d[__i*n_samples+__j] = end_corrected-start_corrected;
                p_matinv_rf_emt_tsc_aux_s[__i*n_samples+__j] = mem_start;
                p_matinv_rf_emt_tsc_aux_e[__i*n_samples+__j] = mem_end;
           }
        }
    }
    else 
    {
        [[maybe_unused]]
        volatile unsigned __int64 warmup_rdtscp_fenced;
        warmup_rdtscp_fenced = gms::common::rdtscp_fenced();
        for(int32_t __i{0}; __i != n_runs; ++__i)
        {
            for(int32_t __j{0}; __j != n_samples; ++__j) 
            {   
                unsigned __int64 start_curr{gms::common::rdtscp_fenced()};       
                gms::math::mat_inv_cholesky_16x16_16xf32<true,true>(matBRe,matBIm,
                                                               matInvBRe,matInvBIm);
                unsigned __int64 end_curr{gms::common::rdtscp_fenced()};
                //remove latency (approximately equal = 42)
                unsigned __int64 start_corrected{start_curr-RDTSCP_LAT};
                p_matinv_rf_emt_s[__i*n_samples+__j] = start_corrected;
                unsigned __int64 end_corrected{end_curr-RDTSCP_LAT};
                p_matinv_rf_emt_e[__i*n_samples+__j] = end_corrected;
                p_matinv_rf_emt_d[__i*n_samples+__j] = end_corrected-start_corrected;
           }
        }
    }
    tid = omp_get_thread_num();
}

#define VEC_PS_LEN 16
#define COMPARE_RESULTS 1

__attribute__((hot))
void test_runner_omp_sections_1st_seq();

void test_runner_omp_sections_1st_seq()
{
    using namespace gms::common;
    constexpr int32_t n_runs{10};
    constexpr int32_t n_samples{50};
    constexpr int32_t n_total{n_runs*n_samples};
    constexpr int32_t nelems{4096};
    __ATTR_ALIGN__(64) float buf_init_re[nelems];
    __ATTR_ALIGN__(64) float buf_init_im[nelems];
    __ATTR_ALIGN__(64) __m512 matBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16];
    __ATTR_ALIGN__(64) __m512 matBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16];
    __ATTR_ALIGN__(64) __m512 matInvBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16];//rf-emt
    __ATTR_ALIGN__(64) __m512 matInvBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16];//rf-emt
    __ATTR_ALIGN__(64) __m512 matInvBReRef[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16];//flexran
    __ATTR_ALIGN__(64) __m512 matInvBImRef[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16];//flexran    
    unsigned __int64 flexran_s[n_total]           = {UINT64_MAX};
    unsigned __int64 flexran_e[n_total]           = {UINT64_MAX}; 
    unsigned __int64 flexran_d[n_total]           = {UINT64_MAX};
    unsigned __int64 rf_emt_s[n_total]            = {UINT64_MAX};
    unsigned __int64 rf_emt_e[n_total]            = {UINT64_MAX};
    unsigned __int64 rf_emt_d[n_total]            = {UINT64_MAX};
    std::uint32_t    flexran_tsc_aux_s[n_total]   = {UINT32_MAX};
    std::uint32_t    flexran_tsc_aux_e[n_total]   = {UINT32_MAX};
    std::uint32_t    rf_emt_tsc_aux_s[n_total]    = {UINT32_MAX};
    std::uint32_t    rf_emt_tsc_aux_e[n_total]    = {UINT32_MAX};
    const float * __restrict ptr_invBRe     = nullptr;
    const float * __restrict ptr_invBIm     = nullptr;
    const float * __restrict ptr_invBReRef  = nullptr;
    const float * __restrict ptr_invBImRef  = nullptr;
    unsigned long long seed_re{0ull};
    unsigned long long seed_im{0ull};
    [[maybe_unused]] std::int32_t printf_ret{};
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
    for(std::int32_t i = 0;i != MAT_SQR_SIZE_16; ++i)
    {   
        const std::int32_t outer_idx = i*MAT_SQR_SIZE_16;
        for(std::int32_t j = 0;j != MAT_SQR_SIZE_16; ++j)  
        {
            const std::int32_t inner_idx = (outer_idx+j)*VEC_PS_LEN;
            register __m512 init_vre = _mm512_load_ps((const float*)&buf_init_re[inner_idx]);                    
            matBRe[i][j]             = init_vre;
            //matInvBReRef[i][j]       = init_vre;
            register __m512 init_vim = _mm512_load_ps((const float*)&buf_init_im[inner_idx]);
            matBIm[i][j]             = init_vim;
            //matInvBImRef[i][j]       = init_vim;
        }
    }
    uint32_t flexran_tid;
    uint32_t rf_emt_tid;
    int32_t setenv_ret;
    setenv_ret = setenv("OMP_PROC_BIND","true",1);
    if(setenv_ret==-1)
    {
        printf_ret = printf("[**ERROR**]: -- setenv reported an error=%d\n",setenv_ret);
    }
    setenv_ret = setenv("OMP_PROC_BIND","spread",1);
    if(setenv_ret==-1)
    {
        printf_ret = printf("[**ERROR**]: -- setenv reported an error=%d\n",setenv_ret);
    }

#pragma omp parallel sections 
{
    #pragma omp section 
    {
        perf_test_mat_inv_cholesky_16x16_flexran(matBRe,matBIm,matInvBReRef,matInvBImRef,
                                                 &flexran_s[0],&flexran_e[0],&flexran_d[0],
                                                 nullptr,nullptr,n_runs,n_samples,flexran_tid);
        print_thread_affinity();
    }

    #pragma omp section 
    {
        perf_test_mat_inv_cholesky_16x16_rf_emt(matBRe,matBIm,matInvBRe,matInvBIm,
                                                &rf_emt_s[0],&rf_emt_e[0],&rf_emt_d[0],
                                                nullptr,nullptr,n_runs,n_samples,rf_emt_tid);
        print_thread_affinity();
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////
    printf_ret = printf(ANSI_COLOR_GREEN "[PERF-TEST]:  matrix_inv_cholesky_16x16(FlexRAN): COMPLETED!!\n");
    printf_ret = printf(ANSI_COLOR_GREEN "[PERF-TEST] -- Executed by Core=%d \n",flexran_tid);
    for(int32_t __i{0}; __i != n_runs; ++__i)
    {
         for(int32_t __j{0}; __j != n_samples; ++__j) 
         {  
            unsigned __int64 s{flexran_s[__i*n_samples+__j]};
            unsigned __int64 e{flexran_e[__i*n_samples+__j]};
            unsigned __int64 d{flexran_d[__i*n_samples+__j]};
            printf(ANSI_COLOR_GREEN "[PMC: RDTSCP] -- Run=%d, start=%llu,end=%llu,delta=%llu\n",__i,s,e,d);
         }
    }
    printf_ret = printf(ANSI_COLOR_GREEN "[PERF-TEST]: matrix_inv_cholesky_16x16(FlexRAN): Finished -- dumping-results" ANSI_RESET_ALL"\n\n");
    
    printf_ret = printf(ANSI_COLOR_WHITE "[PERF-TEST]:  matrix_inv_cholesky_16x16<true,true>(RF-EMT): COMPLETED!!\n");
    printf_ret = printf(ANSI_COLOR_WHITE "[PERF-TEST] -- Executed by Core=%d \n",rf_emt_tid);
    for(int32_t __i{0}; __i != n_runs; ++__i)
    {
         for(int32_t __j{0}; __j != n_samples; ++__j) 
         {  
            unsigned __int64 s{rf_emt_s[__i*n_samples+__j]};
            unsigned __int64 e{rf_emt_e[__i*n_samples+__j]};
            unsigned __int64 d{rf_emt_d[__i*n_samples+__j]};
            printf(ANSI_COLOR_WHITE "[PMC: RDTSCP] -- Run=%d, start=%llu,end=%llu,delta=%llu\n",__i,s,e,d);
         }
    }
    printf_ret = printf(ANSI_COLOR_WHITE "[PERF-TEST]: matrix_inv_cholesky_16x16<true,true>(RF-EMT): Finished -- dumping-results" ANSI_RESET_ALL"\n\n");
#if (COMPARE_RESULTS) == 1
    ptr_invBRe    = reinterpret_cast<const float * __restrict>(&matInvBRe);
    ptr_invBIm    = reinterpret_cast<const float * __restrict>(&matInvBIm);
    ptr_invBReRef = reinterpret_cast<const float * __restrict>(&matInvBReRef);
    ptr_invBImRef = reinterpret_cast<const float * __restrict>(&matInvBImRef);
    volatile std::int32_t bre_err_count = 0;
    volatile std::int32_t bim_err_count = 0;
    for(int32_t __i{0}; __i != nelems; ++__i)
    {  
        const float bre     = ptr_invBRe[__i];
        const float bre_ref = ptr_invBReRef[__i];
        const float bim     = ptr_invBIm[__i];
        const float bim_ref = ptr_invBImRef[__i];
        const bool  bre_cmp = almostEqual2sComplement(bre,bre_ref,9);
        const bool  bim_cmp = almostEqual2sComplement(bim,bim_ref,9);
        if(!bre_cmp)
        {
            ++bre_err_count;
            printf_ret =  printf(ANSI_COLOR_RED "[PERF-TEST]: ***FAILED*** -- (RF-EMT)bre=%.7f,(FlexRAN)bre_ref=%.7f,fail-count=%d\n",bre,bre_ref,bre_err_count);
            printNumber("(RF-EMT)bre",bre,0);
            printNumber("(FlexRAN)bre_ref",bre_ref,0);
        }   
        if(!bim_cmp)
        { 
            ++bim_err_count;
            printf_ret =  printf(ANSI_COLOR_RED "[PERF-TEST]:  ***FAILED*** -- (RF-EMT)bim=%.7f,(FlexRAN)bim_ref=%.7f,fail-count=%d\n",bim,bim_ref,bim_err_count);
            printNumber("(RF-EMT)bim",bim,0);
            printNumber("(FlexRAN)bim_ref",bim_ref,0);
        } 
        if(__builtin_expect(detect_subnormal(bre),0))     { printf_ret = printf(ANSI_COLOR_RED "[PERF-TEST:] -- (RF-EMT)  ***SUBNORMAL*** in I-ch -- %.17f\n",bre);}
        if(__builtin_expect(detect_subnormal(bim),0))     { printf_ret = printf(ANSI_COLOR_RED "[PERF-TEST:] -- (RF-EMT)  ***SUBNORMAL*** in Q-ch -- %.17f\n",bim);}
        if(__builtin_expect(detect_subnormal(bre_ref),0)) { printf_ret = printf(ANSI_COLOR_RED "[PERF-TEST:] -- (FlexRAN) ***SUBNORMAL*** in I-ch -- %.17f\n",bre_ref);}
        if(__builtin_expect(detect_subnormal(bim_ref),0)) { printf_ret = printf(ANSI_COLOR_RED "[PERF-TEST:] -- (FlexRAN) ***SUBNORMAL*** in Q-ch -- %.17f\n",bim_ref);}       
    } 
    if(!bre_err_count) {printf_ret = printf(ANSI_COLOR_GREEN "[PERF-TEST:] -- Validity Test: (I-channel): ***PASSED***\n");}
    if(!bim_err_count) {printf_ret = printf(ANSI_COLOR_GREEN "[PERF-TEST:] -- Validity Test: (Q-channel): ***PASSED***\n");}
#endif 
}

__attribute__((hot))
void test_runner_matrix_inv_cholesky_16x16_flexran_single_thread(const std::int32_t,
                                                                 const std::int32_t);

void test_runner_matrix_inv_cholesky_16x16_flexran_single_thread(const std::int32_t cpu_core,
                                                                 const std::int32_t priority)
{
    using namespace gms::common;
    constexpr int32_t n_runs{10};
    constexpr int32_t n_samples{50};
    constexpr int32_t n_total{n_runs*n_samples};
    constexpr int32_t nelems{4096};
    __ATTR_ALIGN__(64) float buf_init_re[nelems];
    __ATTR_ALIGN__(64) float buf_init_im[nelems];
    __ATTR_ALIGN__(64) __m512 matBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16];
    __ATTR_ALIGN__(64) __m512 matBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16];
    __ATTR_ALIGN__(64) __m512 matInvBReRef[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16];//flexran
    __ATTR_ALIGN__(64) __m512 matInvBImRef[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16];//flexran    
    unsigned __int64 flexran_s[n_total]           = {UINT64_MAX};
    unsigned __int64 flexran_e[n_total]           = {UINT64_MAX}; 
    unsigned __int64 flexran_d[n_total]           = {UINT64_MAX};
    std::uint32_t    flexran_tsc_aux_s[n_total]   = {UINT32_MAX};
    std::uint32_t    flexran_tsc_aux_e[n_total]   = {UINT32_MAX};
    const float * __restrict ptr_invBReRef  = nullptr;
    const float * __restrict ptr_invBImRef  = nullptr;
    unsigned long long seed_re{0ull};
    unsigned long long seed_im{0ull};
    [[maybe_unused]] std::int32_t printf_ret{};
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
    for(std::int32_t i = 0;i != MAT_SQR_SIZE_16; ++i)
    {   
        const std::int32_t outer_idx = i*MAT_SQR_SIZE_16;
        for(std::int32_t j = 0;j != MAT_SQR_SIZE_16; ++j)  
        {
            const std::int32_t inner_idx = (outer_idx+j)*VEC_PS_LEN;
            register __m512 init_vre = _mm512_load_ps((const float*)&buf_init_re[inner_idx]);                    
            matBRe[i][j]             = init_vre;
            //matInvBReRef[i][j]       = init_vre;
            register __m512 init_vim = _mm512_load_ps((const float*)&buf_init_im[inner_idx]);
            matBIm[i][j]             = init_vim;
            //matInvBImRef[i][j]       = init_vim;
        }
    }
    std::uint32_t flexran_tid{};
    std::int32_t setenv_ret{};
    std::int32_t status{};
    setenv_ret = setenv("OMP_PROC_BIND","true",1);
    if(setenv_ret==-1)
    {
        printf_ret = printf("[**ERROR**]: -- setenv reported an error=%d\n",setenv_ret);
    }
    setenv_ret = setenv("OMP_PROC_BIND","spread",1);
    if(setenv_ret==-1)
    {
        printf_ret = printf("[**ERROR**]: -- setenv reported an error=%d\n",setenv_ret);
    }
    status = set_affinity_and_priority(cpu_core,priority);
    if(status > 0)
    { 
       int32_t which = PRIO_PROCESS;
       id_t pid      = getpid();
       errno         = 0;
       int32_t nice  = getpriority(which,pid);
       printf_ret = printf("[***ERROR***]: set_affinity_and_priority status:%d,errno=%d,pid=%d,nice=%d\n",status,errno,pid,nice);
    }

    perf_test_mat_inv_cholesky_16x16_flexran(matBRe,matBIm,matInvBReRef,matInvBImRef,
                                            &flexran_s[0],&flexran_e[0],&flexran_d[0],
                                            nullptr,nullptr,n_runs,n_samples,flexran_tid);
    print_thread_affinity();

//////////////////////////////////////////////////////////////////////////////////////////////////
    printf_ret = printf(ANSI_COLOR_GREEN "[PERF-TEST]:  matrix_inv_cholesky_16x16(FlexRAN): COMPLETED!!\n");
    printf_ret = printf(ANSI_COLOR_GREEN "[PERF-TEST] -- Executed by Core=%d \n",flexran_tid);
    for(int32_t __i{0}; __i != n_runs; ++__i)
    {
         for(int32_t __j{0}; __j != n_samples; ++__j) 
         {  
            unsigned __int64 s{flexran_s[__i*n_samples+__j]};
            unsigned __int64 e{flexran_e[__i*n_samples+__j]};
            unsigned __int64 d{flexran_d[__i*n_samples+__j]};
            printf(ANSI_COLOR_GREEN "[PMC: RDTSCP] -- Run=%d, start=%llu,end=%llu,delta=%llu\n",__i,s,e,d);
         }
    }
    printf_ret = printf(ANSI_COLOR_GREEN "[PERF-TEST]: matrix_inv_cholesky_16x16(FlexRAN): Finished -- dumping-results" ANSI_RESET_ALL"\n\n");
    
#if (COMPARE_RESULTS) == 0
    ptr_invBReRef = reinterpret_cast<const float * __restrict>(&matInvBReRef);
    ptr_invBImRef = reinterpret_cast<const float * __restrict>(&matInvBImRef);
    for(int32_t __i{0}; __i != nelems; ++__i)
    {  
        const float bre_ref = ptr_invBReRef[__i];
        const float bim_ref = ptr_invBImRef[__i];
        printf_ret =  printf(ANSI_COLOR_RED "[PERF-TEST]: -- (FlexRAN)bre_ref=%.7f\n",bre_ref);
        printNumber("(FlexRAN)bre_ref",bre_ref,0);
        printf_ret =  printf(ANSI_COLOR_RED "[PERF-TEST]: -- (FlexRAN)bim_ref=%.7f\n",bim_ref);
        printNumber("(FlexRAN)bim_ref",bim_ref,0);
         
        if(__builtin_expect(detect_subnormal(bre_ref),0)) { printf_ret = printf(ANSI_COLOR_RED "[PERF-TEST:] -- (FlexRAN) ***SUBNORMAL*** in I-ch -- %.17f\n",bre_ref);}
        if(__builtin_expect(detect_subnormal(bim_ref),0)) { printf_ret = printf(ANSI_COLOR_RED "[PERF-TEST:] -- (FlexRAN) ***SUBNORMAL*** in Q-ch -- %.17f\n",bim_ref);}       
    } 
#endif 
}

__attribute__((hot))
void test_runner_matrix_inv_cholesky_16x16_rf_emt_single_thread(const std::int32_t,
                                                                 const std::int32_t);

void test_runner_matrix_inv_cholesky_16x16_rf_emt_single_thread(const std::int32_t cpu_core,
                                                                 const std::int32_t priority)
{
    using namespace gms::common;
    constexpr int32_t n_runs{10};
    constexpr int32_t n_samples{50};
    constexpr int32_t n_total{n_runs*n_samples};
    constexpr int32_t nelems{4096};
    __ATTR_ALIGN__(64) float buf_init_re[nelems];
    __ATTR_ALIGN__(64) float buf_init_im[nelems];
    __ATTR_ALIGN__(64) __m512 matBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16];
    __ATTR_ALIGN__(64) __m512 matBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16];
    __ATTR_ALIGN__(64) __m512 matInvBRe[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16];//flexran
    __ATTR_ALIGN__(64) __m512 matInvBIm[MAT_SQR_SIZE_16][MAT_SQR_SIZE_16];//flexran    
    unsigned __int64 rf_emt_s[n_total]           = {UINT64_MAX};
    unsigned __int64 rf_emt_e[n_total]           = {UINT64_MAX}; 
    unsigned __int64 rf_emt_d[n_total]           = {UINT64_MAX};
    std::uint32_t    rf_emt_tsc_aux_s[n_total]   = {UINT32_MAX};
    std::uint32_t    rf_emt_tsc_aux_e[n_total]   = {UINT32_MAX};
    const float * __restrict ptr_invBRe  = nullptr;
    const float * __restrict ptr_invBIm  = nullptr;
    unsigned long long seed_re{0ull};
    unsigned long long seed_im{0ull};
    [[maybe_unused]] std::int32_t printf_ret{};
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
    for(std::int32_t i = 0;i != MAT_SQR_SIZE_16; ++i)
    {   
        const std::int32_t outer_idx = i*MAT_SQR_SIZE_16;
        for(std::int32_t j = 0;j != MAT_SQR_SIZE_16; ++j)  
        {
            const std::int32_t inner_idx = (outer_idx+j)*VEC_PS_LEN;
            register __m512 init_vre = _mm512_load_ps((const float*)&buf_init_re[inner_idx]);                    
            matBRe[i][j]             = init_vre;
            //matInvBReRef[i][j]       = init_vre;
            register __m512 init_vim = _mm512_load_ps((const float*)&buf_init_im[inner_idx]);
            matBIm[i][j]             = init_vim;
            //matInvBImRef[i][j]       = init_vim;
        }
    }
    std::uint32_t rf_emt_tid{};
    std::int32_t setenv_ret{};
    std::int32_t status{};
    setenv_ret = setenv("OMP_PROC_BIND","true",1);
    if(setenv_ret==-1)
    {
        printf_ret = printf("[**ERROR**]: -- setenv reported an error=%d\n",setenv_ret);
    }
    setenv_ret = setenv("OMP_PROC_BIND","spread",1);
    if(setenv_ret==-1)
    {
        printf_ret = printf("[**ERROR**]: -- setenv reported an error=%d\n",setenv_ret);
    }
    status = set_affinity_and_priority(cpu_core,priority);
    if(status > 0)
    { 
       int32_t which = PRIO_PROCESS;
       id_t pid      = getpid();
       errno         = 0;
       int32_t nice  = getpriority(which,pid);
       printf_ret = printf("[***ERROR***]: set_affinity_and_priority status:%d,errno=%d,pid=%d,nice=%d\n",status,errno,pid,nice);
    }

    perf_test_mat_inv_cholesky_16x16_rf_emt(matBRe,matBIm,matInvBRe,matInvBIm,
                                            &rf_emt_s[0],&rf_emt_e[0],&rf_emt_d[0],
                                            nullptr,nullptr,n_runs,n_samples,rf_emt_tid);
    print_thread_affinity();

//////////////////////////////////////////////////////////////////////////////////////////////////
    printf_ret = printf(ANSI_COLOR_GREEN "[PERF-TEST]:  matrix_inv_cholesky_16x16<true,true>(RF-EMT): COMPLETED!!\n");
    printf_ret = printf(ANSI_COLOR_GREEN "[PERF-TEST] -- Executed by Core=%d \n",rf_emt_tid);
    for(int32_t __i{0}; __i != n_runs; ++__i)
    {
         for(int32_t __j{0}; __j != n_samples; ++__j) 
         {  
            unsigned __int64 s{rf_emt_s[__i*n_samples+__j]};
            unsigned __int64 e{rf_emt_e[__i*n_samples+__j]};
            unsigned __int64 d{rf_emt_d[__i*n_samples+__j]};
            printf(ANSI_COLOR_GREEN "[PMC: RDTSCP] -- Run=%d, start=%llu,end=%llu,delta=%llu\n",__i,s,e,d);
         }
    }
    printf_ret = printf(ANSI_COLOR_GREEN "[PERF-TEST]: matrix_inv_cholesky_16x16<true,true>(RF-EMT): Finished -- dumping-results" ANSI_RESET_ALL"\n\n");
    
#if (COMPARE_RESULTS) == 0
    ptr_invBRe = reinterpret_cast<const float * __restrict>(&matInvBRe);
    ptr_invBIm = reinterpret_cast<const float * __restrict>(&matInvBIm);
    for(int32_t __i{0}; __i != nelems; ++__i)
    {  
        const float bre     = ptr_invBRe[__i];
        const float bim     = ptr_invBIm[__i];
        printf_ret =  printf(ANSI_COLOR_RED "[PERF-TEST]: (RF-EMT)bre=%.7f\n",bre);
        printNumber("(RF-EMT)bre",bre,0);
        printf_ret =  printf(ANSI_COLOR_RED "[PERF-TEST]: (RF-EMT)bim=%.7f\n",bim);
        printNumber("(RF-EMT)bim",bim,0);
        if(__builtin_expect(detect_subnormal(bre),0)) { printf_ret = printf(ANSI_COLOR_RED "[PERF-TEST:] -- (RF-EMT) ***SUBNORMAL*** in I-ch -- %.17f\n",bre);}
        if(__builtin_expect(detect_subnormal(bim),0)) { printf_ret = printf(ANSI_COLOR_RED "[PERF-TEST:] -- (RF-EMT) ***SUBNORMAL*** in Q-ch -- %.17f\n",bim);}       
    } 
#endif 
}

__attribute__((hot))
void test_runner_matrix_inv_cholesky_16x16_flexran_parallel();

void test_runner_matrix_inv_cholesky_16x16_flexran_parallel()
{
    thread_local std::uniform_int_distribution<std::int32_t> rand_thread_prio;
    thread_local std::uint64_t prio_seed{};
    thread_local std::int32_t priority{};
    thread_local std::mt19937 rand_prio_gen;
    int32_t tid{};
#pragma omp parallel private(tid) num_threads(omp_get_max_threads())
{
    rand_thread_prio = std::uniform_int_distribution<std::int32_t>(0,99);
    prio_seed        = __rdtsc();
    rand_prio_gen    = std::mt19937(prio_seed);
    priority         = rand_thread_prio.operator()(rand_prio_gen);
    tid              = omp_get_thread_num();
    test_runner_matrix_inv_cholesky_16x16_flexran_single_thread(tid,priority);
}

}

__attribute__((hot))
void test_runner_matrix_inv_cholesky_16x16_rf_emt_parallel();

void test_runner_matrix_inv_cholesky_16x16_rf_emt_parallel()
{
    thread_local std::uniform_int_distribution<std::int32_t> rand_thread_prio;
    thread_local std::uint64_t prio_seed{};
    thread_local std::int32_t priority{};
    thread_local std::mt19937 rand_prio_gen;
    int32_t tid{};
#pragma omp parallel private(tid) num_threads(omp_get_max_threads())
{
    rand_thread_prio = std::uniform_int_distribution<std::int32_t>(0,99);
    prio_seed        = __rdtsc();
    rand_prio_gen    = std::mt19937(prio_seed);
    priority         = rand_thread_prio.operator()(rand_prio_gen);
    tid              = omp_get_thread_num();
    test_runner_matrix_inv_cholesky_16x16_rf_emt_single_thread(tid,priority);
}

}

int main()
{
    test_runner_omp_sections_1st_seq();
    test_runner_matrix_inv_cholesky_16x16_flexran_single_thread(5,99);
    test_runner_matrix_inv_cholesky_16x16_rf_emt_single_thread(6,99);
    test_runner_matrix_inv_cholesky_16x16_flexran_parallel();
    test_runner_matrix_inv_cholesky_16x16_rf_emt_parallel();
    return 0;
}