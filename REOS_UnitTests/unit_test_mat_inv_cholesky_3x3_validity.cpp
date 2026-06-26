#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <immintrin.h> // __rdtsc()
#include "GMS_malloc.h" // for demangle 
#include "GMS_matrix_inv_cholesky.h"
#include "unit_test_mat_inv_cholesky_defines.h"
/*
   icpc -o unit_test_matrix_inv_cholesky_2x2_validity -fp-model fast=2 -ftz -ggdb -ipo -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_matrix_inv_cholesky.h GMS_matrix_inv_cholesky.cpp unit_test_mat_inv_cholesky_defines.h unit_test_matrix_inv_cholesky_3x3_validity.cpp 
   ASM: 
   icpc -S -fverbose-asm -masm=intel -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_config.h GMS_malloc.h GMS_matrix_inv_cholesky.h GMS_matrix_inv_cholesky.cpp unit_test_mat_inv_cholesky_defines.h unit_test_matrix_inv_cholesky_3x3_validity.cpp 
*/
 
static __m512 constZero = _mm512_setzero_ps();



__attribute__((hot))
__attribute__((noinline))
__attribute__((aligned(32)))
void matrix_inv_cholesky_3x3(__m512 matBRe[N_3][N_3], __m512 matBIm[N_3][N_3],
    __m512 matInvBRe[N_3][N_3], __m512 matInvBIm[N_3][N_3]);

void matrix_inv_cholesky_3x3(__m512 matBRe[N_3][N_3], __m512 matBIm[N_3][N_3],
    __m512 matInvBRe[N_3][N_3], __m512 matInvBIm[N_3][N_3])
{
    // temp matrix and variables for matrix inversion
    __m512 matGRe[N_3][N_3], matGIm[N_3][N_3];
    __m512 matLRe[N_3][N_3], matLIm[N_3][N_3];
    __m512 matD[N_3], matND[N_3];
    __m512 temp0, temp1, temp2;
    int32_t i, j, k;

    /////////////////////////////////// get G, B = G*G', G is a lower triangular matrix
    // Column 0
    GET_G00_REF(matGRe, matBRe, matD, matND);
    GET_G_COL0_REF(matGRe, matGIm, matBRe, matBIm, matD, 1);
    GET_G_COL0_REF(matGRe, matGIm, matBRe, matBIm, matD, 2);

    // Column 1
    GET_G11_REF(matGRe, matGIm, matBRe, matD, matND, temp0);
    GET_G_COL1_REF(matGRe, matGIm, matBRe, matBIm, matD, 2, temp0, temp1);

    // Column 2
    GET_G22_REF(matGRe, matGIm, matBRe, matD, matND, temp0, temp1);


    /////////////////////////////////// get L, L=invG
    // Column 0
    SET_Lii_REF(matLRe, matLIm, matD, 0);
    GET_L_i1i_REF(matLRe, matLIm, matGRe, matGIm, matND, 1, 0);
    GET_L_ji_REF(matLRe, matLIm, matGRe, matGIm, matND, 2, 0, temp0, temp1);

    // Column 1
    SET_Lii_REF(matLRe, matLIm, matD, 1);
    GET_L_i1i_REF(matLRe, matLIm, matGRe, matGIm, matND, 2, 1);

    // Column 2
    SET_Lii_REF(matLRe, matLIm, matD, 2);


    /////////////////////////////////// get invB = L'*L
    for(i = 0; i < N_3; i ++)
    {
        matInvBRe[i][i] = _mm512_mul_ps(matLRe[i][i], matLRe[i][i]);
        for (k = (i+1); k < N_3; k++)
        {
            temp1 = _mm512_add_ps(_mm512_mul_ps(matLRe[k][i], matLRe[k][i]), _mm512_mul_ps(matLIm[k][i], matLIm[k][i]));
            matInvBRe[i][i] = _mm512_add_ps(matInvBRe[i][i], temp1);
        }
        matInvBIm[i][i] = _mm512_setzero_ps();
    }

    for(i = 0; i < N_3; i ++)
    {
        for(j = i+1; j < N_3; j ++)
        {
            matInvBRe[i][j] = _mm512_mul_ps(matLRe[j][i], matLRe[j][j]);
            matInvBIm[i][j] = _mm512_sub_ps(constZero, _mm512_mul_ps(matLIm[j][i], matLRe[j][j]));

            for (k = (j+1); k < N_3; k++)
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

#define VEC_PS_LEN 16

__attribute__((hot))
__attribute__((noinline))
__attribute__((aligned(32)))
void unit_test_matrix_inv_cholesky_3x3_validity();

void unit_test_matrix_inv_cholesky_3x3_validity()
{
    using namespace gms::math;
    constexpr int32_t nelems{144};
    __ATTR_ALIGN__(64) float buf_init_re[nelems];
    __ATTR_ALIGN__(64) float buf_init_im[nelems];
    const float * __restrict ptr_invBRe     = nullptr;
    const float * __restrict ptr_invBIm     = nullptr;
    const float * __restrict ptr_invBReRef  = nullptr;
    const float * __restrict ptr_invBImRef  = nullptr;
    unsigned long long seed_re{0ull};
    unsigned long long seed_im{0ull};
    std::clock_t seed_d{0ull};
    int32_t which(-1);
    printf("[UNIT-TEST]: -- START: mat_inv_cholesky_3x3_16xf32<true,true> vs. matrix_inv_cholesky_3x3 (Validity).\n");
    seed_d = std::clock();
    auto rand_d{std::bind(std::uniform_int_distribution<int32_t>(0,1),std::mt19937(seed_d))};
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
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n", distro_name_re);
                   gms::common::gms_mm_free(distro_name_re);
             }
             else
             {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n", typeid(uniform_distro_re).name());
                   if(distro_name_re != NULL) gms::common::gms_mm_free(distro_name_re);
             }
             std::uniform_real_distribution<float> uniform_distro_im;
             char * distro_name_im{NULL};
             status = 9999;
             distro_name_im = gms::common::demangle(typeid(uniform_distro_im).name(),status);
             if(distro_name_im!= NULL && status == 0)
             {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n", distro_name_im);
                   gms::common::gms_mm_free(distro_name_im);
             }
             else
             {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n", typeid(uniform_distro_im).name());
                   if(distro_name_re != NULL) gms::common::gms_mm_free(distro_name_im);
             }
            __ATTR_ALIGN__(64) __m512 matBRe[3][3];
            __ATTR_ALIGN__(64) __m512 matBIm[3][3];
            __ATTR_ALIGN__(64) __m512 matInvBRe[3][3];
            __ATTR_ALIGN__(64) __m512 matInvBIm[3][3];
            __ATTR_ALIGN__(64) __m512 matInvBReRef[3][3];
            __ATTR_ALIGN__(64) __m512 matInvBImRef[3][3];
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
            for(std::int32_t i = 0;i != MAT_SQR_SIZE_3; ++i)
            {   
                const std::int32_t outer_idx = i*MAT_SQR_SIZE_3;
                for(std::int32_t j = 0;j != MAT_SQR_SIZE_3; ++j)  
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
            printf("[UNIT-TEST]: -- START: mat_inv_cholesky_3x3_16xf32<true,true>.\n");
            gms::math::mat_inv_cholesky_3x3_16xf32<true,false>(matBRe,matBIm,
                                                               matInvBRe,matInvBIm);                                                  
            ptr_invBRe = reinterpret_cast<const float * __restrict>(&matInvBRe);
            ptr_invBIm = reinterpret_cast<const float * __restrict>(&matInvBIm);
            printf("[UNIT-TEST]: -- START: matrix_inv_cholesky_3x3(***FlexRAN(Reference)***).\n");
            matrix_inv_cholesky_3x3(matBRe,matBIm,matInvBReRef,matInvBImRef);
            ptr_invBReRef = reinterpret_cast<const float * __restrict>(&matInvBReRef);
            ptr_invBImRef = reinterpret_cast<const float * __restrict>(&matInvBImRef);
            printf("[UNIT-TEST:] -- Dumping results:\n");
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
               
            } 
            if(!bre_err_count) {printf("[UNIT_TEST:] -- Validity Test: (I-channel): ***PASSED***\n");}
            if(!bim_err_count) {printf("[UNIT_TEST:] -- Validity Test: (Q-channel): ***PASSED***\n");}
        }            
             printf("[UNIT-TEST]: -- END: mat_inv_cholesky_3x3_16xf32<true,true> vs matrix_inv_cholesky_3x3 (Validity).\n");
        
        break;
        case 1 : 
        {
            std::normal_distribution<float> normal_distro_re;
            char * distro_name_re{NULL};
            int32_t status{9999};
            distro_name_re = gms::common::demangle(typeid(normal_distro_re).name(),status);
            if(distro_name_re != NULL && status == 0)
            {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n", distro_name_re);
                   gms::common::gms_mm_free(distro_name_re);
            }
            else
            {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n", typeid(normal_distro_re).name());
                   if(distro_name_re != NULL) gms::common::gms_mm_free(distro_name_re);
            }
            std::normal_distribution<float> normal_distro_im;
            char * distro_name_im{NULL};
            status = 9999;
            distro_name_im = gms::common::demangle(typeid(normal_distro_im).name(),status);
            if(distro_name_im != NULL && status == 0)
            {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n", distro_name_im);
                   gms::common::gms_mm_free(distro_name_im);
            }
            else
            {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n", typeid(normal_distro_im).name());
                   if(distro_name_re != NULL) gms::common::gms_mm_free(distro_name_im);
            }
            __ATTR_ALIGN__(64) __m512 matBRe[3][3];
            __ATTR_ALIGN__(64) __m512 matBIm[3][3];
            __ATTR_ALIGN__(64) __m512 matInvBRe[3][3];
            __ATTR_ALIGN__(64) __m512 matInvBIm[3][3];
            __ATTR_ALIGN__(64) __m512 matInvBReRef[3][3];
            __ATTR_ALIGN__(64) __m512 matInvBImRef[3][3];
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
            for(std::int32_t i = 0;i != MAT_SQR_SIZE_3; ++i)
            {   
                const std::int32_t outer_idx = i*MAT_SQR_SIZE_3;
                for(std::int32_t j = 0;j != MAT_SQR_SIZE_3; ++j)  
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
            printf("[UNIT-TEST]: -- START: mat_inv_cholesky_3x3_16xf32<true,true> ***RF-EMT*** .\n");
            gms::math::mat_inv_cholesky_3x3_16xf32<true,false>(matBRe,matBIm,
                                                               matInvBRe,matInvBIm);
            ptr_invBRe = reinterpret_cast<const float * __restrict>(&matInvBRe);
            ptr_invBIm = reinterpret_cast<const float * __restrict>(&matInvBIm);
            printf("[UNIT-TEST]: -- START: matrix_inv_cholesky_3x3 ***FlexRAN(Reference)***.\n");
            matrix_inv_cholesky_3x3(matBRe,matBIm,matInvBReRef,matInvBImRef);
            ptr_invBReRef = reinterpret_cast<const float * __restrict>(&matInvBReRef);
            ptr_invBImRef = reinterpret_cast<const float * __restrict>(&matInvBImRef);
            printf("[UNIT-TEST:] -- Dumping results:\n");
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
                    [[maybe_unused]] std::int32_t 
                    printf_res_re =  printf("[UNIT-TEST]: ***FAILED*** -- (RF-EMT)bre=%.7f,(FlexRAN)bre_ref=%.7f,fail-count=%d\n",bre,bre_ref,bre_err_count);
                    printNumber("(RF-EMT)bre",bre,0);
                    printNumber("(FlexRAN)bre_ref",bre_ref,0);
                }   
                if(!bim_cmp)
                { 
                    [[maybe_unused]] std::int32_t 
                    printf_res_im =  printf("[UNIT-TEST]:  ***FAILED*** -- (RF-EMT)bim=%.7f,(FlexRAN)bim_ref=%.7f,fail-count=%d\n",bim,bim_ref,bim_err_count);
                    printNumber("(RF-EMT)bim",bim,0);
                    printNumber("(FlexRAN)bim_ref",bim_ref,0);
                } 
            } 
            if(!bre_err_count) {printf("[UNIT_TEST:] -- Validity Test: (I-channel): ***PASSED***\n");}
            if(!bim_err_count) {printf("[UNIT_TEST:] -- Validity Test: (Q-channel): ***PASSED***\n");}
            printf("[UNIT-TEST]: -- END: mat_inv_cholesky_3x3_16xf32<true,true> vs matrix_inv_cholesky_3x3 (Validity).\n");
        }
        break;
        default : 
        printf("[UNIT-TEST:] -- Invalid switch variable=%d\n",which);
                 return;
    }   

}

int main()
{   
    unit_test_matrix_inv_cholesky_3x3_validity();
    return 0;
}

