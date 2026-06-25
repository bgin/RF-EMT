#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <immintrin.h> // __rdtsc()
#include "GMS_malloc.h" // for demangle 
#include "GMS_matrix_inv_cholesky.h"

/*
   icpc -o unit_test_matrix_inv_cholesky_2x2_validity -fp-model fast=2 -ftz -ggdb -ipo -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_matrix_inv_cholesky.h GMS_matrix_inv_cholesky.cpp unit_test_matrix_inv_cholesky_2x2_validity.cpp 
   ASM: 
   icpc -S -fverbose-asm -masm=intel -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_matrix_inv_cholesky.h GMS_matrix_inv_cholesky.cpp unit_test_matrix_inv_cholesky_2x2_validity.cpp 
*/
 
static __m512 constZero = _mm512_setzero_ps();

// complex multiplication: A * A'
#define GET_AxAH_REF(re, im, out)\
{\
    out = _mm512_add_ps(_mm512_mul_ps(re, re), _mm512_mul_ps(im, im));\
}

// complex multiplication: A * B'
#define GET_AxBH_REF(are, aim, bre, bim, outre, outim)\
{\
    outre = _mm512_add_ps(_mm512_mul_ps(are, bre), _mm512_mul_ps(aim, bim));\
    outim = _mm512_sub_ps(_mm512_mul_ps(bre, aim), _mm512_mul_ps(are, bim));\
}

// complex multiplication: A * B
#define GET_AxB_REF(are, aim, bre, bim, outre, outim)\
{\
    outre = _mm512_sub_ps(_mm512_mul_ps(are, bre), _mm512_mul_ps(aim, bim));\
    outim = _mm512_add_ps(_mm512_mul_ps(bre, aim), _mm512_mul_ps(are, bim));\
}

// complex multiplication: A * real(B)
#define GET_AxRealB_REF(are, aim, bre, outre, outim)\
{\
    outre = _mm512_mul_ps(are, bre);\
    outim = _mm512_mul_ps(aim, bre);\
}

#define GET_G00_REF(matGRe, matBRe, matD, matND)\
{\
    matD[0] = _mm512_rsqrt14_ps(matBRe[0][0]);\
    matND[0] = _mm512_sub_ps(constZero, matD[0]);\
}

// get column 0 of matrix G
#define GET_G_COL0_REF(matGRe, matGIm, matBRe, matBIm, matD, j)\
{\
    matGRe[j][0] = _mm512_mul_ps(matBRe[j][0], matD[0]);\
    matGIm[j][0] = _mm512_mul_ps(matBIm[j][0], matD[0]);\
}

// get G11
#define GET_G11_REF(matGRe, matGIm, matBRe, matD, matND, temp0)\
{\
    GET_AxAH_REF(matGRe[1][0], matGIm[1][0], temp0);\
    matGRe[1][1] = _mm512_sub_ps(matBRe[1][1], temp0);\
    matD[1] = _mm512_rsqrt14_ps(matGRe[1][1]);\
    matND[1] = _mm512_sub_ps(constZero, matD[1]);\
}

// set value for Lii, diagonal element
#define SET_Lii_REF(matLRe, matLIm, matD, i)\
{\
    matLRe[i][i] = matD[i];\
    matLIm[i][i] = constZero;\
}

// get element L(i+1, i)
#define GET_L_i1i_REF(matLRe, matLIm, matGRe, matGIm, matND, j, i)\
{\
    GET_AxRealB_REF(matGRe[j][i], matGIm[j][i], matLRe[i][i],  matLRe[j][i], matLIm[j][i]);\
    matLRe[j][i] = _mm512_mul_ps(matLRe[j][i], matND[j]);\
    matLIm[j][i] = _mm512_mul_ps(matLIm[j][i], matND[j]);\
}

#define N_2 2 

__attribute__((hot))
__attribute__((noinline))
__attribute__((aligned(32)))
void matrix_inv_cholesky_2x2(__m512 matBRe[N_2][N_2], __m512 matBIm[N_2][N_2],
    __m512 matInvBRe[N_2][N_2], __m512 matInvBIm[N_2][N_2]);

void matrix_inv_cholesky_2x2(__m512 matBRe[N_2][N_2], __m512 matBIm[N_2][N_2],
    __m512 matInvBRe[N_2][N_2], __m512 matInvBIm[N_2][N_2])
{
    // temp matrix and variables for matrix inversion
    __m512 matGRe[N_2][N_2], matGIm[N_2][N_2];
    __m512 matLRe[N_2][N_2], matLIm[N_2][N_2];
    __m512 matD[N_2], matND[N_2];
    __m512 temp0, temp1;

    /////////////////////////////////// get G, B = G*G', G is a lower triangular matrix
    // Column 0
    GET_G00_REF(matGRe, matBRe, matD, matND);
    GET_G_COL0_REF(matGRe, matGIm, matBRe, matBIm, matD, 1);

    // Column 1
    GET_G11_REF(matGRe, matGIm, matBRe, matD, matND, temp0);


    /////////////////////////////////// get L, L=invG
    // Column 0
    SET_Lii_REF(matLRe, matLIm, matD, 0);
    GET_L_i1i_REF(matLRe, matLIm, matGRe, matGIm, matND, 1, 0);

    // Column 1
    SET_Lii_REF(matLRe, matLIm, matD, 1);


    /////////////////////////////////// get invB = L'*L
    matInvBRe[0][0] = _mm512_mul_ps(matLRe[0][0], matLRe[0][0]);
    temp1 = _mm512_add_ps(_mm512_mul_ps(matLRe[1][0], matLRe[1][0]), _mm512_mul_ps(matLIm[1][0], matLIm[1][0]));
    matInvBRe[0][0] = _mm512_add_ps(matInvBRe[0][0], temp1);
    matInvBIm[0][0] = _mm512_setzero_ps();

    matInvBRe[1][1] = _mm512_mul_ps(matLRe[1][1], matLRe[1][1]);
    matInvBIm[1][1] = _mm512_setzero_ps();

    matInvBRe[0][1] = _mm512_mul_ps(matLRe[1][0], matLRe[1][1]);
    matInvBIm[0][1] = _mm512_sub_ps(constZero, _mm512_mul_ps(matLIm[1][0], matLRe[1][1]));
    matInvBRe[1][0] = matInvBRe[0][1];
    matInvBIm[1][0] = _mm512_sub_ps(constZero, matInvBIm[0][1]);
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
__attribute__((aligned(32)))
void unit_test_matrix_inv_cholesky_2x2_validity();

void unit_test_matrix_inv_cholesky_2x2_validity()
{
    constexpr int32_t nelems{64};
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
    printf("[UNIT-TEST]: -- START: mat_inv_cholesky_2x2_16xf32<true,true> vs. matrix_inv_cholesky_2x2 (Validity).\n");
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
            __ATTR_ALIGN__(64) __m512 matBRe[2][2];
            __ATTR_ALIGN__(64) __m512 matBIm[2][2];
            __ATTR_ALIGN__(64) __m512 matInvBRe[2][2];
            __ATTR_ALIGN__(64) __m512 matInvBIm[2][2];
            __ATTR_ALIGN__(64) __m512 matInvBReRef[2][2];
            __ATTR_ALIGN__(64) __m512 matInvBImRef[2][2];
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
            matBRe[0][0] = _mm512_load_ps((const float*)&buf_init_re[0]);
            matBIm[0][0] = _mm512_load_ps((const float*)&buf_init_im[0]);
            matBRe[1][0] = _mm512_load_ps((const float*)&buf_init_re[16]);
            matBIm[1][0] = _mm512_load_ps((const float*)&buf_init_im[16]);
            matBRe[0][1] = _mm512_load_ps((const float*)&buf_init_re[32]);
            matBIm[0][1] = _mm512_load_ps((const float*)&buf_init_im[32]);
            matBRe[1][1] = _mm512_load_ps((const float*)&buf_init_re[48]);
            matBIm[1][1] = _mm512_load_ps((const float*)&buf_init_im[48]);
            printf("[UNIT-TEST]: -- START: mat_inv_cholesky_2x2_16xf32<true,true>.\n");
            gms::math::mat_inv_cholesky_2x2_16xf32<true,false>(matBRe,matBIm,
                                                               matInvBRe,matInvBIm);                                                  
            ptr_invBRe = reinterpret_cast<const float * __restrict>(&matInvBRe);
            ptr_invBIm = reinterpret_cast<const float * __restrict>(&matInvBIm);
            printf("[UNIT-TEST]: -- START: matrix_inv_cholesky_2x2(***FlexRAN(Reference)***).\n");
            matrix_inv_cholesky_2x2(matBRe,matBIm,matInvBReRef,matInvBImRef);
            ptr_invBReRef = reinterpret_cast<const float * __restrict>(&matInvBReRef);
            ptr_invBImRef = reinterpret_cast<const float * __restrict>(&matInvBImRef);
            printf("[UNIT-TEST:] -- Dumping results:\n");
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
                    printf_res_re =  printf("[UNIT-TEST]: ***FAILED*** -- (RF-EMT)bre=%.7f,(FlexRAN)bre_ref=%.7f\n",bre,bre_ref);
                    printNumber("(RF-EMT)bre",bre,0);
                    printNumber("(FlexRAN)bre_ref",bre_ref,0);
                }   
                if(!bim_cmp)
                { 
                    [[maybe_unused]] std::int32_t 
                    printf_res_im =  printf("[UNIT-TEST]:  ***FAILED*** -- (RF-EMT)bim=%.7f,(FlexRAN)bim_ref=%.7f\n",bim,bim_ref);
                    printNumber("(RF-EMT)bim",bim,0);
                    printNumber("(FlexRAN)bim_ref",bim_ref,0);
                } 
               
            } 
        }            
             printf("[UNIT-TEST]: -- END: mat_inv_cholesky_2x2_16xf32<true,true> vs matrix_inv_cholesky_2x2 (Validity).\n");
        
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
            __ATTR_ALIGN__(64) __m512 matBRe[2][2];
            __ATTR_ALIGN__(64) __m512 matBIm[2][2];
            __ATTR_ALIGN__(64) __m512 matInvBRe[2][2];
            __ATTR_ALIGN__(64) __m512 matInvBIm[2][2];
            __ATTR_ALIGN__(64) __m512 matInvBReRef[2][2];
            __ATTR_ALIGN__(64) __m512 matInvBImRef[2][2];
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
            matBRe[0][0] = _mm512_load_ps((const float*)&buf_init_re[0]);
            matBIm[0][0] = _mm512_load_ps((const float*)&buf_init_im[0]);
            matBRe[1][0] = _mm512_load_ps((const float*)&buf_init_re[16]);
            matBIm[1][0] = _mm512_load_ps((const float*)&buf_init_im[16]);
            matBRe[0][1] = _mm512_load_ps((const float*)&buf_init_re[32]);
            matBIm[0][1] = _mm512_load_ps((const float*)&buf_init_im[32]);
            matBRe[1][1] = _mm512_load_ps((const float*)&buf_init_re[48]);
            matBIm[1][1] = _mm512_load_ps((const float*)&buf_init_im[48]);
            printf("[UNIT-TEST]: -- START: mat_inv_cholesky_2x2_16xf32<true,true> ***RF-EMT*** .\n");
            gms::math::mat_inv_cholesky_2x2_16xf32<true,false>(matBRe,matBIm,
                                                               matInvBRe,matInvBIm);
            ptr_invBRe = reinterpret_cast<const float * __restrict>(&matInvBRe);
            ptr_invBIm = reinterpret_cast<const float * __restrict>(&matInvBIm);
            printf("[UNIT-TEST]: -- START: matrix_inv_cholesky_2x2 ***FlexRAN(Reference)***.\n");
            matrix_inv_cholesky_2x2(matBRe,matBIm,matInvBReRef,matInvBImRef);
            ptr_invBReRef = reinterpret_cast<const float * __restrict>(&matInvBReRef);
            ptr_invBImRef = reinterpret_cast<const float * __restrict>(&matInvBImRef);
            printf("[UNIT-TEST:] -- Dumping results:\n");
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
                    printf_res_re =  printf("[UNIT-TEST]: ***FAILED*** -- (RF-EMT)bre=%.7f,(FlexRAN)bre_ref=%.7f\n",bre,bre_ref);
                    printNumber("(RF-EMT)bre",bre,0);
                    printNumber("(FlexRAN)bre_ref",bre_ref,0);
                }   
                if(!bim_cmp)
                { 
                    [[maybe_unused]] std::int32_t 
                    printf_res_im =  printf("[UNIT-TEST]:  ***FAILED*** -- (RF-EMT)bim=%.7f,(FlexRAN)bim_ref=%.7f\n",bim,bim_ref);
                    printNumber("(RF-EMT)bim",bim,0);
                    printNumber("(FlexRAN)bim_ref",bim_ref,0);
                } 
            } 
            printf("[UNIT-TEST]: -- END: mat_inv_cholesky_2x2_16xf32<true,true> vs matrix_inv_cholesky_2x2 (Validity).\n");
        }
        break;
        default : 
        printf("[UNIT-TEST:] -- Invalid switch variable=%d\n",which);
                 return;
    }   

}

int main()
{   
    unit_test_matrix_inv_cholesky_2x2_validity();
    return 0;
}

