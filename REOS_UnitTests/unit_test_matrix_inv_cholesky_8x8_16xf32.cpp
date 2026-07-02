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

/*
   icpc -o unit_test_matrix_inv_cholesky_8x8_16xf32  -std=c++17 -fp-model fast=2 -ftz -ggdb -ipo -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_matrix_inv_cholesky.h GMS_matrix_inv_cholesky.cpp unit_test_matrix_inv_cholesky_8x8_16xf32.cpp 
   ASM: 
   icpc -S -fverbose-asm -masm=intel -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_matrix_inv_cholesky.h GMS_matrix_inv_cholesky.cpp unit_test_matrix_inv_cholesky_8x8_16xf32.cpp 
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
void unit_test_matrix_inv_cholesky_8x8();

void unit_test_matrix_inv_cholesky_8x8()
{
    using namespace gms::math;
    constexpr int32_t nelems{1024};
    __ATTR_ALIGN__(64) float buf_init_re[nelems];
    __ATTR_ALIGN__(64) float buf_init_im[nelems];
    const float * __restrict ptr_invBRe = nullptr;
    const float * __restrict ptr_invBIm = nullptr;
    unsigned long long seed_re{0ull};
    unsigned long long seed_im{0ull};
    std::clock_t seed_d{0ull};
    int32_t which(-1);
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
            __ATTR_ALIGN__(64) __m512 matBRe[8][8];
            __ATTR_ALIGN__(64) __m512 matBIm[8][8];
            __ATTR_ALIGN__(64) __m512 matInvBRe[8][8];
            __ATTR_ALIGN__(64) __m512 matInvBIm[8][8];
             seed_re = __rdtsc();
             auto rand_re{std::bind(std::uniform_real_distribution<float>(0.0f,1.0f),std::mt19937(seed_re))};
             seed_im = __rdtsc();
             auto rand_im{std::bind(std::uniform_real_distribution<float>(0.0f,1.0f),std::mt19937(seed_im))};
             for(int32_t __i{0}; __i != nelems; ++__i)
             {
                 float rd_re = rand_re();
                 buf_init_re[__i] = rd_re;
                 float rd_im = rand_im();
                 buf_init_im[__i] = rd_im;
             }
            for(std::int32_t i = 0;i != MAT_SQR_SIZE_8; ++i)
            {   
                const std::int32_t outer_idx = i*MAT_SQR_SIZE_8;
                for(std::int32_t j = 0;j != MAT_SQR_SIZE_8; ++j)  
                {
                    const std::int32_t inner_idx = (outer_idx+j)*VEC_PS_LEN;
                    matBRe[i][j] = _mm512_load_ps((const float*)&buf_init_re[inner_idx]);
                    matBIm[i][j] = _mm512_load_ps((const float*)&buf_init_im[inner_idx]);
                }
            }
            printf("[UNIT-TEST]: -- START: mat_inv_cholesky_8x8_16xf32<true,true>.\n");
            gms::math::mat_inv_cholesky_8x8_16xf32<true,true>(matBRe,matBIm,
                                                               matInvBRe,matInvBIm);
            ptr_invBRe = reinterpret_cast<const float * __restrict>(&matInvBRe);
            ptr_invBIm = reinterpret_cast<const float * __restrict>(&matInvBIm);
            printf("[UNIT-TEST:] -- Dumping results:\n");
            for(int32_t __i{0}; __i != nelems; ++__i)
            {

                [[maybe_unused]] std::int32_t 
                printf_res_re =  printf("[UNIT-TEST]: matBRe[%d]=%.7f,matInvBRe[%d]=%.7f\n",__i,buf_init_re[__i],__i,ptr_invBRe[__i]);
                printNumber("matBRe",buf_init_re[__i],0);
                printNumber("matInvBRe",ptr_invBRe[__i],0);
                    
                [[maybe_unused]] std::int32_t 
                printf_res_im =  printf("[UNIT-TEST]: matBIm[%d]=%.7f,matInvBIm[%d]=%.7f\n",__i,buf_init_im[__i],__i,ptr_invBIm[__i]);
                printNumber("matBIm",buf_init_im[__i],0);
                printNumber("matInvBIm",ptr_invBIm[__i],0);
                if(detect_subnormal(ptr_invBRe[__i])) { std::printf("[UNIT-TEST:] ***SUBNORMAL*** in I-ch -- %.17f\n",ptr_invBRe[__i]);}
                if(detect_subnormal(ptr_invBIm[__i])) { std::printf("[UNIT-TEST:] ***SUBNORMAL*** in Q-ch -- %.17f\n",ptr_invBIm[__i]);}
            } 
        }            
             printf("[UNIT-TEST]: -- END: mat_inv_cholesky_8x8_16xf32<true,true>.\n");
        
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
            __ATTR_ALIGN__(64) __m512 matBRe[8][8];
            __ATTR_ALIGN__(64) __m512 matBIm[8][8];
            __ATTR_ALIGN__(64) __m512 matInvBRe[8][8];
            __ATTR_ALIGN__(64) __m512 matInvBIm[8][8];
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
            for(std::int32_t i = 0;i != MAT_SQR_SIZE_8; ++i)
            {   
                const std::int32_t outer_idx = i*MAT_SQR_SIZE_8;
                for(std::int32_t j = 0;j != MAT_SQR_SIZE_8; ++j)  
                {
                    const std::int32_t inner_idx = (outer_idx+j)*VEC_PS_LEN;
                    matBRe[i][j] = _mm512_load_ps((const float*)&buf_init_re[inner_idx]);
                    matBIm[i][j] = _mm512_load_ps((const float*)&buf_init_im[inner_idx]);
                }
            }
            printf("[UNIT-TEST]: -- START: mat_inv_cholesky_8x8_16xf32<true,true>.\n");
            gms::math::mat_inv_cholesky_8x8_16xf32<true,true>(matBRe,matBIm,
                                                               matInvBRe,matInvBIm);
            ptr_invBRe = reinterpret_cast<const float * __restrict>(&matInvBRe);
            ptr_invBIm = reinterpret_cast<const float * __restrict>(&matInvBIm);
            printf("[UNIT-TEST:] -- Dumping results:\n");
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                [[maybe_unused]] std::int32_t 
                printf_res_re =  printf("[UNIT-TEST]: matBRe[%d]=%.7f,matInvBRe[%d]=%.7f\n",__i,buf_init_re[__i],__i,ptr_invBRe[__i]);
                printNumber("matBRe",buf_init_re[__i],0);
                printNumber("matInvBRe",ptr_invBRe[__i],0);
                    
                [[maybe_unused]] std::int32_t 
                printf_res_im =  printf("[UNIT-TEST]: matBIm[%d]=%.7f,matInvBIm[%d]=%.7f\n",__i,buf_init_im[__i],__i,ptr_invBIm[__i]);
                printNumber("matBIm",buf_init_im[__i],0);
                printNumber("matInvBIm",ptr_invBIm[__i],0);
                if(detect_subnormal(ptr_invBRe[__i])) { std::printf("[UNIT-TEST:] ***SUBNORMAL*** in I-ch -- %.17f\n",ptr_invBRe[__i]);}
                if(detect_subnormal(ptr_invBIm[__i])) { std::printf("[UNIT-TEST:] ***SUBNORMAL*** in Q-ch -- %.17f\n",ptr_invBIm[__i]);}
            } 
            printf("[UNIT-TEST]: -- END: mat_inv_cholesky_8x8_16xf32<true,true>.\n");
        }
        break;
        case 2 : 
        {
            __ATTR_ALIGN__(64) __m512 matBRe[8][8];
            __ATTR_ALIGN__(64) __m512 matBIm[8][8];
            __ATTR_ALIGN__(64) __m512 matInvBRe[8][8];
            __ATTR_ALIGN__(64) __m512 matInvBIm[8][8];
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                 float f_i = static_cast<float>(__i);
                 buf_init_re[__i] = buf_init_im[__i] = f_i;
            }
            for(std::int32_t i = 0;i != MAT_SQR_SIZE_8; ++i)
            {   
                const std::int32_t outer_idx = i*MAT_SQR_SIZE_8;
                for(std::int32_t j = 0;j != MAT_SQR_SIZE_8; ++j)  
                {
                    const std::int32_t inner_idx = (outer_idx+j)*VEC_PS_LEN;
                    matBRe[i][j] = _mm512_load_ps((const float*)&buf_init_re[inner_idx]);
                    matBIm[i][j] = _mm512_load_ps((const float*)&buf_init_im[inner_idx]);
                }
            }
            printf("[UNIT-TEST]: -- START: mat_inv_cholesky_8x8_16xf32<true,true>.\n");
            gms::math::mat_inv_cholesky_8x8_16xf32<true,true>(matBRe,matBIm,
                                                               matInvBRe,matInvBIm);
            ptr_invBRe = reinterpret_cast<const float * __restrict>(&matInvBRe);
            ptr_invBIm = reinterpret_cast<const float * __restrict>(&matInvBIm);
            printf("[UNIT-TEST:] -- Dumping results:\n");
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                [[maybe_unused]] std::int32_t 
                printf_res_re =  printf("[UNIT-TEST]: matBRe[%d]=%.7f,matInvBRe[%d]=%.7f\n",__i,buf_init_re[__i],__i,ptr_invBRe[__i]);
                printNumber("matBRe",buf_init_re[__i],0);
                printNumber("matInvBRe",ptr_invBRe[__i],0);
                    
                [[maybe_unused]] std::int32_t 
                printf_res_im =  printf("[UNIT-TEST]: matBIm[%d]=%.7f,matInvBIm[%d]=%.7f\n",__i,buf_init_im[__i],__i,ptr_invBIm[__i]);
                printNumber("matBIm",buf_init_im[__i],0);
                printNumber("matInvBIm",ptr_invBIm[__i],0);
                if(detect_subnormal(ptr_invBRe[__i])) { std::printf("[UNIT-TEST:] ***SUBNORMAL*** in I-ch -- %.17f\n",ptr_invBRe[__i]);}
                if(detect_subnormal(ptr_invBIm[__i])) { std::printf("[UNIT-TEST:] ***SUBNORMAL*** in Q-ch -- %.17f\n",ptr_invBIm[__i]);}
            } 
            printf("[UNIT-TEST]: -- END: mat_inv_cholesky_8x8_16xf32<true,true>.\n");
        }
        break;
        default : 
        printf("[UNIT-TEST:] -- Invalid switch variable=%d\n",which);
                 return;
    }   

}

int main()
{   
    unit_test_matrix_inv_cholesky_8x8();
    return 0;
}

