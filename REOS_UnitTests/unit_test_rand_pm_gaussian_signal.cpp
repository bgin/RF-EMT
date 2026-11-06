#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_rand_pm_gaussian_signal_f32.h"

/*
   icpc -o unit_test_rand_pm_gaussian_signal_f32 -O3 -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_half.h GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_rand_pm_gaussian_signal_f32.h GMS_rand_pm_gaussian_signal_f32.cpp unit_test_rand_pm_gaussian_signal.cpp
   ASM: 
   icpc  -S  -O3 -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GGMS_config.h GMS_malloc.h GMS_half.h GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_rand_pm_gaussian_signal_f32.h GMS_rand_pm_gaussian_signal_f32.cpp unit_test_rand_pm_gaussian_signal.cpp

*/



__ATTR_HOT__ 
__ATTR_ALIGN__(32)
std::int32_t unit_test_rand_pm_gaussian_signal_f32();

std::int32_t unit_test_rand_pm_gaussian_signal_f32()
{
    std::int32_t ret_status = -9999;
    using namespace gms::radiolocation;
    using namespace gms;
    std::int32_t printf_ret = -9999;
    std::int32_t demangle_status = -9999;
    printf_ret =  printf("[UNIT-TEST]: function=%s -- **STARTED**, Date=%s,Time=%s\n", __PRETTY_FUNCTION__,__DATE__,__TIME__);
    printf_ret =  printf("[UNIT-TEST]: [1] -- Initializing Randomly Phase-Modulated Gaussian Signal Generator...\n");
    const std::size_t num_samples   = 1024ull;
    const std::size_t num_gaussians = 64;
    float             tstart        = 0.0f;
    float             tend          = 10.0f;
    float             min_sigma     = 0.05f;
    float             max_sigma     = 0.5f;
    float             max_amplitude = 2.54f;
    rand_pm_gaussian_signal_f32_t gauss_sum_gen(num_samples,num_gaussians,tstart,tend,min_sigma,max_sigma,max_amplitude);
    char * ctor_name{gms::common::demangle(typeid(gauss_sum_gen).name(),demangle_status)};
    if(demangle_status==0 && ctor_name != NULL)
    {
        printf_ret = printf("[UNIT-TEST]: [2]-- Instantiation of object Constructor of type: %s\n", ctor_name);
        gms::common::gms_mm_free(ctor_name);
    }
    else
    {
        printf_ret = printf("[UNIT-TEST]: [2] -- Instantiation of object Constructor of type: %s\n", typeid(gauss_sum_gen).name());
        if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
    } 
    printf_ret =  printf("[UNIT-TEST]: [3] -- Accessing the type=%s internal state...\n",ctor_name);
    gauss_sum_gen.m_signal_samples.info_size_alignment();
    printf_ret =  printf("[UNIT-TEST]: [4] -- Generating Randomly Phase-Modulated Gaussian Signal...\n");
    ret_status = gauss_sum_gen.generate<num_gaussians>();
    if(ret_status<0)
    {
        printf_ret =  printf("[UNIT-TEST]: [4] -- 'generate()' -- ***FAILED***\n");
        if(ret_status == -9999) 
        {
            printf_ret =  printf("[UNIT-TEST]: [4] -- 'generate()' -- ***The Execution Stopped Roughly At The Prologue!!***\n");
            return (ret_status);
        }
        else if(ret_status = -1)
        {
            printf_ret =  printf("[UNIT-TEST]: [4] -- 'generate()' -- ***Size Mismatch Between The Allowed Number of Gaussians And The Automatic Store!!***\n");
            return (ret_status);
        }
    }
    else 
    {
        printf_ret =  printf("[UNIT-TEST]: [4] -- Randomly Phase Modulated Gaussian Signal Is Generated Correctly, Status=%d\n",ret_status);
    }
    printf_ret =  printf("[UNIT-TEST]: [5] -- Displaying The State of the object of type=%s\n\n",ctor_name);
    std::cout << gauss_sum_gen;
    printf("[UNIT_TEST]: function=%s -- **ENDED CORRECTLY**,Date=%s,Time=%s\n", __PRETTY_FUNCTION__,__DATE__,__TIME__);
    return (ret_status);
}


int main()
{   
    std::int32_t generate_status = unit_test_rand_pm_gaussian_signal_f32();
    if(generate_status==0) 
    {
        std::int32_t printf_status = std::printf("[UNIT-TEST]: unit_test_rand_pm_gaussian_signal_f32() -- ***PASSED***,Status=%d\n",generate_status);
    }
    else
    {
        std::int32_t printf_status = std::printf("[UNIT-TEST]: unit_test_rand_pm_gaussian_signal_f32() -- ***FAILED***,Status=%d\n",generate_status);
    }
    return 0;
}