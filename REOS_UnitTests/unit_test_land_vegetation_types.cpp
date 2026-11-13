
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_malloc.h"
#include "GMS_land_vegetation_types.h"
/*
   icpc -o unit_test_land_vegetation_types -O3 -I /home/user/boost_1_85_0 -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_avx512_memset.h GMS_avx512_memset.cpp GMS_half.h GMS_boost_multiarray_types.hpp GMS_module_land_params.h GMS_land_vegeatation_types.h unit_test_land_vegetation_types.cpp
   ASM: 
   icpc  -S  -O3 -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 -I /home/user/boost_1_85_0 GMS_config.h GMS_avx512_memset.h GMS_avx512_memset.cpp GMS_boost_multi_array_wrapper.hpp GMS_module_land_params.h GMS_land_vegeatation_types.h unit_test_land_vegetation_types.cpp
   The new icpc command

    icpc -o unit_test_land_vegetation_types -O3 -I/home/bernard/Software-Development/Cpp-projects/RF_EMT/RF-EMT_unit_tests/Unit_Test_land_vegetation_types -I /home/bernard/Software-Development/Cpp-projects/RF_EMT/RF-EMT_src/boost_1_89_0 -fp-model fast=2 -std=c++17 -ftz -ggdb -Wall -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5 GMS_avx512_memset.cpp  unit_test_land_vegetation_types.cpp
*/


__ATTR_HOT__ 
__ATTR_ALIGN__(32)
std::int32_t unit_test_land_vegetation_types();

std::int32_t unit_test_land_vegetation_types()
{
    char * ctor_name = NULL;
    std::int32_t ret_status = -9999;
    using namespace gms::radiolocation;
    using namespace gms;
    [[maybe_unused]]std::int32_t printf_ret = -9999;
    std::int32_t demangle_status = -9999;
    printf_ret =  printf("[UNIT-TEST]: function=%s -- **STARTED**, Date=%s,Time=%s\n", __PRETTY_FUNCTION__,__DATE__,__TIME__);
    printf_ret =  printf("[UNIT-TEST]: [1] -- Calling Default Constructor...\n");
    vegetation_state_avx512_t veg_state{};
    ctor_name = gms::common::demangle(typeid(veg_state).name(),demangle_status);
    if(demangle_status==0 && ctor_name != NULL)
    {
        printf_ret = printf("[UNIT-TEST]: [1]-- Instantiation of Default Constructor of type: %s\n", ctor_name);
        gms::common::gms_mm_free(ctor_name);
    }
    else
    {
        printf_ret = printf("[UNIT-TEST]: [1] -- Instantiation of Default Constructor of type: 'abi::__cxa_demangle failed'!! %s\n", typeid(veg_state).name());
        if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
    } 
    constexpr std::size_t beg = 0ull;
    constexpr std::size_t end = 1024;
    constexpr std::size_t nband = 6;
    constexpr std::size_t nlevgrnd = 2;
    vegetation_state_avx512_t veg_state2(beg,end,nband,nlevgrnd);
    ctor_name = gms::common::demangle(typeid(veg_state2).name(),demangle_status);
    if(demangle_status==0 && ctor_name != NULL)
    {
        printf_ret = printf("[UNIT-TEST]: [2]-- Instantiation of Four Argument Constructor of type: %s\n", ctor_name);
        gms::common::gms_mm_free(ctor_name);
    }
    else
    {
        printf_ret = printf("[UNIT-TEST]: [2] -- Instantiation of Four Argument Constructor of type: 'abi::__cxa_demangle ***FAILED***!! %s\n", typeid(veg_state2).name());
        if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
    } 
    printf_ret = printf("[UNIT-TEST]: [2]-- Dumping the Four Argument Constructor of type: %s memory layout.\n", ctor_name);
    printf_ret = printf("[UNIT-TEST]: [2]-- =====================================================================================\n");
    std::int32_t ret_status_show_class_memory_layout = veg_state2.show_class_memory_layout();
    if(ret_status_show_class_memory_layout!=0)
    {
        printf_ret = printf("[UNIT-TEST]: [2]-- Dumping the Four Argument Constructor of type: %s memory layout -- ***FAILED***.\n", ctor_name);
        ret_status = ret_status_show_class_memory_layout;
        return ret_status;
    }
    printf_ret =  printf("[UNIT-TEST]: function=%s -- **STARTED**, Date=%s,Time=%s\n", __PRETTY_FUNCTION__,__DATE__,__TIME__);
    ret_status = demangle_status;
    return (ret_status);
}

int main()
{
    std::int32_t generate_status = unit_test_land_vegetation_types();
    if(generate_status==0) 
    {
        std::int32_t printf_status = std::printf("[UNIT-TEST]: unit_test_land_vegetation_types -- ***PASSED***,Status=%d\n",generate_status);
    }
    else
    {
        std::int32_t printf_status = std::printf("[UNIT-TEST]: unit_test_land_vegetation_types -- ***FAILED***,Status=%d\n",generate_status);
    }
    return 0;
}
