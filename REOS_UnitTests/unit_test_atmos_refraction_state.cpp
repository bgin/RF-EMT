#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_atmos_refraction_state.h"

/*
icpc -o unit_test_atmos_refract_state -fopenmp -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_half.h GMS_dyn_array  GMS_sse_memset.h GMS_sse_memset.cpp GMS_atmos_refraction_state.h GMS_atmos_refraction_state.cpp unit_test_atmos_refraction_state.cpp

icpc -S -fverbose-asm -masm=intel -fopenmp -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h GMS_half.h GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_atmos_refraction_state.h GMS_atmos_refraction_state.cpp unit_test_atmos_refraction_state.cpp
*/


__attribute__((hot))
__attribute__((noinline))
void unit_test_ctor_1_omp_alloc();

void unit_test_ctor_1_omp_alloc()
{
    using namespace gms::atmos;
    std::printf("[UNIT-TEST]: %s ---> HAS BEGUN\n", __PRETTY_FUNCTION__);
    constexpr std::size_t lo_bound{1000ull};
    constexpr std::size_t hi_bound{10000ull};
    const std::int32_t    n_params{114};
    std::clock_t seed_alloc;
    atmos_refraction_size_params_t sz_params;
    std::int32_t status;
    seed_alloc = std::clock();
    auto rand_alloc{std::bind(std::uniform_int_distribution<std::size_t>(lo_bound,hi_bound),std::mt19937(seed_alloc))};
    std::size_t * __restrict__ p_sz_params{nullptr};
    p_sz_params = (std::size_t*)&sz_params;
    for(std::int32_t i{0}; i != n_params; ++i) 
    {
        p_sz_params[i] = rand_alloc();
        std::printf("p_sz_params[%d]=%llu\n",i,p_sz_params[i]);
        if(p_sz_params[i]<lo_bound)
        {
            p_sz_params[i] = lo_bound;
        }
        if(p_sz_params[i]>hi_bound)
        {
            p_sz_params[i] = hi_bound;
        }
    }
    atmos_refraction_state_t refraction_state = atmos_refraction_state_t(sz_params,true);
    char * ctor_name{gms::common::demangle(typeid(refraction_state).name(),status)};
    if(status==0 && ctor_name != NULL)
    {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
    }
    else
    {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(refraction_state).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
    }
   
    refraction_state.m_data_1f714.info_size_alignment();
    refraction_state.m_data_beta.info_size_alignment();
    refraction_state.m_data_Bf61b.info_size_alignment();
    refraction_state.m_data_D1.info_size_alignment();
    refraction_state.m_data_d.info_size_alignment();
    refraction_state.m_data_deln0.info_size_alignment();
    refraction_state.m_data_delnA.info_size_alignment();
    refraction_state.m_data_delnb.info_size_alignment();
    refraction_state.m_data_delnc.info_size_alignment();
    refraction_state.m_data_delnh.info_size_alignment();
    refraction_state.m_data_dn0.info_size_alignment();
    refraction_state.m_data_dndr.info_size_alignment();
    std::printf("[UNIT-TEST]: %s ---> HAS ENDED\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_ctor_1_non_omp_alloc();

void unit_test_ctor_1_non_omp_alloc()
{
    using namespace gms::atmos;
    std::printf("[UNIT-TEST]: %s ---> HAS BEGUN\n", __PRETTY_FUNCTION__);
    constexpr std::size_t lo_bound{1000ull};
    constexpr std::size_t hi_bound{10000ull};
    const std::int32_t    n_params{114};
    std::clock_t seed_alloc;
    atmos_refraction_size_params_t sz_params;
    std::int32_t status;
    seed_alloc = std::clock();
    auto rand_alloc{std::bind(std::uniform_int_distribution<std::size_t>(lo_bound,hi_bound),std::mt19937(seed_alloc))};
    std::size_t * __restrict__ p_sz_params{nullptr};
    p_sz_params = (std::size_t*)&sz_params;
    for(std::int32_t i{0}; i != n_params; ++i) 
    {
        p_sz_params[i] = rand_alloc();
        std::printf("p_sz_params[%d]=%llu\n",i,p_sz_params[i]);
        if(p_sz_params[i]<lo_bound)
        {
            p_sz_params[i] = lo_bound;
        }
        if(p_sz_params[i]>hi_bound)
        {
            p_sz_params[i] = hi_bound;
        }
    }
    
    atmos_refraction_state_t refraction_state = atmos_refraction_state_t(sz_params,false);
    char * ctor_name{gms::common::demangle(typeid(refraction_state).name(),status)};
    if(status==0 && ctor_name != NULL)
    {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
    }
    else
    {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(refraction_state).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
    }
    refraction_state.m_data_1f714.info_size_alignment();
    refraction_state.m_data_beta.info_size_alignment();
    refraction_state.m_data_Bf61b.info_size_alignment();
    refraction_state.m_data_D1.info_size_alignment();
    refraction_state.m_data_d.info_size_alignment();
    refraction_state.m_data_deln0.info_size_alignment();
    refraction_state.m_data_delnA.info_size_alignment();
    refraction_state.m_data_delnb.info_size_alignment();
    refraction_state.m_data_delnc.info_size_alignment();
    refraction_state.m_data_delnh.info_size_alignment();
    refraction_state.m_data_dn0.info_size_alignment();
    refraction_state.m_data_dndr.info_size_alignment();
   std::printf("[UNIT-TEST]: %s ---> HAS ENDED\n", __PRETTY_FUNCTION__);
}

int main()
{
    unit_test_ctor_1_non_omp_alloc();
    unit_test_ctor_1_omp_alloc();
    return 0;
}