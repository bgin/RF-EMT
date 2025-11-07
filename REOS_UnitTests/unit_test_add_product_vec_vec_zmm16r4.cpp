
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <immintrin.h>
#include <string>
#include "GMS_malloc.h"
#include "GMS_add_product_vec_zmm16r4.h"

/*
   icpc -o unit_test_add_product_vec_vec_zmm16r4_a  -O3 -fp-model fast=2 -ftz -ggdb -ipo -qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5 -fopenmp \
   GMS_config.h GMS_malloc.h GMS_add_product_vec_zmm16r4.h GMS_add_product_vec_zmm16r4.cpp unit_test_add_product_vec_vec_zmm16r4.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel --O3 qopt-zmm-usage=high -march=skylake-avx512 -mavx512f -fopenmp -falign-functions=32 GMS_config.h GMS_malloc.h GMS_add_product_vec_zmm16r4.h GMS_add_product_vec_zmm16r4.cpp unit_test_add_product_vec_vec_zmm16r4.cpp

*/

int32_t funi_distro_idx = 0;
int32_t fnorm_distr_idx = 0;

__attribute__((hot))
__attribute__((noinline))
__attribute__((aligned(32)))
void unit_test_vaddpv_f32_2x512_u16x_a();

__attribute__((hot))
__attribute__((noinline))
__attribute__((aligned(32)))
void unit_test_vaddpv_f32_2x512_u16x_a()
{
    using namespace gms::common;
    constexpr std::size_t nelems{456748ull};
    constexpr std::size_t nbytes{reinterpret_cast<std::size_t>(sizeof(float)*nelems)};
    constexpr unsigned long long RDTSCP_LAT{42ull};
    std::string fname1 = "UNIT-TEST_uni_distro_vaddpv_f32_2x512_u16x_a_" + std::to_string(funi_distro_idx++)+".csv";
    std::string fname2 = "UNIT-TEST_norm_distro_vaddpv_f32_2x512_u16x_a_" + std::to_string(fnorm_distr_idx++)+".csv";
    float * __restrict__ vec_x{NULL};
    float * __restrict__ vec_y{NULL};
    float * __restrict__ vec_z{NULL};
    FILE  * __restrict__ fp{NULL};
    std::clock_t seed_d{0ull};
    unsigned long long seed_vec_x{0ull};
    unsigned long long seed_vec_y{0ull};
    unsigned long long seed_vec_z{0ull};
    unsigned long long start{0ull};
    unsigned long long end{0ull};
    unsigned long long start_c{0ull};
    unsigned long long end_c{0ull};
    unsigned long long elapsed_tsc{0ull};
    uint32_t tsc_aux{9999};
    float rx,ry,rz;
    int32_t which{-1};
    int32_t idx{};
    bool  vec_x_neq_0{};
    bool  vec_y_neq_0{};
    bool  vec_z_neq_0{};
    seed_d = std::clock();
    auto rand_d{std::bind(std::uniform_int_distribution<int32_t>(0,1),std::mt19937(seed_d))};
    which = rand_d();

    switch (which)
    {
        case 0 :
        {
            std::uniform_real_distribution<float> uniform_distro;
            const char * __restrict__ fout{fname1.c_str()};
            char * distro_name{NULL};
            int32_t status{9999};
            distro_name = demangle(typeid(uniform_distro).name(),status);
            if(distro_name != NULL && status == 0)
            {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
            }
            else
            {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n", typeid(uniform_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
            }
            vec_x = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            vec_x_neq_0 = (vec_x != NULL) && (nbytes > 0ull);
            vec_y = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            vec_y_neq_0 = (vec_y != NULL) && (nbytes > 0ull);
            vec_z = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            vec_z_neq_0 = (vec_z != NULL) && (nbytes > 0ull);
            seed_vec_x = __rdtsc();
            auto rand_vec_x{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                     std::mt19937(seed_vec_x))};
            seed_vec_y = __rdtsc();
            auto rand_vec_y{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                     std::mt19937(seed_vec_y))};
            seed_vec_z = __rdtsc();
            auto rand_vec_z{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                     std::mt19937(seed_vec_z))};                         
            std::size_t __i;
            for(__i = 0; __i != ROUND_TO_EIGHT(nelems,8); __i += 8)
            {
                rx = rand_vec_x();
                vec_x[__i+0ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+0ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+0ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+1ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+1ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+1ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+2ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+2ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+2ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+3ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+3ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+3ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+4ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+4ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+4ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+5ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+5ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+5ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+6ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+6ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+6ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+7ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+7ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+7ull] = rz;
            }
            for(;__i != nelems; ++__i)
            {
                rx = rand_vec_x();
                vec_x[__i] = rx;
                ry = rand_vec_y();
                vec_y[__i] = ry;
                rz = rand_vec_z();
                vec_z[__i] = rz;
            }

            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: vaddpv_f32_2x512_u16x_a\n");
            start   = __rdtscp(&tsc_aux);
            std::int32_t call_state = gms::math::vaddpv_f32_2x512_u16x_a(vec_x,vec_y,vec_z,nelems);
            end     = __rdtscp(&tsc_aux);
            start_c = start-RDTSCP_LAT;
            end_c   = end-RDTSCP_LAT;
            elapsed_tsc = end_c-start_c;
            printf("[UNIT-TEST]: vaddpv_f32_2x512_u16x_a -- TSC=%llu,TSC_AUX=%d,call_state=%d\n",elapsed_tsc,tsc_aux,call_state);
            fp = fopen(fout,"w+");
            if(!fp)
            {
                 std::perror("fopen failed to open a file -- TERMINATING!!");
                 if(vec_y_neq_0) gms_mm_free(vec_y);
                 if(vec_x_neq_0) gms_mm_free(vec_x);
                 if(vec_z_neq_0) gms_mm_free(vec_z);
                 std::terminate();
            }
            
            for(std::size_t __i{0}; __i != nelems; ++__i)
            {
                fprintf(fp,"%.7f\n",vec_z[__i]);
            }
            if(vec_z_neq_0) gms_mm_free(vec_z);
            if(vec_y_neq_0) gms_mm_free(vec_y);
            if(vec_x_neq_0) gms_mm_free(vec_x);
            fclose(fp);
            printf("[UNIT-TEST]: -- END: vaddpv_f32_2x512_u16x_a\n");
        }
        break;
        case 1 :
        {
            std::normal_distribution<float> normal_distro;
            const char * __restrict__ fout{fname2.c_str()};
            char * distro_name{NULL};
            int32_t status{9999};
            distro_name = demangle(typeid(normal_distro).name(),status);
            if(distro_name != NULL && status == 0)
            {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
            }
            else
            {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n", typeid(normal_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
            }

            vec_x = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            vec_y = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            vec_z = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            seed_vec_x = __rdtsc();
            auto rand_vec_x{std::bind(std::normal_distribution<float>(0.0f,0.2f),
                                     std::mt19937(seed_vec_x))};
            seed_vec_y = __rdtsc();
            auto rand_vec_y{std::bind(std::normal_distribution<float>(0.0f,1.0f),
                                     std::mt19937(seed_vec_y))};
            seed_vec_z = __rdtsc();
            auto rand_vec_z{std::bind(std::normal_distribution<float>(0.0f,1.0f),
                                     std::mt19937(seed_vec_z))};
            std::size_t __i;
            for(__i = 0; __i != ROUND_TO_EIGHT(nelems,8); __i += 8)
            {
                rx = rand_vec_x();
                vec_x[__i+0ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+0ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+0ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+1ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+1ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+1ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+2ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+2ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+2ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+3ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+3ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+3ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+4ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+4ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+4ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+5ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+5ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+5ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+6ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+6ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+6ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+7ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+7ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+7ull] = rz;
            }
            for(;__i != nelems; ++__i)
            {
                rx = rand_vec_x();
                vec_x[__i] = rx;
                ry = rand_vec_y();
                vec_y[__i] = ry;
                rz = rand_vec_z();
                vec_z[__i] = rz;
            }

            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: vaddpv_f32_2x512_u16x_a\n");
            start  = __rdtscp(&tsc_aux);
            std::int32_t call_state = gms::math::vaddpv_f32_2x512_u16x_a(vec_x,vec_y,vec_z,nelems);
            end     = __rdtscp(&tsc_aux);
            start_c = start-RDTSCP_LAT;
            end_c   = end-RDTSCP_LAT;
            elapsed_tsc = end_c-start_c;
            printf("[UNIT-TEST]: vaddpv_f32_2x512_u16x_a -- TSC=%llu,TSC_AUX=%d,call_state=%d\n",elapsed_tsc,tsc_aux,call_state);
            fp = fopen(fout,"w+");
            if(!fp)
            {
                 std::perror("fopen failed to open a file -- TERMINATING!!");
                 if(vec_z_neq_0) gms_mm_free(vec_z);
                 if(vec_y_neq_0) gms_mm_free(vec_y);
                 if(vec_x_neq_0) gms_mm_free(vec_x);
                 std::terminate();
            }
            
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                fprintf(fp,"%.7f\n",vec_z[__i]);
            }
            if(vec_z_neq_0) gms_mm_free(vec_z);
            if(vec_y_neq_0) gms_mm_free(vec_y);
            if(vec_x_neq_0) gms_mm_free(vec_x);
            fclose(fp);
            printf("[UNIT-TEST]: -- END: vaddpv_f32_2x512_u16x_a\n");
        }
        break;
        default : 
                 printf("Invalid switch variable=%d\n",which);
                 return;
    } 
}


__attribute__((hot))
__attribute__((noinline))
__attribute__((aligned(32)))
void unit_test_vaddpv_f32_2x512_u10x_a();

__attribute__((hot))
__attribute__((noinline))
__attribute__((aligned(32)))
void unit_test_vaddpv_f32_2x512_u10x_a()
{
    using namespace gms::common;
    constexpr std::size_t nelems{456748ull};
    constexpr std::size_t nbytes{reinterpret_cast<std::size_t>(sizeof(float)*nelems)};
    constexpr unsigned long long RDTSCP_LAT{42ull};
    std::string fname1 = "UNIT-TEST_uni_distro_vaddpv_f32_2x512_u10x_a_" + std::to_string(funi_distro_idx++)+".csv";
    std::string fname2 = "UNIT-TEST_norm_distro_vaddpv_f32_2x512_u10x_a_" + std::to_string(fnorm_distr_idx++)+".csv";
    float * __restrict__ vec_x{NULL};
    float * __restrict__ vec_y{NULL};
    float * __restrict__ vec_z{NULL};
    FILE  * __restrict__ fp{NULL};
    std::clock_t seed_d{0ull};
    unsigned long long seed_vec_x{0ull};
    unsigned long long seed_vec_y{0ull};
    unsigned long long seed_vec_z{0ull};
    unsigned long long start{0ull};
    unsigned long long end{0ull};
    unsigned long long start_c{0ull};
    unsigned long long end_c{0ull};
    unsigned long long elapsed_tsc{0ull};
    uint32_t tsc_aux{9999};
    float rx,ry,rz;
    int32_t which{-1};
    int32_t idx{};
    bool  vec_x_neq_0{};
    bool  vec_y_neq_0{};
    bool  vec_z_neq_0{};
    seed_d = std::clock();
    auto rand_d{std::bind(std::uniform_int_distribution<int32_t>(0,1),std::mt19937(seed_d))};
    which = rand_d();

    switch (which)
    {
        case 0 :
        {
            std::uniform_real_distribution<float> uniform_distro;
            const char * __restrict__ fout{fname1.c_str()};
            char * distro_name{NULL};
            int32_t status{9999};
            distro_name = demangle(typeid(uniform_distro).name(),status);
            if(distro_name != NULL && status == 0)
            {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
            }
            else
            {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n", typeid(uniform_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
            }
            vec_x = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            vec_x_neq_0 = (vec_x != NULL) && (nbytes > 0ull);
            vec_y = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            vec_y_neq_0 = (vec_y != NULL) && (nbytes > 0ull);
            vec_z = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            vec_z_neq_0 = (vec_z != NULL) && (nbytes > 0ull);
            seed_vec_x = __rdtsc();
            auto rand_vec_x{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                     std::mt19937(seed_vec_x))};
            seed_vec_y = __rdtsc();
            auto rand_vec_y{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                     std::mt19937(seed_vec_y))};
            seed_vec_z = __rdtsc();
            auto rand_vec_z{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                     std::mt19937(seed_vec_z))};                         
            std::size_t __i;
            for(__i = 0; __i != ROUND_TO_EIGHT(nelems,8); __i += 8)
            {
                rx = rand_vec_x();
                vec_x[__i+0ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+0ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+0ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+1ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+1ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+1ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+2ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+2ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+2ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+3ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+3ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+3ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+4ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+4ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+4ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+5ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+5ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+5ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+6ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+6ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+6ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+7ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+7ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+7ull] = rz;
            }
            for(;__i != nelems; ++__i)
            {
                rx = rand_vec_x();
                vec_x[__i] = rx;
                ry = rand_vec_y();
                vec_y[__i] = ry;
                rz = rand_vec_z();
                vec_z[__i] = rz;
            }

            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: vaddpv_f32_2x512_u10x_a\n");
            start   = __rdtscp(&tsc_aux);
            std::int32_t call_state = gms::math::vaddpv_f32_2x512_u10x_a(vec_x,vec_y,vec_z,nelems);
            end     = __rdtscp(&tsc_aux);
            start_c = start-RDTSCP_LAT;
            end_c   = end-RDTSCP_LAT;
            elapsed_tsc = end_c-start_c;
            printf("[UNIT-TEST]: vaddpv_f32_2x512_u10x_a -- TSC=%llu,TSC_AUX=%d,call_state=%d\n",elapsed_tsc,tsc_aux,call_state);
            fp = fopen(fout,"w+");
            if(!fp)
            {
                 std::perror("fopen failed to open a file -- TERMINATING!!");
                 if(vec_y_neq_0) gms_mm_free(vec_y);
                 if(vec_x_neq_0) gms_mm_free(vec_x);
                 if(vec_z_neq_0) gms_mm_free(vec_z);
                 std::terminate();
            }
            
            for(std::size_t __i{0}; __i != nelems; ++__i)
            {
                fprintf(fp,"%.7f\n",vec_z[__i]);
            }
            if(vec_z_neq_0) gms_mm_free(vec_z);
            if(vec_y_neq_0) gms_mm_free(vec_y);
            if(vec_x_neq_0) gms_mm_free(vec_x);
            fclose(fp);
            printf("[UNIT-TEST]: -- END: vaddpv_f32_2x512_u10x_a\n");
        }
        break;
        case 1 :
        {
            std::normal_distribution<float> normal_distro;
            const char * __restrict__ fout{fname2.c_str()};
            char * distro_name{NULL};
            int32_t status{9999};
            distro_name = demangle(typeid(normal_distro).name(),status);
            if(distro_name != NULL && status == 0)
            {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
            }
            else
            {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n", typeid(normal_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
            }

            vec_x = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            vec_y = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            vec_z = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            seed_vec_x = __rdtsc();
            auto rand_vec_x{std::bind(std::normal_distribution<float>(0.0f,0.2f),
                                     std::mt19937(seed_vec_x))};
            seed_vec_y = __rdtsc();
            auto rand_vec_y{std::bind(std::normal_distribution<float>(0.0f,1.0f),
                                     std::mt19937(seed_vec_y))};
            seed_vec_z = __rdtsc();
            auto rand_vec_z{std::bind(std::normal_distribution<float>(0.0f,1.0f),
                                     std::mt19937(seed_vec_z))};
            std::size_t __i;
            for(__i = 0; __i != ROUND_TO_EIGHT(nelems,8); __i += 8)
            {
                rx = rand_vec_x();
                vec_x[__i+0ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+0ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+0ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+1ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+1ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+1ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+2ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+2ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+2ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+3ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+3ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+3ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+4ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+4ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+4ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+5ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+5ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+5ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+6ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+6ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+6ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+7ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+7ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+7ull] = rz;
            }
            for(;__i != nelems; ++__i)
            {
                rx = rand_vec_x();
                vec_x[__i] = rx;
                ry = rand_vec_y();
                vec_y[__i] = ry;
                rz = rand_vec_z();
                vec_z[__i] = rz;
            }

            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: vaddpv_f32_2x512_u10x_a\n");
            start  = __rdtscp(&tsc_aux);
            std::int32_t call_state = gms::math::vaddpv_f32_2x512_u10x_a(vec_x,vec_y,vec_z,nelems);
            end     = __rdtscp(&tsc_aux);
            start_c = start-RDTSCP_LAT;
            end_c   = end-RDTSCP_LAT;
            elapsed_tsc = end_c-start_c;
            printf("[UNIT-TEST]: vaddpv_f32_2x512_u10x_a -- TSC=%llu,TSC_AUX=%d,call_state=%d\n",elapsed_tsc,tsc_aux,call_state);
            fp = fopen(fout,"w+");
            if(!fp)
            {
                 std::perror("fopen failed to open a file -- TERMINATING!!");
                 if(vec_z_neq_0) gms_mm_free(vec_z);
                 if(vec_y_neq_0) gms_mm_free(vec_y);
                 if(vec_x_neq_0) gms_mm_free(vec_x);
                 std::terminate();
            }
            
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                fprintf(fp,"%.7f\n",vec_z[__i]);
            }
            if(vec_z_neq_0) gms_mm_free(vec_z);
            if(vec_y_neq_0) gms_mm_free(vec_y);
            if(vec_x_neq_0) gms_mm_free(vec_x);
            fclose(fp);
            printf("[UNIT-TEST]: -- END: vaddpv_f32_2x512_u10x_a\n");
        }
        break;
        default : 
                 printf("Invalid switch variable=%d\n",which);
                 return;
    } 
}


__attribute__((hot))
__attribute__((noinline))
__attribute__((aligned(32)))
void unit_test_vaddpv_f32_2x512_u6x_a();

__attribute__((hot))
__attribute__((noinline))
__attribute__((aligned(32)))
void unit_test_vaddpv_f32_2x512_u6x_a()
{
    using namespace gms::common;
    constexpr std::size_t nelems{456748ull};
    constexpr std::size_t nbytes{reinterpret_cast<std::size_t>(sizeof(float)*nelems)};
    constexpr unsigned long long RDTSCP_LAT{42ull};
    std::string fname1 = "UNIT-TEST_uni_distro_vaddpv_f32_2x512_u6x_a_" + std::to_string(funi_distro_idx++)+".csv";
    std::string fname2 = "UNIT-TEST_norm_distro_vaddpv_f32_2x512_u6x_a_" + std::to_string(fnorm_distr_idx++)+".csv";
    float * __restrict__ vec_x{NULL};
    float * __restrict__ vec_y{NULL};
    float * __restrict__ vec_z{NULL};
    FILE  * __restrict__ fp{NULL};
    std::clock_t seed_d{0ull};
    unsigned long long seed_vec_x{0ull};
    unsigned long long seed_vec_y{0ull};
    unsigned long long seed_vec_z{0ull};
    unsigned long long start{0ull};
    unsigned long long end{0ull};
    unsigned long long start_c{0ull};
    unsigned long long end_c{0ull};
    unsigned long long elapsed_tsc{0ull};
    uint32_t tsc_aux{9999};
    float rx,ry,rz;
    int32_t which{-1};
    int32_t idx{};
    bool  vec_x_neq_0{};
    bool  vec_y_neq_0{};
    bool  vec_z_neq_0{};
    seed_d = std::clock();
    auto rand_d{std::bind(std::uniform_int_distribution<int32_t>(0,1),std::mt19937(seed_d))};
    which = rand_d();

    switch (which)
    {
        case 0 :
        {
            std::uniform_real_distribution<float> uniform_distro;
            const char * __restrict__ fout{fname1.c_str()};
            char * distro_name{NULL};
            int32_t status{9999};
            distro_name = demangle(typeid(uniform_distro).name(),status);
            if(distro_name != NULL && status == 0)
            {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
            }
            else
            {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n", typeid(uniform_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
            }
            vec_x = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            vec_x_neq_0 = (vec_x != NULL) && (nbytes > 0ull);
            vec_y = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            vec_y_neq_0 = (vec_y != NULL) && (nbytes > 0ull);
            vec_z = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            vec_z_neq_0 = (vec_z != NULL) && (nbytes > 0ull);
            seed_vec_x = __rdtsc();
            auto rand_vec_x{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                     std::mt19937(seed_vec_x))};
            seed_vec_y = __rdtsc();
            auto rand_vec_y{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                     std::mt19937(seed_vec_y))};
            seed_vec_z = __rdtsc();
            auto rand_vec_z{std::bind(std::uniform_real_distribution<float>(-3.14159265358979323846264338328f,3.14159265358979323846264338328f),
                                     std::mt19937(seed_vec_z))};                         
            std::size_t __i;
            for(__i = 0; __i != ROUND_TO_EIGHT(nelems,8); __i += 8)
            {
                rx = rand_vec_x();
                vec_x[__i+0ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+0ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+0ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+1ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+1ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+1ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+2ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+2ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+2ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+3ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+3ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+3ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+4ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+4ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+4ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+5ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+5ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+5ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+6ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+6ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+6ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+7ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+7ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+7ull] = rz;
            }
            for(;__i != nelems; ++__i)
            {
                rx = rand_vec_x();
                vec_x[__i] = rx;
                ry = rand_vec_y();
                vec_y[__i] = ry;
                rz = rand_vec_z();
                vec_z[__i] = rz;
            }

            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: vaddpv_f32_2x512_u6x_a\n");
            start   = __rdtscp(&tsc_aux);
            std::int32_t call_state = gms::math::vaddpv_f32_2x512_u6x_a(vec_x,vec_y,vec_z,nelems);
            end     = __rdtscp(&tsc_aux);
            start_c = start-RDTSCP_LAT;
            end_c   = end-RDTSCP_LAT;
            elapsed_tsc = end_c-start_c;
            printf("[UNIT-TEST]: vaddpv_f32_2x512_u6x_a -- TSC=%llu,TSC_AUX=%d,call_state=%d\n",elapsed_tsc,tsc_aux,call_state);
            fp = fopen(fout,"w+");
            if(!fp)
            {
                 std::perror("fopen failed to open a file -- TERMINATING!!");
                 if(vec_y_neq_0) gms_mm_free(vec_y);
                 if(vec_x_neq_0) gms_mm_free(vec_x);
                 if(vec_z_neq_0) gms_mm_free(vec_z);
                 std::terminate();
            }
            
            for(std::size_t __i{0}; __i != nelems; ++__i)
            {
                fprintf(fp,"%.7f\n",vec_z[__i]);
            }
            if(vec_z_neq_0) gms_mm_free(vec_z);
            if(vec_y_neq_0) gms_mm_free(vec_y);
            if(vec_x_neq_0) gms_mm_free(vec_x);
            fclose(fp);
            printf("[UNIT-TEST]: -- END: vaddpv_f32_2x512_u6x_a\n");
        }
        break;
        case 1 :
        {
            std::normal_distribution<float> normal_distro;
            const char * __restrict__ fout{fname2.c_str()};
            char * distro_name{NULL};
            int32_t status{9999};
            distro_name = demangle(typeid(normal_distro).name(),status);
            if(distro_name != NULL && status == 0)
            {
                   printf("[UNIT-TEST]: Instantiated distribution of type: %s\n", distro_name);
                   gms::common::gms_mm_free(distro_name);
            }
            else
            {
                   printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n", typeid(normal_distro).name());
                   if(distro_name != NULL) gms::common::gms_mm_free(distro_name);
            }

            vec_x = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            vec_y = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            vec_z = reinterpret_cast<float * __restrict__>(gms_mm_malloc(nbytes,64ull));
            seed_vec_x = __rdtsc();
            auto rand_vec_x{std::bind(std::normal_distribution<float>(0.0f,0.2f),
                                     std::mt19937(seed_vec_x))};
            seed_vec_y = __rdtsc();
            auto rand_vec_y{std::bind(std::normal_distribution<float>(0.0f,1.0f),
                                     std::mt19937(seed_vec_y))};
            seed_vec_z = __rdtsc();
            auto rand_vec_z{std::bind(std::normal_distribution<float>(0.0f,1.0f),
                                     std::mt19937(seed_vec_z))};
            std::size_t __i;
            for(__i = 0; __i != ROUND_TO_EIGHT(nelems,8); __i += 8)
            {
                rx = rand_vec_x();
                vec_x[__i+0ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+0ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+0ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+1ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+1ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+1ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+2ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+2ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+2ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+3ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+3ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+3ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+4ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+4ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+4ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+5ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+5ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+5ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+6ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+6ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+6ull] = rz;
                rx = rand_vec_x();
                vec_x[__i+7ull] = rx;
                ry = rand_vec_y();
                vec_y[__i+7ull] = ry;
                rz = rand_vec_z();
                vec_z[__i+7ull] = rz;
            }
            for(;__i != nelems; ++__i)
            {
                rx = rand_vec_x();
                vec_x[__i] = rx;
                ry = rand_vec_y();
                vec_y[__i] = ry;
                rz = rand_vec_z();
                vec_z[__i] = rz;
            }

            printf("[UNIT-TEST]: -- End of data initialization\n");
            printf("[UNIT-TEST]: -- START: vaddpv_f32_2x512_u6x_a\n");
            start  = __rdtscp(&tsc_aux);
            std::int32_t call_state = gms::math::vaddpv_f32_2x512_u6x_a(vec_x,vec_y,vec_z,nelems);
            end     = __rdtscp(&tsc_aux);
            start_c = start-RDTSCP_LAT;
            end_c   = end-RDTSCP_LAT;
            elapsed_tsc = end_c-start_c;
            printf("[UNIT-TEST]: vaddpv_f32_2x512_u6x_a -- TSC=%llu,TSC_AUX=%d,call_state=%d\n",elapsed_tsc,tsc_aux,call_state);
            fp = fopen(fout,"w+");
            if(!fp)
            {
                 std::perror("fopen failed to open a file -- TERMINATING!!");
                 if(vec_z_neq_0) gms_mm_free(vec_z);
                 if(vec_y_neq_0) gms_mm_free(vec_y);
                 if(vec_x_neq_0) gms_mm_free(vec_x);
                 std::terminate();
            }
            
            for(int32_t __i{0}; __i != nelems; ++__i)
            {
                fprintf(fp,"%.7f\n",vec_z[__i]);
            }
            if(vec_z_neq_0) gms_mm_free(vec_z);
            if(vec_y_neq_0) gms_mm_free(vec_y);
            if(vec_x_neq_0) gms_mm_free(vec_x);
            fclose(fp);
            printf("[UNIT-TEST]: -- END: vaddpv_f32_2x512_u6x_a\n");
        }
        break;
        default : 
                 printf("Invalid switch variable=%d\n",which);
                 return;
    } 
}



int main()
{  
   unit_test_vaddpv_f32_2x512_u16x_a();
   unit_test_vaddpv_f32_2x512_u10x_a();
   unit_test_vaddpv_f32_2x512_u6x_a();
   return 0;
}