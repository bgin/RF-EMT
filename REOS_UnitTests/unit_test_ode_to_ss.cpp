#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <string>
#include "GMS_ode_to_state_space.h"

/*
   icpc -o unit_test_ode_to_ss -O3 -I /home/bernard/Software-Development/Cpp-projects/RF_EMT/RF-EMT_src/eigen-5.0.0 -fp-model fast=2 -ftz -std=c++17 -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_ode_to_state_space.h unit_test_ode_to_ss.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -O3 -std=c++17 -march=skylake-avx512 -mavx512f  -I /home/bernard/Software-Development/Cpp-projects/RF_EMT/RF-EMT_src/eigen-5.0.0  -falign-functions=32 -fopenmp GMS_config.h GMS_ode_to_state_space.h unit_test_ode_to_ss.cpp

*/

__ATTR_HOT__ 
__ATTR_ALIGN__(32)
void unit_test_ode_to_ss();

void unit_test_ode_to_ss()
{
    using namespace gms::state_space;
    [[nodiscard]] std::int32_t printf_retv;
    printf_retv = std::printf("[UNIT-TEST:] -- of function=%s -- STARTED at day=%s,time=%s!!.\n",__func__,__DATE__,__TIME__);
    printf_retv = std::printf("===============================================================\n");
    printf_retv = std::printf("               ODE_TO_STATE_SPACE class constructors tests\n");
    printf_retv = std::printf("===============================================================\n");

    // Example: Third-order controlled plant ODE
    // y''' + 4y'' + 6y' + 2y = 5u' + 3u
    // Leading coefficient of y''' is 1.
    printf_retv = std::printf("[1] ============== Test Default Constructor =============\n");
    {
        ode_to_state_space_t default_ctor;
        printf_retv = std::printf("m_norder=%llu\n", default_ctor.m_norder);
        printf_retv = std::printf("m_A.rows=%llu\n", default_ctor.m_A.rows());
        printf_retv = std::printf("m_A.cols=%llu\n", default_ctor.m_A.cols());
        printf_retv = std::printf("m_B.rows=%llu\n", default_ctor.m_B.rows());
        printf_retv = std::printf("m_C.rows=%llu\n", default_ctor.m_C.rows());
        printf_retv = std::printf("m_D=%.17f\n",     default_ctor.m_D);
    }

    printf_retv = std::printf("[1] ============== Test MultiArg Constructor =============\n");
    {
        ode_to_state_space_t multiarg_ctor((Eigen::Index)3);
        printf_retv = std::printf("m_norder=%llu\n", multiarg_ctor.m_norder);
        printf_retv = std::printf("m_A.rows=%llu\n", multiarg_ctor.m_A.rows());
        printf_retv = std::printf("m_A.cols=%llu\n", multiarg_ctor.m_A.cols());
        printf_retv = std::printf("m_B.rows=%llu\n", multiarg_ctor.m_B.rows());
        printf_retv = std::printf("m_C.rows=%llu\n", multiarg_ctor.m_C.rows());
        printf_retv = std::printf("m_D=%.17f\n",     multiarg_ctor.m_D);
    }

    printf_retv = std::printf("[1] ============== Test Copy/Move Constructors =============\n");
    {
        ode_to_state_space_t multiarg_ctor((Eigen::Index)3);
        ode_to_state_space_t cpy_ctor(multiarg_ctor);
        printf_retv = std::printf("m_norder=%llu\n", cpy_ctor.m_norder);
        printf_retv = std::printf("m_A.rows=%llu\n", cpy_ctor.m_A.rows());
        printf_retv = std::printf("m_A.cols=%llu\n", cpy_ctor.m_A.cols());
        printf_retv = std::printf("m_B.rows=%llu\n", cpy_ctor.m_B.rows());
        printf_retv = std::printf("m_C.rows=%llu\n", cpy_ctor.m_C.rows());
        printf_retv = std::printf("m_D=%.17f\n",     cpy_ctor.m_D);
        ode_to_state_space_t move_ctor(ode_to_state_space_t((Eigen::Index)3));
        printf_retv = std::printf("m_norder=%llu\n", move_ctor.m_norder);
        printf_retv = std::printf("m_A.rows=%llu\n", move_ctor.m_A.rows());
        printf_retv = std::printf("m_A.cols=%llu\n", move_ctor.m_A.cols());
        printf_retv = std::printf("m_B.rows=%llu\n", move_ctor.m_B.rows());
        printf_retv = std::printf("m_C.rows=%llu\n", move_ctor.m_C.rows());
        printf_retv = std::printf("m_D=%.17f\n",     move_ctor.m_D);
    }

}

int main()
{   
    (void)unit_test_ode_to_ss();
    return 0;
}