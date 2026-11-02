
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_ofdm_signal_modulator.h"

/*
   icpc -o unit_test_ofdm_signal_modulator -O3 -fp-model fast=2 -std=c++17 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_half.h GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_ofdm_signal_modulator.h unit_test_ofdm_signal_modulator.cpp
   ASM: 
   icpc  -S  -O3 -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_half.h GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_cephes_sin_cos.h GMS_ofdm_signal_modulator.h unit_test_ofdm_signal_modulator.cpp

*/



__ATTR_HOT__ 
__ATTR_ALIGN__(32)
bool unit_test_ofdm_signal_modulator_simple();

bool unit_test_ofdm_signal_modulator_simple()
{
    bool ret_status = false;
    using namespace gms::radiolocation;
    using namespace gms;
    [[nodiscard]] std::int32_t printf_ret = -9999;
    printf_ret =  printf("[UNIT-TEST]: function=%s -- **STARTED**\n", __PRETTY_FUNCTION__);
    printf_ret =  printf("[UNIT-TEST]: [1] -- Initializing AVX-512 Manual (partly) Vectorized 5G OFDM Generator Core...\n");
    const std::size_t   USERS   = 600ull;
    const std::size_t   FFT     = 1024ull;
    const std::size_t   CP_LEN  = 72;
    ofdm_signal_modulator_f32_t<PrefetchMode::Enabled> generator(USERS,FFT,CP_LEN);
    printf_ret =  printf("[UNIT-TEST]: [2] -- Mapping user-channels to QAM256 constellation...\n");
    ret_status = generator.generate_channels_randomly(QamScheme::QAM256);
    if(!ret_status) 
    { 
        printf_ret =  printf("[UNIT-TEST]: [2] -- Mapping user-channels to QAM256 constellation: ***FAILED***\n");
        return (ret_status);
    }
    printf_ret =  printf("[UNIT-TEST]: [3] -- Executing IDFT Modulation Core & Cyclic Prefix Insertion...\n");
    ret_status = false;
    ret_status = generator.compute_idft_and_apply_cycle_prefix();
    if(!ret_status) 
    { 
        printf_ret =  printf("[UNIT-TEST]: [3] -- Executing IDFT Modulation Core & Cyclic Prefix Insertion: ***FAILED***\n");
        return (ret_status);
    }
    printf_ret =  printf("[UNIT-TEST]: [4] -- Execution complete! Total Signal Samples Generated: %llu\n",generator.get_total_output_size());
    printf_ret =  printf("[UNIT-TEST]: function=%s -- **ENDED**\n", __PRETTY_FUNCTION__);
    ret_status = true;
    return (ret_status);
}


int main()
{   
    bool test_status = unit_test_ofdm_signal_modulator_simple();
    if(!test_status) {std::printf("[TEST-RUNNER]: unit_test_ofdm_signal_modulator_simple: ***FAILED***!!");}
    return 0;
}