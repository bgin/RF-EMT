/*******************************************************************************
*
*
*  Copyright [2019 - 2023] [Intel Corporation]
* 
*  Licensed under the Apache License, Version 2.0 (the "License");
*  you may not use this file except in compliance with the License.
*  
*  You may obtain a copy of the License at
*  
*     http://www.apache.org/licenses/LICENSE-2.0 
*  
*  Unless required by applicable law or agreed to in writing, software 
*  distributed under the License is distributed on an "AS IS" BASIS, 
*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*  See the License for the specific language governing permissions and 
*  limitations under the License. 
*  
*  SPDX-License-Identifier: Apache-2.0 
*  
* 
*
*******************************************************************************/

/*MIT License
!Copyright (c) 2020 Bernard Gingold
!Permission is hereby granted, free of charge, to any person obtaining a copy
!of this software and associated documentation files (the "Software"), to deal
!in the Software without restriction, including without limitation the rights
!to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!copies of the Software, and to permit persons to whom the Software is
!furnished to do so, subject to the following conditions:
!The above copyright notice and this permission notice shall be included in all
!copies or substantial portions of the Software.
!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!SOFTWARE.
*/

#ifndef __GMS_LDPC_ENCODER_CYCLE_BIT_AVX512_H__
#define __GMS_LDPC_ENCODER_CYCLE_BIT_AVX512_H__ 080720260836

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"

namespace file_info 
{

     static const unsigned int GMS_LDPC_ENCODER_CYCLE_BIT_AVX512_MAJOR = 1;
     static const unsigned int GMS_LDPC_ENCODER_CYCLE_BIT_AVX512_MINOR = 1;
     static const unsigned int GMS_LDPC_ENCODER_CYCLE_BIT_AVX512_MICRO = 0;
     static const unsigned int GMS_LDPC_ENCODER_CYCLE_BIT_AVX512_FULLVER =
       1000U*GMS_LDPC_ENCODER_CYCLE_BIT_AVX512_MAJOR+100U*GMS_LDPC_ENCODER_CYCLE_BIT_AVX512_MINOR+
       10U*GMS_LDPC_ENCODER_CYCLE_BIT_AVX512_MICRO;
     static const char GMS_LDPC_ENCODER_CYCLE_BIT_AVX512_CREATION_DATE[] = "08-07-2026 08:38PM +00200 (WED 08 JUL 2026 GMT+2)";
     static const char GMS_LDPC_ENCODER_CYCLE_BIT_AVX512_BUILD_DATE[]    = __DATE__; 
     static const char GMS_LDPC_ENCODER_CYCLE_BIT_AVX512_BUILD_TIME[]    = __TIME__;
     static const char GMS_LDPC_ENCODER_CYCLE_BIT_AVX512_SYNOPSIS[]      = "LDPC Encoder cycle bit shift helper kernels";

}

namespace gms 
{

namespace l1_phy 
{

__ATTR_ALIGN__(64)
extern std::int32_t g_permute_table_288_to_384[4][32];
__ATTR_ALIGN__(64)
extern std::int16_t g_permute_table_144_to_256[8][32];
__ATTR_ALIGN__(64)
extern std::int16_t g_permute_table_to_128[8][32];

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool reorder_instr_layout>
__ATTR_ALWAYS_INLINE__
static inline 
__m512i 
cycle_bit_lshift_288_to_384_avx512(__m512i data,std::int16_t cycle_lshift,
                                   std::int16_t zc_size,std::int8_t zc_index,
                                   __m512i swap_Idx0) 
{
    __m512i x1,x2,x3,x4;
    const std::int16_t cycle_lshift1 = cycle_lshift >> 5;
    const std::int32_t cycle_lshift2 = cycle_lshift &0x1f;
    const bool cond_reorder_instr_layout = (static_cast<std::int32_t>(reorder_instr_layout) == 
                                            static_cast<std::int32_t>(true));
    _mm_prefetch((const char*)(g_permute_table_288_to_384[zc_index]+cycle_lshift1),_MM_HINT_T0);
    __m512i swap_Idx1 = _mm512_load_si512((void const*)(g_permute_table_288_to_384[zc_index]+cycle_lshift1));
    if constexpr(cond_reorder_instr_layout)
    {
       x1 = _mm512_permutex2var_epi32(data,swap_Idx1,data);
       x3 = _mm512_srli_epi32(x1,cycle_lshift2);
       x2 = _mm512_permutex2var_epi32(x1,swap_Idx0,x1);
       x4 = _mm512_slli_epi32(x2,32-cycle_lshift2);
       x1 = _mm512_or_epi64(x3,x4);
    }
    else 
    {
       x1 = _mm512_permutex2var_epi32(data,swap_Idx1,data);
       x2 = _mm512_permutex2var_epi32(x1,swap_Idx0,x1);
       x3 = _mm512_srli_epi32(x1,cycle_lshift2);
       x4 = _mm512_slli_epi32(x2,32-cycle_lshift2);
       x1 = _mm512_or_epi64(x3,x4);
    }
    return (x1);
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool reorder_instr_layout>
__ATTR_ALWAYS_INLINE__
static inline 
__m512i 
cycle_bit_lshift_144_to_256_avx512(__m512i data,std::int16_t cycle_lshift,
                                   std::int16_t zc_size,std::int8_t zc_index,
                                   __m512i swap_Idx0) 
{
    __m512i x1,x2,x3,x4;
    __m512i swap_Idx1;
    __m256i swap_Idx11;
    std::int32_t cycle_lshift2;
    std::int16_t cycle_lshift1;

    while(cycle_lshift>zc_size)
    {
         cycle_lshift -= zc_size;
    }
    cycle_lshift1 = cycle_lshift >> 4;
    cycle_lshift2 = cycle_lshift & 0xf;

    if(zc_size > 128)
    {
        _mm_prefetch((const char*)(g_permute_table_144_to_256[zc_index]+cycle_lshift1),_MM_HINT_T0);
        swap_Idx11 = _mm256_load_si256((__m256i const*)(g_permute_table_144_to_256[zc_index]+cycle_lshift1));
    }
    else 
    {
        _mm_prefetch((const char*)(g_permute_table_to_128[zc_index]+cycle_lshift1),_MM_HINT_T0);
        swap_Idx11 = _mm256_load_si256((__m256i const*)(g_permute_table_to_128[zc_index]+cycle_lshift1));
    }
    swap_Idx1 = _mm512_broadcast_i32x8(swap_Idx11);
    swap_Idx1 = _mm512_mask_add_epi16(swapIdx1,0xffff0000,swap_Idx1,_mm512_set1_epi16(16));
    const bool cond_reorder_instr_layout = (static_cast<std::int32_t>(reorder_instr_layout) == 
                                            static_cast<std::int32_t>(true));
    if constexpr(cond_reorder_instr_layout)
    {
       x1 = _mm512_permutex2var_epi16(data,swap_Idx1,data);
       x3 = _mm512_srli_epi16(x1,cycle_lshift2);
       x2 = _mm512_permutex2var_epi16(x1,swap_Idx0,x1);
       x4 = _mm512_slli_epi16(x2,16-cycle_lshift2);
       x1 = _mm512_or_epi32(x3,x4);
    }
    else 
    {
       x1 = _mm512_permutex2var_epi16(data,swap_Idx1,data);
       x2 = _mm512_permutex2var_epi16(x1,swap_Idx0,x1);
       x3 = _mm512_srli_epi16(x1,cycle_lshift2);
       x4 = _mm512_slli_epi16(x2,32-cycle_lshift2);
       x1 = _mm512_or_epi32(x3,x4);
    }
    return (x1);
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
__ATTR_ALWAYS_INLINE__
static inline 
__m512i
cycle_bit_lshift_lt_64_avx512(__m512i data,std::int16_t cycle_lshift,
                              std::int16_t zc_size,std::int8_t zc_index,
                              __m512i swap_Idx0) 
{
    __m512i x1,x2,bit_mask;
    std::int32_t e0,e1,e2;

    cycle_lshift %= zc_size;
    if (zc_size >= 64) 
    {
        e0 = 0xFFFFFFFF;
        e1 = 0xFFFFFFFF;
        e2 = (1 << (zc_size - 64)) - 1;
    } else if (zc_size >= 32) 
    {
        e0 = 0xFFFFFFFF;
        e1 = (1 << (zc_size - 32)) - 1;
        e2 = 0;
    } else 
    {      
        e0 = (1 << zc_size) - 1;
        e1 = 0;
        e2 = 0;
    }
    bit_mask = _mm512_set_epi32(0,0,0,0,0,0,0,0,0,0,0,0,0,e2,e1,e0);
    data     = _mm512_and_si512(data,bit_mask);
    x1       = _mm512_srli_epi64(data,cycle_lshift);
    x2       = _mm512_slli_epi64(data,zc_size-cycle_lshift);
    x1       = _mm512_or_epi64(x1,x2);
    x1       = _mm512_and_si512(x1,bit_mask);
    return (x1);
} 

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
__ATTR_ALWAYS_INLINE__
static inline 
__m512i
cycle_bit_lshift_spec_values(__m512i data,std::int16_t cycle_lshift,
                             std::int16_t zc_size,std::int8_t zc_index,
                            __m512i swap_Idx0) 
{
    __m512i x1;
    int zcSizeInBytes = zc_size >> 3;
    std::uint8_t* __restrict__ pBytesData = reinterpret_cast<std::uint8_t* __restrict__>(&data);
    std::uint8_t* __restrict__ pBytesX1   = reinterpret_cast<std::uint8_t* __restrict__>(&x1);
    cycle_lshift %= zc_size;
    std::int32_t cycRightShift = (zc_size - cycle_lshift) % zc_size;
    std::int32_t byteShift1 = cycRightShift >> 3;
    std::int32_t byteShift2 = byteShift1 + 1;
    std::int32_t bitShift = cycRightShift - (byteShift1 << 3);
    for (std::int32_t byteIndex = 0; byteIndex < zcSizeInBytes; byteIndex++)
    {
        std::int32_t topByte    = (byteIndex + zcSizeInBytes - byteShift1) % zcSizeInBytes;
        std::int32_t botByte    = (byteIndex + zcSizeInBytes - byteShift2) % zcSizeInBytes;
        pBytesX1[byteIndex] = (pBytesData[topByte] << bitShift) | (pBytesData[botByte] >> (8 - bitShift));
    }
    for (std::int32_t byteIndex=zcSizeInBytes; byteIndex < 16;byteIndex++){
        pBytesX1[byteIndex] = 0;
    }
    return x1;
}

typedef __m512i (* CYCLE_BIT_LSHIFT_FPTR)(__m512i,std::int16_t,std::int16_t,std::int8_t,__m512i);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool reorder_instr_layout>
__ATTR_ALWAYS_INLINE__
static inline
CYCLE_BIT_LSHIFT_FPTR ldpc_select_lshift_func(std::int16_t zc_size)
{
    if (zc_size >= 288)
        return cycle_bit_lshift_288_to_384_avx512<reorder_instr_layout>;
    else if (zc_size < 64)
        return cycle_bit_lshift_lt_64_avx512;
    else if (zc_size == 72 || zc_size == 88 || zc_size == 104 || zc_size == 120)
        return cycle_bit_lshift_spec_values;
    else 
        return cycle_bit_lshift_144_to_256_avx512<reorder_instr_layout>;
}


} //l1_phy

} //gms

#endif /*__GMS_LDPC_ENCODER_CYCLE_BIT_AVX512_H__*/