/**********************************************************************
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
**********************************************************************/

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

#ifndef __GMS_BIT_REVERSE_SIMD_H__
#define __GMS_BIT_REVERSE_SIMD_H__ 200720261212

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"

namespace file_info 
{

     static const unsigned int GMS_BIT_REVERSE_SIMD_MAJOR = 1;
     static const unsigned int GMS_BIT_REVERSE_SIMD_MINOR = 1;
     static const unsigned int GMS_BIT_REVERSE_SIMD_MICRO = 0;
     static const unsigned int GMS_BIT_REVERSE_SIMD_FULLVER =
       1000U*GMS_BIT_REVERSE_SIMD_MAJOR+100U*GMS_BIT_REVERSE_SIMD_MINOR+
       10U*GMS_BIT_REVERSE_SIMD_MICRO;
     static const char GMS_BIT_REVERSE_SIMD_CREATION_DATE[] = "20-07-2026 12:16PM +00200 (MON 12 JUL 2026 GMT+2)";
     static const char GMS_BIT_REVERSE_SIMD_BUILD_DATE[]    = __DATE__; 
     static const char GMS_BIT_REVERSE_SIMD_BUILD_TIME[]    = __TIME__;
     static const char GMS_BIT_REVERSE_SIMD_SYNOPSIS[]      = "Bit-reversing utilities (scalar and avx/avx512).";

}

namespace gms 
{

namespace l1_phy
{

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=sse
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("sse")
#endif
template<bool use_prefetching>
__ATTR_ALWAYS_INLINE__
static inline 
std::int32_t 
bit_reverse_scalar(std::int8_t * __restrict__ data_inout,
                   const std::int32_t n)  
{
    if(__builtin_expect(n<=0,0)) {return (-1);}
    std::int8_t * __restrict__ p_data_inout{data_inout};
    std::int16_t byte;
    std::int8_t  bit;
    std::int8_t  tmp{};
    const std::int32_t n_rshift3{n>>3};
    constexpr bool cond_use_prefetching = (static_cast<std::int32_t>==
                                          static_cast<std::int32_t>(true));
    for(byte = 0; n_rshift3; ++byte)
    {
        if constexpr(cond_use_prefetching)
        {
            _mm_prefetch((const char*)&p_data_inout[byte*63],_MM_HINT_T0);
        }
        tmp = 0;
        const std::int8_t bval{p_data_inout[byte]};
        for(bit = 0; bit < 8; ++bit)  
        {  
            if(((bval>>bit) & 1) == 1)
            {
                tmp |= (1<<(7-bit));
            } 
        }
        p_data_inout[byte] = tmp;
    }
    return (0);
}

#define REVERSE_BITS_IN_BYTE_AVX(x,res)\
const __m256i hi_consts{_mm256_setr_epi8(0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15,\
                                         0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15)};\
const __m256i lo_consts{_mm256_slli_epi16(hi_consts,4)};\
const __m256i lo_mask{_mm256_set1_epi8(15)};\
const __m256i hi{_mm256_shuffle_epi8(lo_consts,_mm256_and_si256((x),lo_mask))};\
const __m256i lo{_mm256_shuffle_epi8(hi_consts,_mm256_and_si256(_mm256_srli_epi16((x),4),lo_mask))};\
(res) = _mm256_or_si256(lo,hi);



#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=avx
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx")
#endif
template<bool use_prefetching>
__ATTR_ALWAYS_INLINE__
static inline 
std::int32_t
bit_reverse_avx(std::int8_t * __restrict__ data_inout,
                const std::int32_t n_bits)
{
    __m256i vtmp,vtmp_res;
    __m256i      * __restrict__  p_data_inout = reinterpret_cast<__m256i* __restrict__>(&data_inout);
    std::uint8_t * __restrict__  p_tail{nullptr};
    std::uint8_t * __restrict__  p_tmp{nullptr};
    std::uint32_t  bit_div32;
    std::uint32_t  bit_mod32;
    std::uint32_t  byte;
    bit_div32 = n_bits>>8;
    bit_mod32 = (n_bits-(bit_div32<<8)+7)>>3;
    constexpr bool cond_use_prefetching = (static_cast<std::int32_t>==
                                          static_cast<std::int32_t>(true));
    for(std::uint32_t i = 0; i < bit_div32; ++i) 
    {
        if constexpr(cond_use_prefetching)
        {
            _mm_prefetch((const char*)&p_data_inout[i],_MM_HINT_T0);
        }
        vtmp = p_data_inout[i];
        REVERSE_BITS_IN_BYTE_AVX(vtmp,vtmp_res);
        p_data_inout[i] = vtmp_res;
    }
    if(bit_mod32!=0)
    {
       vtmp = *p_data_inout;
       REVERSE_BITS_IN_BYTE_AVX(vtmp,vtmp_res);
       p_tail = reinterpret_cast<std::uint8_t* __restrict__>(&p_data_inout);
       p_tmp  = reinterpret_cast<std::uint8_t* __restrict__>(&vtmp_res);
       for(byte = 0; byte < bit_mod32; ++byte) 
       {
           p_tail[byte] = p_tmp[byte];
       }
    }
}


} //l1_phy

} //gms

#endif /*__GMS_BIT_REVERSE_SIMD_H__*/