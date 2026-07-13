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

#ifndef __GMS_LDPC_ENCODER_5GNR_H__
#define __GMS_LDPC_ENCODER_5GNR_H__ 090720260615

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"

namespace file_info 
{

     static const unsigned int GMS_LDPC_ENCODER_5GNR_MAJOR = 1;
     static const unsigned int GMS_LDPC_ENCODER_5GNR_MINOR = 1;
     static const unsigned int GMS_LDPC_ENCODER_5GNR_MICRO = 0;
     static const unsigned int GMS_LDPC_ENCODER_5GNR_FULLVER =
       1000U*GMS_LDPC_ENCODER_5GNR_MAJOR+100U*GMS_LDPC_ENCODER_5GNR_MINOR+
       10U*GMS_LDPC_ENCODER_5GNR_MICRO;
     static const char GMS_LDPC_ENCODER_5GNR_CREATION_DATE[] = "09-07-2026 06:16PM +00200 (THR 09 JUL 2026 GMT+2)";
     static const char GMS_LDPC_ENCODER_5GNR_BUILD_DATE[]    = __DATE__; 
     static const char GMS_LDPC_ENCODER_5GNR_BUILD_TIME[]    = __TIME__;
     static const char GMS_LDPC_ENCODER_5GNR_SYNOPSIS[]      = "LDPC Encoder foundational data structures.";

}


namespace gms 
{

namespace l1_phy 
{

namespace detail
{

constexpr std::int32_t MAX_CB_BLOCK   = 32;
constexpr std::int32_t PROC_BYTES     = 64;
constexpr std::int32_t WAYS_144to256  = 2;
constexpr std::int32_t WAYS_72to128   = 4;
constexpr std::int32_t WAYS_36to64    = 8;
constexpr std::int32_t WAYS_18to32    = 16;
constexpr std::int32_t WAYS_2to16     = 16;

constexpr std::int32_t I_LS_NUM       = 8;
constexpr std::int32_t ZC_MAX         = 384;
constexpr std::int32_t BG1_COL_TOTAL  = 68;
constexpr std::int32_t BG1_ROW_TOTAL  = 46;
constexpr std::int32_t BG1_COL_INF_NUM= 22;
constexpr std::int32_t BG1_NONZERO_NUM= 307;
constexpr std::int32_t BG2_COL_TOTAL  = 52;
constexpr std::int32_t BG2_ROW_TOTAL  = 42;
constexpr std::int32_t BG2_COL_INF_NUM= 10;
constexpr std::int32_t BG2_NONZERO_NUM= 188;

}

/*
    Based on the original description:
    The bblib_ldpc_encoder_5gnr kernel is a 5G NR LDPC Encoder function.
    It is implemented as defined in TS38212 5.3.2

    __Requirements and Test Coverage:__

    BaseGraph 1(.2). Lifting Factor >176.

    __Algorithm Guidance:__

    The algorithm is implemented as defined in TS38212 5.3.2.
*/

/*
    removed: bblib prefix
    Structure for input parameters in API of LDPC Encoder for 5GNR.
*/
alignas(64) struct gms_ldpc_encoder_5gnr_request_t
{
    std::int8_t * __restrict__ input[detail::MAX_CB_BLOCK];
    std::int32_t  base_graph;
    std::int32_t  n_rows;
    std::uint16_t z_c;
    std::int8_t   num_code_blocks;
};

/*
   structure for outputs of LDPC encoder for 5GNR.
*/
alignas(64) struct gms_ldpc_encoder_5gnr_response_t
{
    std::int8_t * __restrict__ output[detail::MAX_CB_BLOCK];
};

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_prefetching,std::int32_t direct>
__ATTR_ALWAYS_INLINE__
static inline 
void 
adapter_from_288_to_384_avx512(std::int8_t ** __restrict__ p_buff0,
                               std::int8_t *  __restrict__ p_buff1,
                               std::uint16_t zc_size,
                               std::uint32_t cb_len)
                               
{
    __m512i x0;
    std::int8_t * __restrict__ p_buff0_offset = nullptr;
    std::int8_t * __restrict__ p_buff1_offset = nullptr;
    std::uint32_t cb_len_s3{};
    std::int16_t byte_num{};
    __mmask64 mask0{};
    bool cond_use_prefetching{};
    ///////////////////////////////////////////////////////
    cb_len_s3      = cb_len >> 3;
    p_buff0_offset = p_buff0[0];
    byte_num       = zc_size >> 3;
    p_buff1_offset = p_buff1;
    mask0          = (static_cast<__mmask64>(1) << (byte_num))-1;
    cond_use_prefetching = (static_cast<std::int32_t>(use_prefetching) == 
                           static_cast<std::int32_t>(true));
    if constexpr(1 == direct)
    {
        for(std::int16_t i = 0; i < cb_len_s3; i += byte_num)
        { 
            x0 = _mm512_maskz_load_epi8(mask0,p_buff0_offset+i);
            _mm512_store_si512(p_buff1_offset,x0);
            p_buff1_offset += detail::PROC_BYTES;
        }       
    }
    else 
    {
        if constexpr(cond_use_prefetching)
        {
            for(std::int16_t i = 0; i < cb_len_s3; i += byte_num) 
            {
                _mm_prefetch((const char*)p_buff1_offset,_MM_HINT_T0);
                x0 = _mm512_loadu_si512(p_buff1_offset);
                p_buff0_offset += detail::PROC_BYTES;
                _mm512_mask_store_epi8(p_buff0_offset,mask0,x0);
                p_buff1_offset += byte_num;
            }
        }
        else 
        {
            for(std::int16_t i = 0; i < cb_len_s3; i += byte_num) 
            {
                x0 = _mm512_loadu_si512(p_buff1_offset);
                p_buff0_offset += detail::PROC_BYTES;
                _mm512_mask_store_epi8(p_buff0_offset,mask0,x0);
                p_buff1_offset += byte_num;
            }
        }
    }
}

#define BITMASKU8(x) ((1U << (x)) - 1)
#define MIN(a,b) (((a)<(b))?(a):(b))

/*
    No changes were made to this function!!
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
__ATTR_ALWAYS_INLINE__
static inline 
void 
scatter_bytes_slow(std::uint8_t * __restrict dst,
             const std::uint8_t * __restrict src,
             const std::uint32_t num_bits,
             std::uint8_t src_offbits)
{
    while (num_bits != 0) 
    {
        unsigned num_bits_inB = MIN(8, num_bits);
        uint8_t newB;
        if (src_offbits == 0)
            newB = src[0];
        else
            newB = ((src[0] & 0xFF) >> src_offbits) | ((src[1] & 0xFF) << (8 - src_offbits));
        dst[0] = newB & BITMASKU8(num_bits_inB);
        num_bits -= num_bits_inB;
        dst++;
        src++;
    }
}

/*
    No changes were made to this function!!
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
__ATTR_ALWAYS_INLINE__
static inline 
void 
gather_bytes_slow(std::uint8_t * __restrict dst,
                  const std::uint8_t __restrict src,
                  const std::uint32_t num_bits,
                  std::uint8_t dst_offbits)
{
    bool firstByte = true;
    while (num_bits != 0) 
    {
        unsigned num_bits_inB = MIN(8, num_bits);
        uint8_t newB;
        if (dst_offbits == 0) {
            // simple copy
            newB = src[0] & BITMASKU8(num_bits_inB);
            src++;
        }
        else {
            if (firstByte) {
                newB = (dst[0] & BITMASKU8(dst_offbits)) | (src[0] & 0xFF) << dst_offbits;
                num_bits_inB = 8 - dst_offbits;
                firstByte = false;
            }
            else {
                newB = ((src[0] & 0xFF) >> (8 - dst_offbits) | (src[1] & 0xFF) << dst_offbits) & BITMASKU8(num_bits_inB);
                src++;
            }
        }
        dst[0] = newB;
        num_bits -= num_bits_inB;
        dst++;
    }
}

/*
    Minor change was made to this function.
    Removed branch by using if constexpr compile time evaluation
    Desc:
         Plain C and slow implementation of adapter function for scattering/gathering data
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<std::int32_t direct>
__ATTR_ALWAYS_INLINE__
static inline 
void 
adapter_low_speed(std::int8_t ** __restrict__ p_buff0,
                  std::int8_t *  __restrict__ p_buff1,
                  const std::int16_t zc_size,
                  std::uint32_t cb_len)
{
    std::int8_t * __restrict__ p_buff0_offset = &p_buff0[0];
    std::int8_t * __restrict__ p_buff1_offset = p_buff1;
    std::uint32_t cb_len_d_zc_size = cb_len/zc_size;
    std::uint8_t dst_offbits = 0;
    std::uint8_t src_offbits = 0;
    const bool cond_direct   = (1==direct);
    if constexpr(cond_direct)
    {
        for (int16_t i = 0; i < cb_len_d_zc_size; i++) 
        {
            scatter_bytes_slow((uint8_t* __restrict)p_buff1_offset,(uint8_t* __restrict)p_buff0_offset, zc_size, src_offbits);
            uint8_t byteOffset = (src_offbits + zc_size) >> 3;
            src_offbits = (src_offbits + zc_size) - (byteOffset << 3);
            p_buff0_offset = p_buff0_offset + byteOffset;
            p_buff1_offset = p_buff1_offset + detail::PROC_BYTES;
        }
    } 
    else 
    {
        for (int16_t i = 0; i < cb_len_d_zc_size; i++) 
        {
            gather_bytes_slow((uint8_t* __restrict)p_buff0_offset, (uint8_t* __restrict)p_buff1_offset, zc_size, dst_offbits);
            uint8_t byteOffset = (dst_offbits + zc_size) >> 3;
            dst_offbits = (dst_offbits + zc_size) - (byteOffset << 3);
            p_buff0_offset = p_buff0_offset + byteOffset;
            p_buff1_offset = p_buff1_offset + detail::PROC_BYTES;
        }
    }
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_prefetching,std::int32_t direct>
__ATTR_ALWAYS_INLINE__
static inline 
void 
adapter_2ways_144_to_256(std::int8_t ** __restrict pbuff0,
                         std::int8_t  * __restrict pbuff1,
                         std::uint16_t zc_size,
                         const std::uint32_t cb_len)
{

    if(__builtin_expect(zc_size>detail::ZC_MAX,0)) { return;}
    __ATTR_ALIGN__(64)
    int16_t adapterPermuteTableShort[64] = 
    {
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
        32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63
    };
    alignas(64) struct bind_to_cache_line
    {
        __m512i swapIdx0;
        __m512i swapIdx1;
        __m512i x0;
        __m512i x1;
        __m512i x2;
        __m512i x3;
    }hot_data;
    std::int8_t * __restrict__ p_buff0_offset[detail::WAYS_144to256] = {nullptr};
    std::int8_t * __restrict__ p_buff1_offset = nullptr;
    __mmask64 mask0;
    __mmask32 mask1;
    __mmask32 mask2;
    std::uint16_t zc_size_m2,short_num,zc_2way_byte_num;
    //////////////////////////////////////////////////////////////////////////////
    p_buff0_offset[0] = pbuff0[0];
    p_buff0_offset[1] = pbuff0[1];
    p_buff1_offset    = pbuff1;
    short_num         = zc_size>>4;
    zc_2way_byte_num  = zc_size>>2;
    mask0             = (static_cast<__mmask64>(1) << (zc_2way_byte_num))-1;
    if(64 == zc_2way_byte_num)
    {
        mask0 = 0xffffffffffffffff;
    }
    zc_size_m2 = zc_size>>1;
    mask1              = (static_cast<__mmask32>(1) << (short_num))-1;
    const bool cond_direct = (1==direct);
    const bool cond_use_prefetching = static_cast<std::int32_t>(use_prefetching)==
                                      static_cast<std::int32_t>(true);
    if constexpr(cond_direct)
    {
        mask2 = mask1<<16|mask1;
        if constexpr(cond_use_prefetching)
        {
            _mm_prefetch((const char*)&adapterPermuteTableShort[0], _MM_HINT_T0);
            _mm_prefetch((const char*)&adapterPermuteTableShort[16],_MM_HINT_T0);
        }
        hot_data.swapIdx0 = _mm512_load_si512((void const*)&adapterPermuteTableShort[0]);
        hot_data.swapIdx0 = _mm512_mask_load_epi16(hot_data.swapIdx0,mask1<<16,adapterPermuteTableShort+16);
        hot_data.swapIdx1 = _mm512_load_si512((void const*)(adapterPermuteTableShort+short_num));
        hot_data.swapIdx1 = _mm512_mask_load_epi16(hot_data.swapIdx1,mask1<<16,
                                                   (adapterPermuteTableShort+16+short_num));
        for(std::int16_t i = 0; i<cb_len; i += zc_size_m2)
        {
            hot_data.x0 = _mm512_maskz_loadu_epi8(mask0, p_buff0_offset[0]);
            hot_data.x1 = _mm512_maskz_loadu_epi8(mask0, p_buff0_offset[1]);
            p_buff0_offset[0] = p_buff0_offset[0] + zc_2way_byte_num;
            p_buff0_offset[1] = p_buff0_offset[1] + zc_2way_byte_num;
            hot_data.x2 = _mm512_mask_permutex2var_epi16(hot_data.x0, mask2, hot_data.swapIdx0, hot_data.x1);
             _mm512_storeu_si512(p_buff1_offset, hot_data.x2);
            p_buff1_offset = p_buff1_offset + detail::PROC_BYTES;
            hot_data.x3 = _mm512_mask_permutex2var_epi16(hot_data.x0, mask2, hot_data.swapIdx1, hot_data.x1);
            _mm512_storeu_si512(p_buff1_offset, hot_data.x3);
            p_buff1_offset = p_buff1_offset + detail::PROC_BYTES;
        }
    }
    else 
    {
        mask2 = (mask1 << (short_num));
        if constexpr(cond_use_prefetching)
        {
            _mm_prefetch((const char*)&adapterPermuteTableShort[0], _MM_HINT_T0);
            _mm_prefetch((const char*)&adapterPermuteTableShort[16],_MM_HINT_T0);
        }
        hot_data.swapIdx0 = _mm512_load_si512((void const*)&adapterPermuteTableShort[0]);
        hot_data.swapIdx0 = _mm512_mask_load_epi16(hot_data.swapIdx0,mask2,
                                                   ((void const*)(adapterPermuteTableShort+32-short_num)));
        hot_data.swapIdx1 = _mm512_load_si512((void const*)(adapterPermuteTableShort+16));
        hot_data.swapIdx1 = _mm512_mask_load_epi16(hot_data.swapIdx1,mask2,
                                                   ((void const*)(adapterPermuteTableShort+48-short_num)));
        for(std::int16_t i = 0; i<=cb_len; i += zc_size)
        {
            if constexpr(cond_use_prefetching)
            {
                _mm_prefetch((const char*)p_buff1_offset,_MM_HINT_T0);
            }
            hot_data.x0 = _mm512_load_si512(p_buff1_offset);
            if constexpr(cond_use_prefetching)
            {
                _mm_prefetch((const char*)(p_buff1_offset+detail::PROC_BYTES),_MM_HINT_T0);
            }
            p_buff1_offset+=detail::PROC_BYTES;
            hot_data.x1 = _mm512_load_si512(p_buff1_offset);
            if constexpr(cond_use_prefetching)
            {
                _mm_prefetch((const char*)(p_buff1_offset+detail::PROC_BYTES),_MM_HINT_T0);
            }
            p_buff1_offset+=detail::PROC_BYTES;
            hot_data.x2 = _mm512_mask_permutex2var_epi16(hot_data.x0,mask2|mask1,hot_data.swapIdx0,hot_data.x1);
            _mm512_mask_store_epi8(p_buff0_offset[0]+zc_2way_byte_num);
            hot_data.x3 = _mm512_mask_permutex2var_epi16(hot_data.x0,mask2|mask1,hot_data.swapIdx1,hot_data.x1);
            _mm512_mask_store_epi8(p_buff0_offset[1]+zc_2way_byte_num);
            p_buff0_offset[0] += zc_2way_byte_num;
            p_buff0_offset[1] += zc_2way_byte_num;
        }
    }
}

typedef void (*LDPC_ADAPTER_P)(std::int8_t ** __restrict__ , std::int8_t * __restrict__, 
                               std::uint16_t, std::uint32_t , std::int8_t);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_prefetching,std::int32_t direct>
__ATTR_ALWAYS_INLINE__
static inline 
LDPC_ADAPTER_P 
ldpc_select_adapter_func(std::uint16_t zc_size,std::uint8_t num_ways)
{
    if(zc_size<64 || zc_size==72 || zc_size==88 || zc_size==104 || zc_size==120)
       return adapter_low_speed<direct>;
    else if(zc_size>=288 || num_ways==1)
       return adapter_from_288_to_384_avx512<use_prefetching,direct>;
    else 
       return adapter_2ways_144_to_256<use_prefetching,direct>;
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_prefetching,bool reorder_instr_layout>
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
ldpc_encoder_bg1_avx512(std::int8_t * __restrict__,
                        std::int8_t * __restrict__,
                        const std::int16_t shift_matrix,
                        std::int16_t zc_size,std::uint8_t i_ls);

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_prefetching,bool reorder_instr_layout>
__ATTR_HOT__
__ATTR_ALIGN__(32)
void 
ldpc_encoder_bg2_avx512(std::int8_t * __restrict__,
                        std::int8_t * __restrict__,
                        const std::int16_t shift_matrix,
                        std::int16_t zc_size,std::uint8_t i_ls);

/*Top level interface routine*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_prefetching,std::int32_t direct,
         bool reorder_instr_layout>
__ATTR_HOT__
__ATTR_ALIGN__(32)
std::int32_t 
ldpc_encoder_5gnr_iface_avx512(gms_ldpc_encoder_5gnr_request_t * __restrict__ ,
                               gms_ldpc_encoder_5gnr_response_t* __restrict__);

__ATTR_ALIGN__(64)
extern const int16_t Bg1MatrixNumPerCol[detail::BG1_COL_TOTAL];
__ATTR_ALIGN__(64)
extern const int16_t Bg1Address[detail::BG1_NONZERO_NUM];
__ATTR_ALIGN__(64)
extern const int16_t Bg1HShiftMatrix[detail::BG1_NONZERO_NUM*detail::I_LS_NUM];
__ATTR_ALIGN__(64)
extern const int16_t Bg2MatrixNumPerCol[detail::BG2_COL_TOTAL];
__ATTR_ALIGN__(64)
extern const int16_t Bg2Address[detail::BG2_NONZERO_NUM];
__ATTR_ALIGN__(64)
extern const int16_t Bg2HShiftMatrix[detail::BG2_NONZERO_NUM*detail::I_LS_NUM];


} // l1_phy

} // gms

#endif /*__GMS_LDPC_ENCODER_5GNR_H__*/