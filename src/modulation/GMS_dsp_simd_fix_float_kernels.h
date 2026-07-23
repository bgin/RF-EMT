
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
/*
 * Copyright (C) Bernard Gingold, 2020-2026 
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#ifndef __GMS_DSP_SIMD_FIX_FLOAT_KERNELS_H__
#define __GMS_DSP_SIMD_FIX_FLOAT_KERNELS_H__ 210720260642

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"

namespace file_info 
{

     static const unsigned int GMS_DSP_SIMD_FIX_FLOAT_KERNELS_MAJOR = 1;
     static const unsigned int GMS_DSP_SIMD_FIX_FLOAT_KERNELS_MINOR = 1;
     static const unsigned int GMS_DSP_SIMD_FIX_FLOAT_KERNELS_MICRO = 0;
     static const unsigned int GMS_DSP_SIMD_FIX_FLOAT_KERNELS_FULLVER =
       1000U*GMS_DSP_SIMD_FIX_FLOAT_KERNELS_MAJOR+100U*GMS_DSP_SIMD_FIX_FLOAT_KERNELS_MINOR+
       10U*GMS_DSP_SIMD_FIX_FLOAT_KERNELS_MICRO;
     static const char GMS_DSP_SIMD_FIX_FLOAT_KERNELS_CREATION_DATE[] = "21-07-2026 06:42PM +00200 (TUE 21 JUL 2026 GMT+2)";
     static const char GMS_DSP_SIMD_FIX_FLOAT_KERNELS_BUILD_DATE[]    = __DATE__; 
     static const char GMS_DSP_SIMD_FIX_FLOAT_KERNELS_BUILD_TIME[]    = __TIME__;
     static const char GMS_DSP_SIMD_FIX_FLOAT_KERNELS_SYNOPSIS[]      = "DSP Fixed and Floating-Point helper kernels.";

}

namespace gms 
{

namespace l1_phy 
{

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
__ATTR_ALWAYS_INLINE__
static inline 
__m512
negate_zmm16r4_f32(const __m512 v)
{
    return (_mm512_or_ps(_mm512_set1_ps(-0.0f),v));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=avx
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx")
#endif
__ATTR_ALWAYS_INLINE__
static inline 
__m256
negate_ymm8r4_f32(const __m256 v)
{
    return (_mm256_or_ps(_mm256_set1_ps(-0.0f),v)); 
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
negate_zmm32i2_i16(const __m512i v)
{
    return (_mm512_sub_epi16(_mm512_setzero_si512(),_mm512_abs_epi16(v)));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=avx
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx")
#endif
__ATTR_ALWAYS_INLINE__
static inline 
__m256i 
negate_ymm16i2_i16(const __m256i x)
{
    return (_mm256_sub_epi16(_mm256_setzero_si256(),_mm256_abs_epi16(x)));
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
void 
mem_storeu_masked_zmm16r4_f32(__m512 * __restrict__ mem,
                                  const __mmask16 m, 
                                  const __m512 v)
{
    _mm512_mask_storeu_ps(reinterpret_cast<void*>(mem),m,v);
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
void 
mem_store_masked_zmm16r4_f32(__m512 * __restrict__ mem,
                                  const __mmask16 m, 
                                  const __m512 v)
{
    _mm512_mask_store_ps(reinterpret_cast<void*>(mem),m,v);
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=avx
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx")
#endif
__ATTR_ALWAYS_INLINE__
static inline 
void 
mem_store_masked_ymm8r4_f32(__m256 * __restrict__ mem,
                                const __mmask8 m, 
                                const __m256 v)
{
    const std::int32_t msh24{static_cast<std::int32_t>(m<<24)};
    const std::int32_t msh25{static_cast<std::int32_t>(m<<25)};
    const std::int32_t msh26{static_cast<std::int32_t>(m<<26)};
    const std::int32_t msh27{static_cast<std::int32_t>(m<<27)};
    const std::int32_t msh28{static_cast<std::int32_t>(m<<28)};
    const std::int32_t msh29{static_cast<std::int32_t>(m<<29)};
    const std::int32_t msh30{static_cast<std::int32_t>(m<<30)};
    const std::int32_t msh31{static_cast<std::int32_t>(m<<31)};
    const __m256i vmask{_mm256_set_epi32(msh24,msh25,msh26,msh27,
                                         msh28,msh29,msh30,msh31)};
    _mm256_maskstore_ps(reinterpret_cast<float*>(mem),vmask,v);
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
void 
mem_storeu_masked_zmm32i2_i16(__m512i * __restrict__ mem,
                                 const __mmask32 m,
                                 __m512i v)
{
    _mm512_mask_storeu_epi16((short*)mem,m,v);
}



#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=avx
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx")
#endif
__ATTR_ALWAYS_INLINE__
static inline 
void 
mem_storeu_masked_ymm16i2_i16(__m256i * __restrict__ mem,
                                  const __mmask16 m,
                                  __m256i v)
{
    __m256i dst_addr = _mm256_load_si256(reinterpret_cast<__m256i*>(mem));
    const std::int8_t mcnvr = static_cast<std::int8_t>(m);
    const std::int8_t mlsh1 = static_cast<std::int8_t>(m<<1);
    const std::int8_t mlsh2 = static_cast<std::int8_t>(m<<2);
    const std::int8_t mlsh3 = static_cast<std::int8_t>(m<<3);
    const std::int8_t mlsh4 = static_cast<std::int8_t>(m<<4);
    const std::int8_t mlsh5 = static_cast<std::int8_t>(m<<5);
    const std::int8_t mlsh6 = static_cast<std::int8_t>(m<<6);
    const std::int8_t mlsh7 = static_cast<std::int8_t>(m<<7);
    const std::int8_t mrsh1 = static_cast<std::int8_t>(m>>1);
    const std::int8_t mrsh2 = static_cast<std::int8_t>(m>>2);
    const std::int8_t mrsh3 = static_cast<std::int8_t>(m>>3);
    const std::int8_t mrsh4 = static_cast<std::int8_t>(m>>4);
    const std::int8_t mrsh5 = static_cast<std::int8_t>(m>>5);
    const std::int8_t mrsh6 = static_cast<std::int8_t>(m>>6);
    const std::int8_t mrsh7 = static_cast<std::int8_t>(m>>7);
    const std::int8_t mrsh8 = static_cast<std::int8_t>(m>>8);
    
    __m256i mask     = _mm256_set_epi8(mrsh8,mrsh8,mrsh7,mrsh7,
                                       mrsh6,mrsh6,mrsh5,mrsh5,
                                       mrsh4,mrsh4,mrsh3,mrsh3,
                                       mrsh2,mrsh2,mrsh1,mrsh1,
                                       mcnvr,mcnvr,mlsh1,mlsh1,
                                       mlsh2,mlsh2,mlsh3,mlsh3,
                                       mlsh4,mlsh4,mlsh5,mlsh5,
                                       mlsh6,mlsh6,mlsh7,mlsh7);
    __m256i tmp_v{_mm256_blendv_epi8(dst_addr,v,mask)};
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(mem),tmp_v);
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
zmm32i2_mulhrs_zmm32i2_i16(const __m512i x,const __m512i y)
{
   return (_mm512_mulhrs_epi16(x,y));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=avx
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx")
#endif
__ATTR_ALWAYS_INLINE__
static inline 
__m256i 
ymm16i2_mulhrs_ymm16i2_i16(const __m256i x,const __m256i y)
{
   return (_mm256_mulhrs_epi16(x,y));
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
__m512 
zmm16r4_fmadd_zmm16r4_f32(const __m512 x,const __m512 y, 
                              const __m512 z)
{
    return (_mm512_fmadd_ps(x,y,z));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=avx
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx")
#endif
__ATTR_ALWAYS_INLINE__
static inline 
__m256
ymm8r4_fmadd_ymm8r4_f32(const __m256 x,const __m256 y, 
                            const __m256 z)
{
    return (_mm256_fmadd_ps(x,y,z));
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
__m512 
zmm16r4_fmaddsub_zmm16r4_f32(const __m512 x,const __m512 y, 
                                 const __m512 z)
{
    return (_mm512_fmaddsub_ps(x,y,z));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=avx
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx")
#endif
__ATTR_ALWAYS_INLINE__
static inline 
__m256
ymm8r4_fmaddsub_ymm8r4_f32(const __m256 x,const __m256 y, 
                               const __m256 z)
{
    return (_mm256_fmaddsub_ps(x,y,z));
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
zmm32i2_mulhi_zmm32i2_i16(const __m512i x,const __m512i y)
{
    return (_mm512_mulhi_epi16(x,y));
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
zmm32i2_mulhi_zmm32i2_ui16(const __m512i x,const __m512i y)
{
    return (_mm512_mulhi_epu16(x,y));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=avx
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx")
#endif
__ATTR_ALWAYS_INLINE__
static inline 
__m256i
ymm16i2_mulhi_ymm16i2_i16(const __m256i x,const __m256i y)
{
    return (_mm256_mulhi_epi16(x,y));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=avx
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx")
#endif
__ATTR_ALWAYS_INLINE__
static inline 
__m256i
ymm16i2_mulhi_ymm16i2_ui16(const __m256i x,const __m256i y)
{
    return (_mm256_mulhi_epu16(x,y));
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
__m512
zmm16r4_extract_re_f32(const __m512 c)
{
    return (_mm512_moveldup_ps(c));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=avx
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx")
#endif
__ATTR_ALWAYS_INLINE__
static inline
__m256
ymm8r4_extract_re_f32(const __m256 c)
{
    return (_mm256_moveldup_ps(c));
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
zmm32i2_extract_re_i16(const __m512i c)
{
    const __m512i split_idx{_mm512_set_epi8(13,12,13,12,9,8,9,8,5,4,5,4,1,0,1,0,
                                             13,12,13,12,9,8,9,8,5,4,5,4,1,0,1,0,
                                             13,12,13,12,9,8,9,8,5,4,5,4,1,0,1,0,
                                             13,12,13,12,9,8,9,8,5,4,5,4,1,0,1,0)};
    return (_mm512_shuffle_epi8(c,split_idx));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=avx
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx")
#endif
__ATTR_ALWAYS_INLINE__
static inline
__m256i 
ymm16i2_extract_re_i16(const __m256i c)
{
    const __m256i split_idx{_mm256_set_epi8(0,1,0,1,4,5,4,5,8,9,8,9,12,13,12,13,
                                             0,1,0,1,4,5,4,5,8,9,8,9,12,13,12,13)};
    return (_mm256_shuffle_epi8(c,split_idx));
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
__m512
zmm16r4_extract_im_f32(const __m512 c)
{
    return (_mm512_movehdup_ps(c));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=avx
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx")
#endif
__ATTR_ALWAYS_INLINE__
static inline
__m256
ymm8r4_extract_im_f32(const __m256 c)
{
    return (_mm256_movehdup_ps(c));
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
zmm32i2_extract_im_i16(const __m512i c)
{
    const __m512i split_idx{_mm512_set_epi8(15,14,15,14,11,10,11,10,7,6,7,6,3,2,3,2,
                                   15,14,15,14,11,10,11,10,7,6,7,6,3,2,3,2,
                                   15,14,15,14,11,10,11,10,7,6,7,6,3,2,3,2,
                                   15,14,15,14,11,10,11,10,7,6,7,6,3,2,3,2)};
    return (_mm512_shuffle_epi8(c,split_idx));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=avx
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx")
#endif
__ATTR_ALWAYS_INLINE__
static inline
__m256i 
ymm16i2_extract_im_i16(const __m256i c)
{
    const __m256i split_idx{_mm256_set_epi8(2,3,2,3,6,7,6,7,10,11,10,11,14,15,14,15,
                                            2,3,2,3,6,7,6,7,10,11,10,11,14,15,14,15)};
    return (_mm256_shuffle_epi8(c,split_idx));
}


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=avx
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx")
#endif
__ATTR_ALWAYS_INLINE__
static inline
__m256 
ymm8r4_cmplx_sum_unsplitted_f32(const __m256 c)
{
    __m256 t0{_mm256_permute2f128_ps(c,c,0x1)};
    t0 = _mm256_add_ps(c,t0);
    __m256 t1{_mm256_permute_ps(t0,0x4E)};
    t0 = _mm256_add_ps(t1,t0);
    return (t0);
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=avx
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx")
#endif
__ATTR_ALWAYS_INLINE__
static inline
__m256i 
ymm16i2_cmplx_sum_unsplitted_i16(const __m256i c)
{
    __m256i t0{_mm256_permute2f128_si256(c,c,0x1)};
    t0 = _mm256_add_epi16(c,t0);
    __m256i t1{_mm256_shuffle_epi32(t0,0x4e)};
    t0 = _mm256_add_epi16(t0,t1);
    t1 = _mm256_shuffle_epi32(t0,0xb1);
    t0 = _mm256_add_epi16(t0,t1);
    return (t0);
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
zmm32i2_cmplx_sum_unsplitted_i16(const __m512i c)
{
    __m512i t0{_mm512_shuffle_i64x2(c,c,0x4e)};
    t0 = _mm512_add_epi16(c,t0);
    __m512i t1{_mm512_shuffle_i64x2(t0,t0,0xb1)};
    t0 = _mm512_add_epi16(t0,t1);
    t1 = _mm512_shuffle_epi32(t0,(_MM_PERM_ENUM)0x4e);
    t0 = _mm512_add_epi16(t0,t1);
    t1 = _mm512_rol_epi64(t0,32);
    t0 = _mm512_add_epi16(t0,t1);
    return (t0);
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
zmm32i2_split_re_neg_even_i16(const __m512i c)
{
    const __m512i split_idx{_mm512_set_epi8(13,12,13,12,9,8,9,8,5,4,5,4,1,0,1,0,
                                            13,12,13,12,9,8,9,8,5,4,5,4,1,0,1,0,
                                            13,12,13,12,9,8,9,8,5,4,5,4,1,0,1,0,
                                            13,12,13,12,9,8,9,8,5,4,5,4,1,0,1,0)};
    const __m512i t0{_mm512_shuffle_epi8(c,split_idx)};
    return (_mm512_mask_sub_epi16(t0,0xaaaaaaaa,_mm512_setzero_si512(),t0));
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=avx
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx")
#endif
__ATTR_ALWAYS_INLINE__
static inline
__m256i 
ymm16i2_split_re_neg_even_i16(const __m256i c)
{
    const __m256i split_idx{_mm256_set_epi8(0,1,0,1,4,5,4,5,8,9,8,9,12,13,12,13,
                                            0,1,0,1,4,5,4,5,8,9,8,9,12,13,12,13)};
    const __m256i mask_ones{_mm256_set_epi16(1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1)};
    const __m256i t0{_mm256_shuffle_epi8(c,split_idx)};
    return (_mm256_sign_epi16(t0,mask_ones));                                              
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
zmm32i2_split_im_neg_even_i16(const __m512i c)
{
    const __m512i split_idx{_mm512_set_epi8(15,14,15,14,11,10,11,10,7,6,7,6,3,2,3,2,
                                            15,14,15,14,11,10,11,10,7,6,7,6,3,2,3,2,
                                            15,14,15,14,11,10,11,10,7,6,7,6,3,2,3,2,
                                            15,14,15,14,11,10,11,10,7,6,7,6,3,2,3,2)};
    const __m512i t0{_mm512_shuffle_epi8(c,split_idx)};
    return (_mm512_mask_sub_epi16(t0,0xaaaaaaaa,_mm512_setzero_si512(),c));                                       
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=avx
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx")
#endif
__ATTR_ALWAYS_INLINE__
static inline
__m256i 
ymm16i2_split_im_neg_even_i16(const __m256i c)
{
    const __m256i split_idx{_mm256_set_epi8(2,3,2,3,6,7,6,7,10,11,10,11,14,15,14,15,
                                            2,3,2,3,6,7,6,7,10,11,10,11,14,15,14,15)};
    const __m256i mask_ones{_mm256_setr_epi16(1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1)};
    const __m256i t0{_mm256_shuffle_epi8(c,split_idx)};
    return (_mm256_sign_epi16(t0,mask_ones));                                    
}

#if (defined(__INTEL_COMPILER) || defined(__ICC)) && !defined(__GNUC__)

#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
__ATTR_ALWAYS_INLINE__
static inline
__m512i 
zmm32i2_cmul_zmm32i2_i16(const __m512i c1,const __m512i c2)
{
    __ATTR_ALIGN__(64) struct bound_to_cache_line_t
    {
       __m512i splitted_re,splitted_im;
       __m512i neg_imag,ct1;
       __m512i ct2,c3_res;
    }hot_data;
    const __m512i split_idx_re{_mm512_set_epi8(13,12,13,12,9,8,9,8,5,4,5,4,1,0,1,0,
                                             13,12,13,12,9,8,9,8,5,4,5,4,1,0,1,0,
                                             13,12,13,12,9,8,9,8,5,4,5,4,1,0,1,0,
                                             13,12,13,12,9,8,9,8,5,4,5,4,1,0,1,0)};
    hot_data.splitted_re = _mm512_shuffle_epi8(c1,split_idx_re);
    const __m512i split_idx_im{_mm512_set_epi8(15,14,15,14,11,10,11,10,7,6,7,6,3,2,3,2,
                                             15,14,15,14,11,10,11,10,7,6,7,6,3,2,3,2,
                                             15,14,15,14,11,10,11,10,7,6,7,6,3,2,3,2,
                                             15,14,15,14,11,10,11,10,7,6,7,6,3,2,3,2)};
    hot_data.splitted_im = _mm512_shuffle_epi8(c1,split_idx_im);
    const __m512i iq_interleave_idx{_mm512_set_epi8(13,12,15,14,9,8,11,10,5,4,7,6,1,0,3,2,
                                                13,12,15,14,9,8,11,10,5,4,7,6,1,0,3,2,
                                                13,12,15,14,9,8,11,10,5,4,7,6,1,0,3,2,
                                                13,12,15,14,9,8,11,10,5,4,7,6,1,0,3,2)};
    hot_data.ct1 = _mm512_shuffle_epi8(c1,iq_interleave_idx);
    const __m512i neg_imag_idx{_mm512_set_epi16(0x0001,0xffff, 0x0001, 0xffff,
                                                0x0001,0xffff, 0x0001, 0xffff,
                                                0x0001,0xffff, 0x0001, 0xffff,
                                                0x0001,0xffff, 0x0001, 0xffff,
                                                0x0001,0xffff, 0x0001, 0xffff,
                                                0x0001,0xffff, 0x0001, 0xffff,
                                                0x0001,0xffff, 0x0001, 0xffff,
                                                0x0001,0xffff, 0x0001, 0xffff)};
    hot_data.neg_imag = _mm512_mullo_epi16(hot_data.ct1,neg_imag_idx);
    hot_data.ct1 = _mm512_mulhrs_epi16(hot_data.splitted_re,c1);
    hot_data.ct2 = _mm512_mulhrs_epi16(hot_data.splitted_im,hot_data.neg_imag);
    hot_data.c3_res = _mm512_adds_epi16(hot_data.ct1,hot_data.ct2);
    return (hot_data.c3_res);
}

#elif defined(__GNUC__)

#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
template<bool use_software_prefetching>
__ATTR_ALWAYS_INLINE__
static inline
__m512i 
zmm32i2_cmul_zmm32i2_i16(const __m512i c1,const __m512i c2)
{
    __ATTR_ALIGN__(64) struct bound_to_cache_line_t
    {
       __m512i splitted_re,splitted_im;
       __m512i neg_imag,ct1;
       __m512i ct2,c3_res;
    }hot_data;
    __ATTR_ALIGN__(64)
    static constexpr const std::uint16_t buf_neg_imag_idx[32] = 
    {
        0x0001,0xffff, 0x0001, 0xffff,
        0x0001,0xffff, 0x0001, 0xffff,
        0x0001,0xffff, 0x0001, 0xffff,
        0x0001,0xffff, 0x0001, 0xffff,
        0x0001,0xffff, 0x0001, 0xffff,
        0x0001,0xffff, 0x0001, 0xffff,
        0x0001,0xffff, 0x0001, 0xffff,
        0x0001,0xffff, 0x0001, 0xffff
    };
    __ATTR_ALIGN__(64)
    static constexpr const std::uint8_t buf_split_idx_re[64] = 
    {
        13,12,13,12,9,8,9,8,5,4,5,4,1,0,1,0,
        13,12,13,12,9,8,9,8,5,4,5,4,1,0,1,0,
        13,12,13,12,9,8,9,8,5,4,5,4,1,0,1,0,
        13,12,13,12,9,8,9,8,5,4,5,4,1,0,1,0
    };
    __ATTR_ALIGN__(64)
    static constexpr const std::uint8_t buf_split_idx_im[64] = 
    {
        15,14,15,14,11,10,11,10,7,6,7,6,3,2,3,2,
        15,14,15,14,11,10,11,10,7,6,7,6,3,2,3,2,
        15,14,15,14,11,10,11,10,7,6,7,6,3,2,3,2,
        15,14,15,14,11,10,11,10,7,6,7,6,3,2,3,2
    };
    __ATTR_ALIGN__(64)
    static constexpr const std::uint8_t buf_iq_interleave_idx[64] = 
    {
        13,12,15,14,9,8,11,10,5,4,7,6,1,0,3,2,
        13,12,15,14,9,8,11,10,5,4,7,6,1,0,3,2,
        13,12,15,14,9,8,11,10,5,4,7,6,1,0,3,2,
        13,12,15,14,9,8,11,10,5,4,7,6,1,0,3,2
    };
    if constexpr(use_software_prefetching) {_mm_prefetch((const char*)&buf_split_idx_re[0],_MM_HINT_T0);}
    const __m512i split_idx_re{_mm512_load_si512((__m512i*)&buf_split_idx_re[0])};
    hot_data.splitted_re = _mm512_shuffle_epi8(c1,split_idx_re);
    if constexpr(use_software_prefetching) {_mm_prefetch((const char*)&buf_split_idx_im[0],_MM_HINT_T0);}
    const __m512i split_idx_im{_mm512_load_si512((__m512i*)&buf_split_idx_im[0])};
    hot_data.splitted_im = _mm512_shuffle_epi8(c1,split_idx_im);
    if constexpr(use_software_prefetching) {_mm_prefetch((const char*)&buf_iq_interleave_idx[0],_MM_HINT_T0);}
    const __m512i iq_interleave_idx{__mm512_load_si512((__m512i*)&buf_iq_interleave_idx[0])};
    hot_data.ct1 = _mm512_shuffle_epi8(c1,iq_interleave_idx);
    if constexpr(use_software_prefetching) {_mm_prefetch((const char*)&buf_neg_imag_idx[0],_MM_HINT_T0);}
    const __m512i neg_imag_idx{_mm512_load_si512((__m512i*)&buf_neg_imag_idx[0])};
    hot_data.neg_imag = _mm512_mullo_epi16(hot_data.ct1,neg_imag_idx);
    hot_data.ct1 = _mm512_mulhrs_epi16(hot_data.splitted_re,c1);
    hot_data.ct2 = _mm512_mulhrs_epi16(hot_data.splitted_im,hot_data.neg_imag);
    hot_data.c3_res = _mm512_adds_epi16(hot_data.ct1,hot_data.ct2);
    return (hot_data.c3_res);
}
#endif

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_software_prefetch>
__ATTR_ALWAYS_INLINE__
static inline
void 
pack_4x_zmm16r4_vectors(const float * __restrict__ data_in,
                        float * __restrict__ data_out)            
{
    const __m512 * __restrict__ p_data_in = reinterpret_cast<const __m512* __restrict__>(data_in);
    __m512* __restrict__        p_data_out= reinterpret_cast<__m512* __restrict__>(data_out);
    if constexpr(use_software_prefetch)
    {
        _mm_prefetch((const char*)&p_data_in[0],_MM_HINT_T0);
        _mm_prefetch((const char*)&p_data_in[1],_MM_HINT_T0);
        _mm_prefetch((const char*)&p_data_in[2],_MM_HINT_T0);
        _mm_prefetch((const char*)&p_data_in[3],_MM_HINT_T0);
    }
    const __m512i pack_idx_lo{_mm512_setr_epi32(0,4,8,12,1,5,9,13,16,20,24,28,17,21,25,29)};
    const __m512  lo_1234{_mm512_unpacklo_ps(p_data_in[0],p_data_in[1])};
    const __m512  hi_1234{_mm512_unpackhi_ps(p_data_in[0],p_data_in[1])};
    const __m512i pack_idx_hi{_mm512_setr_epi32(2,6,10,14,3,7,11,15,18,22,26,30,19,23,27,31)};
    const __m512  lo_5678{_mm512_unpacklo_ps(p_data_in[2],p_data_in[3])};
    const __m512  hi_5678{_mm512_unpackhi_ps(p_data_in[2],p_data_in[3])};
    p_data_out[0] = _mm512_permutex2var_ps(lo_1234,pack_idx_lo,lo_5678);
    p_data_out[1] = _mm512_permutex2var_ps(lo_1234,pack_idx_hi,lo_5678);
    p_data_out[2] = _mm512_permutex2var_ps(hi_1234,pack_idx_lo,hi_5678);
    p_data_out[3] = _mm512_permutex2var_ps(hi_5678,pack_idx_hi,hi_5678);
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_software_prefetch>
__ATTR_ALWAYS_INLINE__
static inline
void 
pack_8x_zmm16r4_vectors(const float * __restrict__ data_in,
                        float * __restrict__ data_out)
{
     const __m512 * __restrict__ p_data_in = reinterpret_cast<const __m512* __restrict__>(data_in);
    __m512* __restrict__         p_data_out= reinterpret_cast<__m512* __restrict__>(data_out);
    if constexpr(use_software_prefetch)
    {
        _mm_prefetch((const char*)&p_data_in[0],_MM_HINT_T0);
        _mm_prefetch((const char*)&p_data_in[1],_MM_HINT_T0);
    }
    const __m512 lo1{_mm512_unpacklo_ps(p_data_in[0],p_data_in[1])};
    if constexpr(use_software_prefetch)
    {
        _mm_prefetch((const char*)&p_data_in[2],_MM_HINT_T0);
        _mm_prefetch((const char*)&p_data_in[3],_MM_HINT_T0);
    }
    const __m512i split_idx_even{_mm512_setr_epi32(0,8,2,10,1,9,3,11,16,24,18,26,17,25,19,27)};
    const __m512i split_idx_odd{_mm512_setr_epi32(4,12,6,14,5,13,7,15,20,28,22,30,21,29,23,31)};
    const __m512 lo2{_mm512_unpacklo_ps(p_data_in[2],p_data_in[3])};
    const __m512 half04_0{_mm512_unpacklo_ps(lo1,lo2)};
    const __m512 half15_0{_mm512_unpackhi_ps(lo1,lo1)};
    if constexpr(use_software_prefetch)
    {
        _mm_prefetch((const char*)&p_data_in[4],_MM_HINT_T0);
        _mm_prefetch((const char*)&p_data_in[5],_MM_HINT_T0);
    }
    const __m512 lo3{_mm512_unpacklo_ps(p_data_in[4],p_data_in[5])};
    if constexpr(use_software_prefetch)
    {
        _mm_prefetch((const char*)&p_data_in[6],_MM_HINT_T0);
        _mm_prefetch((const char*)&p_data_in[7],_MM_HINT_T0);
    }
    const __m512 lo4{_mm512_unpacklo_ps(p_data_in[6],p_data_in[7])};
    const __m512 half04_1{_mm512_unpacklo_ps(lo3,lo4)};
    p_data_out[0] = _mm512_permutex2var_ps(half04_0,split_idx_even,half04_1);
    p_data_out[4] = _mm512_permutex2var_ps(half04_0,split_idx_odd,half04_1);
    const __m512 half15_1{_mm512_unpackhi_ps(lo3,lo4)};
    p_data_out[1] = _mm512_permutex2var_ps(half15_0,split_idx_even,half15_1);
    p_data_out[5] = _mm512_permutex2var_ps(half15_0,split_idx_odd,half15_1);
    const __m512 up1 = _mm512_unpackhi_ps(p_data_in[0],p_data_in[1]);
    const __m512 up2 = _mm512_unpackhi_ps(p_data_in[2],p_data_in[3]);
    const __m512 up3 = _mm512_unpackhi_ps(p_data_in[4],p_data_in[5]);
    const __m512 up4 = _mm512_unpackhi_ps(p_data_in[6],p_data_in[7]);
    const __m512 half26_0 = _mm512_unpacklo_ps(up1, up2);
    const __m512 half26_1 = _mm512_unpacklo_ps(up3, up4);
    p_data_out[2] = _mm512_permutex2var_ps(half26_0, split_idx_even, half26_1);
    p_data_out[6] = _mm512_permutex2var_ps(half26_0, split_idx_odd, half26_1);
    const __m512 half37_0 = _mm512_unpackhi_ps(up1, up2);
    const __m512 half37_1 = _mm512_unpackhi_ps(up3, up4);
    p_data_out[3] = _mm512_permutex2var_ps(half37_0, split_idx_even, half37_1);
    p_data_out[7] = _mm512_permutex2var_ps(half37_0, split_idx_odd, half37_1);                                                                              
}

#define DO2_UNPACK_LO(x,res0,res1)\
(res0) = _mm512_unpacklo_ps(p_data_in[(x)],p_data_in[(x)+2]);\
(res1) = _mm512_unpacklo_ps(p_data_in[(x)+1],p_data_in[(x)+3]);

#define DO2_UNPACK_HI(x,res0,res1)\
(res0) = _mm512_unpackhi_ps(p_data_in[(x)],p_data_in[(x)+2]);\
(res1) = _mm512_unpackhi_ps(p_data_in[(x)+1],p_data_in[(x)+3]);

#define CREATE_ROWS_TRANSPOSE_ZMM16R4_16X16(src0,src1,src2,src3,dst0,dst1,dst2,dst3)\
{\
temp0 = _mm512_mask_shuffle_f32x4((src0),maskMid,(src1),(src2),_MM_SHUFFLE(0,0,0,0));\
(dst0) = _mm512_permutex2var_ps(temp0,k_idx0,(src3));\
temp1 = _mm512_mask_shuffle_f32x4((src1),maskBack,(src0),(src2),_MM_SHUFFLE(1,1,1,1));\                                    
(dst1) = _mm512_permutex2var_ps(temp1,k_idx1,(src3));\
temp2 = _mm512_mask_shuffle_f32x4((src2),maskFront,(src1),(src0),_MM_SHUFFLE(2,2,2,2));\              
(dst2) = _mm512_permutex2var_ps(temp2,k_idx2,(src3));\
temp3 = _mm512_mask_shuffle_f32x4((src0),maskMid,(src1),(src2),_MM_SHUFFLE(3,3,3,3));\                                        
(dst3) = _mm512_permutex2var_ps(temp3,k_idx3,(src3));\
}
                                       
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_software_prefetch>
__ATTR_ALWAYS_INLINE__
static inline
void 
transpose_zmm16r4_16x16(const float * __restrict__ data_in,
                        float * __restrict__ data_out)
{
    const __m512 * __restrict__ p_data_in{reinterpret_cast<const __m512*__restrict__>(data_in)};
    __m512 * __restrict__ p_data_out{reinterpret_cast<__m512* __restrict__>(data_out)};
    __m512i k_idx0;
    __m512i k_idx1;
    __m512i k_idx2;
    __m512i k_idx3;
    __m512  temp0;
    __m512  temp1;
    __m512  temp2;
    __m512  temp3;
    constexpr __mmask16 maskMid{0b0000111111110000};
    constexpr __mmask16 maskFront{0b1111000011110000};
    constexpr __mmask16 maskBack{0b0000111100001111};
    if constexpr(use_software_prefetching)
    {   
        __ATTR_ALIGN__(64)
        static constexpr const std::int32_t buf_k_idx0[16] = 
        {
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19
        };
        __ATTR_ALIGN__(64)
        static constexpr const std::int32_t buf_k_idx1[16] = 
        {
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 20, 21, 22, 23
        };
        __ATTR_ALIGN__(64)
        static constexpr const std::int32_t buf_k_idx2[16] = 
        {
           12, 13, 14, 15, 4, 5, 6, 7, 8, 9, 10, 11, 24, 25, 26,27
                                             
        };
        __ATTR_ALIGN__(64)
        static constexpr const std::int32_t buf_k_idx3[16] = 
        {
           12, 13, 14, 15, 4, 5, 6, 7, 8, 9, 10, 11, 28, 29, 30,31
                                             
        };
        _mm_prefetch((const char*)&buf_k_idx0[0],_MM_HINT_T0);
        _mm_prefetch((const char*)&buf_k_idx1[0],_MM_HINT_T0);
        _mm_prefetch((const char*)&buf_k_idx2[0],_MM_HINT_T0);
        _mm_prefetch((const char*)&buf_k_idx3[0],_MM_HINT_T0);
        k_idx0 = _mm512_load_si512((__m512i*)&buf_k_idx0[0]);
        k_idx1 = _mm512_load_si512((__m512i*)&buf_k_idx1[0]);
        k_idx2 = _mm512_load_si512((__m512i*)&buf_k_idx2[0]);
        k_idx3 = _mm512_load_si512((__m512i*)&buf_k_idx3[0]);
    }
    else 
    {
        k_idx0 = _mm512_setr_epi32(0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19);
        k_idx1 = _mm512_setr_epi32(0,1,2,3,4,5,6,7,8,9,10,11,20,21,22,23);
        k_idx2 = _mm512_setr_epi32(12,13,14,15,4,5,6,7,8,9,10,11,24,25,26,27);       
        k_idx3 = _mm512_setr_epi32(12,13,14,15,4,5,6,7,8,9,10,11,28,29,30,31);
                                              
    }
    __m512 up_L0_1, up_L0_2;
    DO2_UNPACK_LO(0,up_L0_1,up_L0_2);
    const __m512 R0 = _mm512_unpacklo_ps(up_L0_1,up_L0_2);
    const __m512 R4 = _mm512_unpackhi_ps(up_L0_1,up_L0_2);
    __m512 up_L0_3, up_L0_4;
    DO2_UNPACK_LO(4,up_L0_3,up_L0_4);
    const __m512 R1 = _mm512_unpacklo_ps(up_L0_3,up_L0_4);
    const __m512 R5 = _mm512_unpackhi_ps(up_L0_3,up_L0_4);
    __m512 up_L0_5, up_L0_6;
    DO2_UNPACK_LO(8,up_L0_5,up_L0_6);
    const __m512 R2 = _mm512_unpacklo_ps(up_L0_5,up_L0_6);
    const __m512 R6 = _mm512_unpackhi_ps(up_L0_5,up_L0_6);
    __m512 up_L0_7, up_L0_8;
    DO2_UNPACK_LO(12,up_L0_7,up_L0_8);
    const __m512 R3 = _mm512_unpacklo_ps(up_L0_7,up_L0_8);
    const __m512 R8 = _mm512_unpackhi_ps(up_L0_7,up_L0_8);
    CREATE_ROWS_TRANSPOSE_ZMM16R4_16X16(R0,R1,R2,R3,p_data_out[0],p_data_out[4],p_data_out[8],p_data_out[12]);
    CREATE_ROWS_TRANSPOSE_ZMM16R4_16X16(R4,R5,R6,R7,p_data_out[1],p_data_out[5],p_data_out[9],p_data_out[13]);
    //////////////////////////////////////////////////////////////////////////////////////////
    __m512 up_H0_1, up_H0_2;
    DO2_UNPACK_HI(0,up_H0_1,up_H0_2);
    const __m512 R8 = _mm512_unpacklo_ps(up_H0_1,up_H0_2);
    const __m512 R12 = _mm512_unpackhi_ps(up_H0_1,up_H0_2);
    __m512 up_H0_3, up_H0_4;
    DO2_UNPACK_HI(4,up_H0_3,up_H0_4);
    const __m512 R9 = _mm512_unpacklo_ps(up_H0_3,up_H0_4);
    const __m512 R13 = _mm512_unpackhi_ps(up_H0_3,up_H0_4);
    __m512 up_H0_5, up_H0_6;
    DO2_UNPACK_HI(8,up_H0_5,up_H0_6);
    const __m512 R10 = _mm512_unpacklo_ps(up_H0_5,up_H0_6);
    const __m512 R14 = _mm512_unpackhi_ps(up_H0_5,up_H0_6);
    __m512 up_H0_7, up_H0_8;
    DO2_UNPACK_HI(12,up_H0_7,up_H0_8);
    const __m512 R11 = _mm512_unpacklo_ps(up_H0_7,up_H0_8);
    const __m512 R15 = _mm512_unpackhi_ps(up_H0_7,up_H0_8);
    CREATE_ROWS_TRANSPOSE_ZMM16R4_16X16(R8,R9,R10,R11,p_data_out[2],p_data_out[6],p_data_out[10],p_data_out[14]);
    CREATE_ROWS_TRANSPOSE_ZMM16R4_16X16(R12,R13,R14,R15,p_data_out[3],p_data_out[7],p_data_out[11],p_data_out[15]);
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_software_prefetch>
__ATTR_ALWAYS_INLINE__
static inline
void 
unpack_zmm16r4_16x4(const float * __restrict__ data_in,
                    float * __restrict__ data_out) 
{
    const __m512 * __restrict__ p_data_in{reinterpret_cast<const __m512*__restrict__>(data_in)};
    __m512 * __restrict__ p_data_out{reinterpret_cast<__m512* __restrict__>(data_out)};
    __m512i k_idx_lo;
    __m512i k_idx_hi;
    if constexpr(use_software_prefetch)
    {
        __ATTR_ALIGN__(64)
        static constexpr const std::int32_t buf_idx_lo[16] = 
        {
            0,16,4,20,1,17,5,21,2,18,6,22,3,19,7,23
        };
        __ATTR_ALIGN__(64)
        static constexpr const std::int32_t buf_idx_hi[16] = 
        {
            8,24,12,28,9,25,13,29,10,26,14,30,11,27,15,31                            
        };
        _mm_prefetch((const char*)&buf_idx_lo[0],_MM_HINT_T0);
        _mm_prefetch((const char*)&buf_idx_hi[0],_MM_HINT_T0);
        k_idx_lo = _mm512_load_si512((__m512i*)&buf_idx_lo[0]);
        k_idx_hi = _mm512_load_si512((__m512i*)&buf_idx_hi[0]);
    }
    else 
    {
        k_idx_lo = _mm512_setr_epi32(0,16,4,20,1,17,5,21,2,18,6,22,3,19,7,23);
        k_idx_hi = _mm512_setr_epi32(8,24,12,28,9,25,13,29,10,26,14,30,11,27,15,31);
    }
    const __m512 temp0 = _mm512_permutex2var_ps(p_data_in[0], k_idx_lo, p_data_in[2]);
    const __m512 temp1 = _mm512_permutex2var_ps(p_data_in[1], k_idx_lo, p_data_in[3]);
    p_data_out[0] = _mm512_unpacklo_ps(temp0, temp1);
    p_data_out[1] = _mm512_unpackhi_ps(temp0, temp1);
    const __m512 temp2 = _mm512_permutex2var_ps(p_data_in[0], k_idx_hi, p_data_in[2]);
    const __m512 temp3 = _mm512_permutex2var_ps(p_data_in[1], k_idx_hi, p_data_in[3]);
    p_data_out[2] = _mm512_unpacklo_ps(temp2, temp3);
    p_data_out[3] = _mm512_unpackhi_ps(temp2, temp3);
}

#define CREATE_ROWS_UNPACK_ZMM16R4_16X8(src0,src1,src2,src3,dst0,dst1,dst2,dst3)\
{\
temp0 = _mm512_mask_shuffle_f32x4((src0),maskMid,(src1),(src2),_MM_SHUFFLE(0,0,0,0));\
(dst0) = _mm512_permutex2var_ps(temp0, k_idx0, (src3));\
temp1 = _mm512_mask_shuffle_f32x4((src1),maskBack,(src0),(src2),_MM_SHUFFLE(1,1,1,1));\                                        
(dst1) = _mm512_permutex2var_ps(temp1, k_idx1, (src3));\
temp2 = _mm512_mask_shuffle_f32x4((src2), maskFront,(src1),(src0),_MM_SHUFFLE(2,2,2,2));\              
(dst2) = _mm512_permutex2var_ps(temp2, k_idx2,(src3));\
temp3 = _mm512_mask_shuffle_f32x4((src0), maskMid,(src1),(src2),_MM_SHUFFLE(3,3,3,3));\                                        
(dst3) = _mm512_permutex2var_ps(temp3, k_idx3,(src3));\
}

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx512f")
#endif
template<bool use_software_prefetch>
__ATTR_ALWAYS_INLINE__
static inline
void 
unpack_zmm16r4_16x8(const float * __restrict__ data_in,
                    float * __restrict__ data_out)
{
    const __m512 * __restrict__ p_data_in{reinterpret_cast<const __m512*__restrict__>(data_in)};
    __m512 * __restrict__ p_data_out{reinterpret_cast<__m512* __restrict__>(data_out)};
    __m512 temp0;
    __m512 temp1;
    __m512 temp2;
    __m512 temp3;
    __m512i k_idx0;
    __m512i k_idx1;
    __m512i k_idx2;
    __m512i k_idx3;
    constexpr __mmask16 maskMid = 0b0000111111110000;
    constexpr __mmask16 maskFront = 0b1111000011110000;
    constexpr __mmask16 maskBack = 0b0000111100001111;
    if constexpr(use_software_prefetch)
    {
        __ATTR_ALIGN__(64)
        static constexpr const std::int32_t buf_k_idx0[16] = 
        {
            0, 1, 4, 5, 8, 9, 16, 17, 2, 3, 6, 7, 10, 11, 18, 19
        };
        __ATTR_ALIGN__(64)
        static constexpr const std::int32_t buf_k_idx1[16] = 
        {
            0, 1, 4, 5, 8, 9, 20, 21, 2, 3, 6, 7, 10, 11, 22, 23
        };
        __ATTR_ALIGN__(64)
        static constexpr const std::int32_t buf_k_idx2[16] = 
        {
            12, 13, 4, 5, 8, 9, 24, 25, 14, 15, 6, 7, 10, 11, 26,27                                           
        };
        __ATTR_ALIGN__(64)
        static constexpr const std::int32_t buf_k_idx3[16] = 
        {
            12, 13, 4, 5, 8, 9, 28, 29, 14, 15, 6, 7, 10, 11, 30,31                                         
        };
        _mm_prefetch((const char*)&buf_k_idx0[0],_MM_HINT_T0);
        _mm_prefetch((const char*)&buf_k_idx1[0],_MM_HINT_T0);
        _mm_prefetch((const char*)&buf_k_idx2[0],_MM_HINT_T0);
        _mm_prefetch((const char*)&buf_k_idx3[0],_MM_HINT_T0);
        k_idx0 = _mm512_load_si512((__m512*)&buf_k_idx0[0]);
        k_idx1 = _mm512_load_si512((__m512*)&buf_k_idx1[0]);
        k_idx2 = _mm512_load_si512((__m512*)&buf_k_idx2[0]);
        k_idx3 = _mm512_load_si512((__m512*)&buf_k_idx3[0]);
    }
    else 
    {
        k_idx0 = _mm512_setr_epi32(0,1,4,5,8,9,16,17,2,3,6,7,10,11,18,19);
        k_idx1 = _mm512_setr_epi32(0,1,4,5,8,9,20,21,2,3,6,7,10,11,22,23);
        k_idx2 = _mm512_setr_epi32(12,13,4,5,8,9,24,25,14,15,6,7,10,11,26,27);
        k_idx3 = _mm512_setr_epi32(12,13,4,5,8,9,28,29,14,15,6,7,10,11,30,31);                                                                    
    }
    const __m512 r0 = _mm512_unpacklo_ps(p_data_in[0],p_data_in[1]);
    const __m512 r5 = _mm512_unpackhi_ps(p_data_in[0],p_data_in[1]);
    const __m512 r1 = _mm512_unpacklo_ps(p_data_in[2],p_data_in[3]);
    const __m512 r6 = _mm512_unpackhi_ps(p_data_in[2],p_data_in[3]);
    const __m512 r2 = _mm512_unpacklo_ps(p_data_in[4],p_data_in[5]);
    const __m512 r7 = _mm512_unpackhi_ps(p_data_in[4],p_data_in[5]);
    const __m512 r3 = _mm512_unpacklo_ps(p_data_in[6],p_data_in[7]);
    const __m512 r8 = _mm512_unpackhi_ps(p_data_in[6],p_data_in[7]);
    CREATE_ROWS_UNPACK_ZMM16R4_16X8(r0,r1,r2,r3,p_data_out[0],p_data_out[2],p_data_out[4],p_data_out[6]);
    CREATE_ROWS_UNPACK_ZMM16R4_16X8(r4,r5,r6,r7,p_data_out[1],p_data_out[3],p_data_out[5],p_data_out[7]);
}

} //l1_phy

} //gms



#endif /*__GMS_DSP_SIMD_FIX_FLOAT_KERNELS_H__*/