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
    _mm512_mask_storeu_epi16((reinterpret_cast<void* __restrict__>(mem),m,v));
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
    __m256i tmp_v{_mm256_blendv_epi8(dst_mem,v,mask)};
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
__m512i 
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
    const __m512 split_idx{_mm512_set_epi8(15,14,15,14,11,10,11,10,7,6,7,6,3,2,3,2,
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
    const __m256 t0{_mm256_shuffle_epi8(c,split_idx)};
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
    return (_mm512_mask_sub_epi16(t0,0xaaaaaaaa,_mm512_setzero_si512()));                                       
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

} //l1_phy

} //gms



#endif /*__GMS_DSP_SIMD_FIX_FLOAT_KERNELS_H__*/