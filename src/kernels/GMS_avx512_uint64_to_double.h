


#ifndef __GMS_AVX512_UINT64_TO_DOUBLE_H__
#define __GMS_AVX512_UINT64_TO_DOUBLE_H__ 231020260115

#include <immintrin.h>
#include "GMS_config.h

#if !defined(AVX512_UINT64_TO_DOUBLE_OVERRIDE_COMPILER_CMD_LINE)
#define AVX512_UINT64_TO_DOUBLE_OVERRIDE_COMPILER_CMD_LINE 0
#endif

namespace gms
{

namespace common
{

#if (AVX512_UINT64_TO_DOUBLE_OVERRIDE_COMPILER_CMD_LINE) == 1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=skylake-avx512
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER) || !defined(__ICC))
#pragma GCC optimize("O3")
#pragma GCC target("avx2")
#endif
#endif 
void avx512_uint64_to_double_u8x(const std::uint64_t* __restrict src, double* __restrict dst, std::size_t size) {
    std::size_t i = 0;

   
    std::size_t dst_addr = reinterpret_cast<std::size_t>(dst);
    if (dst_addr % 64 != 0 && size > 0) {
        std::size_t alignment_offset = (64 - (dst_addr % 64)) / sizeof(double);
        std::size_t peel_limit = (alignment_offset < size) ? alignment_offset : size;
        
        for (; i < peel_limit; ++i) {
            dst[i] = static_cast<double>(src[i]);
        }
    }

 
    constexpr std::size_t UNROLL_FACTOR = 8;
    constexpr std::size_t ELEMENTS_PER_REG = 8;
    constexpr std::size_t LOOP_STRIDE = UNROLL_FACTOR * ELEMENTS_PER_REG; // 64 elements

    if (size >= i + LOOP_STRIDE) {
        std::size_t main_limit = size - LOOP_STRIDE;
        for (; i <= main_limit; i += LOOP_STRIDE) {

            __m512i src0 = _mm512_loadu_si512(reinterpret_cast<const __m512i*>(src + i + 0));
            __m512i src1 = _mm512_loadu_si512(reinterpret_cast<const __m512i*>(src + i + 8));
            __m512i src2 = _mm512_loadu_si512(reinterpret_cast<const __m512i*>(src + i + 16));
            __m512i src3 = _mm512_loadu_si512(reinterpret_cast<const __m512i*>(src + i + 24));
            __m512i src4 = _mm512_loadu_si512(reinterpret_cast<const __m512i*>(src + i + 32));
            __m512i src5 = _mm512_loadu_si512(reinterpret_cast<const __m512i*>(src + i + 40));
            __m512i src6 = _mm512_loadu_si512(reinterpret_cast<const __m512i*>(src + i + 48));
            __m512i src7 = _mm512_loadu_si512(reinterpret_cast<const __m512i*>(src + i + 56));

           
            __m512d dst0 = _mm512_cvtepu64_pd(src0);
            __m512d dst1 = _mm512_cvtepu64_pd(src1);
            __m512d dst2 = _mm512_cvtepu64_pd(src2);
            __m512d dst3 = _mm512_cvtepu64_pd(src3);
            __m512d dst4 = _mm512_cvtepu64_pd(src4);
            __m512d dst5 = _mm512_cvtepu64_pd(src5);
            __m512d dst6 = _mm512_cvtepu64_pd(src6);
            __m512d dst7 = _mm512_cvtepu64_pd(src7);

           
            _mm512_store_pd(dst + i + 0, dst0);
            _mm512_store_pd(dst + i + 8, dst1);
            _mm512_store_pd(dst + i + 16, dst2);
            _mm512_store_pd(dst + i + 24, dst3);
            _mm512_store_pd(dst + i + 32, dst4);
            _mm512_store_pd(dst + i + 40, dst5);
            _mm512_store_pd(dst + i + 48, dst6);
            _mm512_store_pd(dst + i + 56, dst7);
        }
    }

    
    for (; i + ELEMENTS_PER_REG <= size; i += ELEMENTS_PER_REG) {
        __m512i src_rem = _mm512_loadu_si512(reinterpret_cast<const __m512i*>(src + i));
        __m512d dst_rem = _mm512_cvtepu64_pd(src_rem);
        _mm512_storeu_pd(dst + i, dst_rem); 
    }

   
    if (i < size) {
        std::size_t remaining = size - i;
       
        __mmask8 mask = (__mmask8)((1U << remaining) - 1U);

        // Masked load, conversion, and store to completely avoid faults and scalar fallback branches
        __m512i src_tail = _mm512_maskz_loadu_epi64(mask, src + i);
        __m512d dst_tail = _mm512_cvtepu64_pd(src_tail);
        _mm512_mask_storeu_pd(dst + i, mask, dst_tail);
    }
}

}

}

#endif /*__GMS_AVX512_UINT64_TO_DOUBLE_H__*/

