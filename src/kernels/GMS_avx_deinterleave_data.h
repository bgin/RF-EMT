/*MIT License
Copyright (c) 2020 Bernard Gingold
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef __GMS_AVX_DEINTERLEAVE_DATA_H__
#define __GMS_AVX_DEINTERLEAVE_DATA_H__ 


namespace file_info 
{

const unsigned int GMS_AVX_DEINTERLEAVE_DATA_MAJOR = 1U;
const unsigned int GMS_AVX_DEINTERLEAVE_DATA_MINOR = 0U;
const unsigned int GMS_AVX_DEINTERLEAVE_DATA_MICRO = 0U;
const unsigned int GMS_AVX_DEINTERLEAVE_DATA_FULLVER =
       1000U*GMS_AVX_DEINTERLEAVE_DATA_MAJOR+
       100U*GMS_AVX_DEINTERLEAVE_DATA_MINOR +
       10U*GMS_AVX_DEINTERLEAVE_DATA_MICRO;
const char * const GMS_AVX_DEINTERLEAVE_DATA_CREATION_DATE = "11-11-2025 01:46 PM +00200 (TUE 11 NOV 2025 GMT+2)";
const char * const GMS_AVX_DEINTERLEAVE_DATA_BUILD_DATE    = __DATE__ ":" __TIME__;
const char * const GMS_AVX_DEINTERLEAVE_DATA_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
const char * const GMS_AVX_DEINTERLEAVE_DATA_DESCRIPTION   = "SIMD AVX data deinterleaving (splitting) intrinsics.";

}

#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"


namespace gms 
{


namespace math 
{

/* splits 16 floats of this layout x0,y0,x1,y1,.....x7,y8
   into two __m256 vectors of layout x0,x1,....x7 and y0,y1,....y7
 */
__ATTR_ALWAYS_INLINE__
static inline
void _mm256_deinterleave_1x16_2x8_ps(__m256 &out_x, __m256 &out_y,
                                 const float * __restrict__ mem_in) 
{
     __m256 v_0_2;
     __m256 v_1_3;
     v_0_2 = _mm256_castps128_ps256(_mm_loadu_ps(&mem_in[0]));
     v_0_2 = _mm256_insertf128_ps(v_0_2,_mm_loadu_ps(&mem_in[8]),1);
     v_1_3 = _mm256_castps128_ps256(_mm_loadu_ps(&mem_in[4]));
     v_1_3 = _mm256_insertf128_ps(v_1_3,_mm_loadu_ps(&mem_in[12]),1);
     out_x = _mm256_shuffle_ps(v_0_2,v_1_3,_MM_SHUFFLE(2,0,2,0));
     out_y = _mm256_shuffle_ps(v_0_2,v_1_3,_MM_SHUFFLE(3,1,3,1));
}

/* splits 16 integers of this layout x0,y0,x1,y1,.....x7,y8
   into two __m256 vectors of layout x0,x1,....x7 and y0,y1,....y7
 */
__ATTR_ALWAYS_INLINE__
static inline
void _mm256_deinterleave_1x16_2x8_epi(__m256i &out_x, __m256i &out_y,
                                      const std::int32_t * __restrict__ mem_in) 
{
     __m256 v_0_2;
     __m256 v_1_3;
     v_0_2 = _mm256_castps128_ps256(_mm_loadu_ps((float*)&mem_in[0]));
     v_0_2 = _mm256_insertf128_ps(v_0_2,_mm_loadu_ps((float*)&mem_in[8]),1);
     v_1_3 = _mm256_castps128_ps256(_mm_loadu_ps((float*)&mem_in[4]));
     v_1_3 = _mm256_insertf128_ps(v_1_3,_mm_loadu_ps((float*)&mem_in[12]),1);
     out_x = _mm256_castps_si256(_mm256_shuffle_ps(v_0_2,v_1_3,_MM_SHUFFLE(2,0,2,0)));
     out_y = _mm256_castps_si256(_mm256_shuffle_ps(v_0_2,v_1_3,_MM_SHUFFLE(3,1,3,1)));
}

}// math


}// gms



























#endif /*__GMS_AVX_DEINTERLEAVE_DATA_H__*/