

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

#ifndef __GMS_AVX512_MEMSET_H__
#define __GMS_AVX512_MEMSET_H__ 100120261244

namespace  file_info
{

    const unsigned int GMS_AVX512_MEMSET_MAJOR = 1U;
    const unsigned int GMS_AVX512_MEMSET_MINOR = 0U;
    const unsigned int GMS_AVX512_MEMSET_MICRO = 0U;
    const unsigned int GMS_AVX512_MEMSET_FULLVER =
      1000U*GMS_AVX512_MEMSET_MAJOR+100U*GMS_AVX512_MEMSET_MINOR+10U*GMS_AVX512_MEMSET_MICRO;
    const char * const GMS_AVX512_MEMSET_CREATE_DATE = "10-01-2026 12:44PM +00200 (SAT 10 JAN 2026 GMT+2)";
    const char * const GMS_AVX512_MEMSET_BUILD_DATE  = __DATE__ ":" __TIME__;
    const char * const GMS_AVX512_MEMSET_AUTHOR      =  "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_AVX512_MEMSET_DESCRIPT    =  "AVX512 memset kernels.";
}

#include <cstdint>
#include "GMS_config.h"

namespace gms
{
namespace common
{


__ATTR_HOT__
__ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
__ATTR_OPTIMIZE_03__
#endif 
void avx512_memset_unroll8x_ps(float * __restrict__ ,const float,std::size_t);


__ATTR_HOT__
__ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
__ATTR_OPTIMIZE_03__
#endif 
void avx512_memset_unroll16x_ps(float * __restrict__ ,const float,std::size_t);


__ATTR_HOT__
__ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
__ATTR_OPTIMIZE_03__
#endif 
void avx512_memset_unroll8x_pd(double * __restrict__ ,const double,std::size_t);


__ATTR_HOT__
__ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
__ATTR_OPTIMIZE_03__
#endif 
void avx512_memset_unroll16x_pd(double * __restrict__,const double,std::size_t);



} // common
} // gms




#endif /*__GMS_AVX512_MEMSET_H__*/