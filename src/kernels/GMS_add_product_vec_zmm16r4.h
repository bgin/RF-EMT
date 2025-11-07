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

#ifndef __GMS_ADD_PRODUCT_VEC_ZMM16R4_H__
#define __GMS_ADD_PRODUCT_VEC_ZMM16R4_H__ 051120250713


namespace file_version 
{

    const unsigned int GMS_ADD_PRODUCT_VEC_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_ADD_PRODUCT_VEC_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_ADD_PRODUCT_VEC_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_ADD_PRODUCT_VEC_ZMM16R4_FULLVER =
      1000U*GMS_ADD_PRODUCT_VEC_ZMM16R4_MAJOR+
      100U*GMS_ADD_PRODUCT_VEC_ZMM16R4_MINOR+
      10U*GMS_ADD_PRODUCT_VEC_ZMM16R4_MICRO;
    const char * const GMS_ADD_PRODUCT_VEC_ZMM16R4_CREATION_DATE = "05-11-2025 07:03 AM +00200 (WED 05 NOV 2025 GMT+2)";
    const char * const GMS_ADD_PRODUCT_VEC_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_ADD_PRODUCT_VEC_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_ADD_PRODUCT_VEC_ZMM16R4_DESCRIPTION   = "AVX512 optimized operations of addition of vector and constant product.";

}

#include <cstdint>
#include "GMS_config.h"

#if !defined(ADD_PRODUCT_VEC_ZMM16R4_ADD_PEEL_LOOP)
#define ADD_PRODUCT_VEC_ZMM16R4_ADD_PEEL_LOOP 1
#endif 

#if !defined(ADD_PRODUCT_VEC_ZMM16R4_SOFT_PREFETCH)
#define ADD_PRODUCT_VEC_ZMM16R4_SOFT_PREFETCH 1 
#endif 


namespace gms 
{

namespace math 
{

         /*
             Accumulates product of vector and scalar
         */
          __ATTR_HOT__
          __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
          __ATTR_OPTIMIZE_03__
#endif 
          std::int32_t 
          vaddpc_f32_1x512_1xf32_u16x_u(const float * __restrict__,
                                        float       * __restrict__,
                                        const float               ,
                                        std::size_t);


          __ATTR_HOT__
          __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
          __ATTR_OPTIMIZE_03__
#endif 
          std::int32_t 
          vaddpc_f32_1x512_1xf32_u16x_a(const float * __restrict__,
                                        float       * __restrict__,
                                        const float               ,
                                        std::size_t);   


          __ATTR_HOT__
          __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
          __ATTR_OPTIMIZE_03__
#endif 
          std::int32_t 
          vaddpc_f32_1x512_1xf32_u10x_u(const float * __restrict__,
                                        float       * __restrict__,
                                        const float               ,
                                        std::size_t);


          __ATTR_HOT__
          __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
          __ATTR_OPTIMIZE_03__
#endif 
          std::int32_t 
          vaddpc_f32_1x512_1xf32_u10x_a(const float * __restrict__,
                                        float       * __restrict__,
                                        const float               ,
                                        std::size_t);  


           __ATTR_HOT__
          __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
          __ATTR_OPTIMIZE_03__
#endif 
          std::int32_t 
          vaddpc_f32_1x512_1xf32_u6x_u(const float * __restrict__,
                                        float       * __restrict__,
                                        const float               ,
                                        std::size_t);


          __ATTR_HOT__
          __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
          __ATTR_OPTIMIZE_03__
#endif 
          std::int32_t 
          vaddpc_f32_1x512_1xf32_u6x_a(const float * __restrict__,
                                        float       * __restrict__,
                                        const float               ,
                                        std::size_t);                  

           /*
             Accumulates product of vector and vector
          */

          __ATTR_HOT__
          __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
          __ATTR_OPTIMIZE_03__
#endif 
          std::int32_t 
          vaddpv_f32_2x512_u16x_u(const float * __restrict__,
                                  const float * __restrict__,
                                  float       * __restrict__,
                                  std::size_t);


          __ATTR_HOT__
          __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
          __ATTR_OPTIMIZE_03__
#endif 
          std::int32_t 
          vaddpv_f32_2x512_u16x_a(const float * __restrict__,
                                  const float * __restrict__,
                                  float       * __restrict__,
                                  std::size_t);   


          __ATTR_HOT__
          __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
          __ATTR_OPTIMIZE_03__
#endif 
          std::int32_t 
          vaddpv_f32_2x512_u10x_u(const float * __restrict__,
                                  const float * __restrict__,
                                  float       * __restrict__,
                                  std::size_t);


          __ATTR_HOT__
          __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
          __ATTR_OPTIMIZE_03__
#endif 
          std::int32_t 
          vaddpv_f32_2x512_u10x_a(const float * __restrict__,
                                  const float * __restrict__,
                                  float       * __restrict__,
                                  std::size_t);


          __ATTR_HOT__
          __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
          __ATTR_OPTIMIZE_03__
#endif 
          std::int32_t 
          vaddpv_f32_2x512_u6x_u(const float * __restrict__,
                                  const float * __restrict__,
                                  float       * __restrict__,
                                  std::size_t);


          __ATTR_HOT__
          __ATTR_ALIGN__(32)
#if defined(__INTEL_COMPILER) || defined(__ICC)
          __ATTR_OPTIMIZE_03__
#endif 
          std::int32_t 
          vaddpv_f32_2x512_u6x_a(const float * __restrict__,
                                  const float * __restrict__,
                                  float       * __restrict__,
                                  std::size_t);

}


}













#endif /*__GMS_ADD_PRODUCT_VEC_ZMM16R4_H__*/