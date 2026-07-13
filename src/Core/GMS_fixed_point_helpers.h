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

#ifndef __GMS_FIXED_POINT_HELPERS_H__
#define __GMS_FIXED_POINT_HELPERS_H__ 130720260749

#include <cstdint>
#include <cmath>
#include "GMS_config.h"

namespace file_info 
{

     static const unsigned int GMS_FIXED_POINT_HELPERS_MAJOR = 1;
     static const unsigned int GMS_FIXED_POINT_HELPERS_MINOR = 1;
     static const unsigned int GMS_FIXED_POINT_HELPERS_MICRO = 0;
     static const unsigned int GMS_FIXED_POINT_HELPERS_FULLVER =
       1000U*GMS_FIXED_POINT_HELPERS_MAJOR+100U*GMS_FIXED_POINT_HELPERS_MINOR+
       10U*GMS_FIXED_POINT_HELPERS_MICRO;
     static const char GMS_FIXED_POINT_HELPERS_CREATION_DATE[] = "13-07-2026 07:48PM +00200 (MON 13 JUL 2026 GMT+2)";
     static const char GMS_FIXED_POINT_HELPERS_BUILD_DATE[]    = __DATE__; 
     static const char GMS_FIXED_POINT_HELPERS_BUILD_TIME[]    = __TIME__;
     static const char GMS_FIXED_POINT_HELPERS_SYNOPSIS[]      = "Simple fixed-point convertion routines.";

}

namespace gms 
{

namespace common 
{

namespace 
{

constexpr std::int32_t FIXED_POINT_11_5_FRACTIONAL_BITS = 5;
constexpr std::int32_t FIXED_POINT_0_8_FRACTIONAL_BITS  = 8;

}

__ATTR_ALWAYS_INLINE__
inline static 
double 
fixed_11_5_to_double(const std::uint16_t value)
{
     return (static_cast<double>(value)/static_cast<double>(1<<FIXED_POINT_11_5_FRACTIONAL_BITS));
}

__ATTR_ALWAYS_INLINE__
inline static 
std::uint16_t 
double_to_fixed_11_5(const double value)
{
    return (static_cast<std::uint16_t>(std::round(value*(1<<FIXED_POINT_11_5_FRACTIONAL_BITS))));
}

__ATTR_ALWAYS_INLINE__
inline static 
float
fixed_11_5_to_float(const std::uint16_t value)
{
     return (static_cast<float>(value)/static_cast<float>(1<<FIXED_POINT_11_5_FRACTIONAL_BITS));
}

__ATTR_ALWAYS_INLINE__
inline static 
std::uint16_t 
float_to_fixed_11_5(const float value)
{
    return (static_cast<std::uint16_t>(std::round(value*(1<<FIXED_POINT_11_5_FRACTIONAL_BITS))));
}

__ATTR_ALWAYS_INLINE__
inline static 
float
fixed_0_8_to_float(const std::uint8_t value)
{
     return (static_cast<float>(value)/static_cast<float>(1<<FIXED_POINT_0_8_FRACTIONAL_BITS));
}

__ATTR_ALWAYS_INLINE__
inline static 
std::uint8_t 
float_to_fixed_0_8(const float value)
{
    return (static_cast<std::uint8_t>(std::round(value*(1<<FIXED_POINT_0_8_FRACTIONAL_BITS))));
}

__ATTR_ALWAYS_INLINE__
inline static 
__m128i 
_mm_vfp32_to_vfixed16(const __m128 vf32, std::int32_t rounding_mode)
{
     return (_mm_cvtps_ph(vf32,rounding_mode));
}

__ATTR_ALWAYS_INLINE__
inline static
__m128 
_mm_vfixed16_to_vfp32(const __m128i vfix16, std::int32_t rounding_mode)
{
     return (_mm_cvtph_ps(vfix16,rounding_mode));
}

__ATTR_ALWAYS_INLINE__
inline static 
__m128i 
_mm256_vfp32_to_vfixed16(const __m256 vf32, std::int32_t rounding_mode)
{
     return (_mm256_cvtps_ph(vf32,rounding_mode));
}

__ATTR_ALWAYS_INLINE__
inline static
__m256
_mm256vfixed16_to_vfp32(const __m128i vfix16, std::int32_t rounding_mode)
{
     return (_mm256_cvtph_ps(vfix16,rounding_mode));
}

__ATTR_ALWAYS_INLINE__
inline static 
__m256i 
_mm512_vfp32_to_vfixed16(const __m512 vf32, std::int32_t rounding_mode)
{
     return (_mm512_cvtps_ph(vf32,rounding_mode));
}

__ATTR_ALWAYS_INLINE__
inline static
__m512
_mm512_vfixed16_to_vfp32(const __m256i vfix16, std::int32_t rounding_mode)
{
     return (_mm512_cvtph_ps(vfix16,rounding_mode));
}


} // common

} //gms

#endif /*__GMS_FIXED_POINT_HELPERS_H__*/