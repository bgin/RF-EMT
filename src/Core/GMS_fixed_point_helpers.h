/* Copyright (C) Bernard Gingold, 2020-2026 
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
float_to_fixed_0_8_sse(const __m128 v)
{
    const __m128 vc256{_mm_set1_ps(256.0f)};
    const __m128 tmp{_mm_mul_ps(v,vc256)};
    return (_mm_cvtepi32_epi8(_mm_cvtps_epi32(tmp)));
}

__ATTR_ALWAYS_INLINE__
inline static 
__m128i
float_to_fixed_0_8_avx(const __m256 v)
{
    const __m256 vc256{_mm256_set1_ps(256.0f)};
    const __m256 tmp{_mm256_mul_ps(v,vc256)};
    return (_mm256_cvtepi32_epi8(_mm256_cvtps_epi32(tmp)));
}

__ATTR_ALWAYS_INLINE__
inline static 
__m128i
float_to_fixed_0_8_avx512(const __m512 v)
{
    const __m512 vc256{_mm512_set1_ps(256.0f)};
    const __m512 tmp{_mm512_mul_ps(v,vc256)};
    return (_mm512_cvtepi32_epi8(_mm512_cvtps_epi32(tmp)));
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
