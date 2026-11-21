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

#ifndef __GMS_REFRACTIVITY_ITUR_P45311_H__
#define __GMS_REFRACTIVITY_ITUR_P45311_H__ 191120261231

namespace file_info {

const unsigned int GMS_REFRACTIVITY_ITUR_P45311_MAJOR = 1;
const unsigned int GMS_REFRACTIVITY_ITUR_P45311_MINOR = 1;
const unsigned int GMS_REFRACTIVITY_ITUR_P45311_MICRO = 0;
const unsigned int GMS_REFRACTIVITY_ITUR_P45311_FULLVER =
    1000U * GMS_REFRACTIVITY_ITUR_P45311_MAJOR + 100U * GMS_REFRACTIVITY_ITUR_P45311_MINOR +
    10U * GMS_REFRACTIVITY_ITUR_P45311_MICRO;
const char *const GMS_REFRACTIVITY_ITUR_P45311_CREATION_DATE =
    "19-11-2026 12:34 +00200 (THR 19 NOV 2026 GMT+2)";
const char *const GMS_REFRACTIVITY_ITUR_P45311_BUILD_DATE = __DATE__ " " __TIME__;
const char *const GMS_REFRACTIVITY_ITUR_P45311_SYNOPSIS =
    "ITU-R P.453-11 - radio refractivity index calculation.";

} // namespace file_info

#include <cstdint>
#include <cmath>
#include "GMS_config.h"

namespace gms 
{

namespace radiolocation
{

/*
Usage: 
const float refraction = water_vapour_pressure_e<float,32ull,1>(1.0f,2.0f,3.0f);
*/

/*
   t ! temperature (C)
   P ! pressure (hPa)
   H ! relative humidity (%)
*/
#pragma omp declare simd simdlen(SimdLen) notinbranch
template<typename T, std::size_t reg_width, std::int32_t water_or_ice, std::size_t SimdLen = reg_width / sizeof(T)>
static inline 
T water_vapour_pressure_e(const T t, const T P, const T H)
{
    static_assert(water_or_ice == 1 || water_or_ice == 2, 
                  "Invalid argument passed: accepting only 1 or 2.");
    static_assert(reg_width == 32ull || reg_width == 64ull,
                  "Invalid argument passed: accepting only 32ull or 64ull.");
    static_assert(SimdLen > 0, "SimdLen must be greater than 0.");

    constexpr T a_wat = static_cast<T>(6.1121);
    constexpr T b_wat = static_cast<T>(18.678);
    constexpr T c_wat = static_cast<T>(257.14);
    constexpr T d_wat = static_cast<T>(234.5);
    //! Ice constants
    constexpr T a_ice = static_cast<T>(6.1125);
    constexpr T b_ice = static_cast<T>(23.036);
    constexpr T c_ice = static_cast<T>(279.82);
    constexpr T d_ice = static_cast<T>(333.7);
    constexpr T C1    = static_cast<T>(1.0);
    constexpr T C0001 = static_cast<T>(0.0001);
    
    T ef_wat, ef_ice, num, den;
    T tt, td, es, e;
    
    if constexpr (water_or_ice == 1)
    {
       tt     = t * t; 
       td     = t / d_wat; 
       num    = (b_wat - td) * t;
       den    = t + c_wat;
       ef_wat = C1 + C0001 * (static_cast<T>(7.2) + P * (static_cast<T>(0.00320) + static_cast<T>(0.0000059) * tt));
       es     = ef_wat * a_wat * std::exp(num / den);
       e      = H * es * static_cast<T>(0.01);
    }
    else if constexpr (water_or_ice == 2)
    {
       tt     = t * t; 
       td     = t / d_ice;  
       num    = (b_ice - td) * t;
       den    = t + c_ice; 
       ef_ice = C1 + C0001 * (static_cast<T>(2.2) + P * (static_cast<T>(0.00382) + static_cast<T>(0.0000064) * tt));
       es     = ef_ice * a_ice * std::exp(num / den); // Note: fixed a_wat -> a_ice here!
       e      = H * es * static_cast<T>(0.01);
    }
    return e;
}

/*
   tc ! temperature (C)
   Tk ! temperatur (K)
   Pd ! dry atmospheric pressure  (hPa)
   P ! pressure (hPa)
   H ! relative humidity (%)
*/
#pragma omp declare simd simdlen(SimdLen) notinbranch
template<typename T, std::size_t reg_width, std::int32_t water_or_ice, std::size_t SimdLen = reg_width / sizeof(T)>
static inline 
T refractivity_index_n(const T tc, const T Tk,
                       const T Pd, const T P,
                       const T H)
{
   float PdT,eT,eTT,e,n_tmp;
   float n;
   e       = water_vapour_pressure_e<T,reg_width,water_or_ice>(tc,P,H);
   PdT     = (static_cast<T>(76.6)*PdT)/Tk; 
   eT      = (static_cast<T>(72.0)*e)/Tk; 
   eTT     = static_cast<T>(3750000.0)*(e/(Tk*Tk)); 
   n_tmp   = PdT+eT+eTT; 
   n       = static_cast<T>(1.0)+n_tmp*static_cast<T>(0.000001);
}




}

}

#endif /*__GMS_REFRACTIVITY_ITUR_P45311_H__*/