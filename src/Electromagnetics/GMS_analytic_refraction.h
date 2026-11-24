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

#ifndef __GMS_ANALYTIC_REFRACTION_H__
#define __GMS_ANALYTIC_REFRACTION_H__ 241120260626

namespace file_info {

const unsigned int GMS_ANALYTIC_REFRACTION_MAJOR = 1;
const unsigned int GMS_ANALYTIC_REFRACTION_MINOR = 1;
const unsigned int GMS_ANALYTIC_REFRACTION_MICRO = 0;
const unsigned int GMS_ANALYTIC_REFRACTION_FULLVER =
    1000U * GMS_ANALYTIC_REFRACTION_MAJOR + 100U * GMS_ANALYTIC_REFRACTION_MINOR +
    10U * GMS_ANALYTIC_REFRACTION_MICRO;
const char *const GMS_ANALYTIC_REFRACTION_CREATION_DATE =
    "24-11-2026 06:27 +00200 (TUE 24 NOV 2026 GMT+2)";
const char *const GMS_ANALYTIC_REFRACTION_BUILD_DATE = __DATE__ " " __TIME__;
const char *const GMS_ANALYTIC_REFRACTION_AUTHOR   = "Programmer:Bernard Gingold, contact: beniekg@gmail.com";
const char *const GMS_ANALYTIC_REFRACTION_SYNOPSIS =
    "Various characteristics and formulae of atmospheric refraction (radio waves and visible light/IR wavelengths)  \
     Based mainly on Колосов М.А. Шабельников А.В. -- Рефракция электромагнитных волн в атмосферах Земли, Венеры и Марса-Советское Радио (1976)";

} // namespace file_info

#include <cstdint>
#include <cmath>
#include "GMS_config.h"

// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (ANALYTIC_REFRACTION_USE_PMC_INSTRUMENTATION)
#define ANALYTIC_REFRACTION_USE_PMC_INSTRUMENTATION 0
#endif 

#if (ANALYTIC_REFRACTION_USE_PMC_INSTRUMENTATION) == 1
#include "GMS_hw_perf_macros.h"

#define PMC_VARS\
    uint64_t prog_counters_start[4] = {};\
    uint64_t prog_counters_stop[4]  = {};\
    uint64_t tsc_start,tsc_stop;\
    uint64_t act_cyc_start,act_cyc_stop;\
    uint64_t ref_cyc_start,ref_cyc_stop;\
    [[maybe_unused]] uint64_t dummy1;\
    [[maybe_unused]] uint64_t dummy2;\
    [[maybe_unused]] uint64_t dummy3;\
    std::int32_t core_counter_width;\
    double utilization,nom_ghz,avg_ghz;
#endif 

namespace gms
{

namespace radiolocation
{

// Reduces any angle in radians to the range [-pi/4, pi/4] for fast convergence
#pragma omp declare simd simdlen(SimdLen)
template<typename T, std::size_t reg_width, std::size_t SimdLen = reg_width / sizeof(T)>
static inline 
T approximated_tan(T x, std::int32_t terms = 10) 
{
    constexpr T PI          = static_cast<T>(3.1415926535897932384626433832795);
    constexpr T HALF_PI     = static_cast<T>(1.5707963267948966192313216916398);
    constexpr T QUARTER_PI  = static_cast<T>(0.7853981633974483096156608458199);
    bool invert{};
    // 1. Reduce x to [-pi/2, pi/2] using periodicity: tan(x + pi) = tan(x)
    x = std::fmod(x, PI);
    if (x > HALF_PI) x -= PI;
    if (x > QUARTER_PI) {
        x = HALF_PI - x;
        invert = true;
    } else if (x < -QUARTER_PI) {
        x = -HALF_PI - x;
        invert = true;
    }

    T sin_val  = x;
    T cos_val  = static_cast<T>(1.0);
    T sin_term = x;
    T cos_term = static_cast<T>(1.0);
    T x_sq     = x * x;

    for (std::int32_t n = 1; n < terms; ++n) {
        const T two_n = static_cast<T>(2*n);
        sin_term  *= -x_sq / (two_n * (two_n + static_cast<T>(1.0)));
        sin_val   += sin_term;
        cos_term *= -x_sq / ((two_n - 1) * two_n);
        cos_val += cos_term;
    }

    T result = sin_val / cos_val;
    return invert ? (static_cast<T>(1.0) / result) : result;
}

/*
    Formula 2.43, page 46
    n    refractive index at dest (observation point)
    n0   refractive index at source
    z    'z' angle at dest
    z0   'z' angle at source
    r    upper limit of integration (radius)
    R0   lower limit of integration (radius)
    phi  'phi' angle at dest
*/
#pragma omp declare simd simdlen(SimdLen) notinbranch
template<typename T, std::int32_t which_tan,std::size_t reg_width, std::size_t SimdLen = reg_width / sizeof(T)>
static inline 
T n_refract_tht_f243(const T n,const T n0,
                     T z,T z0,
                     const T r,const T R0,
                     T phi,T phi0)
{
    constexpr T one = static_cast<T>(1.0);
    T tgz, tgz0, tgphi, tgphi0; 
    T num_d, num_s, den_d, den_s; 
    T rat_s, rat_d;
    T stgz, stgphi, stgz0, stgphi0;
    T n_over_tht;
    if constexpr(which_tan == 1)
    {
       tgz    = std::tan(z);
       stgz   = tgz*tgz;
       tgz0   = std::tan(z0);
       stgz0  = tgz0*tgz0;
       num_d  = n*r*tgz; 
       tgphi  = std::tan(phi);
       stgphi = tgphi*tgphi; 
       tgphi0 = std::tan(phi0);
       stgphi0= tgphi0*tgphi0;
       num_s  = n0*R0*tgz0; 
       den_d  = std::sqrt(one+stgz+stgphi); 
       den_s  = std::sqrt(one+stgz0+stgphi0);
       rat_s  = num_s/den_s; 
       rat_d  = num_d/den_d; 
       n_over_tht = rat_d-rat_s; 
       return (n_over_tht);
    }
    else if constexpr(which_tan == 2)
    {
       tgz    = approximated_tan<T,reg_width>(z);
       stgz   = tgz*tgz;
       tgz0   = approximated_tan<T,reg_width>(z0);
       stgz0  = tgz0*tgz0;
       num_d  = n*r*tgz; 
       tgphi  = approximated_tan<T,reg_width>(phi);
       stgphi = tgphi*tgphi; 
       tgphi0 = approximated_tan<T,reg_width>(phi0);
       stgphi0= tgphi0*tgphi0;
       num_s  = n0*R0*tgz0; 
       den_d  = std::sqrt(one+stgz+stgphi); 
       den_s  = std::sqrt(one+stgz0+stgphi0);
       rat_s  = num_s/den_s; 
       rat_d  = num_d/den_d; 
       n_over_tht = rat_d-rat_s; 
       return (n_over_tht);
    }
}

/*
    Formula 2.43, page 46
    n    refractive index at dest (observation point)
    n0   refractive index at source
    z    'z' angle at dest
    z0   'z' angle at source
    r    upper limit of integration (radius)
    R0   lower limit of integration (radius)
    phi  'phi' angle at dest
*/
#pragma omp declare simd simdlen(SimdLen) notinbranch
template<typename T, std::int32_t which_tan,std::size_t reg_width, std::size_t SimdLen = reg_width / sizeof(T)>
static inline 
T n_refract_phi_f243(const T n,const T n0,
                     T z,T z0,
                     const T r,const T R0,
                     T phi,T phi0)
{
    constexpr T one = static_cast<T>(1.0);
    T n_over_phi; 
    T tgz, tgz0, tgphi, tgphi0; 
    T num_d, num_s, den_d, den_s; 
    T rat_s, rat_d; 
    T stgz, stgphi, stgz0, stgphi0;
    if constexpr(which_tan == 1)
    {
        tgz        = std::tan(z);
        stgz       = tgz*tgz;
        tgz0       = std::tan(z0);
        stgz0      = tgz0*tgz0;
        tgphi      = std::tan(phi);
        stgphi     = tgphi*tgphi; 
        tgphi0     = std::tan(phi0);
        stgphi0    = tgphi0*tgphi0;
        num_d      = n*r*tgphi; 
        num_s      = n0*R0*tgphi0; 
        den_d      = std::sqrt(one+stgz+stgphi); 
        den_s      = std::sqrt(one+stgz0+stgphi0);
        rat_s      = num_s/den_s; 
        rat_d      = num_d/den_d; 
        n_over_phi = rat_d-rat_s;
        return (n_over_phi);
    }
    else if constexpr(which_tan == 2)
    {
        tgz        = approximated_tan<T,reg_width>(z);
        stgz       = tgz*tgz;
        tgz0       = approximated_tan<T,reg_width>(z0);
        stgz0      = tgz0*tgz0;
        tgphi      = approximated_tan<T,reg_width>(phi);
        stgphi     = tgphi*tgphi; 
        tgphi0     = approximated_tan<T,reg_width>(phi0);
        stgphi0    = tgphi0*tgphi0;
        num_d      = n*r*tgphi; 
        num_s      = n0*R0*tgphi0; 
        den_d      = std::sqrt(one+stgz+stgphi); 
        den_s      = std::sqrt(one+stgz0+stgphi0);
        rat_s      = num_s/den_s; 
        rat_d      = num_d/den_d; 
        n_over_phi = rat_d-rat_s;
        return (n_over_phi);
    }
}

/*
    !Радиус кривизны траектории луча, formula 2.51, page: 47
    n refractive index
    z angle
    dndr derivative of refractive index at r
*/
#pragma omp declare simd simdlen(SimdLen) notinbranch
template<typename T ,std::size_t reg_width, std::size_t SimdLen = reg_width / sizeof(T)>
static inline 
T rad_ray_curvature_f251(const T n,const T z,
                         const T dndr) 
{
    T rho,t0,sinz;
    sinz = std::sin(z);
    t0   = -n/sinz; 
    rho  = t0*dndr;
    return (rho);
}

/*
    formula: 2.54, page: 48
    n refractive index
    z angle
    dndr derivative of refractive index at r
*/
#pragma omp declare simd simdlen(SimdLen) notinbranch
template<typename T ,std::size_t reg_width, std::size_t SimdLen = reg_width / sizeof(T)>
static inline 
T key_relative_f254(const T n, const T z,
                    const T dndr)
{
    constexpr T inv_erad = static_cast<T>(0.00015678896205707118218877391);
    T inv_rho,k_rel;
    inv_rho = static_cast<T>(1.0)/rad_ray_curvature_f251<T,reg_width>(n,z,dndr);
    k_rel   = inv_erad*inv_rho;
    return (k_rel);
}

}

}

#endif /*__GMS_ANALYTIC_REFRACTION_H__*/