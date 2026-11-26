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
#include <cstring>
#include <omp.h>
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

namespace 
{

#pragma omp declare simd simdlen(SimdLen)
template<typename T, std::size_t reg_width, std::size_t SimdLen = reg_width / sizeof(T)>
static inline
constexpr T sin_horner_scheme(T x) noexcept 
{
    static_assert(reg_width == 32ull || reg_width == 64ull,
                  "Invalid argument passed: accepting only 32ull or 64ull.");
    static_assert(SimdLen > 0, "SimdLen must be greater than 0.");
    static_assert(std::is_floating_point<T>::value == true, 
                  "Invalid template argument: accepting only float types!!");

    constexpr T two_over_pi = static_cast<T>(0.6366197723675813430755350534900574);
    
    constexpr T pi_over_2_hi = static_cast<T>(1.57079637050628662109375);
    constexpr T pi_over_2_lo = static_cast<T>(-4.37113900018624283e-8);

    T scaled = x * two_over_pi;
    auto k = static_cast<std::int64_t>(scaled + (scaled >= T(0) ? T(0.5) : T(-0.5)));
    T kf = static_cast<T>(k);

    T r = (x - kf * pi_over_2_hi) - kf * pi_over_2_lo;
    T r2 = r * r;

    T result{};
    if ((k & 1) == 0) {
        // Even quadrant: evaluate sin(r) ≈ r + c3*r^3 + c5*r^5 + c7*r^7
        constexpr T c3 = static_cast<T>(-0.166666666666666324348);
        constexpr T c5 = static_cast<T>( 0.00833333333332248946124);
        constexpr T c7 = static_cast<T>(-0.000198412698298579493134);

        result = r * (T(1) + r2 * (c3 + r2 * (c5 + r2 * c7)));
    } else {
        // Odd quadrant: evaluate cos(r) ≈ 1 + c2*r^2 + c4*r^4 + c6*r^6
        constexpr T c2 = static_cast<T>(-0.499999999999999944489);
        constexpr T c4 = static_cast<T>( 0.0416666666666666019037);
        constexpr T c6 = static_cast<T>(-0.00138888888888741095749);

        result = T(1) + r2 * (c2 + r2 * (c4 + r2 * c6));
    }

    return (k & 2) ? -result : result;
}

// Reduces any angle in radians to the range [-pi/4, pi/4] for fast convergence
#pragma omp declare simd simdlen(SimdLen)
template<typename T, std::size_t reg_width, std::size_t SimdLen = reg_width / sizeof(T)>
static inline 
T approximated_tan(T x, std::int32_t terms = 10) noexcept
{
    static_assert(reg_width == 32ull || reg_width == 64ull,
                  "Invalid argument passed: accepting only 32ull or 64ull.");
    static_assert(SimdLen > 0, "SimdLen must be greater than 0.");
    static_assert(std::is_floating_point<T>::value == true, 
                  "Invalid template argument: accepting only float types!!");
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

/* For std::bit_cast substitute
   Adapted from: https://en.cppreference.com/cpp/numeric/bit_cast
*/
template<class To, class From>
std::enable_if_t<
    sizeof(To) == sizeof(From) &&
    std::is_trivially_copyable_v<From> &&
    std::is_trivially_copyable_v<To>,
    To>
// constexpr support needs compiler magic
bit_cast(const From& src) noexcept
{
    static_assert(std::is_trivially_constructible_v<To>,
        "This implementation additionally requires "
        "destination type to be trivially constructible");
    
    To dst;
    std::memcpy(&dst, &src, sizeof(To));
    return dst;
}

#pragma omp declare simd simdlen(SimdLen)
template<typename T, std::size_t reg_width, std::size_t SimdLen = reg_width / sizeof(T)>
static inline
T exp_estrin_scheme(T x) noexcept
{
    
    constexpr T MAX_X = static_cast<T>(709.782712893384);
    constexpr T MIN_X = static_cast<T>(-708.3964185322641);
    x = (x > MAX_X) ? MAX_X : ((x < MIN_X) ? MIN_X : x); // Compiles to branchless MINSD/MAXSD

    constexpr T INV_LN2 = static_cast<T>(1.44269504088896340736);
    constexpr T LN2_HI  = static_cast<T>(0.69314718036912381649); // High bits of ln(2)
    constexpr T LN2_LO  = static_cast<T>(1.90821492927058770002e-10); // Low bits for precision

    T kf = std::round(x * INV_LN2);
    std::int64_t k = static_cast<int64_t>(kf);

    T r = (x - kf * LN2_HI) - kf * LN2_LO;

    constexpr T c0  = static_cast<T>(1.0);
    constexpr T c1  = static_cast<T>(1.0);
    constexpr T c2  = static_cast<T>(0.5);
    constexpr T c3  = static_cast<T>(1.0 / 6.0);
    constexpr T c4  = static_cast<T>(1.0 / 24.0);
    constexpr T c5  = static_cast<T>(1.0 / 120.0);
    constexpr T c6  = static_cast<T>(1.0 / 720.0);
    constexpr T c7  = static_cast<T>(1.0 / 5040.0);
    constexpr T c8  = static_cast<T>(1.0 / 40320.0);
    constexpr T c9  = static_cast<T>(1.0 / 362880.0);
    constexpr T c10 = static_cast<T>(1.0 / 3628800.0);
    constexpr T c11 = static_cast<T>(1.0 / 39916800.0);

    T r2 = r * r;
    T r4 = r2 * r2;
    T r8 = r4 * r4;

    T p01  = std::fma(c1,  r, c0);
    T p23  = std::fma(c3,  r, c2);
    T p45  = std::fma(c5,  r, c4);
    T p67  = std::fma(c7,  r, c6);
    T p89  = std::fma(c9,  r, c8);
    T p1011 = std::fma(c11, r, c10);

    T p03  = std::fma(p23,   r2, p01);
    T p47  = std::fma(p67,   r2, p45);
    T p811 = std::fma(p1011, r2, p89);

    T p07  = std::fma(p47,  r4, p03);
    T poly = std::fma(p811, r8, p07);

    // Double exponent bias is 1023, shifted to bits [62:52]
    std::uint64_t exp_bits = static_cast<uint64_t>(k + 1023) << 52;
    T scale = bit_cast<T>(exp_bits);
    return poly * scale;
}

} // anonymous namespace 

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
                     T phi,T phi0) noexcept
{
    static_assert(reg_width == 32ull || reg_width == 64ull,
                  "Invalid argument passed: accepting only 32ull or 64ull.");
    static_assert(SimdLen > 0, "SimdLen must be greater than 0.");
    static_assert(std::is_floating_point<T>::value == true, 
                  "Invalid template argument: accepting only float types!!");
    static_assert(which_tan ==1 || which_tan ==2,
                  "Invalid template argument: accepting only values: 1 or 2.");
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
                     T phi,T phi0) noexcept
{
    static_assert(reg_width == 32ull || reg_width == 64ull,
                  "Invalid argument passed: accepting only 32ull or 64ull.");
    static_assert(SimdLen > 0, "SimdLen must be greater than 0.");
    static_assert(std::is_floating_point<T>::value == true, 
                  "Invalid template argument: accepting only float types!!");
    static_assert(which_tan ==1 || which_tan ==2,
                  "Invalid template argument: accepting only values: 1 or 2.");
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
                         const T dndr) noexcept
{
    static_assert(reg_width == 32ull || reg_width == 64ull,
                  "Invalid argument passed: accepting only 32ull or 64ull.");
    static_assert(SimdLen > 0, "SimdLen must be greater than 0.");
    static_assert(std::is_floating_point<T>::value == true, 
                  "Invalid template argument: accepting only float types!!");
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
                    const T dndr) noexcept
{   
    static_assert(reg_width == 32ull || reg_width == 64ull,
                  "Invalid argument passed: accepting only 32ull or 64ull.");
    static_assert(SimdLen > 0, "SimdLen must be greater than 0.");
    static_assert(std::is_floating_point<T>::value == true, 
                  "Invalid template argument: accepting only float types!!");
    constexpr T inv_erad = static_cast<T>(0.00015678896205707118218877391);
    T inv_rho,k_rel;
    inv_rho = static_cast<T>(1.0)/rad_ray_curvature_f251<T,reg_width>(n,z,dndr);
    k_rel   = inv_erad*inv_rho;
    return (k_rel);
}

/*
    formula 2.67, page: 52
*/
#pragma omp declare simd simdlen(SimdLen) notinbranch
template<typename T ,std::size_t reg_width, std::size_t SimdLen = reg_width / sizeof(T)>
static inline 
T rho_to_a_f267(const T dndh) noexcept
{
    static_assert(reg_width == 32ull || reg_width == 64ull,
                  "Invalid argument passed: accepting only 32ull or 64ull.");
    static_assert(SimdLen > 0, "SimdLen must be greater than 0.");
    static_assert(std::is_floating_point<T>::value == true, 
                  "Invalid template argument: accepting only float types!!");
    constexpr T inv_erad = static_cast<T>(0.00015678896205707118218877391);
    T R;
    R = inv_erad*dndh;
    return (R);
}

/*
   formula: 1.45, page 29
   dn0  coefficient of refreaction near the Earth surface i.e. dn0 = (240*10e-6->380*10e-6)
   beta coefficient describing the diminishing of 'n' as function of height, i.e. 0.10->0.14 1/km
   h    height 
*/
#pragma omp declare simd simdlen(SimdLen) notinbranch
template<typename T ,std::size_t reg_width, std::int32_t which_exp,std::size_t SimdLen = reg_width / sizeof(T)>
static inline 
T n_avg_h_f145(const T dn0,const T beta,
               const T h) noexcept
{
    static_assert(reg_width == 32ull || reg_width == 64ull,
                  "Invalid argument passed: accepting only 32ull or 64ull.");
    static_assert(SimdLen > 0, "SimdLen must be greater than 0.");
    static_assert(std::is_floating_point<T>::value == true, 
                  "Invalid template argument: accepting only float types!!");
    static_assert(which_exp == 1 || which_exp == 2,
                  "Invalid argument passed: accepting only 1 or 2.");
    T nah,earg;
    earg = beta*h;
    if constexpr(which_exp == 1)
    {
        nah  = static_cast<T>(1.0)+dn0*std::exp(-earg);
    }
    else if constexpr(which_exp == 2)
    {
        nah  = static_cast<T>(1.0)+dn0*exp_estrin_scheme<T,reg_width>(-earg);
    }
    return (nah);
}

/*
   formula 1.46, page: 29
   dn0  ! coefficient of refreaction near the Earth surface i.e. dn0 = (240*10e-6->380*10e-6)
*/
#pragma omp declare simd simdlen(SimdLen) notinbranch
template<typename T ,std::size_t reg_width, std::int32_t which_exp, std::size_t SimdLen = reg_width / sizeof(T)>
static inline 
T approx_beta_coeff_f146(const T dn0) noexcept
{
    static_assert(reg_width == 32ull || reg_width == 64ull,
                  "Invalid argument passed: accepting only 32ull or 64ull.");
    static_assert(SimdLen > 0, "SimdLen must be greater than 0.");
    static_assert(std::is_floating_point<T>::value == true, 
                  "Invalid template argument: accepting only float types!!");
    static_assert(which_exp == 1 || which_exp == 2,
                  "Invalid argument passed: accepting only 1 or 2.");
    T beta,t0,earg;
    t0   = static_cast<T>(0.00000732)/dn0; 
    earg = static_cast<T>(5577.0)*dn0; 
    if constexpr(which_exp == 1)
    {
        beta = t0*std::exp(earg);
    }
    else if constexpr(which_exp) 
    {
        beta = t0*exp_estrin_scheme<T,reg_width>(earg);
    }
    return (beta);
}

/*
      probability integral
*/
#pragma omp declare simd simdlen(SimdLen) notinbranch
template<typename T ,std::size_t reg_width, std::size_t SimdLen = reg_width / sizeof(T)>
static inline 
T prob_integral(const T x) noexcept
{
    static_assert(reg_width == 32ull || reg_width == 64ull,
                  "Invalid argument passed: accepting only 32ull or 64ull.");
    static_assert(SimdLen > 0, "SimdLen must be greater than 0.");
    static_assert(std::is_floating_point<T>::value == true, 
                  "Invalid template argument: accepting only float types!!");
    constexpr T C0707106781186547524400844362105 = static_cast<T>(0.707106781186547524400844362105);
    constexpr T half = static_cast<T>(0.5);
    T result;
    result = half*(static_cast<T>(1.0)*std::erf(x*C0707106781186547524400844362105));
    return (result);
}

/*
   formula 3.37, page: 68
*/
#pragma omp declare simd simdlen(SimdLen) notinbranch
template<typename T ,std::size_t reg_width, std::int32_t which_tan, std::int32_t which_exp,
         std::size_t SimdLen = reg_width / sizeof(T)>
static inline 
T  analytic_solution_L1_f337(const T beta,const T dn0,
                        T z0, const T H) noexcept
{
    static_assert(reg_width == 32ull || reg_width == 64ull,
                  "Invalid argument passed: accepting only 32ull or 64ull.");
    static_assert(SimdLen > 0, "SimdLen must be greater than 0.");
    static_assert(std::is_floating_point<T>::value == true, 
                  "Invalid template argument: accepting only float types!!");
     static_assert(which_tan == 1 || which_tan == 2,
                  "Invalid argument passed: accepting only: which_tan == 1 or which_tan == 2.");
    static_assert(which_exp == 1 || which_exp == 2,
                  "Invalid argument passed: accepting only: which_exp == 1 or which_exp == 2.");
    constexpr T one = static_cast<T>(1.0);
    constexpr T two = static_cast<T>(2.0);
    constexpr T a   = static_cast<T>(6378.0);
    T L1;
    T cosz0,ctgz0,ea1;
    T ea2,exp1,exp2,num2;
    T den2,num1,den1,sdn0;
    T stgz0,rat1,rat2,tgz0;
    ea1   = -two*beta*H;
    ea2   = -beta*H; 
    if constexpr(which_tan == 1)
    {
       tgz0  = std::tan(z0);
    }
    else if constexpr(which_tan == 2)
    {
       tgz0 = approximated_tan<T,reg_width>(z0); 
    }
    ctgz0 = one/tgz0; 
    sdn0  = dn0*dn0; 
    if constexpr(which_exp == 1)
    {
        exp1  = std::exp(ea1);
    }
    else if constexpr(which_exp == 2)
    {
        exp1  = exp_estrin_scheme<T,reg_width>(ea1);
    }
    num1  = beta*a*sdn0*ctgz0;
    cosz0 = std::cos(z0);
    den1  = cosz0*cosz0;
    if constexpr(which_exp == 1)
    { 
        exp2  = std::exp(ea2);
    }
    else if constexpr(which_exp == 2)
    {
        exp2  = exp_estrin_scheme<T,reg_width>(ea2);
    }
    rat1  = num1/den1; 
    stgz0 = two*(tgz0*tgz0); 
    den2  = std::sqrt(one+stgz0*(H/a));
    num2  = exp1-exp2; 
    rat2  = num2/den2; 
    L1    = rat1*rat2;
    return (L1); 
}

/*
    Formula: 3.41, page: 68
*/
#pragma omp declare simd simdlen(SimdLen) notinbranch
template<typename T ,std::size_t reg_width, std::int32_t which_tan, std::int32_t which_exp,
         std::size_t SimdLen = reg_width / sizeof(T)>
static inline 
T analytic_solution_L2_f341(const T dn0,const T beta,
                            T z0, const T H) noexcept
{
    static_assert(reg_width == 32ull || reg_width == 64ull,
                  "Invalid argument passed: accepting only 32ull or 64ull.");
    static_assert(SimdLen > 0, "SimdLen must be greater than 0.");
    static_assert(std::is_floating_point<T>::value == true, 
                  "Invalid template argument: accepting only float types!!");
     static_assert(which_tan == 1 || which_tan == 2,
                  "Invalid argument passed: accepting only: which_tan == 1 or which_tan == 2.");
    static_assert(which_exp == 1 || which_exp == 2,
                  "Invalid argument passed: accepting only: which_exp == 1 or which_exp == 2.");
    constexpr T one = static_cast<T>(1.0);
    constexpr T two = static_cast<T>(2.0);
    constexpr T half= static_cast<T>(0.5);
    constexpr T a   = static_cast<T>(6378.0);
    constexpr T C1253314137315500251207882642406 = static_cast<T>(1.253314137315500251207882642406);
    T L2;
    T sba, ctgz0;
    T sctgz0, tbH, phi1, phi2; 
    T exp1, bactgz0, t0, t1;
    sba    = std::sqrt(beta*a);
    if constexpr(which_tan == 1)
    {
        ctgz0  = one/std::tan(z0);
    }
    else if constexpr(which_tan == 2)
    {
        ctgz0  = one/approximated_tan<T,reg_width>(z0);
    }
    sctgz0 = ctgz0*ctgz0; 
    bactgz0= beta*a*sctgz0; 
    tbH    = two*beta*H; 
    t0     = dn0*std::qrt(beta*a*ctgz0);
    if constexpr(which_exp == 1)
    {
         exp1 = std::exp(sctgz0*half)*C1253314137315500251207882642406;
    }  
    else if constexpr(which_exp == 2)
    {
         exp1 = exp_estrin_scheme<T,reg_width>(sctgz0*half)*C1253314137315500251207882642406;
    }             
    phi1   = prob_integral<T,reg_width>(std::sqrt(bactgz0*tbH));
    phi2   = prob_integral<T,reg_width>(std::sqrt(bactgz0));
    t1     = phi1-phi2; 
    L2     = t0*exp1*t1;
    return (L2);
}

/*
    Formula: 3.42, page: 68
*/
#pragma omp declare simd simdlen(SimdLen) notinbranch
template<typename T ,std::size_t reg_width, std::int32_t which_tan, std::int32_t which_exp,
         std::size_t SimdLen = reg_width / sizeof(T)>
static inline 
T analytic_solution_L3_f342(const T dn0,const T beta,
                            T z0, const T H) noexcept
{
    static_assert(reg_width == 32ull || reg_width == 64ull,
                  "Invalid argument passed: accepting only 32ull or 64ull.");
    static_assert(SimdLen > 0, "SimdLen must be greater than 0.");
    static_assert(std::is_floating_point<T>::value == true, 
                  "Invalid template argument: accepting only float types!!");
    static_assert(which_tan == 1 || which_tan == 2,
                  "Invalid argument passed: accepting only: which_tan == 1 or which_tan == 2.");
    static_assert(which_exp == 1 || which_exp == 2,
                  "Invalid argument passed: accepting only: which_exp == 1 or which_exp == 2.");
    constexpr T one = static_cast<T>(1.0);
    constexpr T two = static_cast<T>(2.0);
    constexpr T four= static_cast<T>(4.5);
    constexpr T a   = static_cast<T>(6378.0);
    constexpr T C1253314137315500251207882642406 = static_cast<T>(1.253314137315500251207882642406);
    T L3;
    T sba, ctgz0;
    T sctgz0, tbH, phi1, phi2; 
    T exp1, bactgz0, t0, t1;
    sba   = std::sqrt(two*beta*a);
    if constexpr(which_tan == 1)
    {
        ctgz0  = one/std::tan(z0);
    }
    else if constexpr(which_tan == 2)
    {
        ctgz0  = one/approximated_tan<T,reg_width>(z0);
    }
    sctgz0 = ctgz0*ctgz0; 
    bactgz0= two*beta*a*sctgz0; 
    tbH    = four*beta*H; 
    t0     = dn0*std::sqrt(beta*a*ctgz0);
    if constexpr(which_exp == 1)
    {
        exp1 = std::exp(sctgz0)*C1253314137315500251207882642406;
    }
    else if constexpr(which_exp == 2)
    {
        exp1 = exp_estrin_scheme<T,reg_width>(sctgz0)*C1253314137315500251207882642406;
    }               
    phi1   = prob_integral<T,reg_width>(std::sqrt(bactgz0+tbH));
    phi2   = prob_integral<T,reg_width>(std::sqrt(bactgz0));
    t1     = phi1-phi2; 
    L3     = t0*exp1*t1;
    return (L3);
}

/*
   Formula: 3.45, computing the angle of Refraction (0 < z0 < 90(deg)) for
   arbitrary value of n(h)
*/
#pragma omp declare simd simdlen(SimdLen) notinbranch
template<typename T ,std::size_t reg_width, std::int32_t which_sin, std::int32_t which_tan,
         std::int32_t which_exp, std::size_t SimdLen = reg_width / sizeof(T)>
static inline 
T refraction_angle_f345(const T n0, const T nh,
                        T z0, const T dn0, const T beta,
                        const T H) noexcept
{
    static_assert(reg_width == 32ull || reg_width == 64ull,
                  "Invalid argument passed: accepting only 32ull or 64ull.");
    static_assert(SimdLen > 0, "SimdLen must be greater than 0.");
    static_assert(std::is_floating_point<T>::value == true, 
                  "Invalid template argument: accepting only float types!!");
    static_assert(which_tan == 1 || which_tan == 2,
                  "Invalid argument passed: accepting only: which_tan == 1 or which_tan == 2.");
    static_assert(which_sin == 1 || which_sin == 2,
                  "Invalid argument passed: accepting only: which_sin == 1 or which_sin == 2.");
    static_assert(which_exp == 1 || which_exp == 2,
                  "Invalid argument passed: accepting only: which_exp == 1 or which_exp == 2.");
    constexpr T a = static_cast<T>(6378.0);
    constexpr T zero = static_cast<T>(0.0);
    constexpr T one = static_cast<T>(1.0);
    T alpha;
    T ctgz0, ln0nh, ssecz0,badn0, ctgzsec0;
    T t0, t1, t2; 
    T L1, L2, L3;
    badn0  = beta*a*dn0;
    L1     = zero;
    if constexpr(which_tan == 1)
    { 
        ctgz0  = one/std::tan(z0);
    }
    else if constexpr(which_tan == 2)
    {
        ctgz0  = one/approximated_tan<T,reg_width>(z0);
    }
    L2       = zero; 
    ln0nh    = std::log(n0/nh); 
    L3       = zero; 
    if constexpr(which_sin == 1)
    {
        t0 = one/std::sin(z0);
    }
    else if constexpr(which_sin == 2)
    {
        t0 = one/sin_horner_scheme<T,reg_width>(z0);
    }
    ssecz0   = t0*t0;
    L1       = analytic_solution_L1_f337<T,reg_width,which_tan,which_exp>(dn0,beta,z0,H);
    t0       = -ctgz0*ln0nh+L1; 
    ctgzsec0 = ctgz0*ssecz0;
    L2       = analytic_solution_L2_f341<T,reg_width,which_tan,which_exp>(dn0,beta,z0,H0);
    t1       = ctgzsec0*L2;
    L3       = analytic_solution_L3_f342<T,reg_width,which_tan,which_exp>(dn0,beta,z0,H);
    t2       = badn0*ctgzsec0*(L3-L2);
    alpha    = t0+t1+t2;
    return (alpha);
}

/*
       ! The angle of arrival close to horizon.
       ! formula 3.51, page: 70
       ! analytic solution L2 for angle near 90 (deg)
*/
#pragma omp declare simd simdlen(SimdLen) notinbranch
template<typename T ,std::size_t reg_width, std::int32_t which_exp, std::int32_t which_tan,
         std::size_t SimdLen = reg_width / sizeof(T)>
static inline 
T analytic_solution_n90_L2_f351(const T dn0, const T beta, T z0) noexcept
{
    static_assert(reg_width == 32ull || reg_width == 64ull,
                  "Invalid argument passed: accepting only 32ull or 64ull.");
    static_assert(SimdLen > 0, "SimdLen must be greater than 0.");
    static_assert(std::is_floating_point<T>::value == true, 
                  "Invalid template argument: accepting only float types!!");
    static_assert(which_tan == 1 || which_tan == 2,
                  "Invalid argument passed: accepting only: which_tan == 1 or which_tan == 2.");
    static_assert(which_exp == 1 || which_exp == 2,
                  "Invalid argument passed: accepting only: which_sin == 1 or which_sin == 2.");
    constexpr T a = static_cast<T>(6378.0);
    constexpr T two = static_cast<T>(2.0);
    constexpr T one = static_cast<T>(1.0);
    constexpr T C1253314137315500251207882642406 = static_cast<T>(1.253314137315500251207882642406);
    constexpr T C0318309886183790671537767526745 = static_cast<T>(0.318309886183790671537767526745);
    T sba, tgz0, ctgz0;
    T earg, exp1, t0, t1, strm;
    T L2;
    sba  = std::sqrt(beta*a);
    if constexpr(which_tan == 1)
    {
        tgz0 = std::tan(z0);
    }
    else if constexpr(which_tan == 2)
    {
        tgz0 = approximated_tan<T,reg_width>(z0);
    }
    ctgz0= one/tgz0;
    earg = beta*a/(two*tgz0*tgz0);
    if constexpr(which_exp == 1)
    {
        exp1 = std::exp(earg);
    }
    else if constexpr(which_exp == 2)
    {
        exp1 = exp_estrin_scheme<T,reg_width>(earg);
    }
    t0   = dn0*(sba/tgz0)*exp1; 
    strm = std::sqrt(two*beta*a*C0318309886183790671537767526745);
    t1   = C1253314137315500251207882642406*(one-strm*ctgz0);                   
    L2   = t0*t1;
    return (L2);
}

/*
       ! The angle of arrival close to horizon.
       ! formula 3.51, page: 70
       ! analytic solution L3 for angle near 90 (deg)
*/
#pragma omp declare simd simdlen(SimdLen) notinbranch
template<typename T ,std::size_t reg_width, std::int32_t which_exp, std::int32_t which_tan,
         std::size_t SimdLen = reg_width / sizeof(T)>
static inline 
T analytic_solution_n90_L3_f351(const T dn0, const T beta, T z0) noexcept
{
    static_assert(reg_width == 32ull || reg_width == 64ull,
                  "Invalid argument passed: accepting only 32ull or 64ull.");
    static_assert(SimdLen > 0, "SimdLen must be greater than 0.");
    static_assert(std::is_floating_point<T>::value == true, 
                  "Invalid template argument: accepting only float types!!");
    static_assert(which_tan == 1 || which_tan == 2,
                  "Invalid argument passed: accepting only: which_tan == 1 or which_tan == 2.");
    static_assert(which_exp == 1 || which_exp == 2,
                  "Invalid argument passed: accepting only: which_exp == 1 or which_sin == 2.");
    constexpr T a = static_cast<T>(6378.0);
    constexpr T two = static_cast<T>(2.0);
    constexpr T one = static_cast<T>(1.0);
    constexpr T four = static_cast<T>(4.0);
    constexpr T C1253314137315500251207882642406 = static_cast<T>(1.253314137315500251207882642406);
    constexpr T C0318309886183790671537767526745 = static_cast<T>(0.318309886183790671537767526745);
    T sba, tgz0, ctgz0;
    T earg, exp1, t0, t1, strm;
    T L3;
    sba  = std::sqrt(two*beta*a);
    if constexpr(which_tan == 1)
    {
        tgz0 = std::tan(z0);
    }
    else if constexpr(which_tan == 2)
    {
        tgz0 = approximated_tan<T,reg_width>(z0);
    }
    ctgz0 = one/tgz0;
    earg = beta*a/(tgz0*tgz0);
    if constexpr(which_exp == 1)
    {
        exp1 = std::exp(earg);
    }
    else if constexpr(which_exp == 2)
    {
        exp1 = exp_estrin_scheme<T,reg_width>(earg);
    }       
    t0   = dn0*(sba/tgz0)*exp1; 
    strm = std::sqrt(four*beta*a*C0318309886183790671537767526745);
    t1   = C1253314137315500251207882642406* (one-strm*ctgz0);          
    L3   = t0*t1;
    return (L3);
}

/*
    ! The angle of arrival close to horizon.
    ! formula 3.51, page: 70
    ! The whole solution for angle alpha near 90 (deg)
*/
#pragma omp declare simd simdlen(SimdLen) notinbranch
template<typename T ,std::size_t reg_width, std::int32_t which_exp, std::int32_t which_tan,
         std::size_t SimdLen = reg_width / sizeof(T)>
static inline 
T  refraction_angle_n90_f351(const T dn0,const T beta,
                             T z0) noexcept
{
    static_assert(reg_width == 32ull || reg_width == 64ull,
                  "Invalid argument passed: accepting only 32ull or 64ull.");
    static_assert(SimdLen > 0, "SimdLen must be greater than 0.");
    static_assert(std::is_floating_point<T>::value == true, 
                  "Invalid template argument: accepting only float types!!");
    static_assert(which_tan == 1 || which_tan == 2,
                  "Invalid argument passed: accepting only: which_tan == 1 or which_tan == 2.");
    static_assert(which_exp == 1 || which_exp == 2,
                  "Invalid argument passed: accepting only: which_exp == 1 or which_sin == 2.");
    constexpr T a   = static_cast<T>(6378.0);
    constexpr T one = static_cast<T>(1.0);
    T ctgz0, badn0, cosz0, scosz0;
    T L2, L3, t0, t1, rat;
    T alpha;
    cosz0 = std::cos(z0);
    badn0 = beta*dn0*a;
    if constexpr(which_tan == 1)
    {
        ctgz0 = one/std::tan(z0);
    }
    else if constexpr(which_tan == 2)
    {
        ctgz0 = one/approximated_tan<T,reg_width>(z0);
    }
    scosz0= cosz0*cosz0; 
    L2  = analytic_solution_n90_L2_f351<T,reg_width,which_exp,which_tan>(dn0,beta,z0);
    rat = ctgz0/scosz0; 
    t0  = -dn0*ctgz0+(one-badn0);
    L3  = analytic_solution_n90_L3_f351<T,reg_width,which_exp,which_tan>(dn0,beta,z0);
    t1  = rat*L2+badn0*rat*L3; 
    alpha = t0*t1;
    return (alpha); 
}

}

}

#endif /*__GMS_ANALYTIC_REFRACTION_H__*/