

/*
    Adapted from Fortran implementation.
    Original author description:
    ! <CONTACT EMAIL="Bruce.Wyman@noaa.gov">
    !   Bruce Wyman
    ! </CONTACT>

       ! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

   ! <OVERVIEW>
       !   Routines for creating Gaussian-shaped land surface topography
       !   for latitude-longitude grids.
   ! </OVERVIEW>
*/
/*
WRF Model Version 3.8.1 (August 12, 2016)
http://wrf-model.org/users/users.php

------------------------
WRF PUBLIC DOMAIN NOTICE
------------------------

WRF was developed at the National Center for Atmospheric Research
(NCAR) which is operated by the University Corporation for
Atmospheric Research (UCAR). NCAR and UCAR make no proprietary
claims, either statutory or otherwise, to this version and
release of WRF and consider WRF to be in the public domain for
use by any person or entity for any purpose without any fee or
charge. UCAR requests that any WRF user include this notice on
any partial or full copies of WRF. WRF is provided on an "AS
IS" basis and any warranties, either express or implied,
including but not limited to implied warranties of
non-infringement, originality, merchantability and fitness for a
particular purpose, are disclaimed. In no event shall
UCAR be liable for any damages, whatsoever, whether direct,
indirect, consequential or special, that arise out of or in
connection with the access, use or performance of WRF, including
infringement actions.
*/
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


#include <algorithm>
#include <cmath>
#include <omp.h>
#include "GMS_idealized_topography.h"
#include "GMS_indices.h"

void gms::radiolocation::gaussian_topography_r4(
    float *__restrict __ATTR_ALIGN__(64) lon, const int32_t nlon,
    float *__restrict __ATTR_ALIGN__(64) lat, const int32_t nlat,
    float *__restrict __ATTR_ALIGN__(64) zsurf, const float height) {

  if (__builtin_expect(lon <= 0, 0) || __builtin_expect(lat <= 0, 0)) {
    return;
  }

  constexpr float tpi = 6.2831853071795864769253f; //! 2.0_sp*pi
  constexpr float dtr = 0.0174532925199432957692f; //! tpi/360.
  constexpr float olon = 90.0f * dtr;
  constexpr float olat = 45.0f * dtr;
  constexpr float wlon = 15.0f * dtr;
  constexpr float wlat = wlon;
  constexpr float rlon = 0.0f;
  constexpr float rlat = 0.0f;
  constexpr float inwlat = 0.066666666666666666666666666667f;
  constexpr float inwlon = inwlat;
  float xx = 0.0f;
  float yy = 0.0f;
#if defined(__INTEL_COMPILER) || defined(__ICC)
  assume_aligned(lon, 64);
  assume_aligned(lat, 64);
  assume_aligned(zsurf, 64)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
  lon = (float *)__builtin_assume_aligned(lon, 64);
  lat = (float *)__builtin_assume_aligned(lat, 64);
  zsurf = (float *)__builtin_assume_aligned(zsurf, 64);
#endif
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#pragma prefetch lat : 0 : 4
#pragma prefetch lat : 1 : 16
#pragma prefetch lon : 0 : 4
#pragma prefetch lon : 1 : 16
#endif
#pragma omp simd simdlen(4) linear(i : 1) private(ilon, dx, xx)
      for (int32_t i = 0; i != nlon; ++i) {
    const float ilon = lon[i];
    float dx = std::abs(ilon - olon);
    dx = std::min(dx, std::abs(dx - tpi));
    xx = std::max(0.0f, dx - rlon) * inwlon;
#pragma omp simd simdlen(4) linear(j : 1) private(ilat, dy, yy)
    for (int32_t j = 0; j != nlat; ++j) {
      const float ilat = lat[j];
      const float dy = std::abs(ilat - olat);
      yy = std::max(0.0f, dy - rlat) * inwlat;
      zsurf[Ix2D(i, nlat, j)] = height * std::exp(-xx * xx - yy * yy);
    }
  }
}

void gms::radiolocation::gaussian_topography_r8(
    double *__restrict __ATTR_ALIGN__(64) lon, const int32_t nlon,
    double *__restrict __ATTR_ALIGN__(64) lat, const int32_t nlat,
    double *__restrict __ATTR_ALIGN__(64) zsurf, const double height) {

  if (__builtin_expect(lon <= 0, 0) || __builtin_expect(lat <= 0, 0)) {
    return;
  }

  constexpr double tpi = 6.2831853071795864769253; //! 2.0_sp*pi
  constexpr double dtr = 0.0174532925199432957692; //! tpi/360.
  constexpr double olon = 90.0 * dtr;
  constexpr double olat = 45.0 * dtr;
  constexpr double wlon = 15.0 * dtr;
  constexpr double wlat = wlon;
  constexpr double rlon = 0.0;
  constexpr double rlat = 0.0;
  constexpr double inwlat = 0.066666666666666666666666666667;
  constexpr double inwlon = inwlat;
  double xx = 0.0;
  double yy = 0.0;
#if defined(__INTEL_COMPILER) || defined(__ICC)
  assume_aligned(lon, 64);
  assume_aligned(lat, 64);
  assume_aligned(zsurf, 64)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
  lon = (double *)__builtin_assume_aligned(lon, 64);
  lat = (double *)__builtin_assume_aligned(lat, 64);
  zsurf = (double *)__builtin_assume_aligned(zsurf, 64);
#endif
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#pragma prefetch lat : 0 : 4
#pragma prefetch lat : 1 : 16
#pragma prefetch lon : 0 : 4
#pragma prefetch lon : 1 : 16
#endif
#pragma omp simd simdlen(8) linear(i : 1) private(ilon, dx, xx)
      for (int32_t i = 0; i != nlon; ++i) {
    const double ilon = lon[i];
    double dx = std::abs(ilon - olon);
    dx = std::min(dx, std::abs(dx - tpi));
    xx = std::max(0.0, dx - rlon) * inwlon;
#pragma omp simd simdlen(8) linear(j : 1) private(ilat, dy, yy)
    for (int32_t j = 0; j != nlat; ++j) {
      const double ilat = lat[j];
      const double dy = std::abs(ilat - olat);
      yy = std::max(0.0, dy - rlat) * inwlat;
      zsurf[Ix2D(i, nlat, j)] = height * std::exp(-xx * xx - yy * yy);
    }
  }
}

void gms::radiolocation::sinusoidal_topography(
    float *__restrict __ATTR_ALIGN__(64) lon, const int32_t nlon,
    float *__restrict __ATTR_ALIGN__(64) lat, const int32_t nlat,
    float *__restrict __ATTR_ALIGN__(64) zsurf, const float height_sin,
    const float m, const float Amp2, const float uneven_fac,
    const float deltalat, const bool uneven_sin) {

  if (__builtin_expect(nlon <= 0, 0) || __builtin_expect(nlat <= 0, 0)) {
    return;
  }
  if (__builtin_expect(height_sin == 0.0f, 0)) {
    return;
  }
  constexpr float pi = 3.1415926535897932384626f;
  constexpr float tpi4 = 2.356194490192344928847f;
  constexpr float spi4 = 5.4977871437821381673096f;
  constexpr float tpi = 6.2831853071795864769253f; //! 2.0*pi
  constexpr float dtr = 0.0174532925199432957692f; //! tpi/360.
  float lat00 = 25.0f + deltalat;
  lat00 = dtr * lat00;
  float lat11 = 65.0f + deltalat;
  lat11 = dtr * lat11;
  float ldiff = lat11 - lat00;
  float ildiff = 1.0f / ldiff;
  float latj, loni;
#if defined(__INTEL_COMPILER) || defined(__ICC)
  assume_aligned(lon, 64);
  assume_aligned(lat, 64);
  assume_aligned(zsurf, 64)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
  lon = (float *)__builtin_assume_aligned(lon, 64);
  lat = (float *)__builtin_assume_aligned(lat, 64);
  zsurf = (float *)__builtin_assume_aligned(zsurf, 64);
#endif
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#pragma prefetch lat : 0 : 4
#pragma prefetch lat : 1 : 16
#pragma prefetch lon : 0 : 4
#pragma prefetch lon : 1 : 16
#endif
      for (int32_t j = 0; j != nlat; ++j) {
    latj = lat[j];
#pragma simd simdlen(4) linear(i : 1) private(loni, t0, t1, t2)
    for (int32_t i = 0; i != nlon; ++i) {
      if (latj > lat00 && latj < lat11) {
        loni = lon[i];
        const float t0 = std::sin((latj - lat00) * pi) * ildiff;
        const float t1 = height_sin * t0 * t0;
        const float t2 = std::cos(m * loni) + Amp2 * std::cos(2.0f * loni);
        zsurf[Ix2D(j, nlon, i)] = t1 * t0;
        if (uneven_sin) {
          /*
             ! uneven_sig = .true.; one mountain taller than the other by
             uneven_fac - half of domain taller; ! this can only be used for
             Amp2=0; m=2; and height_sin > 0.
          */
          if (loni > tpi4 && loni < spi4) {
            zsurf[Ix2D(j, nlon, i)] *= uneven_fac;
          }
        }
      }
    }
  }
}

void gms::radiolocation::sinusoidal_topography(
    double *__restrict __ATTR_ALIGN__(64) lon, const int32_t nlon,
    double *__restrict __ATTR_ALIGN__(64) lat, const int32_t nlat,
    double *__restrict __ATTR_ALIGN__(64) zsurf, const double height_sin,
    const double m, const double Amp2, const double uneven_fac,
    const double deltalat, const bool uneven_sin) {

  if (__builtin_expect(nlon <= 0, 0) || __builtin_expect(nlat <= 0, 0)) {
    return;
  }
  if (__builtin_expect(height_sin == 0.0f, 0)) {
    return;
  }
  constexpr double pi = 3.1415926535897932384626;
  constexpr double tpi4 = 2.356194490192344928847;
  constexpr double spi4 = 5.4977871437821381673096;
  constexpr double tpi = 6.2831853071795864769253; //! 2.0*pi
  constexpr double dtr = 0.0174532925199432957692; //! tpi/360.
  double lat00 = 25.0 + deltalat;
  lat00 = dtr * lat00;
  double lat11 = 65.0 + deltalat;
  lat11 = dtr * lat11;
  double ldiff = lat11 - lat00;
  double ildiff = 1.0 / ldiff;
  double latj, loni;
#if defined(__INTEL_COMPILER) || defined(__ICC)
  assume_aligned(lon, 64);
  assume_aligned(lat, 64);
  assume_aligned(zsurf, 64)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
  lon = (double *)__builtin_assume_aligned(lon, 64);
  lat = (double *)__builtin_assume_aligned(lat, 64);
  zsurf = (double *)__builtin_assume_aligned(zsurf, 64);
#endif
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#pragma prefetch lat : 0 : 4
#pragma prefetch lat : 1 : 16
#pragma prefetch lon : 0 : 4
#pragma prefetch lon : 1 : 16
#endif
      for (int32_t j = 0; j != nlat; ++j) {
    latj = lat[j];
#pragma simd simdlen(8) linear(i : 1) private(loni, t0, t1, t2)
    for (int32_t i = 0; i != nlon; ++i) {
      if (latj > lat00 && latj < lat11) {
        loni = lon[i];
        const double t0 = std::sin((latj - lat00) * pi) * ildiff;
        const double t1 = height_sin * t0 * t0;
        const double t2 = std::cos(m * loni) + Amp2 * std::cos(2.0f * loni);
        zsurf[Ix2D(j, nlon, i)] = t1 * t0;
        if (uneven_sin) {
          /*
             ! uneven_sig = .true.; one mountain taller than the other by
             uneven_fac - half of domain taller; ! this can only be used for
             Amp2=0; m=2; and height_sin > 0.
          */
          if (loni > tpi4 && loni < spi4) {
            zsurf[Ix2D(j, nlon, i)] *= uneven_fac;
          }
        }
      }
    }
  }
}

gms::radiolocation
::terrain_topography_avx512_t
::terrain_topography_avx512_t() noexcept(true)
    : 
m_ids{},m_ide{},m_jds{},m_jde{},m_kds{},m_kde{},
m_ims{},m_ime{},m_jms{},m_jme{},m_kms{},m_kme{},
m_its{},m_ite{},m_jts{},m_jte{},m_kts{},m_kte{},
m_ideal_terrain{},m_spec_bdy_width{},
m_flat_plate_ht{},m_dx{},m_dy{},
m_phb{},m_ph{},m_ibm_ht_u{},m_ibm_ht_v{},m_ibm_ht_w{},
m_ibm_ht_c{},m_ibm_z0 
{ 
  
}

