

/*
    Adapted from Fortran implementation.
    Original author description:
    ! <CONTACT EMAIL="Bruce.Wyman@noaa.gov">
    !   Bruce Wyman
    ! </CONTACT>

       ! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

   ! <OVERVIEW>
       !   Routines for creating Gaussian-shaped land surface topography
       !   for latitude-longitude grm_ids.
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
clam_ims, either statutory or otherwise, to this version and
release of WRF and consm_ider WRF to be in the public domain for
use by any person or entity for any purpose without any fee or
charge. UCAR requests that any WRF user include this notice on
any partial or full copies of WRF. WRF is provm_ided on an "AS
IS" basis and any warranties, either express or implied,
including but not limm_ited to implied warranties of
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
 * along with this program; if not, wrm_ite to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/


#include <algorithm>
#include <cmath>
#include <omp.h>
#include <iostream>
#include "GMS_idealized_topography.h"
#include "GMS_indices.h"

void gms::radiolocation::gaussian_topography_r4(
    float *__restrict __ATTR_ALIGN__(64) lon, const std::int32_t nlon,
    float *__restrict __ATTR_ALIGN__(64) lat, const std::int32_t nlat,
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
  float ilon,ilat;
  float dx,dy;
  std::int32_t i,j;
#if defined(__INTEL_COMPILER) || defined(__ICC)
  __assume_aligned(lon, 64);
  __assume_aligned(lat, 64);
  __assume_aligned(zsurf, 64);
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
  for (i = 0; i != nlon; ++i) {
    ilon = lon[i];
    dx = std::abs(ilon - olon);
    dx = std::min(dx, std::abs(dx - tpi));
    xx = std::max(0.0f, dx - rlon) * inwlon;
#pragma omp simd simdlen(4) linear(j : 1) private(ilat, dy, yy)
    for (j = 0; j != nlat; ++j) {
      ilat = lat[j];
      dy = std::abs(ilat - olat);
      yy = std::max(0.0f, dy - rlat) * inwlat;
      zsurf[Ix2D(i, nlat, j)] = height * std::exp(-xx * xx - yy * yy);
    }
  }
}

void gms::radiolocation::gaussian_topography_r8(
    double *__restrict __ATTR_ALIGN__(64) lon, const std::int32_t nlon,
    double *__restrict __ATTR_ALIGN__(64) lat, const std::int32_t nlat,
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
  double ilon,ilat;
  double dx,dy;
  std::int32_t i,j;
#if defined(__INTEL_COMPILER) || defined(__ICC)
  __assume_aligned(lon, 64);
  __assume_aligned(lat, 64);
  __assume_aligned(zsurf, 64);
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
  for ( i = 0; i != nlon; ++i) {
    ilon = lon[i];
    dx = std::abs(ilon - olon);
    dx = std::min(dx, std::abs(dx - tpi));
    xx = std::max(0.0, dx - rlon) * inwlon;
#pragma omp simd simdlen(8) linear(j : 1) private(ilat, dy, yy)
    for ( j = 0; j != nlat; ++j) {
      ilat = lat[j];
      dy = std::abs(ilat - olat);
      yy = std::max(0.0, dy - rlat) * inwlat;
      zsurf[Ix2D(i, nlat, j)] = height * std::exp(-xx * xx - yy * yy);
    }
  }
}

void gms::radiolocation::sinusoidal_topography(
    float *__restrict __ATTR_ALIGN__(64) lon, const std::int32_t nlon,
    float *__restrict __ATTR_ALIGN__(64) lat, const std::int32_t nlat,
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
  float t0,t1,t2;
  std::int32_t i,j;
#if defined(__INTEL_COMPILER) || defined(__ICC)
  __assume_aligned(lon, 64);
  __assume_aligned(lat, 64);
  __assume_aligned(zsurf, 64);
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
  for (j = 0; j != nlat; ++j) {
    latj = lat[j];
#pragma simd simdlen(4) linear(i : 1) private(loni, t0, t1, t2)
    for ( i = 0; i != nlon; ++i) {
      if (latj > lat00 && latj < lat11) {
        loni = lon[i];
        t0 = std::sin((latj - lat00) * pi) * ildiff;
        t1 = height_sin * t0 * t0;
        t2 = std::cos(m * loni) + Amp2 * std::cos(2.0f * loni);
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
    double *__restrict __ATTR_ALIGN__(64) lon, const std::int32_t nlon,
    double *__restrict __ATTR_ALIGN__(64) lat, const std::int32_t nlat,
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
  double t0,t1,t2;
  std::int32_t i,j;
#if defined(__INTEL_COMPILER) || defined(__ICC)
  __assume_aligned(lon, 64);
  __assume_aligned(lat, 64);
  __assume_aligned(zsurf, 64);
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
  for (j = 0; j != nlat; ++j) {
    latj = lat[j];
#pragma simd simdlen(8) linear(i : 1) private(loni, t0, t1, t2)
    for ( i = 0; i != nlon; ++i) {
      if (latj > lat00 && latj < lat11) {
        loni = lon[i];
        t0 = std::sin((latj - lat00) * pi) * ildiff;
        t1 = height_sin * t0 * t0;
        t2 = std::cos(m * loni) + Amp2 * std::cos(2.0f * loni);
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
::terrain_topography_avx512_t(const terrain_topography_args_t &args) noexcept(false)
    :
m_ids(args.ids),m_ide(args.ide),m_jds(args.jds),
m_jde(args.jde),m_kds(args.kds),m_kde(args.kde),
m_ims(args.ims),m_ime(args.ime),m_jms(args.jms),
m_jme(args.jme),m_kms(args.kms),m_kme(args.kme),
m_its(args.its),m_ite(args.ite),m_jts(args.jts),
m_jte(args.jte),m_kts(args.kts),m_kte(args.kte),
m_ideal_terrain(args.ideal_terrain),
m_spec_bdy_width(args.spec_bdy_width),
m_flat_plate_ht(args.flat_plate_ht),
m_dx(args.dx),m_dy(args.dy),
m_phb(boost::extents[m_jme-m_jms][m_kme-m_kms][m_ime-m_ims]),
m_ph(boost::extents[m_jme-m_jms][m_kme-m_kms][m_ime-m_ims]),
m_ibm_ht_u(boost::extents[m_jme-m_jms][m_ime-m_ims]),
m_ibm_ht_v(boost::extents[m_jme-m_jms][m_ime-m_ims]),
m_ibm_ht_w(boost::extents[m_jme-m_jms][m_ime-m_ims]),
m_ibm_ht_c(boost::extents[m_jme-m_jms][m_ime-m_ims]),
m_ibm_z0(boost::extents[m_jme-m_jms][m_ime-m_ims])
{}


/*
   The dimensions of the terrain topography must be set by the caller
   This setting is not validated insm_ide the body of below stated Constructor.
*/
gms::radiolocation
::terrain_topography_avx512_t
::terrain_topography_avx512_t(const terrain_topography_args_t &args,
                              const common::Aligned64Mat3DF32 &phb,
                              const common::Aligned64Mat3DF32 &ph,
                              const common::Aligned64Mat2DF32 &ibm_ht_u,
                              const common::Aligned64Mat2DF32 &ibm_ht_v,
                              const common::Aligned64Mat2DF32 &ibm_ht_w,
                              const common::Aligned64Mat2DF32 &ibm_ht_c,
                              const common::Aligned64Mat2DF32 &ibm_z0) noexcept(false) 
    : 
m_ids(args.ids),m_ide(args.ide),m_jds(args.jds),
m_jde(args.jde),m_kds(args.kds),m_kde(args.kde),
m_ims(args.ims),m_ime(args.ime),m_jms(args.jms),
m_jme(args.jme),m_kms(args.kms),m_kme(args.kme),
m_its(args.its),m_ite(args.ite),m_jts(args.jts),
m_jte(args.jte),m_kts(args.kts),m_kte(args.kte),
m_ideal_terrain(args.ideal_terrain),
m_spec_bdy_width(args.spec_bdy_width),
m_flat_plate_ht(args.flat_plate_ht),
m_dx(args.dx),m_dy(args.dy),
m_phb(phb),m_ph(ph),m_ibm_ht_u(ibm_ht_u),
m_ibm_ht_v(ibm_ht_v), m_ibm_ht_w(ibm_ht_w),
m_ibm_ht_c(ibm_ht_c), m_ibm_z0(ibm_z0)
{}


/*
   The dimensions of the terrain topography must be set by the caller
   This setting is not validated insm_ide the body of below stated Constructor.
*/
gms::radiolocation
::terrain_topography_avx512_t
::terrain_topography_avx512_t(const terrain_topography_args_t &args,
                              common::Aligned64Mat3DF32 &&phb,
                              common::Aligned64Mat3DF32 &&ph,
                              common::Aligned64Mat2DF32 &&ibm_ht_u,
                              common::Aligned64Mat2DF32 &&ibm_ht_v,
                              common::Aligned64Mat2DF32 &&ibm_ht_w,
                              common::Aligned64Mat2DF32 &&ibm_ht_c,
                              common::Aligned64Mat2DF32 &&ibm_z0) noexcept(true) 
    : 
m_ids(args.ids),m_ide(args.ide),m_jds(args.jds),
m_jde(args.jde),m_kds(args.kds),m_kde(args.kde),
m_ims(args.ims),m_ime(args.ime),m_jms(args.jms),
m_jme(args.jme),m_kms(args.kms),m_kme(args.kme),
m_its(args.its),m_ite(args.ite),m_jts(args.jts),
m_jte(args.jte),m_kts(args.kts),m_kte(args.kte),
m_ideal_terrain(args.ideal_terrain),
m_spec_bdy_width(args.spec_bdy_width),
m_flat_plate_ht(args.flat_plate_ht),
m_dx(args.dx),m_dy(args.dy),
m_phb(std::move(phb)),m_ph(std::move(ph)),m_ibm_ht_u(std::move(ibm_ht_u)),
m_ibm_ht_v(std::move(ibm_ht_v)),m_ibm_ht_w(std::move(ibm_ht_w)),
m_ibm_ht_c(std::move(ibm_ht_c)),m_ibm_z0(std::move(ibm_z0))
{}


// Fortran domain indices: (i,j) maps to multi_array index (i-m_ims,j-m_jms).
std::int32_t 
gms::radiolocation
::terrain_topography_avx512_t
 ::compute_idealized_terrain()
 {
    std::int32_t ret_status = -9999;
    auto U = [&](std::int32_t i, std::int32_t j) -> float& { return m_ibm_ht_u[i-m_ims][j-m_jms]; };
    auto V = [&](std::int32_t i, std::int32_t j) -> float& { return m_ibm_ht_v[i-m_ims][j-m_jms]; };
    auto W = [&](std::int32_t i, std::int32_t j) -> float& { return m_ibm_ht_w[i-m_ims][j-m_jms]; };
    auto C = [&](std::int32_t i, std::int32_t j) -> float& { return m_ibm_ht_c[i-m_ims][j-m_jms]; };
    std::int32_t icm=0, jcm=0, width=0, ista=0, iend=0, jsta=0, jend=0;
    float hm=0, xa=0, offset=0, slope=0, xa1=0, x_dist=0, y_dist=0;
    constexpr float pi = 3.1415926535897932384626433832795;

    if (m_ideal_terrain == 1) {
        for (std::int32_t j=m_jts;j<=m_jte;++j) for (std::int32_t i=m_its;i<=m_ite;++i) {
            if (j!=m_jde) U(i,j)=m_flat_plate_ht;
            if (i!=m_ide) V(i,j)=m_flat_plate_ht;
            if (i!=m_ide && j!=m_jde) W(i,j)=m_flat_plate_ht;
            C(i,j)=m_flat_plate_ht;
        }
    } else if (m_ideal_terrain == 111) {
        for (std::int32_t j=m_jts;j<=m_jte;++j) for (std::int32_t i=m_its;i<=m_ite;++i) {
            if(j!=m_jde) U(i,j)=0; if(i!=m_ide) V(i,j)=0;
            if(i!=m_ide&&j!=m_jde) W(i,j)=0; C(i,j)=0;
        }
        hm=6.0f; width=6;
        ista=(m_ide-m_ids)/2+m_ids-width; iend=(m_ide-m_ids)/2+m_ids+width;
        for(std::int32_t i=m_its;i<=m_ite;++i) if(i>=ista&&i<=iend)
            for(std::int32_t j=m_jts;j<=m_jte;++j) {
                if(j!=m_jde) U(i,j)=hm; if(i!=m_ide) V(i,j)=hm;
                if(i!=m_ide&&j!=m_jde) W(i,j)=hm; C(i,j)=hm;
            }
        jsta=(m_jde-m_jds)/2+m_jds-width; jend=(m_jde-m_jds)/2+m_jds+width;
        for(std::int32_t j=m_jts;j<=m_jte;++j) if(j>=jsta&&j<=jend)
            for(std::int32_t i=m_its;i<=m_ite;++i) {
                if(j!=m_jde) U(i,j)=hm; if(i!=m_ide) V(i,j)=hm;
                if(i!=m_ide&&j!=m_jde) W(i,j)=hm; C(i,j)=hm;
            }
        for(std::int32_t i=m_its;i<=m_ite;++i) for(std::int32_t j=m_jts;j<=m_jte;++j)
            if(i<m_ids+2*m_spec_bdy_width || i>=m_ide-2*m_spec_bdy_width ||
               j<m_jds+2*m_spec_bdy_width || j>=m_jde-2*m_spec_bdy_width) {
                if(i!=m_ide) U(i,j)=0; if(j!=m_jde) V(i,j)=0;
                if(i!=m_ide&&j!=m_jde) W(i,j)=0; C(i,j)=0;
            }
        ista=(m_ide-m_ids)/2+m_ids-width; iend=(m_ide-m_ids)/2+m_ids+width;
        for(std::int32_t i=m_its;i<=m_ite;++i) if(i>=ista&&i<=iend)
            for(std::int32_t j=m_jts;j<=m_jte;++j) {
                jsta=width*(j/width); jend=jsta+width;
                if(j!=m_jde&&i!=ista&&j!=jsta)
                    U(i,j)=hm*(1.0f-std::abs((m_ide-m_ids)/2.0f+m_ids-(i-0.5f))/width)*(1.0f-std::abs(jsta+width/2.0f-j)/(width/2.0f));
                if(i!=m_ide&&i!=ista&&j!=jsta)
                    V(i,j)=hm*(1.0f-std::abs((m_ide-m_ids)/2.0f+m_ids-i)/width)*(1.0f-std::abs(jsta+width/2.0f-(j-0.5f))/(width/2.0f));
                if(j!=m_jde&&i!=ista)
                    W(i,j)=hm*(1.0f-std::abs((m_ide-m_ids)/2.0f+m_ids-i)/width)*(1.0f-std::abs(jsta+width/2.0f-j)/(width/2.0f));
                C(i,j)=hm*(1.0f-std::abs((m_ide-m_ids)/2.0f+m_ids-(i-0.5f))/width)*(1.0f-std::abs(jsta+width/2.0f-(j-0.5f))/(width/2.0f));
            }
        for(std::int32_t i=m_its;i<=m_ite;++i) for(std::int32_t j=m_jts;j<=m_jte;++j)
            if(i<m_ids+m_spec_bdy_width || i>=m_ide-m_spec_bdy_width ||
               j<m_jds+m_spec_bdy_width || j>=m_jde-m_spec_bdy_width) {
                if(i!=m_ide) U(i,j)=0; if(j!=m_jde) V(i,j)=0;
                if(i!=m_ide&&j!=m_jde) W(i,j)=0; C(i,j)=0;
            }
    } else if (m_ideal_terrain == 112 || m_ideal_terrain == 113) {
        for(std::int32_t j=m_jts;j<=m_jte;++j) for(std::int32_t i=m_its;i<=m_ite;++i) {
            if(j!=m_jde) U(i,j)=0; if(i!=m_ide) V(i,j)=0;
            if(i!=m_ide&&j!=m_jde) W(i,j)=0; C(i,j)=0;
        }
        for(std::int32_t j=m_jts;j<=m_jte;++j) for(std::int32_t i=m_its;i<=m_ite;++i)
            if((i/10)%2==1&&(j/10)%2==1) {
                if(j!=m_jde) U(i,j)=25; if(i!=m_ide) V(i,j)=25;
                if(i!=m_ide&&j!=m_jde) W(i,j)=25; C(i,j)=25;
            }
        for(std::int32_t j=m_jts;j<=m_jte;++j) for(std::int32_t i=m_its;i<=m_ite;++i) {
            if(i>=10*(30/10)&&i<=m_ide-10*(30/10)&&j>=10*(34/10)&&j<=m_jde-10*(34/10)) {
                if(j!=m_jde) U(i,j)=0; if(i!=m_ide) V(i,j)=0;
                if(i!=m_ide&&j!=m_jde) W(i,j)=0; C(i,j)=0;
            }
            if(i<=5||i>=m_ide-6||j<=5||j>=m_jde-6) {
                if(j!=m_jde) U(i,j)=0; if(i!=m_ide) V(i,j)=0;
                if(i!=m_ide&&j!=m_jde) W(i,j)=0; C(i,j)=0;
            }
        }
    } else if (m_ideal_terrain==2 || m_ideal_terrain==3 ||
               m_ideal_terrain==33 || m_ideal_terrain==11) {
        hm=250.0f; xa=800.0f; offset=0;
        icm=m_ide/2; jcm=m_jde/2;
        if(m_ideal_terrain==2 || m_ideal_terrain==11) {
            for(std::int32_t j=m_jts;j<=m_jte;++j) for(std::int32_t i=m_its;i<=m_ite;++i) {
                if(j!=m_jde) U(i,j)=offset+hm/(1+std::pow((i-icm-0.5f)/xa,2));
                if(i!=m_ide) V(i,j)=offset+hm/(1+std::pow(float(i-icm)/xa,2));
                if(i!=m_ide&&j!=m_jde) W(i,j)=offset+hm/(1+std::pow(float(i-icm)/xa,2));
                C(i,j)=offset+hm/(1+std::pow((i-icm-0.5f)/xa,2));
                if(j!=m_jde) U(i,j)=offset+hm/(1+std::pow(float(j-jcm)/xa,2));
                if(i!=m_ide) V(i,j)=offset+hm/(1+std::pow((j-jcm-0.5f)/xa,2));
                W(i,j)=offset+hm/(1+std::pow(float(j-jcm)/xa,2));
                C(i,j)=offset+hm/(1+std::pow((j-jcm-0.5f)/xa,2));
                if(m_ideal_terrain==2) {
                    x_dist=(-m_ide/2.0f)*m_dx+i*m_dx;
                    if(i!=m_ide&&j!=m_jde) W(i,j)=offset+hm/(1+std::pow(x_dist/xa,2));
                    if(i!=m_ide) V(i,j)=offset+hm/(1+std::pow(x_dist/xa,2));
                    x_dist=(-m_ide/2.0f)*m_dx-m_dx/2+i*m_dx;
                    if(j!=m_jde) U(i,j)=offset+hm/(1+std::pow(x_dist/xa,2)); C(i,j)=offset+hm/(1+std::pow(x_dist/xa,2));
                    y_dist=(-m_jde/2.0f)*m_dy+j*m_dy;
                    if(i!=m_ide&&j!=m_jde) W(i,j)=offset+hm/(1+std::pow(y_dist/xa,2));
                    if(j!=m_jde) U(i,j)=offset+hm/(1+std::pow(y_dist/xa,2));
                    y_dist=(-m_jde/2.0f)*m_dy-m_dy/2+j*m_dy;
                    if(i!=m_ide) V(i,j)=offset+hm/(1+std::pow(y_dist/xa,2)); C(i,j)=offset+hm/(1+std::pow(y_dist/xa,2));
                } else {
                    x_dist=(-m_jde/2.0f)*m_dy+j*m_dy;
                    if(i!=m_ide&&j!=m_jde) W(i,j)=offset+hm/(1+std::pow(x_dist/xa,2));
                    if(j!=m_jde) U(i,j)=offset+hm/(1+std::pow(x_dist/xa,2));
                    x_dist=(-m_jde/2.0f)*m_dy-m_dy/2+j*m_dy;
                    if(i!=m_ide) V(i,j)=offset+hm/(1+std::pow(x_dist/xa,2)); C(i,j)=offset+hm/(1+std::pow(x_dist/xa,2));
                }
            }
        } else {
            for(std::int32_t j=m_jts;j<=m_jte;++j) for(std::int32_t i=m_its;i<=m_ite;++i) {
                if(j!=m_jde) U(i,j)=hm/(1+std::pow((i-icm-0.5f)/xa,2)+std::pow(float(j-jcm)/xa,2));
                if(i!=m_ide) V(i,j)=hm/(1+std::pow(float(i-icm)/xa,2)+std::pow((j-jcm-0.5f)/xa,2));
                if(i!=m_ide&&j!=m_jde) W(i,j)=hm/(1+std::pow(float(i-icm)/xa,2)+std::pow(float(j-jcm)/xa,2));
                C(i,j)=hm/(1+std::pow((i-icm-0.5f)/xa,2)+std::pow((j-jcm-0.5f)/xa,2));
                x_dist=(-m_ide/2.0f)*m_dx+i*m_dx;
                y_dist=(-m_jde/2.0f)*m_dy+j*m_dy;
                if(i!=m_ide&&j!=m_jde) W(i,j)=offset+hm/(1+std::pow(x_dist/xa,2)+std::pow(y_dist/xa,2));
                y_dist=(-m_jde/2.0f)*m_dy-m_dy/2+j*m_dy;
                if(i!=m_ide) V(i,j)=offset+hm/(1+std::pow(x_dist/xa,2)+std::pow(y_dist/xa,2));
                x_dist=(-m_ide/2.0f)*m_dx-m_dx/2+i*m_dx;
                y_dist=(-m_jde/2.0f)*m_dy+j*m_dy;
                if(j!=m_jde) U(i,j)=offset+hm/(1+std::pow(x_dist/xa,2)+std::pow(y_dist/xa,2));
                y_dist=(-m_jde/2.0f)*m_dy-m_dy/2+j*m_dy;
                C(i,j)=offset+hm/(1+std::pow(x_dist/xa,2)+std::pow(y_dist/xa,2));
            }
        }
    } else if(m_ideal_terrain==21) {
        xa=150; hm=100; offset=0;
        for(std::int32_t j=m_jts;j<=m_jte;++j) for(std::int32_t i=m_its;i<=m_ite;++i) {
            x_dist=(-m_ide/2.0f)*m_dx+i*m_dx;
            if(i!=m_ide&&j!=m_jde) W(i,j)=offset+hm*std::exp(-0.5f*std::pow(x_dist/xa,2));
            if(i!=m_ide) V(i,j)=offset+hm*std::exp(-0.5f*std::pow(x_dist/xa,2));
            x_dist=(-m_ide/2.0f)*m_dx-m_dx/2+i*m_dx;
            if(j!=m_jde) U(i,j)=offset+hm*std::exp(-0.5f*std::pow(x_dist/xa,2));
            C(i,j)=offset+hm*std::exp(-0.5f*std::pow(x_dist/xa,2));
        }
    } else if(m_ideal_terrain==31) {
        xa=150; hm=100; offset=0;
        for(std::int32_t j=m_jts;j<=m_jte;++j) for(std::int32_t i=m_its;i<=m_ite;++i) {
            x_dist=(-m_ide/2.0f)*m_dx+i*m_dx; y_dist=(-m_jde/2.0f)*m_dy+j*m_dy;
            if(i!=m_ide&&j!=m_jde) W(i,j)=hm*std::exp(-0.5f*(std::pow(x_dist/xa,2)+std::pow(y_dist/xa,2)));
            y_dist=(-m_jde/2.0f)*m_dy-m_dy/2+j*m_dy;
            if(i!=m_ide) V(i,j)=hm*std::exp(-0.5f*(std::pow(x_dist/xa,2)+std::pow(y_dist/xa,2)));
            x_dist=(-m_ide/2.0f)*m_dx-m_dx/2+i*m_dx; y_dist=(-m_jde/2.0f)*m_dy+j*m_dy;
            if(j!=m_jde) U(i,j)=hm*std::exp(-0.5f*(std::pow(x_dist/xa,2)+std::pow(y_dist/xa,2)));
            y_dist=(-m_jde/2.0f)*m_dy-m_dy/2+j*m_dy;
            C(i,j)=hm*std::exp(-0.5f*(std::pow(x_dist/xa,2)+std::pow(y_dist/xa,2)));
        }
    } else if(m_ideal_terrain==4) {
        xa=150; hm=50; offset=0;
        for(std::int32_t j=m_jts;j<=m_jte;++j) for(std::int32_t i=m_its;i<=m_ite;++i) {
            if(j!=m_jde) U(i,j)=offset; if(i!=m_ide) V(i,j)=offset;
            if(i!=m_ide&&j!=m_jde) W(i,j)=offset; C(i,j)=offset;
            icm=m_jde/4; x_dist=-icm*m_dy+j*m_dy;
            if(std::abs(x_dist)<=xa) { if(i!=m_ide&&j!=m_jde) W(i,j)=offset+hm; if(i!=m_ide) V(i,j)=offset+hm; }
            x_dist=-icm*m_dy-m_dy/2+j*m_dy;
            if(std::abs(x_dist)<=xa) { if(j!=m_jde) U(i,j)=offset+hm; C(i,j)=offset+hm; }
        }
    } else if(m_ideal_terrain==5) {
        xa=48; hm=60; offset=50; icm=m_ide/2; jcm=m_jde/2;
        for(std::int32_t j=m_jts;j<=m_jte;++j) for(std::int32_t i=m_its;i<=m_ite;++i) {
            x_dist=-icm*m_dx+i*m_dx; y_dist=-jcm*m_dy+j*m_dy;
            if(std::abs(x_dist)<=xa&&std::abs(y_dist)<=xa) { if(i!=m_ide&&j!=m_jde) W(i,j)=offset+hm; } else if(i!=m_ide&&j!=m_jde) W(i,j)=offset;
            y_dist=-jcm*m_dy-m_dy/2+j*m_dy;
            if(std::abs(x_dist)<=xa&&std::abs(y_dist)<=xa) { if(i!=m_ide) V(i,j)=offset+hm; } else if(i!=m_ide) V(i,j)=offset;
            x_dist=-icm*m_dx-m_dx/2+i*m_dx; y_dist=-jcm*m_dy+j*m_dy;
            if(std::abs(x_dist)<=xa&&std::abs(y_dist)<=xa) { if(j!=m_jde) U(i,j)=offset+hm; } else if(j!=m_jde) U(i,j)=offset;
            y_dist=-jcm*m_dy-m_dy/2+j*m_dy;
            if(std::abs(x_dist)<=xa&&std::abs(y_dist)<=xa) C(i,j)=offset+hm; else C(i,j)=offset;
        }
    } else if(m_ideal_terrain==6) {
        hm=1500;
        for(std::int32_t j=m_jts;j<=m_jte;++j) for(std::int32_t i=m_its;i<=m_ite;++i) {
            x_dist=(-m_ide/2.0f)*m_dx+i*m_dx;
            auto h=[&](float x) { float a=std::abs(x); if(a<=500) return 0.0f; if(a<9500) return hm*(0.5f-0.5f*std::cos(pi*(a-500)/9000)); if(a<=10500) return hm; if(a<19500) return hm*(0.5f+0.5f*std::cos(pi*(a-10500)/9000)); return 0.0f; };
            if(i!=m_ide&&j!=m_jde) W(i,j)=h(x_dist); if(i!=m_ide) V(i,j)=h(x_dist);
            x_dist=(-m_ide/2.0f)*m_dx-m_dx/2+i*m_dx;
            if(j!=m_jde) U(i,j)=h(x_dist); C(i,j)=h(x_dist);
        }
    } else if(m_ideal_terrain==7) {
        // The Fortran branch is std::int32_tentionally empty.
        std::cerr << "RF-EMT[compute_idealized_terrain]: ***Exit*** -- reached invalid parameter: m_ideal_terrain=" << m_ideal_terrain << '\n';
        ret_status = -1;
        return (ret_status);
    } else if(m_ideal_terrain==8 || m_ideal_terrain==81) {
        const bool sinusoidal=m_ideal_terrain==81;
        offset=sinusoidal?550.0f:0.0f; slope=0.2126f; xa=1000.0f;
        for(std::int32_t j=m_jts;j<=m_jte;++j) for(std::int32_t i=m_its;i<=m_ite;++i) {
            x_dist=(-m_ide/2)*m_dx-m_dx/2+i*m_dx;
            float hv;
            if(std::abs(x_dist)<(m_ide/2)*m_dx-xa)
                hv=sinusoidal ? offset/2-offset/2*std::cos(2*pi*std::abs(x_dist)/((m_ide-1)*m_dx-2*xa)) : offset+slope*std::abs(x_dist);
            else hv=sinusoidal?offset:offset+slope*((m_ide/2)*m_dx-xa);
            if(i!=m_ide&&j!=m_jde) W(i,j)=hv; if(i!=m_ide) V(i,j)=hv;
            x_dist=(-m_ide/2)*m_dx+(i-1)*m_dx;
            if(std::abs(x_dist)<(m_ide/2)*m_dx-xa)
                hv=sinusoidal ? offset/2-offset/2*std::cos(2*pi*std::abs(x_dist)/((m_ide-1)*m_dx-2*xa)) : offset+slope*std::abs(x_dist);
            else hv=sinusoidal?offset:offset+slope*((m_ide/2)*m_dx-xa);
            if(j!=m_jde) U(i,j)=hv; C(i,j)=hv;
        }
    } else if(m_ideal_terrain==9) {
        hm=3000; xa=25000; xa1=8000;
        for(std::int32_t j=m_jts;j<=m_jte;++j) for(std::int32_t i=m_its;i<=m_ite;++i) {
            x_dist=(-m_ide/2.0f)*m_dx+i*m_dx;
            float h=std::abs(x_dist)<=xa ? hm*std::pow(std::cos(pi*x_dist/(2*xa)),2)*std::pow(std::cos(pi*x_dist/xa1),2) : 0;
            if(i!=m_ide&&j!=m_jde) W(i,j)=h; if(i!=m_ide) V(i,j)=h;
            x_dist=(-m_ide/2.0f)*m_dx-m_dx/2+i*m_dx;
            h=std::abs(x_dist)<=xa ? hm*std::pow(std::cos(pi*x_dist/(2*xa)),2)*std::pow(std::cos(pi*x_dist/xa1),2) : 0;
            if(j!=m_jde) U(i,j)=h; C(i,j)=h;
        }
    } else if(m_ideal_terrain==10) {
        for(std::int32_t j=m_jts;j<=m_jte;++j) for(std::int32_t i=m_its;i<=m_ite;++i) {
            if(j!=m_jde) U(i,j)=offset+U(i,j); if(i!=m_ide) V(i,j)=offset+V(i,j);
            if(i!=m_ide&&j!=m_jde) W(i,j)=offset+W(i,j); C(i,j)=offset+C(i,j);
        }
    } else if(m_ideal_terrain==12) {
        xa=51; hm=60; offset=50; icm=m_ide/2;
        for(std::int32_t j=m_jts;j<=m_jte;++j) for(std::int32_t i=m_its;i<=m_ite;++i) {
            x_dist=-icm*m_dx+i*m_dx;
            if(std::abs(x_dist)<=xa) { if(i!=m_ide&&j!=m_jde) W(i,j)=offset+hm; if(i!=m_ide) V(i,j)=offset+hm; }
            else { if(i!=m_ide&&j!=m_jde) W(i,j)=offset; if(i!=m_ide) V(i,j)=offset; }
            x_dist=-icm*m_dx-m_dx/2+i*m_dx;
            if(std::abs(x_dist)<=xa) { if(j!=m_jde) U(i,j)=offset+hm; C(i,j)=offset+hm; }
            else { if(j!=m_jde) U(i,j)=offset; C(i,j)=offset; }
        }
    } else if(m_ideal_terrain==13 || m_ideal_terrain==14 || m_ideal_terrain==15) {
        std::cerr << "RF-EMT[compute_idealized_terrain]: ***Exit*** -- reached invalid parameter: m_ideal_terrain=" << m_ideal_terrain << '\n';
        ret_status = -2;
        return (ret_status);
    }
    ret_status = 0;
    return (ret_status);
}

std::int32_t 
gms::radiolocation
::terrain_topography_avx512_t
::modify_idealized_terrain()
{
    std::int32_t ret_status = -9999;
    ret_status = compute_idealized_terrain();
    if(ret_status!=0)
    {
      std::cerr << "RF-EMT[modify_idealized_terrain]: ***Exit*** -- received error state provided by the callee: ret_status=" << ret_status << '\n';
      return (ret_status);
    }
    // Map the original Fortran indices to zero-based Boost.MultiArray indices.
    auto U = [&](std::int32_t i, std::int32_t j) -> float& {
        return m_ibm_ht_u[i - m_ims][j - m_jms];
    };
    auto V = [&](std::int32_t i, std::int32_t j) -> float& {
        return m_ibm_ht_v[i - m_ims][j - m_jms];
    };
    auto W = [&](std::int32_t i, std::int32_t j) -> float& {
        return m_ibm_ht_w[i - m_ims][j - m_jms];
    };
    auto C = [&](std::int32_t i, std::int32_t j) -> float& {
        return m_ibm_ht_c[i - m_ims][j - m_jms];
    };
    auto PHB = [&](std::int32_t i, std::int32_t k, std::int32_t j) -> float {
        return m_phb[i - m_ims][k - m_kms][j - m_jms];
    };
    constexpr float gravity = 9.81f;
    // Terrain following the lowest geopotential level.
    auto set_terrain_from_phb = [&]() {
        for (int i = m_its; i <= m_ite; ++i) {
            for (int j = m_jts; j <= m_jte; ++j) {
                W(i, j) = PHB(i, 1, j) / gravity;
                U(i, j) = (PHB(i - 1, 1, j) + PHB(i, 1, j))
                          / (2.0f * gravity);
                V(i, j) = (PHB(i, 1, j - 1) + PHB(i, 1, j))
                          / (2.0f * gravity);
                C(i, j) = (PHB(i - 1, 1, j - 1) + PHB(i - 1, 1, j)
                         + PHB(i, 1, j - 1) + PHB(i, 1, j))
                          / (4.0f * gravity);
            }
        }
    };

    if (m_ideal_terrain == 0) {
        set_terrain_from_phb();
    } else if (m_ideal_terrain == 114) {
        set_terrain_from_phb();

        const int ic = (m_ide - m_ids) / 2 + m_ids;
        const int jc = (m_jde - m_jds) / 2 + m_jds;

        if (((m_its >= ic - 5 || m_ite <= ic + 5) &&
             (m_jts >= jc - 5 || m_jte <= jc + 5))) {
            const int ilo = (m_its > ic - 5) ? m_its : ic - 5;
            const int ihi = (m_ite < ic + 5) ? m_ite : ic + 5;
            const int jlo = (m_jts > jc - 5) ? m_jts : jc - 5;
            const int jhi = (m_jte < jc + 5) ? m_jte : jc + 5;

            for (int i = ilo; i <= ihi; ++i) {
                for (int j = jlo; j <= jhi; ++j) {
                    W(i, j) = 400.0f;
                    U(i, j) = 400.0f;
                    V(i, j) = 400.0f;
                    C(i, j) = 400.0f;
                }
            }
        }
    } else if (m_ideal_terrain == 115) {
        set_terrain_from_phb();

        const float bldx = 75.0f;
        const float bldy = 75.0f;
        const float bldz = 10.0f;
        const float gapx = 150.0f;
        const float gapy = 150.0f;
        const float theta = 10.0f * 2.0f * 4.0f
                          * static_cast<float>(std::atan(1.0)) / 360.0f;

        for (int i = (m_its - 1) * 2; i <= (m_ite - 1) * 2; ++i) {
            for (int j = (m_jts - 1) * 2; j <= (m_jte - 1) * 2; ++j) {
                const float xoff = (j * m_dy / 2.0f)
                                 * std::tan(theta);
                const float bx = std::fmod(
                    i * m_dx / 2.0f + xoff, gapx + bldx);

                if (bx >= (gapx + bldx) / 2.0f - bldx / 2.0f &&
                    bx <= (gapx + bldx) / 2.0f + bldx / 2.0f) {
                    const float yoff = (i * m_dx / 2.0f)
                                     * std::tan(-theta);
                    const float by = std::fmod(
                        j * m_dy / 2.0f + yoff, gapy + bldy);

                    if (by >= (gapy + bldy) / 2.0f - bldy / 2.0f &&
                        by <= (gapy + bldy) / 2.0f + bldy / 2.0f) {
                        const int fi = static_cast<int>(std::floor(i / 2.0f)) + 1;
                        const int fj = static_cast<int>(std::floor(j / 2.0f)) + 1;
                        const bool i_even = std::floor(i / 2.0f) == i / 2.0f;
                        const bool j_even = std::floor(j / 2.0f) == j / 2.0f;

                        if (i_even && j_even) {
                            C(fi, fj) += bldz;
                        } else if (i_even) {
                            U(fi, fj) += bldz;
                        } else if (j_even) {
                            V(fi, fj) += bldz;
                        } else {
                            W(fi, fj) += bldz;
                        }
                    }
                }
            }
        }
    } else if (m_ideal_terrain == 116) {
        set_terrain_from_phb();

        const float bldx = 80.0f;
        const float bldy = 80.0f;
        const float bldz = 10.0f;
        const float gapx = 240.0f;
        const float gapy = 240.0f;
        const float theta = 30.0f; // Source expression uses degrees converted below.
        const float theta_rad = theta * 2.0f * static_cast<float>(std::atan(1.0))
                              / 180.0f;

        for (int j = (m_jts - 1) * 2; j <= (m_jte - 1) * 2 + 1; ++j) {
            const float row = j * m_dy / 2.0f;
            const float row_count = std::floor(row / (gapy + bldy));
            const float by = row - (row_count + 0.5f) * (gapy + bldy);
            const float ymid = (row_count + 0.5f) * (gapy + bldy);

            bool keepgoing = true;
            if (ymid - bldy / 2.0f <=
                    (m_spec_bdy_width - 0.5f) * m_dy ||
                ymid + bldy / 2.0f >=
                    (m_jde - m_spec_bdy_width - 0.5f) * m_dy) {
                keepgoing = false;
            }

            if (by >= -bldy / 2.0f && by <= bldy / 2.0f && keepgoing) {
                const float yoff = row_count * (gapy + bldy);
                const float xoff = yoff * std::tan(-theta_rad);

                for (int i = (m_its - 1) * 2; i <= (m_ite - 1) * 2 + 1; ++i) {
                    const float x = i * m_dx / 2.0f + xoff;
                    const float bx = x - (std::floor(x / (gapx + bldx)) + 0.5f)
                                         * (gapx + bldx);
                    const float xmid = i * m_dx / 2.0f - bx;

                    keepgoing = true;
                    if (xmid - bldx / 2.0f <=
                            (m_spec_bdy_width - 0.5f) * m_dx ||
                        xmid + bldx / 2.0f >=
                            (m_ide - m_spec_bdy_width - 0.5f) * m_dx) {
                        keepgoing = false;
                    }

                    if (bx >= -bldx / 2.0f && bx <= bldx / 2.0f && keepgoing) {
                        const int fi = static_cast<int>(std::floor(i / 2.0f)) + 1;
                        const int fj = static_cast<int>(std::floor(j / 2.0f)) + 1;
                        const bool i_even = std::floor(i / 2.0f) == i / 2.0f;
                        const bool j_even = std::floor(j / 2.0f) == j / 2.0f;

                        if (i_even && j_even) {
                            C(fi, fj) += bldz;
                        } else if (i_even) {
                            U(fi, fj) += bldz;
                        } else if (j_even) {
                            V(fi, fj) += bldz;
                        } else {
                            W(fi, fj) += bldz;
                        }
                    }
                }
            }
        }
    }
    ret_status = 0;
    return (ret_status);
}





