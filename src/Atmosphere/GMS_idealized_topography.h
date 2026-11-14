
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

#ifndef __GMS_IDEALIZED_TOPO_H__
#define __GMS_IDEALIZED_TOPO_H__ 160120220934

namespace file_info {

const unsigned int GMS_IDEALIZED_TOPO_MAJOR = 1;
const unsigned int GMS_IDEALIZED_TOPO_MINOR = 1;
const unsigned int GMS_IDEALIZED_TOPO_MICRO = 0;
const unsigned int GMS_IDEALIZED_TOPO_FULLVER =
    1000U * GMS_IDEALIZED_TOPO_MAJOR + 100U * GMS_IDEALIZED_TOPO_MINOR +
    10U * GMS_IDEALIZED_TOPO_MICRO;
const char *const GMS_IDEALIZED_TOPO_CREATION_DATE =
    "16-01-2022 09:34 +00200 (SAT 16 APR 2022 GMT+2)";
const char *const GMS_IDEALIZED_TOPO_BUILD_DATE = __DATE__ " " __TIME__;
const char *const GMS_IDEALIZED_TOPO_SYNOPSIS =
    "Idealized gaussian and sinusoidal topography simulator."

} // namespace file_info

#include <cstdint>
#include "GMS_config.h"

namespace gms {

namespace radiolocation {

__ATTR_HOT__
__ATTR_ALIGN__(32)
void gaussian_topography_r4(float *__restrict __ATTR_ALIGN__(64) lon,
                            const int32_t nlon,
                            float *__restrict __ATTR_ALIGN__(64) lat,
                            const int32_t nlat,
                            float *__restrict __ATTR_ALIGN__(64) zsurf,
                            const float height);

__ATTR_HOT__
__ATTR_ALIGN__(32)
void gaussian_topography_r8(double *__restrict __ATTR_ALIGN__(64) lon,
                            const int32_t nlon,
                            double *__restrict __ATTR_ALIGN__(64) lat,
                            const int32_t nlat,
                            double *__restrict __ATTR_ALIGN__(64) zsurf,
                            const double height);

__ATTR_HOT__
__ATTR_ALIGN__(32)
void sinusoidal_topography(float *__restrict __ATTR_ALIGN__(64) lon,
                           const int32_t nlon,
                           float *__restrict __ATTR_ALIGN__(64) lat,
                           const int32_t nlat,
                           float *__restrict __ATTR_ALIGN__(64) zsurf,
                           const float height_sin, const float m,
                           const float Amp2, const float uneven_fac,
                           const float deltalat, const bool uneven_sin);

__ATTR_HOT__
__ATTR_ALIGN__(32)
void sinusoidal_topography(double *__restrict __ATTR_ALIGN__(64) lon,
                           const int32_t nlon,
                           double *__restrict __ATTR_ALIGN__(64) lat,
                           const int32_t nlat,
                           double *__restrict __ATTR_ALIGN__(64) zsurf,
                           const double height_sin, const double m,
                           const double Amp2, const double uneven_fac,
                           const double deltalat, const bool uneven_sin);



} // namespace radiolocation

} // namespace gms

#endif /*__GMS_IDEALIZED_TOPO_H__*/