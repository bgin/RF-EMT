
#ifndef __GMS_IDEALIZED_TOPO_H__
#define __GMS_IDEALIZED_TOPO_H__ 160120220934

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