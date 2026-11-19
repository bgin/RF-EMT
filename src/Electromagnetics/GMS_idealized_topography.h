
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

#ifndef __GMS_IDEALIZED_TOPOGRAPHY_H__
#define __GMS_IDEALIZED_TOPOGRAPHY_H__ 160120220934

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
    "Idealized gaussian, sinusoidal and other analytically derived topography simulator.";

} // namespace file_info

#include <cstdint>
#include "GMS_config.h"
#include "GMS_boost_multiarray_types.h"
#include "GMS_avx512_memset.h"
#include "GMS_avx512_memcpy.h"


#if !defined(TERRAIN_TOPOGRAPHY_AVX512_VERBOSE)
#define TERRAIN_TOPOGRAPHY_AVX512_VERBOSE 1
#endif 

#if (TERRAIN_TOPOGRAPHY_AVX512_VERBOSE) == 1
#include <cstdio>
#endif 

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



struct alignas(64) terrain_topography_args_t
{
    std::size_t  ids;
    std::size_t  ide;
    std::size_t  jds;
    std::size_t  jde;
    std::size_t  kds;
    std::size_t  kde;
    std::size_t  ims;
    std::size_t  ime;
    std::size_t  jms;
    std::size_t  jme;
    std::size_t  kms;
    std::size_t  kme;
    std::size_t  its;
    std::size_t  ite;
    std::size_t  jts;
    std::size_t  jte;
    std::size_t  kts;
    std::size_t  kte;
    std::int32_t ideal_terrain;
    std::int32_t spec_bdy_width;
    float        flat_plate_ht;
    float        dx;
    float        dy; 
};

struct alignas(64) terrain_topography_avx512_t 
{
    enum class terrain_types_enum : std::int32_t 
    {
        ideal_terrain               = 0,
        X_m_flat_plate              = 1,
        witch_of_agnesi2D           = 2,
        gaussian_hill2D             = 21,
        witch_of_agnesi3D           = 3,
        gaussian_hill3D             = 31,
        ridge_or_isolated_cube      = 4,
        cube3D                      = 5,
        idealized_trex_valley2D     = 6,
        idealized_trex_valley3D     = 7,
        v_shaped_valley_ridges      = 8,
        sinusoidal_valley_ridges    = 81,
        schar_ideal_advect_test     = 9,
        read_data_1D_ht_array       = 10,
        witch_of_agnesi2D_ydim      = 11,
        djw_okc_read_form_file      = 13,
        rsa_granite_mount_file_read = 14,
        rsa_jb_askervein_file_read  = 15 
    };

    std::size_t  m_ids;
    std::size_t  m_ide;
    std::size_t  m_jds;
    std::size_t  m_jde;
    std::size_t  m_kds;
    std::size_t  m_kde;
    std::size_t  m_ims;
    std::size_t  m_ime;
    std::size_t  m_jms;
    std::size_t  m_jme;
    std::size_t  m_kms;
    std::size_t  m_kme;
    std::size_t  m_its;
    std::size_t  m_ite;
    std::size_t  m_jts;
    std::size_t  m_jte;
    std::size_t  m_kts;
    std::size_t  m_kte;
    std::int32_t m_ideal_terrain;
    std::int32_t m_spec_bdy_width;
    float        m_flat_plate_ht;
    float        m_dx;
    float        m_dy;
    
    common::Aligned64Mat3DF32 m_phb;
    common::Aligned64Mat3DF32 m_ph;
    common::Aligned64Mat2DF32 m_ibm_ht_u;
    common::Aligned64Mat2DF32 m_ibm_ht_v;
    common::Aligned64Mat2DF32 m_ibm_ht_w;
    common::Aligned64Mat2DF32 m_ibm_ht_c;
    common::Aligned64Mat2DF32 m_ibm_z0;

    terrain_topography_avx512_t() noexcept(true) = default;

    terrain_topography_avx512_t(const terrain_topography_args_t &) noexcept(false);

    terrain_topography_avx512_t(const terrain_topography_args_t &,
                                const common::Aligned64Mat3DF32 &,
                                const common::Aligned64Mat3DF32 &,
                                const common::Aligned64Mat2DF32 &,
                                const common::Aligned64Mat2DF32 &,
                                const common::Aligned64Mat2DF32 &,
                                const common::Aligned64Mat2DF32 &,
                                const common::Aligned64Mat2DF32 &) noexcept(false);

    terrain_topography_avx512_t(const terrain_topography_args_t &,
                                common::Aligned64Mat3DF32 &&,
                                common::Aligned64Mat3DF32 &&,
                                common::Aligned64Mat2DF32 &&,
                                common::Aligned64Mat2DF32 &&,
                                common::Aligned64Mat2DF32 &&,
                                common::Aligned64Mat2DF32 &&,
                                common::Aligned64Mat2DF32 &&) noexcept(true);

    terrain_topography_avx512_t(const terrain_topography_avx512_t &) noexcept(false) = default;

    terrain_topography_avx512_t( terrain_topography_avx512_t &&) noexcept(true) = default;

    ~terrain_topography_avx512_t() noexcept(true) = default;

    terrain_topography_avx512_t &
                       operator=(const terrain_topography_avx512_t &) noexcept(false) = default;

    terrain_topography_avx512_t &
                       operator=(terrain_topography_avx512_t &&) noexcept(true) = default; 

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_phb(const float filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        float * __restrict__ ptr = &this->m_phb.data()[0];
        if constexpr(which_unrolling==16)
        {
           common::avx512_memset_unroll16x_ps(&ptr[0],filler,this->m_phb.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
           common::avx512_memset_unroll8x_ps(&ptr[0],filler,this->m_phb.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_ph(const float filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        float * __restrict__ ptr = &this->m_ph.data()[0];
        if constexpr(which_unrolling==16)
        {
           common::avx512_memset_unroll16x_ps(&ptr[0],filler,this->m_ph.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
           common::avx512_memset_unroll8x_ps(&ptr[0],filler,this->m_ph.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_ibm_ht_u(const float filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        float * __restrict__ ptr = &this->m_ibm_ht_u.data()[0];
        if constexpr(which_unrolling==16)
        {
           common::avx512_memset_unroll16x_ps(&ptr[0],filler,this->m_ibm_ht_u.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
           common::avx512_memset_unroll8x_ps(&ptr[0],filler,this->m_ibm_ht_u.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_ibm_ht_v(const float filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        float * __restrict__ ptr = &this->m_ibm_ht_v.data()[0];
        if constexpr(which_unrolling==16)
        {
           common::avx512_memset_unroll16x_ps(&ptr[0],filler,this->m_ibm_ht_v.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
           common::avx512_memset_unroll8x_ps(&ptr[0],filler,this->m_ibm_ht_v.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_ibm_ht_w(const float filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        float * __restrict__ ptr = &this->m_ibm_ht_w.data()[0];
        if constexpr(which_unrolling==16)
        {
           common::avx512_memset_unroll16x_ps(&ptr[0],filler,this->m_ibm_ht_w.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
           common::avx512_memset_unroll8x_ps(&ptr[0],filler,this->m_ibm_ht_w.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_ibm_ht_c(const float filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        float * __restrict__ ptr = &this->m_ibm_ht_c.data()[0];
        if constexpr(which_unrolling==16)
        {
           common::avx512_memset_unroll16x_ps(&ptr[0],filler,this->m_ibm_ht_c.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
           common::avx512_memset_unroll8x_ps(&ptr[0],filler,this->m_ibm_ht_c.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_ibm_z0(const float filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        float * __restrict__ ptr = &this->m_ibm_z0.data()[0];
        if constexpr(which_unrolling==16)
        {
           common::avx512_memset_unroll16x_ps(&ptr[0],filler,this->m_ibm_z0.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
           common::avx512_memset_unroll8x_ps(&ptr[0],filler,this->m_ibm_z0.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    set_phb(const common::Aligned64Mat3DF32 &other_phb) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        if(__builtin_expect(this->m_phb.num_elements()!=other_phb.num_elements(),0)) { ret_status = -1; return (ret_status);}
        double * __restrict__ p_m_phb     = &this->m_phb.data()[0];
        double * __restrict__ p_other_phb = &other_phb.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memcpy_unroll16x_ps(&p_m_phb[0],&other_phb[0],this->m_phb.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memcpy_unroll8x_ps(&p_m_phb[0],&other_phb[0],this->m_phb.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    set_ph(const common::Aligned64Mat3DF32 &other_ph) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        if(__builtin_expect(this->m_ph.num_elements()!=other_ph.num_elements(),0)) { ret_status = -1; return (ret_status);}
        double * __restrict__ p_m_ph     = &this->m_ph.data()[0];
        double * __restrict__ p_other_ph = &other_ph.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memcpy_unroll16x_ps(&p_m_ph[0],&other_ph[0],this->m_ph.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memcpy_unroll8x_ps(&p_m_ph[0],&other_ph[0],this->m_ph.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    set_ibm_ht_u(const common::Aligned64Mat3DF32 &other_ibm_ht_u) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        if(__builtin_expect(this->m_ibm_ht_u.num_elements()!=other_ibm_ht_u.num_elements(),0)) { ret_status = -1; return (ret_status);}
        double * __restrict__ p_m_ibm_ht_u     = &this->m_ibm_ht_u.data()[0];
        double * __restrict__ p_other_ibm_ht_u = &other_ibm_ht_u.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memcpy_unroll16x_ps(&p_m_ibm_ht_u[0],&other_ibm_ht_u[0],this->m_ibm_ht_u.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memcpy_unroll8x_ps(&p_m_ibm_ht_u[0],&other_ibm_ht_u[0],this->m_ibm_ht_u.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    set_ibm_ht_v(const common::Aligned64Mat3DF32 &other_ibm_ht_v) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        if(__builtin_expect(this->m_ibm_ht_v.num_elements()!=other_ibm_ht_v.num_elements(),0)) { ret_status = -1; return (ret_status);}
        double * __restrict__ p_m_ibm_ht_v     = &this->m_ibm_ht_v.data()[0];
        double * __restrict__ p_other_ibm_ht_v = &other_ibm_ht_v.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memcpy_unroll16x_ps(&p_m_ibm_ht_v[0],&other_ibm_ht_v[0],this->m_ibm_ht_v.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memcpy_unroll8x_ps(&p_m_ibm_ht_v[0],&other_ibm_ht_v[0],this->m_ibm_ht_v.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    set_ibm_ht_w(const common::Aligned64Mat3DF32 &other_ibm_ht_w) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        if(__builtin_expect(this->m_ibm_ht_w.num_elements()!=other_ibm_ht_w.num_elements(),0)) { ret_status = -1; return (ret_status);}
        double * __restrict__ p_m_ibm_ht_w     = &this->m_ibm_ht_w.data()[0];
        double * __restrict__ p_other_ibm_ht_w = &other_ibm_ht_w.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memcpy_unroll16x_ps(&p_m_ibm_ht_w[0],&other_ibm_ht_w[0],this->m_ibm_ht_w.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memcpy_unroll8x_ps(&p_m_ibm_ht_w[0],&other_ibm_ht_w[0],this->m_ibm_ht_w.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    set_ibm_ht_c(const common::Aligned64Mat3DF32 &other_ibm_ht_c) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        if(__builtin_expect(this->m_ibm_ht_c.num_elements()!=other_ibm_ht_c.num_elements(),0)) { ret_status = -1; return (ret_status);}
        double * __restrict__ p_m_ibm_ht_c     = &this->m_ibm_ht_c.data()[0];
        double * __restrict__ p_other_ibm_ht_c = &other_ibm_ht_c.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memcpy_unroll16x_ps(&p_m_ibm_ht_c[0],&other_ibm_ht_c[0],this->m_ibm_ht_c.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memcpy_unroll8x_ps(&p_m_ibm_ht_c[0],&other_ibm_ht_c[0],this->m_ibm_ht_c.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    set_ibm_z0(const common::Aligned64Mat3DF32 &other_ibm_z0) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        if(__builtin_expect(this->m_ibm_z0.num_elements()!=other_ibm_z0.num_elements(),0)) { ret_status = -1; return (ret_status);}
        double * __restrict__ p_m_ibm_z0     = &this->m_ibm_z0.data()[0];
        double * __restrict__ p_other_ibm_z0 = &other_ibm_z0.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memcpy_unroll16x_ps(&p_m_ibm_z0[0],&other_ibm_z0[0],this->m_ibm_z0.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memcpy_unroll8x_ps(&p_m_ibm_z0[0],&other_ibm_z0[0],this->m_ibm_z0.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    std::int32_t compute_idealized_terrain();
    
    __ATTR_HOT__
    __ATTR_ALIGN__(32)
    std::int32_t modify_idealied_terrain();

};

} // namespace radiolocation

} // namespace gms

#endif /*__GMS_IDEALIZED_TOPOGRAPHY_H__*/