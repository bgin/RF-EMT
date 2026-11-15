/*
WRF was developed at the National Center for Atmospheric Research (NCAR) which is operated by the University Corporation for Atmospheric Research (UCAR). 
NCAR and UCAR make no proprietary claims, either statutory or otherwise, to this version and release of WRF and consider WRF to be in the public domain 
for use by any person or entity for any purpose without any fee or charge. 
UCAR requests that any WRF user include this notice on any partial or full copies of WRF. 
WRF is provided on an "AS IS" basis and any warranties, either express or implied, including but not limited to implied warranties of non-infringement, 
originality, merchantability and fitness for a particular purpose, are disclaimed. 
In no event shall UCAR be liable for any damages, whatsoever, whether direct, indirect, consequential or special, 
that arise out of or in connection with the access, use or performance of WRF, including infringement actions.
WRF® is a registered trademark of the University Corporation for Atmospheric Research (UCAR).
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

#ifndef __GMS_LAND_VEGETATION_TYPES_H__
#define __GMS_LAND_VEGETATION_TYPES_H__ 101120261316

#include <cstdint>
#include <optional>
#include <utility>
#include <iostream>
#include "GMS_config.h"
#include "GMS_boost_multiarray_types.h"
#include "GMS_module_land_params.h"
#include "GMS_avx512_memset.h"
#include "GMS_avx512_memcpy.h"


namespace file_info 
{

     static const unsigned int GMS_LAND_VEGETATION_TYPES_MAJOR = 1;
     static const unsigned int GMS_LAND_VEGETATION_TYPES_MINOR = 1;
     static const unsigned int GMS_LAND_VEGETATION_TYPES_MICRO = 0;
     static const unsigned int GMS_LAND_VEGETATION_TYPES_FULLVER =
       1000U*GMS_LAND_VEGETATION_TYPES_MAJOR+100U*GMS_LAND_VEGETATION_TYPES_MINOR+
       10U*GMS_LAND_VEGETATION_TYPES_MICRO;
     static const char GMS_LAND_VEGETATION_TYPES_CREATION_DATE[] = "10-11-2025 13:16 +00200 (TUE 10 NOV 2026 GMT+2)";
     static const char GMS_LAND_VEGETATION_TYPES_BUILD_DATE[]    = __DATE__; 
     static const char GMS_LAND_VEGETATION_TYPES_BUILD_TIME[]    = __TIME__;
     static const char GMS_LAND_VEGETATION_TYPES_SYNOPSIS[]      = "Module Land Vegetation Types (based on the WRF 'module_sf_clm.F').";

}

/*
#include <iostream>
#include "boost/multi_array.hpp"

int main() {
    // 1. Get your sizes at runtime (e.g., from user input, file, or calculations)
    int rows, cols, depth;
    std::cout << "Enter dimensions (rows cols depth): ";
    std::cin >> rows >> cols >> depth;

    // 2. Define the array. 
    // The dimension count (3) is compile-time, but the sizes are runtime.
    boost::multi_array<double, 3> my_array(boost::extents[rows][cols][depth]);

    // 3. Use the array as normal
    my_array[0][0][0] = 42.0;
    std::cout << "Stored value: " << my_array[0][0][0] << std::endl;

    return 0;
}

// Create an empty 2D array (0x0)
boost::multi_array<int, 2> dynamic_matrix;

int new_rows = 10;
int new_cols = 20;

// Resize it dynamically at runtime
dynamic_matrix.resize(boost::extents[new_rows][new_cols]);

*/

namespace gms
{

namespace radiolocation
{

using namespace common;

struct alignas(64) vegetation_state_avx512_t
{
    std::size_t m_beg;
    std::size_t m_end;
    std::size_t m_nband; // e.g. X, Ku
    std::size_t m_nlevgrnd;
    Aligned64Vec1DI32 m_frac_veg_nosno;      //fraction of vegetation not covered by snow (0 OR 1) [-]
    Aligned64Vec1DI32 m_frac_veg_nosno_alb;  //!fraction of vegetation not covered by snow (0 OR 1) [-]
    Aligned64Vec1DF64 m_emv;                 //!vegetation emissivity
    Aligned64Vec1DF64 m_z0mv;                //roughness length over vegetation, momentum [m]
    Aligned64Vec1DF64 m_z0hv;                //!roughness length over vegetation, sensible heat [m]
    Aligned64Vec1DF64 m_z0qv;                //!roughness length over vegetation, latent heat [m]
    Aligned64Mat2DF64 m_rootfr;
    Aligned64Mat2DF64 m_rootr;
    Aligned64Mat2DF64 m_rresis;
    Aligned64Vec1DF64 m_dewmx;
    std::optional<Aligned64Vec1DF64> m_rssun;
    std::optional<Aligned64Vec1DF64> m_rssha;
    std::optional<Aligned64Vec1DF64> m_laisun;
    std::optional<Aligned64Vec1DF64> m_laisha;
    Aligned64Vec1DF64 m_btran;
    std::optional<Aligned64Vec1DF64> m_fsun;
    Aligned64Vec1DF64 m_tlai;
    Aligned64Vec1DF64 m_tsai;
    Aligned64Vec1DF64 m_elai;
    Aligned64Vec1DF64 m_esai;
    Aligned64Vec1DF64 m_fwet;
    Aligned64Vec1DF64 m_fdry;
    Aligned64Vec1DF64 m_dt_veg;
    Aligned64Vec1DF64 m_htop;
    Aligned64Vec1DF64 m_hbtop;
    Aligned64Vec1DF64 m_z0m;
    Aligned64Vec1DF64 m_displa;
    std::optional<Aligned64Mat2DF64> m_albd;
    std::optional<Aligned64Mat2DF64> m_albi;
    std::optional<Aligned64Mat2DF64> m_fabd;
    std::optional<Aligned64Mat2DF64> m_fabi;
    std::optional<Aligned64Mat2DF64> m_ftdd;
    std::optional<Aligned64Mat2DF64> m_ftid;
    std::optional<Aligned64Mat2DF64> m_ftii;
    Aligned64Vec1DF64 m_u10;
    Aligned64Vec1DF64 m_ram1;
    Aligned64Vec1DF64 m_fv;
    Aligned64Vec1DF64 m_forc_hgt_u_pft;
    Aligned64Vec1DF64 m_forc_hgt_t_pft;
    Aligned64Vec1DF64 m_forc_hgt_q_pft;
    Aligned64Vec1DF64 m_vds;
    Aligned64Vec1DF64 m_slasun;
    Aligned64Vec1DF64 m_slasha;
    Aligned64Vec1DF64 m_lncsun;
    Aligned64Vec1DF64 m_lncsha;
    Aligned64Vec1DF64 m_vcmxsun;
    Aligned64Vec1DF64 m_vcmxsha;
    // Can be adapted for the Radar Multi-Beam illumination
    Aligned64Vec1DF64 m_gdir;
    Aligned64Mat2DF64 m_omega;
    Aligned64Mat2DF64 m_eff_kid;
    Aligned64Mat2DF64 m_eff_kii;
    Aligned64Mat2DF64 m_sun_faid;
    Aligned64Mat2DF64 m_sun_faii;
    Aligned64Mat2DF64 m_sha_faid;
    Aligned64Mat2DF64 m_sha_faii;
    Aligned64Vec1DF64 m_sandfrac;
    Aligned64Vec1DF64 m_clayfrac;
    Aligned64Vec1DF64 m_soilt;
    Aligned64Vec1DF64 m_smoist;
    Aligned64Vec1DF64 m_sz;
    Aligned64Vec1DF64 m_watopt;
    Aligned64Vec1DF64 m_watdry;
    Aligned64Vec1DF64 m_Wlim;
    Aligned64Vec1DF64 m_tday;
    Aligned64Vec1DF64 m_tdayavg;
    Aligned64Vec1DI32 m_noveg;
    Aligned64Vec1DI32 m_tree;
    Aligned64Vec1DF64 m_dleaf;
    Aligned64Vec1DF64 m_xl;
    // Energy members maybe adequate for the Radar simulation.
    Aligned64Mat2DF64 m_sun_add;          //!sun canopy absorbed direct from direct (W/m**2)
    Aligned64Mat2DF64 m_tot_aid;          //!total canopy absorbed indirect from direct (W/m**2)
    Aligned64Mat2DF64 m_sun_aid;          //!sun canopy absorbed indirect from direct (W/m**2)
    Aligned64Mat2DF64 m_sun_aii;          //!sun canopy absorbed indirect from indirect (W/m**2)
    Aligned64Mat2DF64 m_sha_aid;          //!shade canopy absorbed indirect from direct (W/m**2)
    Aligned64Mat2DF64 m_sha_aii;          //!shade canopy absorbed indirect from indirect (W/m**2)
    Aligned64Mat2DF64 m_sun_atot;         //!sun canopy total absorbed (W/m**2)
    Aligned64Mat2DF64 m_sha_atot;         //!shade canopy total absorbed (W/m**2)
    Aligned64Mat2DF64 m_sun_alf;          //!sun canopy total absorbed by leaves (W/m**2)
    Aligned64Mat2DF64 m_sha_alf;          //!shade canopy total absored by leaves (W/m**2)

    vegetation_state_avx512_t() = default;
    
    vegetation_state_avx512_t(const std::size_t beg,
                              const std::size_t end,
                              const std::size_t nband,
                              const std::size_t nlevgrnd) noexcept(false) 
        :
    m_beg(beg),m_end(end),m_nband(nband),m_nlevgrnd(nlevgrnd),
    m_frac_veg_nosno(boost::extents[m_end-m_beg]),m_frac_veg_nosno_alb(boost::extents[m_end-m_beg]),
    m_emv(boost::extents[m_end-m_beg]),m_z0mv(boost::extents[m_end-m_beg]),m_z0hv(boost::extents[m_end-m_beg]),
    m_z0qv(boost::extents[m_end-m_beg]),m_rootfr(boost::extents[m_nlevgrnd][m_end-m_beg]),
    m_rootr(boost::extents[m_nlevgrnd][m_end-m_beg]),m_rresis(boost::extents[m_nlevgrnd][m_end-m_beg]),
    m_dewmx(boost::extents[m_end-m_beg]),m_rssun(std::in_place,boost::extents[m_end-m_beg]),
    m_rssha(std::in_place,boost::extents[m_end-m_beg]),m_laisun(std::in_place,boost::extents[m_end-m_beg]),
    m_laisha(std::in_place,boost::extents[m_end-m_beg]),m_btran(boost::extents[m_end-m_beg]),
    m_fsun(std::in_place,boost::extents[m_end-m_beg]),m_tlai(boost::extents[m_end-m_beg]),
    m_tsai(boost::extents[m_end-m_beg]),m_elai(boost::extents[m_end-m_beg]),m_esai(boost::extents[m_end-m_beg]),
    m_fwet(boost::extents[m_end-m_beg]),m_fdry(boost::extents[m_end-m_beg]),m_dt_veg(boost::extents[m_end-m_beg]),
    m_htop(boost::extents[m_end-m_beg]),m_hbtop(boost::extents[m_end-m_beg]),m_z0m(boost::extents[m_end-m_beg]),
    m_displa(boost::extents[m_end-m_beg]),m_albd(std::in_place,boost::extents[m_nband][m_end-m_beg]),
    m_albi(std::in_place,boost::extents[m_nband][m_end-m_beg]),m_fabd(std::in_place,boost::extents[m_nband][m_end-m_beg]),
    m_fabi(std::in_place,boost::extents[m_nband][m_end-m_beg]),m_ftdd(std::in_place,boost::extents[m_nband][m_end-m_beg]),
    m_ftid(std::in_place,boost::extents[m_nband][m_end-m_beg]),m_ftii(std::in_place,boost::extents[m_nband][m_end-m_beg]),
    m_u10(boost::extents[m_end-m_beg]),m_ram1(boost::extents[m_end-m_beg]),m_fv(boost::extents[m_end-m_beg]),
    m_forc_hgt_u_pft(boost::extents[m_end-m_beg]),m_forc_hgt_t_pft(boost::extents[m_end-m_beg]),
    m_forc_hgt_q_pft(boost::extents[m_end-m_beg]),m_vds(boost::extents[m_end-m_beg]),
    m_slasun(boost::extents[m_end-m_beg]),m_slasha(boost::extents[m_end-m_beg]),m_lncsun(boost::extents[m_end-m_beg]),
    m_lncsha(boost::extents[m_end-m_beg]),m_vcmxsun(boost::extents[m_end-m_beg]),m_vcmxsha(boost::extents[m_end-m_beg]),
    m_gdir(boost::extents[m_end-m_beg]),m_omega(boost::extents[m_nband][m_end-m_beg]),m_eff_kid(boost::extents[m_nband][m_end-m_beg]),
    m_eff_kii(boost::extents[m_nband][m_end-m_beg]),m_sun_faid(boost::extents[m_nband][m_end-m_beg]),
    m_sun_faii(boost::extents[m_nband][m_end-m_beg]),m_sha_faid(boost::extents[m_nband][m_end-m_beg]),
    m_sha_faii(boost::extents[m_nband][m_end-m_beg]),m_sandfrac(boost::extents[m_end-m_beg]),m_clayfrac(boost::extents[m_end-m_beg]),
    m_soilt(boost::extents[m_end-m_beg]),m_smoist(boost::extents[m_end-m_beg]),m_sz(boost::extents[m_end-m_beg]),
    m_watopt(boost::extents[m_end-m_beg]),m_watdry(boost::extents[m_end-m_beg]),m_Wlim(boost::extents[m_end-m_beg]),
    m_tday(boost::extents[m_end-m_beg]),m_tdayavg(boost::extents[m_end-m_beg]),m_noveg(boost::extents[m_end-m_beg]),
    m_tree(boost::extents[m_end-m_beg]),m_dleaf(boost::extents[m_end-m_beg]),m_xl(boost::extents[m_end-m_beg]),
    m_sun_add(boost::extents[m_nband][m_end-m_beg]),m_tot_aid(boost::extents[m_nband][m_end-m_beg]),
    m_sun_aid(boost::extents[m_nband][m_end-m_beg]),m_sun_aii(boost::extents[m_nband][m_end-m_beg]),
    m_sha_aid(boost::extents[m_nband][m_end-m_beg]),m_sha_aii(boost::extents[m_nband][m_end-m_beg]),
    m_sun_atot(boost::extents[m_nband][m_end-m_beg]),m_sha_atot(boost::extents[m_nband][m_end-m_beg]),
    m_sun_alf(boost::extents[m_nband][m_end-m_beg]),m_sha_alf(boost::extents[m_nband][m_end-m_beg])
    {

    }

    vegetation_state_avx512_t(const vegetation_state_avx512_t &other) noexcept(false) = default;

    vegetation_state_avx512_t(vegetation_state_avx512_t &&other) noexcept(true) = default;

    ~vegetation_state_avx512_t() = default;

    vegetation_state_avx512_t & 
    operator=(const vegetation_state_avx512_t &other) noexcept(false) = default;

    vegetation_state_avx512_t & 
    operator=(vegetation_state_avx512_t &&other) noexcept(true) = default;
    
    template<std::int32_t which_size>
    std::size_t total_size() const 
    {   
        std::size_t row_size = this->m_end-this->m_beg;
        if constexpr(which_size==1)
        {
            return (row_size*this->m_nband);
        }
        else if constexpr(which_size==2)
        {
            return (row_size*this->m_nlevgrnd);
        }
        else if constexpr(which_size==3)
        {
            return (row_size);
        }
    }

    std::int32_t show_class_memory_layout()
    {
        std::int32_t ret_status = -9999;
        std::cout << "alignof(struct vegetation_state_avx512_t) = " << alignof(vegetation_state_avx512_t) << '\n'
                   << "sizeof(struct  vegetation_state_avx512_t) = " << sizeof(vegetation_state_avx512_t)  << '\n'
                   << std::hex << std::showbase        << '\n'
                   << "&this->m_frac_veg_nosno        ="  << (void*)&this->m_frac_veg_nosno     << "\n"
                      "&this->m_frac_veg_nosno_alb    ="  << (void*)&this->m_frac_veg_nosno_alb << "\n"
                      "&this->m_emv                   ="  << (void*)&this->m_emv                << "\n"
                      "&this->m_z0mv                  ="  << (void*)&this->m_z0mv               << "\n"
                      "&this->m_z0hv                  ="  << (void*)&this->m_z0hv               << "\n"
                      "&this->m_z0qv                  ="  << (void*)&this->m_z0qv               << "\n"
                      "&this->m_rootfr                ="  << (void*)&this->m_rootfr             << "\n"
                      "&this->m_rresis                ="  << (void*)&this->m_rresis             << "\n"
                      "&this->m_dewmx                 ="  << (void*)&this->m_dewmx              << "\n"
                      "&this->m_rssun                 ="  << (void*)&this->m_rssun              << "\n"
                      "&this->m_rssha                 ="  << (void*)&this->m_rssha              << "\n"
                      "&this->m_laisun                ="  << (void*)&this->m_laisun             << "\n"
                      "&this->m_laisha                ="  << (void*)&this->m_laisha             << "\n"
                      "&this->m_btran                 ="  << (void*)&this->m_btran              << "\n"
                      "&this->m_fsun                  ="  << (void*)&this->m_fsun               << "\n"
                      "&this->m_tlai                  ="  << (void*)&this->m_tlai               << "\n"
                      "&this->m_tsai                  ="  << (void*)&this->m_tsai               << "\n"
                      "&this->m_elai                  ="  << (void*)&this->m_elai               << "\n"
                      "&this->m_esai                  ="  << (void*)&this->m_esai               << "\n"
                      "&this->m_fwet                  ="  << (void*)&this->m_fwet               << "\n"
                      "&this->m_fdry                  ="  << (void*)&this->m_fdry               << "\n"
                      "&this->m_dt_veg                ="  << (void*)&this->m_dt_veg             << "\n"
                      "&this->m_htop                  ="  << (void*)&this->m_htop               << "\n"
                      "&this->m_hbtop                 ="  << (void*)&this->m_hbtop              << "\n"
                      "&this->m_z0m                   ="  << (void*)&this->m_z0m                << "\n"
                      "&this->m_displa                ="  << (void*)&this->m_displa             << "\n"
                      "&this->m_albd                  ="  << (void*)&this->m_albd               << "\n"
                      "&this->m_albi                  ="  << (void*)&this->m_albi               << "\n"
                      "&this->m_fabd                  ="  << (void*)&this->m_fabd               << "\n"
                      "&this->m_fabi                  ="  << (void*)&this->m_fabi               << "\n"
                      "&this->m_ftdd                  ="  << (void*)&this->m_ftdd               << "\n"
                      "&this->m_ftid                  ="  << (void*)&this->m_ftid               << "\n"
                      "&this->m_ftii                  ="  << (void*)&this->m_ftii               << "\n"
                      "&this->m_u10                   ="  << (void*)&this->m_u10                << "\n"
                      "&this->m_ram1                  ="  << (void*)&this->m_ram1               << "\n"
                      "&this->m_fv                    ="  << (void*)&this->m_fv                 << "\n"
                      "&this->m_forc_hgt_u_pft        ="  << (void*)&this->m_forc_hgt_u_pft     << "\n"
                      "&this->m_forc_hgt_t_pft        ="  << (void*)&this->m_forc_hgt_t_pft     << "\n"
                      "&this->m_forc_hgt_q_pft        ="  << (void*)&this->m_forc_hgt_q_pft     << "\n"
                      "&this->m_vds                   ="  << (void*)&this->m_vds                << "\n"
                      "&this->m_slasun                ="  << (void*)&this->m_slasun             << "\n"
                      "&this->m_slasha                ="  << (void*)&this->m_slasha             << "\n"
                      "&this->m_lncsun                ="  << (void*)&this->m_lncsun             << "\n"
                      "&this->m_lncsha                ="  << (void*)&this->m_lncsha             << "\n"
                      "&this->m_vcmxsun               ="  << (void*)&this->m_vcmxsun            << "\n"
                      "&this->m_vcmxsha               ="  << (void*)&this->m_vcmxsha            << "\n"
                      "&this->m_gdir                  ="  << (void*)&this->m_gdir               << "\n"
                      "&this->m_omega                 ="  << (void*)&this->m_omega              << "\n"
                      "&this->m_eff_kid               ="  << (void*)&this->m_eff_kid            << "\n"
                      "&this->m_eff_kii               ="  << (void*)&this->m_eff_kii            << "\n"
                      "&this->m_sun_faid              ="  << (void*)&this->m_sun_faid           << "\n"
                      "&this->m_sun_faii              ="  << (void*)&this->m_sun_faii           << "\n"
                      "&this->m_sha_faid              ="  << (void*)&this->m_sha_faid           << "\n"
                      "&this->m_sha_faii              ="  << (void*)&this->m_sha_faii           << "\n"
                      "&this->m_sandfrac              ="  << (void*)&this->m_sandfrac           << "\n"
                      "&this->m_clayfrac              ="  << (void*)&this->m_clayfrac           << "\n"
                      "&this->m_soilt                 ="  << (void*)&this->m_soilt              << "\n"
                      "&this->m_smoist                ="  << (void*)&this->m_smoist             << "\n"
                      "&this->m_sz                    ="  << (void*)&this->m_sz                 << "\n"
                      "&this->m_watopt                ="  << (void*)&this->m_watopt             << "\n"
                      "&this->m_watdry                ="  << (void*)&this->m_watdry             << "\n"
                      "&this->m_Wlim                  ="  << (void*)&this->m_Wlim               << "\n"
                      "&this->m_tday                  ="  << (void*)&this->m_tday               << "\n"
                      "&this->m_tdayavg               ="  << (void*)&this->m_tdayavg            << "\n"
                      "&this->m_noveg                 ="  << (void*)&this->m_noveg              << "\n"
                      "&this->m_tree                  ="  << (void*)&this->m_tree               << "\n"
                      "&this->m_dleaf                 ="  << (void*)&this->m_dleaf              << "\n"
                      "&this->m_xl                    ="  << (void*)&this->m_xl                 << "\n"
                      "&this->m_sun_add               ="  << (void*)&this->m_sun_add            << "\n"
                      "&this->m_tot_aid               ="  << (void*)&this->m_tot_aid            << "\n"
                      "&this->m_sun_aid               ="  << (void*)&this->m_sun_aid            << "\n"
                      "&this->m_sun_aii               ="  << (void*)&this->m_sun_aii            << "\n"
                      "&this->m_sha_aid               ="  << (void*)&this->m_sha_aid            << "\n"
                      "&this->m_sha_aii               ="  << (void*)&this->m_sha_aii            << "\n"
                      "&this->m_sun_atot              ="  << (void*)&this->m_sun_atot           << "\n"
                      "&this->m_sha_atot              ="  << (void*)&this->m_sha_atot           << "\n"
                      "&this->m_sun_alf               ="  << (void*)&this->m_sun_alf            << "\n"
                      "&this->m_sha_alf               ="  << (void*)&this->m_sha_alf            << "\n";
        ret_status = 0;
        return (ret_status);    
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_frac_veg_nosno(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_frac_veg_nosno.data()[0];
        if constexpr(which_unrolling==16)
        {
           gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_frac_veg_nosno.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
           gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_frac_veg_nosno.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }
    
    template<std::int32_t which_unrolling>
    std::int32_t 
    init_frac_veg_nosno_alb(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_frac_veg_nosno_alb.data()[0];
        if constexpr(which_unrolling==16)
        {
           gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_frac_veg_nosno_alb.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
           gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_frac_veg_nosno_alb.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_emv(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_emv.data()[0];
        if constexpr(which_unrolling==16)
        {
           gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_emv.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
           gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_emv.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    set_emv(const Aligned64Vec1DF64 &emv) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        if(__builtin_expect(this->m_emv.num_elements()!=emv.num_elements(),0)) { ret_status = -1; return (ret_status);}
        double * __restrict__ this_emv  = &this->m_emv.data()[0];
        double * __restrict__ other_emv = &emv.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memcpy_unroll16x_pd(&this_emv[0],&other_emv[0],this->m_emv.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memcpy_unroll8x_pd(&this_emv[0],&other_emv[0],this->m_emv.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }
    
    template<std::int32_t which_unrolling>
    std::int32_t 
    init_z0mv(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_z0mv.data()[0];
        if constexpr(which_unrolling==16)
        {
           gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_z0mv.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
           gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_z0mv.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    set_z0mv(const Aligned64Vec1DF64 &z0mv) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        if(__builtin_expect(this->m_z0mv.num_elements()!=z0mv.num_elements(),0)) { ret_status = -1; return (ret_status);}
        double * __restrict__ this_z0mv  = &this->m_z0mv.data()[0];
        double * __restrict__ other_z0mv = &z0mv.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memcpy_unroll16x_pd(&this_z0mv[0],&other_z0mv[0],this->m_z0mv.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memcpy_unroll8x_pd(&this_z0mv[0],&other_z0mv[0],this->m_z0mv.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }



    template<std::int32_t which_unrolling>
    std::int32_t 
    init_z0hv(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_z0hv.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_z0hv.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_z0hv.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }
    
    template<std::int32_t which_unrolling>
    std::int32_t 
    init_rootfr(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_rootfr.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_rootfr.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_rootfr.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_rootr(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_rootr.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_rootr.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_rootr.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }
    
    template<std::int32_t which_unrolling>
    std::int32_t 
    init_rresis(const double filler) noexcept(true)
    {   
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_rresis.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_rresis.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_rresis.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }
    
    template<std::int32_t which_unrolling>
    std::int32_t 
    init_dewmx(const double filler) noexcept(true)
    {   
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_dewmx.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_dewmx.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_dewmx.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }
    
    template<std::int32_t which_unrolling>
    std::int32_t 
    init_rssun(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_rssun->data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_rssun->num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_rssun->num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_rssha(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_rssha->data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_rssha->num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_rssha->num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_laisun(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_laisun->data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_laisun->num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_laisun->num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_laisha(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_laisha->data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_laisha->num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_laisha->num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_btran(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_btran.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_btran.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_btran.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_fsun(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_fsun->data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_fsun->num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_fsun->num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_tlai(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_tlai.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_tlai.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_tlai.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_tsai(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_tsai.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_tsai.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_tsai.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_elai(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_elai.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_elai.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_elai.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_esai(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_esai.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_esai.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_esai.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_fwet(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_fwet.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_fwet.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_fwet.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_fdry(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_fdry.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_fdry.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_fdry.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_dt_veg(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_dt_veg.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_dt_veg.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_dt_veg.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_htop(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_htop.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_htop.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_htop.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_hbtop(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_fdry.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_hbtop.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_hbtop.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_z0m(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_z0m.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_z0m.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_z0m.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_displa(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_displa.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_displa.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_displa.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_albd(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_albd->data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_albd->num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_albd->num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_albi(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_albi->data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_albi->num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_albi->num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_fabd(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_fabd->data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_fabd->num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_fabd->num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_fabi(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_fabi->data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_fabi->num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_fabi->num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_ftdd(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_ftdd->data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_ftdd->num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_ftdd->num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_ftid(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_ftid->data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_ftid->num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_ftid->num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_ftii(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_ftii->data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_ftii->num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_ftii->num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_u10(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_u10.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_u10.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_u10.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_ram1(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_ram1.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_ram1.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_ram1.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_fv(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_fv.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_fv.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_fv.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_forc_hgt_u_pft(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_forc_hgt_u_pft.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_forc_hgt_u_pft.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_forc_hgt_u_pft.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_forc_hgt_t_pft(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_forc_hgt_t_pft.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_forc_hgt_t_pft.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_forc_hgt_t_pft.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_forc_hgt_q_pft(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_forc_hgt_q_pft.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_forc_hgt_q_pft.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_forc_hgt_q_pft.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_vds(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_vds.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_vds.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_vds.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

     template<std::int32_t which_unrolling>
    std::int32_t 
    init_slasun(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_slasun.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_slasun.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_slasun.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_slasha(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_slasha.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_slasha.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_slasha.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_lncsun(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_lncsun.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_lncsun.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_lncsun.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_lncsha(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_lncsha.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_lncsha.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_lncsha.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_vcmxsun(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_vcmxsun.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_vcmxsun.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_vcmxsun.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_vcmxsha(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_vcmxsha.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_vcmxsha.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_vcmxsha.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_gdir(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_gdir.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_gdir.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_gdir.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_omega(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_omega.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_omega.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_omega.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_eff_kid(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_eff_kid.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_eff_kid.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_eff_kid.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_eff_kii(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_eff_kii.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_eff_kii.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_eff_kii.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_sun_faid(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_sun_faid.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_sun_faid.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_sun_faid.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_sun_faii(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_sun_faii.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_sun_faii.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_sun_faii.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_sun_sandfrac(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_sandfrac.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_sandfrac.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_sandfrac.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_clayfrac(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_clayfrac.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_clayfrac.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_clayfrac.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_solit(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_soilt.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_soilt.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_soilt.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_smoist(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_smoist.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_smoist.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_smoist.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_sz(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_sz.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_sz.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_sz.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_watopt(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_watdry.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_watopt.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_watdry.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_Wlim(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_Wlim.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_Wlim.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_Wlim.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_tday(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_tday.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_tday.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_tday.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_tdayavg(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_tdayavg.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_tdayavg.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_tdayavg.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_noveg(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_noveg.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_noveg.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_noveg.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_tree(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_tree.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_tree.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_tree.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

     template<std::int32_t which_unrolling>
    std::int32_t 
    init_dleaf(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_dleaf.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_dleaf.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_dleaf.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_xl(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_xl.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_xl.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_xl.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_sun_add(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_sun_add.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_sun_add.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_sun_add.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_tot_aid(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_tot_aid.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_tot_aid.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_tot_aid.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_sun_aid(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_sun_aid.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_sun_aid.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_sun_aid.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_sun_aii(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_sun_aii.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_sun_aii.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_sun_aii.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_sha_aid(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_sha_aid.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_sha_aid.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_sha_aid.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_sha_aii(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_sha_aii.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_sha_aii.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_sha_aii.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_sun_atot(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_sun_atot.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_sun_atot.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_sun_atot.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_sha_atot(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_sha_atot.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_sha_atot.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_sha_atot.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_sun_alf(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_sun_alf.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_sun_alf.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_sun_alf.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

    template<std::int32_t which_unrolling>
    std::int32_t 
    init_sha_alf(const double filler) noexcept(true)
    {
        static_assert((which_unrolling==16) || (which_unrolling==8), "The template argument 'which_unrolling' is ***Incorrect***");
        std::int32_t ret_status = -9999;
        double * __restrict__ ptr = &this->m_sha_alf.data()[0];
        if constexpr(which_unrolling==16)
        {
            gms::common::avx512_memset_unroll16x_pd(&ptr[0],filler,this->m_sha_alf.num_elements());
        }
        else if constexpr(which_unrolling==8)
        {
            gms::common::avx512_memset_unroll8x_pd(&ptr[0],filler,this->m_sha_alf.num_elements());
        }
        ret_status = 0;
        return (ret_status);
    }

};

}

}
#endif /*__GMS_LAND_VEGETATION_TYPES_H__*/