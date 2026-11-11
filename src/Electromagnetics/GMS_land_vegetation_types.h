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
#include "GMS_config.h"
#include "GMS_boost_multiarray_types.hpp"
#include "GMS_module_land_params.h"


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
    Aligned64Vec1DF64 m_omega;
    Aligned64Vec1DF64 m_eff_kid;
    Aligned64Vec1DF64 m_eff_kii;
    Aligned64Vec1DF64 m_sun_faid;
    Aligned64Vec1DF64 m_sun_faii;
    Aligned64Vec1DF64 m_sha_faid;
    Aligned64Vec1DF64 m_sha_faii;
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

    vegetation_state_avx512_t() noexcept
        : 
    m_beg{},m_end{},m_nband{},m_nlevgrnd{},
    m_frac_veg_nosno{},m_frac_veg_nosno_alb{},
    m_emv{},m_z0mv{},m_z0hv{},m_z0qv{},m_rootfr{},
    m_rootr{},m_rresis{},m_dewmx{},m_rssun(),
    m_rssha(),m_laisun(),m_laisha(),m_btran{},
    m_fsun(),m_tlai{},m_tsai{},m_elai{},m_esai{},
    m_fwet{},m_fdry{},m_dt_veg{},m_htop{},m_hbtop{},
    m_z0m{},m_displa{},m_albd(),m_albi(),m_fabd(),
    m_fabi(),m_ftdd(),m_ftid(),m_ftii(),m_u10{},
    m_ram1{},m_fv{},m_forc_hgt_u_pft{},
    m_forc_hgt_t_pft{},m_forc_hgt_q_pft{},m_vds{},
    m_slasun{},m_slasha{},m_lncsun{},m_lncsha{},
    m_vcmxsun{},m_vcmxsha{},m_gdir{},m_omega{},
    m_eff_kid{},m_eff_kii{},m_sun_faid{},m_sun_faii{},
    m_sha_faid{},m_sha_faii{},m_sandfrac{},m_clayfrac{},
    m_soilt{},m_smoist{},m_sz{},m_watopt{},m_watdry{},
    m_Wlim{},m_tday{},m_tdayavg{},m_noveg{},m_tree{},
    m_dleaf{},m_xl{},m_sun_add{},m_tot_aid{},
    m_sun_aid{},m_sun_aii{},m_sha_aid{},m_sha_aii{},
    m_sun_atot{},m_sha_atot{},m_sun_alf{},m_sha_alf{}
    {

    }
    
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

    vegetation_state_avx512_t(const vegetation_state_avx512_t &other) noexcept(false) 
        :
    m_beg(other.beg),m_end(other.end),m_nband(other.nband),m_nlevgrnd(other.nlevgrnd),
    m_frac_veg_nosno(other.m_frac_veg_nosno),m_frac_veg_nosno_alb(other.m_frac_veg_nosno_alb),
    m_emv(other.m_emv),m_z0mv(other.m_z0mv),m_z0hv(other.m_z0hv),
    m_z0qv(other.m_z0qv),m_rootfr(other.m_rootfr),
    m_rootr(other.m_rootr),m_rresis(other.m_rresis),
    m_dewmx(other.m_dewmx),m_rssun(other.m_rssun),
    m_rssha(other.m_rssha),m_laisun(other.m_laisun),
    m_laisha(other.m_laisha),m_btran(other.m_btran),
    m_fsun(other.m_fsun),m_tlai(other.m_tlai),
    m_tsai(other.m_tsai),m_elai(other.m_elai),m_esai(other.m_esai),
    m_fwet(other.m_fwet),m_fdry(other.m_fdry),m_dt_veg(other.m_dt_veg),
    m_htop(other.m_htop),m_hbtop(other.m_hbtop),m_z0m(other.m_z0m),
    m_displa(other.m_displa),m_albd(other.m_albd),
    m_albi(other.m_albi),m_fabd(other.m_fabd),
    m_fabi(other.m_fabi),m_ftdd(other.m_ftdd),
    m_ftid(other.m_ftid),m_ftii(other.m_ftii),
    m_u10(other.m_u10),m_ram1(other.m_ram1),m_fv(other.m_fv),
    m_forc_hgt_u_pft(other.m_forc_hgt_u_pft),m_forc_hgt_t_pft(other.m_forc_hgt_t_pft),
    m_forc_hgt_q_pft(other.m_forc_hgt_q_pft),m_vds(other.m_vds),
    m_slasun(other.m_slasun),m_slasha(other.m_slasha),m_lncsun(other.m_lncsun),
    m_lncsha(other.m_lncsha),m_vcmxsun(other.m_vcmxsun),m_vcmxsha(other.m_vcmxsha),
    m_gdir(other.m_gdir),m_omega(other.m_omega),m_eff_kid(other.m_eff_kid),
    m_eff_kii(other.m_eff_kii),m_sun_faid(other.m_sun_faid),
    m_sun_faii(other.m_sun_faii),m_sha_faid(other.m_sha_faid),
    m_sha_faii(other.m_sha_faii),m_sandfrac(other.m_sandfrac),m_clayfrac(other.m_clayfrac),
    m_soilt(other.m_soilt),m_smoist(other.m_smoist),m_sz(other.m_sz),
    m_watopt(other.m_watopt),m_watdry(other.m_watdry),m_Wlim(other.m_Wlim),
    m_tday(other.m_tday),m_tdayavg(other.m_tdayavg),m_noveg(other.m_noveg),
    m_tree(other.m_tree),m_dleaf(other.m_dleaf),m_xl(other.m_xl),
    m_sun_add(other.m_sun_add),m_tot_aid(other.m_tot_aid),
    m_sun_aid(other.m_sun_aid),m_sun_aii(other.m_sun_aii),
    m_sha_aid(other.m_sha_aid),m_sha_aii(other.m_sha_aii),
    m_sun_atot(other.m_sun_atot),m_sha_atot(other.m_sha_atot),
    m_sun_alf(other.m_sun_alf),m_sha_alf(other.m_sha_alf)
    {

    }

};

}

}
#endif /*__GMS_LAND_VEGETATION_TYPES_H__*/