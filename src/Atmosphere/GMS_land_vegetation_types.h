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
    m_beg(other.m_beg),m_end(other.m_end),m_nband(other.m_nband),m_nlevgrnd(other.m_nlevgrnd),
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

     vegetation_state_avx512_t(vegetation_state_avx512_t &&other) noexcept(true) 
        :
    m_beg(other.m_beg),m_end(other.m_end),m_nband(other.m_nband),m_nlevgrnd(other.m_nlevgrnd),
    m_frac_veg_nosno(std::move(other.m_frac_veg_nosno)),m_frac_veg_nosno_alb(std::move(other.m_frac_veg_nosno_alb)),
    m_emv(std::move(other.m_emv)),m_z0mv(std::move(other.m_z0mv)),m_z0hv(std::move(other.m_z0hv)),
    m_z0qv(std::move(other.m_z0qv)),m_rootfr(std::move(other.m_rootfr)),
    m_rootr(std::move(other.m_rootr)),m_rresis(std::move(other.m_rresis)),
    m_dewmx(std::move(other.m_dewmx)),m_rssun(std::move(other.m_rssun)),
    m_rssha(std::move(other.m_rssha)),m_laisun(std::move(other.m_laisun)),
    m_laisha(std::move(other.m_laisha)),m_btran(std::move(other.m_btran)),
    m_fsun(std::move(other.m_fsun)),m_tlai(std::move(other.m_tlai)),
    m_tsai(std::move(other.m_tsai)),m_elai(std::move(other.m_elai)),m_esai(std::move(other.m_esai)),
    m_fwet(std::move(other.m_fwet)),m_fdry(std::move(other.m_fdry)),m_dt_veg(std::move(other.m_dt_veg)),
    m_htop(std::move(other.m_htop)),m_hbtop(std::move(other.m_hbtop)),m_z0m(std::move(other.m_z0m)),
    m_displa(std::move(other.m_displa)),m_albd(std::move(other.m_albd)),
    m_albi(std::move(other.m_albi)),m_fabd(std::move(other.m_fabd)),
    m_fabi(std::move(other.m_fabi)),m_ftdd(std::move(other.m_ftdd)),
    m_ftid(std::move(other.m_ftid)),m_ftii(std::move(other.m_ftii)),
    m_u10(std::move(other.m_u10)),m_ram1(std::move(other.m_ram1)),m_fv(std::move(other.m_fv)),
    m_forc_hgt_u_pft(std::move(other.m_forc_hgt_u_pft)),m_forc_hgt_t_pft(std::move(other.m_forc_hgt_t_pft)),
    m_forc_hgt_q_pft(std::move(other.m_forc_hgt_q_pft)),m_vds(std::move(other.m_vds)),
    m_slasun(std::move(other.m_slasun)),m_slasha(std::move(other.m_slasha)),m_lncsun(std::move(other.m_lncsun)),
    m_lncsha(std::move(other.m_lncsha)),m_vcmxsun(std::move(other.m_vcmxsun)),m_vcmxsha(std::move(other.m_vcmxsha)),
    m_gdir(std::move(other.m_gdir)),m_omega(std::move(other.m_omega)),m_eff_kid(std::move(other.m_eff_kid)),
    m_eff_kii(std::move(other.m_eff_kii)),m_sun_faid(std::move(other.m_sun_faid)),
    m_sun_faii(std::move(other.m_sun_faii)),m_sha_faid(std::move(other.m_sha_faid)),
    m_sha_faii(std::move(other.m_sha_faii)),m_sandfrac(std::move(other.m_sandfrac)),m_clayfrac(std::move(other.m_clayfrac)),
    m_soilt(std::move(other.m_soilt)),m_smoist(std::move(other.m_smoist)),m_sz(std::move(other.m_sz)),
    m_watopt(std::move(other.m_watopt)),m_watdry(std::move(other.m_watdry)),m_Wlim(std::move(other.m_Wlim)),
    m_tday(std::move(other.m_tday)),m_tdayavg(std::move(other.m_tdayavg)),m_noveg(std::move(other.m_noveg)),
    m_tree(std::move(other.m_tree)),m_dleaf(std::move(other.m_dleaf)),m_xl(std::move(other.m_xl)),
    m_sun_add(std::move(other.m_sun_add)),m_tot_aid(std::move(other.m_tot_aid)),
    m_sun_aid(std::move(other.m_sun_aid)),m_sun_aii(std::move(other.m_sun_aii)),
    m_sha_aid(std::move(other.m_sha_aid)),m_sha_aii(std::move(other.m_sha_aii)),
    m_sun_atot(std::move(other.m_sun_atot)),m_sha_atot(std::move(other.m_sha_atot)),
    m_sun_alf(std::move(other.m_sun_alf)),m_sha_alf(std::move(other.m_sha_alf))
    {

    }

    ~vegetation_state_avx512_t() noexcept(true) 
    {

    }

    friend void swap(vegetation_state_avx512_t &lhs, vegetation_state_avx512_t &rhs) noexcept(true) 
    {
        using std::swap;
        swap(lhs.m_beg,rhs.m_beg);
        swap(lhs.m_end,rhs.m_end);
        swap(lhs.m_nband,rhs.m_nband);
        swap(lhs.m_nlevgrnd,rhs.m_nlevgrnd);
        swap(lhs.m_frac_veg_nosno,rhs.m_frac_veg_nosno);
        swap(lhs.m_frac_veg_nosno_alb,rhs.m_frac_veg_nosno_alb);
        swap(lhs.m_emv,rhs.m_emv);
        swap(lhs.m_z0mv,rhs.m_z0mv);
        swap(lhs.m_z0hv,rhs.m_z0hv);
        swap(lhs.m_z0qv,rhs.m_z0qv);
        swap(lhs.m_rootfr,rhs.m_rootfr);
        swap(lhs.m_rresis,rhs.m_rresis);
        swap(lhs.m_dewmx,rhs.m_dewmx);
        swap(lhs.m_rssun,rhs.m_rssun);
        swap(lhs.m_rssha,rhs.m_rssha);
        swap(lhs.m_laisun,rhs.m_laisun);
        swap(lhs.m_laisha,lhs.m_laisha);
        swap(lhs.m_btran,rhs.m_btran);
        swap(lhs.m_fsun,rhs.m_fsun);
        swap(lhs.m_tlai,rhs.m_tlai);
        swap(lhs.m_tsai,rhs.m_tsai);
        swap(lhs.m_elai,rhs.m_elai);
        swap(lhs.m_esai,rhs.m_esai);
        swap(lhs.m_fwet,rhs.m_fwet);
        swap(lhs.m_fdry,rhs.m_fdry);
        swap(lhs.m_dt_veg,rhs.m_dt_veg);
        swap(lhs.m_htop,rhs.m_htop);
        swap(lhs.m_hbtop,rhs.m_hbtop);
        swap(lhs.m_z0m,rhs.m_z0m);
        swap(lhs.m_displa,rhs.m_displa);
        swap(lhs.m_albd,rhs.m_albd);
        swap(lhs.m_albi,rhs.m_albi);
        swap(lhs.m_fabd,rhs.m_fabd);
        swap(lhs.m_fabi,rhs.m_fabi);
        swap(lhs.m_ftdd,rhs.m_ftdd);
        swap(lhs.m_ftid,rhs.m_ftid);
        swap(lhs.m_ftii,rhs.m_ftii);
        swap(lhs.m_u10,rhs.m_u10);
        swap(lhs.m_ram1,rhs.m_ram1);
        swap(lhs.m_fv,rhs.m_fv);
        swap(lhs.m_forc_hgt_u_pft,rhs.m_forc_hgt_u_pft);
        swap(lhs.m_forc_hgt_t_pft,rhs.m_forc_hgt_t_pft);
        swap(lhs.m_forc_hgt_q_pft,rhs.m_forc_hgt_q_pft);
        swap(lhs.m_vds,rhs.m_vds);
        swap(lhs.m_slasun,rhs.m_slasun);
        swap(lhs.m_slasha,rhs.m_slasha);
        swap(lhs.m_lncsun,rhs.m_lncsun);
        swap(lhs.m_lncsha,rhs.m_lncsha);
        swap(lhs.m_vcmxsun,rhs.m_vcmxsun);
        swap(lhs.m_vcmxsha,rhs.m_vcmxsha);
        swap(lhs.m_gdir,rhs.m_gdir);
        swap(lhs.m_omega,rhs.m_omega);
        swap(lhs.m_eff_kid,rhs.m_eff_kid);
        swap(lhs.m_eff_kii,rhs.m_eff_kii);
        swap(lhs.m_sun_faid,rhs.m_sun_faid);
        swap(lhs.m_sun_faii,rhs.m_sun_faii);
        swap(lhs.m_sha_faid,rhs.m_sha_faid);
        swap(lhs.m_sha_faii,rhs.m_sha_faii);
        swap(lhs.m_sandfrac,rhs.m_sandfrac);
        swap(lhs.m_clayfrac,rhs.m_clayfrac);
        swap(lhs.m_soilt,rhs.m_soilt);
        swap(lhs.m_smoist,rhs.m_smoist);
        swap(lhs.m_sz,rhs.m_sz);
        swap(lhs.m_watopt,rhs.m_watopt);
        swap(lhs.m_watdry,rhs.m_watdry);
        swap(lhs.m_Wlim,rhs.m_Wlim);
        swap(lhs.m_tday,rhs.m_tday);
        swap(lhs.m_tdayavg,rhs.m_tdayavg);
        swap(lhs.m_noveg,rhs.m_noveg);
        swap(lhs.m_tree,rhs.m_tree);
        swap(lhs.m_dleaf,rhs.m_dleaf);
        swap(lhs.m_xl,rhs.m_xl);
        swap(lhs.m_sun_add,rhs.m_sun_add);
        swap(lhs.m_tot_aid,rhs.m_tot_aid);
        swap(lhs.m_sun_aid,rhs.m_sun_aid);
        swap(lhs.m_sun_aii,rhs.m_sun_aii);
        swap(lhs.m_sha_aid,rhs.m_sha_aid);
        swap(lhs.m_sha_aii,rhs.m_sha_aii);
        swap(lhs.m_sun_atot,rhs.m_sun_atot);
        swap(lhs.m_sha_atot,rhs.m_sha_atot);
        swap(lhs.m_sun_alf,rhs.m_sun_alf);
        swap(lhs.m_sha_alf,rhs.m_sha_alf);
    }

    vegetation_state_avx512_t & 
    operator=(const vegetation_state_avx512_t &other) noexcept(false)
    {
          if(__builtin_expect(this==&other,0)) { return (*this);}
          vegetation_state_avx512_t tmp{other};
          swap(tmp,*this);
          return (*this);
    }

    vegetation_state_avx512_t & 
    operator=(vegetation_state_avx512_t &&other) noexcept(true)
    {
          if(__builtin_expect(this==&other,0)) { return (*this);}
          swap(*this,other);
          return (*this);
    }

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