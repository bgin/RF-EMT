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

#ifndef __GMS_LAND_VEGETATION_H__
#define __GMS_LAND_VEGETATION_H__ 101120261316

#include <cstdint>
#include <optional>
#include <utility>
#include "GMS_config.h"
#include "GMS_dyn_array.h"
#include "GMS_module_land_params.h"
#include "GMS_avx512_memset.h"


namespace file_info 
{

     static const unsigned int GMS_LAND_VEGETATION_MAJOR = 1;
     static const unsigned int GMS_LAND_VEGETATION_MINOR = 0;
     static const unsigned int GMS_LAND_VEGETATION_MICRO = 0;
     static const unsigned int GMS_LAND_VEGETATION_FULLVER =
       1000U*GMS_LAND_VEGETATION_MAJOR+100U*GMS_LAND_VEGETATION_MINOR+
       10U*GMS_LAND_VEGETATION_MICRO;
     static const char GMS_LAND_VEGETATION_CREATION_DATE[] = "10-11-2025 13:16 +00200 (TUE 10 NOV 2026 GMT+2)";
     static const char GMS_LAND_VEGETATION_BUILD_DATE[]    = __DATE__; 
     static const char GMS_LAND_VEGETATION_BUILD_TIME[]    = __TIME__;
     static const char GMS_LAND_VEGETATION_SYNOPSIS[]      = "Module Land Vegetation Types (based on the WRF 'module_sf_clm.F').";

}

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
    darray_r4_t m_frac_veg_nosno;      //fraction of vegetation not covered by snow (0 OR 1) [-]
    darray_r4_t m_frac_veg_nosno_alb;  //!fraction of vegetation not covered by snow (0 OR 1) [-]
    darray_r8_t m_emv;                 //!vegetation emissivity
    darray_r8_t m_z0mv;                //roughness length over vegetation, momentum [m]
    darray_r8_t m_z0hv;                //!roughness length over vegetation, sensible heat [m]
    darray_r8_t m_z0qv;                //!roughness length over vegetation, latent heat [m]
    darray_r8_t m_rootfr;
    darray_r8_t m_rootr;
    darray_r8_t m_rresis;
    darray_r8_t m_dewmx;
    std::optional<darray_r8_t> m_rssun;
    std::optional<darray_r8_t> m_rssha;
    std::optional<darray_r8_t> m_laisun;
    std::optional<darray_r8_t> m_laisha;
    darray_r8_t m_btran;
    std::optional<darray_r8_t> m_fsun;
    darray_r8_t m_tlai;
    darray_r8_t m_tsai;
    darray_r8_t m_elai;
    darray_r8_t m_esai;
    darray_r8_t m_fwet;
    darray_r8_t m_fdry;
    darray_r8_t m_dt_veg;
    darray_r8_t m_htop;
    darray_r8_t m_hbtop;
    darray_r8_t m_z0m;
    darray_r8_t m_displa;
    std::optional<darray_r8_t> m_albd;
    std::optional<darray_r8_t> m_albi;
    std::optional<darray_r8_t> m_fabd;
    std::optional<darray_r8_t> m_fabi;
    std::optional<darray_r8_t> m_ftdd;
    std::optional<darray_r8_t> m_ftid;
    std::optional<darray_r8_t> m_ftii;
    darray_r8_t m_u10;
    darray_r8_t m_ram1;
    darray_r8_t m_fv;
    darray_r8_t m_forc_hgt_u_pft;
    darray_r8_t m_forc_hgt_t_pft;
    darray_r8_t m_forc_hgt_q_pft;
    darray_r8_t m_vds;
    darray_r8_t m_slasun;
    darray_r8_t m_slasha;
    darray_r8_t m_lncsun;
    darray_r8_t m_lncsha;
    darray_r8_t m_vcmxsun;
    darray_r8_t m_vcmxsha;
    // Can be adapted for the Radar Multi-Beam illumination
    darray_r8_t m_gdir;
    darray_r8_t m_omega;
    darray_r8_t m_eff_kid;
    darray_r8_t m_eff_kii;
    darray_r8_t m_sun_faid;
    darray_r8_t m_sun_faii;
    darray_r8_t m_sha_faid;
    darray_r8_t m_sha_faii;
    darray_r8_t m_sandfrac;
    darray_r8_t m_clayfrac;
    darray_r8_t m_soilt;
    darray_r8_t m_smoist;
    darray_r8_t m_sz;
    darray_r8_t m_watopt;
    darray_r8_t m_watdry;
    darray_r8_t m_Wlim;
    darray_r8_t m_tday;
    darray_r8_t m_tdayavg;
    darray_r4_t m_noveg;
    darray_r4_t m_tree;
    darray_r8_t m_dleaf;
    darray_r8_t m_xl;
    // Energy members maybe adequate for the Radar simulation.
    darray_r8_t m_sun_add;          //!sun canopy absorbed direct from direct (W/m**2)
    darray_r8_t m_tot_aid;          //!total canopy absorbed indirect from direct (W/m**2)
    darray_r8_t m_sun_aid;          //!sun canopy absorbed indirect from direct (W/m**2)
    darray_r8_t m_sun_aii;          //!sun canopy absorbed indirect from indirect (W/m**2)
    darray_r8_t m_sha_aid;          //!shade canopy absorbed indirect from direct (W/m**2)
    darray_r8_t m_sha_aii;          //!shade canopy absorbed indirect from indirect (W/m**2)
    darray_r8_t m_sun_atot;         //!sun canopy total absorbed (W/m**2)
    darray_r8_t m_sha_atot;         //!shade canopy total absorbed (W/m**2)
    darray_r8_t m_sun_alf;          //!sun canopy total absorbed by leaves (W/m**2)
    darray_r8_t m_sha_alf;          //!shade canopy total absored by leaves (W/m**2)

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
    m_frac_veg_nosno(m_end-m_beg),m_frac_veg_nosno_alb(m_end-m_beg),
    m_emv(m_end-m_beg),m_z0mv(m_end-m_beg),m_z0hv(m_end-m_beg),
    m_z0qv(m_end-m_beg),m_rootfr(m_nlevgrnd*(m_end-m_beg)),
    m_rootr(m_nlevgrnd*(m_end-m_beg)),m_rresis(m_nlevgrnd*(m_end-m_beg)),
    m_dewmx(m_end-m_beg),m_rssun(std::in_place,m_end-m_beg),
    m_rssha(std::in_place,m_end-m_beg),m_laisun(std::in_place,m_end-m_beg),
    m_laisha(std::in_place,m_end-m_beg),m_btran(m_end-m_beg),
    m_fsun(std::in_place,m_end-m_beg),m_tlai(m_end-m_beg),
    m_tsai(m_end-m_beg),m_elai(m_end-m_beg),m_esai(m_end-m_beg),
    m_fwet(m_end-m_beg),m_fdry(m_end-m_beg),m_dt_veg(m_end-m_beg),
    m_htop(m_end-m_beg),m_hbtop(m_end-m_beg),m_z0m(m_end-m_beg),
    m_displa(m_end-m_beg),m_albd(std::in_place,m_nband*(m_end-m_beg)),
    m_albi(std::in_place,m_nband*(m_end-m_beg)),m_fabd(std::in_place,m_nband*(m_end-m_beg)),
    m_fabi(std::in_place,m_nband*(m_end-m_beg)),m_ftdd(std::in_place,m_nband*(m_end-m_beg)),
    m_ftid(std::in_place,m_nband*(m_end-m_beg)),m_ftii(std::in_place,m_nband*(m_end-m_beg)),
    m_u10(m_end-m_beg),m_ram1(m_end-m_beg),m_fv(m_end-m_beg),
    m_forc_hgt_u_pft(m_end-m_beg),m_forc_hgt_t_pft(m_end-m_beg),
    m_forc_hgt_q_pft(m_end-m_beg),m_vds(m_end-m_beg),
    m_slasun(m_end-m_beg),m_slasha(m_end-m_beg),m_lncsun(m_end-m_beg),
    m_lncsha(m_end-m_beg),m_vcmxsun(m_end-m_beg),m_vcmxsha(m_end-m_beg),
    m_gdir(m_end-m_beg),m_omega(m_nband*(m_end-m_beg)),m_eff_kid(m_nband*(m_end-m_beg)),
    m_eff_kii(m_nband*(m_end-m_beg)),m_sun_faid(m_nband*(m_end-m_beg)),
    m_sun_faii(m_nband*(m_end-m_beg)),m_sha_faid(m_nband*(m_end-m_beg)),
    m_sha_faii(m_nband*(m_end-m_beg)),m_sandfrac(m_end-m_beg),m_clayfrac(m_end-m_beg),
    m_soilt(m_end-m_beg),m_smoist(m_end-m_beg),m_sz(m_end-m_beg),
    m_watopt(m_end-m_beg),m_watdry(m_end-m_beg),m_Wlim(m_end-m_beg),
    m_tday(m_end-m_beg),m_tdayavg(m_end-m_beg),m_noveg(m_end-m_beg),
    m_tree(m_end-m_beg),m_dleaf(m_end-m_beg),m_xl(m_end-m_beg),
    m_sun_add(m_nband*(m_end-m_beg)),m_tot_aid(m_nband*(m_end-m_beg)),
    m_sun_aid(m_nband*(m_end-m_beg)),m_sun_aii(m_nband*(m_end-m_beg)),
    m_sha_aid(m_nband*(m_end-m_beg)),m_sha_aii(m_nband*(m_end-m_beg)),
    m_sun_atot(m_nband*(m_end-m_beg)),m_sha_atot(m_nband*(m_end-m_beg)),
    m_sun_alf(m_nband*(m_end-m_beg)),m_sha_alf(m_nband*(m_end-m_beg))
    {

    }


};

}

}

#endif /*__GMS_LAND_VEGETATION_H__*/