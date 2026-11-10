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

struct alignas(64) pft_state_type_avx512_t
{
    
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
};

}

}
#endif /*__GMS_LAND_VEGETATION_TYPES_H__*/