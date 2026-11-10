
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

#ifndef __GMS_MODULE_LAND_PARAMS_H__
#define __GMS_MODULE_LAND_PARAMS_H__ 101120261026

#include <cstdint>
#include "GMS_config.h"

namespace file_info 
{

     static const unsigned int GMS_MODULE_LAND_PARAMS_MAJOR = 1;
     static const unsigned int GMS_MODULE_LAND_PARAMS_MINOR = 1;
     static const unsigned int GMS_MODULE_LAND_PARAMS_MICRO = 0;
     static const unsigned int GMS_MODULE_LAND_PARAMS_FULLVER =
       1000U*GMS_MODULE_LAND_PARAMS_MAJOR+100U*GMS_MODULE_LAND_PARAMS_MINOR+
       10U*GMS_MODULE_LAND_PARAMS_MICRO;
     static const char GMS_MODULE_LAND_PARAMS_CREATION_DATE[] = "10-11-2025 10:44AM +00200 (TUE 10 NOV 2026 GMT+2)";
     static const char GMS_MODULE_LAND_PARAMS_BUILD_DATE[]    = __DATE__; 
     static const char GMS_MODULE_LAND_PARAMS_BUILD_TIME[]    = __TIME__;
     static const char GMS_MODULE_LAND_PARAMS_SYNOPSIS[]      = "Module Land constant parameters (based on the WRF 'module_sf_clm.F').";

}

namespace gms
{


namespace radiolocation 
{


namespace mod_land_varpar
{

constexpr const static std::int32_t lsmlon = 1; // maximum number of longitude points on lsm grid
constexpr const static std::int32_t lsmlat = 1; // number of latitude points on lsm grid

// Define the number of levels
constexpr const static std::int32_t nlevsoi = 10; // ! number of soil layers
constexpr const static std::int32_t nlevlak = 10; // number of lake layers
constexpr const static std::int32_t nlevsno = 5;  // maximum number of snow layers
constexpr const static std::int32_t nlevgrnd = 10; // number of ground layers (includes lower layers that are hydrologically inactive)
constexpr const static std::int32_t nlevurb  = nlevgrnd; //number of urban layers (must equal nlevgrnd right now)

// Miscellanous defines
constexpr const static std::int32_t numwat  = 5; //number of water types (soil, ice, 2 lakes, wetland)
constexpr const static std::int32_t numrad  = 2; // number of solar radiation bands: vis, nir
constexpr const static std::int32_t numsolar = 2; // number of solar type bands: direct, diffuse
constexpr const static std::int32_t ndst    = 4; // number of dust size classes (BGC only)
constexpr const static std::int32_t dst_src_nbr = 3; // number of size distns in src soil (BGC only)
constexpr const static std::int32_t sz_nbr  = 200; // number of sub-grid bins in large bin of dust size distribution (BGC only)
constexpr const static std::int32_t numcol  = 8;  // number of soil color types
constexpr const static std::int32_t rtmlon  = 720;// number of rtm longitudes
constexpr const static std::int32_t rtmlat  = 320;// number of rtm latitudes
#if defined(CROP) 
constexpr const static std::int32_t numpft  = 20;// number of plant types
constexpr const static std::int32_t numcft  = 6; // actual # of crops
constexpr const static std::int32_t numveg  = 16; //number of veg types (without specific crop)
#else 
constexpr const static std::int32_t numpft  = 16; //actual # of pfts (without bare)
constexpr const static std::int32_t numcft  = 2;  //actual # of crops
constexpr const static std::int32_t numveg  = numpft; // number of veg types (without specific crop)
#endif 
constexpr const static std::int32_t maxpatch_urb = 1;  //the current coupling not include urban
#if defined(CROP)
constexpr const static std::int32_t maxpatch_cft = 4; 
#else 
constexpr const static std::int32_t maxpatch_cft = 2;
#endif
constexpr const static std::int32_t maxpatch_pft   = 4;
constexpr const static std::int32_t npatch_urban   = maxpatch_pft + 1;
constexpr const static std::int32_t npatch_lake    = npatch_urban + maxpatch_urb;
constexpr const static std::int32_t npatch_wet     = npatch_lake  + 1;
constexpr const static std::int32_t npatch_glacier = npatch_wet   + 1;
constexpr const static std::int32_t npatch_crop    = npatch_glacier + maxpatch_cft;
constexpr const static std::int32_t maxpatch       = npatch_crop;

}

}

}

#endif /*__GMS_MODULE_LAND_PARAMS_H__*/