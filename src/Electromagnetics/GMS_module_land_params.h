
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
constexpr const static std::int32_t maxpatch_cft = 2; //
#endif
constexpr const static std::int32_t maxpatch_pft   = 4;
constexpr const static std::int32_t npatch_urban   = maxpatch_pft + 1;
constexpr const static std::int32_t npatch_lake    = npatch_urban + maxpatch_urb;
constexpr const static std::int32_t npatch_wet     = npatch_lake  + 1;
constexpr const static std::int32_t npatch_glacier = npatch_wet   + 1;
constexpr const static std::int32_t npatch_crop    = npatch_glacier + maxpatch_cft;
constexpr const static std::int32_t maxpatch       = npatch_crop;
#if (defined CROP)
constexpr const static std::int32_t max_pft_per_gcell = numpft+1 + 3 + maxpatch_urb;
#else
constexpr const static std::int32_t max_pft_per_gcell = numpft+1+3+maxpatch_urb+numcft;
#endif
constexpr const static std::int32_t max_pft_per_lu    = std::max(numpft+1, std::max(numcft, maxpatch_urb));
constexpr const static std::int32_t max_pft_per_col   = std::max(numpft+1, std::max(numcft, maxpatch_urb));
std::int32_t num_landcover_types{}; 
} // mod_land_varpar

namespace mod_land_vapar
{

//std::int32_t i;
// Initialize physical constants
constexpr const static double cday   = 86400.0;                     //!sec in calendar day ~ sec
constexpr const static std::int32_t idx_Mie_snw_mx = 1471;           //! number of effective radius indices used in Mie lookup table [idx]
constexpr const static std::int32_t idx_T_max      = 11;            // ! maxiumum temperature index used in aging lookup table [idx]
constexpr const static std::int32_t idx_Tgrd_max   = 31;            // ! maxiumum temperature gradient index used in aging lookup table [idx]
constexpr const static std::int32_t idx_rhos_max   = 8;             // ! maxiumum snow density index used in aging lookup table [idx]
constexpr const static std::int32_t numrad_snw     = 5;             //  ! number of spectral bands used in snow model [nbr]
constexpr const static double pie    = 3.141592653589793;           //! pi
constexpr const static double rpi    = 3.141592653589793;           //! pi
constexpr const static double grav   = 9.80616;                     //!gravity constant [m/s2]
constexpr const static double sb     = 5.67e-8;                     //!stefan-boltzmann constant  [W/m2/K4]
constexpr const static double vkc    = 0.4;                         //!von Karman constant [-]
constexpr const static double rgas   = 8314.468;                    //!Universal gas constant ~ J/K/kmole
constexpr const static double rwat   = 461.5046;                    //!gas constant for water vapor [J/(kg K)]
constexpr const static double rair   = 287.0423;                    //!gas constant for dry air [J/kg/K]
constexpr const static double roverg = 47062.73;                    //!Rw/g constant = (8.3144/0.018)/(9.80616)*1000. mm/K
constexpr const static double cpliq  = 4.188e3;                     //!Specific heat of water [J/kg-K]
constexpr const static double cpice  = 2.11727e3;                   //!Specific heat of ice [J/kg-K]
constexpr const static double cpair  = 1.00464e3;                   //!specific heat of dry air [J/kg/K]
constexpr const static double hvap   = 2.501e6;                     //!Latent heat of evap for water [J/kg]
constexpr const static double hfus   = 3.337e5;                     //!Latent heat of fusion for ice [J/kg]
constexpr const static double hsub   = 2.501e6+3.337e5;             //!Latent heat of sublimation    [J/kg]
constexpr const static double denh2o = 1.000e3;                     //!density of liquid water [kg/m3]
constexpr const static double denice = 0.917e3;                     //!density of ice [kg/m3]
constexpr const static double tkair  = 0.023;                       //!thermal conductivity of air   [W/m/k]
constexpr const static double tkice  = 2.290;                       //!thermal conductivity of ice   [W/m/k]
constexpr const static double tkwat  = 0.6;                         // !thermal conductivity of water [W/m/k]
constexpr const static double tfrz   = 273.16;                      //!freezing temperature [K]
constexpr const static double tcrit  = 2.5;                         //!critical temperature to determine rain or snow
constexpr const static double po2    = 0.209;                       //!constant atmospheric partial pressure  O2 (mol/mol)
constexpr const static double pco2   = 355.e-06;                    //!constant atmospheric partial pressure CO2 (mol/mol)
constexpr const static double pstd   = 101325.0;                    //!standard pressure ~ pascals
constexpr const static double bdsno  = 250.;                        //!bulk density snow (kg/m**3)
constexpr const static double re     = 6.37122e6*0.001;             //!radius of earth (km)
constexpr const static double secspday= 86400.0;                    //! Seconds per day
constexpr const static double spval = 1.e36;                        //! special value for real data
constexpr const static std::int32_t ispval = -9999;                 //! special value for int data
constexpr const static double alpha_aero = 1.0;                     //!constant for aerodynamic parameter weighting
constexpr const static double tlsai_crit = 2.0;                     //!critical value of elai+esai for which aerodynamic parameters are maximum
constexpr const static double watmin = 0.01;                        //!minimum soil moisture (mm)
constexpr const static double zlnd = 0.01;                          //!Roughness length for soil [m]
constexpr const static double zsno = 0.0024;                        //!Roughness length for snow [m]
constexpr const static double csoilc = 0.004;                       //!Drag coefficient for soil under canopy [-]
constexpr const static double capr   = 0.34;                        //!Tuning factor to turn first layer T into surface T  
constexpr const static double cnfac  = 0.5;                         //!Crank Nicholson factor between 0 and 1
constexpr const static double ssi    = 0.033;                       //!Irreducible water saturation of snow
constexpr const static double wimp   = 0.05;                        //!Water impremeable if porosity less than wimp
constexpr const static double pondmx = 10.0;                        //!Ponding depth (mm)
constexpr const static double pondmx_urban = 1.0;                   //!Ponding depth for urban roof and impervious road (mm)
constexpr const static double o2_molar_const = 0.209;               //!constant atmospheric O2 molar ratio (mol/mol)
constexpr const static double maxwattabfract = 1.0;                 //!Max water table fraction for landunit that is not wet or ice
#if (defined C13)

constexpr const static double preind_atm_del13c = -6.0   //! preindustrial value for atmospheric del13C
constexpr const static double preind_atm_ratio = SHR_CONST_PDB + (preind_atm_del13c * SHR_CONST_PDB)/1000.0  //! 13C/12C
constexpr const static double c13ratio = preind_atm_ratio/(1.0+preind_atm_ratio) //! 13C/(12+13)C preind atmosphere
#endif
constexpr const static double ht_efficiency_factor = 0.75;          //!efficiency factor for urban heating (-)
constexpr const static double ac_efficiency_factor = 0.25;          //!efficiency factor for urban air conditioning (-)
constexpr const static double ht_wasteheat_factor = 1.0/ht_efficiency_factor;  //!wasteheat factor for urban heating (-)
constexpr const static double ac_wasteheat_factor = 1.0/ac_efficiency_factor;  //!wasteheat factor for urban air conditioning (-)
constexpr const static double wasteheat_limit = 100.               //!limit on wasteheat (W/m2)

/*
  !------------------------------------------------------------------
  ! Initialize water type constants
  !------------------------------------------------------------------

  ! "water" types 
  !   1     soil
  !   2     land ice (glacier)
  !   3     deep lake
  !   4     shallow lake
  !   5     wetland: swamp, marsh, etc
*/
constexpr const static std::int32_t istsoil = 1;  //!soil         "water" type
constexpr const static std::int32_t istice  = 2;  //!land ice     "water" type
constexpr const static std::int32_t istdlak = 3;  //!deep lake    "water" type
constexpr const static std::int32_t istslak = 4;  //!shallow lake "water" type
constexpr const static std::int32_t istwet  = 5;  //!wetland      "water" type
constexpr const static std::int32_t isturb  = 6;  //!urban        landunit type
#ifdef (CROP)
constexpr const static std::int32_t istcrop = 7;  //!crop         landunit type
#endif
constexpr const static std::int32_t icol_roof        = 61;
constexpr const static std::int32_t icol_sunwall     = 62;
constexpr const static std::int32_t icol_shadewall   = 63;
constexpr const static std::int32_t icol_road_imperv = 64;
constexpr const static std::int32_t icol_road_perv   = 65;

}// mod_land_varpar

}

}

#endif /*__GMS_MODULE_LAND_PARAMS_H__*/