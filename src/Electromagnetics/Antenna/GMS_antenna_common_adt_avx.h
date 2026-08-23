
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

#ifndef __GMS_ANTENNA_COMMON_ADT_AVX_H__
#define __GMS_ANTENNA_COMMON_ADT_AVX_H__ 171220230546



namespace file_info {

     const unsigned int GMS_ANTENNA_COMMON_ADT_AVX_MAJOR = 1;
     const unsigned int GMS_ANTENNA_COMMON_ADT_AVX_MINOR = 0;
     const unsigned int GMS_ANTENNA_COMMON_ADT_AVX_MICRO = 0;
     const unsigned int GMS_ANTENNA_COMMON_ADT_AVX_FULLVER =
       1000U*GMS_ANTENNA_COMMON_ADT_AVX_MAJOR+100U*GMS_ANTENNA_COMMON_ADT_AVX_MINOR+
       10U*GMS_ANTENNA_COMMON_ADT_AVX_MICRO;
     static const char GMS_ANTENNA_COMMON_ADT_AVX_CREATION_DATE[] = "17-12-2023 05:46 +00200 (SUN 17 DEC 2023 GMT+2)";
     static const char GMS_ANTENNA_COMMON_ADT_AVX_BUILD_DATE[]    = __DATE__;
     static const char GMS_ANTENNA_COMMON_ADT_AVX_BUILD_TIME[]    = __TIME__;
     static const char GMS_ANTENNA_COMMON_ADT_AVX_SYNOPSIS[]      = "Antenna model common abstract data types -- AVX based.";

}


/*
 Purpose:
 !                        Derived data types for 'antenna_sensor' module implementation.
 !                        Various characteristics of different antenna types  
 !                        Based mainly on book titled (rus):          
 !                        Проектирование антенно фидерных устройств. Жук М.С. Молочков Ю.Б
*/


#include <cstdint>
#include "GMS_config"
#include "GMS_dyn_containers_avx.hpp"



namespace gms {


          namespace  radiolocation {
          


           
                  


              
           
            struct Ff239R4x8x_t {
                      // Psi function's (phi,theta) values
                      std::size_t        nph; // number of Psi function's phi values  -- 1st dimension
                      std::size_t        nth; // number of Psi function's theta values -- 2nd dimension
                      DC1D_m256_t f239;
           };
           
           // ! Elementary electric dipoles (radiation patterns) (2.40)
           
            struct  Ff240R4x8x_t {
                      // F(theta) values per each dipole
                      std::size_t        ndip; // number of dipoles
                      std::size_t        nth; // number of  F(theta) values -- 2nd dimension
                      DC1D_m256_t f240;
           };
              

            //! Sinusoidal current distribution (2.43)
           
            struct  If243C4x8x_t {
                      // F(theta) values per each dipole
                      std::size_t        nz; // number of  F(theta) values -- 2nd dimension
                      DC1D_xmm4c4_t f243;
           };
           
         
         //  Radiation pattern of similiar EM radiators (2.96) 
          
          struct  Ff296R4x8x_t {
                      // F(theta) values per each dipole
                      std::size_t        nth; // number of theta values (2.96)
                      std::size_t        nph; //  number of phi values (2.96) -- 2nd dimension
                      DC1D_m256_t f296;
          };     
           
          // !Linear phase error values apperture edge (2.98)
         
          struct  Ff298R4x8x_t {
                      // F(theta) values per each dipole
                      std::size_t        nph; //number of linear phase error values apperture edge (2.98)
                      DC1D_m256_t f298;
         };  
           
           
           // !Radiation pattern including linear phase error term (2.100)  
          
          struct  Ff2100R4x8x_t {
                      // F(theta) values per each dipole
                      std::size_t        nph; //number of values for radiation pattern linear phase error (2.100)
                      DC1D_m256_t f2100;                    
           };   
           
           //  !Radiation pattern including quadratic phase error term (2.102)
          
          struct  Ff2102R4x8x_t {
                      // F(theta) values per each dipole
                     std::size_t        nph; //number of values for radiation pattern quadratic phase error (2.102)
                     DC1D_m256_t f2102;
          };  
           
           // Radiation pattern cubic phase error and 
           // cosinusoidal amplitude distribution (2.107) 
          
          struct  Ff2107R4x8x_t {
                      // F(theta) values per each dipole
                      std::size_t        nph; //number of values for radiation pattern quadratic phase error
                                              // ! and cosinusoidal amplitude distribution (2.105)
                      DC1D_m256_t f2107;
          };

           
           // Average of radiation pattern of 2D (single) array (2.110)
          
          struct  Ff2110R4x8x_t {
                      // F(theta) values per each dipole
                      std::size_t        nph; // number of phi values (2.110)
                      std::size_t        nth; //  number of theta values (2.110) -- 2nd dimension
                      DC1D_m256_t          f2110;
         };    
           
           //  ! Average of radiation pattern of 1D (single) array (2.110a)
          
          struct  Ff2110aR4x8x_t {
                      // F(theta) values per each dipole
                      std::size_t        nth; //number of thet values for radiation pattern 
                      DC1D_m256_t          f2100a;
         };  
           
           // ! Power-averaged of radiation pattern of 2D (single) array (2.111)
          
          struct  Ff2111R4x8x_t {
                      // F(theta) values per each dipole
                     
                      std::size_t        nph; // number of phi values (2.111)
                      std::size_t        nth; //  number of theta values (2.111) -- 2nd dimension
                      DC1D_m256_t          f2111;
           };   
           
           
          // ! Power-average of radiation pattern of 1D (single) array (2.111a) 
          
          struct Ff2111aR4x8x_t {
          
                      std::size_t        nth;   //! number of theta values 1D array (2.111a)
                      DC1D_m256_t          f2111a;
          }; 
          
          
         // Phi values of coefficient of directional pattern (2.127)
         
         struct DPf2127R4x8x_t {
         
                      std::size_t       nph;  // ! number of phi values of coefficient of directional pattern (2.127)
                      std::size_t       nth;  // number of theta values of coefficient of directional pattern (2.127)
                      std::size_t       nae;  //  number of radiating elements or discrete antennas (2.127)
                      DC1D_m256_t         f2127;
         };
         
         
         // Values of real parts of m-th and n-th antenna impedance (2.143)
         
         struct RMNf2143R4x8x_t {
                      
                      std::size_t       nrmn; // number of real parts of m-th and n-th antenna impedance (2.143)
                      DC1D_m256_t         f2143;
         };
         
         
         // Values of imaginary parts of m-th and n-th antenna impedance (2.144)
         
         struct IMNf2144R4x8x_t {
         
                      std::size_t      nimn; // number of imaginary parts of m-th and n-th antenna impedance (2.144)
                      DC1D_m256_t        f2144;
         };
         
         
         // Values of mutual impedance of two antennas as a function of their distance (2.145)
         
         struct R12f2145R4x8x_t {
         
                      std::size_t      nr12; // number of values of mutual impedance of two antennas as a function of their distance (2.145)
                      DC1D_m256_t        f2145;
         };
         
         
         // Values of real parts of m-th and n-th antenna impedance (2.148)
         
         struct XMNf2148R4x8x_t {
                      
                      std::size_t      nmn; // number of real parts of m-th and n-th antenna impedance (2.148)
                      DC1D_m256_t        f2148;
         };
         
         
         // Theta values for complex radiating pattern (2.149)
         
         struct RPf2149R4x2x_t {
                    
                      std::size_t      nth;  // number of theta values for complex radiating pattern (2.149)
                      DC1D_m256_t        f2149;   
         };
         
         
         // Values of mutual impedance (real part) of two antennas as an 
                 //                 ! function of complex radiation pattern (2.150)
         
         struct R12f2150R4x2x_t {
                     
                      std::size_t      nr;    // ! number of values of mutual impedance (real part) of two antennas as an 
                                              // ! function of complex radiation pattern (2.150)
                      DC1D_m256_t        rf2150;
         };
         
         
         // Values of mutual impedance (imaginary part) of two antennas as an 
                              //    ! function of complex radiation pattern (2.150)
         
         struct X12f2150R4x8x_t {
         
                      std::size_t      nx;  // !number of values of mutual impedance (imaginary part) of two antennas as an 
                                            // ! function of complex radiation pattern (2.150)
                      DC1D_m256_t        xf2150;
         };
         
         
         //  The values 'height' of an antenna (EM-meaning) (2.153)
         
         struct HGf2135R4x8x_t {
         
                       std::size_t     nh;  //number of 'height' values of antenna (EM-meaning)  (2.153)
                       DC1D_m256_t       f2135;
         };
         
         
         // The values 'height' of an antenna (EM-meaning) symmetric vibrator (2.154)
         
         struct HSf2154R4x8x_t {
         
                       std::size_t     nh;    // number of 'height' values of antenna (EM-meaning) symmetric vibrator (2.154)
                       DC1D_m256_t       f2154;
         };
         
         
         //  The  area values as (function of an area) of an 
                                 // ! antenna (EM-meaning) general case (2.159)
         
         struct Af2159R4x8x_t {
               
                        std::size_t    na;    //  ! number of area values (function of an area) of an 
                                  //! antenna (EM-meaning) general case (2.159)
                        DC1D_m256_t      f2159;
         };
         
         
         // The  area values as (function of an area) of an 
                                  //! antenna (EM-meaning) a very narrow beam (2.160) 
         
         struct Af2160R4x8x_t {
         
                        std::size_t    na;     //  ! number of area values (function of an area) of an 
                                               // ! antenna (EM-meaning) a very narrow beam (2.160)
                        DC1D_m256_t      f2160;
         };
         
         
         // The  area values as (function of an area) of an 
                                 // ! antenna (EM-meaning) a sine-symmetric apperture (2.161)
         
         struct Af2161R4x8x_t {
         
                        std::size_t    na;    // ! number of area values (function of an area) of an 
                                              // ! antenna (EM-meaning) a sine-symmetric aperture (2.161)
                        DC1D_m256_t      f2161;
         };
         
         
         //! The  area values as (function of an area) of an 
                                 // ! antenna (EM-meaning) coaxial to Electric field tangent to apperture (2.162)
         
         struct Af2162R4x8x_t {
         
                        std::size_t    na;   // ! number of area values (function of an area) of an 
                                             //! antenna (EM-meaning) coaxial orientation of E-field tangent to apperture (2.162)
                        DC1D_m256_t      f2162;
         };
         
         
           //  The values of complex Fresnel coefficients, vertical polarization (2.169)
         
         struct Rvf2169R4x8x_t {
         
                         std::size_t   nr;   //number of complex Fresnel coefficients, vertical polarization (2.169)
                         DC2D_ymm8c4_t f2169;
         };
         
         
         //  The values of complex Fresnel coefficients, horizontal polarization (2.170)
         
         struct Rvf2170R4x8x_t {
         
                         std::size_t   nr;   //number of complex Fresnel coefficients, horizontal polarization (2.169)
                         DC2D_ymm8c4_t f2170;
         };
         
         
         // ! The values of Electrical field, 
                                  //! vertical polarization, receiption point (2.172)
         
         struct Evf2172R4x8x_t {
         
                         std::size_t   nr;   // ! number of values of Electrical field, 
                                             // ! vertical polarization, receiption point (2.172)
                         DC2D_ymm8c4_t  f2172;
         };
         
         
         //  ! The values of Electrical field, 
                                  //! horizontal polarization, receiption point (2.173)
         
         struct Evf2173R4x8x_t {
         
                         std::size_t   nr;   // ! number of values of Electrical field, 
                                             // ! horizontal polarization, receiption point (2.173)
                         DC2D_ymm8c4_t f2173;
         };
         
         
         // ! Internal antenna noise temperature values (2.179)
         
         struct Tf2179R4x8x_t {
                          
                         std::size_t  nt; // ! number of values of Electrical field, 
                                 // ! horizontal polarization, receiption point (2.173)
                         DC1D_m256_t    f2179;
         };
         
         
         //  External antenna noise temperature values (2.180) 
         
         struct Tf2180R4x8x_t {
                       
                         std::size_t  nt; // number of external to antenna noise temperature values (2.180)
                         DC1D_m256_t    f2180;
         };
         
         
         //  The values of noise atmosphere temperature (2.181)
         
         struct Tf2181R4x8x_t {
                          
                         std::size_t  nt;
                         DC1D_m256_t    f2181;
         };
         
         
         // The values of total antenna noise temperature (2.182)
         
         struct Tf2182R4x8x_t {
                       
                         std::size_t  nt;
                         DC1D_m256_t    f2182;
         };
         
         
         // The values of optical curve length for 'n' curves (2.186)
         
         struct Qf2186R4x8x_t {
                      
                         std::size_t  nq;
                         DC1D_m256_t    f2186;
         };
         
         
         // The values of 'x,y,z' (eikonal) coordinate (2.187)
         
         struct L2187R4x8x_t {
                        // ! number of values of 'x' (eikonal) coordinate (2.187)
                        // ! number of values of 'y' (eikonal) coordinate (2.187)
                        // ! number of values of 'z' (eikonal) coordinate (2.187)
                         std::size_t nLx;
                         std::size_t nLy;
                         std::size_t nLz;
                         DC1D_m256_t   Lx;
                         DC1D_m256_t   Ly;
                         DC1D_m256_t   Lz;
         };
         // 

       } // radiolocation

} // gms














#endif /*__GMS_ANTENNA_COMMON_ADT_AVX_H__*/
