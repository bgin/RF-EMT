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

#ifndef __GMS_ANTENNA_COMMON_ADT_H__
#define __GMS_ANTENNA_COMMON_ADT_H__ 201120231117


namespace file_info {

     const unsigned int GMS_ANTENNA_COMMON_ADT_MAJOR = 2;
     const unsigned int GMS_ANTENNA_COMMON_ADT_MINOR = 0;
     const unsigned int GMS_ANTENNA_COMMON_ADT_MICRO = 0;
     const unsigned int GMS_ANTENNA_COMMON_ADT_FULLVER =
       1000U*GMS_ANTENNA_COMMON_ADT_MAJOR+100U*GMS_ANTENNA_COMMON_ADT_MINOR+
       10U*GMS_ANTENNA_COMMON_ADT_MICRO;
     static const char GMS_ANTENNA_COMMON_ADT_CREATION_DATE[] = "20-11-2023 11:17 +00200 (SUN 20 NOV 2023 GMT+2)";
     static const char GMS_ANTENNA_COMMON_ADT_BUILD_DATE[]    = __DATE__; 
     static const char GMS_ANTENNA_COMMON_ADT_BUILD_TIME[]    = __TIME__;
     static const char GMS_ANTENNA_COMMON_ADT_SYNOPSIS[]      = "Antenna model common abstract data types.";

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
#include "GMS_dyn_containers.hpp"



namespace gms {


          namespace  radiolocation {
          


           
                  


              
            template<typename T>
            struct Ff239R1x_t {
                      // Psi function's (phi,theta) values
                      std::size_t        nph; // number of Psi function's phi values  -- 1st dimension
                      std::size_t        nth; // number of Psi function's theta values -- 2nd dimension
                      DC1D_t<T> f239;
           };
           
           // ! Elementary electric dipoles (radiation patterns) (2.40)
            template<typename T>
            struct  Ff240R1x_t {
                      // F(theta) values per each dipole
                      std::size_t        ndip; // number of dipoles
                      std::size_t        nth; // number of  F(theta) values -- 2nd dimension
                      DC1D_t<T> f240;
           };
              

            //! Sinusoidal current distribution (2.43)
            template<typename T>
            struct  If243C1x_t {
                      // F(theta) values per each dipole
                      std::size_t        nz; // number of  F(theta) values -- 2nd dimension
                      DC2D_t<T> f243;
           };
           
         
         //  Radiation pattern of similiar EM radiators (2.96) 
          template<typename T,
          struct  Ff296R1x_t {
                      // F(theta) values per each dipole
                      std::size_t        nth; // number of theta values (2.96)
                      std::size_t        nph; //  number of phi values (2.96) -- 2nd dimension
                      DC1D_t<T> f296;
          };     
           
          // !Linear phase error values apperture edge (2.98)
          template<typename T,
          struct  Ff298R1x_t {
                      // F(theta) values per each dipole
                      std::size_t        nph; //number of linear phase error values apperture edge (2.98)
                      DC1D_t<T> f298;
         };   
           
           
           // !Radiation pattern including linear phase error term (2.100)  
          template<typename T>
          struct  Ff2100R1x_t {
                      // F(theta) values per each dipole
                      std::size_t        nph; //number of values for radiation pattern linear phase error (2.100)
                      DC1D_t<T> f2100;                    
           };   
           
           //  !Radiation pattern including quadratic phase error term (2.102)
          template<typename T>
          struct  Ff2102R1x_t {
                      // F(theta) values per each dipole
                     std::size_t        nph; //number of values for radiation pattern quadratic phase error (2.102)
                     DC1D_t<T> f2102;
          };  
           
           // Radiation pattern cubic phase error and 
           // cosinusoidal amplitude distribution (2.107) 
          template<typename T>
          struct  Ff2107R1x_t {
                      // F(theta) values per each dipole
                      std::size_t        nph; //number of values for radiation pattern quadratic phase error
                                              // ! and cosinusoidal amplitude distribution (2.105)
                      DC1D_t<T> f2107;
          };

           
           // Average of radiation pattern of 2D (single) array (2.110)
          template<typename T>
          struct  Ff2110R1x_t {
                      // F(theta) values per each dipole
                      std::size_t        nph; // number of phi values (2.110)
                      std::size_t        nth; //  number of theta values (2.110) -- 2nd dimension
                      DC1D_t<T>          f2110;
         };    
           
           //  ! Average of radiation pattern of 1D (single) array (2.110a)
          template<typename T>
          struct  Ff2110aR1x_t {
                      // F(theta) values per each dipole
                      std::size_t        nth; //number of thet values for radiation pattern 
                      DC1D_t<T>          f2100a;
         };  
           
           // ! Power-averaged of radiation pattern of 2D (single) array (2.111)
          template<typename T>
          struct  Ff2111R1x_t {
                      // F(theta) values per each dipole
                     
                      std::size_t        nph; // number of phi values (2.111)
                      std::size_t        nth; //  number of theta values (2.111) -- 2nd dimension
                      DC1D_t<T>          f2111;
           };   
           
           
          // ! Power-average of radiation pattern of 1D (single) array (2.111a) 
          template<typename T>
          struct Ff2111aR1x_t {
          
                      std::size_t        nth;   //! number of theta values 1D array (2.111a)
                      DC1D_t<T>          f2111a;
          }; 
          
          
         // Phi values of coefficient of directional pattern (2.127)
         template<typename T>
         struct DPf2127R1x_t {
         
                      std::size_t       nph;  // ! number of phi values of coefficient of directional pattern (2.127)
                      std::size_t       nth;  // number of theta values of coefficient of directional pattern (2.127)
                      std::size_t       nae;  //  number of radiating elements or discrete antennas (2.127)
                      DC1D_t<T>         f2127;
         };
         
         
         // Values of real parts of m-th and n-th antenna impedance (2.143)
         template<typename T>
         struct RMNf2143R1x_t {
                      
                      std::size_t       nrmn; // number of real parts of m-th and n-th antenna impedance (2.143)
                      DC1D_t<T>         f2143;
         };
         
         
         // Values of imaginary parts of m-th and n-th antenna impedance (2.144)
         template<typename T>
         struct IMNf2144R1x_t {
         
                      std::size_t      nimn; // number of imaginary parts of m-th and n-th antenna impedance (2.144)
                      DC1D_t<T>        f2144;
         };
         
         
         // Values of mutual impedance of two antennas as a function of their distance (2.145)
         template<typename T>
         struct R12f2145R1x_t {
         
                      std::size_t      nr12; // number of values of mutual impedance of two antennas as a function of their distance (2.145)
                      DC1D_t<T>        f2145;
         };
         
         
         // Values of real parts of m-th and n-th antenna impedance (2.148)
         template<typename T>
         struct XMNf2148R1x_t {
                      
                      std::size_t      nmn; // number of real parts of m-th and n-th antenna impedance (2.148)
                      DC1D_t<T>        f2148;
         };
         
         
         // Theta values for complex radiating pattern (2.149)
         template<typename T>
         struct RPf2149R1x_t {
                    
                      std::size_t      nth;  // number of theta values for complex radiating pattern (2.149)
                      DC1D_t<T>        f2149;   
         };
         
         
         // Values of mutual impedance (real part) of two antennas as an 
                 //                 ! function of complex radiation pattern (2.150)
         template<typename T>
         struct R12f2150R1x_t {
                     
                      std::size_t      nr;    // ! number of values of mutual impedance (real part) of two antennas as an 
                                              // ! function of complex radiation pattern (2.150)
                      DC1D_t<T>        rf2150;
         };
         
         
         // Values of mutual impedance (imaginary part) of two antennas as an 
                              //    ! function of complex radiation pattern (2.150)
         template<typename T>
         struct X12f2150R1x_t {
         
                      std::size_t      nx;  // !number of values of mutual impedance (imaginary part) of two antennas as an 
                                            // ! function of complex radiation pattern (2.150)
                      DC1D_t<T>        xf2150;
         };
         
         
         //  The values 'height' of an antenna (EM-meaning) (2.153)
         template<typename T>
         struct HGf2135R1x_t {
         
                       std::size_t     nh;  //number of 'height' values of antenna (EM-meaning)  (2.153)
                       DC1D_t<T>       f2135;
         };
         
         
         // The values 'height' of an antenna (EM-meaning) symmetric vibrator (2.154)
         template<typename T>
         struct HSf2154R1x_t {
         
                       std::size_t     nh;    // number of 'height' values of antenna (EM-meaning) symmetric vibrator (2.154)
                       DC1D_t<T>       f2154;
         };
         
         
         //  The  area values as (function of an area) of an 
                                 // ! antenna (EM-meaning) general case (2.159)
         template<typename T>
         struct Af2159R1x_t {
               
                        std::size_t    na;    //  ! number of area values (function of an area) of an 
                                  //! antenna (EM-meaning) general case (2.159)
                        DC1D_t<T>      f2159;
         };
         
         
         // The  area values as (function of an area) of an 
                                  //! antenna (EM-meaning) a very narrow beam (2.160) 
         template<typename T>
         struct Af2160R1x_t {
         
                        std::size_t    na;     //  ! number of area values (function of an area) of an 
                                               // ! antenna (EM-meaning) a very narrow beam (2.160)
                        DC1D_t<T>      f2160;
         };
         
         
         // The  area values as (function of an area) of an 
                                 // ! antenna (EM-meaning) a sine-symmetric apperture (2.161)
         template<typename T>
         struct Af2161R1x_t {
         
                        std::size_t    na;    // ! number of area values (function of an area) of an 
                                              // ! antenna (EM-meaning) a sine-symmetric aperture (2.161)
                        DC1D_t<T>      f2161;
         };
         
         
         //! The  area values as (function of an area) of an 
                                 // ! antenna (EM-meaning) coaxial to Electric field tangent to apperture (2.162)
         template<typename T>
         struct Af2162R1x_t {
         
                        std::size_t    na;   // ! number of area values (function of an area) of an 
                                             //! antenna (EM-meaning) coaxial orientation of E-field tangent to apperture (2.162)
                        DC1D_t<T>      f2162;
         };
         
         
           //  The values of complex Fresnel coefficients, vertical polarization (2.169)
         template<typename T>
         struct Rvf2169R1x_t {
         
                         std::size_t   nr;   //number of complex Fresnel coefficients, vertical polarization (2.169)
                         DC2D_t<T>     f2169;
         };
         
         
         //  The values of complex Fresnel coefficients, horizontal polarization (2.170)
         template<typename T>
         struct Rvf2170R1x_t {
         
                         std::size_t   nr;   //number of complex Fresnel coefficients, horizontal polarization (2.169)
                         DC2D_t<T>     f2170;
         };
         
         
         // ! The values of Electrical field, 
                                  //! vertical polarization, receiption point (2.172)
         template<typename T>
         struct Evf2172R1x_t {
         
                         std::size_t   nr;   // ! number of values of Electrical field, 
                                             // ! vertical polarization, receiption point (2.172)
                         DC2D_t<T>     f2172;
         };
         
         
         //  ! The values of Electrical field, 
                                  //! horizontal polarization, receiption point (2.173)
         template<typename T>
         struct Evf2173R1x_t {
         
                         std::size_t   nr;   // ! number of values of Electrical field, 
                                             // ! horizontal polarization, receiption point (2.173)
                         DC2D_t<T>     f2173;
         };
         
         
         // ! Internal antenna noise temperature values (2.179)
         template<typename T>
         struct Tf2179R1x_t {
                          
                         std::size_t  nt; // ! number of values of Electrical field, 
                                 // ! horizontal polarization, receiption point (2.173)
                         DC1D_t<T>    f2179;
         };
         
         
         //  External antenna noise temperature values (2.180) 
         template<typename T>
         struct Tf2180R1x_t {
                       
                         std::size_t  nt; // number of external to antenna noise temperature values (2.180)
                         DC1D_t<T>    f2180;
         };
         
         
         //  The values of noise atmosphere temperature (2.181)
         template<typename T>
         struct Tf2181R1x_t {
                          
                         std::size_t  nt;
                         DC1D_t<T>    f2181;
         };
         
         
         // The values of total antenna noise temperature (2.182)
         template<typename T>
         struct Tf2182R1x_t {
                       
                         std::size_t  nt;
                         DC1D_t<T>    f2182;
         };
         
         
         // The values of optical curve length for 'n' curves (2.186)
         template<typename T>
         struct Qf2186R1x_t {
                      
                         std::size_t  nq;
                         DC1D_t<T>    f2186;
         };
         
         
         // The values of 'x,y,z' (eikonal) coordinate (2.187)
         template<typename T>
         struct L2187R1x_t {
                        // ! number of values of 'x' (eikonal) coordinate (2.187)
                        // ! number of values of 'y' (eikonal) coordinate (2.187)
                        // ! number of values of 'z' (eikonal) coordinate (2.187)
                         std::size_t nLx;
                         std::size_t nLy;
                         std::size_t nLz;
                         DC1D_t<T>   Lx;
                         DC1D_t<T>   Ly;
                         DC1D_t<T>   Lz;
         };
         // 

       } // radiolocation

} // gms














#endif /*__GMS_ANTENNA_COMMON_ADT_H__*/
