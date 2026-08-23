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

#ifndef __GMS_DIPOLE_ANT_TYPES_AVX_H__
#define __GMS_DIPOLE_ANT_TYPES_AVX_H__ 231220230916


namespace file_info {

     const unsigned int GMS_DIPOLE_ANT_TYPES_AVX_MAJOR = 1;
     const unsigned int GMS_DIPOLE_ANT_TYPES_AVX_MINOR = 0;
     const unsigned int GMS_DIPOLE_ANT_TYPES_AVX_MICRO = 0;
     const unsigned int GMS_DIPOLE_ANT_TYPES_AVX_FULLVER =
       1000U*GMS_DIPOLE_ANT_TYPES_AVX_MAJOR+100U*GMS_DIPOLE_ANT_TYPES_AVX_MINOR+
       10U*GMS_DIPOLE_ANT_TYPES_AVX_MICRO;
     const char * const GMS_DIPOLE_ANT_TYPES_AVX_CREATION_DATE = "23-12-2023 09:16 +00200 (SAT 23 DEC 2023 GMT+2)";
     const char * const GMS_DIPOLE_ANT_TYPES_AVX_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_DIPOLE_ANT_TYPES_AVX_SYNOPSIS      = "Dipole antenna model abstract data types - AVX based."

}


/*
 Purpose:
 !                        Dipole antenna abstract data types representing various characteristics
 !                        of [mainly] real Radiation Patterns of dipole antennae.
 !                        Various characteristics of different antenna types  
 !                        Based mainly on book titled (rus):          
 !                        Проектирование антенно фидерных устройств. Жук М.С. Молочков Ю.Б
*/


#include <cstdint>
#include "GMS_config"
#include "GMS_dyn_containers_avx.hpp"


namespace gms {

           namespace radiolocation {
     
     
                 //  ! First dimension -- number of dipole radiating elements in the dipole array!!
             
                 //! 'z' values of current distribution of symmetric
                                     //! dipole (3.1)
                 
                 struct Izf31C4x8x_t {
                         
                        std::size_t nIz; // ! number of 'z' values of current distribution of symmetric
                                         // ! dipole (3.1)
                        std::size_t nd;  //! number of dipoles in dipole array (3.1)
                        DC2D_ymm8c4_t   Iz;
                 };  
                 
                 
                 // ! 'z' values of current distribution of symmetric
                              // ! very thin dipole (3.4)                      
                 
                 struct Izf34C4x8x_t {
                        
                        std::size_t nIz; // ! number of 'z' values of current distribution of symmetric
                                         // ! very thin dipole (3.4)
                        std::size_t nd;  // ! number of dipoles in dipole array (3.4)
                        DC2D_ymm8c4_t  Iz;
                 };
                 
                 
                // 'theta' values of sinusoidal current distribution
                              // ! for dipole radiation pattern (3.5)   
                 
                struct Ftf35R4x8x_t {    
                   
                         std::size_t nF; //  ! number of 'theta' values of sinusoidal current distribution
                                         //  ! for dipole radiation pattern (3.5)
                         std::size_t nd; //  ! number of dipoles in dipole array (3.5)
                         DC1D_m256_t Ft;
                };
                
                
                // ! 'D' values of directed coefficients (KND, rus.)
                               //! for symmetric dipole (3.6)
                 
                struct Df36R4x8x_t {
                          
                          std::size_t nD; //! number of 'D' values of directed coefficients (KND, rus.)
                                          // ! for symmetric dipole (3.6)
                          std::size_t nd; // ! number of dipoles in dipole array (3.6)
                          DC1D_m256_t Dc;
                };
                
                
                // ! 'R' values of symmetric dipole impedance (function of length)
                              // ! (3.8)
                
                struct Rf38R4x8x_t {
                            
                          std::size_t nR; //! number of 'R' values of symmetric dipole impedance (function of length)
                                          //! (3.8)
                          std::size_t nd; //! number of dipoles in dipole array (3.8) 
                          DC1D_m256_t Rd;
                };
                
                
                // ! 'X' values of symmetric dipole reactive 
                               //! impedance (function of length) (3.9)
                
                struct Xf39R4x8x_t {
                          
                          std::size_t nX; //! number of 'X' values of symmetric dipole reactive 
                                          //! impedance (function of length) (3.9)
                          std::size_t nd; // ! number of dipoles in dipole array (3.9)
                          DC1D_m256_t   Xd;
                };
                
                
                //  ! Values for dipole impedance 'R' (function of length) (3.11)
                
                struct Rf311R4x8x_t {
                           
                          std::size_t nR; // ! number of values for dipole impedance 'R' (function of length) (3.11)
                          std::size_t nd; // ! number of dipoles in dipole array (3.11)
                         DC1D_m256_t   Rd;
                };
                
                
                // Values for dipole impedance 'X' (function of length) (3.11)
                
                struct Xf311R4x8x_t {
                          
                          std::size_t nX; // ! number of values for dipole impedance 'R' (function of length) (3.11)
                          std::size_t nd; // ! number of dipoles in dipole array (3.11)
                          DC1D_m256_t   Xd;
                };
                
                
                //  ! Values for short dipole impedance 'Z' (total value) (3.12)   
                
                struct Zf312R4x8x_t {
                          
                          std::size_t nZ; // ! number of values for short dipole impedance 'Z' (total value) (3.12)    
                          std::size_t nd; // number of dipoles in dipole array (3.12)  
                          DC1D_m256_t   Zd;
                };
                
                
                //  Values for 'wave' dipole impedance (function of length) (3.13)  
                
                struct rf313R4x8x_t {
                          
                          std::size_t nr; //! number of values for 'wave' dipole impedance (function of length) (3.13)
                          std::size_t nd; // ! number of dipoles in dipole array (3.13)
                          DC1D_m256_t   rd;
                };
                
                
                //  Values for 'R' length-invariant dipole impedance (3.14)
                
                struct Rf314R4x8x_t {
                          
                          std::size_t nR; // number of values for 'R' length-invariant dipole impedance (3.14)
                          std::size_t nd; // number of dipoles in dipole array (3.14)
                         DC1D_m256_t   R;
                };
                
                
                // Values for 'X' length-invariant dipole impedance (3.15)
                
                struct Xf315R4x8x_t {
                         
                          std::size_t nX; // number of values for 'X' length-invariant dipole impedance (3.15)
                          std::size_t nd; // number of dipoles in dipole array (3.15)
                          DC1D_m256_t   X;
                };
                
                
                // Beta ratio values (part of 3.15,3.14 formulae) (3.16)
                 
                struct Bf316R4x8x_t {
                      
                          std::size_t nB; // number of beta ratio values (part of 3.15,3.14 formulae) (3.16)
                          std::size_t nd; // number of dipoles in dipole array (3.16)
                          DC1D_m256_t   Br;   
                };
                
                
                // Values of scattering coefficient 'Gamma' (3.21)    
                
                struct Gf321R4x8x_t {
                       
                          std::size_t nG; // number of values of scattering coefficient 'Gamma' (3.21)
                          std::size_t nd; // number of dipoles in dipole array (3.21)
                          DC1D_m256_t   Gs;
                };
                
                
                // Values of radiation pattern of charged dipole (3.27)
                
                struct Ftf327R4x8x_t {
                          
                          std::size_t nFt; // number of values of radiation pattern of charged dipole (3.27)
                          std::size_t nd;  // number of dipoles in dipole array (3.27)
                          DC1D_m256_t   Ft;
                };
                
                
                // Values for 'R' active impedance of charged vibrator (3.28)
                 
                struct Rf328R4x8x_t {
                          
                          std::size_t nR; //  number of values for 'R' active impedance of charged vibrator (3.28)
                          std::size_t nd; //  number of dipoles in dipole array (3.28)
                          DC1D_m256_t   Rd; 
                };
                
                
                // Values of ingress impedance of thin biconic dipole (3.31) 
                
                struct rf331R4x8x_t {
                        
                          std::size_t nr; // number of values of ingress impedance of thin biconic dipole (3.31)
                          std::size_t nd; //  number of dipoles in dipole array (3.31)
                          DC1D_m256_t   rd;
                };
                
                
                // Values of total impedance of dipole (i.e. active,reactive) (3.33)  
                
                struct Zf333R4x8x_t {
                           
                           std::size_t nZ; // number of values of total impedance of dipole (active,reactive) (3.33)
                           std::size_t nd; //  number of dipoles in dipole array (3.33)
                           DC1D_m256_t   Zd;
                };
                
                
                //  Values of active impedance component of dipole (3.34, a part of 3.33)
                
                struct Xf334R4x8x_t {
                       
                           std::size_t nX; // number of values of reactive impedance component of dipole (3.34, a part of 3.33)
                           DC1D_m256_t   Xd;
                };
                
                
                // Values of active impedance component of dipole (3.34, a part of 3.33)
                
                struct Rf334R4x8x_t {
                           
                           std::size_t nR; // number of values of active impedance component of dipole (3.34, a part of 3.33)
                           DC1D_m256_t   Rd;
                };
                
                
                // Values of input impedance of the biconic dipole (3.35)
                
                struct Zf335R4x8x_t {
                           
                           std::size_t nZ; //  number of values of input impedance of the biconic dipole (3.35)
                           std::size_t nd; //   number of dipoles in dipole array (3.35)
                           DC1D_m256_t   Zd;
                };
                
                
                // !  ! Values horizontal-component [phi,delta, radiation pattern] 
                // ! of symmetric horizontal dipole (plane-parallel) (3.45)
                
                struct Ff345R4x8x_t {
                           
                           std::size_t np; //  number of phi values [radiation pattern] of symmetric horizontal dipole (plane-parallel) (3.45)
                           std::size_t nt; //  number of delta values [radiation pattern] of symmetric horizontal dipole (plane-parallel) (3.45)
                           std::size_t nd; //  number of dipoles in dipole array (3.45)
                           DC1D_m256_t   Ft;
                };
                
                
                // !  ! Values vertical-component [phi,delta, radiation pattern] 
                //! of symmetric horizontal dipole (plane-parallel) (3.46)
                 
                struct Ff346R4x8x_t {
                           
                           std::size_t np; // number of phi values horizontal-part [radiation pattern] of symmetric horizontal dipole (plane-parallel) (3.46)
                           std::size_t nt; // number of delta values horizontal-part [radiation pattern] of symmetric horizontal dipole (plane-parallel) (3.46)
                           std::size_t nd; // number of dipoles in dipole array (3.46)
                           DC1D_m256_t   Ft;
                };
                
                
                // Values of vertical-component of radiation pattern [ideally conducting earth] (3.49)
                
                struct Ff34950R4x8x_t {
                            
                            std::size_t nv; // number of values of vertical-component of radiation pattern [ideally conducting earth] (3.49)
                            std::size_t nh; // number of values of horizontal-component of radiation pattern [ideally conducting earth] (3.50)
                            std::size_t nd; // number of dipoles in dipole array (3.49,3.50)
                            DC1D_m256_t   Fv;
                            DC1D_m256_t   Fh;
                };
                
                
                // Values of vertical-component of radiation pattern of vertical symmetric dipole (3.52)
                
                struct Ff352R4x8x_t {
                            
                            std::size_t nF; // number of values of vertical-component of radiation pattern of vertical symmetric dipole (3.52)
                            std::size_t nd; // number of dipoles in dipole array (3.52)
                            DC1D_m256_t   Ft;
                };
                
                
                // Values of vertical-component of radiation pattern for assymetrical cylindrical dipole (3.66)
                
                struct Ff366R4x8x_t {
                            
                            std::size_t nF; // number of values of vertical-component of radiation pattern for assymetrical cylindrical dipole (3.66)
                            std::size_t nd; // number of dipoles in dipole array (3.66)
                            DC1D_m256_t   Ft;
                };
                
                
                // Azimuthal 'phi' and 'theta' coordinate values of radiation pattern for V-antenna (3.97)
                
                struct Ff397R4x8x_t {
                            
                            std::size_t np; // number of azimuthal 'phi' coordinate values of radiation pattern for V-antenna (3.97)
                            std::size_t nt; // number of azimuthal 'theta' coordinate values of radiation pattern for V-antenna (3.97)
                            std::size_t nd; // number of dipoles in dipole array (3.97)
                            DC1D_m256_t   Fpt;
                };
                
                
                // Meridional 'phi' and 'theta' coordinate values of radiation pattern for V-antenna (3.98)
                
                struct Ff398R4x8x_t {
                            
                            std::size_t np; // number of meridional 'phi' coordinate values of radiation pattern for V-antenna (3.98)
                            std::size_t nt; // number of meridional 'theta' coordinate values of radiation pattern for V-antenna (3.98)
                            std::size_t nd; // number of dipoles in dipole array (3.98)
                            DC1D_m256_t   Fpt;
                };
                
                
                // Impedance values of V-antenna located in the free space (3.99)
                
                struct Rf399R4x8x_t {
                            
                            std::size_t nR; // Impedance values of V-antenna located in the free space (3.99)
                            std::size_t nd; // number of dipoles in dipole array (3.99)
                            DC1D_m256_t   Rd;
                };
                
                
                // Phi and theta values of radiation pattern for V-antenna (Pistolkors antenna) (3.100) 
                
                struct Ff3100R4x8x_t {
                            
                            std::size_t np; // number of phi values of radiation pattern for V-antenna (Pistolkors antenna) (3.100)
                            std::size_t nt; // number of theta values of radiation pattern for V-antenna (Pistolkors antenna) (3.100)
                            std::size_t nd; // number of dipoles in dipole array (3.100)
                            DC1D_m256_t   Fpt;
                };
                
                
                // Phi values [horizontal plane] of radiation pattern 
                             //  ! for V-antenna (Pistolkors antenna) (3.101) 
                
                struct Ff3101R4x8x_t {
                            
                            std::size_t np; // ! number of phi values [horizontal plane] of radiation pattern 
                                            // ! for V-antenna (Pistolkors antenna) (3.101)
                            std::size_t nd; // number of dipoles in dipole array (3.101) 
                            DC1D_m256_t   Fp;
                };
                
                
                // Values (W/cm) of critical voltage for dipole antenna (3.104)
                
                struct Ef3104R4x8x_t {
                            
                            std::size_t nE; // number of values (W/cm) of critical voltage for dipole antenna (3.104)
                            DC1D_m256_t   E;
                };
                
                
                // Values for mutual impedance of two dipole antennae (3.110)
                
                struct Zf3110R4x8x_t {
                            
                            std::size_t nZ; // number of values for mutual impedance of two dipole antennae (3.110)
                            DC1D_m256_t   Z;
                };
                
                
                // Phi values of radiation pattern for two-dipole antenna (horizontal-plane) (3.122)
                
                struct Ff3122R4x8x_t {
                             
                            std::size_t np; // number of phi values of radiation pattern for two-dipole antenna (horizontal plane)(3.122)
                            DC1D_m256_t   Fp;
                };
                
                
                // Phi values of radiation pattern for two-dipole antenna 
                              // ! plane-perpendicular to antennae axis(3.123)
                
                struct Ff3123R4x8x_t {
                            
                            std::size_t np; // ! number of phi values of radiation pattern for two-dipole antenna 
                                            // ! plane-perpendicular to antennae axis(3.123)
                            DC1D_m256_t   Fp;
                };
                
                
                // ! Theta values of radiation pattern for 2D array of dipoles (horizontal plane) (3.126)
                
                struct Ff3126R4x8x_t {
                            
                            std::size_t nt; //  number of theta values of radiation pattern for 2D array of dipoles (horizontal plane) (3.126)
                            DC1D_m256_t   Ft;
                };
                
                
                 // ! Theta values of radiation pattern for 2D array of dipoles (vertical plane) (3.127)
                
                struct Ff3127R4x8x_t {
                            
                            std::size_t nt; //  number of theta values of radiation pattern for 2D array of dipoles (vertical plane) (3.127)
                            DC1D_m256_t   Ft;
                };
                
                
                //  Theta values (horizontal-plane) of radiation pattern for aperiodic reflector (3.136)
                
                struct Ff3136R4x8x_t {
                            
                            std::size_t nt; // number of theta values (horizontal-plane) of radiation pattern for aperiodic reflector (3.136)
                            DC1D_m256_t   Ft;
                };
                
                
                // Delta values (vertical-plane) of radiation pattern for aperiodic reflector (3.137)
                
                struct Ff3137R4x8x_t {
                         
                            std::size_t nd; // number of delta values (vertical-plane) of radiation pattern for aperiodic reflector (3.137)
                            DC1D_m256_t   Fd;  
                           
                };
      }



}


















#endif /*__GMS_DIPOLE_ANT_TYPES_AVX_H__*/
