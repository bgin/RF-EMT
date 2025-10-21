#ifndef _ATM_REFRACTIVEINDEX_H
#define _ATM_REFRACTIVEINDEX_H
/*******************************************************************************
 * ALMA - Atacama Large Millimiter Array
 * (c) Instituto de Estructura de la Materia, 2009
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 *
 * "@(#) $Id: ATMRefractiveIndex.h Exp $"
 *
 * who       when      what
 * --------  --------  ----------------------------------------------
 * pardo     24/03/09  created
 */

#ifndef __cplusplus
#error This is a C++ include file and cannot be used from plain C
#endif

#include "ATMCommon.h"
#include <complex>



ATM_NAMESPACE_BEGIN

/*! \brief RefractiveIndex allows to extract absorption and Phase coefficient(s) at
 *  a given frequency and P/T/gas densities.
 *
 *  \anchor definitions Definitions to intepret the outputs:
 *    - Phase and amplitude response of a plane wave propagating a distance \f$z\f$ at frequency \f$\nu\f$:
 *      \f$E(z)=e^{ikz(1+N)}E(0)\f$
 *    - \f$k=2\pi\nu/c\f$ is the <b>free space wave number</b>.
 *    - \f$N=(n-1)=(N_r+iN_i)\f$ is the <b>Refractivity</b> (\f$n\f$ is the <b>Refraction Index</b>). <br>
 *    - \f$N=N_{g1}+N_{g2}+N_{g3}+ ...\f$ (\f$gi\f$ are atmospheric gases or account for other mechanisms
 *      that modify the refractivity such as <br> collision-induced absorption by \f$O_2-N_2\f$, \f$N_2-N_2\f$,
 *      \f$O_2-O_2\f$, \f$O_2-H_2O\f$ or \f$N_2-H_2O\f$ pairs.
 *    - For a given gas \f$g\f$:
 *       -# \f$N_g=(N_{rg}+iN_{ig})\f$ is the <b>Refractivity</b> of that gas.
 *       -# \f$N_g/\rho_g\f$, where \f$\rho_g\f$ is the <b>number density</b> of gas \f$g\f$, is the <b>Specific Refractivity</b> of that gas.
 *       -# \f$2\pi\nu N_{rg}/c=\phi_g\f$ \f$(rad\cdot m^{-1})\f$ is the <b>Phase Dispersion Coefficient</b> of gas \f$g\f$
 *       -# \f$2\pi\nu N_{ig}/c=\kappa_g\f$ \f$(m^{-1})\f$ is the <b>Absorption Coefficient</b> of gas \f$g\f$
 *       -# \f$2\pi\nu N_{rg}/(c\rho_g)=\phi_g/\rho_{g}\f$ \f$(rad\cdot m^2)\f$ is the <b>Specific Phase Dispersion Coefficient</b> of gas \f$g\f$
 *       -# \f$2\pi\nu N_{ig}/(c\rho_g)=\kappa_g/\rho_{g}\f$ \f$(m^2)\f$ is the <b>Specific Absorption Coefficient</b> of gas \f$g\f$.
 *       -# \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ is the result of <b>getRefractivity</b> operators for gas \f$g\f$.
 *       -# \f$[2\pi\nu/(c\rho_g)]\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^2,m^2)\f$ is the result og <b>getSpecificRefractivity</b> operators for gas \f$g\f$.
 *    - Therefore in order to obtain Phase Dispersion and Absorption Coefficients:
 *	 -# Absorption Coefficient (\f$m^{-1}\f$) = Imaginary part of the result of the <b>getRefractivity</b> operator.
 *       -# Phase Dispersion Coefficient (\f$rad\cdot m^{-1}\f$) = Real part of the result of the <b>getRefractivity</b> operator.
 *       -# Specific Absorption Coefficient (\f$m^2\f$) = Imaginary part of the result of the <b>getSpecificRefractivity</b> operator.
 *       -# Specific Phase Dispersion Coefficient (\f$rad\cdot m^2\f$) = Real part of the result of the <b>getSpecificRefractivity</b> operator.
 *
 *    - <b>species \f$(g)\f$ codes: <b>
 *       -# \f$^{16}O^{16}O\f$  
 *       -# \f$^{16}O^{16}O vib\f$ 
 *       -# \f$^{16}O^{18}O\f$  
 *       -# \f$^{16}O^{17}O\f$  
 *       -# \f$CO \f$    
 *       -# \f$N_2O\f$  
 *       -# \f$NO_2\f$   
 *       -# \f$SO_2\f$   
 *       -# \f$CNTH2O \f$
 *       -# \f$CNTDRY \f$      
 *       -# \f$HH^{16}O \f$     
 *       -# \f$HH^{16}O v2 \f$
 *       -# \f$HH^{18}O \f$     
 *       -# \f$HH^{17}O  \f$  
 *       -# \f$HDO\f$   
 *       -# \f$^{16}O^{16}O^{16}O\f$     
 *       -# \f$^{16}O^{16}O^{16}O v2\f$ 
 *       -# \f$^{16}O^{16}O^{16}O v1\f$
 *       -# \f$^{16}O^{16}O^{16}O v3 \f$ 
 *       -# \f$^{16}O^{16}O^{18}O   \f$  
 *       -# \f$^{16}O^{16}O^{17}O   \f$  
 *       -# \f$^{16}O^{18}O^{16}O  \f$       
 *       -# \f$^{16}O^{17}O^{16}O \f$
 */
class RefractiveIndex
{
public:

  //@{

  /** The constructor has no arguments */
  RefractiveIndex();

  virtual ~RefractiveIndex();

  //@}

  //@{


  
  /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for \f$g=O_2\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, 
      and <b>frequency</b> in GHz. */
  std::complex<double> getRefractivity_o2(double temperature,double pressure, double wvpressure,double frequency);

  /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for \f$g=O_2\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, 
      <b>frequency</b> in GHz, <b>width</b> (channel width around <b>frequency</b>) in GHz, 
      and <b>n</b> (number of frequency points for averaging within <b>width</b>). */
  std::complex<double> getRefractivity_o2(double temperature,double pressure,double wvpressure,
				     double frequency,double width,unsigned int n);
  


  /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for \f$g=H_2O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, 
      and <b>frequency</b> in GHz. */
  std::complex<double> getRefractivity_h2o(double temperature, double pressure, double wvpressure, double frequency);

  /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for \f$g=H_2O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, 
      <b>frequency</b> in GHz, <b>width</b> (channel width around <b>frequency</b>) in GHz, 
      and <b>n</b> (number of frequency points for averaging within <b>width</b>). */
  std::complex<double> getRefractivity_h2o(double temperature,double pressure,double wvpressure,
				      double frequency,double width,unsigned int n);







  /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for \f$g=O_3\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, and <b>frequency</b> in GHz, */
  std::complex<double> getSpecificRefractivity_o3(double temperature,double pressure,double frequency);

  /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for \f$g=O_3\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, and <b>n</b> (number of frequency points for averaging within <b>width</b>) */
  std::complex<double> getSpecificRefractivity_o3(double temperature,double pressure,double frequency,
					     double width,unsigned int n);


  /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for \f$g=O_3\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      and <b>numberdensity</b> in molecules of \f$O_3\f$ m\f$^{-3}\f$ */
  inline std::complex<double> getRefractivity_o3(double temperature, double pressure, double frequency, double numberdensity)
    {return getSpecificRefractivity_o3(temperature, pressure, frequency) * numberdensity;}


  /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for \f$g=O_3\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, <b>n</b> (number of frequency points for averaging within <b>width</b>), 
      and <b>numberdensity</b> in molecules of \f$O_3\f$ m\f$^{-3}\f$ */
  std::complex<double> getRefractivity_o3(double temperature,double pressure,double frequency,
				     double width,unsigned int n,double numberdensity)
    {return getSpecificRefractivity_o3(temperature, pressure, frequency, width, n) * numberdensity;}



  /************************** 16o16o (species 1) **************************************************************/


  /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for \f$^{16}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, and <b>frequency</b> in GHz. */
  inline std::complex<double> getSpecificRefractivity_16o16o(double temperature,double pressure,double wvpressure,double frequency)
    {unsigned int species=1; return mkSpecificRefractivity(species, temperature, pressure, wvpressure, frequency);}

  /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for \f$^{16}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, and <b>n</b> (number of frequency points for averaging within <b>width</b>). */
  inline std::complex<double> getSpecificRefractivity_16o16o(double temperature,double pressure,double wvpressure,double frequency,double width,unsigned int n)
    {unsigned int species=1; return mkSpecificRefractivity(species, temperature, pressure, wvpressure, frequency, width, n);}



  /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for \f$^{16}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, <b>frequency</b> in GHz, 
      and <b>numberdensity</b> in molecules of \f$^{16}O^{16}O\f$ m\f$^{-3}\f$. */
  inline std::complex<double> getRefractivity_16o16o(double temperature,double pressure,double wvpressure,double frequency,double numberdensity)
    {return getSpecificRefractivity_16o16o(temperature, pressure, wvpressure, frequency)* numberdensity;}

  /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for \f$^{16}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, <b>n</b> (number of frequency points for averaging within <b>width</b>), 
      and <b>numberdensity</b> in molecules of \f$^{16}O^{16}O\f$ m\f$^{-3}\f$. */
  std::complex<double> getRefractivity_16o16o(double temperature,double pressure,double wvpressure,double frequency,double width,unsigned int n,double numberdensity)
    {return getSpecificRefractivity_16o16o(temperature, pressure, wvpressure, frequency, width, n)* numberdensity;}

  /*************************************************************************************************************/


  /************************** 16o16o_vib (species 2) ***********************************************************/


  /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for the first vibrationally excited state 
      of \f$^{16}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, and <b>frequency</b> in GHz. */
  inline std::complex<double> getSpecificRefractivity_16o16o_vib(double temperature,double pressure,double wvpressure,double frequency)
    {unsigned int species=2; return mkSpecificRefractivity(species, temperature, pressure, wvpressure, frequency);}

  /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for the first vibrationally excited state 
      of \f$^{16}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, and <b>n</b> (number of frequency points for averaging within <b>width</b>). */
  inline std::complex<double> getSpecificRefractivity_16o16o_vib(double temperature,double pressure,double wvpressure,double frequency,double width,unsigned int n)
    {unsigned int species=2; return mkSpecificRefractivity(species, temperature, pressure, wvpressure, frequency, width, n);}



   /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for the first vibrationally excited state 
       of \f$^{16}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, <b>frequency</b> in GHz, 
      and <b>numberdensity</b> in molecules of vibrationally excited \f$^{16}O^{16}O\f$ m\f$^{-3}\f$. */
 inline std::complex<double> getRefractivity_16o16o_vib(double temperature,double pressure,double wvpressure,double frequency,double numberdensity)
    {return getSpecificRefractivity_16o16o_vib(temperature, pressure, wvpressure, frequency)* numberdensity;}

  /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for the first vibrationally excited state 
      of \f$^{16}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, <b>n</b> (number of frequency points for averaging within <b>width</b>), 
      and <b>numberdensity</b> in molecules of vibrationally excited \f$^{16}O^{16}O\f$ m\f$^{-3}\f$. */
 std::complex<double> getRefractivity_16o16o_vib(double temperature,double pressure,double wvpressure,double frequency,double width,unsigned int n,double numberdensity)
    {return getSpecificRefractivity_16o16o_vib(temperature, pressure, wvpressure, frequency, width, n)* numberdensity;}

  /*************************************************************************************************************/


  /************************** 16o18o (species 3) ***************************************************************/


   /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for \f$^{16}O^{18}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, and <b>frequency</b> in GHz. */
 inline std::complex<double> getSpecificRefractivity_16o18o(double temperature,double pressure,double wvpressure,double frequency)
    {unsigned int species=3; return mkSpecificRefractivity(species, temperature, pressure, wvpressure, frequency);}

  /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for \f$^{16}O^{18}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, and <b>n</b> (number of frequency points for averaging within <b>width</b>). */
  inline std::complex<double> getSpecificRefractivity_16o18o(double temperature,double pressure,double wvpressure,double frequency,double width,unsigned int n)
    {unsigned int species=3; return mkSpecificRefractivity(species, temperature, pressure, wvpressure, frequency, width, n);}



   /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for \f$^{16}O^{18}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, <b>frequency</b> in GHz, 
      and <b>numberdensity</b> in molecules of \f$^{16}O^{18}O\f$ m\f$^{-3}\f$. */
 inline std::complex<double> getRefractivity_16o18o(double temperature,double pressure,double wvpressure,double frequency,double numberdensity)
    {return getSpecificRefractivity_16o18o(temperature, pressure, wvpressure, frequency)* numberdensity;}
 
  /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for \f$^{16}O^{18}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, <b>n</b> (number of frequency points for averaging within <b>width</b>), 
      and <b>numberdensity</b> in molecules of \f$^{16}O^{18}O\f$ m\f$^{-3}\f$. */
 std::complex<double> getRefractivity_16o18o(double temperature,double pressure,double wvpressure,double frequency,double width,unsigned int n,double numberdensity)
    {return getSpecificRefractivity_16o18o(temperature, pressure, wvpressure, frequency, width, n)* numberdensity;}

  /*************************************************************************************************************/

  
  /************************** 16o17o (species 4) ***************************************************************/


   /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for \f$^{16}O^{17}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, and <b>frequency</b> in GHz. */
  inline std::complex<double> getSpecificRefractivity_16o17o(double temperature,double pressure,double wvpressure,double frequency)
    {unsigned int species=4; return mkSpecificRefractivity(species, temperature, pressure, wvpressure, frequency);}
 
  /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for \f$^{16}O^{17}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, and <b>n</b> (number of frequency points for averaging within <b>width</b>). */
 inline std::complex<double> getSpecificRefractivity_16o17o(double temperature,double pressure,double wvpressure,double frequency,double width,unsigned int n)
    {unsigned int species=4; return mkSpecificRefractivity(species, temperature, pressure, wvpressure, frequency, width, n);}



   /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for \f$^{16}O^{17}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, <b>frequency</b> in GHz, 
      and <b>numberdensity</b> in molecules of \f$^{16}O^{17}O\f$ m\f$^{-3}\f$. */
  inline std::complex<double> getRefractivity_16o17o(double temperature,double pressure,double wvpressure,double frequency,double numberdensity)
    {return getSpecificRefractivity_16o17o(temperature, pressure, wvpressure, frequency)* numberdensity;}

 /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for \f$^{16}O^{17}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, <b>n</b> (number of frequency points for averaging within <b>width</b>), 
      and <b>numberdensity</b> in molecules of \f$^{16}O^{17}O\f$ m\f$^{-3}\f$. */
  std::complex<double> getRefractivity_16o17o(double temperature,double pressure,double wvpressure,double frequency,double width,unsigned int n,double numberdensity)
    {return getSpecificRefractivity_16o17o(temperature, pressure, wvpressure, frequency, width, n)* numberdensity;}

  /*************************************************************************************************************/


 
  /************************** co (species 5) *******************************************************************/


   /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for \f$CO\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, and <b>frequency</b> in GHz. */
  inline std::complex<double> getSpecificRefractivity_co(double temperature,double pressure,double frequency)
    {unsigned int species=5; return mkSpecificRefractivity(species, temperature, pressure, frequency);}
  
  /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for \f$CO\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, and <b>n</b> (number of frequency points for averaging within <b>width</b>). */
  std::complex<double> getSpecificRefractivity_co(double temperature,double pressure,double frequency,double width,unsigned int n)
    {unsigned int species=5; return mkSpecificRefractivity(species, temperature, pressure, frequency, width, n);}
  
  

    /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for \f$CO\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      and <b>numberdensity</b> in molecules of \f$CO\f$ m\f$^{-3}\f$. */
 inline std::complex<double> getRefractivity_co(double temperature,double pressure,double frequency,double numberdensity)
    {return getSpecificRefractivity_co(temperature, pressure, frequency)* numberdensity;}
  
  /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for \f$CO\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, <b>n</b> (number of frequency points for averaging within <b>width</b>), 
      and <b>numberdensity</b> in molecules of \f$CO\f$ m\f$^{-3}\f$. */
 std::complex<double> getRefractivity_co(double temperature,double pressure,double frequency,double width,unsigned int n,double numberdensity)
    {return getSpecificRefractivity_co(temperature, pressure, frequency, width, n)* numberdensity;}
  
  /*************************************************************************************************************/



  /************************** n2o (species 6) *******************************************************************/


   /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for \f$N_2O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, and <b>frequency</b> in GHz. */
  inline std::complex<double> getSpecificRefractivity_n2o(double temperature,double pressure,double frequency)
    {unsigned int species=6; return mkSpecificRefractivity(species, temperature, pressure, frequency);}
  
   /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for \f$N_2O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, and <b>n</b> (number of frequency points for averaging within <b>width</b>). */
 std::complex<double> getSpecificRefractivity_n2o(double temperature,double pressure,double frequency,double width,unsigned int n)
    {unsigned int species=6; return mkSpecificRefractivity(species, temperature, pressure, frequency, width, n);}
  
  
    /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for \f$N_2O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      and <b>numberdensity</b> in molecules of \f$N_2O\f$ m\f$^{-3}\f$. */
  inline std::complex<double> getRefractivity_n2o(double temperature,double pressure,double frequency,double numberdensity)
    {return getSpecificRefractivity_no2(temperature, pressure, frequency)* numberdensity;}
  
   /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for \f$N_2O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, <b>n</b> (number of frequency points for averaging within <b>width</b>), 
      and <b>numberdensity</b> in molecules of \f$N_2O\f$ m\f$^{-3}\f$. */
 std::complex<double> getRefractivity_n2o(double temperature,double pressure,double frequency,double width,unsigned int n,double numberdensity)
    {return getSpecificRefractivity_no2(temperature, pressure, frequency, width, n)* numberdensity;}
  
  /*************************************************************************************************************/



  /************************** no2 (species 7) ******************************************************************/


   /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for \f$NO_2\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, and <b>frequency</b> in GHz. */
  inline std::complex<double> getSpecificRefractivity_no2(double temperature,double pressure,double frequency)
    {unsigned int species=7; return mkSpecificRefractivity(species, temperature, pressure, frequency);}
  
    /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for \f$NO_2\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, and <b>n</b> (number of frequency points for averaging within <b>width</b>). */
 std::complex<double> getSpecificRefractivity_no2(double temperature,double pressure,double frequency,double width,unsigned int n)
    {unsigned int species=7; return mkSpecificRefractivity(species, temperature, pressure, frequency, width, n);}
  
  

     /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for \f$NO_2\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      and <b>numberdensity</b> in molecules of \f$NO_2\f$ m\f$^{-3}\f$. */
 inline std::complex<double> getRefractivity_no2(double temperature,double pressure,double frequency,double numberdensity)
    {return getSpecificRefractivity_no2(temperature, pressure, frequency)* numberdensity;}
  
    /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for \f$NO_2\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, <b>n</b> (number of frequency points for averaging within <b>width</b>), 
      and <b>numberdensity</b> in molecules of \f$NO_2\f$ m\f$^{-3}\f$. */
 std::complex<double> getRefractivity_no2(double temperature,double pressure,double frequency,double width,unsigned int n,double numberdensity)
    {return getSpecificRefractivity_no2(temperature, pressure, frequency, width, n)* numberdensity;}
  
  /*************************************************************************************************************/


  /************************** so2 (species 8) ******************************************************************/


   /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for \f$SO_2\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, and <b>frequency</b> in GHz. */
  inline std::complex<double> getSpecificRefractivity_so2(double temperature,double pressure,double frequency)
    {unsigned int species=8; return mkSpecificRefractivity(species, temperature, pressure, frequency);}
  
    /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for \f$SO_2\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, and <b>n</b> (number of frequency points for averaging within <b>width</b>). */
  std::complex<double> getSpecificRefractivity_so2(double temperature,double pressure,double frequency,double width,unsigned int n)
    {unsigned int species=8; return mkSpecificRefractivity(species, temperature, pressure, frequency, width, n);}
  
  
     /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for \f$SO_2\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      and <b>numberdensity</b> in molecules of \f$SO_2\f$ m\f$^{-3}\f$. */
  inline std::complex<double> getRefractivity_so2(double temperature,double pressure,double frequency,double numberdensity)
    {return getSpecificRefractivity_so2(temperature, pressure, frequency)* numberdensity;}
  
     /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for \f$SO_2\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, <b>n</b> (number of frequency points for averaging within <b>width</b>), 
      and <b>numberdensity</b> in molecules of \f$SO_2\f$ m\f$^{-3}\f$. */
 std::complex<double> getRefractivity_so2(double temperature,double pressure,double frequency,double width,unsigned int n,double numberdensity)
    {return getSpecificRefractivity_so2(temperature, pressure, frequency, width, n)* numberdensity;}
  
  /*************************************************************************************************************/

 

  /************************** cnth2o (species 9) ***************************************************************/

  /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ corresponding to the contribution of 
      the "wet" collision induced (\f$O_2-H_2O\f$ and \f$N_2-H_2O\f$) processes (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, and <b>frequency</b> in GHz. */
  inline std::complex<double> getSpecificRefractivity_cnth2o(double temperature,double pressure,double wvpressure,double frequency)
    {unsigned int species=9; return mkSpecificRefractivity(species, temperature, pressure, wvpressure, frequency);}

  /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ corresponding to the contribution of 
      the "wet" collision induced (\f$O_2-H_2O\f$ and \f$N_2-H_2O\f$) processes (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, and <b>n</b> (number of frequency points for averaging within <b>width</b>). */
  std::complex<double> getSpecificRefractivity_cnth2o(double temperature,double pressure,double wvpressure,double frequency,double width,unsigned int n)
    {unsigned int species=9; return mkSpecificRefractivity(species, temperature, pressure, wvpressure, frequency, width, n);}

  /*************************************************************************************************************/



  /************************** cntdry (species 10) **************************************************************/


  /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ corresponding to the contribution of 
      the "dry" collision induced (\f$O_2-O_2\f$, \f$N_2-N_2\f$ and \f$N_2-O_2\f$) processes (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, and <b>frequency</b> in GHz. */
  inline std::complex<double> getSpecificRefractivity_cntdry(double temperature,double pressure,double wvpressure,double frequency)
    {unsigned int species=10; return mkSpecificRefractivity(species, temperature, pressure, wvpressure, frequency);}

  /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ corresponding to the contribution of 
      the "dry" collision induced (\f$O_2-O_2\f$, \f$N_2-N_2\f$ and \f$N_2-O_2\f$) processes (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, and <b>n</b> (number of frequency points for averaging within <b>width</b>). */  
  std::complex<double> getSpecificRefractivity_cntdry(double temperature,double pressure,double wvpressure,double frequency,double width,unsigned int n)
    {unsigned int species=10; return mkSpecificRefractivity(species, temperature, pressure, wvpressure, frequency, width, n);}

  /*************************************************************************************************************/



  /************************** hh16o (species 11) ***************************************************************/

 
  /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for \f$H_2^{16}\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, and <b>frequency</b> in GHz. */
  inline std::complex<double> getSpecificRefractivity_hh16o(double temperature,double pressure,double wvpressure,double frequency)
    {unsigned int species=11; return mkSpecificRefractivity(species, temperature, pressure, wvpressure, frequency);}

  /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for \f$H_2^{16}\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, and <b>n</b> (number of frequency points for averaging within <b>width</b>). */
  std::complex<double> getSpecificRefractivity_hh16o(double temperature,double pressure,double wvpressure,double frequency,double width,unsigned int n)
    {unsigned int species=11; return mkSpecificRefractivity(species, temperature, pressure, wvpressure, frequency, width, n);}

  /*************************************************************************************************************/


  /************************** hh16o_v2 (species 12) *************************************************************/

  /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for the v2 vibrational state of \f$H_2^{16}\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, and <b>frequency</b> in GHz. */
  inline std::complex<double> getSpecificRefractivity_hh16o_v2(double temperature,double pressure,double wvpressure,double frequency)
    {unsigned int species=12; return mkSpecificRefractivity(species, temperature, pressure, wvpressure, frequency);}

  /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for the v2 vibrational state of \f$H_2^{16}\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, and <b>n</b> (number of frequency points for averaging within <b>width</b>). */
  std::complex<double> getSpecificRefractivity_hh16o_v2(double temperature,double pressure,double wvpressure,double frequency,double width,unsigned int n)
    {unsigned int species=12; return mkSpecificRefractivity(species, temperature, pressure, wvpressure, frequency, width, n);}

  /*************************************************************************************************************/


  /************************** hh18o (species 13) ***************************************************************/

  /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for \f$H_2^{18}\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, and <b>frequency</b> in GHz. */
  inline std::complex<double> getSpecificRefractivity_hh18o(double temperature,double pressure,double wvpressure,double frequency)
    {unsigned int species=13; return mkSpecificRefractivity(species, temperature, pressure, wvpressure, frequency);}

  /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for \f$H_2^{18}\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, and <b>n</b> (number of frequency points for averaging within <b>width</b>). */
  std::complex<double> getSpecificRefractivity_hh18o(double temperature,double pressure,double wvpressure,double frequency,double width,unsigned int n)
    {unsigned int species=13; return mkSpecificRefractivity(species, temperature, pressure, wvpressure, frequency, width, n);}

  
  /*************************************************************************************************************/


  /************************** hh17o (species 14) ***************************************************************/

  /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for \f$H_2^{17}\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, and <b>frequency</b> in GHz. */
  inline std::complex<double> getSpecificRefractivity_hh17o(double temperature,double pressure,double wvpressure,double frequency)
    {unsigned int species=14; return mkSpecificRefractivity(species, temperature, pressure, wvpressure, frequency);}

  /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for \f$H_2^{17}\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, and <b>n</b> (number of frequency points for averaging within <b>width</b>). */
  std::complex<double> getSpecificRefractivity_hh17o(double temperature,double pressure,double wvpressure,double frequency,double width,unsigned int n)
    {unsigned int species=14; return mkSpecificRefractivity(species, temperature, pressure, wvpressure, frequency, width, n);}

  
  /*************************************************************************************************************/


  /************************** hdo (species 15) *****************************************************************/

  /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for \f$HDO\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, and <b>frequency</b> in GHz. */
  inline std::complex<double> getSpecificRefractivity_hdo(double temperature,double pressure,double wvpressure,double frequency)
    {unsigned int species=15; return mkSpecificRefractivity(species, temperature, pressure, wvpressure, frequency);}

   /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for \f$HDO\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>wvpressure</b> (water vapor partial pressure) in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, and <b>n</b> (number of frequency points for averaging within <b>width</b>). */
 std::complex<double> getSpecificRefractivity_hdo(double temperature,double pressure,double wvpressure,double frequency,double width,unsigned int n)
    {unsigned int species=15; return mkSpecificRefractivity(species, temperature, pressure, wvpressure, frequency, width, n);}

  
  /*************************************************************************************************************/



  /************************** 16o16o16o (species 16) ************************************************************/


   /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for \f$^{16}O^{16}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, and <b>frequency</b> in GHz. */
  inline std::complex<double> getSpecificRefractivity_16o16o16o(double temperature,double pressure,double frequency)
    {unsigned int species=16; return mkSpecificRefractivity(species, temperature, pressure, frequency);}
  
     /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for \f$^{16}O^{16}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, and <b>n</b> (number of frequency points for averaging within <b>width</b>). */
 std::complex<double> getSpecificRefractivity_16o16o16o(double temperature,double pressure,double frequency,double width,unsigned int n)
    {unsigned int species=16; return mkSpecificRefractivity(species, temperature, pressure, frequency, width, n);}
  
  

     /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for \f$^{16}O^{16}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      and <b>numberdensity</b> in molecules of \f$O_3\f$ m\f$^{-3}\f$. */
 inline std::complex<double> getRefractivity_16o16o16o(double temperature,double pressure,double frequency,double numberdensity)
    {return getSpecificRefractivity_16o16o16o(temperature, pressure, frequency)* numberdensity;}
  
  
     /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for \f$^{16}O^{16}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, <b>n</b> (number of frequency points for averaging within <b>width</b>), 
      and <b>numberdensity</b> in molecules of \f$O_3\f$ m\f$^{-3}\f$. */
  std::complex<double> getRefractivity_16o16o16o(double temperature,double pressure,double frequency,double width,unsigned int n,double numberdensity)
    {return getSpecificRefractivity_16o16o16o(temperature, pressure, frequency, width, n)* numberdensity;}
  
  /*************************************************************************************************************/

 

  /************************** 16o16o16o v2 (species 17) *********************************************************/
  

   /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for vibrationally excited (v2 state)
       \f$^{16}O^{16}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, and <b>frequency</b> in GHz. */
  inline std::complex<double> getSpecificRefractivity_16o16o16o_v2(double temperature,double pressure,double frequency)
    {unsigned int species=17; return mkSpecificRefractivity(species, temperature, pressure, frequency);}
  
     /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for vibrationally excited (v2 state)
	 \f$^{16}O^{16}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, and <b>n</b> (number of frequency points for averaging within <b>width</b>). */
  std::complex<double> getSpecificRefractivity_16o16o16o_v2(double temperature,double pressure,double frequency,double width,unsigned int n)
    {unsigned int species=17; return mkSpecificRefractivity(species, temperature, pressure, frequency, width, n);}
  

     /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for vibrationally excited (v2 state)
	 \f$^{16}O^{16}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      and <b>numberdensity</b> in molecules of \f$^{16}O^{16}O^{16}O\f$  \f$v=2\f$    m\f$^{-3}\f$. */
  inline std::complex<double> getRefractivity_16o16o16o_v2(double temperature,double pressure,double frequency,double numberdensity)
    {return getSpecificRefractivity_16o16o16o_v2(temperature, pressure, frequency)* numberdensity;}

     /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for vibrationally excited (v2 state)
	 \f$^{16}O^{16}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, <b>n</b> (number of frequency points for averaging within <b>width</b>), 
      and <b>numberdensity</b> in molecules of \f$^{16}O^{16}O^{16}O\f$  \f$v=2\f$   m\f$^{-3}\f$. */
  std::complex<double> getRefractivity_16o16o16o_v2(double temperature,double pressure,double frequency,double width,unsigned int n,double numberdensity)
    {return getSpecificRefractivity_16o16o16o_v2(temperature, pressure, frequency, width, n)* numberdensity;}

  /*************************************************************************************************************/



  /************************** 16o16o16o v1 (species 18) *********************************************************/
  

   /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for vibrationally excited (v1 state)
       \f$^{16}O^{16}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, and <b>frequency</b> in GHz. */
  inline std::complex<double> getSpecificRefractivity_16o16o16o_v1(double temperature,double pressure,double frequency)
    {unsigned int species=18; return mkSpecificRefractivity(species, temperature, pressure, frequency);}

     /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for vibrationally excited (v1 state)
	 \f$^{16}O^{16}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, and <b>n</b> (number of frequency points for averaging within <b>width</b>). */
  std::complex<double> getSpecificRefractivity_16o16o16o_v1(double temperature,double pressure,double frequency,double width,unsigned int n)
    {unsigned int species=18; return mkSpecificRefractivity(species, temperature, pressure, frequency, width, n);}


     /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for vibrationally excited (v1 state)
	 \f$^{16}O^{16}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      and <b>numberdensity</b> in molecules of \f$^{16}O^{16}O^{16}O\f$  \f$v=1\f$   m\f$^{-3}\f$. */
  inline std::complex<double> getRefractivity_16o16o16o_v1(double temperature,double pressure,double frequency,double numberdensity)
    {return getSpecificRefractivity_16o16o16o_v1(temperature, pressure, frequency)* numberdensity;}

     /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for vibrationally excited (v1 state)
	 \f$^{16}O^{16}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, <b>n</b> (number of frequency points for averaging within <b>width</b>), 
      and <b>numberdensity</b> in molecules of \f$^{16}O^{16}O^{16}O\f$ \f$v=1\f$  m\f$^{-3}\f$. */
  std::complex<double> getRefractivity_16o16o16o_v1(double temperature,double pressure,double frequency,double width,unsigned int n,double numberdensity)
    {return getSpecificRefractivity_16o16o16o_v1(temperature, pressure, frequency, width, n)* numberdensity;}

  /*************************************************************************************************************/


  
  /************************** 16o16o16o v3 (species 19) *********************************************************/


   /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for vibrationally excited (v3 state)
       \f$^{16}O^{16}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, and <b>frequency</b> in GHz. */
  inline std::complex<double> getSpecificRefractivity_16o16o16o_v3(double temperature,double pressure,double frequency)
    {unsigned int species=19; return mkSpecificRefractivity(species, temperature, pressure, frequency);}

     /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for vibrationally excited (v3 state)
	 \f$^{16}O^{16}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, and <b>n</b> (number of frequency points for averaging within <b>width</b>). */
  std::complex<double> getSpecificRefractivity_16o16o16o_v3(double temperature,double pressure,double frequency,double width,unsigned int n)
    {unsigned int species=19; return mkSpecificRefractivity(species, temperature, pressure, frequency, width, n);}

     /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for vibrationally excited (v3 state)
	 \f$^{16}O^{16}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      and <b>numberdensity</b> in molecules of \f$^{16}O^{16}O^{16}O\f$ \f$v=3\f$  m\f$^{-3}\f$. */
  inline std::complex<double> getRefractivity_16o16o16o_v3(double temperature,double pressure,double frequency,double numberdensity)
    {return getSpecificRefractivity_16o16o16o_v3(temperature, pressure, frequency)* numberdensity;}

     /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for vibrationally excited (v3 state)
	 \f$^{16}O^{16}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, <b>n</b> (number of frequency points for averaging within <b>width</b>), 
      and <b>numberdensity</b> in molecules of \f$^{16}O^{16}O^{16}O\f$ \f$v=3\f$  m\f$^{-3}\f$. */
  std::complex<double> getRefractivity_16o16o16o_v3(double temperature,double pressure,double frequency,double width,unsigned int n,double numberdensity)
    {return getSpecificRefractivity_16o16o16o_v3(temperature, pressure, frequency, width, n)* numberdensity;}

  /*************************************************************************************************************/



  /************************** 16o16o18o    (species 20) *********************************************************/


   /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for 
       \f$^{16}O^{16}O^{18}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, and <b>frequency</b> in GHz. */
  inline std::complex<double> getSpecificRefractivity_16o16o18o(double temperature,double pressure,double frequency)
      {unsigned int species=20; return mkSpecificRefractivity(species, temperature, pressure, frequency);}

     /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for 
	 \f$^{16}O^{16}O^{18}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, and <b>n</b> (number of frequency points for averaging within <b>width</b>). */
  std::complex<double> getSpecificRefractivity_16o16o18o(double temperature,double pressure,double frequency,double width,unsigned int n)
    {unsigned int species=20; return mkSpecificRefractivity(species, temperature, pressure, frequency, width, n);}


     /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for 
	 \f$^{16}O^{16}O^{18}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      and <b>numberdensity</b> in molecules of \f$^{16}O^{16}O^{18}O\f$ m\f$^{-3}\f$. */
  inline std::complex<double> getRefractivity_16o16o18o(double temperature,double pressure,double frequency,double numberdensity)
    {return getSpecificRefractivity_16o16o18o(temperature, pressure, frequency)* numberdensity;}

     /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for 
	 \f$^{16}O^{16}O^{18}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, <b>n</b> (number of frequency points for averaging within <b>width</b>), 
      and <b>numberdensity</b> in molecules of \f$^{16}O^{16}O^{18}O\f$ m\f$^{-3}\f$. */
  std::complex<double> getRefractivity_16o16o18o(double temperature,double pressure,double frequency,double width,unsigned int n,double numberdensity)
    {return getSpecificRefractivity_16o16o18o(temperature, pressure, frequency, width, n)* numberdensity;}

  /*************************************************************************************************************/



  /************************** 16o16o17o    (species 21) *********************************************************/

   /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for 
       \f$^{16}O^{16}O^{17}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, and <b>frequency</b> in GHz. */
  inline std::complex<double> getSpecificRefractivity_16o16o17o(double temperature,double pressure,double frequency)
      {unsigned int species=21; return mkSpecificRefractivity(species, temperature, pressure, frequency);}

     /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for 
	 \f$^{16}O^{16}O^{17}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, and <b>n</b> (number of frequency points for averaging within <b>width</b>). */
  std::complex<double> getSpecificRefractivity_16o16o17o(double temperature,double pressure,double frequency,double width,unsigned int n)
    {unsigned int species=21; return mkSpecificRefractivity(species, temperature, pressure, frequency, width, n);}


     /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for 
	 \f$^{16}O^{16}O^{17}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      and <b>numberdensity</b> in molecules of \f$^{16}O^{16}O^{17}O\f$ m\f$^{-3}\f$. */
  inline std::complex<double> getRefractivity_16o16o17o(double temperature,double pressure,double frequency,double numberdensity)
    {return getSpecificRefractivity_16o16o17o(temperature, pressure, frequency)* numberdensity;}

     /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for 
	 \f$^{16}O^{16}O^{17}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, <b>n</b> (number of frequency points for averaging within <b>width</b>), 
      and <b>numberdensity</b> in molecules of \f$^{16}O^{16}O^{17}O\f$ m\f$^{-3}\f$. */
  std::complex<double> getRefractivity_16o16o17o(double temperature,double pressure,double frequency,double width,unsigned int n,double numberdensity)
    {return getSpecificRefractivity_16o16o17o(temperature, pressure, frequency, width, n)* numberdensity;}

  /*************************************************************************************************************/



  /************************** 16o18o16o    (species 22) *********************************************************/

   /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for 
       \f$^{16}O^{18}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, and <b>frequency</b> in GHz. */
  inline std::complex<double> getSpecificRefractivity_16o18o16o(double temperature,double pressure,double frequency)
      {unsigned int species=22; return mkSpecificRefractivity(species, temperature, pressure, frequency);}

     /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for 
	 \f$^{16}O^{18}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, and <b>n</b> (number of frequency points for averaging within <b>width</b>). */
  std::complex<double> getSpecificRefractivity_16o18o16o(double temperature,double pressure,double frequency,double width,unsigned int n)
    {unsigned int species=22; return mkSpecificRefractivity(species, temperature, pressure, frequency, width, n);}

     /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for 
	 \f$^{16}O^{18}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      and <b>numberdensity</b> in molecules of \f$^{16}O^{18}O^{16}O\f$ m\f$^{-3}\f$. */
  inline std::complex<double> getRefractivity_16o18o16o(double temperature,double pressure,double frequency,double numberdensity)
    {return getSpecificRefractivity_16o18o16o(temperature, pressure, frequency)* numberdensity;}

      /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for 
	 \f$^{16}O^{18}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, <b>n</b> (number of frequency points for averaging within <b>width</b>), 
      and <b>numberdensity</b> in molecules of \f$^{16}O^{18}O^{16}O\f$ m\f$^{-3}\f$. */
 std::complex<double> getRefractivity_16o18o16o(double temperature,double pressure,double frequency,double width,unsigned int n,double numberdensity)
    {return getSpecificRefractivity_16o18o16o(temperature, pressure, frequency, width, n)* numberdensity;}

  /*************************************************************************************************************/



  /************************** 16o17o16o    (species 23) *********************************************************/

   /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for 
       \f$^{16}O^{17}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, and <b>frequency</b> in GHz. */
  inline std::complex<double> getSpecificRefractivity_16o17o16o(double temperature,double pressure,double frequency)
      {unsigned int species=22; return mkSpecificRefractivity(species, temperature, pressure, frequency);}

     /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for 
	 \f$^{16}O^{17}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, and <b>n</b> (number of frequency points for averaging within <b>width</b>). */
  std::complex<double> getSpecificRefractivity_16o17o16o(double temperature,double pressure,double frequency,double width,unsigned int n)
    {unsigned int species=22; return mkSpecificRefractivity(species, temperature, pressure, frequency, width, n);}

     /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for 
	 \f$^{16}O^{17}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      and <b>numberdensity</b> in molecules of \f$^{16}O^{17}O^{16}O\f$ m\f$^{-3}\f$. */
  inline std::complex<double> getRefractivity_16o17o16o(double temperature,double pressure,double frequency,double numberdensity)
    {return getSpecificRefractivity_16o17o16o(temperature, pressure, frequency)* numberdensity;}

      /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for 
	 \f$^{16}O^{17}O^{16}O\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      <b>width</b> (channel width around <b>frequency</b>) in GHz, <b>n</b> (number of frequency points for averaging within <b>width</b>), 
      and <b>numberdensity</b> in molecules of \f$^{16}O^{17}O^{16}O\f$ m\f$^{-3}\f$. */
  std::complex<double> getRefractivity_16o17o16o(double temperature,double pressure,double frequency,double width,unsigned int n,double numberdensity)
    {return getSpecificRefractivity_16o17o16o(temperature, pressure, frequency, width, n)* numberdensity;}

  /*************************************************************************************************************/

     /** It returns \f$(2\pi\nu/c)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{-1},m^{-1})\f$ for 
	 \f$species\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, <b>frequency</b> in GHz, 
      and <b>numberdensity</b> in molecules of \f$species\f$ m\f$^{-3}\f$.  */
  std::complex<double> getRefractivity(unsigned int species,double temperature,double pressure,double frequency,double numberdensity)
    {return getSpecificRefractivity(species, temperature, pressure, frequency)* numberdensity;}

   /** It returns \f$(2\pi\nu/c\rho_g)\cdot(N_{rg}+iN_{ig})\f$ with units \f$(rad\cdot m^{2},m^{2})\f$ for 
       \f$species\f$ (see \ref definitions) <br>
      The parameters are <b>temperature</b> in K, <b>pressure</b> in hPa, and <b>frequency</b> in GHz.   */
  std::complex<double> getSpecificRefractivity(unsigned int species,double temperature,double pressure,double frequency)
    {return mkSpecificRefractivity(species, temperature, pressure, frequency);}



  //@}

private:

  std::complex<double> mkSpecificRefractivity(unsigned int species,    // species --> 1 to 23
					 double temperature,                 
					 double pressure,
					 double wvpressure,
					 double frequency);

  std::complex<double> mkSpecificRefractivity(unsigned int species, double temperature, double pressure, double frequency){return mkSpecificRefractivity(species, temperature, pressure, double(0.0), frequency);}

  std::complex<double> mkSpecificRefractivity(unsigned int species,    // species --> 1 to 23
					 double temperature,                 
					 double pressure,
					 double wvpressure,
					 double frequency,
                                         double width,
                                         unsigned int n);

  std::complex<double> mkSpecificRefractivity(unsigned int species, double temperature, double pressure, double frequency, double width, unsigned int n){return mkSpecificRefractivity(species, temperature, pressure, double(0.0), frequency, width, n);}

  std::complex<double> mkSpecificRefractivity_16o16o(double temperature,          /// 1
                                                double pressure,
                                                double wvpressure,
                                                double frequency);
  std::complex<double> mkSpecificRefractivity_16o16o_vib(double temperature,      /// 2
                                                    double pressure,
                                                    double wvpressure,
                                                    double frequency);
  std::complex<double> mkSpecificRefractivity_16o18o(double temperature,          /// 3
                                                double pressure,
                                                double wvpressure,
                                                double frequency);
  std::complex<double> mkSpecificRefractivity_16o17o(double temperature,          /// 4
                                                double pressure,
                                                double wvpressure,
                                                double frequency);
  std::complex<double> mkSpecificRefractivity_co(double temperature,              /// 5
                                            double pressure,
                                            double frequency);
  std::complex<double> mkSpecificRefractivity_n2o(double temperature,             /// 6
                                             double pressure,
                                             double frequency);
  std::complex<double> mkSpecificRefractivity_no2(double temperature,             /// 7
                                             double pressure,
                                             double frequency);
  std::complex<double> mkSpecificRefractivity_so2(double temperature,             /// 8
                                             double pressure,
                                             double frequency);
  std::complex<double> mkSpecificRefractivity_cnth2o(double temperature,          /// 9
                                                double pressure,
                                                double wvpressure,
                                                double frequency);
  std::complex<double> mkSpecificRefractivity_cntdry(double temperature,          /// 10
                                                double pressure,
                                                double wvpressure,
                                                double frequency);
  std::complex<double> mkSpecificRefractivity_hh16o(double temperature,           /// 11
                                               double pressure,
                                               double wvpressure,
                                               double frequency);
  std::complex<double> mkSpecificRefractivity_hh16o_v2(double temperature,        /// 12
                                                  double pressure,
                                                  double wvpressure,
                                                  double frequency);
  std::complex<double> mkSpecificRefractivity_hh18o(double temperature,           /// 13
                                               double pressure,
                                               double wvpressure,
                                               double frequency);
  std::complex<double> mkSpecificRefractivity_hh17o(double temperature,           /// 14
                                               double pressure,
                                               double wvpressure,
                                               double frequency);
  std::complex<double> mkSpecificRefractivity_hdo(double temperature,             /// 15
                                             double pressure,
                                             //double wvpressure,            // unused parameter
                                             double frequency);
  std::complex<double> mkSpecificRefractivity_16o16o16o(double temperature,       /// 16
                                                   double pressure,
                                                   double frequency);
  std::complex<double> mkSpecificRefractivity_16o16o16o_v2(double temperature,    /// 17
                                                   double pressure,
                                                   double frequency);
  std::complex<double> mkSpecificRefractivity_16o16o16o_v1(double temperature,    /// 18
                                                   double pressure,
                                                   double frequency);
  std::complex<double> mkSpecificRefractivity_16o16o16o_v3(double temperature,    /// 19
                                                   double pressure,
                                                   double frequency);
  std::complex<double> mkSpecificRefractivity_16o16o18o(double temperature,       /// 20
                                                   double pressure,
                                                   double frequency);
  std::complex<double> mkSpecificRefractivity_16o16o17o(double temperature,       /// 21
                                                   double pressure,
                                                   double frequency);
  std::complex<double> mkSpecificRefractivity_16o18o16o(double temperature,       /// 22
                                                   double pressure,
                                                   double frequency);
  std::complex<double> mkSpecificRefractivity_16o17o16o(double temperature,       /// 23
                                                   double pressure,
                                                   double frequency);

  unsigned int vpIndex(double nu);
  double linebroadening(double frequency,
                        double temperature,
                        double pressure,
                        double mmol,
                        double dv0_lines,
                        double texp_lines);
  double linebroadening_water(double frequency,
                              double temperature,
                              double pressure,
                              double wvpressure,
                              double ensanche1,
                              double ensanche2,
                              double ensanche3,
                              double ensanche4);
  double linebroadening_hh18o_hh17o(double temperature,
                                    double pressure,
                                    double ph2o,
                                    double dv0,
                                    double dvlm,
                                    double temp_exp);
  double linebroadening_o2(double frequency,
                           double temperature,
                           double pressure,
                           double ph2o,
                           double mmol,
                           double ensanche1,
                           double ensanche2);
  double interf_o2(double temperature,
                   double pressure,
                   double ensanche3,
                   double ensanche4);
  std::complex<double> lineshape(double frequency,
                            double linefreq,
                            double linebroad,
                            double interf);

}; // class RefractiveIndex

ATM_NAMESPACE_END

#endif /*!_ATM_REFRACTIVEINDEX_H*/
