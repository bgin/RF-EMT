#ifndef _ATM_PROFILE_H
#define _ATM_PROFILE_H
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
 * "@(#) $Id: ATMProfile.h Exp $"
 *
 * who       when      what
 * --------  --------  ----------------------------------------------
 * pardo     24/03/09  created
 */

#ifndef __cplusplus
#error This is a C++ include file and cannot be used from plain C
#endif

#include "ATMCommon.h"
#include "ATMHumidity.h"
#include "ATMLength.h"
#include "ATMMassDensity.h"
#include "ATMNumberDensity.h"
#include "ATMPressure.h"
#include "ATMTemperature.h"
#include "ATMEnumerations.h"

#include <string>
#include <vector>

using std::string;
using std::vector;


ATM_NAMESPACE_BEGIN

/** \brief Class for an atmospheric profile object.
 *
 *         An  atmospheric profile is composed of 4 quantities as a function of altitude z:
 *         - the layer thickness
 *         - the pressure P
 *         - the temperature T and
 *         - the gas densities for H2O, O3, CO, N2O, NO2, SO2.<br>
 *
 *   This object is needed  for computing the absorption and
 *   phase coefficients, as well as for performing radiative transfer
 *   calculations (only layer thickness/T are needed).
 *
 *   This class builds an atmospheric profile that can be used to
 *   calculate absorption and phase coefficients, as well as to
 *   perform forward and/or retrieval radiative transfer calculations.
 *   It is composed of a set of parameters needed to build a layer thickness/P/T/gas densities
 *   densities profile from simple parameters currently available at observatories
 *   (from weather stations for example) using functions from
 *   the ATM library. The set of input parameters consists of the
 *   pressure P, the temperature T and the relative humidity at the ground, the altitude of
 *   the site, the tropospheric temperature lapse rate,... The profile is built as:
 *   thickness of the considered atmospheric layers above the site, and mean
 *   P,T,H2O,O3,CO,N2O, NO2, SO2 in
 *   them. The total number of atmospheric layers in the particular profile
 *   is also available (a negative value indicates an error).<br>
 *   The zenith column  of water vapor can be calculated
 *   by simply integrating the H2O profile. \f[ column\;\;H_{2}O\;\; (mm) =\sum_{layers} \Delta z \cdot [ H_{2}O ] \f]
 */
class AtmProfile
{
public:
  /** A constructor of an empty profile with n layers, that can be filled up later. */
  AtmProfile(unsigned int n);

  /** A long constructor of the atmospheric profile from the basic set of parameters described above.
      Please note that this constructor assumes that the &altitude of the antenna is the SAME of the
      weather station that provides: &groundPressure, &groundTemperature, and &relativeHumidity */
  AtmProfile(const Length &altitude,
             const Pressure &groundPressure,
             const Temperature &groundTemperature,
             double tropoLapseRate,
             const Humidity &relativeHumidity,
             const Length &wvScaleHeight,
             const Pressure &pressureStep,
             double pressureStepFactor,
             const Length &topAtmProfile,
             unsigned int atmType); //Atmospheretype atmType);


  /** A long constructor of the atmospheric profile from the basic set of parameters described above which,
      in addition, includes user-defined temperature profile  */
  AtmProfile(const Length &altitude,
             const Pressure &groundPressure,
             const Temperature &groundTemperature,
             double tropoLapseRate,
             const Humidity &relativeHumidity,
             const Length &wvScaleHeight,
             const Pressure &pressureStep,
             double pressureStepFactor,
             const Length &topAtmProfile,
             unsigned int atmType,
	     const vector<Length> &v_layerBoundaries,
	     const vector<Temperature> &v_layerTemperature);



  /** A long constructor of the atmospheric profile from the basic set of parameters described above. */
  /*
   AtmProfile(Length altitude, Pressure    groundPressure,
   Temperature groundTemperature, double  tropoLapseRate,
   Humidity relativeHumidity, Length wvScaleHeight,
   Pressure pressureStep, double pressureStepFactor,
   Length topAtmProfile, int typeAtm);

   */

  /** A short constructor of the atmospheric profile. With respect to the long constructor, this one considers
   the following (recomended) default values: pressureStep = 10 mb,
   pressureStepFactor = 1.2, topAtmProfile = 48.0 km */
  AtmProfile(const Length &altitude,
             const Pressure &groundPressure,
             const Temperature &groundTemperature,
             double tropoLapseRate,
             const Humidity &relativeHumidity,
             const Length &wvScaleHeight,
             unsigned int atmType); //Atmospheretype atmType);

  /*AtmProfile(const Length &altitude,
             const Pressure &groundPressure,
             const Temperature &groundTemperature,
	     const Humidity &relativeHumidity,
	     const vector<Length> &v_layerTop,
	     const vector<Pressure> &v_layerTopPressure,
	     const vector<Temperature> &v_layerTopTemperature,
	     const vector<Humidity> &v_layerTopHumidity); */

  /** The user provides his own atmospheric profile (basic one: four vectors for layer thickness in m,
   average pressure in each layer in mb, average temperature in each layer in K, and average water
   vapor density in each layer in kg/m**3). QUESTION: SHOULD O3, CO, N2O, NO2, SO2 PROFILES BE FILLED UP
   INTERNALLY FROM A STANDARD ATMOSPHERE OR LEFT ITS ABUNDANCES EQUAL TO ZERO ? */
  AtmProfile(const Length &altitude,
	     const vector<Length> &v_layerThickness,
             const vector<Pressure> &v_layerPressure,
             const vector<Temperature> &v_layerTemperature,
             const vector<MassDensity> &v_layerWaterVapor);
  AtmProfile(const vector<Length> &v_layerBoundaries,
             const vector<Pressure> &v_layerPressure,
             const vector<Temperature> &v_layerTemperature,
             const vector<MassDensity> &v_layerWaterVapor);

  AtmProfile(const Length &altitude,
	     const vector<Length> &v_layerThickness,
             const vector<Pressure> &v_layerPressure,
             const vector<Temperature> &v_layerTemperature,
             const vector<NumberDensity> &v_layerWaterVapor);
  AtmProfile(const vector<Length> &v_layerBoundaries,
             const vector<Pressure> &v_layerPressure,
             const vector<Temperature> &v_layerTemperature,
             const vector<NumberDensity> &v_layerWaterVapor);

  /** The user provides his own atmospheric profile (in this case five vectors for layer thickness in m,
   average pressure in each layer in mb, average temperature in each layer in K, average water vapor
   density in each layer in kg/m**3, and average ozone number density in each layer in molecules/m**3)
   QUESTION: SHOULD CO and N2O, NO2, SO2 PROFILES BE FILLED UP
   INTERNALLY FROM A STANDARD ATMOSPHERE OR LEFT ITS ABUNDANCES EQUAL TO ZERO ? */

  AtmProfile(const Length &altitude,
	     const vector<Length> &v_layerThickness,
             const vector<Pressure> &v_layerPressure,
             const vector<Temperature> &v_layerTemperature,
             const vector<MassDensity> &v_layerWaterVapor,
             const vector<NumberDensity> &v_layerO3);
  AtmProfile(const Length &altitude,
	     const vector<Length> &v_layerThickness,
             const vector<Pressure> &v_layerPressure,
             const vector<Temperature> &v_layerTemperature,
             const vector<NumberDensity> &v_layerWaterVapor,
             const vector<NumberDensity> &v_layerO3);

  /** The user provides his own atmospheric profile (in this case seven vectors for layer thickness in m,
   average pressure in each layer in mb, average temperature in each layer in K, average water vapor
   density in each layer in kg/m**3, average ozone number density in each layer in molecules/m**3,
   average CO number density in each layer in molecules/m**3, average N2O, NO2, SO2 number density in each
   layer in molecules/m**3) */
  AtmProfile(const Length &altitude,
	     const vector<Length> &v_layerThickness,
             const vector<Pressure> &v_layerPressure,
             const vector<Temperature> &v_layerTemperature,
             const vector<MassDensity> &v_layerWaterVapor,
             const vector<NumberDensity> &v_layerO3,
             const vector<NumberDensity> &v_layerCO,
             const vector<NumberDensity> &v_layerN2O,
             const vector<NumberDensity> &v_layerNO2,
             const vector<NumberDensity> &v_layerSO2);

  AtmProfile(const Length &altitude,
	     const vector<Length> &v_layerThickness,
             const vector<Pressure> &v_layerPressure,
             const vector<Temperature> &v_layerTemperature,
             const vector<NumberDensity> &v_layerWaterVapor,
             const vector<NumberDensity> &v_layerO3,
             const vector<NumberDensity> &v_layerCO,
             const vector<NumberDensity> &v_layerN2O,
             const vector<NumberDensity> &v_layerNO2,
             const vector<NumberDensity> &v_layerSO2);

  AtmProfile(const AtmProfile &a); // copy constructor

  virtual ~AtmProfile() {}

  /** Setter to update the AtmProfile if some basic atmospheric parameter has changed.
   * @pre   an atmospheric profile already exists
   * @param altitude          the new altitude, a Length
   * @param groundPressure    the Pressure at the ground level
   * @param groundTemperature the Temperature at the ground level
   * @param tropoLapseRate    the tropospheric lapse rate
   * @param relativeHumidity  the relative Humidity
   * @param wvScaleHeight     the scale height of the water vapor, a Length
   * @return true if the atmospheric profile has been update, else false because the basic parameters have not changed
   * @post   the atmospheric profile has been updated, only if at least one of the basic parameters has changed
   *
   * \note  there is an overriding on this method in the sub-class AbsorptionPhaseProfile which in turn
   *        has an  overriding in its WaterVaporRetrieval sub-class. Hence this method must not be overloaded
   *        in this AtmProfile class.
   */
  bool setBasicAtmosphericParameters(const Length &altitude,
                                     const Pressure &groundPressure,
                                     const Temperature &groundTemperature,
                                     double tropoLapseRate,
                                     const Humidity &relativeHumidity,
                                     const Length &wvScaleHeight);

  /** Accessor to the type of current atmosphere **/
  string getAtmosphereType() const;

  /** Accessor to the type of atmosphere specified by the number**/
  static string getAtmosphereType(unsigned int typeAtm);

  /** Accessor to the current Ground Temperature used in the object */
  Temperature getGroundTemperature() const { return groundTemperature_; }

  /** Accessor to the current Tropospheric Lapse Rate used in the object (temperature/length units) */
  double getTropoLapseRate() const { return tropoLapseRate_; }

  /** Accessor to the current Ground Pressure used in the object (pressure units) */
  Pressure getGroundPressure() const { return groundPressure_; }

  /** Accessor to the current Tropopause Temperature used in the object */
  Temperature getTropopauseTemperature() const { return tropoTemperature_; }

  /** Accessor to the current Ground Relative Humidity in the object (humidity units) */
  Humidity getRelativeHumidity() const { return relativeHumidity_; }

  /** Accessor to the current Water Vapor Scale Height in the object (length units) */
  Length getWvScaleHeight() const { return wvScaleHeight_; }

  /** Accessor to the current Primary Pressure Step in the object. The Primary Pressure Step
   (pressure units) is used to define the thickness of the first layer in the profile.
   Pressure difference between the boundaries of first layer will be equal to the
   Primary Pressure Step. */
  Pressure getPressureStep() const { return pressureStep_; }

  /** Accessor to the current Pressure_Step_Factor in the object. The Pressure_Step_Factor
   (no units) is the Pressure step change between consecutive layers when moving upwards.
   Pressure difference between the boundaries of the (n+1)-th layer (DP_n+1) will
   be DP_n**DP1. */
  Pressure getPressureStepFactor() const { return pressureStepFactor_; }

  /** Accessor to the ground altitude of site (length units) */
  Length getAltitude() const { return altitude_; }

  /** Alternative accessor to the ground altitude of site (length units) */
  Length getGroundAltitude() const { return altitude_; }

   /** setter for the ground altitude of site (length units). Careful! It will remove or add layers if necessary.
       Ground values of T/P/h would change as well */
  void setAltitude(const Length &groundaltitude);

  /** Alternative setter for the ground altitude of site (length units). Careful! It will remove or add layers if necessary.
       Ground values of T/P/h would change as well */
  // void setGroundAltitude(const Length &groundaltitude);

  /** Accessor to the altitude of the tropopause (length units) */
  Length getTropopauseAltitude() const { return tropoAltitude_; }

  /** Accessor to the Maximum allowed altitude for the Atmospheric Profile above the site (length units) */
  Length getTopAtmProfile() const { return topAtmProfile_; }

  /** Accessor to the number of layers of the atmospheric profile */
  unsigned int getNumLayer() const { return numLayer_; }

  /** Method to access the Temperature Profile  */
  vector<Temperature> getTemperatureProfile() const;

  /**
   * Method to access the average Temperature in layer i (thickness of layers in ThicknessProfile)
   * @exception AtmException if the layer is not valid.
   **/
  Temperature getLayerTemperature(unsigned int i) const;

  /**
   * Method to access the Temperature at bottom of layer i (thickness of layers in ThicknessProfile)
   * @exception AtmException if the layer is not valid.
   **/
  Temperature getLayerBottomTemperature(unsigned int i) const;

  /**
   * Method to access the Temperature at top of layer i (thickness of layers in ThicknessProfile)
   * @exception AtmException if the layer is not valid.
   **/
  Temperature getLayerTopTemperature(unsigned int i) const;

  /** Setter for the average Temperature in layer i (allows to touch one layer each time once a profile has been defined) */
  void setLayerTemperature(unsigned int i, const Temperature &layerTemperature);

  /** Method to retrieve the layer thickness from site altitude upwards.
   *  Use Altitude to + ThicknessProfile to know the vertical grid. */
  vector<Length> getThicknessProfile() const;

  /**
   * Method to access the layer thickness of layer i
   * @exception AtmException if the layer is not valid.
   */
  Length getLayerThickness(unsigned int i) const;

  /** Setter for the thickness of layer i (allows to touch one layer each time once a profile has been defined). We do
   *  not advise to use this one unless you change P and T accordingly */
  void setLayerThickness(unsigned int i, const Length &layerThickness);
  //void setLayerThickness(const Length &layerThickness, unsigned int i) { setLayerThickness(i, layerThickness); }

  /**
   * Method to access the Bottom Height of layer i above the Ground
   * @exception AtmException if the layer is not valid.
   */
  Length getLayerBottomHeightAboveGround(unsigned int i) const;

  /**
   * Method to access the Bottom Height of layer i above the Sea Level
   * @exception AtmException if the layer is not valid.
   */
  Length getLayerBottomHeightAboveSeaLevel(unsigned int i) const;


  /**
   * Method to access the Top Height of layer i above the Ground
   * @exception AtmException if the layer is not valid.
   */
  Length getLayerTopHeightAboveGround(unsigned int i) const;

  /**
   * Method to access the Top Height of layer i above the Sea Level
   * @exception AtmException if the layer is not valid.
   */
  Length getLayerTopHeightAboveSeaLevel(unsigned int i) const;


  /** Function to retrieve Average Water vapor density in a given layer in kg/m**3
   *  (thickness of layers in ThicknessProfile)
   * @exception AtmException if the layer is not valid.
   */
  MassDensity getLayerWaterVaporMassDensity(unsigned int i) const;
  MassDensity getLayerBottomWaterVaporMassDensity(unsigned int i) const;
  MassDensity getLayerTopWaterVaporMassDensity(unsigned int i) const;
  NumberDensity getLayerWaterVaporNumberDensity(unsigned int i) const;
  NumberDensity getLayerBottomWaterVaporNumberDensity(unsigned int i) const;
  NumberDensity getLayerTopWaterVaporNumberDensity(unsigned int i) const;

  /** Setter for the average Water vapor density in layer i in kg/m**3 (allows to touch one layer each
   *  time once a profile has been defined) */
  void setLayerWaterVaporMassDensity(unsigned int i, const MassDensity &layerWaterVapor);
  //void setLayerWaterVaporMassDensity(const MassDensity &layerWaterVapor, unsigned int i) { setLayerWaterVaporMassDensity(i, layerWaterVapor); }
  void setLayerWaterVaporNumberDensity(unsigned int i, const NumberDensity &layerWaterVapor);
  //void setLayerWaterVapor(const NumberDensity &layerWaterVapor, unsigned int i) { setLayerWaterVapor(i, layerWaterVapor); }

  /** Method to get the Pressure Profile */
  vector<Pressure> getPressureProfile() const;

  /** Method to access the average Pressure in layer i
   * @exception AtmException if the layer is not valid.
   */
  Pressure getLayerPressure(unsigned int i) const;

  /** Method to access the Pressure at bottom of layer i
   * @exception AtmException if the layer is not valid.
   */
  Pressure getLayerBottomPressure(unsigned int i) const;

  /** Method to access the Pressure at top of layer i
   * @exception AtmException if the layer is not valid.
   */
  Pressure getLayerTopPressure(unsigned int i) const;

  /** Setter for the average Pressure in layer i (allows to touch one layer each
   *  time once a profile has been defined) */
  void setLayerPressure(unsigned int i, const Pressure &layerPressure) { v_layerPressure_[i] = layerPressure.get(Pressure::UnitMilliBar); }
  //void setLayerPressure(const Pressure &layerPressure, unsigned int i) { setLayerPressure(i, layerPressure); }

  /** Function to retrieve CO density in a given layer (thickness of layers
   *  in ThicknessProfile)  */
  NumberDensity getLayerCO(unsigned int i) const { return NumberDensity(v_layerCO_[i], NumberDensity::UnitInverseCubicMeter); }
  /** Setter for the average number density of CO in layer i in molecules/m**3 (allows to touch one layer each
   *  time once a profile has been defined) */
  void setLayerCO(unsigned int i, const NumberDensity &layerCO) { v_layerCO_[i] = layerCO.get(NumberDensity::UnitInverseCubicMeter); }
  //void setLayerCO(const NumberDensity &layerCO, unsigned int i) { setLayerCO(i, layerCO); }

  /** Function to retrieve O3 density in a given layer (thickness of layers
   *  in ThicknessProfile) */
  NumberDensity getLayerO3(unsigned int i) const { return NumberDensity(v_layerO3_[i], NumberDensity::UnitInverseCubicMeter); }
  /** Setter for the average number density of O3 in layer i in molecules/m**3 (allows to touch one layer each
   *  time once a profile has been defined) */
  void setLayerO3(unsigned int i, const NumberDensity &layerO3) { v_layerO3_[i] = layerO3.get(NumberDensity::UnitInverseCubicMeter); }
  //void setLayerO3(const NumberDensity &layerO3, unsigned int i) { setLayerO3(i, layerO3); }

  /** Function to retrieve N2O density in a given layer (thickness of layers
   *  in ThicknessProfile)   */
  NumberDensity getLayerN2O(unsigned int i) const { return NumberDensity(v_layerN2O_[i], NumberDensity::UnitInverseCubicMeter); }
  /** Setter for the average number density of N2O in layer i in molecules/m**3 (allows to touch one layer each
   *  time once a profile has been defined) */
  void setLayerN2O(unsigned int i, const NumberDensity &layerN2O) { v_layerN2O_[i] = layerN2O.get(NumberDensity::UnitInverseCubicMeter); }
  //void setLayerN2O(const NumberDensity &layerN2O, unsigned int i) { setLayerN2O(i, layerN2O); }

  /** Function to retrieve NO2 density in a given layer (thickness of layers
   *  in ThicknessProfile)   */
  NumberDensity getLayerNO2(unsigned int i) const { return NumberDensity(v_layerNO2_[i], NumberDensity::UnitInverseCubicMeter); }
  /** Setter for the average number density of NO2 in layer i in molecules/m**3 (allows to touch one layer each
   *  time once a profile has been defined) */
  void setLayerNO2(unsigned int i, const NumberDensity &layerNO2) { v_layerNO2_[i] = layerNO2.get(NumberDensity::UnitInverseCubicMeter); }
  //void setLayerNO2(const NumberDensity &layerNO2, unsigned int i) { setLayerNO2(i, layerNO2); }

  /** Function to retrieve SO2 density in a given layer (thickness of layers
   *  in ThicknessProfile)   */
  NumberDensity getLayerSO2(unsigned int i) const { return NumberDensity(v_layerSO2_[i], NumberDensity::UnitInverseCubicMeter); }
  /** Setter for the average number density of SO2 in layer i in molecules/m**3 (allows to touch one layer each
   *  time once a profile has been defined) */
  void setLayerSO2(unsigned int i, const NumberDensity &layerSO2) { v_layerSO2_[i] = layerSO2.get(NumberDensity::UnitInverseCubicMeter); }
  //void setLayerSO2(const NumberDensity &layerSO2, unsigned int i) { setLayerSO2(i, layerSO2); }

  void setBasicAtmosphericParameterThresholds(const Length &altitudeThreshold,
					     const Pressure &groundPressureThreshold,
					     const Temperature &groundTemperatureThreshold,
					     double tropoLapseRateThreshold,
					     const Humidity &relativeHumidityThreshold,
					     const Length &wvScaleHeightThreshold);


  /** Method to get the zenith column  of water vapor. It is computed by
   *  simply integrating the H2O profile:
   *  \f$ column\;\;H_{2}O\;\; (mm) =\sum_{layers} \Delta z \cdot [ H_{2}O ] \f$
   */
  Length getGroundWH2O() const;

  // Thresholds
  Length      getAltitudeThreshold()   const {return altitudeThreshold_;};
  Pressure    getGroundPressureThreshold() const {return groundPressureThreshold_;};
  Temperature getGroundTemperatureThreshold() const {return groundTemperatureThreshold_;};
  double      getTropoLapseRateThreshold() const {return tropoLapseRateThreshold_;};
  Humidity    getRelativeHumidityThreshold() const {return relativeHumidityThreshold_;};
  Length      getWvScaleHeightThreshold() const {return wvScaleHeightThreshold_;};

  void        setAltitudeThreshold(const Length &altitudeThreshold) {altitudeThreshold_=altitudeThreshold;};
  void        setGroundPressureThreshold(const Pressure &groundPressureThreshold) {groundPressureThreshold_=groundPressureThreshold;};
  void        setGroundTemperatureThreshold(const Temperature &groundTemperatureThreshold) {groundTemperatureThreshold_=groundTemperatureThreshold;};


  //@}

protected:
  unsigned int typeAtm_; //!< 1: tropical, 2: midlatSummer, 3: midlatWinter, 4: subarcticSummer, 5: subarcticWinter
  Temperature groundTemperature_; //!< Ambient temperature at the site (K)
  double tropoLapseRate_; //!< tropospheric lapse rate in K/km
  Temperature tropoTemperature_; //!< Temperature at the tropopause
  unsigned int tropoLayer_; //!< Layer where tropopause starts
  Length tropoAltitude_; //!< Altitude where tropopause starts
  Pressure groundPressure_; //!< Ground pressure at the site
  Humidity relativeHumidity_; /** Relative humidity at the site (%)
   used only to make an estimate
   of the water vapor column, first guess) */
  Length wvScaleHeight_; //!< scale height of water vapor distribution (km)
  Pressure pressureStep_; //!< Pressure basic step (mb)
  double pressureStepFactor_; /** Multiplicative factor for presure steps.
   Example of pressure parameters:
   P_ground=550; DP: 10; DP1:
   1.2 ==> The pressure levels will
   then be 550, 560, 572, 586.4, .... */

  Length altitude_; //!< Altitude of the site (km)
  Length topAtmProfile_; //!< Top of atmospheric profile (km)
  unsigned int numLayer_; //!< Total number of layers in the output	atmospheric profiles
  double fractionLast_;  //!< Fraction of last layer needed for some calculations
  bool newBasicParam_;
  vector<double> v_layerThickness_; //!< Thickness of layer (m)
  vector<double> v_layerTemperature_; //!< Temp. of layers (K)
  vector<double> v_layerTemperature0_; //!< Temp. at bottom of layers (K)
  vector<double> v_layerTemperature1_; //!< Temp. at top of layer (K)
  vector<double> v_layerWaterVapor_; //!< Average water vapor kg/m**3 in the layer
  vector<double> v_layerWaterVapor0_; //!< Water vapor kg/m**3 at bottom of layer
  vector<double> v_layerWaterVapor1_; //!< Water vapor kg/m**3 at top of layer
  vector<double> v_layerPressure_; //!< Pressure of layers (mb)
  vector<double> v_layerPressure0_; //!< Pressure at bottom of layers (mb)
  vector<double> v_layerPressure1_; //!< Pressure at top of layers (mb)
  vector<double> v_layerCO_; //!< CO in molecules per m**3
  vector<double> v_layerO3_; //!< O3 in molecules per m**3
  vector<double> v_layerN2O_; //!< N2O in molecules per m**3
  vector<double> v_layerNO2_; //!< NO2 in molecules per m**3
  vector<double> v_layerSO2_; //!< SO2 in molecules per m**3
  Length altitudeThreshold_;
  Pressure groundPressureThreshold_;
  Temperature groundTemperatureThreshold_;
  double tropoLapseRateThreshold_;
  Humidity relativeHumidityThreshold_;
  Length wvScaleHeightThreshold_;

  unsigned int ier_;

  /** Default constructor (required if copy constructor in derived classes) */
  AtmProfile() {}

  /*
   * \fn Method to build the profile,
   */
  unsigned int mkAtmProfile(); /** returns error code: <0 unsuccessful           */

  /** Method to update an atmospheric profile based on one or more new basic parameter(s)
   * @param altitude          the new altitude, a Length
   * @param groundPressure    the Pressure at the ground level
   * @param groundTemperature the Temperature at the ground level
   * @param tropoLapseRate    the tropospheric lapse rate
   * @param relativeHumidity  the relative Humidity
   * @param wvScaleHeight     the scale height of the water vapor, a Length
   * @return true if the atmospheric profile has been update, else false because the basic parameters have not changed
   * @post   the atmospheric profile has been updated only if at least one of the basic parameters has changed
   */
  bool updateAtmProfile(const Length &altitude,
                        const Pressure &groundPressure,
                        const Temperature &groundTemperature,
                        double tropoLapseRate,
                        const Humidity &relativeHumidity,
                        const Length &wvScaleHeight);
  void initBasicAtmosphericParameterThresholds();

private:
  MassDensity rwat(const Temperature &t, const Humidity &rh, const Pressure &p) const;
  Humidity rwat_inv(const Temperature &tt, const MassDensity &dd, const Pressure &pp) const;
  vector<NumberDensity> st76(const Length &ha, unsigned int tip) const;
  double poli2(double ha,
               double x1,
               double x2,
               double x3,
               double y1,
               double y2,
               double y3) const;
  /*
   // static Pressure pressureStep_default(  10.0,Pressure::UnitMilliBar);
   //  static double pSd; // = 10.0;  // pressureStep_default in mb
   static Pressure pressureStep_default_;


   static double pressureStepFactor_default_; // =    1.2;
   //static double pSFd; // = 1.2;



   static Length  topAtmProfile_default_; //(  48.0,Length::UnitKiloMeter);
   //static tAPd; // = 48.0;
   */

  /** Default value of Atmospheric type (to reproduce behavior above the tropopause) */
  //static Atmospheretype atmType_default =  tropical;

}; // class AtmProfile

ATM_NAMESPACE_END
/*
  Pressure AtmProfile::pressureStep_default_(Pressure(1.2,Pressure::UnitMilliBar));
  double AtmProfile::pressureStepFactor_default_(1.2);
  Length AtmProfile::topAtmProfile_default_(Length(48,Length::UnitKiloMeter));
*/

/** \page AtmProfile_example Example with the AtmProfile class.
 *  Our main function always starts like this:
 *  \skip #include "AtmTypeName.h"
 *  \until {
 *
 *  First we set a constant that will be used in the calculations:
 *  \skip double
 *  \until 0.04799274551
 *  Then we set the parameters necessary to create an atmospheric profile (if you do not want to enter the whole profile yourself):
 *  \skip  Atmospheretype
 *  \until 1.2
 *  The above input information is printed on the screen
 *  \skip cout
 *  \until K/km
 *  With all this information an AtmProfile object can now be created:
 *  \skipline AtmProfile
 *  The AtmProfile object "myProflie" is accessed and partially shown on the screen
 *  \skip cout
 *  \until scale height */


#endif /*!_ATM_PROFILE_H*/
