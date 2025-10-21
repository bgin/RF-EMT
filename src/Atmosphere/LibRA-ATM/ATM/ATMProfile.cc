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
 * "@(#) $Id: ATMProfile.cpp Exp $"
 *
 * who       when      what
 * --------  --------  ----------------------------------------------
 * pardo     24/03/09  created
 */

#include "ATMProfile.h"
#include "ATMException.h"
#include <iostream>
#include <math.h>
#include <sstream>



ATM_NAMESPACE_BEGIN

AtmProfile::AtmProfile(unsigned int n)
{
  numLayer_ = n;
  initBasicAtmosphericParameterThresholds();
  for(unsigned int i = 0; i < numLayer_; ++i) {
    v_layerO3_.push_back(0.0);
    v_layerCO_.push_back(0.0);
    v_layerN2O_.push_back(0.0);
    v_layerNO2_.push_back(0.0);
    v_layerSO2_.push_back(0.0);
    v_layerThickness_.push_back(0.0);
    v_layerTemperature_.push_back(0.0);
    v_layerTemperature0_.push_back(0.0);
    v_layerTemperature1_.push_back(0.0);
    v_layerPressure_.push_back(0.0);
    v_layerPressure0_.push_back(0.0);
    v_layerPressure1_.push_back(0.0);
    v_layerWaterVapor_.push_back(0.0);
    v_layerWaterVapor0_.push_back(0.0);
    v_layerWaterVapor1_.push_back(0.0);
  }
}



AtmProfile::AtmProfile(const Length &altitude,
                       const Pressure &groundPressure,
                       const Temperature &groundTemperature,
                       double tropoLapseRate,
                       const Humidity &relativeHumidity,
                       const Length &wvScaleHeight,
                       const Pressure &pressureStep,
                       double pressureStepFactor,
                       const Length &topAtmProfile,
                       unsigned int atmType) :
      //			 Atmospheretype atmType):
      typeAtm_(atmType), groundTemperature_(groundTemperature),
      tropoLapseRate_(tropoLapseRate), groundPressure_(groundPressure),
      relativeHumidity_(relativeHumidity), wvScaleHeight_(wvScaleHeight),
      pressureStep_(pressureStep), pressureStepFactor_(pressureStepFactor),
      altitude_(altitude), topAtmProfile_(topAtmProfile)
{
  numLayer_ = 0;
  numLayer_ = mkAtmProfile();
  // std::cout << "ad numLayer_=" << numLayer_ << std::endl;
  initBasicAtmosphericParameterThresholds();
  newBasicParam_ = true;
}

AtmProfile::AtmProfile(const Length &altitude,
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
		       const vector<Temperature> &v_layerTemperature):
      typeAtm_(atmType), groundTemperature_(groundTemperature),
      tropoLapseRate_(tropoLapseRate), groundPressure_(groundPressure),
      relativeHumidity_(relativeHumidity), wvScaleHeight_(wvScaleHeight),
      pressureStep_(pressureStep), pressureStepFactor_(pressureStepFactor),
      altitude_(altitude), topAtmProfile_(topAtmProfile)
{
  numLayer_ = 0;
  numLayer_ = mkAtmProfile();
  initBasicAtmosphericParameterThresholds();
  newBasicParam_ = true;
  unsigned int nL1 = v_layerBoundaries.size();
  unsigned int nL2 = v_layerTemperature.size();
  if(nL1 == nL2 ) {
    double h=altitude_.get(Length::UnitMeter);
    double h0;
    double h1;
    double counter;
    double avT;
    for(unsigned int n = 0; n < numLayer_; n++) {
      h0=h;
      h1 = h0 + v_layerThickness_[n];
      counter = 0.0;
      avT = 0.0;
      for(unsigned int m = 0; m < nL1; m++) {
	if( h0 <= v_layerBoundaries[m].get(Length::UnitMeter) && h1 >= v_layerBoundaries[m].get(Length::UnitMeter) ){
	  //	  std::cout << "n=" << n << " h0=" << h0 << " h1=" << h1 << " v_layerBoundaries[" << m << "]=" <<  v_layerBoundaries[m].get(Length::UnitMeter) << std::endl;
	  avT = avT + v_layerTemperature[m].get(Temperature::UnitKelvin);
	  counter = counter + 1.0;
	}
      }
      if(avT > 0.0){
	// std::cout << "layer" << n << "old average temperature:" << v_layerTemperature_[n] << std::endl;
	v_layerTemperature_[n] = avT/counter;
	// std::cout << "layer" << n << "new average temperature:" << v_layerTemperature_[n] << std::endl;
      }
      h=h1;
    }
  }
}

AtmProfile::AtmProfile(const Length &altitude,
                       const Pressure &groundPressure,
                       const Temperature &groundTemperature,
                       double tropoLapseRate,
                       const Humidity &relativeHumidity,
                       const Length &wvScaleHeight,
                       unsigned int atmType) :
      //			 Atmospheretype atmType):
      typeAtm_(atmType), groundTemperature_(groundTemperature),
      tropoLapseRate_(tropoLapseRate), groundPressure_(groundPressure),
      relativeHumidity_(relativeHumidity), wvScaleHeight_(wvScaleHeight),
      pressureStep_(10.0, Pressure::UnitMilliBar), pressureStepFactor_(1.2), altitude_(altitude),
      topAtmProfile_(48.0, Length::UnitKiloMeter)
{
  numLayer_ = 0;
  numLayer_ = mkAtmProfile();
  initBasicAtmosphericParameterThresholds();
  newBasicParam_ = true;
}

AtmProfile::AtmProfile(const Length &altitude,
		       const vector<Length> &v_layerThickness,
                       const vector<Pressure> &v_layerPressure,
                       const vector<Temperature> &v_layerTemperature,
                       const vector<MassDensity> &v_layerWaterVapor)
{
  newBasicParam_ = true;
  unsigned int nL1 = v_layerThickness.size();
  unsigned int nL2 = v_layerPressure.size();
  unsigned int nL3 = v_layerTemperature.size();
  unsigned int nL4 = v_layerWaterVapor.size();

  if(nL1 + 1 == nL2 && nL2 == nL3 && nL3 == nL4) {
    numLayer_ = nL1;
    for(unsigned int n = 0; n < numLayer_; n++) {
      v_layerO3_.push_back(0.0);
      v_layerCO_.push_back(0.0);
      v_layerN2O_.push_back(0.0);
      v_layerNO2_.push_back(0.0);
      v_layerSO2_.push_back(0.0);
      v_layerThickness_.push_back(v_layerThickness[n].get(Length::UnitMeter));
      v_layerTemperature_.push_back((v_layerTemperature[n].get(Temperature::UnitKelvin) + v_layerTemperature[n+1].get(Temperature::UnitKelvin))/2.0);
      v_layerTemperature0_.push_back(v_layerTemperature[n].get(Temperature::UnitKelvin));
      v_layerTemperature1_.push_back(v_layerTemperature[n+1].get(Temperature::UnitKelvin));
      v_layerPressure_.push_back(exp((log(v_layerPressure[n].get(Pressure::UnitMilliBar))+ log(v_layerPressure[n + 1].get(Pressure::UnitMilliBar)))/2.0));
      v_layerPressure0_.push_back(v_layerPressure[n].get(Pressure::UnitMilliBar));
      v_layerPressure1_.push_back(v_layerPressure[n+1].get(Pressure::UnitMilliBar));
      v_layerWaterVapor_.push_back(exp((log(v_layerWaterVapor[n].get(MassDensity::UnitKiloGramPerCubicMeter))+ log(v_layerWaterVapor[n + 1].get(MassDensity::UnitKiloGramPerCubicMeter)))/2.0));
      v_layerWaterVapor0_.push_back(v_layerWaterVapor[n].get(MassDensity::UnitKiloGramPerCubicMeter));
      v_layerWaterVapor1_.push_back(v_layerWaterVapor[n+1].get(MassDensity::UnitKiloGramPerCubicMeter));
    }
  } else {
    numLayer_ = 0;
  }
  initBasicAtmosphericParameterThresholds();
}

AtmProfile::AtmProfile(const vector<Length> &v_layerBoundaries,
                       const vector<Pressure> &v_layerPressure,
                       const vector<Temperature> &v_layerTemperature,
                       const vector<MassDensity> &v_layerWaterVapor)
{
  newBasicParam_ = true;
  unsigned int nL1 = v_layerBoundaries.size();
  unsigned int nL2 = v_layerPressure.size();
  unsigned int nL3 = v_layerTemperature.size();
  unsigned int nL4 = v_layerWaterVapor.size();

  if(nL1 == nL2 && nL2 == nL3 && nL3 == nL4) {
    numLayer_ = nL1 - 1;
    altitude_ = v_layerBoundaries[0];
    for(unsigned int n = 0; n < numLayer_; n++) {
      v_layerO3_.push_back(0.0);
      v_layerCO_.push_back(0.0);
      v_layerN2O_.push_back(0.0);
      v_layerNO2_.push_back(0.0);
      v_layerSO2_.push_back(0.0);
      v_layerThickness_.push_back(v_layerBoundaries[n+1].get(Length::UnitMeter) - v_layerBoundaries[n].get(Length::UnitMeter));
      v_layerTemperature_.push_back((v_layerTemperature[n].get(Temperature::UnitKelvin) + v_layerTemperature[n+1].get(Temperature::UnitKelvin))/2.0);
      v_layerTemperature0_.push_back(v_layerTemperature[n].get(Temperature::UnitKelvin));
      v_layerTemperature1_.push_back(v_layerTemperature[n+1].get(Temperature::UnitKelvin));
      v_layerPressure_.push_back(exp((log(v_layerPressure[n].get(Pressure::UnitMilliBar))+ log(v_layerPressure[n + 1].get(Pressure::UnitMilliBar)))/2.0));
      v_layerPressure0_.push_back(v_layerPressure[n].get(Pressure::UnitMilliBar));
      v_layerPressure1_.push_back(v_layerPressure[n+1].get(Pressure::UnitMilliBar));
      v_layerWaterVapor_.push_back(exp((log(v_layerWaterVapor[n].get(MassDensity::UnitKiloGramPerCubicMeter))+ log(v_layerWaterVapor[n + 1].get(MassDensity::UnitKiloGramPerCubicMeter)))/2.0));
      v_layerWaterVapor0_.push_back(v_layerWaterVapor[n].get(MassDensity::UnitKiloGramPerCubicMeter));
      v_layerWaterVapor1_.push_back(v_layerWaterVapor[n+1].get(MassDensity::UnitKiloGramPerCubicMeter));
    }
  } else {
    numLayer_ = 0;
  }
  initBasicAtmosphericParameterThresholds();
}

AtmProfile::AtmProfile(const Length &altitude,
		       const vector<Length> &v_layerThickness,
                       const vector<Pressure> &v_layerPressure,
                       const vector<Temperature> &v_layerTemperature,
                       const vector<NumberDensity> &v_layerWaterVapor)
{
  newBasicParam_ = true;
  unsigned int nL1 = v_layerThickness.size();
  unsigned int nL2 = v_layerPressure.size();
  unsigned int nL3 = v_layerTemperature.size();
  unsigned int nL4 = v_layerWaterVapor.size();

  if(nL1 + 1 == nL2 && nL2 == nL3 && nL3 == nL4) {
    numLayer_ = nL1;
    for(unsigned int n = 0; n < numLayer_; n++) {
      v_layerO3_.push_back(0.0);
      v_layerCO_.push_back(0.0);
      v_layerN2O_.push_back(0.0);
      v_layerNO2_.push_back(0.0);
      v_layerSO2_.push_back(0.0);
      v_layerThickness_.push_back(v_layerThickness[n].get(Length::UnitMeter));
      v_layerTemperature_.push_back((v_layerTemperature[n].get(Temperature::UnitKelvin) + v_layerTemperature[n+1].get(Temperature::UnitKelvin))/2.0);
      v_layerTemperature0_.push_back(v_layerTemperature[n].get(Temperature::UnitKelvin));
      v_layerTemperature1_.push_back(v_layerTemperature[n+1].get(Temperature::UnitKelvin));
      v_layerPressure_.push_back(exp((log(v_layerPressure[n].get(Pressure::UnitMilliBar))+ log(v_layerPressure[n + 1].get(Pressure::UnitMilliBar)))/2.0));
      v_layerPressure0_.push_back(v_layerPressure[n].get(Pressure::UnitMilliBar));
      v_layerPressure1_.push_back(v_layerPressure[n+1].get(Pressure::UnitMilliBar));
      v_layerWaterVapor_.push_back( (exp((log(v_layerWaterVapor[n].get(NumberDensity::UnitInverseCubicMeter))+ log(v_layerWaterVapor[n + 1].get(NumberDensity::UnitInverseCubicMeter)))/2.0))* 18.0
				    / (1000.0 * 6.023e23));
      v_layerWaterVapor0_.push_back(v_layerWaterVapor[n].get(NumberDensity::UnitInverseCubicMeter) * 18.0
          / (1000.0 * 6.023e23));
      v_layerWaterVapor1_.push_back(v_layerWaterVapor[n+1].get(NumberDensity::UnitInverseCubicMeter) * 18.0
          / (1000.0 * 6.023e23));
    }
  } else {
    numLayer_ = 0;
  }
  initBasicAtmosphericParameterThresholds();
}

AtmProfile::AtmProfile(const vector<Length> &v_layerBoundaries,
                       const vector<Pressure> &v_layerPressure,
                       const vector<Temperature> &v_layerTemperature,
                       const vector<NumberDensity> &v_layerWaterVapor)
{
  newBasicParam_ = true;
  unsigned int nL1 = v_layerBoundaries.size();
  unsigned int nL2 = v_layerPressure.size();
  unsigned int nL3 = v_layerTemperature.size();
  unsigned int nL4 = v_layerWaterVapor.size();

  if(nL1 == nL2 && nL2 == nL3 && nL3 == nL4) {
    numLayer_ = nL1 - 1;
    altitude_ = v_layerBoundaries[0];
    for(unsigned int n = 0; n < numLayer_; n++) {
      v_layerO3_.push_back(0.0);
      v_layerCO_.push_back(0.0);
      v_layerN2O_.push_back(0.0);
      v_layerNO2_.push_back(0.0);
      v_layerSO2_.push_back(0.0);
      v_layerThickness_.push_back(v_layerBoundaries[n+1].get(Length::UnitMeter) - v_layerBoundaries[n].get(Length::UnitMeter));
      v_layerTemperature_.push_back((v_layerTemperature[n].get(Temperature::UnitKelvin) + v_layerTemperature[n+1].get(Temperature::UnitKelvin))/2.0);
      v_layerTemperature0_.push_back(v_layerTemperature[n].get(Temperature::UnitKelvin));
      v_layerTemperature1_.push_back(v_layerTemperature[n+1].get(Temperature::UnitKelvin));
      v_layerPressure_.push_back(exp((log(v_layerPressure[n].get(Pressure::UnitMilliBar))+ log(v_layerPressure[n + 1].get(Pressure::UnitMilliBar)))/2.0));
      v_layerPressure0_.push_back(v_layerPressure[n].get(Pressure::UnitMilliBar));
      v_layerPressure1_.push_back(v_layerPressure[n+1].get(Pressure::UnitMilliBar));
      v_layerWaterVapor_.push_back( (exp((log(v_layerWaterVapor[n].get(NumberDensity::UnitInverseCubicMeter))+ log(v_layerWaterVapor[n + 1].get(NumberDensity::UnitInverseCubicMeter)))/2.0))* 18.0
				     / (1000.0 * 6.023e23));
      v_layerWaterVapor0_.push_back(v_layerWaterVapor[n].get(NumberDensity::UnitInverseCubicMeter) * 18.0
          / (1000.0 * 6.023e23));
      v_layerWaterVapor1_.push_back(v_layerWaterVapor[n+1].get(NumberDensity::UnitInverseCubicMeter) * 18.0
          / (1000.0 * 6.023e23));
    }
  } else {
    numLayer_ = 0;
  }
  initBasicAtmosphericParameterThresholds();
}

AtmProfile::AtmProfile(const Length &altitude,
		       const vector<Length> &v_layerThickness,
                       const vector<Pressure> &v_layerPressure,
                       const vector<Temperature> &v_layerTemperature,
                       const vector<MassDensity> &v_layerWaterVapor,
                       const vector<NumberDensity> &v_layerO3)
{
  newBasicParam_ = true;
  unsigned int nL1 = v_layerThickness.size();
  unsigned int nL2 = v_layerPressure.size();
  unsigned int nL3 = v_layerTemperature.size();
  unsigned int nL4 = v_layerWaterVapor.size();
  unsigned int nL5 = v_layerO3.size();

  if(nL1 * 1 == nL2 && nL2 == nL3 && nL3 == nL4 && nL4 == nL5) {
    numLayer_ = nL1;
    for(unsigned int n = 0; n < numLayer_; n++) {
      v_layerO3_.push_back(v_layerO3[n].get(NumberDensity::UnitInverseCubicMeter));
      v_layerCO_.push_back(0.0);
      v_layerN2O_.push_back(0.0);
      v_layerNO2_.push_back(0.0);
      v_layerSO2_.push_back(0.0);
      v_layerThickness_.push_back(v_layerThickness[n].get(Length::UnitMeter));
      v_layerTemperature_.push_back((v_layerTemperature[n].get(Temperature::UnitKelvin) + v_layerTemperature[n+1].get(Temperature::UnitKelvin))/2.0);
      v_layerTemperature0_.push_back(v_layerTemperature[n].get(Temperature::UnitKelvin));
      v_layerTemperature1_.push_back(v_layerTemperature[n+1].get(Temperature::UnitKelvin));
      v_layerPressure_.push_back(exp((log(v_layerPressure[n].get(Pressure::UnitMilliBar))+ log(v_layerPressure[n + 1].get(Pressure::UnitMilliBar)))/2.0));
      v_layerPressure0_.push_back(v_layerPressure[n].get(Pressure::UnitMilliBar));
      v_layerPressure1_.push_back(v_layerPressure[n+1].get(Pressure::UnitMilliBar));
      v_layerWaterVapor_.push_back(exp((log(v_layerWaterVapor[n].get(MassDensity::UnitKiloGramPerCubicMeter))+ log(v_layerWaterVapor[n + 1].get(MassDensity::UnitKiloGramPerCubicMeter)))/2.0));
      v_layerWaterVapor0_.push_back(v_layerWaterVapor[n].get(MassDensity::UnitKiloGramPerCubicMeter));
      v_layerWaterVapor1_.push_back(v_layerWaterVapor[n+1].get(MassDensity::UnitKiloGramPerCubicMeter));
    }
  } else {
    numLayer_ = 0;
  }
  initBasicAtmosphericParameterThresholds();
}

AtmProfile::AtmProfile(const Length &altitude,
		       const vector<Length> &v_layerThickness,
                       const vector<Pressure> &v_layerPressure,
                       const vector<Temperature> &v_layerTemperature,
                       const vector<NumberDensity> &v_layerWaterVapor,
                       const vector<NumberDensity> &v_layerO3)
{
  newBasicParam_ = true;
  unsigned int nL1 = v_layerThickness.size();
  unsigned int nL2 = v_layerPressure.size();
  unsigned int nL3 = v_layerTemperature.size();
  unsigned int nL4 = v_layerWaterVapor.size();
  unsigned int nL5 = v_layerO3.size();

  if(nL1 + 1 == nL2 && nL2 == nL3 && nL3 == nL4 && nL4 == nL5) {
    numLayer_ = nL1;
    for(unsigned int n = 0; n < numLayer_; n++) {
      v_layerO3_.push_back(v_layerO3[n].get(NumberDensity::UnitInverseCubicMeter));
      v_layerCO_.push_back(0.0);
      v_layerN2O_.push_back(0.0);
      v_layerNO2_.push_back(0.0);
      v_layerSO2_.push_back(0.0);
      v_layerThickness_.push_back(v_layerThickness[n].get(Length::UnitMeter));
      v_layerTemperature_.push_back((v_layerTemperature[n].get(Temperature::UnitKelvin) + v_layerTemperature[n+1].get(Temperature::UnitKelvin))/2.0);
      v_layerTemperature0_.push_back(v_layerTemperature[n].get(Temperature::UnitKelvin));
      v_layerTemperature1_.push_back(v_layerTemperature[n+1].get(Temperature::UnitKelvin));
      v_layerPressure_.push_back(exp((log(v_layerPressure[n].get(Pressure::UnitMilliBar))+ log(v_layerPressure[n + 1].get(Pressure::UnitMilliBar)))/2.0));
      v_layerPressure0_.push_back(v_layerPressure[n].get(Pressure::UnitMilliBar));
      v_layerPressure1_.push_back(v_layerPressure[n+1].get(Pressure::UnitMilliBar));
      v_layerWaterVapor_.push_back( (exp((log(v_layerWaterVapor[n].get(NumberDensity::UnitInverseCubicMeter))+ log(v_layerWaterVapor[n + 1].get(NumberDensity::UnitInverseCubicMeter)))/2.0))* 18.0
				    / (1000.0 * 6.023e23));
      v_layerWaterVapor0_.push_back(v_layerWaterVapor[n].get(NumberDensity::UnitInverseCubicMeter) * 18.0
          / (1000.0 * 6.023e23));
      v_layerWaterVapor1_.push_back(v_layerWaterVapor[n+1].get(NumberDensity::UnitInverseCubicMeter) * 18.0
          / (1000.0 * 6.023e23));
    }
  } else {
    numLayer_ = 0;
  }
  initBasicAtmosphericParameterThresholds();
}

AtmProfile::AtmProfile(const Length &altitude,
		       const vector<Length> &v_layerThickness,
                       const vector<Pressure> &v_layerPressure,
                       const vector<Temperature> &v_layerTemperature,
                       const vector<MassDensity> &v_layerWaterVapor,
                       const vector<NumberDensity> &v_layerO3,
                       const vector<NumberDensity> &v_layerCO,
                       const vector<NumberDensity> &v_layerN2O,
                       const vector<NumberDensity> &v_layerNO2,
                       const vector<NumberDensity> &v_layerSO2)
{
  newBasicParam_ = true;
  unsigned int nL1 = v_layerThickness.size();
  unsigned int nL2 = v_layerPressure.size();
  unsigned int nL3 = v_layerTemperature.size();
  unsigned int nL4 = v_layerWaterVapor.size();
  unsigned int nL5 = v_layerO3.size();
  unsigned int nL6 = v_layerCO.size();
  unsigned int nL7 = v_layerN2O.size();
  unsigned int nL8 = v_layerNO2.size();
  unsigned int nL9 = v_layerSO2.size();

  if(nL1 + 1 == nL2 && nL2 == nL3 && nL3 == nL4 && nL4 == nL5 && nL5 == nL6 && nL6
      == nL7 && nL7 == nL8 && nL8 == nL9) {
    numLayer_ = nL1;
    for(unsigned int n = 0; n < numLayer_; n++) {
      v_layerO3_.push_back(v_layerO3[n].get(NumberDensity::UnitInverseCubicMeter));
      v_layerCO_.push_back(v_layerCO[n].get(NumberDensity::UnitInverseCubicMeter));
      v_layerN2O_.push_back(v_layerN2O[n].get(NumberDensity::UnitInverseCubicMeter));
      v_layerNO2_.push_back(v_layerNO2[n].get(NumberDensity::UnitInverseCubicMeter));
      v_layerSO2_.push_back(v_layerSO2[n].get(NumberDensity::UnitInverseCubicMeter));
      v_layerThickness_.push_back(v_layerThickness[n].get(Length::UnitMeter));
      v_layerTemperature_.push_back((v_layerTemperature[n].get(Temperature::UnitKelvin) + v_layerTemperature[n+1].get(Temperature::UnitKelvin))/2.0);
      v_layerTemperature0_.push_back(v_layerTemperature[n].get(Temperature::UnitKelvin));
      v_layerTemperature1_.push_back(v_layerTemperature[n+1].get(Temperature::UnitKelvin));
      v_layerPressure_.push_back(exp((log(v_layerPressure[n].get(Pressure::UnitMilliBar))+ log(v_layerPressure[n + 1].get(Pressure::UnitMilliBar)))/2.0));
      v_layerPressure0_.push_back(v_layerPressure[n].get(Pressure::UnitMilliBar));
      v_layerPressure1_.push_back(v_layerPressure[n+1].get(Pressure::UnitMilliBar));
      v_layerWaterVapor_.push_back(exp((log(v_layerWaterVapor[n].get(MassDensity::UnitKiloGramPerCubicMeter))+ log(v_layerWaterVapor[n + 1].get(MassDensity::UnitKiloGramPerCubicMeter)))/2.0));
      v_layerWaterVapor0_.push_back(v_layerWaterVapor[n].get(MassDensity::UnitKiloGramPerCubicMeter));
      v_layerWaterVapor1_.push_back(v_layerWaterVapor[n+1].get(MassDensity::UnitKiloGramPerCubicMeter));
    }
  } else {
    numLayer_ = 0;
  }
  initBasicAtmosphericParameterThresholds();
}

AtmProfile::AtmProfile(const Length &altitude,
		       const vector<Length> &v_layerThickness,
                       const vector<Pressure> &v_layerPressure,
                       const vector<Temperature> &v_layerTemperature,
                       const vector<NumberDensity> &v_layerWaterVapor,
                       const vector<NumberDensity> &v_layerO3,
                       const vector<NumberDensity> &v_layerCO,
                       const vector<NumberDensity> &v_layerN2O,
                       const vector<NumberDensity> &v_layerNO2,
                       const vector<NumberDensity> &v_layerSO2)
{
  newBasicParam_ = true;
  unsigned int nL1 = v_layerThickness.size();
  unsigned int nL2 = v_layerPressure.size();
  unsigned int nL3 = v_layerTemperature.size();
  unsigned int nL4 = v_layerWaterVapor.size();
  unsigned int nL5 = v_layerO3.size();
  unsigned int nL6 = v_layerCO.size();
  unsigned int nL7 = v_layerN2O.size();
  unsigned int nL8 = v_layerNO2.size();
  unsigned int nL9 = v_layerSO2.size();

  if(nL1 + 1 == nL2 && nL2 == nL3 && nL3 == nL4 && nL4 == nL5 && nL5 == nL6 && nL6
      == nL7 && nL7 == nL8 && nL8 == nL9) {
    numLayer_ = nL1;
    for(unsigned int n = 0; n < numLayer_; n++) {
      v_layerO3_.push_back(v_layerO3[n].get(NumberDensity::UnitInverseCubicMeter));
      v_layerCO_.push_back(v_layerCO[n].get(NumberDensity::UnitInverseCubicMeter));
      v_layerN2O_.push_back(v_layerN2O[n].get(NumberDensity::UnitInverseCubicMeter));
      v_layerNO2_.push_back(v_layerNO2[n].get(NumberDensity::UnitInverseCubicMeter));
      v_layerSO2_.push_back(v_layerSO2[n].get(NumberDensity::UnitInverseCubicMeter));
      v_layerThickness_.push_back(v_layerThickness[n].get(Length::UnitMeter));
      v_layerTemperature_.push_back((v_layerTemperature[n].get(Temperature::UnitKelvin) + v_layerTemperature[n+1].get(Temperature::UnitKelvin))/2.0);
      v_layerTemperature0_.push_back(v_layerTemperature[n].get(Temperature::UnitKelvin));
      v_layerTemperature1_.push_back(v_layerTemperature[n+1].get(Temperature::UnitKelvin));
      v_layerPressure_.push_back(exp((log(v_layerPressure[n].get(Pressure::UnitMilliBar))+ log(v_layerPressure[n + 1].get(Pressure::UnitMilliBar)))/2.0));
      v_layerPressure0_.push_back(v_layerPressure[n].get(Pressure::UnitMilliBar));
      v_layerPressure1_.push_back(v_layerPressure[n+1].get(Pressure::UnitMilliBar));
      v_layerWaterVapor_.push_back( (exp((log(v_layerWaterVapor[n].get(NumberDensity::UnitInverseCubicMeter))+ log(v_layerWaterVapor[n + 1].get(NumberDensity::UnitInverseCubicMeter)))/2.0))* 18.0
				    / (1000.0 * 6.023e23));
      v_layerWaterVapor0_.push_back(v_layerWaterVapor[n].get(NumberDensity::UnitInverseCubicMeter) * 18.0
          / (1000.0 * 6.023e23));
      v_layerWaterVapor1_.push_back(v_layerWaterVapor[n+1].get(NumberDensity::UnitInverseCubicMeter) * 18.0
          / (1000.0 * 6.023e23));
    }
  } else {
    numLayer_ = 0;
  }
  initBasicAtmosphericParameterThresholds();
}

AtmProfile::AtmProfile(const AtmProfile &a)
{ //:AtmType(a.type_){
  // std::cout<<"AtmProfile copy constructor"<<endl;  COMMENTED OUT BY JUAN MAY/16/2005
  typeAtm_ = a.typeAtm_;
  groundTemperature_ = a.groundTemperature_;
  tropoLapseRate_ = a.tropoLapseRate_;
  groundPressure_ = a.groundPressure_;
  relativeHumidity_ = a.relativeHumidity_;
  wvScaleHeight_ = a.wvScaleHeight_;
  pressureStep_ = a.pressureStep_;
  pressureStepFactor_ = a.pressureStepFactor_;
  altitude_ = a.altitude_;
  topAtmProfile_ = a.topAtmProfile_;
  numLayer_ = a.numLayer_;
  newBasicParam_ = a.newBasicParam_;
  v_layerThickness_.reserve(numLayer_);
  v_layerPressure_.reserve(numLayer_);
  v_layerPressure0_.reserve(numLayer_);
  v_layerPressure1_.reserve(numLayer_);
  v_layerTemperature_.reserve(numLayer_);
  v_layerTemperature0_.reserve(numLayer_);
  v_layerTemperature1_.reserve(numLayer_);
  v_layerWaterVapor_.reserve(numLayer_);
  v_layerWaterVapor0_.reserve(numLayer_);
  v_layerWaterVapor1_.reserve(numLayer_);
  v_layerCO_.reserve(numLayer_);
  v_layerO3_.reserve(numLayer_);
  v_layerN2O_.reserve(numLayer_);
  v_layerNO2_.reserve(numLayer_);
  v_layerSO2_.reserve(numLayer_);
  // std::cout << "numLayer_=" << numLayer_ << std::endl;  COMMENTED OUT BY JUAN MAY/16/2005
  for(unsigned int n = 0; n < numLayer_; n++) {
    v_layerThickness_.push_back(a.v_layerThickness_[n]);
    v_layerTemperature_.push_back(a.v_layerTemperature_[n]);
    v_layerTemperature0_.push_back(a.v_layerTemperature0_[n]);
    v_layerTemperature1_.push_back(a.v_layerTemperature1_[n]);
    //cout << "n=" << n << std::endl;
    v_layerWaterVapor_.push_back(a.v_layerWaterVapor_[n]);
    v_layerWaterVapor0_.push_back(a.v_layerWaterVapor0_[n]);
    v_layerWaterVapor1_.push_back(a.v_layerWaterVapor1_[n]);
    v_layerPressure_.push_back(a.v_layerPressure_[n]);
    v_layerPressure0_.push_back(a.v_layerPressure1_[n]);
    v_layerPressure1_.push_back(a.v_layerPressure1_[n]);
    v_layerCO_.push_back(a.v_layerCO_[n]);
    v_layerO3_.push_back(a.v_layerO3_[n]);
    v_layerN2O_.push_back(a.v_layerN2O_[n]);
    v_layerNO2_.push_back(a.v_layerNO2_[n]);
    v_layerSO2_.push_back(a.v_layerSO2_[n]);
  }
  altitudeThreshold_ = a.altitudeThreshold_;
  groundPressureThreshold_ = a.groundPressureThreshold_;
  groundTemperatureThreshold_ = a.groundTemperatureThreshold_;
  tropoLapseRateThreshold_ = a.tropoLapseRateThreshold_;
  relativeHumidityThreshold_ = a.relativeHumidityThreshold_;
  wvScaleHeightThreshold_ = a.wvScaleHeightThreshold_;
}
void AtmProfile::setBasicAtmosphericParameterThresholds(const Length &altitudeThreshold,
							const Pressure &groundPressureThreshold,
							const Temperature &groundTemperatureThreshold,
							double tropoLapseRateThreshold,
							const Humidity &relativeHumidityThreshold,
							const Length &wvScaleHeightThreshold)
{
  altitudeThreshold_ = altitudeThreshold;
  groundPressureThreshold_ = groundPressureThreshold;
  groundTemperatureThreshold_ = groundTemperatureThreshold;
  tropoLapseRateThreshold_ = tropoLapseRateThreshold;
  relativeHumidityThreshold_ = relativeHumidityThreshold;
  wvScaleHeightThreshold_ =wvScaleHeightThreshold;
}

void AtmProfile::initBasicAtmosphericParameterThresholds()
{
  altitudeThreshold_ = Length(1.0,Length::UnitMeter);
  groundPressureThreshold_ =  Pressure(99.,Pressure::UnitPascal);  // DB 2014-03-06 : Choose 99 Pascal instead of 100 Pascal because the threshold must be lower than the value of delta_pressure used in getAverageDispersiveDryPathLength_GroundPressureDerivative()
  groundTemperatureThreshold_ =  Temperature(0.3,Temperature::UnitKelvin);
  tropoLapseRateThreshold_ = 0.01;
  relativeHumidityThreshold_ = Humidity(100.,Percent::UnitPercent);
  wvScaleHeightThreshold_ = Length(20.,Length::UnitMeter);
}

bool AtmProfile::updateAtmProfile(const Length &altitude,
                                  const Pressure &groundPressure,
                                  const Temperature &groundTemperature,
                                  double tropoLapseRate,
                                  const Humidity &relativeHumidity,
                                  const Length &wvScaleHeight)
{

  /* TODO A faire: pour decider s'il faut recalculer le profile on devrait plutot donner des seuils, eg
   if(fabs(altitude_.get()-altitude.get())>0.1)mkNewProfile=true;
   */

  unsigned int numLayer;
  bool mkNewProfile = false;

//   if(altitude_.get() != altitude.get()) mkNewProfile = true; //if(mkNewProfile)cout<<"altitude has changed"          <<endl;
//   if(groundPressure_.get() != groundPressure.get()) mkNewProfile = true; //if(mkNewProfile)cout<<"ground pressure has changed"   <<endl;
//   if(groundTemperature_.get() != groundTemperature.get()) mkNewProfile = true; //if(mkNewProfile)cout<<"ground temperature has changed"<<endl;
//   if(wvScaleHeight_.get() != wvScaleHeight.get()) mkNewProfile = true; //if(mkNewProfile)cout<<"wv scale height has changed"   <<endl;
//   if(tropoLapseRate_ != tropoLapseRate) mkNewProfile = true; //if(mkNewProfile)cout<<"tropo lapse rate has changed"  <<endl;
//   if(relativeHumidity_.get() != relativeHumidity.get()) mkNewProfile = true; //if(mkNewProfile)cout<<"relative humidity has changed" <<endl;
  //! all these thresholds shoudl be parametrized.
  if(fabs(altitude_.get()-altitude.get())> altitudeThreshold_.get()) { // Length(0.1,Length::UnitMeter).get())  {
    mkNewProfile = true;
    //cout<<"altitude has changed"          <<endl;
  }
  if(fabs(groundPressure_.get()-groundPressure.get())> groundPressureThreshold_.get()) { // Pressure(100.,Pressure::UnitPascal).get()) {
    mkNewProfile = true;
    //cout<<"ground pressure has changed"   <<endl;
  }
  if(fabs(groundTemperature_.get()-groundTemperature.get()) > groundTemperatureThreshold_.get()) { //  Temperature(0.3,Temperature::UnitKelvin).get()) {
    mkNewProfile = true;
    //cout<<"ground temperature has changed"<<endl;
  }
  if(fabs(wvScaleHeight_.get()-wvScaleHeight.get()) > wvScaleHeightThreshold_.get()) { // Length(20.,Length::UnitMeter).get() ) {
    mkNewProfile = true;
    //cout<<"wv scale height has changed"   <<endl;
  }
  if(fabs(tropoLapseRate_-tropoLapseRate) > tropoLapseRateThreshold_) { // 0.01) {
    mkNewProfile = true;
    //cout<<"tropo lapse rate has changed"  <<endl;
  }
  // we use WVR...
  if(fabs(relativeHumidity_.get()-relativeHumidity.get())> relativeHumidityThreshold_.get()) { // Humidity(100.,Percent::UnitPercent).get()) {
    mkNewProfile = true;
    //cout<<"relative humidity has changed" <<endl;
  }
  if(mkNewProfile) {
    newBasicParam_ = true;
    altitude_ = altitude;
    groundPressure_ = groundPressure;
    groundTemperature_ = groundTemperature;
    tropoLapseRate_ = tropoLapseRate;
    relativeHumidity_ = relativeHumidity;
    wvScaleHeight_ = wvScaleHeight;
    numLayer = mkAtmProfile();
    numLayer_ = numLayer;
    //      std::cout << "There are new basic parameters, with " << numLayer_ << " layers " << std::endl;
  } else {
    numLayer = getNumLayer();
    numLayer_ = numLayer;
  }

  return mkNewProfile;
}

// Note that this setBasicAtmosphericParameters will be overrided by the subclasses.
bool AtmProfile::setBasicAtmosphericParameters(const Length &altitude,
                                               const Pressure &groundPressure,
                                               const Temperature &groundTemperature,
                                               double tropoLapseRate,
                                               const Humidity &relativeHumidity,
                                               const Length &wvScaleHeight)
{
  return updateAtmProfile(altitude,
                          groundPressure,
                          groundTemperature,
                          tropoLapseRate,
                          relativeHumidity,
                          wvScaleHeight);
}

string AtmProfile::getAtmosphereType() const
{
  return getAtmosphereType(typeAtm_);
}

string AtmProfile::getAtmosphereType(unsigned int typeAtm)
{
    string type;
    string typeNames[] = { "TROPICAL", "MIDLATSUMMER", "MIDLATWINTER",
                           "SUBARTSUMMER", "SUBARTWINTER","US_ST76"};

    if(typeAtm < typeATM_end) {
      type = typeNames[typeAtm-1];
    } else {
      type = "DEFAULT";
    }
    return type;
}

vector<Temperature> AtmProfile::getTemperatureProfile() const
{
  vector<Temperature> t;
  t.reserve(v_layerTemperature_.size());
  for(unsigned int i = 0; i < v_layerTemperature_.size(); i++) {
    Temperature tt(v_layerTemperature_[i], Temperature::UnitKelvin);
    t.push_back(tt);
  }
  return t;
}

Temperature AtmProfile::getLayerTemperature(unsigned int i) const
{
  /*if(i > v_layerTemperature_.size() - 1) {
    Temperature t(-999.0, Temperature::UnitKelvin);
    return t;
  } else {
    Temperature t(v_layerTemperature_[i], Temperature::UnitKelvin);
    return t;
  }*/
  if(i > v_layerTemperature_.size() - 1) {
    std::ostringstream oss;
    oss << "Not a valid layer: " << i;
    throw AtmException(ATM_EXCEPTION_ARGS(oss.str().c_str()));
  }
  return Temperature(v_layerTemperature_[i], Temperature::UnitKelvin);
}

Temperature AtmProfile::getLayerBottomTemperature(unsigned int i) const
{
  if(i > v_layerTemperature0_.size() - 1) {
    std::ostringstream oss;
    oss << "Not a valid layer: " << i;
    throw AtmException(ATM_EXCEPTION_ARGS(oss.str().c_str()));
  }
  return Temperature(v_layerTemperature0_[i], Temperature::UnitKelvin);
}

Temperature AtmProfile::getLayerTopTemperature(unsigned int i) const
{
  if(i > v_layerTemperature1_.size() - 1) {
    std::ostringstream oss;
    oss << "Not a valid layer: " << i;
    throw AtmException(ATM_EXCEPTION_ARGS(oss.str().c_str()));
  }
  return Temperature(v_layerTemperature1_[i], Temperature::UnitKelvin);
}


void AtmProfile::setAltitude(const Length &groundaltitude)
{

  if (groundaltitude <= altitude_){

    // std::cout << "extrapolar para abajo" << std::endl;

    int nextralayers = int( 0.50001+(altitude_-groundaltitude).get(Length::UnitMeter)/v_layerThickness_[0]);
    if (nextralayers == 0) {nextralayers=1;}
    //    std::cout << "aaa=" << ( 0.50001+(altitude_-groundaltitude).get(Length::UnitMeter)/v_layerThickness_[0]) << std::endl;
    double newThickness = (altitude_-groundaltitude).get(Length::UnitMeter) / nextralayers;

    //  std::cout << "Number of extra layers: " << nextralayers << " of thickness = " <<  newThickness << " m"  << std::endl;
    //  std::cout << "nextralayers-1 " << nextralayers-1 << std::endl;

    for(int i=nextralayers; i>0; i--){
      v_layerThickness_.push_back(newThickness);
      v_layerTemperature_.push_back(v_layerTemperature_[0]);
      v_layerTemperature0_.push_back(v_layerTemperature0_[0]);
      v_layerTemperature1_.push_back(v_layerTemperature1_[0]);
      v_layerWaterVapor_.push_back(v_layerWaterVapor_[0]);
      v_layerWaterVapor0_.push_back(v_layerWaterVapor0_[0]);
      v_layerWaterVapor1_.push_back(v_layerWaterVapor1_[0]);
      v_layerPressure_.push_back(v_layerPressure_[0]);
      v_layerPressure0_.push_back(v_layerPressure0_[0]);
      v_layerPressure1_.push_back(v_layerPressure1_[0]);
      v_layerCO_.push_back(v_layerCO_[0]);
      v_layerO3_.push_back(v_layerO3_[0]);
      v_layerN2O_.push_back(v_layerN2O_[0]);
      v_layerNO2_.push_back(v_layerNO2_[0]);
      v_layerSO2_.push_back(v_layerSO2_[0]);
    }

    for(int i=v_layerThickness_.size()-1; i>nextralayers-1; i--){
      v_layerThickness_[i] = v_layerThickness_[i-nextralayers];
      v_layerTemperature_[i] = v_layerTemperature_[i-nextralayers];
      v_layerTemperature0_[i] = v_layerTemperature0_[i-nextralayers];
      v_layerTemperature1_[i] = v_layerTemperature1_[i-nextralayers];
      v_layerWaterVapor_[i] = v_layerWaterVapor_[i-nextralayers];
      v_layerWaterVapor0_[i] = v_layerWaterVapor0_[i-nextralayers];
      v_layerWaterVapor1_[i] = v_layerWaterVapor1_[i-nextralayers];
      v_layerPressure_[i] = v_layerPressure_[i-nextralayers];
      v_layerPressure0_[i] = v_layerPressure0_[i-nextralayers];
      v_layerPressure1_[i] = v_layerPressure1_[i-nextralayers];
      v_layerCO_[i] = v_layerCO_[i-nextralayers];
      v_layerO3_[i] = v_layerO3_[i-nextralayers];
      v_layerN2O_[i] = v_layerN2O_[i-nextralayers];
      v_layerNO2_[i] = v_layerNO2_[i-nextralayers];
      v_layerSO2_[i] = v_layerSO2_[i-nextralayers];
    }

    //    std::cout << "nextralayers=" << nextralayers << std::endl;

    for(int i=nextralayers; i>0 ; i--){

      v_layerThickness_[i-1] = newThickness;
      v_layerTemperature1_[i-1] = v_layerTemperature0_[i];
      v_layerTemperature0_[i-1] = v_layerTemperature0_[i]-tropoLapseRate_*0.001*v_layerThickness_[i-1];
      v_layerTemperature_[i-1] = (v_layerTemperature0_[i]+v_layerTemperature1_[i])/2.0;

      v_layerPressure1_[i-1] = v_layerPressure0_[i];
      v_layerPressure0_[i-1] = v_layerPressure0_[i]* exp(-0.0341695 *
					 (-newThickness)/ v_layerTemperature_[i-1] );
      v_layerPressure_[i-1] = exp((log(v_layerPressure0_[i-1])+ log(v_layerPressure1_[i-1]))/2.0);

      v_layerWaterVapor1_[i-1] = v_layerWaterVapor0_[i];
      v_layerWaterVapor0_[i-1] = v_layerWaterVapor0_[i]*(v_layerPressure0_[i-1]/v_layerPressure0_[i]);
      v_layerWaterVapor_[i-1] = exp((log(v_layerWaterVapor0_[i-1])+ log(v_layerWaterVapor1_[i-1]))/2.0);

      v_layerCO_[i-1] = v_layerCO_[i]*(v_layerPressure0_[i-1]/v_layerPressure0_[i]);
      v_layerO3_[i-1] = v_layerO3_[i]*(v_layerPressure0_[i-1]/v_layerPressure0_[i]);
      v_layerN2O_[i-1] = v_layerN2O_[i]*(v_layerPressure0_[i-1]/v_layerPressure0_[i]);
      v_layerNO2_[i-1] = v_layerNO2_[i]*(v_layerPressure0_[i-1]/v_layerPressure0_[i]);
      v_layerSO2_[i-1] = v_layerSO2_[i]*(v_layerPressure0_[i-1]/v_layerPressure0_[i]);

    }

    tropoLayer_ = tropoLayer_ + (v_layerThickness_.size()-numLayer_);
    numLayer_=v_layerThickness_.size();

  }else{
    //    std::cout << "extrapolar para arriba" << std::endl;
    double extraheight = fabs(groundaltitude.get(Length::UnitMeter)-altitude_.get(Length::UnitMeter));
    double cumulheight = 0.0;

    // std::cout << "before v_layerTemperature_.size()=" << v_layerTemperature_.size() << std::endl;

    for(unsigned int i=0; i<numLayer_; i++){
      cumulheight = cumulheight +  v_layerThickness_[i];
      // std::cout << i << " cumulheight = " << cumulheight <<  " extraheight= " << extraheight << std::endl;
      if (cumulheight>=extraheight) {

	double remainingheight = fabs(cumulheight-extraheight);
	// std::cout << "remainingheight = " <<  remainingheight << " m" << std::endl;
	v_layerThickness_.erase(v_layerThickness_.begin(),v_layerThickness_.begin()+i);
	v_layerThickness_[0] = remainingheight;


	v_layerTemperature0_.erase(v_layerTemperature0_.begin(),v_layerTemperature0_.begin()+i);
      	v_layerTemperature1_.erase(v_layerTemperature1_.begin(),v_layerTemperature1_.begin()+i);
	v_layerTemperature_.erase(v_layerTemperature_.begin(),v_layerTemperature_.begin()+i);
        v_layerTemperature0_[0] = v_layerTemperature0_[1]-tropoLapseRate_*0.001*v_layerThickness_[0];
        v_layerTemperature_[0] = (v_layerTemperature0_[0]+v_layerTemperature1_[0])/2.0;

	v_layerPressure0_.erase(v_layerPressure0_.begin(),v_layerPressure0_.begin()+i);
	// std::cout << "v_layerPressure1_[0]=" << v_layerPressure1_[0] << " remainingheight=" << remainingheight << std::endl;
	v_layerPressure1_.erase(v_layerPressure1_.begin(),v_layerPressure1_.begin()+i);
	v_layerPressure_.erase(v_layerPressure_.begin(),v_layerPressure_.begin()+i);
	v_layerPressure0_[0] = v_layerPressure1_[0]* exp(-0.0341695 *
					 (-remainingheight)/ v_layerTemperature_[0] );
	v_layerPressure_[0] = exp((log(v_layerPressure0_[0])+ log(v_layerPressure1_[0]))/2.0);

	v_layerWaterVapor0_.erase(v_layerWaterVapor0_.begin(),v_layerWaterVapor0_.begin()+i);
	v_layerWaterVapor1_.erase(v_layerWaterVapor1_.begin(),v_layerWaterVapor1_.begin()+i);
	v_layerWaterVapor_.erase(v_layerWaterVapor_.begin(),v_layerWaterVapor_.begin()+i);
	v_layerWaterVapor0_[0] = v_layerWaterVapor1_[0]*(v_layerPressure0_[0]/v_layerPressure1_[0]);
	v_layerWaterVapor_[0] = exp((log(v_layerWaterVapor0_[0])+ log(v_layerWaterVapor1_[0]))/2.0);

	v_layerCO_.erase(v_layerCO_.begin(),v_layerCO_.begin()+i);
	v_layerO3_.erase(v_layerO3_.begin(),v_layerO3_.begin()+i);
	v_layerN2O_.erase(v_layerN2O_.begin(),v_layerN2O_.begin()+i);
	v_layerNO2_.erase(v_layerNO2_.begin(),v_layerNO2_.begin()+i);
	v_layerSO2_.erase(v_layerSO2_.begin(),v_layerSO2_.begin()+i);

	tropoLayer_ = tropoLayer_ - i;

	// std::cout << "after v_layerTemperature_.size()=" << v_layerTemperature_.size() << std::endl;
	break;
      }
    }

    numLayer_ = v_layerThickness_.size();

  }

  altitude_ = groundaltitude;
  groundTemperature_ = v_layerTemperature0_[0];

  //  std::cout << "v_layerPressure0_[0]=" << v_layerPressure0_[0] << std::endl;

  groundPressure_ = v_layerPressure0_[0]*100.0;  // default units for Pressure are Pascals

  /*
  const Temperature Tdif(((true_antenna_altitude.get(Length::UnitKiloMeter)-altitude.get(Length::UnitKiloMeter))*tropoLapseRate),Temperature::UnitKelvin);
  groundTemperature_ = groundTemperature+Tdif;
  groundPressure_ = groundPressure * exp(-0.0341695 * pow((6.371/(6.371+altitude.get(Length::UnitKiloMeter))),2.0)
					 * (true_antenna_altitude.get(Length::UnitMeter)-altitude.get(Length::UnitMeter))
					 / ((groundTemperature.get(Temperature::UnitKelvin)+Tdif.get(Temperature::UnitKelvin))/2.0) );
  */



}

void AtmProfile::setLayerTemperature(unsigned int i, const Temperature &layerTemperature)
{
  if(i < v_layerTemperature_.size()) {
    v_layerTemperature_[i] = layerTemperature.get(Temperature::UnitKelvin);
  }
}

vector<Length> AtmProfile::getThicknessProfile() const
{
  vector<Length> l;
  l.reserve(v_layerThickness_.size());
  for(unsigned int i = 0; i < v_layerThickness_.size(); i++) {
    Length ll(v_layerThickness_[i], Length::UnitMeter);
    l.push_back(ll);
  }
  return l;
}

Length AtmProfile::getLayerThickness(unsigned int i) const
{
  /*if(i > v_layerThickness_.size() - 1) {
    Length l(-999.0, Length::UnitMeter);
    return l;
  } else {
    Length l(v_layerThickness_[i], Length::UnitMeter);
    return l;
  }*/
  if(i > v_layerThickness_.size() - 1) {
    std::ostringstream oss;
    oss << "Not a valid layer: " << i;
    throw AtmException(ATM_EXCEPTION_ARGS(oss.str().c_str()));
  }
  return Length(v_layerThickness_[i], Length::UnitMeter);
}

Length AtmProfile::getLayerBottomHeightAboveGround(unsigned int i) const
{
  if(i > v_layerThickness_.size() - 1) {
    std::ostringstream oss;
    oss << "Not a valid layer: " << i;
    throw AtmException(ATM_EXCEPTION_ARGS(oss.str().c_str()));
  }
  double h=0.0;
  for(unsigned int j = 0; j < i; j++) {
    h = h + v_layerThickness_[j];
  }
  return Length(h, Length::UnitMeter);
}

Length AtmProfile::getLayerTopHeightAboveGround(unsigned int i) const
{
  if(i > v_layerThickness_.size() - 1) {
    std::ostringstream oss;
    oss << "Not a valid layer: " << i;
    throw AtmException(ATM_EXCEPTION_ARGS(oss.str().c_str()));
  }
  double h=0.0;
  for(unsigned int j = 0; j < i+1; j++) {
    h = h + v_layerThickness_[j];
  }
  return Length(h, Length::UnitMeter);
}

Length AtmProfile::getLayerBottomHeightAboveSeaLevel(unsigned int i) const
{
  if(i > v_layerThickness_.size() - 1) {
    std::ostringstream oss;
    oss << "Not a valid layer: " << i;
    throw AtmException(ATM_EXCEPTION_ARGS(oss.str().c_str()));
  }
  double h=altitude_.get(Length::UnitMeter);
  for(unsigned int j = 0; j < i; j++) {
    h = h + v_layerThickness_[j];
  }
  return Length(h, Length::UnitMeter);
}

Length AtmProfile::getLayerTopHeightAboveSeaLevel(unsigned int i) const
{
  if(i > v_layerThickness_.size() - 1) {
    std::ostringstream oss;
    oss << "Not a valid layer: " << i;
    throw AtmException(ATM_EXCEPTION_ARGS(oss.str().c_str()));
  }
  double h=altitude_.get(Length::UnitMeter);
  for(unsigned int j = 0; j < i+1; j++) {
    h = h + v_layerThickness_[j];
  }
  return Length(h, Length::UnitMeter);
}


void AtmProfile::setLayerThickness(unsigned int i, const Length &layerThickness)
{
  if(i < v_layerThickness_.size()) {
    v_layerThickness_[i] = layerThickness.get(Length::UnitMeter);
  }
}

MassDensity AtmProfile::getLayerWaterVaporMassDensity(unsigned int i) const
{
  /*if(i > v_layerWaterVapor_.size() - 1) {
    MassDensity m(-999.0, MassDensity::UnitKiloGramPerCubicMeter);
    return m;
  } else {
    MassDensity m(v_layerWaterVapor_[i], MassDensity::UnitKiloGramPerCubicMeter);
    return m;
  }*/
  if(i > v_layerWaterVapor_.size() - 1) {
    std::ostringstream oss;
    oss << "Not a valid layer: " << i;
    throw AtmException(ATM_EXCEPTION_ARGS(oss.str().c_str()));
  }
  return MassDensity(v_layerWaterVapor_[i], MassDensity::UnitKiloGramPerCubicMeter);
}
MassDensity AtmProfile::getLayerBottomWaterVaporMassDensity(unsigned int i) const
{
  if(i > v_layerWaterVapor0_.size() - 1) {
    std::ostringstream oss;
    oss << "Not a valid layer: " << i;
    throw AtmException(ATM_EXCEPTION_ARGS(oss.str().c_str()));
  }
  return MassDensity(v_layerWaterVapor0_[i], MassDensity::UnitKiloGramPerCubicMeter);
}
MassDensity AtmProfile::getLayerTopWaterVaporMassDensity(unsigned int i) const
{
  if(i > v_layerWaterVapor1_.size() - 1) {
    std::ostringstream oss;
    oss << "Not a valid layer: " << i;
    throw AtmException(ATM_EXCEPTION_ARGS(oss.str().c_str()));
  }
  return MassDensity(v_layerWaterVapor1_[i], MassDensity::UnitKiloGramPerCubicMeter);
}

NumberDensity AtmProfile::getLayerWaterVaporNumberDensity(unsigned int i) const
{
  /*if(i > v_layerWaterVapor_.size() - 1) {
    NumberDensity m(-999.0, NumberDensity::UnitInverseCubicMeter);
    return m;
  } else {
    NumberDensity
        m(v_layerWaterVapor_[i] * 6.023e23 * 1000.0 / 18.0, NumberDensity::UnitInverseCubicMeter);
    return m;
  }*/
  if(i > v_layerWaterVapor_.size() - 1) {
    std::ostringstream oss;
    oss << "Not a valid layer: " << i;
    throw AtmException(ATM_EXCEPTION_ARGS(oss.str().c_str()));
  }
  return NumberDensity(v_layerWaterVapor_[i] * 6.023e23 * 1000.0 / 18.0, NumberDensity::UnitInverseCubicMeter);
}
NumberDensity AtmProfile::getLayerBottomWaterVaporNumberDensity(unsigned int i) const
{
  if(i > v_layerWaterVapor0_.size() - 1) {
    std::ostringstream oss;
    oss << "Not a valid layer: " << i;
    throw AtmException(ATM_EXCEPTION_ARGS(oss.str().c_str()));
  }
  return NumberDensity(v_layerWaterVapor0_[i] * 6.023e23 * 1000.0 / 18.0, NumberDensity::UnitInverseCubicMeter);
}
NumberDensity AtmProfile::getLayerTopWaterVaporNumberDensity(unsigned int i) const
{
  if(i > v_layerWaterVapor1_.size() - 1) {
    std::ostringstream oss;
    oss << "Not a valid layer: " << i;
    throw AtmException(ATM_EXCEPTION_ARGS(oss.str().c_str()));
  }
  return NumberDensity(v_layerWaterVapor1_[i] * 6.023e23 * 1000.0 / 18.0, NumberDensity::UnitInverseCubicMeter);
}

void AtmProfile::setLayerWaterVaporMassDensity(unsigned int i, const MassDensity &layerWaterVapor)
{
  if(i <= v_layerWaterVapor_.size() - 1) {
    v_layerWaterVapor_[i] = layerWaterVapor.get(MassDensity::UnitKiloGramPerCubicMeter);
  }
}

void AtmProfile::setLayerWaterVaporNumberDensity(unsigned int i, const NumberDensity &layerWaterVapor)
{
  if(i <= v_layerWaterVapor_.size() - 1) {
    v_layerWaterVapor_[i] = layerWaterVapor.get(NumberDensity::UnitInverseCubicMeter) * 18.0 / (6.023e23 * 1000.0);
  }
}

vector<Pressure> AtmProfile::getPressureProfile() const
{
  vector<Pressure> p;
  p.reserve(v_layerPressure_.size());
  for(unsigned int i = 0; i < v_layerPressure_.size(); i++) {
    Pressure pp(v_layerPressure_[i], Pressure::UnitMilliBar);
    p.push_back(pp);
  }
  return p;
}

Pressure AtmProfile::getLayerPressure(unsigned int i) const
{
  /*if(i > v_layerPressure_.size() - 1) {
    Pressure p(-999.0, Pressure::UnitMilliBar);
    return p;
  } else {
    Pressure p(v_layerPressure_[i], Pressure::UnitMilliBar);
    return p;
  }*/
  if(i > v_layerPressure_.size() - 1) {
    std::ostringstream oss;
    oss << "Not a valid layer: " << i;
    throw AtmException(ATM_EXCEPTION_ARGS(oss.str().c_str()));
  }
  return Pressure(v_layerPressure_[i], Pressure::UnitMilliBar);
}

Pressure AtmProfile::getLayerBottomPressure(unsigned int i) const
{
  if(i > v_layerPressure0_.size() - 1) {
    std::ostringstream oss;
    oss << "Not a valid layer: " << i;
    throw AtmException(ATM_EXCEPTION_ARGS(oss.str().c_str()));
  }
  return Pressure(v_layerPressure0_[i], Pressure::UnitMilliBar);
}

Pressure AtmProfile::getLayerTopPressure(unsigned int i) const
{
  if(i > v_layerPressure1_.size() - 1) {
    std::ostringstream oss;
    oss << "Not a valid layer: " << i;
    throw AtmException(ATM_EXCEPTION_ARGS(oss.str().c_str()));
  }
  return Pressure(v_layerPressure1_[i], Pressure::UnitMilliBar);
}


Length AtmProfile::getGroundWH2O() const
{
  double wm = 0;
  for(unsigned int j = 0; j < numLayer_; j++) {
    wm = wm + v_layerWaterVapor_[j] * v_layerThickness_[j]; // kg/m**2 or mm (from m*kg/m**3 IS units)
  }
  wm = wm * 1e-3; // (pasar de mm a m)
  Length wh2o(wm);
  return wh2o;
}

MassDensity AtmProfile::rwat(const Temperature &tt, const Humidity &rh, const Pressure &pp) const
{
  double t = tt.get(Temperature::UnitKelvin);
  double p = pp.get(Pressure::UnitMilliBar);
  double u = rh.get(Percent::UnitPercent);
  double e, es, rwat0;

  if(p <= 0 || t <= 0 || u <= 0) {
    return MassDensity(0.0, MassDensity::UnitGramPerCubicMeter);
  } else {
    es = 6.105 * exp(25.22 / t * (t - 273.0) - 5.31 * log(t / 273.0));
    e = 1.0 - (1.0 - u / 100.0) * es / p;
    e = es * u / 100.0 / e;
    rwat0 = e * 216.502 / t; //(en g/m*3)
  }
  return MassDensity(rwat0, MassDensity::UnitGramPerCubicMeter);
}

Humidity AtmProfile::rwat_inv(const Temperature &tt, const MassDensity &dd, const Pressure &pp) const
{
  double p = pp.get(Pressure::UnitMilliBar);
  double t = tt.get(Temperature::UnitKelvin);
  double r = dd.get(MassDensity::UnitGramPerCubicMeter);
  double es, e, rinv;

  if(p <= 0 || t <= 0 || r <= 0) {
    rinv = 0.0;
  } else {
    es = 6.105 * exp(25.22 / t * (t - 273.0) - 5.31 * log(t / 273.0));
    e = r * t / 216.502;
    rinv = 100 * (e * (p - es) / (es * (p - e)));
    if(rinv < 0 && p < 3) rinv = 0.0;
  }
  return Humidity(rinv, Percent::UnitPercent);
}

vector<NumberDensity> AtmProfile::st76(const Length &h, unsigned int tip) const
{
  unsigned int i1, i2, i3, i_layer;
  double x1, x2, x3, d;
  vector<NumberDensity> minorden;
  NumberDensity o3den, n2oden, coden, no2den, so2den;
  static const double avogad = 6.022045E+23;
  static const double airmwt = 28.964;
  //    static const double h2omwt=18.015;
  double ha = h.get(Length::UnitKiloMeter);


  static const double
      alt[50] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
		  13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0,
		  24.0, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5, 40.0, 42.5, 45.0, 47.5,
		  50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0,
		  105.0, 110.0, 115.0, 120.0 };     // REFERENCE LEVELS IN KM

  static const double
      ozone[6][50] = { { 2.869E-02, 3.150E-02, 3.342E-02, 3.504E-02, 3.561E-02, 3.767E-02, 3.989E-02, 4.223E-02,
			 4.471E-02, 5.000E-02, 5.595E-02, 6.613E-02, 7.815E-02, 9.289E-02, 1.050E-01, 1.256E-01,
			 1.444E-01, 2.500E-01, 5.000E-01, 9.500E-01, 1.400E+00, 1.800E+00, 2.400E+00, 3.400E+00,
			 4.300E+00, 5.400E+00, 7.800E+00, 9.300E+00, 9.850E+00, 9.700E+00, 8.800E+00, 7.500E+00,
			 5.900E+00, 4.500E+00, 3.450E+00, 2.800E+00, 1.800E+00, 1.100E+00, 6.500E-01, 3.000E-01,
			 1.800E-01, 3.300E-01, 5.000E-01, 5.200E-01, 5.000E-01, 4.000E-01, 2.000E-01, 5.000E-02,
			 5.000E-03, 5.000E-04 },
		       { 3.017E-02, 3.337E-02, 3.694E-02, 4.222E-02, 4.821E-02, 5.512E-02, 6.408E-02, 7.764E-02,
			 9.126E-02, 1.111E-01, 1.304E-01, 1.793E-01, 2.230E-01, 3.000E-01, 4.400E-01, 5.000E-01,
			 6.000E-01, 7.000E-01, 1.000E+00, 1.500E+00, 2.000E+00, 2.400E+00, 2.900E+00, 3.400E+00,
			 4.000E+00, 4.800E+00, 6.000E+00, 7.000E+00, 8.100E+00, 8.900E+00, 8.700E+00, 7.550E+00,
			 5.900E+00, 4.500E+00, 3.500E+00, 2.800E+00, 1.800E+00, 1.300E+00, 8.000E-01, 4.000E-01,
			 1.900E-01, 2.000E-01, 5.700E-01, 7.500E-01, 7.000E-01, 4.000E-01, 2.000E-01, 5.000E-02,
			 5.000E-03, 5.000E-04 },
		       { 2.778E-02, 2.800E-02, 2.849E-02, 3.200E-02, 3.567E-02, 4.720E-02, 5.837E-02, 7.891E-02,
			 1.039E-01, 1.567E-01, 2.370E-01, 3.624E-01, 5.232E-01, 7.036E-01, 8.000E-01, 9.000E-01,
			 1.100E+00, 1.400E+00, 1.800E+00, 2.300E+00, 2.900E+00, 3.500E+00, 3.900E+00, 4.300E+00,
			 4.700E+00, 5.100E+00, 5.600E+00, 6.100E+00, 6.800E+00, 7.100E+00, 7.200E+00, 6.900E+00,
			 5.900E+00, 4.600E+00, 3.700E+00, 2.750E+00, 1.700E+00, 1.000E-00, 5.500E-01, 3.200E-01,
			 2.500E-01, 2.300E-01, 5.500E-01, 8.000E-01, 8.000E-01, 4.000E-01, 2.000E-01, 5.000E-02,
			 5.000E-03, 5.000E-04 },
		       { 2.412E-02, 2.940E-02, 3.379E-02, 3.887E-02, 4.478E-02, 5.328E-02, 6.564E-02, 7.738E-02,
			 9.114E-02, 1.420E-01, 1.890E-01, 3.050E-01, 4.100E-01, 5.000E-01, 6.000E-01, 7.000E-01,
			 8.500E-01, 1.000E+00, 1.300E+00, 1.700E+00, 2.100E+00, 2.700E+00, 3.300E+00, 3.700E+00,
			 4.200E+00, 4.500E+00, 5.300E+00, 5.700E+00, 6.900E+00, 7.700E+00, 7.800E+00, 7.000E+00,
			 5.400E+00, 4.200E+00, 3.200E+00, 2.500E+00, 1.700E+00, 1.200E+00, 8.000E-01, 4.000E-01,
			 2.000E-01, 1.800E-01, 6.500E-01, 9.000E-01, 8.000E-01, 4.000E-01, 2.000E-01, 5.000E-02,
			 5.000E-03, 5.000E-04 },
		       { 1.802E-02, 2.072E-02, 2.336E-02, 2.767E-02, 3.253E-02, 3.801E-02, 4.446E-02, 7.252E-02,
			 1.040E-01, 2.100E-01, 3.000E-01, 3.500E-01, 4.000E-01, 6.500E-01, 9.000E-01, 1.200E+00,
			 1.500E+00, 1.900E+00, 2.450E+00, 3.100E+00, 3.700E+00, 4.000E+00, 4.200E+00, 4.500E+00,
			 4.600E+00, 4.700E+00, 4.900E+00, 5.400E+00, 5.900E+00, 6.200E+00, 6.250E+00, 5.900E+00,
			 5.100E+00, 4.100E+00, 3.000E+00, 2.600E+00, 1.600E+00, 9.500E-01, 6.500E-01, 5.000E-01,
			 3.300E-01, 1.300E-01, 7.500E-01, 8.000E-01, 8.000E-01, 4.000E-01, 2.000E-01, 5.000E-02,
			 5.000E-03, 5.000E-04 },
		       { 2.660E-02, 2.931E-02, 3.237E-02, 3.318E-02, 3.387E-02, 3.768E-02, 4.112E-02, 5.009E-02,
			 5.966E-02, 9.168E-02, 1.313E-01, 2.149E-01, 3.095E-01, 3.846E-01, 5.030E-01, 6.505E-01,
			 8.701E-01, 1.187E+00, 1.587E+00, 2.030E+00, 2.579E+00, 3.028E+00, 3.647E+00, 4.168E+00,
			 4.627E+00, 5.118E+00, 5.803E+00, 6.553E+00, 7.373E+00, 7.837E+00, 7.800E+00, 7.300E+00,
			 6.200E+00, 5.250E+00, 4.100E+00, 3.100E+00, 1.800E+00, 1.100E+00, 7.000E-01, 3.000E-01,
			 2.500E-01, 3.000E-01, 5.000E-01, 7.000E-01, 7.000E-01, 4.000E-01, 2.000E-01, 5.000E-02,
			 5.000E-03, 5.000E-04 } };

  static const double
      n2o[6][50] = { { 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01,
		       3.200E-01, 3.195E-01, 3.179E-01, 3.140E-01, 3.095E-01, 3.048E-01, 2.999E-01, 2.944E-01,
		       2.877E-01, 2.783E-01, 2.671E-01, 2.527E-01, 2.365E-01, 2.194E-01, 2.051E-01, 1.967E-01,
		       1.875E-01, 1.756E-01, 1.588E-01, 1.416E-01, 1.165E-01, 9.275E-02, 6.693E-02, 4.513E-02,
		       2.751E-02, 1.591E-02, 9.378E-03, 4.752E-03, 3.000E-03, 2.065E-03, 1.507E-03, 1.149E-03,
		       8.890E-04, 7.056E-04, 5.716E-04, 4.708E-04, 3.932E-04, 3.323E-04, 2.837E-04, 2.443E-04,
		       2.120E-04, 1.851E-04 },
		     { 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01,
		       3.195E-01, 3.163E-01, 3.096E-01, 2.989E-01, 2.936E-01, 2.860E-01, 2.800E-01, 2.724E-01,
		       2.611E-01, 2.421E-01, 2.174E-01, 1.843E-01, 1.607E-01, 1.323E-01, 1.146E-01, 1.035E-01,
		       9.622E-02, 8.958E-02, 8.006E-02, 6.698E-02, 4.958E-02, 3.695E-02, 2.519E-02, 1.736E-02,
		       1.158E-02, 7.665E-03, 5.321E-03, 3.215E-03, 2.030E-03, 1.397E-03, 1.020E-03, 7.772E-04,
		       6.257E-04, 5.166E-04, 4.352E-04, 3.727E-04, 3.237E-04, 2.844E-04, 2.524E-04, 2.260E-04,
		       2.039E-04, 1.851E-04 },
		     { 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01,
		       3.195E-01, 3.163E-01, 3.096E-01, 2.989E-01, 2.936E-01, 2.860E-01, 2.800E-01, 2.724E-01,
		       2.611E-01, 2.421E-01, 2.174E-01, 1.843E-01, 1.621E-01, 1.362E-01, 1.230E-01, 1.124E-01,
		       1.048E-01, 9.661E-02, 8.693E-02, 7.524E-02, 6.126E-02, 5.116E-02, 3.968E-02, 2.995E-02,
		       2.080E-02, 1.311E-02, 8.071E-03, 4.164E-03, 2.629E-03, 1.809E-03, 1.321E-03, 1.007E-03,
		       7.883E-04, 6.333E-04, 5.194E-04, 4.333E-04, 3.666E-04, 3.140E-04, 2.717E-04, 2.373E-04,
		       2.089E-04, 1.851E-04 },
		     { 3.100E-01, 3.100E-01, 3.100E-01, 3.100E-01, 3.079E-01, 3.024E-01, 2.906E-01, 2.822E-01,
		       2.759E-01, 2.703E-01, 2.651E-01, 2.600E-01, 2.549E-01, 2.494E-01, 2.433E-01, 2.355E-01,
		       2.282E-01, 2.179E-01, 2.035E-01, 1.817E-01, 1.567E-01, 1.350E-01, 1.218E-01, 1.102E-01,
		       9.893E-02, 8.775E-02, 7.327E-02, 5.941E-02, 4.154E-02, 3.032E-02, 1.949E-02, 1.274E-02,
		       9.001E-03, 6.286E-03, 4.558E-03, 2.795E-03, 1.765E-03, 1.214E-03, 8.866E-04, 6.756E-04,
		       5.538E-04, 4.649E-04, 3.979E-04, 3.459E-04, 3.047E-04, 2.713E-04, 2.439E-04, 2.210E-04,
		       2.017E-04, 1.851E-04 },
		     { 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01,
		       3.195E-01, 3.163E-01, 3.096E-01, 2.989E-01, 2.936E-01, 2.860E-01, 2.800E-01, 2.724E-01,
		       2.611E-01, 2.421E-01, 2.174E-01, 1.843E-01, 1.621E-01, 1.362E-01, 1.230E-01, 1.122E-01,
		       1.043E-01, 9.570E-02, 8.598E-02, 7.314E-02, 5.710E-02, 4.670E-02, 3.439E-02, 2.471E-02,
		       1.631E-02, 1.066E-02, 7.064E-03, 3.972E-03, 2.508E-03, 1.726E-03, 1.260E-03, 9.602E-04,
		       7.554E-04, 6.097E-04, 5.024E-04, 4.210E-04, 3.579E-04, 3.080E-04, 2.678E-04, 2.350E-04,
		       2.079E-04, 1.851E-04 },
		     { 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01,
		       3.200E-01, 3.195E-01, 3.179E-01, 3.140E-01, 3.095E-01, 3.048E-01, 2.999E-01, 2.944E-01,
		       2.877E-01, 2.783E-01, 2.671E-01, 2.527E-01, 2.365E-01, 2.194E-01, 2.051E-01, 1.967E-01,
		       1.875E-01, 1.756E-01, 1.588E-01, 1.416E-01, 1.165E-01, 9.275E-02, 6.693E-02, 4.513E-02,
		       2.751E-02, 1.591E-02, 9.378E-03, 4.752E-03, 3.000E-03, 2.065E-03, 1.507E-03, 1.149E-03,
		       8.890E-04, 7.056E-04, 5.716E-04, 4.708E-04, 3.932E-04, 3.323E-04, 2.837E-04, 2.443E-04,
		       2.120E-04, 1.851E-04 } };

  static const double
      co[6][50] = { { 1.500E-01, 1.450E-01, 1.399E-01, 1.349E-01, 1.312E-01, 1.303E-01, 1.288E-01, 1.247E-01,
		      1.185E-01, 1.094E-01, 9.962E-02, 8.964E-02, 7.814E-02, 6.374E-02, 5.025E-02, 3.941E-02,
		      3.069E-02, 2.489E-02, 1.966E-02, 1.549E-02, 1.331E-02, 1.232E-02, 1.232E-02, 1.307E-02,
		      1.400E-02, 1.521E-02, 1.722E-02, 1.995E-02, 2.266E-02, 2.487E-02, 2.738E-02, 3.098E-02,
		      3.510E-02, 3.987E-02, 4.482E-02, 5.092E-02, 5.985E-02, 6.960E-02, 9.188E-02, 1.938E-01,
		      5.688E-01, 1.549E+00, 3.849E+00, 6.590E+00, 1.044E+01, 1.705E+01, 2.471E+01, 3.358E+01,
		      4.148E+01, 5.000E+01 },
		    { 1.500E-01, 1.450E-01, 1.399E-01, 1.349E-01, 1.312E-01, 1.303E-01, 1.288E-01, 1.247E-01,
		      1.185E-01, 1.094E-01, 9.962E-02, 8.964E-02, 7.814E-02, 6.374E-02, 5.025E-02, 3.941E-02,
		      3.069E-02, 2.489E-02, 1.966E-02, 1.549E-02, 1.331E-02, 1.232E-02, 1.232E-02, 1.307E-02,
		      1.400E-02, 1.521E-02, 1.722E-02, 1.995E-02, 2.266E-02, 2.487E-02, 2.716E-02, 2.962E-02,
		      3.138E-02, 3.307E-02, 3.487E-02, 3.645E-02, 3.923E-02, 4.673E-02, 6.404E-02, 1.177E-01,
		      2.935E-01, 6.815E-01, 1.465E+00, 2.849E+00, 5.166E+00, 1.008E+01, 1.865E+01, 2.863E+01,
		      3.890E+01, 5.000E+01 },
		    { 1.500E-01, 1.450E-01, 1.399E-01, 1.349E-01, 1.312E-01, 1.303E-01, 1.288E-01, 1.247E-01,
		      1.185E-01, 1.094E-01, 9.962E-02, 8.964E-02, 7.814E-02, 6.374E-02, 5.025E-02, 3.941E-02,
		      3.069E-02, 2.489E-02, 1.966E-02, 1.549E-02, 1.331E-02, 1.232E-02, 1.232E-02, 1.307E-02,
		      1.400E-02, 1.498E-02, 1.598E-02, 1.710E-02, 1.850E-02, 1.997E-02, 2.147E-02, 2.331E-02,
		      2.622E-02, 3.057E-02, 3.803E-02, 6.245E-02, 1.480E-01, 2.926E-01, 5.586E-01, 1.078E+00,
		      1.897E+00, 2.960E+00, 4.526E+00, 6.862E+00, 1.054E+01, 1.709E+01, 2.473E+01, 3.359E+01,
		      4.149E+01, 5.000E+01 },
		    { 1.500E-01, 1.450E-01, 1.399E-01, 1.349E-01, 1.312E-01, 1.303E-01, 1.288E-01, 1.247E-01,
		      1.185E-01, 1.094E-01, 9.962E-02, 8.964E-02, 7.814E-02, 6.374E-02, 5.025E-02, 3.941E-02,
		      3.069E-02, 2.489E-02, 1.966E-02, 1.549E-02, 1.331E-02, 1.232E-02, 1.232E-02, 1.307E-02,
		      1.400E-02, 1.510E-02, 1.649E-02, 1.808E-02, 1.997E-02, 2.183E-02, 2.343E-02, 2.496E-02,
		      2.647E-02, 2.809E-02, 2.999E-02, 3.220E-02, 3.650E-02, 4.589E-02, 6.375E-02, 1.176E-01,
		      3.033E-01, 7.894E-01, 1.823E+00, 3.402E+00, 5.916E+00, 1.043E+01, 1.881E+01, 2.869E+01,
		      3.892E+01, 5.000E+01 },
		    { 1.500E-01, 1.450E-01, 1.399E-01, 1.349E-01, 1.312E-01, 1.303E-01, 1.288E-01, 1.247E-01,
		      1.185E-01, 1.094E-01, 9.962E-02, 8.964E-02, 7.814E-02, 6.374E-02, 5.025E-02, 3.941E-02,
		      3.069E-02, 2.489E-02, 1.966E-02, 1.549E-02, 1.331E-02, 1.232E-02, 1.232E-02, 1.307E-02,
		      1.400E-02, 1.521E-02, 1.722E-02, 2.037E-02, 2.486E-02, 3.168E-02, 4.429E-02, 6.472E-02,
		      1.041E-01, 1.507E-01, 2.163E-01, 3.141E-01, 4.842E-01, 7.147E-01, 1.067E+00, 1.516E+00,
		      2.166E+00, 3.060E+00, 4.564E+00, 6.877E+00, 1.055E+01, 1.710E+01, 2.473E+01, 3.359E+01,
		      4.149E+01, 5.000E+01 },
		    { 1.500E-01, 1.450E-01, 1.399E-01, 1.349E-01, 1.312E-01, 1.303E-01, 1.288E-01, 1.247E-01,
		      1.185E-01, 1.094E-01, 9.962E-02, 8.964E-02, 7.814E-02, 6.374E-02, 5.025E-02, 3.941E-02,
		      3.069E-02, 2.489E-02, 1.966E-02, 1.549E-02, 1.331E-02, 1.232E-02, 1.232E-02, 1.307E-02,
		      1.400E-02, 1.498E-02, 1.598E-02, 1.710E-02, 1.850E-02, 2.009E-02, 2.220E-02, 2.497E-02,
		      2.824E-02, 3.241E-02, 3.717E-02, 4.597E-02, 6.639E-02, 1.073E-01, 1.862E-01, 3.059E-01,
		      6.375E-01, 1.497E+00, 3.239E+00, 5.843E+00, 1.013E+01, 1.692E+01, 2.467E+01, 3.356E+01,
		      4.148E+01, 5.000E+01 } };

  static const double
      den[6][50] = { { 2.450E+19, 2.231E+19, 2.028E+19, 1.827E+19, 1.656E+19, 1.499E+19, 1.353E+19, 1.218E+19,
		       1.095E+19, 9.789E+18, 8.747E+18, 7.780E+18, 6.904E+18, 6.079E+18, 5.377E+18, 4.697E+18,
		       4.084E+18, 3.486E+18, 2.877E+18, 2.381E+18, 1.981E+18, 1.651E+18, 1.381E+18, 1.169E+18,
		       9.920E+17, 8.413E+17, 5.629E+17, 3.807E+17, 2.598E+17, 1.789E+17, 1.243E+17, 8.703E+16,
		       6.147E+16, 4.352E+16, 3.119E+16, 2.291E+16, 1.255E+16, 6.844E+15, 3.716E+15, 1.920E+15,
		       9.338E+14, 4.314E+14, 1.801E+14, 7.043E+13, 2.706E+13, 1.098E+13, 4.445E+12, 1.941E+12,
		       8.706E+11, 4.225E+11 },
		     { .496E+19, 2.257E+19, 2.038E+19, 1.843E+19, 1.666E+19, 1.503E+19, 1.351E+19, 1.212E+19,
		       1.086E+19, 9.716E+18, 8.656E+18, 7.698E+18, 6.814E+18, 6.012E+18, 5.141E+18, 4.368E+18,
		       3.730E+18, 3.192E+18, 2.715E+18, 2.312E+18, 1.967E+18, 1.677E+18, 1.429E+18, 1.223E+18,
		       1.042E+18, 8.919E+17, 6.050E+17, 4.094E+17, 2.820E+17, 1.927E+17, 1.338E+17, 9.373E+16,
		       6.624E+16, 4.726E+16, 3.398E+16, 2.500E+16, 1.386E+16, 7.668E+15, 4.196E+15, 2.227E+15,
		       1.109E+15, 4.996E+14, 1.967E+14, 7.204E+13, 2.541E+13, 9.816E+12, 3.816E+12, 1.688E+12,
		       8.145E+11, 4.330E+11 },
		     { .711E+19, 2.420E+19, 2.158E+19, 1.922E+19, 1.724E+19, 1.542E+19, 1.376E+19, 1.225E+19,
		       1.086E+19, 9.612E+18, 8.472E+18, 7.271E+18, 6.237E+18, 5.351E+18, 4.588E+18, 3.931E+18,
		       3.368E+18, 2.886E+18, 2.473E+18, 2.115E+18, 1.809E+18, 1.543E+18, 1.317E+18, 1.125E+18,
		       9.633E+17, 8.218E+17, 5.536E+17, 3.701E+17, 2.486E+17, 1.647E+17, 1.108E+17, 7.540E+16,
		       5.202E+16, 3.617E+16, 2.570E+16, 1.863E+16, 1.007E+16, 5.433E+15, 2.858E+15, 1.477E+15,
		       7.301E+14, 3.553E+14, 1.654E+14, 7.194E+13, 3.052E+13, 1.351E+13, 6.114E+12, 2.952E+12,
		       1.479E+12, 7.836E+11 },
		     { .549E+19, 2.305E+19, 2.080E+19, 1.873E+19, 1.682E+19, 1.508E+19, 1.357E+19, 1.216E+19,
		       1.088E+19, 9.701E+18, 8.616E+18, 7.402E+18, 6.363E+18, 5.471E+18, 4.699E+18, 4.055E+18,
		       3.476E+18, 2.987E+18, 2.568E+18, 2.208E+18, 1.899E+18, 1.632E+18, 1.403E+18, 1.207E+18,
		       1.033E+18, 8.834E+17, 6.034E+17, 4.131E+17, 2.839E+17, 1.938E+17, 1.344E+17, 9.402E+16,
		       6.670E+16, 4.821E+16, 3.516E+16, 2.581E+16, 1.421E+16, 7.946E+15, 4.445E+15, 2.376E+15,
		       1.198E+15, 5.311E+14, 2.022E+14, 7.221E+13, 2.484E+13, 9.441E+12, 3.624E+12, 1.610E+12,
		       7.951E+11, 4.311E+11 },
		     { .855E+19, 2.484E+19, 2.202E+19, 1.950E+19, 1.736E+19, 1.552E+19, 1.383E+19, 1.229E+19,
		       1.087E+19, 9.440E+18, 8.069E+18, 6.898E+18, 5.893E+18, 5.039E+18, 4.308E+18, 3.681E+18,
		       3.156E+18, 2.704E+18, 2.316E+18, 1.982E+18, 1.697E+18, 1.451E+18, 1.241E+18, 1.061E+18,
		       9.065E+17, 7.742E+17, 5.134E+17, 3.423E+17, 2.292E+17, 1.533E+17, 1.025E+17, 6.927E+16,
		       4.726E+16, 3.266E+16, 2.261E+16, 1.599E+16, 8.364E+15, 4.478E+15, 2.305E+15, 1.181E+15,
		       6.176E+14, 3.127E+14, 1.531E+14, 7.244E+13, 3.116E+13, 1.403E+13, 6.412E+12, 3.099E+12,
		       1.507E+12, 7.814E+11 },
		     { .548E+19, 2.313E+19, 2.094E+19, 1.891E+19, 1.704E+19, 1.532E+19, 1.373E+19, 1.228E+19,
		       1.094E+19, 9.719E+18, 8.602E+18, 7.589E+18, 6.489E+18, 5.546E+18, 4.739E+18, 4.050E+18,
		       3.462E+18, 2.960E+18, 2.530E+18, 2.163E+18, 1.849E+18, 1.575E+18, 1.342E+18, 1.144E+18,
		       9.765E+17, 8.337E+17, 5.640E+17, 3.830E+17, 2.524E+17, 1.761E+17, 1.238E+17, 8.310E+16,
		       5.803E+16, 4.090E+16, 2.920E+16, 2.136E+16, 1.181E+16, 6.426E+15, 3.386E+15, 1.723E+15,
		       8.347E+14, 3.832E+14, 1.711E+14, 7.136E+13, 2.924E+13, 1.189E+13, 5.033E+12, 2.144E+12,
		       9.688E+11, 5.114E+11 } };

 static const double so2[50] = { 3.00E-04,  2.74E-04,  2.36E-04,  1.90E-04,  1.46E-04,
				 1.18E-04,  9.71E-05,  8.30E-05,  7.21E-05,  6.56E-05,
				 6.08E-05,  5.79E-05,  5.60E-05,  5.59E-05,  5.64E-05,
				 5.75E-05,  5.75E-05,  5.37E-05,  4.78E-05,  3.97E-05,
				 3.19E-05,  2.67E-05,  2.28E-05,  2.07E-05,  1.90E-05,
				 1.75E-05,  1.54E-05,  1.34E-05,  1.21E-05,  1.16E-05,
				 1.21E-05,  1.36E-05,  1.65E-05,  2.10E-05,  2.77E-05,
				 3.56E-05,  4.59E-05,  5.15E-05,  5.11E-05,  4.32E-05,
				 2.83E-05,  1.33E-05,  5.56E-06,  2.24E-06,  8.96E-07,
				 3.58E-07,  1.43E-07,  5.73E-08,  2.29E-08,  9.17E-09};

 static const double no2[50] = { 2.30E-05,  2.30E-05,  2.30E-05,  2.30E-05,  2.30E-05,
				 2.30E-05,  2.30E-05,  2.30E-05,  2.30E-05,  2.32E-05,
				 2.38E-05,  2.62E-05,  3.15E-05,  4.45E-05,  7.48E-05,
				 1.71E-04,  3.19E-04,  5.19E-04,  7.71E-04,  1.06E-03,
				 1.39E-03,  1.76E-03,  2.16E-03,  2.58E-03,  3.06E-03,
				 3.74E-03,  4.81E-03,  6.16E-03,  7.21E-03,  7.28E-03,
				 6.26E-03,  4.03E-03,  2.17E-03,  1.15E-03,  6.66E-04,
				 4.43E-04,  3.39E-04,  2.85E-04,  2.53E-04,  2.31E-04,
				 2.15E-04,  2.02E-04,  1.92E-04,  1.83E-04,  1.76E-04,
				 1.70E-04,  1.64E-04,  1.59E-04,  1.55E-04,  1.51E-04};

 //   --------------------------------------------[ altitude for interpolation ]

  if((ha < 0.0) || (ha > 120.0)) {

    o3den = NumberDensity(0.0, NumberDensity::UnitInverseCubicMeter);
    n2oden = NumberDensity(0.0, NumberDensity::UnitInverseCubicMeter);
    coden = NumberDensity(0.0, NumberDensity::UnitInverseCubicMeter);
    no2den = NumberDensity(0.0, NumberDensity::UnitInverseCubicMeter);
    so2den = NumberDensity(0.0, NumberDensity::UnitInverseCubicMeter);

  } else {

    i1 = 0; i2 = 0; i3 = 0;
    x1 = 0.0; x2 = 0.0; x3 = 0.0;

    for(i_layer = 0; i_layer < 50; i_layer++) {

      if(ha < alt[i_layer]) {

        if(i_layer == 49) {

          i1 = i_layer - 2; i2 = i_layer - 1; i3 = i_layer;

        } else {

          if(i_layer == 0) {

            i1 = i_layer; i2 = i_layer + 1; i3 = i_layer + 2;

          } else {

            i1 = i_layer - 1; i2 = i_layer; i3 = i_layer + 1;

          }

        }

        x1 = alt[i1]; x2 = alt[i2]; x3 = alt[i3];
        goto calc;
      }
    }

    calc:

    if(x1 == 0.0 && x2 == 0.0 && x3 == 0.0) {

      o3den = NumberDensity(0.0, NumberDensity::UnitInverseCubicMeter);
      n2oden = NumberDensity(0.0, NumberDensity::UnitInverseCubicMeter);
      coden = NumberDensity(0.0, NumberDensity::UnitInverseCubicMeter);
      no2den = NumberDensity(0.0, NumberDensity::UnitInverseCubicMeter);
      so2den = NumberDensity(0.0, NumberDensity::UnitInverseCubicMeter);

    } else {

      //     ----- ---------------------------------------[ pressure ]

      //	p=poli2(ha,x1,x2,x3,p[tip][i1],p[tip][i2],p[tip][i3]);

      //     --------------------------------------------[ TEMPERATURE ]

      //	t=poli2(ha,x1,x2,x3,t[tip][i1],t[tip][i2],t[tip][i3]);

      //     --------------------------------------------[ DENSITY ]

      d = poli2(ha,x1,x2,x3,den[tip - 1][i1],den[tip - 1][i2],den[tip - 1][i3])
	* airmwt * 1e6 / avogad; // en g/m**3

      //   D=A+B*HA+C*HA2    en cm**-3, pero dividido por 1E19   D=D*airmwt*1e6/avogad  // en g/m**3

      //     --------------------------------------------[ H2O ]

      //    rr=poli2(ha,x1,x2,x3,agua[tip][i1],agua[tip][i2],agua[tip][i3]);   // RR=A+B*HA+C*HA2  en ppmv
      //    rr=rr*1e-6*(h2omwt/airmwt)*d;  // en g/m**3

      //     --------------------------------------------[ OZONE ]

      o3den= NumberDensity(poli2(ha, x1, x2, x3, ozone[tip - 1][i1], ozone[tip- 1][i2], ozone[tip - 1][i3])
			   * 1e-12 * d * avogad / airmwt,NumberDensity::UnitInverseCubicCentiMeter);
      if(i2 < 30){         // FIX 15/11/2018
	o3den=o3den*0.82;  // FIX 15/11/2018
      }else{               // FIX 15/11/2018
	o3den=o3den*1.65;  // FIX 15/11/2018
      }                    // FIX 15/11/2018

      //   std::cout << "ha o3den " << ha << " " << o3den.get() << std::endl;

      //  OZONO=A+B*HA+C*HA2     // en ppmv

      //     --------------------------------------------[ N2O ]

      n2oden = NumberDensity(poli2(ha, x1, x2, x3, n2o[tip - 1][i1], n2o[tip - 1][i2], n2o[tip - 1][i3])
			     * 1e-12 * d * avogad / airmwt, NumberDensity::UnitInverseCubicCentiMeter);
      //  N2O=A+B*HA+C*HA2       // en ppmv

      //    --------------------------------------------[ NO2 ]

      no2den = NumberDensity(poli2(ha, x1, x2, x3, no2[i1], no2[i2], no2[i3]) * 1e-12 * d * avogad / airmwt, NumberDensity::UnitInverseCubicCentiMeter);

      //    --------------------------------------------[ SO2 ]

      so2den = NumberDensity(poli2(ha, x1, x2, x3, no2[i1], so2[i2], so2[i3]) * 1e-12 * d * avogad / airmwt, NumberDensity::UnitInverseCubicCentiMeter);

      //    --------------------------------------------[ CO ]

      coden = NumberDensity(poli2(ha,x1,x2,x3,co[tip-1][i1],co[tip-1][i2],co[tip-1][i3])
			    * 1e-12 * d * avogad / airmwt,NumberDensity::UnitInverseCubicCentiMeter);

    }

  }

  minorden.push_back(o3den);
  minorden.push_back(n2oden);
  minorden.push_back(coden);
  minorden.push_back(no2den);
  minorden.push_back(so2den);

  return minorden;
}

double AtmProfile::poli2(double ha,
                         double x1,
                         double x2,
                         double x3,
                         double y1,
                         double y2,
                         double y3) const
{
  double a, b, c;

  c = (y3 - y2) * (x2 - x1) - (y2 - y1) * (x3 - x2);
  b = (x2 - x1) * (x3 * x3 - x2 * x2) - (x2 * x2 - x1 * x1) * (x3 - x2);
  c = c / b;
  b = (y2 - y1) - c * (x2 * x2 - x1 * x1);
  b = b / (x2 - x1);
  a = y1 - c * x1 * x1 - b * x1;

  return a + b * ha + c * pow(ha, 2);
}

unsigned int AtmProfile::mkAtmProfile()
{
  static const double
      hx[20] = { 9.225, 10.225, 11.225, 12.850, 14.850, 16.850, 18.850, 22.600, 26.600, 30.600,
		 34.850, 40.850, 46.850, 52.850, 58.850, 65.100, 73.100, 81.100, 89.100, 95.600 };

  static const double
      px[6][20] = { { 0.3190E+03, 0.2768E+03, 0.2391E+03, 0.1864E+03, 0.1354E+03, 0.9613E+02, 0.6833E+02,
		      0.3726E+02, 0.2023E+02, 0.1121E+02, 0.6142E+01, 0.2732E+01, 0.1260E+01, 0.6042E+00,
		      0.2798E+00, 0.1202E+00, 0.3600E-01, 0.9162E-02, 0.2076E-02, 0.6374E-03 },
		    { 0.3139E+03, 0.2721E+03, 0.2350E+03, 0.1833E+03, 0.1332E+03, 0.9726E+02, 0.7115E+02,
		      0.3992E+02, 0.2185E+02, 0.1216E+02, 0.6680E+01, 0.2985E+01, 0.1400E+01, 0.6780E+00,
		      0.3178E+00, 0.1380E+00, 0.4163E-01, 0.9881E-02, 0.2010E-02, 0.5804E-03 },
		    { 0.2892E+03, 0.2480E+03, 0.2124E+03, 0.1649E+03, 0.1206E+03, 0.8816E+02, 0.6433E+02,
		      0.3558E+02, 0.1901E+02, 0.1014E+02, 0.5316E+01, 0.2255E+01, 0.1022E+01, 0.4814E+00,
		      0.2206E+00, 0.9455E-01, 0.3000E-01, 0.8729E-02, 0.2332E-02, 0.8164E-03 },
		    { 0.3006E+03, 0.2587E+03, 0.2223E+03, 0.1739E+03, 0.1288E+03, 0.9495E+02, 0.7018E+02,
		      0.3983E+02, 0.2199E+02, 0.1232E+02, 0.6771E+01, 0.3056E+01, 0.1452E+01, 0.7051E+00,
		      0.3353E+00, 0.1459E+00, 0.4431E-01, 0.1024E-01, 0.1985E-02, 0.5627E-03 },
		    { 0.2731E+03, 0.2335E+03, 0.1995E+03, 0.1546E+03, 0.1130E+03, 0.8252E+02, 0.6017E+02,
		      0.3314E+02, 0.1751E+02, 0.9306E+01, 0.4826E+01, 0.1988E+01, 0.8645E+00, 0.4000E+00,
		      0.1819E+00, 0.7866E-01, 0.2639E-01, 0.8264E-02, 0.2364E-02, 0.8439E-03 },
                    { 0.2978E+03, 0.2559E+03, 0.2191E+03, 0.1698E+03, 0.1240E+03, 0.9060E+02, 0.6621E+02,
                      0.3688E+02, 0.1999E+02, 0.1087E+02, 0.5862E+01, 0.2565E+01, 0.1182E+01, 0.5572E+00,
                      0.2551E+00, 0.1074E+00, 0.3224E-01, 0.8697E-02, 0.2158E-02, 0.6851E-03 } };

  static const double
      tx[6][20] = { { 0.2421E+03, 0.2354E+03, 0.2286E+03, 0.2180E+03, 0.2046E+03, 0.1951E+03, 0.2021E+03,
		      0.2160E+03, 0.2250E+03, 0.2336E+03, 0.2428E+03, 0.2558E+03, 0.2686E+03, 0.2667E+03,
		      0.2560E+03, 0.2357E+03, 0.2083E+03, 0.1826E+03, 0.1767E+03, 0.1841E+03 },
		    { 0.2403E+03, 0.2338E+03, 0.2273E+03, 0.2167E+03, 0.2157E+03, 0.2157E+03, 0.2177E+03,
		      0.2223E+03, 0.2271E+03, 0.2349E+03, 0.2448E+03, 0.2596E+03, 0.2740E+03, 0.2727E+03,
		      0.2603E+03, 0.2394E+03, 0.2045E+03, 0.1715E+03, 0.1644E+03, 0.1785E+03 },
		    { 0.2239E+03, 0.2196E+03, 0.2191E+03, 0.2183E+03, 0.2173E+03, 0.2163E+03, 0.2153E+03,
		      0.2152E+03, 0.2152E+03, 0.2178E+03, 0.2275E+03, 0.2458E+03, 0.2637E+03, 0.2633E+03,
		      0.2531E+03, 0.2407E+03, 0.2243E+03, 0.2072E+03, 0.1991E+03, 0.2091E+03 },
		    { 0.2301E+03, 0.2252E+03, 0.2252E+03, 0.2252E+03, 0.2252E+03, 0.2252E+03, 0.2252E+03,
		      0.2252E+03, 0.2298E+03, 0.2361E+03, 0.2468E+03, 0.2648E+03, 0.2756E+03, 0.2763E+03,
		      0.2662E+03, 0.2392E+03, 0.2023E+03, 0.1681E+03, 0.1609E+03, 0.1770E+03 },
		    { 0.2172E+03, 0.2172E+03, 0.2172E+03, 0.2172E+03, 0.2172E+03, 0.2161E+03, 0.2149E+03,
		      0.2126E+03, 0.2127E+03, 0.2165E+03, 0.2222E+03, 0.2368E+03, 0.2516E+03, 0.2601E+03,
		      0.2523E+03, 0.2486E+03, 0.2388E+03, 0.2215E+03, 0.2033E+03, 0.2113E+03 },
                    { 228.260, 221.838, 216.778, 216.700, 216.700, 216.700, 216.700, 219.200, 223.136,
		      227.340, 236.110, 252.746, 268.936, 265.057, 250.174, 233.026, 212.656, 196.466,
                      187.260, 189.204}   };

  double T_ground = groundTemperature_.get(Temperature::UnitKelvin); //  std::cout<<"T_ground: " << T_ground <<"K"<<endl;
  double P_ground = groundPressure_.get(Pressure::UnitMilliBar); //  std::cout<<"P_ground: " << P_ground <<"mb"<<endl;
  double rh = relativeHumidity_.get(Percent::UnitPercent); //  std::cout<<"rh:       " << rh <<"%"<<endl;
  double h0 = wvScaleHeight_.get(Length::UnitKiloMeter); //  std::cout<<"h0:       " << h0 <<"km"<<endl;
  double dp = pressureStep_.get(Pressure::UnitMilliBar); //  std::cout<<"dp:       " << dp <<"mb"<<endl;
  double alti = altitude_.get(Length::UnitKiloMeter); //  std::cout<<"alti:     " << alti <<"km"<<endl;
  double atmh = topAtmProfile_.get(Length::UnitKiloMeter); //  std::cout<<"atmh:     " << atmh <<"km"<<endl;
  double dp1 = pressureStepFactor_; //  std::cout<<"dp1:      " << dp1 <<" "<<endl;
  double dt = tropoLapseRate_; // TODO implementer des unites (K/km) ici localement
  double prLimit;
  double minmin;

  if(typeAtm_ == 1) {
    prLimit = 100.0;   // 100.0 230.2;
  } else if(typeAtm_ == 2) {
    prLimit = 100.0;   // 198.0;
  } else if(typeAtm_ == 3) {
    prLimit = 100.0;   // 300.0
  } else if(typeAtm_ == 4) {
    prLimit = 100.0;   // 311.0
  } else if(typeAtm_ == 5) {
    prLimit = 100.0;   // 332.0
  } else if(typeAtm_ == 6) {
    prLimit = 100.0;   // 332.0
  } else {
    typeAtm_ = 6;      // Default: US St76 Atm.
    prLimit = 100.0;   // 250.0;
  }

  unsigned int npp = 0; // number of layers initialized

  double rt = 6371.2E+0; // Earth radius in km
  double g0 = 9.80665E+0; // Earth gravity at the surface  (m/s**2)

  // static bool first = true;  // [-Wunused_but_set_variable]

  unsigned int i;
  unsigned int i0;
  unsigned int j;
  unsigned int k;
  // double wh2o, wh2o0;  // [-Wunused_but_set_variable]
  double wgr0;
  double g = g0;
  double www, altura;
  // double abun_ozono, abun_n2o, abun_co;
  NumberDensity ozono, n2o, co;
  double wgr, dh;
// double humrel; // [-Wunused_but_set_variable]
  // bool   errcode;

  //int    nmaxLayers=40;  // FV peut etre devrions nos avoir un garde-fou au cas ou le nb de couches serait stupidement trop grand

  vector<double> v_layerPressure;
  vector<double> v_layerPressure0;
  vector<double> v_layerPressure1;
  vector<double> v_layerTemperature;
  vector<double> v_layerTemperature0;
  vector<double> v_layerTemperature1;
  vector<double> v_layerThickness;
  vector<double> v_layerWaterVapor;
  vector<double> v_layerWaterVapor0;
  vector<double> v_layerWaterVapor1;
  vector<double> v_layerO3;
  vector<double> v_layerCO;
  vector<double> v_layerN2O;
  vector<double> v_layerNO2;
  vector<double> v_layerSO2;

  vector<double> v_layerO3_aux;
  vector<double> v_layerCO_aux;
  vector<double> v_layerN2O_aux;
  vector<double> v_layerNO2_aux;
  vector<double> v_layerSO2_aux;
  vector<double> v_layerPressure_aux;
  vector<double> v_layerPressure0_aux;
  vector<double> v_layerPressure1_aux;
  vector<double> v_layerTemperature_aux;
  vector<double> v_layerTemperature0_aux;
  vector<double> v_layerTemperature1_aux;
  vector<double> v_layerThickness_aux;
  vector<double> v_layerWaterVapor_aux;
  vector<double> v_layerWaterVapor0_aux;
  vector<double> v_layerWaterVapor1_aux;
  vector<NumberDensity> minorden;

  v_layerPressure.push_back(P_ground);
  v_layerThickness.push_back(alti * 1000);
  v_layerTemperature.push_back(T_ground);

  wgr = rwat(Temperature(v_layerTemperature[0], Temperature::UnitKelvin),
             Humidity(rh, Percent::UnitPercent),
             Pressure(v_layerPressure[0], Pressure::UnitMilliBar)).get(MassDensity::UnitGramPerCubicMeter);

  // Absolute Humidity in gr/m**3 at altitude alti
  wgr0 = wgr * exp(alti / h0); // h0 in km ==> wgr0 would be the absolute humidity (gr/m**3) at sea level
  // wh2o0 = wgr0 * h0; // wh2o0 is the zenith column of water above sea level (kgr/m**2 or mm) // [-Wunused_but_set_variable]
  // wh2o = wh2o0 * exp(-alti / h0); // wh2o is the zenith column of water above alti (kgr/m**2 or mm) // [-Wunused_but_set_variable]

  v_layerWaterVapor.push_back(wgr); // in gr/m**3

  i0 = 0;
  i = 0;
  j = 0;
  npp = 0;

  //  std::cout << "layer " << i << " v_layerThickness[" << i  << "]=" << v_layerThickness[i] << " v_layerPressure[" << i << "]=" << v_layerPressure[i] << std::endl;

  bool control = true;
  //  bool control2 = true;
  while(true) {
    i++;

    //    std::cout << "layer: " << i << " " << v_layerPressure[i - 1] - dp * pow(dp1, i - 1) << " " << prLimit << std::endl;

    if( (v_layerPressure[i - 1] - dp * pow(dp1, i - 1) <= prLimit ) ) { //&& (v_layerTemperature[i - 1] <= tx[typeAtm_ - 1][0]) ) {  //  &&(fabs(prLimit-v_layerPressure[i - 1]-dp*pow(dp1,i-1))>=20.0) ) {

      if(control) {

	//std::cout << "capa " << i-1 << " temperature=" << v_layerTemperature[i - 1] << " K " <<  " tx[0]=" << tx[typeAtm_ - 1][0] << std::endl;

	//	std::cout <<  "prLimit=" << prLimit << std::endl;
	// std::cout << "aaa " <<  v_layerPressure[i - 1] - dp * pow(dp1, i - 1) << std::endl;

	/*      for(k = 0; k < 20; k++) {
          if(control2) {
            if(v_layerPressure[i - 1] - dp * pow(dp1, i - 1) >= px[typeAtm_ - 1][k]) {
              j = k;
	      std::cout << "P=" << v_layerPressure[i - 1] - dp * pow(dp1, i - 1)
			<< " prLimit=" << prLimit << " px[" << typeAtm_ - 1 << "][" << k << "]="
			<< px[typeAtm_ - 1][k] << std::endl;
              control2 = false;
            }
          }
	  } */

	minmin = 20000.0;

	for(k = 0; k < 20; k++) {

	  if( (fabs(v_layerPressure[i - 1] - dp * pow(dp1, i - 1) > 1.05*px[typeAtm_ - 1][k])) &&
	      (fabs(v_layerPressure[i - 1] - dp * pow(dp1, i - 1) - px[typeAtm_ - 1][k])) <= minmin ) {

	    j = k;

	    /*	    std::cout << "P=" << v_layerPressure[i - 1] - dp * pow(dp1, i - 1)
		      << " prLimit=" << prLimit << " px[" << typeAtm_ - 1 << "][" << k << "]="
		      << px[typeAtm_ - 1][k] << std::endl; */
	    minmin = fabs(v_layerPressure[i - 1] - dp * pow(dp1, i - 1) - px[typeAtm_ - 1][k]);
	    /* std::cout << " minmin=" << minmin << std::endl; */

	  }

	}


      }

      //      std::cout << "i,j,v_layerPressure.size()-1=" << i << "," << j << "," << v_layerPressure.size() - 1 << std::endl;

      if(i0 == 0) i0 = i - 1;

      //      std::cout << "i0,j " << i0 << "," << j << std::endl;


      if(i < v_layerPressure.size() - 1) {
        if(control) {
          int j0 = (j>0) ? j-1:0;
          v_layerPressure[i] = px[typeAtm_ - 1][j0];
          v_layerTemperature[i] = tx[typeAtm_ - 1][j0];

	  std::cout << "v_layerTemperature[i]=" << v_layerTemperature[i] << std::endl;
	  //	    - tx[typeAtm_ - 1][0] + v_layerTemperature[i0];   COMMENTED OUT 31/5/2017

	  www = v_layerWaterVapor[i - 1] / 1000.0;
          dh = 288.6948 * v_layerTemperature[i - 1] * (1.0 + 0.61 * www / 1000.0)
              * log(v_layerPressure[i - 1] / v_layerPressure[i]) / g;
          v_layerThickness[i] = v_layerThickness[i - 1] + dh;
        } else {
          j++;
          v_layerPressure[i] = px[typeAtm_ - 1][j - 1];
          v_layerTemperature[i] = tx[typeAtm_ - 1][j - 1];
	    //  - tx[typeAtm_ - 1][0] + v_layerTemperature[i0];   COMMENTED OUT 31/5/2017
	  v_layerThickness[i] = (hx[j] - hx[j - 1]) * 1000.0 + v_layerThickness[i - 1];
        }
	//        std::cout << "layer " << i << " j=" << j << " v_layerThickness[" << i << "]="  << v_layerThickness[i] << " v_layerPressure[" << i << "]=" << v_layerPressure[i]  << std::endl;
        v_layerWaterVapor[i] = wgr0 * exp(-v_layerThickness[i] / (1000.0 * h0));
      } else {
        if(control) {
          int j0 = (j>0) ? j-1:0;
          v_layerPressure.push_back(px[typeAtm_ - 1][j0]);
          v_layerTemperature.push_back(tx[typeAtm_ - 1][j0]);
	  //- tx[typeAtm_ - 1][0] + v_layerTemperature[i0]);   COMMENTED OUT 31/5/2017

	  www = v_layerWaterVapor[i - 1] / 1000.0;
          dh = 288.6948 * v_layerTemperature[i - 1] * (1.0 + 0.61 * www / 1000.0)
              * log(v_layerPressure[i - 1] / v_layerPressure[i]) / g;
          v_layerThickness.push_back(v_layerThickness[i - 1] + dh);
        } else {
          j++;
          v_layerPressure.push_back(px[typeAtm_ - 1][j - 1]);
          v_layerTemperature.push_back(tx[typeAtm_ - 1][j - 1]);
	    //  - tx[typeAtm_ - 1][0] + v_layerTemperature[i0]);   COMMENTED OUT 31/5/2017
          v_layerThickness.push_back((hx[j] - hx[j - 1]) * 1000.0 + v_layerThickness[i - 1]);
        }
        v_layerWaterVapor.push_back(wgr0 * exp(-v_layerThickness[i] / (1000.0 * h0)));
      }
      //  std::cout << "layer " << i << " j=" << j << " v_layerThickness[" << i << "]="  << v_layerThickness[i] << " v_layerPressure[" << i << "]=" << v_layerPressure[i]  << std::endl;
     if(control) {
        tropoLayer_ = i - 1;
	//	std::cout << "tropoLayer_=" << tropoLayer_ << std::endl;
        control = false;
      }
      //      std::cout << "2 layer " << i << " alt at end of layer = " <<  v_layerThickness[i] << " atmh=" << atmh*1000 << std::endl;
      if(v_layerThickness[i] > (atmh * 1000.0)) break;
    } else {

      //      std::cout << "i,j,v_layerPressure.size()-1=" << i << "," << j << "," << v_layerPressure.size() - 1 << std::endl;

      if(i > v_layerPressure.size() - 1) {
        v_layerPressure.push_back(v_layerPressure[i - 1] - dp * pow(dp1, i - 1));

        www = v_layerWaterVapor[i - 1] / 1000.0; // in kg/m**3
        g = g0 * pow(1. + ((v_layerThickness[i - 1] / 1000.0)) / rt, -2); // gravity corrected for the height
        dh = 288.6948 * v_layerTemperature[i - 1] * (1.0 + 0.61 * www / 1000.0)
            * log(v_layerPressure[i - 1] / v_layerPressure[i]) / g;
        v_layerThickness.push_back(v_layerThickness[i - 1] + dh);
	//std::cout << "layer " << i << " v_layerThickness[" << i << "]=" << v_layerThickness[i]  << " v_layerPressure[" << i << "]=" << v_layerPressure[i] << std::endl;
        v_layerTemperature.push_back(v_layerTemperature[i - 1] + dt * dh / 1000.0);
        v_layerWaterVapor.push_back(wgr0 * exp(-v_layerThickness[i] / (1000.0 * h0))); //r[i] in kgr/(m**2*1000m) [gr/m**3]
      } else {
        v_layerPressure[i] = v_layerPressure[i - 1] - dp * pow(dp1, i - 1);
        www = v_layerWaterVapor[i - 1] / 1000.0; // in kg/m**3
        g = g0 * pow(1. + ((v_layerThickness[i - 1] / 1000.0)) / rt, -2); // gravity corrected for the height
        dh = 288.6948 * v_layerTemperature[i - 1] * (1.0 + 0.61 * www / 1000.0)
            * log(v_layerPressure[i - 1] / v_layerPressure[i]) / g;
        v_layerThickness[i] = v_layerThickness[i - 1] + dh;
	//std::cout << "layer " << i << " v_layerThickness[" << i << "]=" << v_layerThickness[i]  << " v_layerPressure[" << i << "]=" << v_layerPressure[i] << std::endl;
        v_layerTemperature[i] = v_layerTemperature[i - 1] + dt * dh / 1000.0;
        v_layerWaterVapor[i] = wgr0 * exp(-v_layerThickness[i] / (1000.0 * h0)); //r[i] in kgr/(m**2*1000m) [gr/m**3]
      }

      //std::cout << "dh=" << dh << std::endl;
       if( v_layerTemperature[i] <= tx[typeAtm_ - 1][0] && v_layerPressure[i] <= px[typeAtm_ - 1][0] )
	{
	  // std::cout << "v_layerThickness[i]=" << v_layerThickness[i] <<  " " << v_layerThickness[i-1] <<std::endl;
	  // std::cout << "dh=" << dh << std::endl;
	  // std::cout << "old prLimit =" << prLimit << std::endl;
	  prLimit=v_layerPressure[i];  //14NOV2017
      	  // std::cout << "new prLimit =" << prLimit << std::endl;
	}

//       humrel = rwat_inv(Temperature(v_layerTemperature[i], Temperature::UnitKelvin),
//                         MassDensity(v_layerWaterVapor[i], MassDensity::UnitGramPerCubicMeter),
//                         Pressure(v_layerPressure[i], Pressure::UnitMilliBar)).get(Percent::UnitPercent);


      /*	cout << "layer " << i
       << " P " << v_layerPressure[i] << " prLimit " << prLimit
       << " T " << v_layerTemperature[i]
       << " Alti (layer top in m) " << v_layerThickness[i]
       << " WaterVapor " << v_layerWaterVapor[i] << std::endl;  	*/
      //      std::cout << "1 layer " << i << " alt at end of layer = " <<  v_layerThickness[i] << " atmh=" << atmh*1000 << std::endl;
      if(v_layerThickness[i] > (atmh * 1000.0)) break;
    }
    //    std::cout << "v_layerPressure[" << i << "]" << v_layerPressure[i] << "  v_layerTemperature[" << i << "] " << v_layerTemperature[i] << std::endl;
  }

  npp = i - 1;

  // std::cout << "npp=" << npp << std::endl;
  // std::cout << "tropoLayer=" << tropoLayer_ << std::endl;
  tropoTemperature_ = Temperature(v_layerTemperature[tropoLayer_], Temperature::UnitKelvin);
  tropoAltitude_ = Length(v_layerThickness[tropoLayer_], Length::UnitMeter);
  // std::cout << "tropoAltitude=" << tropoAltitude_.get(Length::UnitKiloMeter) << " km" << std::endl;
  // std::cout << "tropoTemperature=" << tropoTemperature_.get(Temperature::UnitKelvin) << " K" << std::endl;
  // std::cout << "ground Altitude=" << altitude_.get(Length::UnitKiloMeter) << " km" << std::endl;
  // std::cout << "ground Temperature=" << groundTemperature_.get(Temperature::UnitKelvin) << " K" << std::endl;
  // std::cout << "Calculated Lapse Rate=" << (tropoTemperature_.get(Temperature::UnitKelvin)-groundTemperature_.get(Temperature::UnitKelvin))/(tropoAltitude_.get(Length::UnitKiloMeter)-altitude_.get(Length::UnitKiloMeter)) << " K/km" << std::endl;


  altura = alti;

  /*
   if(first){
   v_layerO3_.reserve(npp);
   v_layerCO_.reserve(npp);
   v_layerN2O_.reserve(npp);
   v_layerNO2_.reserve(npp);
   v_layerSO2_.reserve(npp);
   }
   */


  for(unsigned int jj = 0; jj < npp; jj++) {

    v_layerTemperature0.push_back(v_layerTemperature[jj]);
    v_layerTemperature1.push_back(v_layerTemperature[jj + 1]);
    v_layerPressure0.push_back(v_layerPressure[jj]);
    v_layerPressure1.push_back(v_layerPressure[jj + 1]);
    v_layerWaterVapor0.push_back(1.0E-3*v_layerWaterVapor[jj]);
    v_layerWaterVapor1.push_back(1.0E-3*v_layerWaterVapor[jj + 1]);

      }

  for(j = 0; j < npp; j++) {

    v_layerThickness[j] = (v_layerThickness[j + 1] - v_layerThickness[j]); // in m
    altura = altura + v_layerThickness[j] / 2.0E3; // in km
    v_layerTemperature[j] = (v_layerTemperature[j + 1] + v_layerTemperature[j])            // AVERAGING MAY/30/2017
        / 2.; // in K
    v_layerPressure[j] = exp((log(v_layerPressure[j + 1])
        + log(v_layerPressure[j])) / 2.0); // in mb
    v_layerWaterVapor[j] = 1.0E-3 * exp((log(v_layerWaterVapor[j + 1]) + // in kg/m**3
        log(v_layerWaterVapor[j])) / 2.0); //  1.0E-3 ?


    //      std::cout << "type_=" << type_ << std::endl;
    //      std::cout << "typeAtm_=" << typeAtm_ << std::endl;

    unsigned int atmType = typeAtm_; // conversion in int

    //      std::cout << "going to minorden with atmType=" << atmType << std::endl;

    minorden = st76(Length(altura, Length::UnitKiloMeter), atmType);

    //      std::cout << "Ozone: " << abun_ozono << "  " << ozono.get(NumberDensity::UnitInverseCubicCentiMeter) << std::endl;
    // std::cout << "N2O  : " << abun_n2o << "  " << n2o.get(NumberDensity::UnitInverseCubicCentiMeter) << std::endl;
    // std::cout << "CO   : " << abun_co  << "  " << co.get(NumberDensity::UnitInverseCubicCentiMeter) << std::endl;


    /* std::cout << j << " " << v_layerO3_.size()-1 << std::endl; */

    /*    if(j>v_layerO3_.size()-1){
     v_layerO3_.push_back(  1.E6*abun_ozono); // in m**-3
     v_layerCO_.push_back(  1.E6*abun_co   ); // in m**-3
     v_layerN2O_.push_back( 1.E6*abun_n2o  ); // in m**-3
     v_layerNO2_.push_back( 1.E6*abun_no2  ); // in m**-3
     v_layerSO2_.push_back( 1.E6*abun_so2  ); // in m**-3
     std::cout << "uno" << std::endl;
     }else{ */

    /*    v_layerO3_[j]  = 1.E6*abun_ozono;        // in m**-3
     v_layerCO_[j]  = 1.E6*abun_co;           // in m**-3
     v_layerN2O_[j] = 1.E6*abun_n2o;          // in m**-3
     v_layerNO2_[j] = 1.E6*abun_no2;          // in m**-3
     v_layerSO2_[j] = 1.E6*abun_so2;          // in m**-3 */

    // v_layerO3.push_back(1.E6*abun_ozono);        // in m**-3
    // v_layerCO.push_back(1.E6*abun_co);           // in m**-3
    // v_layerN2O.push_back(1.E6*abun_n2o);          // in m**-3
    // v_layerNO2.push_back(1.E6*abun_no2);          // in m**-3
    // v_layerSO2.push_back(1.E6*abun_so2);          // in m**-3

    v_layerO3.push_back(1.E6 * minorden[0].get(NumberDensity::UnitInverseCubicCentiMeter)); // in m**-3
    v_layerCO.push_back(1.E6 * minorden[2].get(NumberDensity::UnitInverseCubicCentiMeter)); // in m**-3
    v_layerN2O.push_back(1.E6 * minorden[1].get(NumberDensity::UnitInverseCubicCentiMeter)); // in m**-3
    v_layerNO2.push_back(1.E6 * minorden[3].get(NumberDensity::UnitInverseCubicCentiMeter)); // in m**-3
    v_layerSO2.push_back(1.E6 * minorden[4].get(NumberDensity::UnitInverseCubicCentiMeter)); // in m**-3



    /*     } */

    altura = altura + v_layerThickness[j] / 2.0E3;

    /* std::cout << "j=" << j << "v_layerThickness.size()=" << v_layerThickness.size() <<  std::endl; */

  }

  /*  std::cout << "j=" << j << " v_layerThickness_aux.size()= " << v_layerThickness_aux.size() <<  std::endl; */

  /* if(v_layerThickness.size()>npp){
   for(j=npp; j<v_layerThickness.size(); j++){
   delete v_layerThickness[j];
   }
   } */

  for(j = 0; j < npp; j++) {

    v_layerPressure_aux.push_back(v_layerPressure[j]);
    v_layerPressure0_aux.push_back(v_layerPressure0[j]);
    v_layerPressure1_aux.push_back(v_layerPressure1[j]);
    v_layerTemperature_aux.push_back(v_layerTemperature[j]);
    v_layerTemperature0_aux.push_back(v_layerTemperature0[j]);
    v_layerTemperature1_aux.push_back(v_layerTemperature1[j]);
    v_layerThickness_aux.push_back(v_layerThickness[j]);
    v_layerWaterVapor_aux.push_back(v_layerWaterVapor[j]);
    v_layerWaterVapor0_aux.push_back(v_layerWaterVapor0[j]);
    v_layerWaterVapor1_aux.push_back(v_layerWaterVapor1[j]);
    v_layerO3_aux.push_back(v_layerO3[j]);
    v_layerCO_aux.push_back(v_layerCO[j]);
    v_layerN2O_aux.push_back(v_layerN2O[j]);
    v_layerNO2_aux.push_back(v_layerNO2[j]);
    v_layerSO2_aux.push_back(v_layerSO2[j]);

  }

  if(j > v_layerPressure_.size() - 1) { // ?????
    v_layerPressure_.reserve(npp);
    v_layerPressure0_.reserve(npp);
    v_layerPressure1_.reserve(npp);
    v_layerTemperature_.reserve(npp);
    v_layerTemperature0_.reserve(npp);
    v_layerTemperature1_.reserve(npp);
    v_layerThickness_.reserve(npp);
    v_layerWaterVapor_.reserve(npp);
    v_layerWaterVapor0_.reserve(npp);
    v_layerWaterVapor1_.reserve(npp);
    v_layerO3_.reserve(npp);
    v_layerCO_.reserve(npp);
    v_layerN2O_.reserve(npp);
    v_layerNO2_.reserve(npp);
    v_layerSO2_.reserve(npp);
  }

  v_layerPressure_ = v_layerPressure_aux;
  v_layerPressure0_ = v_layerPressure0_aux;
  v_layerPressure1_ = v_layerPressure1_aux;
  v_layerTemperature_ = v_layerTemperature_aux;
  v_layerTemperature0_ = v_layerTemperature0_aux;
  v_layerTemperature1_ = v_layerTemperature1_aux;
  v_layerThickness_ = v_layerThickness_aux;
  v_layerWaterVapor_ = v_layerWaterVapor_aux;
  v_layerWaterVapor0_ = v_layerWaterVapor0_aux;
  v_layerWaterVapor1_ = v_layerWaterVapor1_aux;
  v_layerO3_ = v_layerO3_aux;
  v_layerCO_ = v_layerCO_aux;
  v_layerN2O_ = v_layerN2O_aux;
  v_layerNO2_ = v_layerNO2_aux;
  v_layerSO2_ = v_layerSO2_aux;

  // first = false; // ?????    [-Wunused_but_set_variable]

  return npp;
}

ATM_NAMESPACE_END

