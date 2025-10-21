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
 * "@(#) $Id: ATMSkyStatus.cpp Exp $"
 *
 * who       when      what
 * --------  --------  ----------------------------------------------
 * pardo     24/03/09  created
 */

#include "ATMSkyStatus.h"

#include <iostream>
#include <math.h>

ATM_NAMESPACE_BEGIN

// Constructors

SkyStatus::SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile) :
  RefractiveIndexProfile(refractiveIndexProfile), airMass_(1.0),
      skyBackgroundTemperature_(2.73, Temperature::UnitKelvin)
{

  iniSkyStatus();

}

SkyStatus::SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
                     double airMass) :
  RefractiveIndexProfile(refractiveIndexProfile), airMass_(airMass),
      skyBackgroundTemperature_(2.73, Temperature::UnitKelvin)
{

  iniSkyStatus();

}

SkyStatus::SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
                     const Temperature &temperatureBackground) :
  RefractiveIndexProfile(refractiveIndexProfile), airMass_(1.0),
      skyBackgroundTemperature_(temperatureBackground)
{

  iniSkyStatus();

}

SkyStatus::SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
                     const Length &wh2o) :
  RefractiveIndexProfile(refractiveIndexProfile), airMass_(1.0),
      skyBackgroundTemperature_(2.73, Temperature::UnitKelvin), wh2o_user_(wh2o)
{

  iniSkyStatus();

}

SkyStatus::SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
                     const Temperature &temperatureBackground,
                     double airMass) :
  RefractiveIndexProfile(refractiveIndexProfile), airMass_(airMass),
      skyBackgroundTemperature_(temperatureBackground)
{

  iniSkyStatus();

}

SkyStatus::SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
                     double airMass,
                     const Temperature &temperatureBackground) :
  RefractiveIndexProfile(refractiveIndexProfile), airMass_(airMass),
      skyBackgroundTemperature_(temperatureBackground)
{

  iniSkyStatus();

}

SkyStatus::SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
                     const Length &wh2o,
                     double airMass) :
  RefractiveIndexProfile(refractiveIndexProfile), airMass_(airMass),
      skyBackgroundTemperature_(2.73, Temperature::UnitKelvin), wh2o_user_(wh2o)
{

  iniSkyStatus();

}

SkyStatus::SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
                     double airMass,
                     const Length &wh2o) :
  RefractiveIndexProfile(refractiveIndexProfile), airMass_(airMass),
      skyBackgroundTemperature_(2.73, Temperature::UnitKelvin), wh2o_user_(wh2o)
{

  iniSkyStatus();

}

SkyStatus::SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
                     const Length &wh2o,
                     const Temperature &temperatureBackground) :
  RefractiveIndexProfile(refractiveIndexProfile), airMass_(1.0),
      skyBackgroundTemperature_(temperatureBackground), wh2o_user_(wh2o)
{

  iniSkyStatus();

}

SkyStatus::SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
                     const Temperature &temperatureBackground,
                     const Length &wh2o) :
  RefractiveIndexProfile(refractiveIndexProfile), airMass_(1.0),
      skyBackgroundTemperature_(temperatureBackground), wh2o_user_(wh2o)
{

  iniSkyStatus();

}

SkyStatus::SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
                     const Temperature &temperatureBackground,
                     double airMass,
                     const Length &wh2o) :
  RefractiveIndexProfile(refractiveIndexProfile), airMass_(airMass),
      skyBackgroundTemperature_(temperatureBackground), wh2o_user_(wh2o)
{

  iniSkyStatus();

}

SkyStatus::SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
                     const Temperature &temperatureBackground,
                     const Length &wh2o,
                     double airMass) :
  RefractiveIndexProfile(refractiveIndexProfile), airMass_(airMass),
      skyBackgroundTemperature_(temperatureBackground), wh2o_user_(wh2o)
{

  iniSkyStatus();

}

SkyStatus::SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
                     double airMass,
                     const Temperature &temperatureBackground,
                     const Length &wh2o) :
  RefractiveIndexProfile(refractiveIndexProfile), airMass_(airMass),
      skyBackgroundTemperature_(temperatureBackground), wh2o_user_(wh2o)
{

  iniSkyStatus();

}

SkyStatus::SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
                     double airMass,
                     const Length &wh2o,
                     const Temperature &temperatureBackground) :
  RefractiveIndexProfile(refractiveIndexProfile), airMass_(airMass),
      skyBackgroundTemperature_(temperatureBackground), wh2o_user_(wh2o)
{

  iniSkyStatus();

}

SkyStatus::SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
                     const Length &wh2o,
                     const Temperature &temperatureBackground,
                     double airMass) :
  RefractiveIndexProfile(refractiveIndexProfile), airMass_(airMass),
      skyBackgroundTemperature_(temperatureBackground), wh2o_user_(wh2o)
{

  iniSkyStatus();

}

SkyStatus::SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
                     const Length &wh2o,
                     double airMass,
                     const Temperature &temperatureBackground) :
  RefractiveIndexProfile(refractiveIndexProfile), airMass_(airMass),
      skyBackgroundTemperature_(temperatureBackground), wh2o_user_(wh2o)
{

  iniSkyStatus();

}

SkyStatus::SkyStatus(const SkyStatus & a) : RefractiveIndexProfile(a)
{

  // 2015-03-02 (DB) :
  // Use constructor of base class RefractiveIndexProfile()
  // and comment the following code because it is executed in RefractiveIndexProfile() constructor

//   groundTemperature_ = a.groundTemperature_;
//   tropoLapseRate_ = a.tropoLapseRate_;
//   groundPressure_ = a.groundPressure_;
//   relativeHumidity_ = a.relativeHumidity_;
//   wvScaleHeight_ = a.wvScaleHeight_;
//   pressureStep_ = a.pressureStep_;
//   pressureStepFactor_ = a.pressureStepFactor_;
//   altitude_ = a.altitude_;
//   topAtmProfile_ = a.topAtmProfile_;
//   numLayer_ = a.numLayer_;
//   newBasicParam_ = a.newBasicParam_;
//   v_layerThickness_.reserve(numLayer_);
//   v_layerTemperature_.reserve(numLayer_);
//   v_layerWaterVapor_.reserve(numLayer_);
//   v_layerCO_.reserve(numLayer_);
//   v_layerO3_.reserve(numLayer_);
//   v_layerN2O_.reserve(numLayer_);
//   v_layerNO2_.reserve(numLayer_);
//   v_layerSO2_.reserve(numLayer_);
//   for(unsigned int n = 0; n < numLayer_; n++) {
//     v_layerThickness_.push_back(a.v_layerThickness_[n]);
//     v_layerTemperature_.push_back(a.v_layerTemperature_[n]);
//     //v_layerDeltaT_.push_back(a.v_layerDeltaT_[n]);
//     //cout << "n=" << n << endl;
//     v_layerWaterVapor_.push_back(a.v_layerWaterVapor_[n]);
//     v_layerPressure_.push_back(a.v_layerPressure_[n]);
//     v_layerCO_.push_back(a.v_layerCO_[n]);
//     v_layerO3_.push_back(a.v_layerO3_[n]);
//     v_layerN2O_.push_back(a.v_layerN2O_[n]);
//     v_layerNO2_.push_back(a.v_layerNO2_[n]);
//     v_layerSO2_.push_back(a.v_layerSO2_[n]);
//   }

//   // level Spectral Grid
//   freqUnits_ = a.freqUnits_;
//   v_chanFreq_ = a.v_chanFreq_;

//   v_numChan_ = a.v_numChan_;
//   v_refChan_ = a.v_refChan_;
//   v_refFreq_ = a.v_refFreq_;
//   v_chanSep_ = a.v_chanSep_;
//   v_maxFreq_ = a.v_maxFreq_;
//   v_minFreq_ = a.v_minFreq_;
//   v_intermediateFrequency_ = a.v_intermediateFrequency_;
//   v_loFreq_ = a.v_loFreq_;

//   v_sidebandSide_ = a.v_sidebandSide_;
//   v_sidebandType_ = a.v_sidebandType_;

//   vv_assocSpwId_ = a.vv_assocSpwId_;
//   vv_assocNature_ = a.vv_assocNature_;

//   v_transfertId_ = a.v_transfertId_;

//   // level Absorption Profile
//   vv_N_H2OLinesPtr_.reserve(a.v_chanFreq_.size());
//   vv_N_H2OContPtr_.reserve(a.v_chanFreq_.size());
//   vv_N_O2LinesPtr_.reserve(a.v_chanFreq_.size());
//   vv_N_DryContPtr_.reserve(a.v_chanFreq_.size());
//   vv_N_O3LinesPtr_.reserve(a.v_chanFreq_.size());
//   vv_N_COLinesPtr_.reserve(a.v_chanFreq_.size());
//   vv_N_N2OLinesPtr_.reserve(a.v_chanFreq_.size());
//   vv_N_NO2LinesPtr_.reserve(a.v_chanFreq_.size());
//   vv_N_SO2LinesPtr_.reserve(a.v_chanFreq_.size());

//   vector<complex<double> >* v_N_H2OLinesPtr;
//   vector<complex<double> >* v_N_H2OContPtr;
//   vector<complex<double> >* v_N_O2LinesPtr;
//   vector<complex<double> >* v_N_DryContPtr;
//   vector<complex<double> >* v_N_O3LinesPtr;
//   vector<complex<double> >* v_N_COLinesPtr;
//   vector<complex<double> >* v_N_N2OLinesPtr;
//   vector<complex<double> >* v_N_NO2LinesPtr;
//   vector<complex<double> >* v_N_SO2LinesPtr;

//   for(unsigned int nc = 0; nc < v_chanFreq_.size(); nc++) {

//     v_N_H2OLinesPtr = new vector<complex<double> > ;
//     v_N_H2OLinesPtr->reserve(numLayer_);
//     v_N_H2OContPtr = new vector<complex<double> > ;
//     v_N_H2OContPtr->reserve(numLayer_);
//     v_N_O2LinesPtr = new vector<complex<double> > ;
//     v_N_O2LinesPtr->reserve(numLayer_);
//     v_N_DryContPtr = new vector<complex<double> > ;
//     v_N_DryContPtr->reserve(numLayer_);
//     v_N_O3LinesPtr = new vector<complex<double> > ;
//     v_N_O3LinesPtr->reserve(numLayer_);
//     v_N_COLinesPtr = new vector<complex<double> > ;
//     v_N_COLinesPtr->reserve(numLayer_);
//     v_N_N2OLinesPtr = new vector<complex<double> > ;
//     v_N_N2OLinesPtr->reserve(numLayer_);
//     v_N_NO2LinesPtr = new vector<complex<double> > ;
//     v_N_NO2LinesPtr->reserve(numLayer_);
//     v_N_SO2LinesPtr = new vector<complex<double> > ;
//     v_N_SO2LinesPtr->reserve(numLayer_);

//     for(unsigned int n = 0; n < numLayer_; n++) {

//       // cout << "numLayer_=" << nc << " " << n << endl; // COMMENTED OUT BY JUAN MAY/16/2005

//       v_N_H2OLinesPtr->push_back(a.vv_N_H2OLinesPtr_[nc]->at(n));
//       v_N_H2OContPtr->push_back(a.vv_N_H2OContPtr_[nc]->at(n));
//       v_N_O2LinesPtr->push_back(a.vv_N_O2LinesPtr_[nc]->at(n));
//       v_N_DryContPtr->push_back(a.vv_N_DryContPtr_[nc]->at(n));
//       v_N_O3LinesPtr->push_back(a.vv_N_O3LinesPtr_[nc]->at(n));
//       v_N_COLinesPtr->push_back(a.vv_N_COLinesPtr_[nc]->at(n));
//       v_N_N2OLinesPtr->push_back(a.vv_N_N2OLinesPtr_[nc]->at(n));
//       v_N_NO2LinesPtr->push_back(a.vv_N_NO2LinesPtr_[nc]->at(n));
//       v_N_SO2LinesPtr->push_back(a.vv_N_SO2LinesPtr_[nc]->at(n));

//     }

//     vv_N_H2OLinesPtr_.push_back(v_N_H2OLinesPtr);
//     vv_N_H2OContPtr_.push_back(v_N_H2OContPtr);
//     vv_N_O2LinesPtr_.push_back(v_N_O2LinesPtr);
//     vv_N_DryContPtr_.push_back(v_N_DryContPtr);
//     vv_N_O3LinesPtr_.push_back(v_N_O3LinesPtr);
//     vv_N_COLinesPtr_.push_back(v_N_COLinesPtr);
//     vv_N_N2OLinesPtr_.push_back(v_N_N2OLinesPtr);
//     vv_N_NO2LinesPtr_.push_back(v_N_NO2LinesPtr);
//     vv_N_SO2LinesPtr_.push_back(v_N_SO2LinesPtr);

//   }

  // level Atm Radiance

  airMass_ = a.airMass_;
  skyBackgroundTemperature_ = a.skyBackgroundTemperature_;
  wh2o_user_ = a.wh2o_user_;

}

SkyStatus::~SkyStatus()
{
  rmSkyStatus();
}

bool SkyStatus::setBasicAtmosphericParameters(const Length &altitude,
                                              const Pressure &groundPressure,
                                              const Temperature &groundTemperature,
                                              double tropoLapseRate,
                                              const Humidity &relativeHumidity,
                                              const Length &wvScaleHeight)
{
  bool update = updateProfilesAndRadiance(altitude,
                                          groundPressure,
                                          groundTemperature,
                                          tropoLapseRate,
                                          relativeHumidity,
                                          wvScaleHeight);
  return update;
}

bool SkyStatus::setBasicAtmosphericParameters(const Length &altitude,
                                              const Length &wvScaleHeight)
{
  bool update = updateProfilesAndRadiance(altitude,
                                          groundPressure_,
                                          groundTemperature_,
                                          tropoLapseRate_,
                                          relativeHumidity_,
                                          wvScaleHeight);
  return update;
}

bool SkyStatus::setBasicAtmosphericParameters(const Length &altitude)
{
  bool update = updateProfilesAndRadiance(altitude,
                                          groundPressure_,
                                          groundTemperature_,
                                          tropoLapseRate_,
                                          relativeHumidity_,
                                          wvScaleHeight_);
  return update;
}

bool SkyStatus::setBasicAtmosphericParameters(const Temperature &groundTemperature)
{
  bool update = updateProfilesAndRadiance(altitude_,
                                          groundPressure_,
                                          groundTemperature,
                                          tropoLapseRate_,
                                          relativeHumidity_,
                                          wvScaleHeight_);
  return update;
}

bool SkyStatus::setBasicAtmosphericParameters(const Pressure &groundPressure)
{

  bool update = updateProfilesAndRadiance(altitude_,
                                          groundPressure,
                                          groundTemperature_,
                                          tropoLapseRate_,
                                          relativeHumidity_,
                                          wvScaleHeight_);

  return update;
}

bool SkyStatus::setBasicAtmosphericParameters(const Humidity &relativeHumidity)
{
  bool update = updateProfilesAndRadiance(altitude_,
                                          groundPressure_,
                                          groundTemperature_,
                                          tropoLapseRate_,
                                          relativeHumidity,
                                          wvScaleHeight_);
  return update;
}

bool SkyStatus::setBasicAtmosphericParameters(double tropoLapseRate)
{
  bool update = updateProfilesAndRadiance(altitude_,
                                          groundPressure_,
                                          groundTemperature_,
                                          tropoLapseRate,
                                          relativeHumidity_,
                                          wvScaleHeight_);
  return update;
}

bool SkyStatus::setBasicAtmosphericParameters(const Length &altitude,
                                              const Temperature &groundTemperature)
{
  bool update = updateProfilesAndRadiance(altitude,
                                          groundPressure_,
                                          groundTemperature,
                                          tropoLapseRate_,
                                          relativeHumidity_,
                                          wvScaleHeight_);
  return update;
}

bool SkyStatus::setBasicAtmosphericParameters(const Length &altitude,
                                              const Pressure &groundPressure)
{
  bool update = updateProfilesAndRadiance(altitude,
                                          groundPressure,
                                          groundTemperature_,
                                          tropoLapseRate_,
                                          relativeHumidity_,
                                          wvScaleHeight_);
  return update;
}

bool SkyStatus::setBasicAtmosphericParameters(const Length &altitude,
                                              const Humidity &relativeHumidity)
{
  bool update = updateProfilesAndRadiance(altitude,
                                          groundPressure_,
                                          groundTemperature_,
                                          tropoLapseRate_,
                                          relativeHumidity,
                                          wvScaleHeight_);
  return update;
}

bool SkyStatus::setBasicAtmosphericParameters(const Length &altitude,
                                              double tropoLapseRate)
{
  bool update = updateProfilesAndRadiance(altitude,
                                          groundPressure_,
                                          groundTemperature_,
                                          tropoLapseRate,
                                          relativeHumidity_,
                                          wvScaleHeight_);
  return update;
}

bool SkyStatus::setBasicAtmosphericParameters(const Temperature &groundTemperature,
                                              const Pressure &groundPressure)
{
  bool update = updateProfilesAndRadiance(altitude_,
                                          groundPressure,
                                          groundTemperature,
                                          tropoLapseRate_,
                                          relativeHumidity_,
                                          wvScaleHeight_);
  return update;
}
bool SkyStatus::setBasicAtmosphericParameters(const Pressure &groundPressure,
                                              const Temperature &groundTemperature)
{
  bool update = updateProfilesAndRadiance(altitude_,
                                          groundPressure,
                                          groundTemperature,
                                          tropoLapseRate_,
                                          relativeHumidity_,
                                          wvScaleHeight_);
  return update;
}

bool SkyStatus::setBasicAtmosphericParameters(const Temperature &groundTemperature,
                                              const Humidity &relativeHumidity)
{
  bool update = updateProfilesAndRadiance(altitude_,
                                          groundPressure_,
                                          groundTemperature,
                                          tropoLapseRate_,
                                          relativeHumidity,
                                          wvScaleHeight_);
  return update;
}
bool SkyStatus::setBasicAtmosphericParameters(const Humidity &relativeHumidity,
                                              const Temperature &groundTemperature)
{
  bool update = updateProfilesAndRadiance(altitude_,
                                          groundPressure_,
                                          groundTemperature,
                                          tropoLapseRate_,
                                          relativeHumidity,
                                          wvScaleHeight_);
  return update;
}

bool SkyStatus::setBasicAtmosphericParameters(const Temperature &groundTemperature,
                                              double tropoLapseRate)
{
  bool update = updateProfilesAndRadiance(altitude_,
                                          groundPressure_,
                                          groundTemperature,
                                          tropoLapseRate,
                                          relativeHumidity_,
                                          wvScaleHeight_);
  return update;
}
bool SkyStatus::setBasicAtmosphericParameters(double tropoLapseRate,
                                              const Temperature &groundTemperature)
{
  bool update = updateProfilesAndRadiance(altitude_,
                                          groundPressure_,
                                          groundTemperature,
                                          tropoLapseRate,
                                          relativeHumidity_,
                                          wvScaleHeight_);
  return update;
}

bool SkyStatus::setBasicAtmosphericParameters(const Temperature &groundTemperature,
                                              const Length &wvScaleHeight)
{
  bool update = updateProfilesAndRadiance(altitude_,
                                          groundPressure_,
                                          groundTemperature,
                                          tropoLapseRate_,
                                          relativeHumidity_,
                                          wvScaleHeight);
  return update;
}

bool SkyStatus::setBasicAtmosphericParameters(const Pressure &groundPressure,
                                              const Humidity &relativeHumidity)
{
  bool update = updateProfilesAndRadiance(altitude_,
                                          groundPressure,
                                          groundTemperature_,
                                          tropoLapseRate_,
                                          relativeHumidity,
                                          wvScaleHeight_);
  return update;
}
bool SkyStatus::setBasicAtmosphericParameters(const Humidity &relativeHumidity,
                                              const Pressure &groundPressure)
{
  bool update = updateProfilesAndRadiance(altitude_,
                                          groundPressure,
                                          groundTemperature_,
                                          tropoLapseRate_,
                                          relativeHumidity,
                                          wvScaleHeight_);
  return update;
}

bool SkyStatus::setBasicAtmosphericParameters(const Pressure &groundPressure,
                                              double tropoLapseRate)
{
  bool update = updateProfilesAndRadiance(altitude_,
                                          groundPressure,
                                          groundTemperature_,
                                          tropoLapseRate,
                                          relativeHumidity_,
                                          wvScaleHeight_);
  return update;
}
bool SkyStatus::setBasicAtmosphericParameters(double tropoLapseRate,
                                              const Pressure &groundPressure)
{
  bool update = updateProfilesAndRadiance(altitude_,
                                          groundPressure,
                                          groundTemperature_,
                                          tropoLapseRate,
                                          relativeHumidity_,
                                          wvScaleHeight_);
  return update;
}

bool SkyStatus::setBasicAtmosphericParameters(const Pressure &groundPressure,
                                              const Length &wvScaleHeight)
{
  bool update = updateProfilesAndRadiance(altitude_,
                                          groundPressure,
                                          groundTemperature_,
                                          tropoLapseRate_,
                                          relativeHumidity_,
                                          wvScaleHeight);
  return update;
}

bool SkyStatus::setBasicAtmosphericParameters(const Humidity &relativeHumidity,
                                              double tropoLapseRate)
{
  bool update = updateProfilesAndRadiance(altitude_,
                                          groundPressure_,
                                          groundTemperature_,
                                          tropoLapseRate,
                                          relativeHumidity,
                                          wvScaleHeight_);
  return update;
}
bool SkyStatus::setBasicAtmosphericParameters(double tropoLapseRate,
                                              const Humidity &relativeHumidity)
{
  bool update = updateProfilesAndRadiance(altitude_,
                                          groundPressure_,
                                          groundTemperature_,
                                          tropoLapseRate,
                                          relativeHumidity,
                                          wvScaleHeight_);
  return update;
}

bool SkyStatus::setBasicAtmosphericParameters(const Humidity &relativeHumidity,
                                              const Length &wvScaleHeight)
{
  bool update = updateProfilesAndRadiance(altitude_,
                                          groundPressure_,
                                          groundTemperature_,
                                          tropoLapseRate_,
                                          relativeHumidity,
                                          wvScaleHeight);
  return update;
}

bool SkyStatus::setBasicAtmosphericParameters(double tropoLapseRate,
                                              const Length &wvScaleHeight)
{
  bool update = updateProfilesAndRadiance(altitude_,
                                          groundPressure_,
                                          groundTemperature_,
                                          tropoLapseRate,
                                          relativeHumidity_,
                                          wvScaleHeight);
  return update;
}

bool SkyStatus::updateProfilesAndRadiance(const Length &altitude,
                                          const Pressure &groundPressure,
                                          const Temperature &groundTemperature,
                                          double tropoLapseRate,
                                          const Humidity &relativeHumidity,
                                          const Length &wvScaleHeight)
{

  bool updated = false;

  bool mkNewAbsPhaseProfile = updateRefractiveIndexProfile(altitude,
                                                           groundPressure,
                                                           groundTemperature,
                                                           tropoLapseRate,
                                                           relativeHumidity,
                                                           wvScaleHeight);

  if(mkNewAbsPhaseProfile) {
    updated = true;
  }
  return updated;
}

Opacity SkyStatus::getH2OLinesOpacity(unsigned int nc)
{
  if(!chanIndexIsValid(nc)) return (double) -999.0;
  double kv = 0;
  for(unsigned int j = 0; j < numLayer_; j++) {
    kv = kv + imag(vv_N_H2OLinesPtr_[nc]->at(j)) * v_layerThickness_[j];
  }
  return ((getUserWH2O().get()) / (getGroundWH2O().get())) * kv;
}

Opacity SkyStatus::getH2OLinesOpacityUpTo(unsigned int nc, Length refalti)
{
  unsigned int ires; unsigned int numlayerold; Length alti;
  Opacity opacityout0; Opacity opacityout1; Opacity zeroOp(0.0,Opacity::UnitNeper);
  double fractionLast; double g1; double g2;

  if(refalti.get(Length::UnitKiloMeter) <= altitude_.get(Length::UnitKiloMeter)) {
    return zeroOp;
  }else{
      fractionLast = 1.0; numlayerold = numLayer_;
      g1=getGroundWH2O().get();
      opacityout0=getH2OLinesOpacity(nc); ires=numlayerold-1; alti=altitude_;
      for(unsigned int i=0; i<numLayer_; i++){
	if(alti.get(Length::UnitKiloMeter) < refalti.get(Length::UnitKiloMeter) &&  (alti.get(Length::UnitKiloMeter)+v_layerThickness_[i]/1000.0) >= refalti.get(Length::UnitKiloMeter))
	  { ires=i; fractionLast = (refalti.get(Length::UnitMeter)-alti.get(Length::UnitMeter))/v_layerThickness_[i]; }
	alti = alti + Length(v_layerThickness_[i],Length::UnitMeter);
      }
      numLayer_ = ires; g2=getGroundWH2O().get();
      opacityout0=getH2OLinesOpacity(nc)*(g2/g1);
      numLayer_ = ires+1; g2=getGroundWH2O().get();
      opacityout1=getH2OLinesOpacity(nc)*(g2/g1);
      numLayer_ = numlayerold;
      return opacityout0+(opacityout1-opacityout0)*fractionLast;
  }
}

/*
Opacity SkyStatus::getTotalOpacityUpTo(unsigned int nc, Length refalti)    //15NOV2017
  {
    unsigned int ires; unsigned int numlayerold; Length alti;
    Opacity opacityout; Opacity opacityout0;
    Opacity opacityout1; Opacity zeroOp(0.0,Opacity::UnitNeper);
    double fractionLast; double g1; double g2;

    if(refalti.get(Length::UnitKiloMeter) <= altitude_.get(Length::UnitKiloMeter)) {
      return zeroOp;
    }else{
      fractionLast = 1.0; numlayerold = numLayer_;
      g1=getGroundWH2O().get(); ires = numlayerold-1; alti=altitude_;
      for(unsigned int i=0; i<numLayer_; i++){
	if(alti.get(Length::UnitKiloMeter) < refalti.get(Length::UnitKiloMeter) &&  (alti.get(Length::UnitKiloMeter)+v_layerThickness_[i]/1000.0) >= refalti.get(Length::UnitKiloMeter))
	  { ires=i; fractionLast = (refalti.get(Length::UnitMeter)-alti.get(Length::UnitMeter))/v_layerThickness_[i]; }
	alti = alti + Length(v_layerThickness_[i],Length::UnitMeter);
      }
      numLayer_ = ires+1; g2=getGroundWH2O().get();
      // std::cout << "En getTotalOpacityUpTo numlayerold numLayer_ fractionLast " << numlayerold << " " << numLayer_ << " " << fractionLast << std::endl;
      // std::cout << "calculando opacityout1" << std::endl;
      opacityout1 = getWetOpacity(nc)*(g2/g1) + getDryOpacity(nc);
      // std::cout << "opacityout1= " << opacityout1.get("np") << std::endl;
      // std::cout << "calculando opacityout0" << std::endl;
      // std::cout << "ires=" << ires << std::endl;
      if(ires == 0){
	opacityout0 = zeroOp;
      }else{
	numLayer_ = ires; g2=getGroundWH2O().get();
	opacityout0 = getWetOpacity(nc)*(g2/g1) + getDryOpacity(nc);
      }
      //    std::cout << "opacityout0= " << opacityout0.get("np") << std::endl;
      opacityout = opacityout0 + (opacityout1-opacityout0)*fractionLast;
      numLayer_ = numlayerold;
      return opacityout;
    }
  }
*/

Opacity SkyStatus::getH2OContOpacity(unsigned int nc)
{
  if(!chanIndexIsValid(nc)) return (double) -999.0;
  double kv = 0;
  for(unsigned int j = 0; j < numLayer_; j++) {
    kv = kv + imag(vv_N_H2OContPtr_[nc]->at(j)) * v_layerThickness_[j];
  }
  return ((getUserWH2O().get()) / (getGroundWH2O().get())) * kv;
}

Opacity SkyStatus::getH2OContOpacityUpTo(unsigned int nc, Length refalti)
{
  unsigned int ires; unsigned int numlayerold; Length alti;
  Opacity opacityout0; Opacity opacityout1; Opacity zeroOp(0.0,Opacity::UnitNeper);
  double fractionLast; double g1; double g2;


  if(refalti.get(Length::UnitKiloMeter) <= altitude_.get(Length::UnitKiloMeter)) {
    return zeroOp;
  }else{
      fractionLast = 1.0; numlayerold = numLayer_;
      g1=getGroundWH2O().get();
      opacityout0=getH2OContOpacity(nc); ires=numlayerold-1; alti=altitude_;
      for(unsigned int i=0; i<numLayer_; i++){
	if(alti.get(Length::UnitKiloMeter) < refalti.get(Length::UnitKiloMeter) &&  (alti.get(Length::UnitKiloMeter)+v_layerThickness_[i]/1000.0) >= refalti.get(Length::UnitKiloMeter))
	  { ires=i; fractionLast = (refalti.get(Length::UnitMeter)-alti.get(Length::UnitMeter))/v_layerThickness_[i]; }
	alti = alti + Length(v_layerThickness_[i],Length::UnitMeter);
      }
      numLayer_ = ires; g2=getGroundWH2O().get();
      opacityout0=getH2OContOpacity(nc)*(g2/g1);
      numLayer_ = ires+1; g2=getGroundWH2O().get();
      opacityout1=getH2OContOpacity(nc)*(g2/g1);
      numLayer_ = numlayerold;
      return opacityout0+(opacityout1-opacityout0)*fractionLast;
  }
}

Angle SkyStatus::getDispersiveH2OPhaseDelay(unsigned int nc)
{
  if(!chanIndexIsValid(nc)) {
    Angle aa(0.0, Angle::UnitDegree);
    return aa;
  }
  double kv = 0;
  for(unsigned int j = 0; j < numLayer_; j++) {
    kv = kv + real(vv_N_H2OLinesPtr_[nc]->at(j)) * v_layerThickness_[j];
  }
  Angle aa(((getUserWH2O().get()) / (getGroundWH2O().get())) * kv * 57.29578,
           Angle::UnitDegree);
  return aa;
}

Length SkyStatus::getDispersiveH2OPathLength(unsigned int nc)
{
  if(!chanIndexIsValid(nc)) {
    Length ll(0.0, Length::UnitMilliMeter);
    return ll;
  }
  double wavelength = 299792458.0 / v_chanFreq_[nc]; // in m
  Length ll((wavelength / 360.0) * getDispersiveH2OPhaseDelay(nc).get(Angle::UnitDegree),
            Length::UnitMeter);
  return ll;
}

Angle SkyStatus::getNonDispersiveH2OPhaseDelay(unsigned int nc)
{
  double kv = 0;
  if(!chanIndexIsValid(nc)) {
    Angle aa(0.0, Angle::UnitDegree);
    return aa;
  }
  for(unsigned int j = 0; j < numLayer_; j++) {
    kv = kv + real(vv_N_H2OContPtr_[nc]->at(j)) * v_layerThickness_[j];
  }
  Angle aa(((getUserWH2O().get()) / (getGroundWH2O().get())) * kv * 57.29578,
           Angle::UnitDegree);
  return aa;
}

Length SkyStatus::getNonDispersiveH2OPathLength(unsigned int nc)
{
  if(!chanIndexIsValid(nc)) {
    Length ll(0.0, Length::UnitMilliMeter);
    return ll;
  }
  double wavelength = 299792458.0 / v_chanFreq_[nc]; // in m
  Length
      ll((wavelength / 360.0) * getNonDispersiveH2OPhaseDelay(nc).get(Angle::UnitDegree),
         Length::UnitMeter);
  return ll;
}

Angle SkyStatus::getAverageDispersiveH2OPhaseDelay(unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) {
    Angle aa(-999.0, Angle::UnitDegree);
    return aa;
  }
  double av = 0.0;
  for(unsigned int i = 0; i < getNumChan(spwid); i++) {
    av = av + getDispersiveH2OPhaseDelay(v_transfertId_[spwid] + i).get(Angle::UnitDegree);
  }
  av = av / getNumChan(spwid);
  Angle average(av, Angle::UnitDegree);
  return average;
}

Length SkyStatus::getAverageDispersiveH2OPathLength(unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) {
    Length ll(0.0, Length::UnitMilliMeter);
    return ll;
  }
  double av = 0.0;
  for(unsigned int i = 0; i < getNumChan(spwid); i++) {
    av = av + getDispersiveH2OPathLength(v_transfertId_[spwid] + i).get(Length::UnitMilliMeter);
  }
  av = av / getNumChan(spwid);
  Length average(av, Length::UnitMilliMeter);
  return average;
}

Angle SkyStatus::getAverageNonDispersiveH2OPhaseDelay(unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) {
    Angle aa(0.0, Angle::UnitDegree);
    return aa;
  }
  double av = 0.0;
  for(unsigned int i = 0; i < getNumChan(spwid); i++) {
    av = av
        + getNonDispersiveH2OPhaseDelay(v_transfertId_[spwid] + i).get(Angle::UnitDegree);
  }
  av = av / getNumChan(spwid);
  Angle average(av, Angle::UnitDegree);
  return average;
}

Length SkyStatus::getAverageNonDispersiveH2OPathLength(unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) {
    Length ll(0.0, Length::UnitMilliMeter);
    return ll;
  }
  double av = 0.0;
  for(unsigned int i = 0; i < getNumChan(spwid); i++) {
    av = av
        + getNonDispersiveH2OPathLength(v_transfertId_[spwid] + i).get(Length::UnitMilliMeter);
  }
  av = av / getNumChan(spwid);
  Length average(av, Length::UnitMilliMeter);
  return average;
}

Temperature SkyStatus::getAverageTebbSky(unsigned int spwid,
                                         const Length &wh2o,
                                         double airmass,
                                         double skycoupling,
                                         const Temperature &Tspill)
{
  Temperature tt(-999, Temperature::UnitKelvin);
  if(!spwidAndIndexAreValid(spwid, 0)) {
    return tt;
  }
  if(wh2o.get() < 0.0) {
    return tt;
  }
  // if(skycoupling<0.0 || skycoupling>1.0){return tt;}
  if(airmass < 1.0) {
    return tt;
  }
  if(Tspill.get(Temperature::UnitKelvin) < 0.0 || Tspill.get(Temperature::UnitKelvin) > 350.0) {
    return tt;
  }
  return RT(((wh2o.get()) / (getGroundWH2O().get())),
            skycoupling,
            Tspill.get(Temperature::UnitKelvin),
            airmass,
            spwid);
}

Temperature SkyStatus::getAverageTebbSky(unsigned int spwid,
                                         const Length &wh2o,
                                         double airmass,
                                         double skycoupling,
                                         double signalgain,     // adition
                                         const Temperature &Tspill)
{
  Temperature tt(-999, Temperature::UnitKelvin);
  if(!spwidAndIndexAreValid(spwid, 0)) {
    return tt;
  }
  if(wh2o.get() < 0.0) {
    return tt;
  }
  // if(skycoupling<0.0 || skycoupling>1.0){return tt;}
  if(airmass < 1.0) {
    return tt;
  }
  if(Tspill.get(Temperature::UnitKelvin) < 0.0 || Tspill.get(Temperature::UnitKelvin) > 350.0) {
    return tt;
  }
  return signalgain*RT(((wh2o.get()) / (getGroundWH2O().get())),
            skycoupling,
            Tspill.get(Temperature::UnitKelvin),
            airmass,
            spwid)+
    +(1-signalgain)*RT(((wh2o.get()) / (getGroundWH2O().get())),
            skycoupling,
            Tspill.get(Temperature::UnitKelvin),
            airmass,
            getAssocSpwId(spwid)[0]);
}

Temperature SkyStatus::getTebbSky(unsigned int spwid,
                                  unsigned int nc,
                                  const Length &wh2o,
                                  double airmass,
                                  double skycoupling,
                                  const Temperature &Tspill)
{
  Temperature tt(-999, Temperature::UnitKelvin);
  if(!spwidAndIndexAreValid(spwid, nc)) {
    return tt;
  }
  if(wh2o.get() < 0.0) {
    return tt;
  }
  if(skycoupling < 0.0 || skycoupling > 1.0) {
    return tt;
  }
  if(airmass < 1.0) {
    return tt;
  }
  if(Tspill.get(Temperature::UnitKelvin) < 0.0 || Tspill.get(Temperature::UnitKelvin) > 350.0) {
    return tt;
  }
  return RT(((wh2o.get()) / (getGroundWH2O().get())),
            skycoupling,
            Tspill.get(Temperature::UnitKelvin),
            airmass,
            spwid,
            nc);
}


Temperature SkyStatus::getAverageTrjSky(unsigned int spwid,
                                         const Length &wh2o,
                                         double airmass,
                                         double skycoupling,
                                         const Temperature &Tspill)
{
  Temperature tt(-999, Temperature::UnitKelvin);
  if(!spwidAndIndexAreValid(spwid, 0)) {
    return tt;
  }
  if(wh2o.get() < 0.0) {
    return tt;
  }
  // if(skycoupling<0.0 || skycoupling>1.0){return tt;}
  if(airmass < 1.0) {
    return tt;
  }
  if(Tspill.get(Temperature::UnitKelvin) < 0.0 || Tspill.get(Temperature::UnitKelvin) > 350.0) {
    return tt;
  }
  return RTRJ(((wh2o.get()) / (getGroundWH2O().get())),
            skycoupling,
            Tspill.get(Temperature::UnitKelvin),
            airmass,
            spwid);
}

Temperature SkyStatus::getAverageTrjSky(unsigned int spwid,
                                         const Length &wh2o,
                                         double airmass,
                                         double skycoupling,
                                         double signalgain,     // adition
                                         const Temperature &Tspill)
{
  Temperature tt(-999, Temperature::UnitKelvin);
  if(!spwidAndIndexAreValid(spwid, 0)) {
    return tt;
  }
  if(wh2o.get() < 0.0) {
    return tt;
  }
  // if(skycoupling<0.0 || skycoupling>1.0){return tt;}
  if(airmass < 1.0) {
    return tt;
  }
  if(Tspill.get(Temperature::UnitKelvin) < 0.0 || Tspill.get(Temperature::UnitKelvin) > 350.0) {
    return tt;
  }
  return signalgain*RTRJ(((wh2o.get()) / (getGroundWH2O().get())),
            skycoupling,
            Tspill.get(Temperature::UnitKelvin),
            airmass,
            spwid)+
    +(1-signalgain)*RTRJ(((wh2o.get()) / (getGroundWH2O().get())),
            skycoupling,
            Tspill.get(Temperature::UnitKelvin),
            airmass,
            getAssocSpwId(spwid)[0]);
}

Temperature SkyStatus::getTrjSky(unsigned int spwid,
                                  unsigned int nc,
                                  const Length &wh2o,
                                  double airmass,
                                  double skycoupling,
                                  const Temperature &Tspill)
{
  Temperature tt(-999, Temperature::UnitKelvin);
  if(!spwidAndIndexAreValid(spwid, nc)) {
    return tt;
  }
  if(wh2o.get() < 0.0) {
    return tt;
  }
  if(skycoupling < 0.0 || skycoupling > 1.0) {
    return tt;
  }
  if(airmass < 1.0) {
    return tt;
  }
  if(Tspill.get(Temperature::UnitKelvin) < 0.0 || Tspill.get(Temperature::UnitKelvin) > 350.0) {
    return tt;
  }
  return RTRJ(((wh2o.get()) / (getGroundWH2O().get())),
            skycoupling,
            Tspill.get(Temperature::UnitKelvin),
            airmass,
            spwid,
            nc);
}



Angle SkyStatus::getDispersiveH2OPhaseDelay(unsigned int spwid, unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) {
    Angle aa(0.0, Angle::UnitDegree);
    return aa;
  }
  return getDispersiveH2OPhaseDelay(v_transfertId_[spwid] + nc);
}

Length SkyStatus::getDispersiveH2OPathLength(unsigned int spwid,
                                             unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) {
    Length ll(0.0, Length::UnitMilliMeter);
    return ll;
  }
  return getDispersiveH2OPathLength(v_transfertId_[spwid] + nc);
}

Angle SkyStatus::getNonDispersiveH2OPhaseDelay(unsigned int spwid,
                                               unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) {
    Angle aa(0.0, Angle::UnitDegree);
    return aa;
  }
  return getNonDispersiveH2OPhaseDelay(v_transfertId_[spwid] + nc);
}

Length SkyStatus::getNonDispersiveH2OPathLength(unsigned int spwid,
                                                unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) return (double) 0.0;
  return getNonDispersiveH2OPathLength(v_transfertId_[spwid] + nc);
}

Length SkyStatus::WaterVaporRetrieval_fromFTS(unsigned int spwId,
                                              const vector<double> &v_transmission,
                                              const Frequency &f1,
                                              const Frequency &f2)
{
  if(f1.get() > f2.get()) {
    return Length(-999, Length::UnitMilliMeter);
  }
  if(v_transmission.size() == getSpectralWindow(spwId).size()) {
    return mkWaterVaporRetrieval_fromFTS(spwId,
                                         v_transmission,
                                         // getAirMass(),    // unused parameter
                                         f1,
                                         f2);
  } else {
    return Length(-999.0, Length::UnitMilliMeter);
  }
}

  Length SkyStatus::WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                                 const vector<Temperature> &v_tebb,
                                                 const vector<double> &skycoupling,
                                                 const vector<Temperature> &tspill)
  {
    vector<vector<double> > spwId_filters;
    vector<double> spwId_filter;
    for(unsigned int i = 0; i < spwId.size(); i++) {
      for(unsigned int n = 0; n < v_numChan_[spwId[i]]; n++) {
        spwId_filter.push_back(1.0);
      }
      spwId_filters.push_back(spwId_filter);
      spwId_filter.clear();
    }
    return WaterVaporRetrieval_fromTEBB(spwId,
                                        v_tebb,
                                        spwId_filters,
                                        skycoupling,
                                        tspill);
  }

  Length SkyStatus::WaterVaporRetrieval_fromTEBB(unsigned int spwId,
                                                 const vector<Temperature> &v_tebb,
                                                 double skycoupling,
                                                 const Temperature &tspill)
  {
    vector<double> spwId_filter;
    for(unsigned int n = 0; n < v_numChan_[spwId]; n++) {
      spwId_filter.push_back(1.0);
    }
    return WaterVaporRetrieval_fromTEBB(spwId,
                                        v_tebb,
                                        spwId_filter,
                                        skycoupling,
                                        tspill);
  }


  Length SkyStatus::WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
						 const vector< vector<Temperature> > &vv_tebb,
						 const vector<double> &skycoupling,
						 const vector<Temperature> &tspill)
  {
    vector<vector<double> > spwId_filters;
    vector<double> spwId_filter;
    for(unsigned int i = 0; i < spwId.size(); i++) {
      for(unsigned int n = 0; n < v_numChan_[spwId[i]]; n++) {
        spwId_filter.push_back(1.0);
      }
      spwId_filters.push_back(spwId_filter);
      spwId_filter.clear();
    }
    return WaterVaporRetrieval_fromTEBB(spwId,
                                        vv_tebb,
                                        spwId_filters,
                                        skycoupling,
                                        tspill);
  }

Length SkyStatus::WaterVaporRetrieval_fromTEBB(unsigned int spwId,
                                               const Percent &signalGain,
                                               const vector<Temperature> &v_tebb,
                                               double skycoupling,
                                               const Temperature &tspill)
{
  vector<double> spwId_filter;
  for(unsigned int n = 0; n < v_numChan_[spwId]; n++) {
    spwId_filter.push_back(1.0);
  }
  return WaterVaporRetrieval_fromTEBB(spwId,
                                      signalGain,
                                      v_tebb,
                                      spwId_filter,
                                      skycoupling,
                                      tspill);
}
Length SkyStatus::WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                               const vector<Percent> &signalGain,
                                               const vector<vector<Temperature> > &vv_tebb,
                                               const vector<double> &skycoupling,
                                               const vector<Temperature> &tspill)
{
  vector<vector<double> > spwId_filters;
  vector<double> spwId_filter;
  for(unsigned int i = 0; i < spwId.size(); i++) {
    for(unsigned int n = 0; n < v_numChan_[spwId[i]]; n++) {
      spwId_filter.push_back(1.0);
    }
    spwId_filters.push_back(spwId_filter);
    spwId_filter.clear();
  }
  return WaterVaporRetrieval_fromTEBB(spwId,
                                      signalGain,
                                      vv_tebb,
                                      spwId_filters,
                                      skycoupling,
                                      tspill);
}
Length SkyStatus::WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                               const vector<Percent> &signalGain,
                                               const vector<Temperature> &v_tebb,
                                               const vector<double> &skycoupling,
                                               const vector<Temperature> &tspill)
{
  vector<vector<double> > spwId_filters;
  vector<double> spwId_filter;
  for(unsigned int i = 0; i < spwId.size(); i++) {
    for(unsigned int n = 0; n < v_numChan_[spwId[i]]; n++) {
      spwId_filter.push_back(1.0);
    }
    spwId_filters.push_back(spwId_filter);
    spwId_filter.clear();
  }
  return WaterVaporRetrieval_fromTEBB(spwId,
                                      signalGain,
                                      v_tebb,
                                      spwId_filters,
                                      skycoupling,
                                      tspill);
}

Length SkyStatus::WaterVaporRetrieval_fromTEBB(unsigned int spwId,
                                               const vector<Temperature> &v_tebb,
                                               const vector<double> &spwId_filter,
                                               double skycoupling,
                                               const Temperature &tspill)
{
  if(v_tebb.size() == getSpectralWindow(spwId).size()) {
    return mkWaterVaporRetrieval_fromTEBB(spwId,
                                          Percent(100.0, Percent::UnitPercent),
                                          v_tebb,
                                          getAirMass(),
                                          spwId_filter,
                                          skycoupling,
                                          tspill);
  } else {
    return Length(-999.0, Length::UnitMilliMeter);
  }
}
Length SkyStatus::WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                               const vector<vector<Temperature> > &vv_tebb,
                                               const vector<vector<double> > &spwId_filters,
                                               const vector<double> &skycoupling,
                                               const vector<Temperature> &tspill)
{
  vector<Percent> signalGain;
  for(unsigned int i = 0; i < spwId.size(); i++) {
    signalGain.push_back(Percent(100.0, Percent::UnitPercent));
  }
  return mkWaterVaporRetrieval_fromTEBB(spwId,
                                        signalGain,
                                        vv_tebb,
                                        getAirMass(),
                                        spwId_filters,
                                        skycoupling,
                                        tspill);
}
Length SkyStatus::WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                               const vector<Temperature> &v_tebb,
                                               const vector<vector<double> > &spwId_filters,
                                               const vector<double> &skycoupling,
                                               const vector<Temperature> &tspill)
{
  vector<Percent> signalGain;
  for(unsigned int i = 0; i < spwId.size(); i++) {
    signalGain.push_back(Percent(100.0, Percent::UnitPercent));
  }
  return mkWaterVaporRetrieval_fromTEBB(spwId,
                                        signalGain,
                                        v_tebb,
                                        getAirMass(),
                                        spwId_filters,
                                        skycoupling,
                                        tspill);
}

Length SkyStatus::WaterVaporRetrieval_fromTEBB(unsigned int spwId,
                                               const Percent &signalGain,
                                               const vector<Temperature> &v_tebb,
                                               const vector<double> &spwId_filter,
                                               double skycoupling,
                                               const Temperature &tspill)
{
  if(v_tebb.size() == getSpectralWindow(spwId).size()) {
    return mkWaterVaporRetrieval_fromTEBB(spwId,
                                          signalGain,
                                          v_tebb,
                                          getAirMass(),
                                          spwId_filter,
                                          skycoupling,
                                          tspill);
  } else {
    return Length(-999.0, Length::UnitMilliMeter);
  }
}
Length SkyStatus::WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                               const vector<Percent> &signalGain,
                                               const vector<vector<Temperature> > &vv_tebb,
                                               const vector<vector<double> > &spwId_filters,
                                               const vector<double> &skycoupling,
                                               const vector<Temperature> &tspill)
{
  return mkWaterVaporRetrieval_fromTEBB(spwId,
                                        signalGain,
                                        vv_tebb,
                                        getAirMass(),
                                        spwId_filters,
                                        skycoupling,
                                        tspill);
}
Length SkyStatus::WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                               const vector<Percent> &signalGain,
                                               const vector<Temperature> &v_tebb,
                                               const vector<vector<double> > &spwId_filters,
                                               const vector<double> &skycoupling,
                                               const vector<Temperature> &tspill)
{
  return mkWaterVaporRetrieval_fromTEBB(spwId,
                                        signalGain,
                                        v_tebb,
                                        getAirMass(),
                                        spwId_filters,
                                        skycoupling,
                                        tspill);
}

Length SkyStatus::WaterVaporRetrieval_fromTEBB(unsigned int spwId,
                                               const vector<Temperature> &v_tebb,
                                               double airmass,
                                               double skycoupling,
                                               const Temperature &tspill)
{
  vector<double> spwId_filter;
  for(unsigned int n = 0; n < v_numChan_[spwId]; n++) {
    spwId_filter.push_back(1.0);
  }
  return WaterVaporRetrieval_fromTEBB(spwId,
                                      v_tebb,
                                      spwId_filter,
                                      airmass,
                                      skycoupling,
                                      tspill);
}
Length SkyStatus::WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                               const vector<vector<Temperature> > &vv_tebb,
                                               double airmass,
                                               const vector<double> &skycoupling,
                                               const vector<Temperature> &tspill)
{
  vector<vector<double> > spwId_filters;
  vector<double> spwId_filter;
  for(unsigned int i = 0; i < spwId.size(); i++) {
    for(unsigned int n = 0; n < v_numChan_[spwId[i]]; n++) {
      spwId_filter.push_back(1.0);
    }
    spwId_filters.push_back(spwId_filter);
    spwId_filter.clear();
  }
  return WaterVaporRetrieval_fromTEBB(spwId,
                                      vv_tebb,
                                      spwId_filters,
                                      airmass,
                                      skycoupling,
                                      tspill);
}
Length SkyStatus::WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                               const vector<Temperature> &v_tebb,
                                               double airmass,
                                               const vector<double> &skycoupling,
                                               const vector<Temperature> &tspill)
{
  vector<vector<double> > spwId_filters;
  vector<double> spwId_filter;
  for(unsigned int i = 0; i < spwId.size(); i++) {
    for(unsigned int n = 0; n < v_numChan_[spwId[i]]; n++) {
      spwId_filter.push_back(1.0);
    }
    spwId_filters.push_back(spwId_filter);
    spwId_filter.clear();
  }
  return WaterVaporRetrieval_fromTEBB(spwId,
                                      v_tebb,
                                      spwId_filters,
                                      airmass,
                                      skycoupling,
                                      tspill);
}

Length SkyStatus::WaterVaporRetrieval_fromTEBB(unsigned int spwId,
                                               const Percent &signalGain,
                                               const vector<Temperature> &v_tebb,
                                               double airmass,
                                               double skycoupling,
                                               const Temperature &tspill)
{
  vector<double> spwId_filter;
  for(unsigned int n = 0; n < v_numChan_[spwId]; n++) {
    spwId_filter.push_back(1.0);
  }
  return WaterVaporRetrieval_fromTEBB(spwId,
                                      signalGain,
                                      v_tebb,
                                      spwId_filter,
                                      airmass,
                                      skycoupling,
                                      tspill);
}
Length SkyStatus::WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                               const vector<Percent> &signalGain,
                                               const vector<vector<Temperature> > &vv_tebb,
                                               double airmass,
                                               const vector<double> &skycoupling,
                                               const vector<Temperature> &tspill)
{
  vector<vector<double> > spwId_filters;
  vector<double> spwId_filter;
  for(unsigned int i = 0; i < spwId.size(); i++) {
    for(unsigned int n = 0; n < v_numChan_[spwId[i]]; n++) {
      spwId_filter.push_back(1.0);
    }
    spwId_filters.push_back(spwId_filter);
    spwId_filter.clear();
  }
  return WaterVaporRetrieval_fromTEBB(spwId,
                                      signalGain,
                                      vv_tebb,
                                      spwId_filters,
                                      airmass,
                                      skycoupling,
                                      tspill);
}

  Length SkyStatus::WaterVaporRetrieval_fromTEBB(unsigned int spwId,
						 const Percent &signalGain,
						 const Temperature &tebb,
						 double airmass,
						 double skycoupling,
						 const Temperature &tspill)
  {
    vector<unsigned int> spwIdv;
    vector<Percent> signalGainv;
    vector<Temperature> tebbv;
    vector<double> skycouplingv;
    vector<Temperature> tspillv;
    spwIdv.push_back(spwId);
    signalGainv.push_back(signalGain);
    tebbv.push_back(tebb);
    skycouplingv.push_back(skycoupling);
    tspillv.push_back(tspill);
    return WaterVaporRetrieval_fromTEBB(spwIdv,signalGainv,tebbv,airmass,skycouplingv,tspillv);
  }



Length SkyStatus::WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                               const vector<Percent> &signalGain,
                                               const vector<Temperature> &v_tebb,
                                               double airmass,
                                               const vector<double> &skycoupling,
                                               const vector<Temperature> &tspill)
{
  vector<vector<double> > spwId_filters;
  vector<double> spwId_filter;
  for(unsigned int i = 0; i < spwId.size(); i++) {
    for(unsigned int n = 0; n < v_numChan_[spwId[i]]; n++) {
      spwId_filter.push_back(1.0);
    }
    spwId_filters.push_back(spwId_filter);
    spwId_filter.clear();
  }
  return WaterVaporRetrieval_fromTEBB(spwId,
                                      signalGain,
                                      v_tebb,
                                      spwId_filters,
                                      airmass,
                                      skycoupling,
                                      tspill);
}

Length SkyStatus::WaterVaporRetrieval_fromTEBB(unsigned int spwId,
                                               const vector<Temperature> &v_tebb,
                                               const vector<double> &spwId_filter,
                                               double airmass,
                                               double skycoupling,
                                               const Temperature &tspill)
{
  if(v_tebb.size() == getSpectralWindow(spwId).size()) {
    return mkWaterVaporRetrieval_fromTEBB(spwId,
                                          Percent(100.0, Percent::UnitPercent),
                                          v_tebb,
                                          airmass,
                                          spwId_filter,
                                          skycoupling,
                                          tspill);
  } else {
    return Length(-999.0, Length::UnitMilliMeter);
  }
}
Length SkyStatus::WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                               const vector<vector<Temperature> > &vv_tebb,
                                               const vector<vector<double> > &spwId_filters,
                                               double airmass,
                                               const vector<double> &skycoupling,
                                               const vector<Temperature> &tspill)
{
  for(unsigned int j = 0; j < spwId.size(); j++) {
    if(vv_tebb[j].size() != getSpectralWindow(spwId[j]).size()) {
      return Length(-999.0, Length::UnitMilliMeter);
    }
  }
  vector<Percent> signalGain;
  for(unsigned int i = 0; i < spwId.size(); i++) {
    signalGain.push_back(Percent(100.0, Percent::UnitPercent));
  }
  return mkWaterVaporRetrieval_fromTEBB(spwId,
                                        signalGain,
                                        vv_tebb,
                                        airmass,
                                        spwId_filters,
                                        skycoupling,
                                        tspill);
}
Length SkyStatus::WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                               const vector<Temperature> &v_tebb,
                                               const vector<vector<double> > &spwId_filters,
                                               double airmass,
                                               const vector<double> &skycoupling,
                                               const vector<Temperature> &tspill)
{
  vector<Percent> signalGain;
  for(unsigned int i = 0; i < spwId.size(); i++) {
    signalGain.push_back(Percent(100.0, Percent::UnitPercent));
  }
  return mkWaterVaporRetrieval_fromTEBB(spwId,
                                        signalGain,
                                        v_tebb,
                                        airmass,
                                        spwId_filters,
                                        skycoupling,
                                        tspill);
}

Length SkyStatus::WaterVaporRetrieval_fromTEBB(unsigned int spwId,
                                               const Percent &signalGain,
                                               const vector<Temperature> &v_tebb,
                                               const vector<double> &spwId_filter,
                                               double airmass,
                                               double skycoupling,
                                               const Temperature &tspill)
{
  if(v_tebb.size() == getSpectralWindow(spwId).size()) {
    return mkWaterVaporRetrieval_fromTEBB(spwId,
                                          signalGain,
                                          v_tebb,
                                          airmass,
                                          spwId_filter,
                                          skycoupling,
                                          tspill);
  } else {
    return Length(-999.0, Length::UnitMilliMeter);
  }
}

Length SkyStatus::WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                               const vector<Percent> &signalGain,
                                               const vector<vector<Temperature> > &vv_tebb,
                                               const vector<vector<double> > &spwId_filter,
                                               double airmass,
                                               const vector<double> &skycoupling,
                                               const vector<Temperature> &tspill)
{

  if(spwId.size() != signalGain.size()) {
    return Length(-999.0, Length::UnitMilliMeter);
  }
  for(unsigned int j = 0; j < spwId.size(); j++) {
    if(vv_tebb[j].size() != getSpectralWindow(spwId[j]).size()) {
      return Length(-999.0, Length::UnitMilliMeter);
    }
  }
  if(spwId.size() != spwId_filter.size()) {
    return Length(-999.0, Length::UnitMilliMeter);
  }
  if(spwId.size() != skycoupling.size()) {
    return Length(-999.0, Length::UnitMilliMeter);
  }
  if(spwId.size() != tspill.size()) {
    return Length(-999.0, Length::UnitMilliMeter);
  }

  return mkWaterVaporRetrieval_fromTEBB(spwId,
                                        signalGain,
                                        vv_tebb,
                                        airmass,
                                        spwId_filter,
                                        skycoupling,
                                        tspill);
}

Length SkyStatus::WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                               const vector<Percent> &signalGain,
                                               const vector<Temperature> &v_tebb,
                                               const vector<vector<double> > &spwId_filter,
                                               double airmass,
                                               const vector<double> &skycoupling,
                                               const vector<Temperature> &tspill)
{
  if(spwId.size() != signalGain.size()) {
    return Length(-999.0, Length::UnitMilliMeter);
  }
  if(spwId.size() != v_tebb.size()) {
    return Length(-999.0, Length::UnitMilliMeter);
  }
  if(spwId.size() != spwId_filter.size()) {
    return Length(-999.0, Length::UnitMilliMeter);
  }
  if(spwId.size() != skycoupling.size()) {
    return Length(-999.0, Length::UnitMilliMeter);
  }
  if(spwId.size() != tspill.size()) {
    return Length(-999.0, Length::UnitMilliMeter);
  }

  return mkWaterVaporRetrieval_fromTEBB(spwId,
                                        signalGain,
                                        v_tebb,
                                        airmass,
                                        spwId_filter,
                                        skycoupling,
                                        tspill);
}

Length SkyStatus::WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                               const Percent &signalGain,
                                               const vector<Temperature> &v_tebb,
                                               const vector<vector<double> > &spwId_filters,
                                               double airmass,
                                               double skycoupling,
                                               const Temperature &tspill)
{
  vector<Percent> v_signalGain;
  vector<double> v_skycoupling;
  vector<Temperature> v_tspill;
  v_signalGain.reserve(spwId.size());
  v_skycoupling.reserve(spwId.size());
  v_tspill.reserve(spwId.size());

  for(unsigned int j = 0; j < spwId.size(); j++) {
    v_signalGain.push_back(signalGain);
    v_skycoupling.push_back(skycoupling);
    v_tspill.push_back(tspill);
  }
  return WaterVaporRetrieval_fromTEBB(spwId,
                                      v_signalGain,
                                      v_tebb,
                                      spwId_filters,
                                      airmass,
                                      v_skycoupling,
                                      v_tspill);
}
Length SkyStatus::WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                               const Percent &signalGain,
                                               const vector<Temperature> &v_tebb,
                                               double airmass,
                                               double skycoupling,
                                               const Temperature &tspill)
{
  vector<vector<double> > spwId_filters;
  vector<double> spwId_filter;
  for(unsigned int i = 0; i < spwId.size(); i++) {
    for(unsigned int n = 0; n < v_numChan_[spwId[i]]; n++) {
      spwId_filter.push_back(1.0);
    }
    spwId_filters.push_back(spwId_filter);
    spwId_filter.clear();
  }
  return WaterVaporRetrieval_fromTEBB(spwId,
                                      signalGain,
                                      v_tebb,
                                      spwId_filters,
                                      airmass,
                                      skycoupling,
                                      tspill);
}

double SkyStatus::SkyCouplingRetrieval_fromTEBB(unsigned int spwId,
                                                const vector<Temperature> &v_tebb,
                                                double skycoupling,
                                                const Temperature &tspill)
{
  vector<double> spwId_filter;
  for(unsigned int n = 0; n < v_numChan_[spwId]; n++) {
    spwId_filter.push_back(1.0);
  }
  return SkyCouplingRetrieval_fromTEBB(spwId,
                                       v_tebb,
                                       spwId_filter,
                                       skycoupling,
                                       tspill);
}

double SkyStatus::SkyCouplingRetrieval_fromTEBB(unsigned int spwId,
                                                const vector<Temperature> &v_tebb,
                                                const vector<double> &spwId_filter,
                                                double skycoupling,
                                                const Temperature &tspill)
{
  if(v_tebb.size() == getSpectralWindow(spwId).size()) {
    return mkSkyCouplingRetrieval_fromTEBB(spwId,
                                           Percent(100, Percent::UnitPercent),
                                           v_tebb,
                                           getAirMass(),
                                           spwId_filter,
                                           skycoupling,
                                           tspill);
  } else {
    return -999.0;
  }
}

double SkyStatus::SkyCouplingRetrieval_fromTEBB(unsigned int spwId,
                                                const vector<Temperature> &v_tebb,
                                                double airmass,
                                                double skycoupling,
                                                const Temperature &tspill)
{
  vector<double> spwId_filter;
  for(unsigned int n = 0; n < v_numChan_[spwId]; n++) {
    spwId_filter.push_back(1.0);
  }
  return SkyCouplingRetrieval_fromTEBB(spwId,
                                       v_tebb,
                                       spwId_filter,
                                       airmass,
                                       skycoupling,
                                       tspill);
}

double SkyStatus::SkyCouplingRetrieval_fromTEBB(unsigned int spwId,
                                                const vector<Temperature> &v_tebb,
                                                const vector<double> &spwId_filter,
                                                double airmass,
                                                double skycoupling,
                                                const Temperature &tspill)
{
  if(v_tebb.size() == getSpectralWindow(spwId).size()) {
    return mkSkyCouplingRetrieval_fromTEBB(spwId,
                                           Percent(100, Percent::UnitPercent),
                                           v_tebb,
                                           airmass,
                                           spwId_filter,
                                           skycoupling,
                                           tspill);
  } else {
    return -999.0;
  }
}

double SkyStatus::getSigmaTransmissionFit(unsigned int spwId,
                                          const vector<double> &v_transmission,
                                          double airm,
                                          const Frequency &f1,
                                          const Frequency &f2)
{
  if(f1.get() > f2.get()) {
    return -999.0;
  }
  if(v_transmission.size() == getSpectralWindow(spwId).size()) {
    double rms = 0.0;
    unsigned int num = 0;

    for(unsigned int i = 0; i < v_transmission.size(); i++) {
      if((getSpectralWindow(spwId)[i] * 1E-09 >= f1.get(Frequency::UnitGigaHertz)
          && getSpectralWindow(spwId)[i] * 1E-09 <= f2.get(Frequency::UnitGigaHertz))) {
        num++;
        rms = rms + pow((v_transmission[i] - exp(-airm
						 * ( (getDryContOpacity(spwId, i).get()+getO2LinesOpacity(spwId, i).get()+0.65*getO3LinesOpacity(spwId, i).get()  )  // getDryOpacity(spwId, i).get()
                                + getWetOpacity(spwId, i).get()))),
                        2);
      }
    }
    rms = sqrt(rms / num);
    return rms;
  } else {
    return -999.0;
  }
}

Temperature SkyStatus::getSigmaFit(unsigned int spwId,
                                   const vector<Temperature> &v_tebbspec,
                                   const Length &wh2o,
                                   double airmass,
                                   double skyCoupling,
                                   const Temperature &Tspill)
{
  Temperature ttt(-999, Temperature::UnitKelvin);
  if(!spwidAndIndexAreValid(spwId, 0)) {
    return ttt;
  }
  if(v_tebbspec.size() != getSpectralWindow(spwId).size()) {
    return ttt;
  }
  if(wh2o.get(Length::UnitMilliMeter) < 0.0) {
    return ttt;
  }
  if(skyCoupling < 0.0 || skyCoupling > 1.0) {
    return ttt;
  }
  if(airmass < 1.0) {
    return ttt;
  }
  if(Tspill.get(Temperature::UnitKelvin) < 0.0 || Tspill.get(Temperature::UnitKelvin) > 350.0) {
    return ttt;
  }
  double rms = 0.0;
  unsigned int num = 0;
  for(unsigned int i = 0; i < v_tebbspec.size(); i++) {
    if(v_tebbspec[i].get() > 0.0) {
      num++;
      rms = rms + pow((v_tebbspec[i].get(Temperature::UnitKelvin) - getTebbSky(spwId,
                                                           i,
                                                           wh2o,
                                                           airmass,
                                                           skyCoupling,
                                                           Tspill).get(Temperature::UnitKelvin)),
                      2);
    }
  }
  rms = sqrt(rms / num);
  return Temperature(rms, Temperature::UnitKelvin);
}

Length SkyStatus::mkWaterVaporRetrieval_fromFTS(unsigned int spwId,
                                                const vector<double> &measuredSkyTransmission,
                                                //double airm,   // unused parameter
                                                const Frequency &fre1,
                                                const Frequency &fre2)
{
  double pfit_wh2o;
  double deltaa = 0.02;
  double sig_fit = -999.0;
  double eps = 0.01;
  vector<double> transmission_fit;
  transmission_fit.reserve(measuredSkyTransmission.size()); // allows resizing

  unsigned int num;
  double flamda;
  unsigned int niter = 20;
  double alpha;
  double beta;
  double array;
  double f1;
  double psave;
  double f2;
  double deriv;
  double chisq1;
  double chisqr;
  double pfit_wh2o_b;
  double res;
  Length wh2o_retrieved(-999.0, Length::UnitMilliMeter);
  Length werr(-888, Length::UnitMilliMeter);
  double sigma_fit_transm0;
  Length sigma_wh2o;

  num = 0;

  flamda = 0.001;
  pfit_wh2o = 1.0; // (getUserWH2O().get(Length::UnitMilliMeter))/(getGroundWH2O().get(Length::UnitMilliMeter));

  unsigned int nl = 0;

  if(fre1.get(Frequency::UnitGigaHertz) < 0) {
    nl = getSpectralWindow(spwId).size();
  } else {
    for(unsigned int i = 0; i < getSpectralWindow(spwId).size(); i++) {
      if(getSpectralWindow(spwId)[i] * 1E-09 >= fre1.get(Frequency::UnitGigaHertz)
          && getSpectralWindow(spwId)[i] * 1E-09 <= fre2.get(Frequency::UnitGigaHertz)) {
        nl = nl + 1;
      }
    }
  }

  for(unsigned int kite = 0; kite < niter; kite++) {

    num = num + 1;

    beta = 0.0;
    alpha = 0.0;

    //    for(unsigned int i=0; i<getSpectralWindow(spwId).size(); i++){

    for(unsigned int i = 0; i < getSpectralWindow(spwId).size(); i++) {

      if(nl == getSpectralWindow(spwId).size() || (getSpectralWindow(spwId)[i]
          * 1E-09 >= fre1.get(Frequency::UnitGigaHertz) && getSpectralWindow(spwId)[i] * 1E-09
          <= fre2.get(Frequency::UnitGigaHertz))) {
	//
        transmission_fit[i] = exp(-( (getDryContOpacity(spwId, i).get()+getO2LinesOpacity(spwId, i).get()+0.65*getO3LinesOpacity(spwId, i).get() )                                              //getDryOpacity(spwId, i).get()
				     + pfit_wh2o * getWetOpacity(spwId, i).get() ) );
	//
        f1 = transmission_fit[i];
        psave = pfit_wh2o;
        pfit_wh2o = pfit_wh2o + deltaa;
        f2 = exp(-( (getDryContOpacity(spwId, i).get()+getO2LinesOpacity(spwId, i).get()+0.65*getO3LinesOpacity(spwId, i).get() )     //getDryOpacity(spwId, i).get()
		   + pfit_wh2o* getWetOpacity(spwId, i).get()));
        deriv = (f2 - f1) / deltaa;
        pfit_wh2o = psave;
        beta = beta + (measuredSkyTransmission[i] - transmission_fit[i])
            * deriv;
        alpha = alpha + deriv * deriv;
      }
    }

    chisq1 = 0;
    for(unsigned int i = 0; i < getSpectralWindow(spwId).size(); i++) {
      if(nl == getSpectralWindow(spwId).size() || (getSpectralWindow(spwId)[i]
          * 1E-09 >= fre1.get(Frequency::UnitGigaHertz) && getSpectralWindow(spwId)[i] * 1E-09
          <= fre2.get(Frequency::UnitGigaHertz))) {
        res = -transmission_fit[i] + measuredSkyTransmission[i];
        chisq1 = chisq1 + res * res;
      }
    }
    if(nl > 1) {
      chisq1 = chisq1 / (nl - 1);
    }

    adjust: array = 1.0 / (1.0 + flamda);
    pfit_wh2o_b = pfit_wh2o;
    pfit_wh2o_b = pfit_wh2o_b + beta * array / alpha;
    if(pfit_wh2o_b < 0.0) pfit_wh2o_b = 0.9 * pfit_wh2o;

    chisqr = 0;
    for(unsigned int i = 0; i < getSpectralWindow(spwId).size(); i++) {
      if(nl == getSpectralWindow(spwId).size() || (getSpectralWindow(spwId)[i]
          * 1E-09 >= fre1.get(Frequency::UnitGigaHertz) && getSpectralWindow(spwId)[i] * 1E-09
          <= fre2.get(Frequency::UnitGigaHertz))) {
        transmission_fit[i] = exp(-( ( getDryContOpacity(spwId, i).get()+getO2LinesOpacity(spwId, i).get()+0.65*getO3LinesOpacity(spwId, i).get() )             //      (getDryOpacity(spwId, i).get())
				     + pfit_wh2o_b * getWetOpacity(spwId, i).get()));
        res = -transmission_fit[i] + measuredSkyTransmission[i];
        chisqr = chisqr + res * res;
      }
    }

    if(nl > 1) {
      chisqr = chisqr / (nl - 1);
    }

    if(fabs(chisq1 - chisqr) > 0.001) {
      if(chisq1 < chisqr) {
        flamda = flamda * 10.0;
        goto adjust;
      }
    }

    flamda = flamda / 10.0;
    sig_fit = sqrt(chisqr);
    pfit_wh2o = pfit_wh2o_b;
    sigma_wh2o = Length(sqrt(array / alpha) * sig_fit * pfit_wh2o
        * (getUserWH2O().get()), Length::UnitMilliMeter);

    if(fabs(sqrt(chisq1) - sqrt(chisqr)) < eps) {

      /*	for(unsigned int i=0; i<getSpectralWindow(spwId).size(); i++){
       cout << getSpectralWindow(spwId)[i]*1E-09 << "  "  << measuredSkyTransmission[i] << "  "  << transmission_fit[i] << endl;
       } */

      sigma_fit_transm0 = sig_fit;

      wh2o_retrieved = Length(pfit_wh2o * getUserWH2O().get(Length::UnitMilliMeter), Length::UnitMilliMeter);

      goto salir;

    }

  }

  wh2o_retrieved = werr; // Extra error code, fit not reached after 20 iterations
  sigma_fit_transm0 = -888.0; // Extra error code, fit not reached after 20 iterations
  sigma_wh2o = werr; // Extra error code, fit not reached after 20 iterations

  salir:

  sigma_transmission_FTSfit_ = sigma_fit_transm0;

  if(wh2o_retrieved.get() > 0.0) {
    wh2o_user_ = wh2o_retrieved;
  }
  return wh2o_retrieved;

}

Length SkyStatus::mkWaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                                 const vector<Percent> &signalGain,
                                                 const vector<Temperature> &measuredAverageSkyTEBB,
                                                 double airm,
                                                 const vector<vector<double> > &spwId_filter,
                                                 const vector<double> &skycoupling,
                                                 const vector<Temperature> &tspill)
{

  double pfit_wh2o;
  double deltaa = 0.02;
  double sig_fit = -999.0;
  double eps = 0.001;
  vector<Temperature> v_tebb_fit;

  unsigned int num;
  double flamda;
  unsigned int niter = 20;
  double alpha;
  double beta;
  double array;
  double f1;
  double psave;
  double f2;
  double deriv;
  double chisq1;
  double chisqr;
  double pfit_wh2o_b;
  double res;
  Length wh2o_retrieved(-999.0, Length::UnitMilliMeter);
  Length werr(-888, Length::UnitMilliMeter);
  double sigma_fit_transm0;
  Length sigma_wh2o;
  v_tebb_fit.reserve(measuredAverageSkyTEBB.size());

  num = 0;
  flamda = 0.001;

  pfit_wh2o = (getUserWH2O().get(Length::UnitMilliMeter)) / (getGroundWH2O().get(Length::UnitMilliMeter));

  for(unsigned int kite = 0; kite < niter; kite++) {

    num = num + 1;

    beta = 0.0;
    alpha = 0.0;

    for(unsigned int j = 0; j < spwId.size(); j++) {


      f1 = RT(pfit_wh2o,
	      skycoupling[j],
	      tspill[j].get(Temperature::UnitKelvin),
	      airm,
	      spwId[j],
	      spwId_filter[j],
	      signalGain[j]);    // if signalGain[j] < 1.0 then there is an image side
                                 // band and it is correctly taken into account in RT

      v_tebb_fit[j] = Temperature(f1, Temperature::UnitKelvin);

      psave = pfit_wh2o;
      pfit_wh2o = pfit_wh2o + deltaa;

      f2 = RT(pfit_wh2o,
              skycoupling[j],
              tspill[j].get(Temperature::UnitKelvin),
              airm,
              spwId[j],
              spwId_filter[j],
              signalGain[j]);

      deriv = (f2 - f1) / deltaa;

      pfit_wh2o = psave;
      beta = beta + ((measuredAverageSkyTEBB[j]).get(Temperature::UnitKelvin) - f1) * deriv;
      alpha = alpha + deriv * deriv;

    }

    chisq1 = 0;

    for(unsigned int j = 0; j < measuredAverageSkyTEBB.size(); j++) {

      res = -v_tebb_fit[j].get(Temperature::UnitKelvin) + (measuredAverageSkyTEBB[j]).get(Temperature::UnitKelvin);
      chisq1 = chisq1 + res * res;

    }

    if(measuredAverageSkyTEBB.size() > 1) {
      chisq1 = chisq1 / (measuredAverageSkyTEBB.size() - 1);
    }

    adjust: array = 1.0 / (1.0 + flamda);
    pfit_wh2o_b = pfit_wh2o;
    pfit_wh2o_b = pfit_wh2o_b + beta * array / alpha;
    if(pfit_wh2o_b < 0.0) pfit_wh2o_b = 0.9 * pfit_wh2o;

    chisqr = 0;

    for(unsigned int j = 0; j < spwId.size(); j++) {

      v_tebb_fit[j] = Temperature(RT(pfit_wh2o_b,
                                     skycoupling[j],
                                     tspill[j].get(Temperature::UnitKelvin),
                                     airm,
                                     spwId[j],
                                     spwId_filter[j],
                                     signalGain[j]), Temperature::UnitKelvin);

      res = -v_tebb_fit[j].get(Temperature::UnitKelvin) + (measuredAverageSkyTEBB[j]).get(Temperature::UnitKelvin);
      chisqr = chisqr + res * res;

    }

    if(spwId.size() > 1) {
      chisqr = chisqr / (spwId.size() - 1);
    }

    if(fabs(chisq1 - chisqr) > 0.001) {
      if(chisq1 < chisqr) {
        flamda = flamda * 10.0;
        goto adjust;
      }
    }

    flamda = flamda / 10.0;
    sig_fit = sqrt(chisqr);
    pfit_wh2o = pfit_wh2o_b;
    sigma_wh2o = Length(sqrt(array / alpha) * sig_fit * pfit_wh2o
        * (getGroundWH2O().get(Length::UnitMilliMeter)), Length::UnitMilliMeter);

    if(fabs(sqrt(chisq1) - sqrt(chisqr)) < eps) {

      sigma_fit_transm0 = sig_fit;

      wh2o_retrieved = Length(pfit_wh2o * getGroundWH2O().get(Length::UnitMilliMeter), Length::UnitMilliMeter);

      goto salir;

    }

  }

  wh2o_retrieved = werr; // Extra error code, fit not reached after 20 iterations
  sigma_fit_transm0 = -888.0; // Extra error code, fit not reached after 20 iterations
  sigma_wh2o = werr; // Extra error code, fit not reached after 20 iterations

  salir:

  sigma_TEBBfit_ = Temperature(sigma_fit_transm0, Temperature::UnitKelvin);
  if(wh2o_retrieved.get() > 0.0) {
    wh2o_user_ = wh2o_retrieved;
  }
  return wh2o_retrieved;

}

Length SkyStatus::mkWaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                                 const vector<Percent> &signalGain,
                                                 const vector<vector<Temperature> > &measuredSkyTEBB,
                                                 double airm,
                                                 const vector<vector<double> > &spwId_filter,
                                                 const vector<double> &skycoupling,
                                                 const vector<Temperature> &tspill)
{

  double pfit_wh2o;
  double deltaa = 0.02;
  double sig_fit = -999.0;
  double eps = 0.01;
  vector<vector<Temperature> > vv_tebb_fit;
  vector<Temperature> v_tebb_fit;

  unsigned int num;
  double flamda;
  unsigned int niter = 20;
  double alpha;
  double beta;
  double array;
  double f1;
  double psave;
  double f2;
  double deriv;
  double chisq1;
  double chisqr;
  double pfit_wh2o_b;
  double res;
  Length wh2o_retrieved(-999.0, Length::UnitMilliMeter);
  Length werr(-888, Length::UnitMilliMeter);
  double sigma_fit_transm0;
  Length sigma_wh2o;
  vector<double> spwIdNorm;
  vector<double> validchannels;
  spwIdNorm.reserve(spwId.size());
  validchannels.reserve(spwId.size());
  vv_tebb_fit.reserve(spwId.size());

  for(unsigned int j = 0; j < spwId.size(); j++) {
    spwIdNorm[j] = 0.0;
    validchannels[j] = 0.0;
    for(unsigned int i = 0; i < getSpectralWindow(spwId[j]).size(); i++) {
      if(spwId_filter[j][i] > 0) {
        spwIdNorm[j] = spwIdNorm[j] + spwId_filter[j][i];
        validchannels[j] = validchannels[j] + 1.0;
      }
    }
  }

  num = 0;
  flamda = 0.001;
  pfit_wh2o = (getUserWH2O().get(Length::UnitMilliMeter)) / (getGroundWH2O().get(Length::UnitMilliMeter));

  for(unsigned int kite = 0; kite < niter; kite++) {

    num = num + 1;

    beta = 0.0;
    alpha = 0.0;

    for(unsigned int j = 0; j < spwId.size(); j++) {

      v_tebb_fit.clear();

      for(unsigned int i = 0; i < getSpectralWindow(spwId[j]).size(); i++) {

        if(spwId_filter[j][i] > 0) {

          if(signalGain[j].get() < 1.0) {
            f1 = RT(pfit_wh2o,
                    skycoupling[j],
                    tspill[j].get(Temperature::UnitKelvin),
                    airm,
                    spwId[j],
                    i) * signalGain[j].get() + RT(pfit_wh2o,
                                                  skycoupling[j],
                                                  tspill[j].get(Temperature::UnitKelvin),
                                                  airm,
                                                  getAssocSpwId(spwId[j])[0],
                                                  i)
                * (1 - signalGain[j].get());
          } else {
            f1 = RT(pfit_wh2o,
                    skycoupling[j],
                    tspill[j].get(Temperature::UnitKelvin),
                    airm,
                    spwId[j],
                    i);
          }

          psave = pfit_wh2o;
          pfit_wh2o = pfit_wh2o + deltaa;

          v_tebb_fit.push_back(Temperature(f1, Temperature::UnitKelvin));

          if(signalGain[j].get() < 1.0) {
            f2 = (RT(pfit_wh2o,
                     skycoupling[j],
                     tspill[j].get(Temperature::UnitKelvin),
                     airm,
                     spwId[j],
                     i) * signalGain[j].get() + RT(pfit_wh2o,
                                                   skycoupling[j],
                                                   tspill[j].get(Temperature::UnitKelvin),
                                                   airm,
                                                   getAssocSpwId(spwId[j])[0],
                                                   i) * (1
                - signalGain[j].get()));
          } else {
            f2 = (RT(pfit_wh2o,
                     skycoupling[j],
                     tspill[j].get(Temperature::UnitKelvin),
                     airm,
                     spwId[j],
                     i) * signalGain[j].get());
          }

          f2 = f2 * spwId_filter[j][i] * validchannels[j] / spwIdNorm[j];
          f1 = f1 * spwId_filter[j][i] * validchannels[j] / spwIdNorm[j];
          deriv = (f2 - f1) / deltaa;

          pfit_wh2o = psave;
          beta = beta + ((measuredSkyTEBB[j][i]).get(Temperature::UnitKelvin) - f1) * deriv; //*spwId_filter[j][i]/spwIdNorm[j];
          alpha = alpha + deriv * deriv;

        } else {

          v_tebb_fit.push_back(Temperature(-999, Temperature::UnitKelvin));

        }

      }

      vv_tebb_fit.push_back(v_tebb_fit);

    }

    chisq1 = 0;

    for(unsigned int j = 0; j < spwId.size(); j++) {

      for(unsigned int i = 0; i < getSpectralWindow(spwId[j]).size(); i++) {
        if(spwId_filter[j][i] > 0) {
          res
              = (-vv_tebb_fit[j][i].get(Temperature::UnitKelvin) + (measuredSkyTEBB[j][i]).get(Temperature::UnitKelvin))
                  * spwId_filter[j][i] * validchannels[j] / spwIdNorm[j];
          chisq1 = chisq1 + res * res;
        }
      }

    }

    if(spwId.size() > 1) {
      chisq1 = chisq1 / (spwId.size() - 1);
    }

    adjust: array = 1.0 / (1.0 + flamda);
    pfit_wh2o_b = pfit_wh2o;
    pfit_wh2o_b = pfit_wh2o_b + beta * array / alpha;
    if(pfit_wh2o_b < 0.0) pfit_wh2o_b = 0.9 * pfit_wh2o;

    chisqr = 0;

    for(unsigned int j = 0; j < spwId.size(); j++) {

      for(unsigned int i = 0; i < getSpectralWindow(spwId[j]).size(); i++) {

        if(spwId_filter[j][i] > 0) {

          if(signalGain[j].get() < 1.0) {
            vv_tebb_fit[j][i] = Temperature(RT(pfit_wh2o_b,
                                               skycoupling[j],
                                               tspill[j].get(Temperature::UnitKelvin),
                                               airm,
                                               spwId[j],
                                               i) * signalGain[j].get()
                + RT(pfit_wh2o_b,
                     skycoupling[j],
                     tspill[j].get(Temperature::UnitKelvin),
                     airm,
                     getAssocSpwId(spwId[j])[0],
                     i) * (1 - signalGain[j].get()), Temperature::UnitKelvin);
          } else {
            vv_tebb_fit[j][i] = Temperature(RT(pfit_wh2o_b,
                                               skycoupling[j],
                                               tspill[j].get(Temperature::UnitKelvin),
                                               airm,
                                               spwId[j],
                                               i) * signalGain[j].get(), Temperature::UnitKelvin);
          }

          res
              = (-vv_tebb_fit[j][i].get(Temperature::UnitKelvin) + (measuredSkyTEBB[j][i]).get(Temperature::UnitKelvin))
                  * spwId_filter[j][i] * validchannels[j] / spwIdNorm[j];

          chisqr = chisqr + res * res;
        }
      }

    }

    if(spwId.size() > 1) {
      chisqr = chisqr / (spwId.size() - 1);
    }

    if(fabs(chisq1 - chisqr) > 0.001) {
      if(chisq1 < chisqr) {
        flamda = flamda * 10.0;
        goto adjust;
      }
    }

    flamda = flamda / 10.0;
    sig_fit = sqrt(chisqr);
    pfit_wh2o = pfit_wh2o_b;
    sigma_wh2o = Length(sqrt(array / alpha) * sig_fit * pfit_wh2o
        * (getGroundWH2O().get(Length::UnitMilliMeter)), Length::UnitMilliMeter);

    if(fabs(sqrt(chisq1) - sqrt(chisqr)) < eps) {

      sigma_fit_transm0 = sig_fit;

      wh2o_retrieved = Length(pfit_wh2o * getGroundWH2O().get(Length::UnitMilliMeter), Length::UnitMilliMeter);

      goto salir;

    }

  }

  wh2o_retrieved = werr; // Extra error code, fit not reached after 20 iterations
  sigma_fit_transm0 = -888.0; // Extra error code, fit not reached after 20 iterations
  sigma_wh2o = werr; // Extra error code, fit not reached after 20 iterations

  salir:

  sigma_TEBBfit_ = Temperature(sigma_fit_transm0, Temperature::UnitKelvin);
  if(wh2o_retrieved.get() > 0.0) {
    wh2o_user_ = wh2o_retrieved;
  }
  return wh2o_retrieved;

}

double SkyStatus::mkSkyCouplingRetrieval_fromTEBB(unsigned int spwId,
                                                  const Percent &signalGain,
                                                  const vector<Temperature> &measuredSkyTEBB,
                                                  double airm,
                                                  const vector<double> &spwId_filter,
                                                  double skycoupling,
                                                  const Temperature &tspill)
{

  // double pfit_wh2o;  // [-Wunused_but_set_variable]
  double pfit_skycoupling;
  double pfit_skycoupling_b;
  double deltaa = 0.02;
  // double sig_fit = -999.0;  // [-Wunused_but_set_variable]
  double eps = 0.01;
  vector<Temperature> tebb_fit;
  tebb_fit.reserve(measuredSkyTEBB.size()); // allows resizing

  unsigned int num;
  double flamda;
  unsigned int niter = 20;
  double alpha;
  double beta;
  double array;
  double f1;
  double psave;
  double f2;
  double deriv;
  double chisq1;
  double chisqr;
  double res;
  Length wh2o_retrieved(-999.0, Length::UnitMilliMeter);
  Length werr(-888, Length::UnitMilliMeter);
  // double sigma_fit_transm0; // [-Wunused_but_set_variable]
  Length sigma_wh2o;

  num = 0;

  flamda = 0.001;
  // pfit_wh2o = (getUserWH2O().get(Length::UnitMilliMeter)) / (getGroundWH2O().get(Length::UnitMilliMeter)); // [-Wunused_but_set_variable]
  pfit_skycoupling = 1.0;

  for(unsigned int kite = 0; kite < niter; kite++) {

    num = num + 1;

    beta = 0.0;
    alpha = 0.0;

    //    for(unsigned int i=0; i<getSpectralWindow(spwId).size(); i++){

    mkWaterVaporRetrieval_fromTEBB(spwId,
                                   signalGain,
                                   measuredSkyTEBB,
                                   airm,
                                   spwId_filter,
                                   pfit_skycoupling * skycoupling,
                                   tspill);
    f1 = sigma_TEBBfit_.get(Temperature::UnitKelvin);
    // cout << "pfit_skycoupling =" << pfit_skycoupling << "  f1=" << f1 << "  wh2o_user_=" << wh2o_user_.get(Length::UnitMilliMeter) << " mm" << endl;
    psave = pfit_skycoupling;
    pfit_skycoupling = pfit_skycoupling + deltaa;
    mkWaterVaporRetrieval_fromTEBB(spwId,
                                   signalGain,
                                   measuredSkyTEBB,
                                   airm,
                                   spwId_filter,
                                   pfit_skycoupling * skycoupling,
                                   tspill);
    f2 = sigma_TEBBfit_.get(Temperature::UnitKelvin);
    // cout << "pfit_skycoupling =" << pfit_skycoupling << "  f2=" << f2 << "  wh2o_user_=" << wh2o_user_.get(Length::UnitMilliMeter) << " mm" << endl;
    deriv = (f2 - f1) / deltaa;
    pfit_skycoupling = psave;
    beta = beta - f1 * deriv;
    alpha = alpha + deriv * deriv;

    chisq1 = 0;
    res = f1;
    chisq1 = chisq1 + res * res;

    adjust: array = 1.0 / (1.0 + flamda);
    pfit_skycoupling_b = pfit_skycoupling;
    pfit_skycoupling_b = pfit_skycoupling_b + beta * array / alpha;
    if(pfit_skycoupling_b < 0.0) pfit_skycoupling_b = 0.9 * pfit_skycoupling;
    chisqr = 0;
    mkWaterVaporRetrieval_fromTEBB(spwId,
                                   signalGain,
                                   measuredSkyTEBB,
                                   airm,
                                   spwId_filter,
                                   pfit_skycoupling_b * skycoupling,
                                   tspill);
    res = sigma_TEBBfit_.get(Temperature::UnitKelvin);
    // cout << "pfit_skycoupling =" << pfit_skycoupling_b << "  res=" << res << "  wh2o_user_=" << wh2o_user_.get(Length::UnitMilliMeter) << " mm" << endl;
    chisqr = chisqr + res * res;

    if(fabs(chisq1 - chisqr) > 0.001) {
      if(chisq1 < chisqr) {
        flamda = flamda * 10.0;
        goto adjust;
      }
    }

    flamda = flamda / 10.0;
    // sig_fit = sqrt(chisqr); // [-Wunused_but_set_variable]
    pfit_skycoupling = pfit_skycoupling_b;

    if(fabs(sqrt(chisq1) - sqrt(chisqr)) < eps) {
      goto salir;
    }

  }

  wh2o_retrieved = werr; // Extra error code, fit not reached after 20 iterations
  //sigma_fit_transm0 = -888.0; // Extra error code, fit not reached after 20 iterations [-Wunused_but_set_variable]
  sigma_wh2o = werr; // Extra error code, fit not reached after 20 iterations

  salir: return pfit_skycoupling * skycoupling;

}


  double SkyStatus::RT(double pfit_wh2o,
		       double skycoupling,
		       double tspill,
		       double airmass,
		       unsigned int spwid,
		       const vector<double> &spwId_filter,
		       const Percent &signalgain)
  {

    double tebb_channel = 0.0;
    double rtr;
    double norm = 0.0;

    for(unsigned int n = 0; n < v_numChan_[spwid]; n++) {
      if(spwId_filter[n] > 0) {
        norm = norm + spwId_filter[n];
      }
    }

    if(norm == 0.0) {
      return norm;
    }

    for(unsigned int n = 0; n < v_numChan_[spwid]; n++) {

      if(spwId_filter[n] > 0) {

        if(signalgain.get() < 1.0) {
          rtr = RT(pfit_wh2o, skycoupling, tspill, airmass, spwid, n)
              * signalgain.get() + RT(pfit_wh2o,
                                      skycoupling,
                                      tspill,
                                      airmass,
                                      getAssocSpwId(spwid)[0],
                                      n) * (1.0 - signalgain.get());
        } else {
          rtr = RT(pfit_wh2o, skycoupling, tspill, airmass, spwid, n);
        }
        tebb_channel = tebb_channel + rtr * spwId_filter[n] / norm;
      }
    }

    return tebb_channel;
  }



double SkyStatus::RT(double pfit_wh2o,
                     double skycoupling,
                     double tspill,
                     double airm,
                     unsigned int spwid,
                     unsigned int nc)
{

  double radiance;
  double singlefreq;
  // double chanwidth;  // [-Wunused_but_set_variable]
  double tebb;
  double h_div_k = 0.04799274551; /* plank=6.6262e-34,boltz=1.3806E-23 */
  double kv;
  double tau_layer;
  double tbgr = skyBackgroundTemperature_.get(Temperature::UnitKelvin);
  double ratioWater = pfit_wh2o;

  singlefreq = getChanFreq(spwid, nc).get(Frequency::UnitGigaHertz);
  // chanwidth = getChanWidth(spwid, nc).get(Frequency::UnitGigaHertz); // [-Wunused_but_set_variable]
  // cout << "Chan freq. =" << singlefreq << " GHz" << endl;
  // cout << "Chan width =" << chanwidth << " GHz" << endl;

  tebb=0.0;
  kv = 0.0;
  radiance = 0.0;

  for(unsigned int i = 0; i < getNumLayer(); i++) {

    tau_layer = ((getAbsTotalWet(spwid, nc, i).get()) * ratioWater+ getAbsTotalDry(spwid, nc, i).get()) * getLayerThickness(i).get();

    // cout << i << "  " <<  getAbsTotalWet(spwid, nc, i).get() << "  " << getAbsTotalDry(spwid, nc, i).get() << endl;

    radiance = radiance + (1.0 / (exp(h_div_k * singlefreq/ getLayerTemperature(i).get()) - 1.0)) * exp(-kv * airm) * (1.0- exp(-airm * tau_layer));

    kv = kv + tau_layer;

  }

  radiance = skycoupling * (radiance + (1.0 / (exp(h_div_k * singlefreq / tbgr)- 1.0)) * exp(-kv * airm)) + (1.0 / (exp(h_div_k * singlefreq / tspill)- 1.0)) * (1 - skycoupling);

  tebb = h_div_k * singlefreq / log(1 + (1 / radiance));
  // tebb = tebb+  h_div_k * singlefreq * radiance;
  // cout << "singlefreq = " << singlefreq <<  " total opacity = " << kv << " tebb = " << tebb << endl;

  return tebb;

}


double SkyStatus::RTRJ(double pfit_wh2o,
                     double skycoupling,
                     double tspill,
                     double airm,
                     unsigned int spwid,
                     unsigned int nc)
{

  double radiance;
  double singlefreq;
  // double chanwidth;  // [-Wunused_but_set_variable]
  double trj;
  double h_div_k = 0.04799274551; /* plank=6.6262e-34,boltz=1.3806E-23 */
  double kv;
  double tau_layer;
  double tbgr = skyBackgroundTemperature_.get(Temperature::UnitKelvin);
  double ratioWater = pfit_wh2o;

  singlefreq = getChanFreq(spwid, nc).get(Frequency::UnitGigaHertz);
  // chanwidth = getChanWidth(spwid, nc).get(Frequency::UnitGigaHertz); // [-Wunused_but_set_variable]
  trj=0.0;
  kv = 0.0;
  radiance = 0.0;

  for(unsigned int i = 0; i < getNumLayer(); i++) {

    tau_layer = ((getAbsTotalWet(spwid, nc, i).get()) * ratioWater+ getAbsTotalDry(spwid, nc, i).get()) * getLayerThickness(i).get();

    radiance = radiance + (1.0 / (exp(h_div_k * singlefreq/ getLayerTemperature(i).get()) - 1.0)) * exp(-kv * airm) * (1.0- exp(-airm * tau_layer));

    kv = kv + tau_layer;

  }

  radiance = skycoupling * (radiance + (1.0 / (exp(h_div_k * singlefreq / tbgr)- 1.0)) * exp(-kv * airm)) + (1.0 / (exp(h_div_k * singlefreq / tspill)- 1.0)) * (1 - skycoupling);

  trj = h_div_k * singlefreq * radiance;

  return trj;

}



void SkyStatus::iniSkyStatus()
{

  Length wh2o_default(1, Length::UnitMilliMeter);
  Length wh2o_default_neg(-999, Length::UnitMilliMeter);
  Temperature temp_default_neg(-999, Temperature::UnitKelvin);

  if(wh2o_user_.get() <= 0.0 || wh2o_user_.get() > (getGroundWH2O().get())
      * (200 / (getRelativeHumidity().get(Percent::UnitPercent)))) {
    wh2o_user_ = wh2o_default;
  }

}

double SkyStatus::getAverageNonDispersiveDryPathLength_GroundPressureDerivative(unsigned int spwid)
{
  Pressure ref = getGroundPressure();
  Length a =
      RefractiveIndexProfile::getAverageNonDispersiveDryPathLength(spwid);
  setBasicAtmosphericParameters(ref + Pressure(1.0, Pressure::UnitMilliBar));
  Length b =
      RefractiveIndexProfile::getAverageNonDispersiveDryPathLength(spwid);
  setBasicAtmosphericParameters(ref);
  return b.get(Length::UnitMicrons) - a.get(Length::UnitMicrons);
}

double SkyStatus::getAverageNonDispersiveDryPathLength_GroundTemperatureDerivative(unsigned int spwid)
{
  Temperature ref = getGroundTemperature();
  double oldLapseRate = getTropoLapseRate();
  Length a =
      RefractiveIndexProfile::getAverageNonDispersiveDryPathLength(spwid);
  double newLapseRate =
      ((getTropopauseTemperature() - ref - Temperature(1.0, Temperature::UnitKelvin)).get(Temperature::UnitKelvin))
          / ((getTropopauseAltitude() - getAltitude()).get(Length::UnitKiloMeter));
  setBasicAtmosphericParameters(ref + Temperature(1.0, Temperature::UnitKelvin), newLapseRate);
  Length b =
      RefractiveIndexProfile::getAverageNonDispersiveDryPathLength(spwid);
  setBasicAtmosphericParameters(ref, oldLapseRate);
  return b.get(Length::UnitMicrons) - a.get(Length::UnitMicrons);
}

double SkyStatus::getAverageDispersiveDryPathLength_GroundPressureDerivative(unsigned int spwid)
{
  Pressure ref = getGroundPressure();
  Length a = RefractiveIndexProfile::getAverageDispersiveDryPathLength(spwid);
  // scanf("%d",&e);
  setBasicAtmosphericParameters(ref + Pressure(1.0, Pressure::UnitMilliBar));
  Length b = RefractiveIndexProfile::getAverageDispersiveDryPathLength(spwid);
  setBasicAtmosphericParameters(ref);
  return b.get(Length::UnitMicrons) - a.get(Length::UnitMicrons);
}

double SkyStatus::getAverageDispersiveDryPathLength_GroundTemperatureDerivative(unsigned int spwid)
{
  Temperature ref = getGroundTemperature();
  double oldLapseRate = getTropoLapseRate();
  Length a = RefractiveIndexProfile::getAverageDispersiveDryPathLength(spwid);
  double newLapseRate =
      ((getTropopauseTemperature() - ref - Temperature(1.0, Temperature::UnitKelvin)).get(Temperature::UnitKelvin))
          / ((getTropopauseAltitude() - getAltitude()).get(Length::UnitKiloMeter));
  setBasicAtmosphericParameters(ref + Temperature(1.0, Temperature::UnitKelvin), newLapseRate);
  Length b = RefractiveIndexProfile::getAverageDispersiveDryPathLength(spwid);
  setBasicAtmosphericParameters(ref, oldLapseRate);
  return b.get(Length::UnitMicrons) - a.get(Length::UnitMicrons);
}

Temperature SkyStatus::getWVRAverageSigmaTskyFit(const vector<WVRMeasurement> &RadiometerData,
                                                 unsigned int n,
                                                 unsigned int m)
{
  double sigma = 0.0;
  double tr;
  Temperature sigmaT;
  if(m < n) {
    return Temperature(-999, Temperature::UnitKelvin);
  }
  for(unsigned int i = n; i < m; i++) {
    tr = RadiometerData[i].getSigmaFit().get(Temperature::UnitKelvin);
    if(tr < 0) {
      return Temperature(-999, Temperature::UnitKelvin);
    }
    sigma = sigma + tr * tr;
  }
  if(m == n) {
    sigmaT = RadiometerData[n].getSigmaFit();
  } else {
    sigma = sqrt(sigma / (m - n));
    sigmaT = Temperature(sigma, Temperature::UnitKelvin);
  }
  return sigmaT;
}
Temperature SkyStatus::getWVRSigmaChannelTskyFit(const vector<WVRMeasurement> &RadiometerData,
						 unsigned int ichan,
                                                 unsigned int n,
                                                 unsigned int m)
{
  double sigma = 0.0;
  double dtr;
  Temperature sigmaT;
  if(m <= n) {
    return Temperature(-999, Temperature::UnitKelvin);
  }
  for(unsigned int i = n; i < m; i++) {
    dtr = RadiometerData[i].getmeasuredSkyBrightness()[ichan].get(Temperature::UnitKelvin)
      -RadiometerData[i].getfittedSkyBrightness()[ichan].get(Temperature::UnitKelvin);
    sigma = sigma + dtr * dtr;
  }
  sigma = sqrt(sigma / (m - n));
  sigmaT = Temperature(sigma, Temperature::UnitKelvin);
  return sigmaT;
}

WVRMeasurement SkyStatus::mkWaterVaporRetrieval_fromWVR(const vector<Temperature> &measuredSkyBrightnessVector,
                                                        const vector<unsigned int> &IdChannels,
                                                        const vector<double> &skyCoupling,
                                                        const vector<Percent> &signalGain,
                                                        const Temperature &spilloverTemperature,
                                                        const Angle &elevation)
{
  double tspill = spilloverTemperature.get(Temperature::UnitKelvin);
  double pfit_wh2o;
  double deltaa = 0.02;
  double sig_fit = -999.0;
  double eps = 0.01;
  vector<double> tebb_fit;
  tebb_fit.reserve(measuredSkyBrightnessVector.size());
  double airm = 1.0 / sin((3.1415926 * elevation.get(Angle::UnitDegree)) / 180.0);
  unsigned int num;
  double flamda;
  unsigned int niter = 20;
  double alpha;
  double beta;
  double array;
  double f1;
  double psave;
  double f2;
  double deriv;
  double chisq1;
  double chisqr;
  double pfit_wh2o_b;
  double res;
  Length wh2o_retrieved(-999.0, Length::UnitMilliMeter);
  Length werr(-888, Length::UnitMilliMeter);
  Temperature sigma_fit_temp0;
  Temperature t_astro;
  Length sigma_wh2o;

  num = 0;

  flamda = 0.001;

  pfit_wh2o = (getUserWH2O().get(Length::UnitMilliMeter)) / (getGroundWH2O().get(Length::UnitMilliMeter));

  //    cout << "pfit_wh2o=" << pfit_wh2o << endl;

  for(unsigned int kite = 0; kite < niter; kite++) {

    num = num + 1;

    beta = 0.0;
    alpha = 0.0;

    for(unsigned int i = 0; i < IdChannels.size(); i++) {

      tebb_fit[i] = RT(pfit_wh2o,
                       skyCoupling[i],
                       tspill,
                       airm,
                       IdChannels[i],
                       signalGain[i]);
      // cout << i << " " << tebb_fit[i] << endl;


      f1 = tebb_fit[i];
      psave = pfit_wh2o;
      pfit_wh2o = pfit_wh2o + deltaa;
      f2 = RT(pfit_wh2o,
              skyCoupling[i],
              tspill,
              airm,
              IdChannels[i],
              signalGain[i]);
      deriv = (f2 - f1) / deltaa;
      pfit_wh2o = psave;
      beta = beta + (measuredSkyBrightnessVector[i].get(Temperature::UnitKelvin) - tebb_fit[i])
          * deriv;

      alpha = alpha + deriv * deriv;

    }

    chisq1 = 0;
    for(unsigned int i = 0; i < measuredSkyBrightnessVector.size(); i++) {
      res = -tebb_fit[i] + measuredSkyBrightnessVector[i].get(Temperature::UnitKelvin);
      chisq1 = chisq1 + res * res;
    }
    if(measuredSkyBrightnessVector.size() > 1) {
      chisq1 = chisq1 / (measuredSkyBrightnessVector.size() - 1);
    }

    adjust: array = 1.0 / (1.0 + flamda);
    pfit_wh2o_b = pfit_wh2o;
    pfit_wh2o_b = pfit_wh2o_b + beta * array / alpha;
    if(pfit_wh2o_b < 0.0) pfit_wh2o_b = 0.9 * pfit_wh2o;

    for(unsigned int i = 0; i < IdChannels.size(); i++) {

      tebb_fit[i] = RT(pfit_wh2o_b,
                       skyCoupling[i],
                       tspill,
                       airm,
                       IdChannels[i],
                       signalGain[i]);

    }

    chisqr = 0;
    for(unsigned int i = 0; i < IdChannels.size(); i++) {
      res = -tebb_fit[i] + measuredSkyBrightnessVector[i].get(Temperature::UnitKelvin);
      chisqr = chisqr + res * res;
    }
    if(IdChannels.size() > 1) {
      chisqr = chisqr / (IdChannels.size() - 1);
    }

    if(fabs(chisq1 - chisqr) > 0.001) {
      if(chisq1 < chisqr) {
        flamda = flamda * 10.0;
        goto adjust;
      }
    }

    flamda = flamda / 10.0;
    sig_fit = sqrt(chisqr);
    pfit_wh2o = pfit_wh2o_b;
    sigma_wh2o = Length(sqrt(array / alpha) * sig_fit * pfit_wh2o
        * (getGroundWH2O().get(Length::UnitMilliMeter)), Length::UnitMilliMeter);

    if(fabs(sqrt(chisq1) - sqrt(chisqr)) < eps) {

      sigma_fit_temp0 = Temperature(sig_fit, Temperature::UnitKelvin);

      wh2o_retrieved = Length(pfit_wh2o * getGroundWH2O().get(Length::UnitMilliMeter), Length::UnitMilliMeter);

      goto salir;

    }

  }

  wh2o_retrieved = werr; // Extra error code, fit not reached after 20 iterations
  sigma_fit_temp0 = Temperature(sig_fit, Temperature::UnitKelvin); // Extra error code, fit not reached after 20 iterations
  sigma_wh2o = werr; // Extra error code, fit not reached after 20 iterations

  salir:

  vector<Temperature> ttt;

  for(unsigned int i = 0; i < IdChannels.size(); i++) {
    ttt.push_back(Temperature(tebb_fit[i], Temperature::UnitKelvin));
  }

  if(wh2o_retrieved.get() > 0.0) {
    wh2o_user_ = wh2o_retrieved;
  }
  return WVRMeasurement(elevation,
                        measuredSkyBrightnessVector,
                        ttt,
                        wh2o_retrieved,
                        sigma_fit_temp0);

}

void SkyStatus::WaterVaporRetrieval_fromWVR(vector<WVRMeasurement> &RadiometerData,
                                            unsigned int n,
                                            unsigned int m)
{

  for(unsigned int i = n; i < m; i++) {

    WaterVaporRetrieval_fromWVR(RadiometerData[i]);

  }

}

void SkyStatus::WaterVaporRetrieval_fromWVR(WVRMeasurement &RadiometerData)
{

  WVRMeasurement RadiometerData_withRetrieval;

  //    cout << waterVaporRadiometer_.getIdChannels().size() << endl;
  //    cout << RadiometerData.getmeasuredSkyBrightness()[0].get(Temperature::UnitKelvin) << " K" << endl;
  //    cout << waterVaporRadiometer_.getIdChannels()[1] << endl;
  //    cout << getAssocSpwIds(waterVaporRadiometer_.getIdChannels())[1] << endl;


  //    cout << "zz=" << waterVaporRadiometer_.getIdChannels().size() << endl;

  RadiometerData_withRetrieval
      = mkWaterVaporRetrieval_fromWVR(RadiometerData.getmeasuredSkyBrightness(),
                                      waterVaporRadiometer_.getIdChannels(),
                                      waterVaporRadiometer_.getSkyCoupling(),
                                      waterVaporRadiometer_.getsignalGain(),
                                      waterVaporRadiometer_.getSpilloverTemperature(),
                                      RadiometerData.getElevation());

  // cout << "_fromWVR Sky Coupling = " <<  waterVaporRadiometer_.getSkyCoupling()[0] << endl;
  // cout << "Signal Gain = " << waterVaporRadiometer_.getsignalGain()[0].get(Percent::UnitPercent) << " %" <<  endl;
  // cout << "Spillover Temp. = " << waterVaporRadiometer_.getSpilloverTemperature().get(Temperature::UnitKelvin) << " K" << endl;
  // cout << "Elevation = " << RadiometerData.getElevation().get(Angle::UnitDegree) << endl;
  // cout << "PWV=" << RadiometerData_withRetrieval.getretrievedWaterVaporColumn().get(Length::UnitMilliMeter) << " mm" << endl;

  RadiometerData.setretrievedWaterVaporColumn(RadiometerData_withRetrieval.getretrievedWaterVaporColumn());
  RadiometerData.setfittedSkyBrightness(RadiometerData_withRetrieval.getfittedSkyBrightness());
  RadiometerData.setSigmaFit(RadiometerData_withRetrieval.getSigmaFit());

}

void SkyStatus::updateSkyCoupling_fromWVR(vector<WVRMeasurement> &RadiometerData,
                                          unsigned int n,
                                          unsigned int m)
{
  double pfit;
  double deltaa = 0.02;
  // double sig_fit = -999.0; // [-Wunused_but_set_variable]
  double eps = 0.01;

  unsigned int num;
  double flamda;
  unsigned int niter = 20;
  double alpha;
  double beta;
  double array;
  double f1;
  double psave;
  double f2;
  double deriv;
  double chisq1;
  double chisqr;
  double pfit_b;
  double res;

  num = 0;

  flamda = 0.001;
  pfit = 0.5;

  // Find the maximum of skycoupling on the WVR channels
  // This value is used to assure that the found skycoupling does not exceed 1
  double maxCoupling =0.;
  for (unsigned int i=0; i<waterVaporRadiometer_.getSkyCoupling().size(); i++)
    if (waterVaporRadiometer_.getSkyCoupling()[i]>maxCoupling)
      maxCoupling = waterVaporRadiometer_.getSkyCoupling()[i];


  for(unsigned int kite = 0; kite < niter; kite++) {

    num = num + 1;
    beta = 0.0;
    alpha = 0.0;

    if (pfit*maxCoupling>1.5)
      pfit = 1.-deltaa;
    f1 = sigmaSkyCouplingRetrieval_fromWVR(pfit,
                                           waterVaporRadiometer_,
                                           RadiometerData,
                                           n,
                                           m);
    psave = pfit;
    pfit = pfit + deltaa;
    f2 = sigmaSkyCouplingRetrieval_fromWVR(pfit,
                                           waterVaporRadiometer_,
                                           RadiometerData,
                                           n,
                                           m);
    deriv = (f2 - f1) / deltaa;
    pfit = psave;
    beta = beta - f1 * deriv;
    alpha = alpha + deriv * deriv;
    chisq1 = f1 * f1;

    adjust: array = 1.0 / (1.0 + flamda);
    pfit_b = pfit;
    pfit_b = pfit_b + beta * array / alpha;

    chisqr = 0.;

    if(pfit_b < 0.0) {
      pfit_b = 0.9 * pfit;
    }
    if (pfit_b*maxCoupling>1.5)
      pfit_b = 1;

    res = sigmaSkyCouplingRetrieval_fromWVR(pfit_b,
                                            waterVaporRadiometer_,
                                            RadiometerData,
                                            n,
                                            m);
    chisqr = chisqr + res * res;

    if(fabs(chisq1 - chisqr) > 0.001) {
      if(chisq1 < chisqr) {
        flamda = flamda * 10.0;
        goto adjust;
      }
    }

    flamda = flamda / 10.0;
    // sig_fit = sqrt(chisqr); // [-Wunused_but_set_variable]
    pfit = pfit_b;

    if(fabs(sqrt(chisq1) - sqrt(chisqr)) < eps) {
      goto salir;
    }

  }

  salir: waterVaporRadiometer_.multiplySkyCoupling(pfit);

  //    cout << "pfit=" << pfit << "  sky_coupling: " << waterVaporRadiometer_.getSkyCoupling()[0]   << endl;


}
void SkyStatus::updateSkyCouplingChannel_fromWVR(vector<WVRMeasurement> &RadiometerData,
						 unsigned int ichan,
						 unsigned int n,
						 unsigned int m)
{
  double pfit;
  double deltaa = 0.02;
  // double sig_fit = -999.0;    // [-Wunused_but_set_variable]
  double eps = 0.01;

  unsigned int num;
  double flamda;
  unsigned int niter = 20;
  double alpha;
  double beta;
  double array;
  double f1;
  double psave;
  double f2;
  double deriv;
  double chisq1;
  double chisqr;
  double pfit_b;
  double res;

  num = 0;

  flamda = 0.001;
  pfit = 1.00;

  // This value is used to assure that the found skycoupling does not exceed 1
  double maxCoupling = waterVaporRadiometer_.getSkyCoupling()[ichan];


  for(unsigned int kite = 0; kite < niter; kite++) {

    num = num + 1;
    beta = 0.0;
    alpha = 0.0;

    if (pfit*maxCoupling>1)
      pfit = 1.-deltaa;
    f1 = sigmaSkyCouplingChannelRetrieval_fromWVR(pfit,
						  waterVaporRadiometer_,
						  RadiometerData,
						  ichan,
						  n,
						  m);
    psave = pfit;
    pfit = pfit + deltaa;
    f2 = sigmaSkyCouplingChannelRetrieval_fromWVR(pfit,
						  waterVaporRadiometer_,
						  RadiometerData,
						  ichan,
						  n,
						  m);
    deriv = (f2 - f1) / deltaa;
    pfit = psave;
    beta = beta - f1 * deriv;
    alpha = alpha + deriv * deriv;
    chisq1 = f1 * f1;

    adjust: array = 1.0 / (1.0 + flamda);
    pfit_b = pfit;
    pfit_b = pfit_b + beta * array / alpha;

    chisqr = 0.;

    if(pfit_b < 0.0) {
      pfit_b = 0.9 * pfit;
    }
    if (pfit_b*maxCoupling>1)
      pfit_b = 1/maxCoupling;

    res = sigmaSkyCouplingChannelRetrieval_fromWVR(pfit_b,
						   waterVaporRadiometer_,
						   RadiometerData,
						   ichan,
						   n,
						   m);
    chisqr = chisqr + res * res;

    if(fabs(chisq1 - chisqr) > 0.001) {
      if(chisq1 < chisqr) {
        flamda = flamda * 10.0;
        goto adjust;
      }
    }

    flamda = flamda / 10.0;
    // sig_fit = sqrt(chisqr);   // [-Wunused_but_set_variable]
    pfit = pfit_b;

    if(fabs(sqrt(chisq1) - sqrt(chisqr)) < eps) {
      goto salir;
    }

  }

 salir: waterVaporRadiometer_.multiplySkyCouplingChannel(ichan, pfit);

  //    cout << "pfit=" << pfit << "  sky_coupling: " << waterVaporRadiometer_.getSkyCoupling()[0]   << endl;


}

double SkyStatus::sigmaSkyCouplingRetrieval_fromWVR(double par_fit,
                                                    const WaterVaporRadiometer &wvr,
                                                    vector<WVRMeasurement> &RadiometerData,
                                                    unsigned int n,
                                                    unsigned int m)
{

  vector<double> skyCoupling = wvr.getSkyCoupling();

  for(unsigned int i = 0; i < skyCoupling.size(); i++) {
    skyCoupling[i] = skyCoupling[i] * par_fit;
  }

  WVRMeasurement RadiometerData_withRetrieval;

  for(unsigned int i = n; i < m; i++) {

    RadiometerData_withRetrieval
        = mkWaterVaporRetrieval_fromWVR(RadiometerData[i].getmeasuredSkyBrightness(),
                                        wvr.getIdChannels(),
                                        skyCoupling,
                                        wvr.getsignalGain(),
                                        wvr.getSpilloverTemperature(),
                                        RadiometerData[i].getElevation());

    RadiometerData[i].setretrievedWaterVaporColumn(RadiometerData_withRetrieval.getretrievedWaterVaporColumn());
    RadiometerData[i].setfittedSkyBrightness(RadiometerData_withRetrieval.getfittedSkyBrightness());
    RadiometerData[i].setSigmaFit(RadiometerData_withRetrieval.getSigmaFit());

  }

  return getWVRAverageSigmaTskyFit(RadiometerData, n, m).get(Temperature::UnitKelvin);

}

double SkyStatus::sigmaSkyCouplingChannelRetrieval_fromWVR(double par_fit,
							   const WaterVaporRadiometer &wvr,
							   vector<WVRMeasurement> &RadiometerData,
							   unsigned int ichan,
							   unsigned int n,
							   unsigned int m)
{

  vector<double> skyCoupling = wvr.getSkyCoupling();

  //for(unsigned int i = 0; i < skyCoupling.size(); i++) {
  skyCoupling[ichan] = skyCoupling[ichan] * par_fit;
    //}

  WVRMeasurement RadiometerData_withRetrieval;

  for(unsigned int i = n; i < m; i++) {

    RadiometerData_withRetrieval
        = mkWaterVaporRetrieval_fromWVR(RadiometerData[i].getmeasuredSkyBrightness(),
                                        wvr.getIdChannels(),
                                        skyCoupling,
                                        wvr.getsignalGain(),
                                        wvr.getSpilloverTemperature(),
                                        RadiometerData[i].getElevation());

    RadiometerData[i].setretrievedWaterVaporColumn(RadiometerData_withRetrieval.getretrievedWaterVaporColumn());
    RadiometerData[i].setfittedSkyBrightness(RadiometerData_withRetrieval.getfittedSkyBrightness());
    RadiometerData[i].setSigmaFit(RadiometerData_withRetrieval.getSigmaFit());

  }

  return getWVRAverageSigmaTskyFit(RadiometerData, n, m).get(Temperature::UnitKelvin);
  // return getWVRSigmaChannelTskyFit(RadiometerData, ichan, n, m).get(Temperature::UnitKelvin);

}

void SkyStatus::rmSkyStatus()
{

  /* CODE À ÉCRIRE */

}
ATM_NAMESPACE_END

