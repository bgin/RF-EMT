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
 * "@(#) $Id: ATMRefractiveIndexProfile.cpp Exp $"
 *
 * who       when      what
 * --------  --------  ----------------------------------------------
 * pardo     24/03/09  created
 */

#include "ATMRefractiveIndexProfile.h"

#include <iostream>
#include <math.h>
#include <string>
#include <vector>



ATM_NAMESPACE_BEGIN

// Constructors

RefractiveIndexProfile::RefractiveIndexProfile(const Frequency &freq,
                                               const AtmProfile &atmProfile) :
  AtmProfile(atmProfile), SpectralGrid(freq)
{
  mkRefractiveIndexProfile();
}

RefractiveIndexProfile::RefractiveIndexProfile(const SpectralGrid &spectralGrid,
                                               const AtmProfile &atmProfile) :
  AtmProfile(atmProfile), SpectralGrid(spectralGrid)
{
  mkRefractiveIndexProfile();
}

RefractiveIndexProfile::RefractiveIndexProfile(const RefractiveIndexProfile & a) : AtmProfile(a), SpectralGrid(a)
{
  //   std::cout<<"Enter RefractiveIndexProfile copy constructor version Fri May 20 00:59:47 CEST 2005"<<endl;

  // level AtmProfile

  // type_ = a.type_;
  // prLimit_ = a.prLimit_;
  /*
   v_hx_.reserve(a.v_hx_.size());
   v_px_.reserve(a.v_px_.size());
   v_tx_.reserve(a.v_tx_.size());
   for(unsigned int n=0; n<a.v_hx_.size(); n++){
   v_hx_.push_back(a.v_hx_[n]);
   v_px_.push_back(a.v_px_[n]);
   v_tx_.push_back(a.v_tx_[n]);
   }
   */

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

  // copy thresholds
  altitudeThreshold_ = a.altitudeThreshold_;
  groundPressureThreshold_ = a.groundPressureThreshold_;
  groundTemperatureThreshold_ = a.groundTemperatureThreshold_;
  tropoLapseRateThreshold_ = a.tropoLapseRateThreshold_;
  relativeHumidityThreshold_ = a.relativeHumidityThreshold_;
  wvScaleHeightThreshold_ = a.wvScaleHeightThreshold_;

  v_layerThickness_.reserve(numLayer_);
  v_layerTemperature_.reserve(numLayer_);
  v_layerWaterVapor_.reserve(numLayer_);
  v_layerCO_.reserve(numLayer_);
  v_layerO3_.reserve(numLayer_);
  v_layerN2O_.reserve(numLayer_);
  v_layerNO2_.reserve(numLayer_);
  v_layerSO2_.reserve(numLayer_);


  for(unsigned int n = 0; n < numLayer_; n++) {
    v_layerThickness_.push_back(a.v_layerThickness_[n]);
    v_layerTemperature_.push_back(a.v_layerTemperature_[n]);
    //cout << "n=" << n << std::endl;
    v_layerWaterVapor_.push_back(a.v_layerWaterVapor_[n]);
    v_layerPressure_.push_back(a.v_layerPressure_[n]);
    v_layerCO_.push_back(a.v_layerCO_[n]);
    v_layerO3_.push_back(a.v_layerO3_[n]);
    v_layerN2O_.push_back(a.v_layerN2O_[n]);
    v_layerNO2_.push_back(a.v_layerNO2_[n]);
    v_layerSO2_.push_back(a.v_layerSO2_[n]);
  }

  // level Spectral Grid
  freqUnits_ = a.freqUnits_;
  v_chanFreq_ = a.v_chanFreq_;

  v_numChan_ = a.v_numChan_;
  v_refChan_ = a.v_refChan_;
  v_refFreq_ = a.v_refFreq_;
  v_chanSep_ = a.v_chanSep_;
  v_maxFreq_ = a.v_maxFreq_;
  v_minFreq_ = a.v_minFreq_;
  v_intermediateFrequency_ = a.v_intermediateFrequency_;
  v_loFreq_ = a.v_loFreq_;

  v_sidebandSide_ = a.v_sidebandSide_;
  v_sidebandType_ = a.v_sidebandType_;

  vv_assocSpwId_ = a.vv_assocSpwId_;
  vv_assocNature_ = a.vv_assocNature_;

  v_transfertId_ = a.v_transfertId_;

  // level Absorption Profile
  vv_N_H2OLinesPtr_.reserve(a.v_chanFreq_.size());
  vv_N_H2OContPtr_.reserve(a.v_chanFreq_.size());
  vv_N_O2LinesPtr_.reserve(a.v_chanFreq_.size());
  vv_N_DryContPtr_.reserve(a.v_chanFreq_.size());
  vv_N_O3LinesPtr_.reserve(a.v_chanFreq_.size());
  vv_N_COLinesPtr_.reserve(a.v_chanFreq_.size());
  vv_N_N2OLinesPtr_.reserve(a.v_chanFreq_.size());
  vv_N_NO2LinesPtr_.reserve(a.v_chanFreq_.size());
  vv_N_SO2LinesPtr_.reserve(a.v_chanFreq_.size());

  std::vector<std::complex<double> >* v_N_H2OLinesPtr;
  std::vector<std::complex<double> >* v_N_H2OContPtr;
  std::vector<std::complex<double> >* v_N_O2LinesPtr;
  std::vector<std::complex<double> >* v_N_DryContPtr;
  std::vector<std::complex<double> >* v_N_O3LinesPtr;
  std::vector<std::complex<double> >* v_N_COLinesPtr;
  std::vector<std::complex<double> >* v_N_N2OLinesPtr;
  std::vector<std::complex<double> >* v_N_NO2LinesPtr;
  std::vector<std::complex<double> >* v_N_SO2LinesPtr;

  for(unsigned int nc = 0; nc < v_chanFreq_.size(); nc++) {

    v_N_H2OLinesPtr = new std::vector<std::complex<double> > ;
    v_N_H2OLinesPtr->reserve(numLayer_);
    v_N_H2OContPtr = new std::vector<std::complex<double> > ;
    v_N_H2OContPtr->reserve(numLayer_);
    v_N_O2LinesPtr = new std::vector<std::complex<double> > ;
    v_N_O2LinesPtr->reserve(numLayer_);
    v_N_DryContPtr = new std::vector<std::complex<double> > ;
    v_N_DryContPtr->reserve(numLayer_);
    v_N_O3LinesPtr = new std::vector<std::complex<double> > ;
    v_N_O3LinesPtr->reserve(numLayer_);
    v_N_COLinesPtr = new std::vector<std::complex<double> > ;
    v_N_COLinesPtr->reserve(numLayer_);
    v_N_N2OLinesPtr = new std::vector<std::complex<double> > ;
    v_N_N2OLinesPtr->reserve(numLayer_);
    v_N_NO2LinesPtr = new std::vector<std::complex<double> > ;
    v_N_NO2LinesPtr->reserve(numLayer_);
    v_N_SO2LinesPtr = new std::vector<std::complex<double> > ;
    v_N_SO2LinesPtr->reserve(numLayer_);

    for(unsigned int n = 0; n < numLayer_; n++) {

      // std::cout << "numLayer_=" << nc << " " << n << std::endl; // COMMENTED OUT BY JUAN MAY/16/2005

      v_N_H2OLinesPtr->push_back(a.vv_N_H2OLinesPtr_[nc]->at(n));
      v_N_H2OContPtr->push_back(a.vv_N_H2OContPtr_[nc]->at(n));
      v_N_O2LinesPtr->push_back(a.vv_N_O2LinesPtr_[nc]->at(n));
      v_N_DryContPtr->push_back(a.vv_N_DryContPtr_[nc]->at(n));
      v_N_O3LinesPtr->push_back(a.vv_N_O3LinesPtr_[nc]->at(n));
      v_N_COLinesPtr->push_back(a.vv_N_COLinesPtr_[nc]->at(n));
      v_N_N2OLinesPtr->push_back(a.vv_N_N2OLinesPtr_[nc]->at(n));
      v_N_NO2LinesPtr->push_back(a.vv_N_NO2LinesPtr_[nc]->at(n));
      v_N_SO2LinesPtr->push_back(a.vv_N_SO2LinesPtr_[nc]->at(n));

    }

    vv_N_H2OLinesPtr_.push_back(v_N_H2OLinesPtr);
    vv_N_H2OContPtr_.push_back(v_N_H2OContPtr);
    vv_N_O2LinesPtr_.push_back(v_N_O2LinesPtr);
    vv_N_DryContPtr_.push_back(v_N_DryContPtr);
    vv_N_O3LinesPtr_.push_back(v_N_O3LinesPtr);
    vv_N_COLinesPtr_.push_back(v_N_COLinesPtr);
    vv_N_N2OLinesPtr_.push_back(v_N_N2OLinesPtr);
    vv_N_NO2LinesPtr_.push_back(v_N_NO2LinesPtr);
    vv_N_SO2LinesPtr_.push_back(v_N_SO2LinesPtr);

  }

}

RefractiveIndexProfile::RefractiveIndexProfile()
{
}

RefractiveIndexProfile::~RefractiveIndexProfile()
{
  rmRefractiveIndexProfile();
}

void RefractiveIndexProfile::rmRefractiveIndexProfile()
{
  // for every frequency channel delete the pointer to the absorption profile
  for(unsigned int nc = 0; nc < v_chanFreq_.size(); nc++) {
    delete vv_N_H2OLinesPtr_[nc];
    delete vv_N_H2OContPtr_[nc];
    delete vv_N_O2LinesPtr_[nc];
    delete vv_N_DryContPtr_[nc];
    delete vv_N_O3LinesPtr_[nc];
    delete vv_N_COLinesPtr_[nc];
    delete vv_N_N2OLinesPtr_[nc];
    delete vv_N_NO2LinesPtr_[nc];
    delete vv_N_SO2LinesPtr_[nc];
  }
}

bool RefractiveIndexProfile::updateRefractiveIndexProfile(const Length &altitude,
                                                          const Pressure &groundPressure,
                                                          const Temperature &groundTemperature,
                                                          double tropoLapseRate,
                                                          const Humidity &relativeHumidity,
                                                          const Length &wvScaleHeight)
{

  bool updated = false;
  bool mkNewAtmProfile = updateAtmProfile(altitude,
                                          groundPressure,
                                          groundTemperature,
                                          tropoLapseRate,
                                          relativeHumidity,
                                          wvScaleHeight);
  unsigned int numLayer = getNumLayer();

  if(vv_N_H2OLinesPtr_.size() < v_chanFreq_.size()) {
    mkNewAtmProfile = true;
    std::cout << " RefractiveIndexProfile: number of spectral windows has increased"
        << std::endl;
  }
  if(mkNewAtmProfile) {
    if(numLayer) {
      mkRefractiveIndexProfile();
      updated = true;
    } else {
      std::cout << " RefractiveIndexProfile: ERROR:  getNumLayer() returns 0"
          << std::endl;
    }
  }
  return updated;
}

// NB: this interface is required because the sub-class .... overrides this method.
bool RefractiveIndexProfile::setBasicAtmosphericParameters(const Length &altitude,
                                                           const Pressure &groundPressure,
                                                           const Temperature &groundTemperature,
                                                           double tropoLapseRate,
                                                           const Humidity &relativeHumidity,
                                                           const Length &wvScaleHeight)
{
  bool updated = updateRefractiveIndexProfile(altitude,
                                              groundPressure,
                                              groundTemperature,
                                              tropoLapseRate,
                                              relativeHumidity,
                                              wvScaleHeight);
  return updated;
}

void RefractiveIndexProfile::mkRefractiveIndexProfile()
{

  //    static const double abun_18o=0.0020439;
  //    static const double abun_17o=0.0003750;
  //    static const double abun_D=0.000298444;
  //    static const double o2_mixing_ratio=0.2092;
  //    static const double mmol_h2o=18.005059688;  //   20*0.0020439+19*(0.0003750+2*0.000298444)+18*(1-0.0020439-0.0003750-2*0.000298444)

  // static bool first = true;   // [-Wunused_but_set_variable]

  double abun_O3, abun_CO, abun_N2O, abun_NO2, abun_SO2;
  double wvt, wv;
  // double t; // [-Wunused_but_set_variable]
  double nu;
  // double nu2, nu_pi; // [-Wunused_but_set_variable]
  // double width;
  // unsigned int npoints;
  RefractiveIndex atm;
  //    double sumAbsO3Lines1, sumAbsCOLines1, sumAbsN2OLines1, sumAbsNO2Lines1, sumAbsSO2Lines1;


  //TODO we will have to put numLayer_ and v_chanFreq_.size() const
  //we do not want to resize! ==> pas de setter pour SpectralGrid

  //cout << "vv_N_H2OLinesPtr_.size()=" << vv_N_H2OLinesPtr_.size() << std::endl;
  if(vv_N_H2OLinesPtr_.size() == 0) { // first time
    vv_N_H2OLinesPtr_.reserve(v_chanFreq_.size());
    vv_N_H2OContPtr_.reserve(v_chanFreq_.size());
    vv_N_O2LinesPtr_.reserve(v_chanFreq_.size());
    vv_N_DryContPtr_.reserve(v_chanFreq_.size());
    vv_N_O3LinesPtr_.reserve(v_chanFreq_.size());
    vv_N_COLinesPtr_.reserve(v_chanFreq_.size());
    vv_N_N2OLinesPtr_.reserve(v_chanFreq_.size());
    vv_N_NO2LinesPtr_.reserve(v_chanFreq_.size());
    vv_N_SO2LinesPtr_.reserve(v_chanFreq_.size());
  } else {
    if(vv_N_H2OLinesPtr_.size() == v_chanFreq_.size()) // there are new basic param
    rmRefractiveIndexProfile(); // delete all the layer profiles for all the frequencies
  }

  std::vector<std::complex<double> >* v_N_H2OLinesPtr;
  std::vector<std::complex<double> >* v_N_H2OContPtr;
  std::vector<std::complex<double> >* v_N_O2LinesPtr;
  std::vector<std::complex<double> >* v_N_DryContPtr;
  std::vector<std::complex<double> >* v_N_O3LinesPtr;
  std::vector<std::complex<double> >* v_N_COLinesPtr;
  std::vector<std::complex<double> >* v_N_N2OLinesPtr;
  std::vector<std::complex<double> >* v_N_NO2LinesPtr;
  std::vector<std::complex<double> >* v_N_SO2LinesPtr;

  // std::cout << "v_chanFreq_.size()=" << v_chanFreq_.size() << std::endl;
  // std::cout << "numLayer_=" << numLayer_ << std::endl;
  // std::cout << "v_chanFreq_[0]=" << v_chanFreq_[0] << std::endl;
  // check if new spectral windows have been added
  unsigned int ncmin;
  /*  std::cout << "vv_N_H2OLinesPtr_.size()="<<vv_N_H2OLinesPtr_.size()<<endl; */
  ncmin = vv_N_H2OLinesPtr_.size(); // will be > 0 if spectral window(s) have been added
  if(newBasicParam_) ncmin = 0;

  //    std::cout << "ncmin=" << ncmin << std::endl;

  for(unsigned int nc = ncmin; nc < v_chanFreq_.size(); nc++) {

    v_N_H2OLinesPtr = new std::vector<std::complex<double> > ;
    v_N_H2OContPtr = new std::vector<std::complex<double> > ;
    v_N_O2LinesPtr = new std::vector<std::complex<double> > ;
    v_N_DryContPtr = new std::vector<std::complex<double> > ;
    v_N_O3LinesPtr = new std::vector<std::complex<double> > ;
    v_N_COLinesPtr = new std::vector<std::complex<double> > ;
    v_N_N2OLinesPtr = new std::vector<std::complex<double> > ;
    v_N_NO2LinesPtr = new std::vector<std::complex<double> > ;
    v_N_SO2LinesPtr = new std::vector<std::complex<double> > ;
    v_N_H2OLinesPtr->reserve(numLayer_);
    v_N_H2OContPtr->reserve(numLayer_);
    v_N_O2LinesPtr->reserve(numLayer_);
    v_N_DryContPtr->reserve(numLayer_);
    v_N_O3LinesPtr->reserve(numLayer_);
    v_N_COLinesPtr->reserve(numLayer_);
    v_N_N2OLinesPtr->reserve(numLayer_);
    v_N_NO2LinesPtr->reserve(numLayer_);
    v_N_SO2LinesPtr->reserve(numLayer_);

    nu = 1.0E-9 * v_chanFreq_[nc]; // ATM uses GHz units

    // std::cout << "freq. points =" << v_chanFreq_.size() << std::endl;

    /*       TO BE IMPLEMENTED IN NEXT RELEASE

    if (v_chanFreq_.size()>1){
      if(nc==0){
	width = fabs(v_chanFreq_[nc+1]-v_chanFreq_[nc])*1e-9;       // width en GHz para ATM
	npoints=(unsigned int)round(width*100);                     // One point every 10 MHz
      }else{
	if(nc==v_chanFreq_.size()-1){
	  width = fabs(v_chanFreq_[nc]-v_chanFreq_[nc-1])*1e-9;     // width en GHz para ATM
	  npoints=(unsigned int)round(width*100);                   // One point every 10 MHz
	}else{
	  width = fabs((v_chanFreq_[nc+1]-v_chanFreq_[nc-1])/2.0)*1e-9;    // width en GHz para ATM
	  npoints=(unsigned int)round(width*100);                          // One point every 10 MHz
	}
      }
    }else{
      width = 0.001;      // default width = 1 MHz = 0.001 GHz
      npoints=1;
    }

    if(npoints==0){npoints=1;}

    */



    // std::cout << "nc =" << nc << " nu=" << nu << " width=" << width << " GHz    npoints=" << npoints << std::endl;

    // nu2 = nu * nu;      // [-Wunused_but_set_variable]
    // nu_pi = nu / M_PI;    // [-Wunused_but_set_variable]

    for(unsigned int j = 0; j < numLayer_; j++) {

      wv = v_layerWaterVapor_[j] * 1000.0; // se multiplica por 10**3 por cuestión de unidades en las rutinas fortran.
      wvt = wv * v_layerTemperature_[j] / 217.0; // v_layerWaterVapor_[j] está en kg/m**3
      // t = v_layerTemperature_[j] / 300.0;    // [-Wunused_but_set_variable]


      // std::cout <<"ATMRefractiveIndexProfile: " << v_layerTemperature_[j] << " K " << v_layerPressure_[j] << " mb "  << nu << " GHz " << std::endl;
      // std::cout <<"ATMRefractiveIndexProfile: O2" <<  atm.getRefractivity_o2(v_layerTemperature_[j],v_layerPressure_[j],wvt,nu) << std::endl;
      // std::cout << "ATMRefractiveIndexProfile: O2" << atm.getRefractivity_o2(v_layerTemperature_[j],v_layerPressure_[j],wvt,nu,width,npoints) << std::endl;
      // std::cout << "ATMRefractiveIndexProfile: H2O" << atm.getRefractivity_h2o(v_layerTemperature_[j],v_layerPressure_[j],wvt,nu) << std::endl;
      // std::cout << "ATMRefractiveIndexProfile: H2O" << atm.getRefractivity_h2o(v_layerTemperature_[j],v_layerPressure_[j],wvt,nu,width,npoints) << std::endl;
      // std::cout <<"ATMRefractiveIndexProfile: O3" <<  atm.getRefractivity_o3(v_layerTemperature_[j],v_layerPressure_[j],nu,v_layerO3_[j]) << std::endl;
      // std::cout << "ATMRefractiveIndexProfile: O3" << atm.getRefractivity_o3(v_layerTemperature_[j],v_layerPressure_[j],nu,width,npoints,v_layerO3_[j]) << std::endl;
      // std::cout <<"ATMRefractiveIndexProfile: CO" <<  atm.getSpecificRefractivity_co(v_layerTemperature_[j],v_layerPressure_[j],nu) << std::endl;
      // std::cout << "ATMRefractiveIndexProfile: CO" << atm.getSpecificRefractivity_co(v_layerTemperature_[j],v_layerPressure_[j],nu,width,npoints) << std::endl;

      v_N_O2LinesPtr->push_back(atm.getRefractivity_o2(v_layerTemperature_[j],
                                                       v_layerPressure_[j],
                                                       wvt,
                                                       nu));     // ,width,npoints)); TO BE IMPLEMENTED IN NEXT RELEASE

      std::complex<double> cont_h2o =
        atm.getSpecificRefractivity_cnth2o(v_layerTemperature_[j],
                                           v_layerPressure_[j],
                                           wvt,
                                           nu);                   // ,width,npoints); TO BE IMPLEMENTED IN NEXT RELEASE
      std::complex<double> cont_dry =
          atm.getSpecificRefractivity_cntdry(v_layerTemperature_[j],
                                             v_layerPressure_[j],
                                             wvt,
                                             nu);                  // ,width,npoints); TO BE IMPLEMENTED IN NEXT RELEASE

      v_N_H2OContPtr->push_back(cont_h2o);
      v_N_DryContPtr->push_back(cont_dry);

      if(v_layerWaterVapor_[j] > 0) {

        v_N_H2OLinesPtr->push_back(atm.getRefractivity_h2o(v_layerTemperature_[j],
                                                           v_layerPressure_[j],
                                                           wvt,
                                                           nu)); // ,width,npoints)); TO BE IMPLEMENTED IN NEXT RELEASE
      } else {
        v_N_H2OLinesPtr->push_back(0.0);
      }

      //	if(v_layerO3_[j]<0.0||j==10){cout << "v_layerO3_[" << j << "]=" << v_layerO3_[j] << std::endl;}

      if(v_layerO3_[j] > 0) {

        abun_O3 = v_layerO3_[j] * 1E-6;
        v_N_O3LinesPtr->push_back(atm.getRefractivity_o3(v_layerTemperature_[j],
                                                         v_layerPressure_[j],
                                                         nu,      // width,npoints, TO BE IMPLEMENTED IN NEXT RELEASE
                                                         abun_O3 * 1e6));
      } else {
        v_N_O3LinesPtr->push_back(0.0);
      }

      if(v_layerCO_[j] > 0) {
        abun_CO = v_layerCO_[j] * 1E-6; // in cm^-3
        v_N_COLinesPtr->push_back(atm.getSpecificRefractivity_co(v_layerTemperature_[j],
								 v_layerPressure_[j],
								 nu)            // ,width,npoints) TO BE IMPLEMENTED IN NEXT RELEASE
				  * abun_CO * 1e6); // m^2 * m^-3 = m^-1
        } else {
        v_N_COLinesPtr->push_back(0.0);
      }

      if(v_layerN2O_[j] > 0) {
        abun_N2O = v_layerN2O_[j] * 1E-6;
        v_N_N2OLinesPtr->push_back(atm.getSpecificRefractivity_n2o(v_layerTemperature_[j],
								   v_layerPressure_[j],
								   nu)             // ,width,npoints) TO BE IMPLEMENTED IN NEXT RELEASE
				   * abun_N2O * 1e6); // m^2 * m^-3 = m^-1
        } else {
        v_N_N2OLinesPtr->push_back(0.0);
      }

      if(v_layerNO2_[j] > 0) {
        abun_NO2 = v_layerNO2_[j] * 1E-6;
        v_N_NO2LinesPtr->push_back(atm.getSpecificRefractivity_no2(v_layerTemperature_[j],
								   v_layerPressure_[j],
								   nu)             // ,width,npoints) TO BE IMPLEMENTED IN NEXT RELEASE
				   * abun_NO2 * 1e6); // m^2 * m^-3 = m^-1
      } else {
        v_N_NO2LinesPtr->push_back(0.0);
      }

      if(v_layerSO2_[j] > 0) {
        abun_SO2 = v_layerSO2_[j] * 1E-6;
        v_N_SO2LinesPtr->push_back(atm.getSpecificRefractivity_so2(v_layerTemperature_[j],
								   v_layerPressure_[j],
								   nu)            // ,width,npoints) TO BE IMPLEMENTED IN NEXT RELEASE
				   * abun_SO2 * 1e6); // m^2 * m^-3 = m^-1
      } else {
        v_N_SO2LinesPtr->push_back(0.0);
      }
    }

    // if(vv_N_H2OLinesPtr_.size() == 0) first = true;  // [-Wunused_but_set_variable]

    if(vv_N_H2OLinesPtr_.size() < v_chanFreq_.size()) {
      vv_N_H2OLinesPtr_.push_back(v_N_H2OLinesPtr);
      vv_N_H2OContPtr_.push_back(v_N_H2OContPtr);
      vv_N_O2LinesPtr_.push_back(v_N_O2LinesPtr);
      vv_N_DryContPtr_.push_back(v_N_DryContPtr);
      vv_N_O3LinesPtr_.push_back(v_N_O3LinesPtr);
      vv_N_COLinesPtr_.push_back(v_N_COLinesPtr);
      vv_N_N2OLinesPtr_.push_back(v_N_N2OLinesPtr);
      vv_N_NO2LinesPtr_.push_back(v_N_NO2LinesPtr);
      vv_N_SO2LinesPtr_.push_back(v_N_SO2LinesPtr);
    } else {
      vv_N_H2OLinesPtr_[nc] = v_N_H2OLinesPtr;
      vv_N_H2OContPtr_[nc] = v_N_H2OContPtr;
      vv_N_O2LinesPtr_[nc] = v_N_O2LinesPtr;
      vv_N_DryContPtr_[nc] = v_N_DryContPtr;
      vv_N_O3LinesPtr_[nc] = v_N_O3LinesPtr;
      vv_N_COLinesPtr_[nc] = v_N_COLinesPtr;
      vv_N_N2OLinesPtr_[nc] = v_N_N2OLinesPtr;
      vv_N_NO2LinesPtr_[nc] = v_N_NO2LinesPtr;
      vv_N_SO2LinesPtr_[nc] = v_N_SO2LinesPtr;
    }

  }

  newBasicParam_ = false;
  // first = false;  // [-Wunused_but_set_variable]
}


Opacity RefractiveIndexProfile::getDryOpacityUpTo(unsigned int nc, Length refalti)
{
  unsigned int ires; unsigned int numlayerold; Length alti;  double fractionLast;
  Opacity opacityout0; Opacity opacityout1; Opacity zeroOp(0.0,Opacity::UnitNeper);

  if(refalti.get(Length::UnitKiloMeter) <= altitude_.get(Length::UnitKiloMeter)) {
    return zeroOp;
  }else{
      fractionLast = 1.0; numlayerold = numLayer_;
      opacityout0=getDryOpacity(nc); ires=numlayerold-1; alti=altitude_;
      for(unsigned int i=0; i<numLayer_; i++){
	if(alti.get(Length::UnitKiloMeter) < refalti.get(Length::UnitKiloMeter) &&  (alti.get(Length::UnitKiloMeter)+v_layerThickness_[i]/1000.0) >= refalti.get(Length::UnitKiloMeter))
	  { ires=i; fractionLast = (refalti.get(Length::UnitMeter)-alti.get(Length::UnitMeter))/v_layerThickness_[i]; }
	alti = alti + Length(v_layerThickness_[i],Length::UnitMeter);
      }
      numLayer_ = ires;
      opacityout0=getDryOpacity(nc);
      numLayer_ = ires+1;
      opacityout1=getDryOpacity(nc);
      numLayer_ = numlayerold;
      return opacityout0+(opacityout1-opacityout0)*fractionLast;
  }
}
Opacity RefractiveIndexProfile::getDryOpacity(unsigned int nc)
{
  if(!chanIndexIsValid(nc)) return Opacity(-999.0);
  double kv = 0;
  for(unsigned int j = 0; j < numLayer_; j++) {
    kv = kv + imag(vv_N_O2LinesPtr_[nc]->at(j) + vv_N_DryContPtr_[nc]->at(j)
		   + vv_N_O3LinesPtr_[nc]->at(j)  + vv_N_COLinesPtr_[nc]->at(j)
		   + vv_N_N2OLinesPtr_[nc]->at(j) + vv_N_NO2LinesPtr_[nc]->at(j)
		   + vv_N_SO2LinesPtr_[nc]->at(j)) * v_layerThickness_[j];
  }
  return Opacity(kv);
}


Opacity RefractiveIndexProfile::getAverageDryOpacity(unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) return Opacity(-999.0);
  Opacity totalaverage;
  totalaverage = Opacity(0.0, Opacity::UnitNeper);
  for(unsigned int nc = 0; nc < getNumChan(spwid); nc++) {
    totalaverage = totalaverage + getDryOpacity(spwid, nc);
  }
  totalaverage = totalaverage / getNumChan(spwid);
  return totalaverage;
}

Opacity RefractiveIndexProfile::getAverageO2LinesOpacity(unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) return Opacity(-999.0);
  Opacity totalaverage;
  totalaverage = Opacity(0.0, Opacity::UnitNeper);
  for(unsigned int nc = 0; nc < getNumChan(spwid); nc++) {
    totalaverage = totalaverage + getO2LinesOpacity(spwid, nc);
  }
  totalaverage = totalaverage / getNumChan(spwid);
  return totalaverage;
}

Opacity RefractiveIndexProfile::getAverageO3LinesOpacity(unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) return Opacity(-999.0);
  Opacity totalaverage;
  totalaverage = Opacity(0.0, Opacity::UnitNeper);
  for(unsigned int nc = 0; nc < getNumChan(spwid); nc++) {
    /*      std::cout << " Freq = " << getChanFreq(spwid,nc).get(Frequency::UnitGigaHertz)
     << " O3 opacity = " << getO3LinesOpacity(spwid,nc).get(Opacity::UnitNeper)
     << " O3 pathlength = " << getO3LinesPathLength(spwid,nc).get(Length::Microns)
     << " O2 opacity = " << getO2LinesOpacity(spwid,nc).get(Opacity::UnitNeper)
     << " O2 pathlength = " << getO2LinesPathLength(spwid,nc).get(Length::Microns)
     << std::endl; */
    totalaverage = totalaverage + getO3LinesOpacity(spwid, nc);
  }
  totalaverage = totalaverage / getNumChan(spwid);
  return totalaverage;
}

Opacity RefractiveIndexProfile::getAverageN2OLinesOpacity(unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) return Opacity(-999.0);
  Opacity totalaverage;
  totalaverage = Opacity(0.0, Opacity::UnitNeper);
  for(unsigned int nc = 0; nc < getNumChan(spwid); nc++) {
    totalaverage = totalaverage + getN2OLinesOpacity(spwid, nc);
  }
  totalaverage = totalaverage / getNumChan(spwid);
  return totalaverage;
}

Opacity RefractiveIndexProfile::getAverageNO2LinesOpacity(unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) return Opacity(-999.0);
  Opacity totalaverage;
  totalaverage = Opacity(0.0, Opacity::UnitNeper);
  for(unsigned int nc = 0; nc < getNumChan(spwid); nc++) {
    totalaverage = totalaverage + getNO2LinesOpacity(spwid, nc);
  }
  totalaverage = totalaverage / getNumChan(spwid);
  return totalaverage;
}

Opacity RefractiveIndexProfile::getAverageSO2LinesOpacity(unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) return Opacity(-999.0);
  Opacity totalaverage;
  totalaverage = Opacity(0.0, Opacity::UnitNeper);
  for(unsigned int nc = 0; nc < getNumChan(spwid); nc++) {
    totalaverage = totalaverage + getSO2LinesOpacity(spwid, nc);
  }
  totalaverage = totalaverage / getNumChan(spwid);
  return totalaverage;
}




Opacity RefractiveIndexProfile::getAverageCOLinesOpacity(unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) return Opacity(-999.0);
  Opacity totalaverage;
  totalaverage = Opacity(0.0, Opacity::UnitNeper);
  for(unsigned int nc = 0; nc < getNumChan(spwid); nc++) {
    totalaverage = totalaverage + getCOLinesOpacity(spwid, nc);
  }
  totalaverage = totalaverage / getNumChan(spwid);
  return totalaverage;
}

Opacity RefractiveIndexProfile::getAverageDryContOpacity(unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) return Opacity(-999.0);
  Opacity totalaverage;
  totalaverage = Opacity(0.0, Opacity::UnitNeper);
  for(unsigned int nc = 0; nc < getNumChan(spwid); nc++) {
    totalaverage = totalaverage + getDryContOpacity(spwid, nc);
  }
  totalaverage = totalaverage / getNumChan(spwid);
  return totalaverage;
}

Opacity RefractiveIndexProfile::getDryContOpacity()
{
  return getDryContOpacity(0);
}

Opacity RefractiveIndexProfile::getDryContOpacity(unsigned int nc)
{
  if(!chanIndexIsValid(nc)) return Opacity(-999.0);
  double kv = 0;
  for(unsigned int j = 0; j < numLayer_; j++) {
    kv = kv + imag(vv_N_DryContPtr_[nc]->at(j)) * v_layerThickness_[j];
  }
  return Opacity(kv);
}

Opacity RefractiveIndexProfile::getDryContOpacity(unsigned int spwid,
                                                  unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) return Opacity(-999.0);
  return getDryContOpacity(v_transfertId_[spwid] + nc);
}

Opacity RefractiveIndexProfile::getO2LinesOpacity()
{
  return getO2LinesOpacity(0);
}

Opacity RefractiveIndexProfile::getO2LinesOpacity(unsigned int nc)
{
  if(!chanIndexIsValid(nc)) return Opacity(-999.0);
  double kv = 0;
  for(unsigned int j = 0; j < numLayer_; j++) {
    kv = kv + imag(vv_N_O2LinesPtr_[nc]->at(j)) * v_layerThickness_[j];
  }
  return Opacity(kv);
}

Opacity RefractiveIndexProfile::getO2LinesOpacity(unsigned int spwid,
                                                  unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) return Opacity(-999.0);
  return getO2LinesOpacity(v_transfertId_[spwid] + nc);
}

Opacity RefractiveIndexProfile::getCOLinesOpacity()
{
  return getCOLinesOpacity(0);
}

Opacity RefractiveIndexProfile::getCOLinesOpacity(unsigned int nc)
{
  if(!chanIndexIsValid(nc)) return Opacity(-999.0);
  double kv = 0;
  for(unsigned int j = 0; j < numLayer_; j++) {
    kv = kv + imag(vv_N_COLinesPtr_[nc]->at(j)) * v_layerThickness_[j];
  }
  return Opacity(kv);
}

Opacity RefractiveIndexProfile::getCOLinesOpacity(unsigned int spwid,
                                                  unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) return Opacity(-999.0);
  return getCOLinesOpacity(v_transfertId_[spwid] + nc);
}

Opacity RefractiveIndexProfile::getN2OLinesOpacity()
{
  return getN2OLinesOpacity(0);
}

Opacity RefractiveIndexProfile::getNO2LinesOpacity()
{
  return getNO2LinesOpacity(0);
}

Opacity RefractiveIndexProfile::getSO2LinesOpacity()
{
  return getSO2LinesOpacity(0);
}

Opacity RefractiveIndexProfile::getN2OLinesOpacity(unsigned int nc)
{
  if(!chanIndexIsValid(nc)) return Opacity(-999.0);
  double kv = 0;
  for(unsigned int j = 0; j < numLayer_; j++) {
    kv = kv + imag(vv_N_N2OLinesPtr_[nc]->at(j)) * v_layerThickness_[j];
  }
  return Opacity(kv);
}

Opacity RefractiveIndexProfile::getN2OLinesOpacity(unsigned int spwid,
                                                   unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) return Opacity(-999.0);
  return getN2OLinesOpacity(v_transfertId_[spwid] + nc);
}

Opacity RefractiveIndexProfile::getNO2LinesOpacity(unsigned int nc)
{
  if(!chanIndexIsValid(nc)) return Opacity(-999.0);
  double kv = 0;
  for(unsigned int j = 0; j < numLayer_; j++) {
    kv = kv + imag(vv_N_NO2LinesPtr_[nc]->at(j)) * v_layerThickness_[j];
  }
  return Opacity(kv);
}

Opacity RefractiveIndexProfile::getNO2LinesOpacity(unsigned int spwid,
                                                   unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) return Opacity(-999.0);
  return getNO2LinesOpacity(v_transfertId_[spwid] + nc);
}


Opacity RefractiveIndexProfile::getSO2LinesOpacity(unsigned int nc)
{
  if(!chanIndexIsValid(nc)) return Opacity(-999.0);
  double kv = 0;
  for(unsigned int j = 0; j < numLayer_; j++) {
    kv = kv + imag(vv_N_SO2LinesPtr_[nc]->at(j)) * v_layerThickness_[j];
  }
  return Opacity(kv);
}

Opacity RefractiveIndexProfile::getSO2LinesOpacity(unsigned int spwid,
                                                   unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) return Opacity(-999.0);
  return getSO2LinesOpacity(v_transfertId_[spwid] + nc);
}


Opacity RefractiveIndexProfile::getO3LinesOpacity()
{
  return getO3LinesOpacity(0);
}

Opacity RefractiveIndexProfile::getO3LinesOpacity(unsigned int nc)
{
  if(!chanIndexIsValid(nc)) return Opacity(-999.0);
  double kv = 0;
  for(unsigned int j = 0; j < numLayer_; j++) {
    kv = kv + imag(vv_N_O3LinesPtr_[nc]->at(j)) * v_layerThickness_[j];
  }
  return Opacity(kv);
}

Opacity RefractiveIndexProfile::getO3LinesOpacity(unsigned int spwid,
                                                  unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) return Opacity(-999.0);
  return getO3LinesOpacity(v_transfertId_[spwid] + nc);
}

Opacity RefractiveIndexProfile::getWetOpacity(const Length &integratedwatercolumn)
{
//   std::cout << "1 integratedwatercolumn.get()="  << integratedwatercolumn.get() << std::endl;
//   std::cout << "2 getGroundWH2O().get()="   << getGroundWH2O().get() << std::endl;
//   std::cout << "3 getWetOpacity()="  << std::endl;
//   std::cout << "4 " << std::endl;
  return getWetOpacity(getGroundWH2O(),0)*(integratedwatercolumn.get()/getGroundWH2O().get());

  // 2010_SEP02: return getWetOpacity(integratedwatercolumn,0)*(integratedwatercolumn.get()/getGroundWH2O().get());
}

Opacity RefractiveIndexProfile::getWetOpacity(const Length &integratedwatercolumn,
                                              unsigned int nc)
{
  if(!chanIndexIsValid(nc)) return Opacity(-999.0);
  double kv = 0;
  /*  std::cout<<"nc="<<nc<<endl; */
  for(unsigned int j = 0; j < numLayer_; j++) {
    kv = kv + imag(vv_N_H2OLinesPtr_[nc]->at(j) + vv_N_H2OContPtr_[nc]->at(j))
        * v_layerThickness_[j];

  }
  return Opacity(kv*(integratedwatercolumn.get()/getGroundWH2O().get()));
}

Opacity RefractiveIndexProfile::getWetOpacity(const Length & integratedwatercolumn,
                                              unsigned int spwid,
                                              unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) return Opacity(-999.0);
  return getWetOpacity(integratedwatercolumn,v_transfertId_[spwid] + nc);
}

Opacity RefractiveIndexProfile::getAverageWetOpacity(const Length &integratedwatercolumn,
                                                     unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) return Opacity(-999.0);
  Opacity totalaverage;
  totalaverage = Opacity(0.0, Opacity::UnitNeper);
  for(unsigned int nc = 0; nc < getNumChan(spwid); nc++) {
    totalaverage = totalaverage + getWetOpacity(integratedwatercolumn, spwid, nc);
  }
  totalaverage = totalaverage / getNumChan(spwid);
  return totalaverage;
}

Opacity RefractiveIndexProfile::getH2OLinesOpacity(const Length &integratedwatercolumn)
{
  return getH2OLinesOpacity(integratedwatercolumn,0);
}

Opacity RefractiveIndexProfile::getH2OLinesOpacity(const Length &integratedwatercolumn,
                                                   unsigned int nc)
{
  if(!chanIndexIsValid(nc)) return Opacity(-999.0);
  double kv = 0;
  for(unsigned int j = 0; j < numLayer_; j++) {
    /*    std::cout <<"j="<<j<<" abs H2O Lines ="<<vv_N_H2OLinesPtr_[nc]->at(j) <<endl; */
    kv = kv + imag(vv_N_H2OLinesPtr_[nc]->at(j)) * v_layerThickness_[j];
  }
  return Opacity(kv*(integratedwatercolumn.get()/getGroundWH2O().get()));
}

Opacity RefractiveIndexProfile::getH2OLinesOpacity(const Length &integratedwatercolumn,
                                                   unsigned int spwid,
                                                   unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) return Opacity(-999.0);
  return getH2OLinesOpacity(integratedwatercolumn,v_transfertId_[spwid] + nc);
}

Opacity RefractiveIndexProfile::getAverageH2OLinesOpacity(const Length &integratedwatercolumn,
                                                          unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) return Opacity(-999.0);
  Opacity totalaverage;
  totalaverage = Opacity(0.0, Opacity::UnitNeper);
  for(unsigned int nc = 0; nc < getNumChan(spwid); nc++) {
    totalaverage = totalaverage + getH2OLinesOpacity(integratedwatercolumn,spwid, nc);
  }
  totalaverage = totalaverage / getNumChan(spwid);
  return totalaverage;
}


Opacity RefractiveIndexProfile::getH2OContOpacity(const Length &integratedwatercolumn)
{
  return getH2OContOpacity(integratedwatercolumn,0);
}


Opacity RefractiveIndexProfile::getH2OContOpacity(const Length &integratedwatercolumn,
                                                  unsigned int nc)
{
  if(!chanIndexIsValid(nc)) return Opacity(-999.0);
  double kv = 0;
  for(unsigned int j = 0; j < numLayer_; j++) {
    kv = kv + imag(vv_N_H2OContPtr_[nc]->at(j)) * v_layerThickness_[j];
  }
  return Opacity(kv*(integratedwatercolumn.get()/getGroundWH2O().get()));
}


Opacity RefractiveIndexProfile::getH2OContOpacity(const Length &integratedwatercolumn,
                                                  unsigned int spwid,
                                                  unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) return Opacity(-999.0);
  return getH2OContOpacity(integratedwatercolumn,v_transfertId_[spwid] + nc);
}



Opacity RefractiveIndexProfile::getAverageH2OContOpacity(const Length &integratedwatercolumn,
                                                         unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) return Opacity(-999.0);
  Opacity totalaverage;
  totalaverage = Opacity(0.0, Opacity::UnitNeper);
  for(unsigned int nc = 0; nc < getNumChan(spwid); nc++) {
    totalaverage = totalaverage + getH2OContOpacity(integratedwatercolumn,spwid, nc);
  }
  totalaverage = totalaverage / getNumChan(spwid);
  return totalaverage;
}

Angle RefractiveIndexProfile::getDispersiveH2OPhaseDelay(const Length &integratedwatercolumn)
{
  return getDispersiveH2OPhaseDelay(integratedwatercolumn,0);
}

Length RefractiveIndexProfile::getDispersiveH2OPathLength(const Length &integratedwatercolumn)
{
  return getDispersiveH2OPathLength(integratedwatercolumn,0);
}

Angle RefractiveIndexProfile::getDispersiveH2OPhaseDelay(const Length &integratedwatercolumn,
                                                         unsigned int nc)
{
  if(!chanIndexIsValid(nc)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  double kv = 0;
  for(unsigned int j = 0; j < numLayer_; j++) {
    kv = kv + real(vv_N_H2OLinesPtr_[nc]->at(j)) * v_layerThickness_[j];
  }
  Angle aa(kv*(integratedwatercolumn.get()/getGroundWH2O().get())* 57.29578, Angle::UnitDegree);
  return aa;
}

Length RefractiveIndexProfile::getDispersiveH2OPathLength(const Length &integratedwatercolumn,
                                                          unsigned int nc)
{
  if(!chanIndexIsValid(nc)) {
    return Length(-999.0, Length::UnitMeter);
  }
  double wavelength = 299792458.0 / v_chanFreq_[nc]; // in m
  Length ll((wavelength / 360.0) * getDispersiveH2OPhaseDelay(integratedwatercolumn,nc).get(Angle::UnitDegree),
            Length::UnitMeter);
  return ll;
}

Angle RefractiveIndexProfile::getDispersiveH2OPhaseDelay(const Length &integratedwatercolumn,
                                                         unsigned int spwid,
                                                         unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  return getDispersiveH2OPhaseDelay(integratedwatercolumn,v_transfertId_[spwid] + nc);
}

Angle RefractiveIndexProfile::getAverageDispersiveH2OPhaseDelay(const Length &integratedwatercolumn,
                                                                unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  double av = 0.0;
  for(unsigned int i = 0; i < getNumChan(spwid); i++) {
    av = av + getDispersiveH2OPhaseDelay(integratedwatercolumn,v_transfertId_[spwid] + i).get(Angle::UnitDegree);
  }
  av = av / getNumChan(spwid);
  Angle average(av, Angle::UnitDegree);
  return average;
}

Length RefractiveIndexProfile::getDispersiveH2OPathLength(const Length &integratedwatercolumn,
                                                          unsigned int spwid,
                                                          unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) {
    return Length(-999.0, Length::UnitMeter);
  }
  return getDispersiveH2OPathLength(integratedwatercolumn,v_transfertId_[spwid] + nc);
}

Length RefractiveIndexProfile::getAverageDispersiveH2OPathLength(const Length &integratedwatercolumn,
                                                                 unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) {
    return Length(-999.0, Length::UnitMeter);
  }
  double av = 0.0;
  for(unsigned int i = 0; i < getNumChan(spwid); i++) {
    av = av + getDispersiveH2OPathLength(integratedwatercolumn,v_transfertId_[spwid] + i).get(Length::UnitMilliMeter);
  }
  av = av / getNumChan(spwid);
  Length average(av, Length::UnitMilliMeter);
  return average;
}

Angle RefractiveIndexProfile::getNonDispersiveDryPhaseDelay()
{
  return getNonDispersiveDryPhaseDelay(0);
}

Length RefractiveIndexProfile::getNonDispersiveDryPathLength()
{
  return getNonDispersiveDryPathLength(0);
}

Angle RefractiveIndexProfile::getDispersiveDryPhaseDelay()
{
  return getDispersiveDryPhaseDelay(0);
}

Length RefractiveIndexProfile::getDispersiveDryPathLength()
{
  return getDispersiveDryPathLength(0);
}

Angle RefractiveIndexProfile::getNonDispersiveDryPhaseDelay(unsigned int nc)
{
  if(!chanIndexIsValid(nc)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  double kv = 0;
  for(unsigned int j = 0; j < numLayer_; j++) {
    kv = kv + real(vv_N_DryContPtr_[nc]->at(j)) * v_layerThickness_[j];
  }
  Angle aa(kv * 57.29578, Angle::UnitDegree);
  return aa;
}

Angle RefractiveIndexProfile::getDispersiveDryPhaseDelay(unsigned int nc)
{
  //    std::cout << "getO2LinesPhaseDelay(" << nc << ")=" << getO2LinesPhaseDelay(nc).get(Angle::UnitDegree)  << std::endl;
  // std::cout << "getO3LinesPhaseDelay(" << nc << ")=" << getO3LinesPhaseDelay(nc).get(Angle::UnitDegree) << std::endl;
  //    std::cout << "getN2OLinesPhaseDelay(" << nc << ")=" << getN2OLinesPhaseDelay(nc).get(Angle::UnitDegree) << std::endl;
  //    std::cout << "getNO2LinesPhaseDelay(" << nc << ")=" << getNO2LinesPhaseDelay(nc).get(Angle::UnitDegree) << std::endl;
  //    std::cout << "getSO2LinesPhaseDelay(" << nc << ")=" << getSO2LinesPhaseDelay(nc).get(Angle::UnitDegree) << std::endl;
  //    std::cout << "getCOLinesPhaseDelay(" << nc << ")=" << getCOLinesPhaseDelay(nc).get(Angle::UnitDegree) << std::endl;
  return getO2LinesPhaseDelay(nc) + getO3LinesPhaseDelay(nc)
      + getN2OLinesPhaseDelay(nc) + getCOLinesPhaseDelay(nc)
      + getNO2LinesPhaseDelay(nc) + getSO2LinesPhaseDelay(nc);
}

Length RefractiveIndexProfile::getNonDispersiveDryPathLength(unsigned int nc)
{
  if(!chanIndexIsValid(nc)) {
    return Length(-999.0, Length::UnitMeter);
  }
  double wavelength = 299792458.0 / v_chanFreq_[nc]; // in m
  Length
      ll((wavelength / 360.0) * getNonDispersiveDryPhaseDelay(nc).get(Angle::UnitDegree),
         Length::UnitMeter);
  return ll;
}

Angle RefractiveIndexProfile::getNonDispersiveDryPhaseDelay(unsigned int spwid,
                                                            unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  return getNonDispersiveDryPhaseDelay(v_transfertId_[spwid] + nc);
}

Length RefractiveIndexProfile::getDispersiveDryPathLength(unsigned int nc)
{
  if(!chanIndexIsValid(nc)) {
    return Length(-999.0, Length::UnitMeter);
  }
  double wavelength = 299792458.0 / v_chanFreq_[nc]; // in m
  Length ll((wavelength / 360.0) * getDispersiveDryPhaseDelay(nc).get(Angle::UnitDegree),
            Length::UnitMeter);
  return ll;
}

Angle RefractiveIndexProfile::getDispersiveDryPhaseDelay(unsigned int spwid,
                                                         unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  return getDispersiveDryPhaseDelay(v_transfertId_[spwid] + nc);
}

Angle RefractiveIndexProfile::getAverageNonDispersiveDryPhaseDelay(unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  double av = 0.0;
  for(unsigned int i = 0; i < getNumChan(spwid); i++) {
    av = av
        + getNonDispersiveDryPhaseDelay(v_transfertId_[spwid] + i).get(Angle::UnitDegree);
  }
  av = av / getNumChan(spwid);
  Angle average(av, Angle::UnitDegree);
  return average;
}

Angle RefractiveIndexProfile::getAverageDispersiveDryPhaseDelay(unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  double av = 0.0;
  for(unsigned int i = 0; i < getNumChan(spwid); i++) {
    av = av + getDispersiveDryPhaseDelay(v_transfertId_[spwid] + i).get(Angle::UnitDegree);
  }
  av = av / getNumChan(spwid);
  Angle average(av, Angle::UnitDegree);
  return average;
}

Length RefractiveIndexProfile::getNonDispersiveDryPathLength(unsigned int spwid,
                                                             unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) {
    return Length(-999.0, Length::UnitMeter);
  }
  return getNonDispersiveDryPathLength(v_transfertId_[spwid] + nc);
}

Length RefractiveIndexProfile::getDispersiveDryPathLength(unsigned int spwid,
                                                          unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) {
    return Length(-999.0, Length::UnitMeter);
  }
  return getDispersiveDryPathLength(v_transfertId_[spwid] + nc);
}

Length RefractiveIndexProfile::getAverageNonDispersiveDryPathLength(unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) {
    return Length(-999.0, Length::UnitMeter);
  }
  double av = 0.0;
  for(unsigned int i = 0; i < getNumChan(spwid); i++) {
    av = av
        + getNonDispersiveDryPathLength(v_transfertId_[spwid] + i).get(Length::UnitMilliMeter);
  }
  av = av / getNumChan(spwid);
  Length average(av, Length::UnitMilliMeter);
  return average;
}

Length RefractiveIndexProfile::getAverageDispersiveDryPathLength(unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) {
    return Length(-999.0, Length::UnitMeter);
  }
  double av = 0.0;
  for(unsigned int i = 0; i < getNumChan(spwid); i++) {
    av = av + getDispersiveDryPathLength(v_transfertId_[spwid] + i).get(Length::UnitMilliMeter);
  }

  av = av / getNumChan(spwid);
  Length average(av, Length::UnitMilliMeter);
  return average;
}

Angle RefractiveIndexProfile::getO2LinesPhaseDelay()
{
  return getO2LinesPhaseDelay(0);
}

Length RefractiveIndexProfile::getO2LinesPathLength()
{
  return getO2LinesPathLength(0);
}

Angle RefractiveIndexProfile::getO2LinesPhaseDelay(unsigned int nc)
{
  if(!chanIndexIsValid(nc)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  double kv = 0;
  for(unsigned int j = 0; j < numLayer_; j++) {
    kv = kv + real(vv_N_O2LinesPtr_[nc]->at(j)) * v_layerThickness_[j];
  }
  Angle aa(kv * 57.29578, Angle::UnitDegree);
  return aa;
}

Length RefractiveIndexProfile::getO2LinesPathLength(unsigned int nc)
{
  if(!chanIndexIsValid(nc)) {
    return Length(-999.0, Length::UnitMeter);
  }
  double wavelength = 299792458.0 / v_chanFreq_[nc]; // in m
  Length ll((wavelength / 360.0) * getO2LinesPhaseDelay(nc).get(Angle::UnitDegree), Length::UnitMeter);
  return ll;
}

Angle RefractiveIndexProfile::getO2LinesPhaseDelay(unsigned int spwid,
                                                   unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  return getO2LinesPhaseDelay(v_transfertId_[spwid] + nc);
}

Angle RefractiveIndexProfile::getAverageO2LinesPhaseDelay(unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  double av = 0.0;
  for(unsigned int i = 0; i < getNumChan(spwid); i++) {
    av = av + getO2LinesPhaseDelay(v_transfertId_[spwid] + i).get(Angle::UnitDegree);
  }
  av = av / getNumChan(spwid);
  Angle average(av, Angle::UnitDegree);
  return average;
}

Length RefractiveIndexProfile::getO2LinesPathLength(unsigned int spwid,
                                                    unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) {
    return Length(-999.0, Length::UnitMeter);
  }
  return getO2LinesPathLength(v_transfertId_[spwid] + nc);
}

Length RefractiveIndexProfile::getAverageO2LinesPathLength(unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) {
    return Length(-999.0, Length::UnitMeter);
  }
  double av = 0.0;
  for(unsigned int i = 0; i < getNumChan(spwid); i++) {
    av = av + getO2LinesPathLength(v_transfertId_[spwid] + i).get(Length::UnitMilliMeter);
  }
  av = av / getNumChan(spwid);
  Length average(av, Length::UnitMilliMeter);
  return average;
}

Angle RefractiveIndexProfile::getO3LinesPhaseDelay()
{
  return getO3LinesPhaseDelay(0);
}

Length RefractiveIndexProfile::getO3LinesPathLength()
{
  return getO3LinesPathLength(0);
}

Angle RefractiveIndexProfile::getO3LinesPhaseDelay(unsigned int nc)
{
  if(!chanIndexIsValid(nc)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  double kv = 0;

  //    if(nc=66){cout << "vv_N_O3LinesPtr_" << ".size()=" << vv_N_O3LinesPtr_.size() << std::endl;}

  for(unsigned int j = 0; j < numLayer_; j++) {
    /* if(nc=66){
     std::cout << "j=" << j << " vv_N_O3LinesPtr_[" << nc << "]->at(" << j << ")="  << vv_N_O3LinesPtr_[nc]->at(j) << std::endl;
     } */
    kv = kv + real(vv_N_O3LinesPtr_[nc]->at(j)) * v_layerThickness_[j];
  }
  Angle aa(kv * 57.29578, Angle::UnitDegree);
  return aa;
}

Length RefractiveIndexProfile::getO3LinesPathLength(unsigned int nc)
{
  if(!chanIndexIsValid(nc)) {
    return Length(-999.0, Length::UnitMeter);
  }
  double wavelength = 299792458.0 / v_chanFreq_[nc]; // in m
  Length ll((wavelength / 360.0) * getO3LinesPhaseDelay(nc).get(Angle::UnitDegree), Length::UnitMeter);
  return ll;
}

Angle RefractiveIndexProfile::getO3LinesPhaseDelay(unsigned int spwid,
                                                   unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  return getO3LinesPhaseDelay(v_transfertId_[spwid] + nc);
}

Angle RefractiveIndexProfile::getAverageO3LinesPhaseDelay(unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  double av = 0.0;
  for(unsigned int i = 0; i < getNumChan(spwid); i++) {
    av = av + getO3LinesPhaseDelay(v_transfertId_[spwid] + i).get(Angle::UnitDegree);
  }
  av = av / getNumChan(spwid);
  Angle average(av, Angle::UnitDegree);
  return average;
}

Length RefractiveIndexProfile::getO3LinesPathLength(unsigned int spwid,
                                                    unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) {
    return Length(-999.0, Length::UnitMeter);
  }
  return getO3LinesPathLength(v_transfertId_[spwid] + nc);
}

Length RefractiveIndexProfile::getAverageO3LinesPathLength(unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) {
    return Length(-999.0, Length::UnitMeter);
  }
  double av = 0.0;
  for(unsigned int i = 0; i < getNumChan(spwid); i++) {
    av = av + getO3LinesPathLength(v_transfertId_[spwid] + i).get(Length::UnitMilliMeter);
  }
  av = av / getNumChan(spwid);
  Length average(av, Length::UnitMilliMeter);
  return average;
}

Angle RefractiveIndexProfile::getCOLinesPhaseDelay()
{
  return getCOLinesPhaseDelay(0);
}

Length RefractiveIndexProfile::getCOLinesPathLength()
{
  return getCOLinesPathLength(0);
}

Angle RefractiveIndexProfile::getCOLinesPhaseDelay(unsigned int nc)
{
  if(!chanIndexIsValid(nc)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  double kv = 0;
  for(unsigned int j = 0; j < numLayer_; j++) {
    kv = kv + real(vv_N_COLinesPtr_[nc]->at(j)) * v_layerThickness_[j];
  }
  Angle aa(kv * 57.29578, Angle::UnitDegree);
  return aa;
}

Length RefractiveIndexProfile::getCOLinesPathLength(unsigned int nc)
{
  if(!chanIndexIsValid(nc)) {
    return Length(-999.0, Length::UnitMeter);
  }
  double wavelength = 299792458.0 / v_chanFreq_[nc]; // in m
  Length ll((wavelength / 360.0) * getCOLinesPhaseDelay(nc).get(Angle::UnitDegree), Length::UnitMeter);
  return ll;
}

Angle RefractiveIndexProfile::getCOLinesPhaseDelay(unsigned int spwid,
                                                   unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  return getCOLinesPhaseDelay(v_transfertId_[spwid] + nc);
}

Angle RefractiveIndexProfile::getAverageCOLinesPhaseDelay(unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  double av = 0.0;
  for(unsigned int i = 0; i < getNumChan(spwid); i++) {
    av = av + getCOLinesPhaseDelay(v_transfertId_[spwid] + i).get(Angle::UnitDegree);
  }
  av = av / getNumChan(spwid);
  Angle average(av, Angle::UnitDegree);
  return average;
}

Length RefractiveIndexProfile::getCOLinesPathLength(unsigned int spwid,
                                                    unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) {
    return Length(-999.0, Length::UnitMeter);
  }
  return getCOLinesPathLength(v_transfertId_[spwid] + nc);
}

Length RefractiveIndexProfile::getAverageCOLinesPathLength(unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) {
    return Length(-999.0, Length::UnitMeter);
  }
  double av = 0.0;
  for(unsigned int i = 0; i < getNumChan(spwid); i++) {
    av = av + getCOLinesPathLength(v_transfertId_[spwid] + i).get(Length::UnitMilliMeter);
  }
  av = av / getNumChan(spwid);
  Length average(av, Length::UnitMilliMeter);
  return average;
}

Angle RefractiveIndexProfile::getN2OLinesPhaseDelay()
{
  return getN2OLinesPhaseDelay(0);
}

Length RefractiveIndexProfile::getN2OLinesPathLength()
{
  return getN2OLinesPathLength(0);
}

Angle RefractiveIndexProfile::getN2OLinesPhaseDelay(unsigned int nc)
{
  if(!chanIndexIsValid(nc)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  double kv = 0;
  for(unsigned int j = 0; j < numLayer_; j++) {
    kv = kv + real(vv_N_N2OLinesPtr_[nc]->at(j)) * v_layerThickness_[j];
  }
  Angle aa(kv * 57.29578, Angle::UnitDegree);
  return aa;
}

Length RefractiveIndexProfile::getN2OLinesPathLength(unsigned int nc)
{
  if(!chanIndexIsValid(nc)) {
    return Length(-999.0, Length::UnitMeter);
  }
  double wavelength = 299792458.0 / v_chanFreq_[nc]; // in m
  Length ll((wavelength / 360.0) * getN2OLinesPhaseDelay(nc).get(Angle::UnitDegree), Length::UnitMeter);
  return ll;
}

Angle RefractiveIndexProfile::getN2OLinesPhaseDelay(unsigned int spwid,
                                                    unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  return getN2OLinesPhaseDelay(v_transfertId_[spwid] + nc);
}

Angle RefractiveIndexProfile::getAverageN2OLinesPhaseDelay(unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  double av = 0.0;
  for(unsigned int i = 0; i < getNumChan(spwid); i++) {
    av = av + getN2OLinesPhaseDelay(v_transfertId_[spwid] + i).get(Angle::UnitDegree);
  }
  av = av / getNumChan(spwid);
  Angle average(av, Angle::UnitDegree);
  return average;
}

Length RefractiveIndexProfile::getN2OLinesPathLength(unsigned int spwid,
                                                     unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) {
    return Length(-999.0, Length::UnitMeter);
  }
  return getN2OLinesPathLength(v_transfertId_[spwid] + nc);
}

Length RefractiveIndexProfile::getAverageN2OLinesPathLength(unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) {
    return Length(-999.0, Length::UnitMeter);
  }
  double av = 0.0;
  for(unsigned int i = 0; i < getNumChan(spwid); i++) {
    av = av + getN2OLinesPathLength(v_transfertId_[spwid] + i).get(Length::UnitMilliMeter);
  }
  av = av / getNumChan(spwid);
  Length average(av, Length::UnitMilliMeter);
  return average;
}






Angle RefractiveIndexProfile::getNO2LinesPhaseDelay()
{
  return getNO2LinesPhaseDelay(0);
}

Length RefractiveIndexProfile::getNO2LinesPathLength()
{
  return getNO2LinesPathLength(0);
}

Angle RefractiveIndexProfile::getNO2LinesPhaseDelay(unsigned int nc)
{
  if(!chanIndexIsValid(nc)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  double kv = 0;
  for(unsigned int j = 0; j < numLayer_; j++) {
    kv = kv + real(vv_N_NO2LinesPtr_[nc]->at(j)) * v_layerThickness_[j];
  }
  Angle aa(kv * 57.29578, Angle::UnitDegree);
  return aa;
}

Length RefractiveIndexProfile::getNO2LinesPathLength(unsigned int nc)
{
  if(!chanIndexIsValid(nc)) {
    return Length(-999.0, Length::UnitMeter);
  }
  double wavelength = 299792458.0 / v_chanFreq_[nc]; // in m
  Length ll((wavelength / 360.0) * getNO2LinesPhaseDelay(nc).get(Angle::UnitDegree), Length::UnitMeter);
  return ll;
}

Angle RefractiveIndexProfile::getNO2LinesPhaseDelay(unsigned int spwid,
                                                    unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  return getNO2LinesPhaseDelay(v_transfertId_[spwid] + nc);
}

Angle RefractiveIndexProfile::getAverageNO2LinesPhaseDelay(unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  double av = 0.0;
  for(unsigned int i = 0; i < getNumChan(spwid); i++) {
    av = av + getNO2LinesPhaseDelay(v_transfertId_[spwid] + i).get(Angle::UnitDegree);
  }
  av = av / getNumChan(spwid);
  Angle average(av, Angle::UnitDegree);
  return average;
}

Length RefractiveIndexProfile::getNO2LinesPathLength(unsigned int spwid,
                                                     unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) {
    return Length(-999.0, Length::UnitMeter);
  }
  return getNO2LinesPathLength(v_transfertId_[spwid] + nc);
}

Length RefractiveIndexProfile::getAverageNO2LinesPathLength(unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) {
    return Length(-999.0, Length::UnitMeter);
  }
  double av = 0.0;
  for(unsigned int i = 0; i < getNumChan(spwid); i++) {
    av = av + getNO2LinesPathLength(v_transfertId_[spwid] + i).get(Length::UnitMilliMeter);
  }
  av = av / getNumChan(spwid);
  Length average(av, Length::UnitMilliMeter);
  return average;
}








Angle RefractiveIndexProfile::getSO2LinesPhaseDelay()
{
  return getSO2LinesPhaseDelay(0);
}

Length RefractiveIndexProfile::getSO2LinesPathLength()
{
  return getSO2LinesPathLength(0);
}

Angle RefractiveIndexProfile::getSO2LinesPhaseDelay(unsigned int nc)
{
  if(!chanIndexIsValid(nc)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  double kv = 0;
  for(unsigned int j = 0; j < numLayer_; j++) {
    kv = kv + real(vv_N_SO2LinesPtr_[nc]->at(j)) * v_layerThickness_[j];
  }
  Angle aa(kv * 57.29578, Angle::UnitDegree);
  return aa;
}

Length RefractiveIndexProfile::getSO2LinesPathLength(unsigned int nc)
{
  if(!chanIndexIsValid(nc)) {
    return Length(-999.0, Length::UnitMeter);
  }
  double wavelength = 299792458.0 / v_chanFreq_[nc]; // in m
  Length ll((wavelength / 360.0) * getSO2LinesPhaseDelay(nc).get(Angle::UnitDegree), Length::UnitMeter);
  return ll;
}

Angle RefractiveIndexProfile::getSO2LinesPhaseDelay(unsigned int spwid,
                                                    unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  return getSO2LinesPhaseDelay(v_transfertId_[spwid] + nc);
}

Angle RefractiveIndexProfile::getAverageSO2LinesPhaseDelay(unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  double av = 0.0;
  for(unsigned int i = 0; i < getNumChan(spwid); i++) {
    av = av + getSO2LinesPhaseDelay(v_transfertId_[spwid] + i).get(Angle::UnitDegree);
  }
  av = av / getNumChan(spwid);
  Angle average(av, Angle::UnitDegree);
  return average;
}

Length RefractiveIndexProfile::getSO2LinesPathLength(unsigned int spwid,
                                                     unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) {
    return Length(-999.0, Length::UnitMeter);
  }
  return getSO2LinesPathLength(v_transfertId_[spwid] + nc);
}

Length RefractiveIndexProfile::getAverageSO2LinesPathLength(unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) {
    return Length(-999.0, Length::UnitMeter);
  }
  double av = 0.0;
  for(unsigned int i = 0; i < getNumChan(spwid); i++) {
    av = av + getSO2LinesPathLength(v_transfertId_[spwid] + i).get(Length::UnitMilliMeter);
  }
  av = av / getNumChan(spwid);
  Length average(av, Length::UnitMilliMeter);
  return average;
}


Angle RefractiveIndexProfile::getNonDispersiveH2OPhaseDelay(const Length &integratedwatercolumn)
{
  return getNonDispersiveH2OPhaseDelay(integratedwatercolumn,0);
}

Length RefractiveIndexProfile::getNonDispersiveH2OPathLength(const Length &integratedwatercolumn)
{
  return getNonDispersiveH2OPathLength(integratedwatercolumn,0);
}

Angle RefractiveIndexProfile::getNonDispersiveH2OPhaseDelay(const Length &integratedwatercolumn,
                                                            unsigned int nc)
{
  double kv = 0;
  if(!chanIndexIsValid(nc)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  for(unsigned int j = 0; j < numLayer_; j++) {
    kv = kv + real(vv_N_H2OContPtr_[nc]->at(j)) * v_layerThickness_[j];
  }
  Angle aa(kv*(integratedwatercolumn.get()/getGroundWH2O().get())* 57.29578, Angle::UnitDegree);
  return aa;
}

Length RefractiveIndexProfile::getNonDispersiveH2OPathLength(const Length &integratedwatercolumn,
                                                             unsigned int nc)
{
  if(!chanIndexIsValid(nc)) {
    return Length(-999.0, Length::UnitMeter);
  }
  double wavelength = 299792458.0 / v_chanFreq_[nc]; // in m
  Length ll((wavelength / 360.0) * getNonDispersiveH2OPhaseDelay(integratedwatercolumn,nc).get(Angle::UnitDegree),Length::UnitMeter);
  return ll;
}

Angle RefractiveIndexProfile::getNonDispersiveH2OPhaseDelay(const Length &integratedwatercolumn,
                                                            unsigned int spwid,
                                                            unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  return getNonDispersiveH2OPhaseDelay(integratedwatercolumn,v_transfertId_[spwid] + nc);
}

Angle RefractiveIndexProfile::getAverageNonDispersiveH2OPhaseDelay(const Length &integratedwatercolumn,
                                                                   unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) {
    return Angle(-999.0, Angle::UnitDegree);
  }
  double av = 0.0;
  for(unsigned int i = 0; i < getNumChan(spwid); i++) {
    av = av
        + getNonDispersiveH2OPhaseDelay(v_transfertId_[spwid] + i).get(Angle::UnitDegree);
  }
  av = av / getNumChan(spwid);
  Angle average(av*(integratedwatercolumn.get()/getGroundWH2O().get()), Angle::UnitDegree);
  return average;
}

Length RefractiveIndexProfile::getNonDispersiveH2OPathLength(const Length &integratedwatercolumn,
                                                             unsigned int spwid,
                                                             unsigned int nc)
{
  if(!spwidAndIndexAreValid(spwid, nc)) {
    return Length(-999.0);
  }
  return getNonDispersiveH2OPathLength(integratedwatercolumn,v_transfertId_[spwid] + nc);
}

Length RefractiveIndexProfile::getAverageNonDispersiveH2OPathLength(const Length &integratedwatercolumn,
                                                                    unsigned int spwid)
{
  if(!spwidAndIndexAreValid(spwid, 0)) {
    return Length(-999.0);
  }
  double av = 0.0;
  for(unsigned int i = 0; i < getNumChan(spwid); i++) {
    av = av
      + getNonDispersiveH2OPathLength(integratedwatercolumn,v_transfertId_[spwid] + i).get(Length::UnitMeter);
  }
  av = av / getNumChan(spwid);
  Length average(av, Length::UnitMeter);
  return average;
}

// NB: the function chanIndexIsValid will be overrided by ....
bool RefractiveIndexProfile::chanIndexIsValid(unsigned int nc)
{
  if(nc < vv_N_H2OLinesPtr_.size()) return true;
  if(nc < v_chanFreq_.size()) {
    std::cout
        << " RefractiveIndexProfile: Requested index in a new spectral window ==> update profile"
        << std::endl;
    mkRefractiveIndexProfile();
    // std::cout << "...and we return" << std::endl;
    return true;
  }
  std::cout << " RefractiveIndexProfile: ERROR: Invalid channel frequency index"
      << std::endl;
  return false;
}

// NB: the function spwidAndIndexAreValid will be overrided by ...
bool RefractiveIndexProfile::spwidAndIndexAreValid(unsigned int spwid,
                                                   unsigned int idx)
{

  if(spwid > getNumSpectralWindow() - 1) {
    std::cout
        << " RefractiveIndexProfile: ERROR: spectral window identifier out of range "
        << std::endl;
    return false;
  }
  if(idx > getNumChan(spwid) - 1) {
    std::cout << " RefractiveIndexProfile: ERROR: channel index out of range "
        << std::endl;
    return false;
  }
  unsigned int nc = v_transfertId_[spwid] + idx;
  bool valid = chanIndexIsValid(nc);
  return valid;
}

void RefractiveIndexProfile::updateNewSpectralWindows()
{
  mkRefractiveIndexProfile();
}

ATM_NAMESPACE_END
