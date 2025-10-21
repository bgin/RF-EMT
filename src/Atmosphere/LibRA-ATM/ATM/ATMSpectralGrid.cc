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
 * "@(#) $Id: ATMSpectralGrid.cpp Exp $"
 *
 * who       when      what
 * --------  --------  ----------------------------------------------
 * pardo     24/03/09  created
 */

#include "ATMSpectralGrid.h"

#include <iostream>
#include <limits>
#include <math.h>
#include <vector>
#include <string>



ATM_NAMESPACE_BEGIN

// public methods:
// constructors
SpectralGrid::SpectralGrid(const Frequency &oneFreq)
{
  v_chanFreq_.reserve(1);
  freqUnits_ = Frequency::UnitHertz;
  v_transfertId_.resize(0); // not sure this is necessary!
  unsigned int numChan = 1;
  unsigned int refChan = 0;
  Frequency chanSep(0.0);
  add(numChan, refChan, oneFreq, chanSep);
  std::vector<unsigned int> v_dummyInt;
  vv_assocSpwId_.push_back(v_dummyInt); // put an empty std::vector
  std::vector<std::string> v_dummyString;
  vv_assocNature_.push_back(v_dummyString); // put an empty std::vector
}

SpectralGrid::SpectralGrid(unsigned int numChan,
                           unsigned int refChan,
                           const Frequency &refFreq,
                           const Frequency &chanSep)
{
  //  cout << " SpectralGrid constructor" << endl;
  v_chanFreq_.reserve(numChan);
  freqUnits_ = Frequency::UnitHertz;
  v_transfertId_.resize(0); // not sure this is necessary!
  add(numChan, refChan, refFreq, chanSep);
  std::vector<unsigned int> v_dummyInt;
  vv_assocSpwId_.push_back(v_dummyInt); // put an empty std::vector
  std::vector<std::string> v_dummyString;
  vv_assocNature_.push_back(v_dummyString); // put an empty std::vector
}

SpectralGrid::SpectralGrid(unsigned int numChan,
                           unsigned int refChan,
                           const Frequency &refFreq,
                           const Frequency &chanSep,
                           const Frequency &intermediateFreq,
                           const SidebandSide &sbSide,
                           const SidebandType &sbType)
{
  freqUnits_ = Frequency::UnitHertz;
  v_transfertId_.resize(0); // not sure this is necessary!
  v_chanFreq_.reserve(numChan);
  /* cout << " la" << endl; */
  add(numChan, refChan, refFreq, chanSep, intermediateFreq, sbSide, sbType);
}

SpectralGrid::SpectralGrid(unsigned int numChan,
                           unsigned int refChan,
                           double* chanFreq,
                           Frequency::Units freqUnits)
{
  v_chanFreq_.reserve(numChan);
  freqUnits_ = Frequency::UnitHertz;
  v_transfertId_.resize(0); // not sure this is necessary!
  add(numChan, refChan, chanFreq, freqUnits);
}

void SpectralGrid::add(unsigned int numChan,
                       unsigned int refChan,
                       const Frequency &refFreq,
                       const Frequency &chanSep,
                       const Frequency &intermediateFreq,
                       const SidebandSide &sbSide,
                       const SidebandType &sbType)
{


  double chSep;
  std::vector<std::string> v_assocNature;
  std::vector<unsigned int> v_assocSpwId;

  unsigned int spwId = v_transfertId_.size();

  if(sbSide == LSB) { // LSB tuning
    // the LSB:
    add(numChan, refChan, refFreq, chanSep); // LSB
    v_sidebandSide_[spwId] = LSB;
    v_sidebandType_[spwId] = sbType;
    double loFreq    = refFreq.get() + intermediateFreq.get(); // store loFreq for USB
    v_loFreq_[spwId] = loFreq;
    v_assocSpwId.push_back(v_numChan_.size());
    vv_assocSpwId_[vv_assocSpwId_.size() - 1] = v_assocSpwId;
    v_assocNature.push_back("USB");
    vv_assocNature_[vv_assocNature_.size() - 1] = v_assocNature;

    // the USB:
    spwId = v_transfertId_.size();
    double refFreqUSB = refFreq.get() + 2.*intermediateFreq.get();  // fix refFreq in the image band (refChan is unchanged)
    chSep = -chanSep.get();
    add(numChan, refChan, Frequency(refFreqUSB), Frequency(chSep));

    v_sidebandSide_[spwId] = USB;
    v_sidebandType_[spwId] = sbType;
    v_loFreq_[spwId] = loFreq;

    v_assocSpwId[0] = v_numChan_.size() - 2;
    vv_assocSpwId_[vv_assocSpwId_.size() - 1] = v_assocSpwId;
    v_assocNature[0] = "LSB";
    vv_assocNature_[vv_assocNature_.size() - 1] = v_assocNature;

  } else { // USB tuning
    // the USB:
    add(numChan, refChan, refFreq, chanSep);

    v_sidebandSide_[spwId] = USB;
    v_sidebandType_[spwId] = sbType;
    double loFreq    = refFreq.get() - intermediateFreq.get();
    v_loFreq_[spwId] = loFreq;

    v_assocSpwId.push_back(v_numChan_.size());
    vv_assocSpwId_[vv_assocSpwId_.size() - 1] = v_assocSpwId;
    v_assocNature.push_back("LSB");
    vv_assocNature_[vv_assocNature_.size() - 1] = v_assocNature;

    // the LSB:
    spwId = v_transfertId_.size();
    double refFreqLSB = refFreq.get() - 2.*intermediateFreq.get();  // fix refFreq in the image band (refChan is unchanged)
    chSep = -chanSep.get();
    add(numChan, refChan, Frequency(refFreqLSB), Frequency(chSep));

    v_sidebandSide_[spwId] = LSB;
    v_sidebandType_[spwId] = sbType;
    v_loFreq_[spwId] = loFreq;

    v_assocSpwId[0] = v_numChan_.size() - 2;
    vv_assocSpwId_[vv_assocSpwId_.size() - 1] = v_assocSpwId;
    v_assocNature[0] = "USB";
    vv_assocNature_[vv_assocNature_.size() - 1] = v_assocNature;
  }
}

unsigned int SpectralGrid::add(unsigned int numChan,
                               unsigned int refChan,
                               const Frequency &refFreq,
                               const Frequency &chanSep)
{
  freqUnits_ = Frequency::UnitHertz;

  unsigned int spwId = v_transfertId_.size();
  v_loFreq_.push_back(refFreq.get());

  if(spwId == 0) {
    v_transfertId_.push_back(0);
  } else {
    v_transfertId_.push_back(v_transfertId_[spwId - 1] + v_numChan_[spwId - 1]);
  }

  v_numChan_.push_back(numChan);
  v_refChan_.push_back(refChan);

  v_refFreq_.push_back(refFreq.get(freqUnits_));
  v_chanSep_.push_back(chanSep.get(freqUnits_));

  double* chanFreq = new double[numChan];
  {
	// CAS-9762 refChan should not be 1-base
	// ICT-9490	Change refChan convention so that it is 0-based
    // double freqOffset = v_refFreq_[spwId] - v_chanSep_[spwId]
    //     * (double) (v_refChan_[spwId] - 1.);
    double freqOffset = v_refFreq_[spwId] - v_chanSep_[spwId]
        * (double) (v_refChan_[spwId]);

    for(unsigned int i = 0; i < numChan; i++) {
      chanFreq[i] = freqOffset + (double) i * v_chanSep_[spwId];
    }
    appendChanFreq(numChan, chanFreq);

    if(numChan > 1) {
      if(chanFreq[0] > chanFreq[1]) {
        v_minFreq_.push_back(chanFreq[numChan - 1]);
        v_maxFreq_.push_back(chanFreq[0]);
      } else {
        v_minFreq_.push_back(chanFreq[0]);
        v_maxFreq_.push_back(chanFreq[numChan - 1]);
      }
    } else {
      v_minFreq_.push_back(chanFreq[0]);
      v_maxFreq_.push_back(chanFreq[0]);
    }

    v_sidebandSide_ .push_back(NOSB);
    v_sidebandType_ .push_back(NOTYPE);
    v_intermediateFrequency_.push_back(0.0);
  }
  std::vector<unsigned int> v_dummyAssoc;
  vv_assocSpwId_.push_back(v_dummyAssoc);
  std::vector<std::string> v_dummyNature;
  vv_assocNature_.push_back(v_dummyNature);

  delete [] chanFreq;
  return spwId;
}

void SpectralGrid::appendChanFreq(unsigned int numChan, double* chanFreq)
{

  //    unsigned int k=v_chanFreq_.size();
  for(unsigned int i = 0; i < numChan; i++) {
    v_chanFreq_.push_back(chanFreq[i]); // cout << i << "v_chanFreq_="<<v_chanFreq_[k+i]<<endl;
  }
}

void SpectralGrid::appendChanFreq(unsigned int numChan, const std::vector<double> &chanFreq)
{

  //    unsigned int k=v_chanFreq_.size();
  for(unsigned int i = 0; i < numChan; i++) {
    v_chanFreq_.push_back(chanFreq[i]); // cout << i << "v_chanFreq_="<<v_chanFreq_[k+i]<<endl;
  }
}

unsigned int SpectralGrid::add(unsigned int numChan,
                               unsigned int refChan,
                               double* chanFreq,
                               Frequency::Units freqUnits)
{
  double fact = 1.0;
  if(freqUnits == Frequency::UnitGigaHertz) fact = 1.0E9;
  if(freqUnits == Frequency::UnitMegaHertz) fact = 1.0E6;
  if(freqUnits == Frequency::UnitKiloHertz) fact = 1.0E3;

  unsigned int spwId = v_transfertId_.size();
  if(spwId == 0) {
    v_transfertId_.push_back(0);
  } else {
    v_transfertId_.push_back(v_transfertId_[spwId - 1] + v_numChan_[spwId - 1]);
  }

  v_numChan_.push_back(numChan);
  v_refChan_.push_back(refChan);

  bool regular = true;
  double minFreq = 1.E30;
  double maxFreq = 0;
  double chanSep = 0;
  if(numChan > 1) chanSep = fact * (chanFreq[1] - chanFreq[0]);

  chanFreq[0] = fact * chanFreq[0];
  for(unsigned int i = 1; i < numChan; i++) {
    chanFreq[i] = fact * chanFreq[i];
    if(fabs(chanFreq[i] - chanFreq[i - 1] - chanSep) > 1.0E-12) regular = false;
    if(chanFreq[i] < minFreq) minFreq = chanFreq[i];
    if(chanFreq[i] > maxFreq) maxFreq = chanFreq[i];
  }
  appendChanFreq(numChan, chanFreq);

  v_refFreq_.push_back(chanFreq[refChan - 1]);
  if(regular) {
    v_chanSep_.push_back(chanSep);
  } else {
    v_chanSep_.push_back(0);
  }
  v_sidebandSide_ .push_back(NOSB);
  v_sidebandType_ .push_back(NOTYPE);
  v_intermediateFrequency_.push_back(0.0);

  return spwId;
}

SpectralGrid::SpectralGrid(unsigned int numChan,
                           double refFreq,
                           double* chanFreq,
                           Frequency::Units freqUnits)
{
  v_chanFreq_.reserve(numChan);
  freqUnits_ = Frequency::UnitHertz;
  v_transfertId_.resize(0); // not sure this is necessary!
  add(numChan, refFreq, chanFreq, freqUnits);
  std::vector<unsigned int> v_dummyInt;
  vv_assocSpwId_.push_back(v_dummyInt); // put an empty std::vector
  std::vector<std::string> v_dummyString;
  vv_assocNature_.push_back(v_dummyString); // put an empty std::vector
}

SpectralGrid::SpectralGrid(double refFreq,
                           const std::vector<double> &chanFreq,
                           Frequency::Units freqUnits)
{
  v_chanFreq_.reserve(chanFreq.size());
  freqUnits_ = Frequency::UnitHertz;
  v_transfertId_.resize(0); // not sure this is necessary!
  add(chanFreq.size(), refFreq, chanFreq, freqUnits);
  std::vector<unsigned int> v_dummyInt;
  vv_assocSpwId_.push_back(v_dummyInt); // put an empty std::vector
  std::vector<std::string> v_dummyString;
  vv_assocNature_.push_back(v_dummyString); // put an empty std::vector
}

SpectralGrid::SpectralGrid(const std::vector<double> &chanFreq, Frequency::Units freqUnits)
{
  v_chanFreq_.reserve(chanFreq.size());
  freqUnits_ = Frequency::UnitHertz;
  v_transfertId_.resize(0); // not sure this is necessary!
  double refFreq = (Frequency(chanFreq[0], freqUnits)).get(Frequency::UnitHertz); // We take the frequency of the first channel as
  // reference frequency because it has not been specified
  add(chanFreq.size(), refFreq, chanFreq, freqUnits);
  std::vector<unsigned int> v_dummyInt;
  vv_assocSpwId_.push_back(v_dummyInt); // put an empty std::vector
  std::vector<std::string> v_dummyString;
  vv_assocNature_.push_back(v_dummyString); // put an empty std::vector
}

SpectralGrid::SpectralGrid(const std::vector<Frequency> &chanFreq)
{
  v_chanFreq_.reserve(chanFreq.size());
  freqUnits_ = Frequency::UnitHertz;
  v_transfertId_.resize(0); // not sure this is necessary!
  double refFreq = chanFreq[0].get(Frequency::UnitHertz); // We take the frequency of the first channel as
  // reference frequency because it has not been specified
  std::vector<double> chanFreq_double;
  for(unsigned int i = 0; i < chanFreq.size(); i++) {
    chanFreq_double.push_back(chanFreq[i].get(Frequency::UnitGigaHertz));
  }

  add(chanFreq.size(), refFreq, chanFreq_double, Frequency::UnitGigaHertz);
  std::vector<unsigned int> v_dummyInt;
  vv_assocSpwId_.push_back(v_dummyInt); // put an empty std::vector
  std::vector<std::string> v_dummyString;
  vv_assocNature_.push_back(v_dummyString); // put an empty std::vector
}

unsigned int SpectralGrid::add(unsigned int numChan,
                               double refFreq,
                               double* chanFreq,
                               Frequency::Units freqUnits)
{

  bool regular = true;
  double fact = 1.0;
  if(freqUnits == Frequency::UnitGigaHertz) fact = 1.0E9;
  if(freqUnits == Frequency::UnitMegaHertz) fact = 1.0E6;
  if(freqUnits == Frequency::UnitKiloHertz) fact = 1.0E3;

  freqUnits_ = Frequency::UnitHertz;

  unsigned int spwId = v_transfertId_.size();
  if(spwId == 0) {
    v_transfertId_.push_back(0);
  } else {
    v_transfertId_.push_back(v_transfertId_[spwId - 1] + v_numChan_[spwId - 1]);
  }

  v_numChan_.push_back(numChan);
  v_refFreq_.push_back(fact * refFreq);

  double chanSep = fact * (chanFreq[1] - chanFreq[0]);
  double minFreq = 1.E30;
  double maxFreq = 0;

  chanFreq[0] = fact * chanFreq[0];
  for(unsigned int i = 1; i < numChan; i++) {
    chanFreq[i] = fact * chanFreq[i];
    if(fabs(chanFreq[i] - chanFreq[i - 1] - chanSep) > 1.0E-12) regular = false;
    if(chanFreq[i] < minFreq) minFreq = chanFreq[i];
    if(chanFreq[i] > maxFreq) maxFreq = chanFreq[i];
  }
  appendChanFreq(numChan, chanFreq);
  v_minFreq_.push_back(minFreq);
  v_maxFreq_.push_back(maxFreq);

  if(numChan > 1) {
    if(regular) {
      v_refChan_.push_back((unsigned int) (1. + (refFreq - v_chanFreq_[0]
          + 1.E-12) / chanSep));
      v_chanSep_.push_back(chanSep);
    } else {
      v_refChan_.push_back(0);
      v_chanSep_.push_back(0.0);
    }
  } else {
    v_refChan_.push_back(0);
    v_chanSep_.push_back(0.0);
  }
  v_sidebandSide_ .push_back(NOSB);
  v_sidebandType_ .push_back(NOTYPE);
  v_intermediateFrequency_.push_back(0.0);

  return spwId;
}

  unsigned int SpectralGrid::add(unsigned int numChan,
                                 double refFreq,
                                 const std::vector<double> &chanFreq,
                                 Frequency::Units freqUnits)
{

  bool regular = true;
  double fact = 1.0;
  if(freqUnits == Frequency::UnitGigaHertz) fact = 1.0E9;
  if(freqUnits == Frequency::UnitMegaHertz) fact = 1.0E6;
  if(freqUnits == Frequency::UnitKiloHertz) fact = 1.0E3;

  freqUnits_ =Frequency::UnitHertz;

  unsigned int spwId = v_transfertId_.size();
  if(spwId == 0) {
    v_transfertId_.push_back(0);
  } else {
    v_transfertId_.push_back(v_transfertId_[spwId - 1] + v_numChan_[spwId - 1]);
  }

  v_numChan_.push_back(numChan);
  v_refFreq_.push_back(fact * refFreq);

  double chanSep = fact * (chanFreq[1] - chanFreq[0]);

  std::vector<double> chanFreqHz(numChan);
  chanFreqHz[0] = fact * chanFreq[0];
  double minFreq = chanFreqHz[0];
  double maxFreq = chanFreqHz[0];

  for(unsigned int i = 1; i < numChan; i++) {
    chanFreqHz[i] = fact * chanFreq[i];
    if(fabs(chanFreqHz[i] - chanFreqHz[i - 1] - chanSep) > 1.0E-12) regular = false;
    if(chanFreqHz[i] < minFreq) minFreq = chanFreqHz[i];
    if(chanFreqHz[i] > maxFreq) maxFreq = chanFreqHz[i];
  }
  appendChanFreq(numChan, chanFreqHz);
  v_minFreq_.push_back(minFreq);
  v_maxFreq_.push_back(maxFreq);

  if(numChan > 1) {
    if(regular) {
      v_refChan_.push_back((unsigned int) (1. + (refFreq - v_chanFreq_[0]
          + 1.E-12) / chanSep));
      v_chanSep_.push_back(chanSep);
    } else {
      v_refChan_.push_back(0);
      v_chanSep_.push_back(0.0);
    }
  } else {
    v_refChan_.push_back(0);
    v_chanSep_.push_back(0.0);
  }
  v_sidebandSide_ .push_back(NOSB);
  v_sidebandType_ .push_back(NOTYPE);
  v_intermediateFrequency_.push_back(0.0);

  return spwId;
}

SpectralGrid::SpectralGrid(const SpectralGrid & a)
{
  // cout <<  " SpectralGrid copy constructor" << endl;
  freqUnits_ = a.freqUnits_;
  v_chanFreq_ = a.v_chanFreq_;

  v_numChan_ = a.v_numChan_; // cout << numChan_ << endl;
  v_refChan_ = a.v_refChan_; // cout << refChan_ << endl;
  v_refFreq_ = a.v_refFreq_; // cout << refChan_ << endl;
  v_chanSep_ = a.v_chanSep_; // cout << chanSep_ << endl;
  v_maxFreq_ = a.v_maxFreq_; // cout << maxFreq_ << endl;
  v_minFreq_ = a.v_minFreq_; // cout << minFreq_ << endl;
  v_intermediateFrequency_ = a.v_intermediateFrequency_;
  v_loFreq_ = a.v_loFreq_;
  v_sidebandSide_ = a.v_sidebandSide_;
  v_sidebandType_ = a.v_sidebandType_;
  vv_assocSpwId_ = a.vv_assocSpwId_;
  vv_assocNature_ = a.vv_assocNature_;
  v_transfertId_ = a.v_transfertId_;

  // cout << "v_chanFreq_.size()=" << v_chanFreq_.size() << endl;
}

SpectralGrid::SpectralGrid()
{
}

// destructor
SpectralGrid::~SpectralGrid()
{
}

bool SpectralGrid::wrongSpwId(unsigned int spwId) const
{
  if(spwId > (v_transfertId_.size() - 1)) {
    std::cout << " SpectralGrid: ERROR: " << spwId
        << " is a wrong spectral window identifier" << std::endl;
    return (bool) true;
  }
  return (bool) false;
}

// accessors and utilities:
unsigned int SpectralGrid::getNumSpectralWindow() const
{
  return v_transfertId_.size();
}
unsigned int SpectralGrid::getNumChan() const
{
  return v_numChan_[0];
}
unsigned int SpectralGrid::getNumChan(unsigned int spwId) const
{
  if(wrongSpwId(spwId)) return 0;
  return v_numChan_[spwId];
}

unsigned int SpectralGrid::getRefChan() const
{
  return v_refChan_[0];
}
unsigned int SpectralGrid::getRefChan(unsigned int spwId) const
{
  if(wrongSpwId(spwId)) return 32767;
  return v_refChan_[spwId];
}

Frequency SpectralGrid::getRefFreq() const
{
  return Frequency(v_refFreq_[0], Frequency::UnitHertz);
}
Frequency SpectralGrid::getRefFreq(unsigned int spwId) const
{
  if(wrongSpwId(spwId)) return 32767.;
  return Frequency(v_refFreq_[spwId], Frequency::UnitHertz);
}

Frequency SpectralGrid::getChanSep() const
{
  return Frequency(v_chanSep_[0], Frequency::UnitHertz);
}
Frequency SpectralGrid::getChanSep(unsigned int spwId) const
{
  if(wrongSpwId(spwId)) return 32767.;
  return Frequency(v_chanSep_[spwId], Frequency::UnitHertz);
}

Frequency SpectralGrid::getChanFreq(unsigned int i) const
{
  return Frequency(v_chanFreq_[i], Frequency::UnitHertz);
}

Frequency SpectralGrid::getChanWidth(unsigned int i) const
{
 if(i == 0){
   return getChanFreq(i+1)-getChanFreq(i);
 }else{
   return getChanFreq(i)-getChanFreq(i-1);
 }
}

Frequency SpectralGrid::getChanFreq(unsigned int spwId, unsigned int chanIdx) const
{
  if(wrongSpwId(spwId)) return 32767.;
  return Frequency(v_chanFreq_[v_transfertId_[spwId] + chanIdx], Frequency::UnitHertz);
}

Frequency SpectralGrid::getChanWidth(unsigned int spwId, unsigned int chanIdx) const
{
  if(wrongSpwId(spwId)) return 32767.;

  unsigned int banda=spwId;
  unsigned int canalmasuno=chanIdx+1;
  unsigned int canal=chanIdx;
  unsigned int canalmenosuno=chanIdx-1;
  //  cout << "banda,canal-1,canal,canal+1= " << banda << " " << canalmenosuno << " " << canal << " " << canalmasuno << endl;

  if(chanIdx == 0){
    return getChanFreq(spwId,canalmasuno)-getChanFreq(spwId,canal);
  }else{
    //  cout << "ChanFreq(" << banda << "," << canal << ")=" << getChanFreq(banda,canal).get(Frequency::UnitGigaHertz) << endl;
    //  cout << "ChanFreq(" << banda << "," << canalmenosuno << ")=" << getChanFreq(banda,canalmenosuno).get(Frequency::UnitGigaHertz) << endl;
    return getChanFreq(banda,canal)-getChanFreq(banda,canalmenosuno);
  }
}

std::vector<double> SpectralGrid::getSbChanFreq(unsigned int spwId,
                                           unsigned int chanIdx,
                                           const std::string &units) const
{
  std::vector<double> v_dummyVector;
  if(wrongSpwId(spwId)) return v_dummyVector;
  v_dummyVector.push_back(getChanFreq(spwId, chanIdx).get(units));
  for(unsigned int n = 0; n < vv_assocNature_[spwId].size(); n++) {
    if(vv_assocNature_[spwId][n] == "USB" || vv_assocNature_[spwId][n] == "LSB") {
      unsigned int assocSpwId = vv_assocSpwId_[spwId][n];
      v_dummyVector.push_back(getChanFreq(assocSpwId, chanIdx).get(units));
    }
  }
  return v_dummyVector;
}

std::vector<double> SpectralGrid::getSpectralWindow(unsigned int spwId) const
{
  std::vector<double> v_chanFreq;
  if(wrongSpwId(spwId)) return v_chanFreq;
  v_chanFreq.reserve(v_numChan_[spwId]);
  for(unsigned int n = 0; n < v_numChan_[spwId]; n++)
    v_chanFreq.push_back(v_chanFreq_[v_transfertId_[spwId] + n]);
  return v_chanFreq;
}

Frequency SpectralGrid::getMinFreq() const
{
  return Frequency(v_minFreq_[0], Frequency::UnitHertz);
}
Frequency SpectralGrid::getMinFreq(unsigned int spwId) const
{
  if(wrongSpwId(spwId)) return 32767.;
  return Frequency(v_minFreq_[spwId], Frequency::UnitHertz);
}

Frequency SpectralGrid::getMaxFreq() const
{
  return Frequency(v_maxFreq_[0], Frequency::UnitHertz);
}
Frequency SpectralGrid::getMaxFreq(unsigned int spwId) const
{
  if(wrongSpwId(spwId)) return 32767.;
  return Frequency(v_maxFreq_[spwId], Frequency::UnitHertz);
}

double SpectralGrid::getChanNum(double freq) const
{
  if(v_numChan_[0] == 1) return 1;
  if(v_chanSep_[0] == 0.0) { // irregular grid, look for the nearest channel
    double sep = 1.E30;
    int k = -1;
    for(unsigned int i = 0; i < v_numChan_[0]; i++) {
      if(sep > fabs(v_chanFreq_[v_transfertId_[0] + i] - freq)) {
        sep = fabs(v_chanFreq_[v_transfertId_[0] + i] - freq);
        k = i;
      }
    }
    return (double) 1 - v_refChan_[0] + k; // channel the nearest
  } else { // regular spectral grid
    return (freq - v_refFreq_[0]) / v_chanSep_[0];
  }
}
double SpectralGrid::getChanNum(unsigned int spwId, double freq) const
{
  if(wrongSpwId(spwId)) return 32767.;
  if(v_numChan_[spwId] == 1) return 1;
  if(v_chanSep_[spwId] == 0.0) { // irregular grid, look for the nearest channel
    double sep = 1.E30;
    int k = -1;
    for(unsigned int i = 0; i < v_numChan_[spwId]; i++) {
      if(sep > fabs(v_chanFreq_[v_transfertId_[spwId] + i] - freq)) {
        sep = fabs(v_chanFreq_[v_transfertId_[spwId] + i] - freq);
        k = i;
      }
    }
    return (double) 1 - v_refChan_[spwId] + k; // channel the nearest
  } else { // regular spectral grid
    return (freq - v_refFreq_[spwId]) / v_chanSep_[spwId];
  }
}

Frequency SpectralGrid::getBandwidth() const
{
  return Frequency(v_maxFreq_[0] - v_minFreq_[0], Frequency::UnitHertz);
}

Frequency SpectralGrid::getBandwidth(unsigned int spwId) const
{
  if(wrongSpwId(spwId)) return 32767.;
  return Frequency(v_maxFreq_[spwId] - v_minFreq_[spwId], Frequency::UnitHertz);
}

bool SpectralGrid::isRegular() const
{
  if(v_chanSep_[0] == 0.0) return false;
  return true;
}

bool SpectralGrid::isRegular(unsigned int spwId) const
{
  if(wrongSpwId(spwId)) return false;
  if(v_chanSep_[spwId] == 0.0) return false;
  return true;
}

std::string SpectralGrid::getSidebandSide(unsigned int spwId) const
{
  if(!wrongSpwId(spwId)) {
    if(vv_assocSpwId_[spwId].size() == 0) {
      /* cout << "WARNING: the spectral window with the identifier "<< spwId
       << " has no associated spectral window "<< endl; */
      return "";
    }
    if(v_sidebandSide_[spwId] == NOSB) return "NoSB";
    if(v_sidebandSide_[spwId] == LSB) return "LSB";
    if(v_sidebandSide_[spwId] == USB) return "USB";
  }
  return "";
}

std::string SpectralGrid::getSidebandType(unsigned int spwId) const
{
  if(!wrongSpwId(spwId)) {
    if(vv_assocSpwId_[spwId].size() == 0) {
      /* cout << "WARNING: the spectral window with the identifier "<< spwId
       << " has no associated spectral window "<< endl; */
    }
    return "";
    if(v_sidebandType_[spwId] == DSB) return " DSB";
    if(v_sidebandType_[spwId] == SSB) return " SSB";
    if(v_sidebandType_[spwId] == TWOSB) return " 2SB";
  }
  return "";
}

std::string SpectralGrid::getSideband(unsigned int spwId) const
{
  if(!wrongSpwId(spwId)) {
    if(vv_assocSpwId_[spwId].size() == 0) {
      /* cout << "WARNING: the spectral window with the identifier "<< spwId
       << " has no associated spectral window "<< endl; */
    }
    return "";
    if(getSidebandSide(spwId) == "NoSB") {
      return getSidebandSide(spwId);
    } else {
      std::string sbTypeSide = getSidebandSide(spwId) + " with type ";
      return (sbTypeSide + getSidebandType(spwId));
    }
  }
  return "";
}

std::vector<std::string> SpectralGrid::getAssocNature(unsigned int spwId) const
{
  if(!wrongSpwId(spwId)) {
    if(vv_assocNature_[spwId].size() == 0) {
      /* cout << "WARNING: the spectral window with the identifier "<< spwId
       << " has no associated spectral window "<< endl; */
    }
    return vv_assocNature_[spwId];
  }
  std::vector<std::string> v_dummyVector;
  return v_dummyVector;
}

std::vector<unsigned int> SpectralGrid::getAssocSpwId(unsigned int spwId) const
{
  if(!wrongSpwId(spwId)) {
    if(vv_assocSpwId_[spwId].size() == 0) {
      /* cout << "WARNING: the spectral window with the identifier "<< spwId
       << " has no associated spectral window "<< endl; */
      /* cout << "vv_assocSpwId_[" <<spwId<<"]=" <<  vv_assocSpwId_[spwId][0] << endl; */
    }
    return vv_assocSpwId_[spwId];

  }
  std::vector<unsigned int> v_dummyVector;
  return v_dummyVector;
}

std::vector<unsigned int> SpectralGrid::getAssocSpwIds(const std::vector<unsigned int> &spwIds) const
{

  unsigned int spwId;
  std::vector<unsigned int> assoc_spwIds;

  for(unsigned int n = 0; n < spwIds.size(); n++) {

    spwId = spwIds[n];

    if(!wrongSpwId(spwId)) {
      if(vv_assocSpwId_[spwId].size() == 0) {
        /* cout << "WARNING: the spectral window with the identifier "<< spwId
         << " has no associated spectral window "<< endl; */
        /* cout << "vv_assocSpwId_[" <<spwId<<"]=" <<  vv_assocSpwId_[spwId][0] << endl; */
      }

      assoc_spwIds.push_back((vv_assocSpwId_[spwId])[0]);

    } else {
      assoc_spwIds.push_back(spwId);
    }

  }

  return assoc_spwIds;

}

double SpectralGrid::getLoFrequency() const
{
  return v_loFreq_[0];
}

double SpectralGrid::getLoFrequency(unsigned int spwId) const
{
  if(wrongSpwId(spwId)) return 32767.;
  return v_loFreq_[spwId];
}

ATM_NAMESPACE_END

