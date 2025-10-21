#ifndef _ATM_REFRACTIVEINDEXPROFILE_H
#define _ATM_REFRACTIVEINDEXPROFILE_H
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
 * "@(#) $Id: ATMRefractiveIndexProfile.h Exp $"
 *
 * who       when      what
 * --------  --------  ----------------------------------------------
 * pardo     24/03/09  created
 */

#ifndef __cplusplus
#error This is a C++ include file and cannot be used from plain C
#endif

#include "ATMAngle.h"
#include "ATMCommon.h"
#include "ATMInverseLength.h"
#include "ATMOpacity.h"
#include "ATMProfile.h"
#include "ATMSpectralGrid.h"
#include "ATMRefractiveIndex.h"

#include <complex>

ATM_NAMESPACE_BEGIN

/**  \brief Profile of the absorption and Phase coefficient(s) at given frequency(ies) for an
 *   atmospheric profile (P/T/gas densities).
 *
 *   Using the ATM library the vertical profile of the absorption
 *   coefficient corresponding to an atmospheric P/T/gas profile are build
 *   The absorption coefficient in each layer can be used
 *   later for radiative transfer purposes. The atmospheric P/T/gas profile is inherited
 *   from an object of type AtmProfile. The absorption coefficient profile
 *   is provided separately for O2 lines, H2O lines, CO lines, O3 lines, N2O lines, NO2 lines,
 *   SO2 lines dry and wet continuum. The total absorption coefficient profile (addition of
 *   those just described) is also available, since it is the only profile really needed
 *   for the WaterVaporRetrieval class which derives from this
 *   RefractiveIndexProfile class.
 */
class RefractiveIndexProfile: public AtmProfile, public SpectralGrid
{
public:

  //@{
  /** A full constructor for the case of a profile of absorption coefficients
   *  at a single frequency
   */
  RefractiveIndexProfile(const Frequency &frequency,
                         const AtmProfile &atmProfile);

  /** A full constructor for the case of a profile of absorption coefficients
   *  for a set of frequency points.
   */
  RefractiveIndexProfile(const SpectralGrid &spectralGrid,
                         const AtmProfile &atmProfile);

  /** A copy constructor for deep copy
   */
  RefractiveIndexProfile(const RefractiveIndexProfile &);

  RefractiveIndexProfile();

  virtual ~RefractiveIndexProfile();

  //@}

  //@{
  /** Setter to update the AtmProfile and RefractiveIndexProfile if some basic atmospheric parameter has changed.
   * @pre   an RefractiveIndexProfile already exists
   * @param altitude          the new altitude, a Length
   * @param groundPressure    the Pressure at the ground level
   * @param groundTemperature the Temperature at the ground level
   * @param tropoLapseRate    the tropospheric lapse rate
   * @param relativeHumidity  the relative Humidity
   * @param wvScaleHeight     the scale height of the water vapor, a Length
   * @return true if the atmospheric profile has been updated, else false because the basic parameters have not changed
   * @post   the atmospheric profile, and accordingly the absorption phase profile, have been updated, unless none
   *         of all the input basic parameters correspond to new values.
   *
   * \note  there is an overriding on this method in the sub-class in WaterVaporRadiance sub-class. Hence this method
   *        must not be overloaded in this RefractiveIndexProfile class.
   */
  bool setBasicAtmosphericParameters(const Length &altitude,
                                     const Pressure &groundPressure,
                                     const Temperature &groundTemperature,
                                     double tropoLapseRate,
                                     const Humidity &relativeHumidity,
                                     const Length &wvScaleHeight);
  //@}

  //@{

  /** Setter of new spectral windows (adds frequencies to the SpectralGrid and calculates the corresponding absorption coefficients) */
  inline void addNewSpectralWindow(unsigned int numChan,
                                   unsigned int refChan,
                                   const Frequency &refFreq,
                                   const Frequency &chanSep)
  {
    SpectralGrid::add(numChan, refChan, refFreq, chanSep);
    mkRefractiveIndexProfile();
  }
  /** Setter of new spectral windows (adds frequencies to the SpectralGrid and calculates the corresponding absorption coefficients) */
  inline void addNewSpectralWindow(unsigned int numChan,
                                   unsigned int refChan,
                                   double* chanFreq,
                                   Frequency::Units freqUnits)
  {
    SpectralGrid::add(numChan, refChan, chanFreq, freqUnits);
    mkRefractiveIndexProfile();
  }
  /** Setter of new spectral windows (adds frequencies to the SpectralGrid and calculates the corresponding absorption coefficients) */
  inline void addNewSpectralWindow(unsigned int numChan,
                                   double refFreq,
                                   double* chanFreq,
                                   Frequency::Units freqUnits)
  {
    SpectralGrid::add(numChan, refFreq, chanFreq, freqUnits);
    mkRefractiveIndexProfile();
  }
  /** Setter of new spectral windows (adds frequencies to the SpectralGrid and calculates the corresponding absorption coefficients) */
  inline void addNewSpectralWindow(unsigned int numChan,
                                   double refFreq,
                                   const vector<double> &chanFreq,
                                   Frequency::Units freqUnits)
  {
    SpectralGrid::add(numChan, refFreq, chanFreq, freqUnits);
    mkRefractiveIndexProfile();
  }
  /** Setter of new spectral windows (adds frequencies to the SpectralGrid and calculates the corresponding absorption coefficients) */
  inline void addNewSpectralWindow(unsigned int numChan,
                                   unsigned int refChan,
                                   const Frequency &refFreq,
                                   const Frequency &chanSep,
                                   const Frequency &intermediateFreq,
                                   const SidebandSide &sbSide,
                                   const SidebandType &sbType)
  {
    SpectralGrid::add(numChan,
                      refChan,
                      refFreq,
                      chanSep,
                      intermediateFreq,
                      sbSide,
                      sbType);
    mkRefractiveIndexProfile();
  }
  inline void addNewSpectralWindow(const vector<Frequency> &chanFreq)
  {
    SpectralGrid::add(chanFreq);
    mkRefractiveIndexProfile();
  }

  /** Accessor to total number of Frequency points in the SpectralGrid object */
  inline unsigned int getNumIndividualFrequencies() const
  {
    return v_chanFreq_.size();
  }

  /** Accessor to the single frequency (or the frequency of the first grid point
   *  in case of a spectrum) in Hz (SI default unit)*/
  inline double getFrequency() const
  {
    return getChanFreq(0).get();
  }

  /** Accessor to the single frequency (or the frequency of the first grid point
   *  in case of a spectrum) in specified unit
   */
  inline double getFrequency(Frequency::Units freqUnits) const
  {
    return getChanFreq(0).get(freqUnits);
  }
  //@}

  //@{
  /** Accessor to get H2O lines Absorption Coefficient at layer nl, for single frequency RefractiveIndexProfile object */
  InverseLength getAbsH2OLines(unsigned int nl) const
  {
    return InverseLength(imag((vv_N_H2OLinesPtr_[0]->at(nl))), InverseLength::UnitInverseMeter);
  }
  /** Accessor to get H2O lines Absorption Coefficient at layer nl and frequency channel nf, for RefractiveIndexProfile object with a spectral grid */
  InverseLength getAbsH2OLines(unsigned int nf, unsigned int nl) const
  {
    return InverseLength(imag((vv_N_H2OLinesPtr_[nf]->at(nl))), InverseLength::UnitInverseMeter);
  }
  /** Accessor to get H2O Continuum Absorption Coefficient at layer nl, spectral window spwid and channel nf */
  InverseLength getAbsH2OLines(unsigned int spwid,
                               unsigned int nf,
                               unsigned int nl) const
  {
    unsigned int j = v_transfertId_[spwid] + nf;
    return InverseLength(imag((vv_N_H2OLinesPtr_[j]->at(nl))), InverseLength::UnitInverseMeter);
  }

  /** Accessor to get H2O Continuum Absorption Coefficient at layer nl, for single frequency RefractiveIndexProfile object */
  InverseLength getAbsH2OCont(unsigned int nl) const
  {
    return InverseLength(imag((vv_N_H2OContPtr_[0]->at(nl))), InverseLength::UnitInverseMeter);
  }
  /** Accessor to get H2O Continuum Absorption Coefficient at layer nl and frequency channel nf, for RefractiveIndexProfile object with a spectral grid */
  InverseLength getAbsH2OCont(unsigned int nf, unsigned int nl) const
  {
    return InverseLength(imag((vv_N_H2OContPtr_[nf]->at(nl))), InverseLength::UnitInverseMeter);
  }
  /** Accessor to get H2O Continuum Absorption Coefficient at layer nl, spectral window spwid and channel nf */
  InverseLength getAbsH2OCont(unsigned int spwid,
                              unsigned int nf,
                              unsigned int nl) const
  {
    unsigned int j = v_transfertId_[spwid] + nf;
    return InverseLength(imag((vv_N_H2OContPtr_[j]->at(nl))), InverseLength::UnitInverseMeter);
  }

  /** Function to retrieve O2 lines Absorption Coefficient at layer nl, for single frequency RefractiveIndexProfile object */
  InverseLength getAbsO2Lines(unsigned int nl) const
  {
    return InverseLength(imag((vv_N_O2LinesPtr_[0]->at(nl))), InverseLength::UnitInverseMeter);
  }
  /** Function to retrieve O2 lines Absorption Coefficient at layer nl and frequency channel nf, for RefractiveIndexProfile object with a spectral grid */
  InverseLength getAbsO2Lines(unsigned int nf, unsigned int nl) const
  {
    return InverseLength(imag((vv_N_O2LinesPtr_[nf]->at(nl))), InverseLength::UnitInverseMeter);
  }
  /** Function to retrieve O2 lines Absorption Coefficient at layer nl, spectral window spwid and channel nf */
  InverseLength getAbsO2Lines(unsigned int spwid,
                              unsigned int nf,
                              unsigned int nl) const
  {
    unsigned int j = v_transfertId_[spwid] + nf;
    return InverseLength(imag((vv_N_O2LinesPtr_[j]->at(nl))), InverseLength::UnitInverseMeter);
  }

  /** Function to retrieve Dry continuum Absorption Coefficient at layer nl, for single frequency RefractiveIndexProfile object */
  InverseLength getAbsDryCont(unsigned int nl) const
  {
    return InverseLength(imag((vv_N_DryContPtr_[0]->at(nl))), InverseLength::UnitInverseMeter);
  }
  /** Function to retrieve Dry continuum Absorption Coefficient at layer nl and frequency channel nf, for RefractiveIndexProfile object with a spectral grid */
  InverseLength getAbsDryCont(unsigned int nf, unsigned int nl) const
  {
    return InverseLength(imag((vv_N_DryContPtr_[nf]->at(nl))), InverseLength::UnitInverseMeter);
  }
  /** Function to retrieve Dry continuum Absorption Coefficient at layer nl, spectral window spwid and channel nf */
  InverseLength getAbsDryCont(unsigned int spwid,
                              unsigned int nf,
                              unsigned int nl) const
  {
    unsigned int j = v_transfertId_[spwid] + nf;
    return InverseLength(imag((vv_N_DryContPtr_[j]->at(nl))), InverseLength::UnitInverseMeter);
  }

  /** Function to retrieve O3 lines Absorption Coefficient at layer nl, for single frequency RefractiveIndexProfile object */
  InverseLength getAbsO3Lines(unsigned int nl) const
  {
    return InverseLength(imag((vv_N_O3LinesPtr_[0]->at(nl))), InverseLength::UnitInverseMeter);
  }
  /** Function to retrieve O3 lines Absorption Coefficient at layer nl and frequency channel nf, for RefractiveIndexProfile object with a spectral grid */
  InverseLength getAbsO3Lines(unsigned int nf, unsigned int nl) const
  {
    return InverseLength(imag((vv_N_O3LinesPtr_[nf]->at(nl))), InverseLength::UnitInverseMeter);
  }
  /** Function to retrieve O3 lines Absorption Coefficient at layer nl, spectral window spwid and channel nf */
  InverseLength getAbsO3Lines(unsigned int spwid,
                              unsigned int nf,
                              unsigned int nl) const
  {
    unsigned int j = v_transfertId_[spwid] + nf;
    return InverseLength(imag((vv_N_O3LinesPtr_[j]->at(nl))), InverseLength::UnitInverseMeter);
  }

  /** Function to retrieve CO lines Absorption Coefficient at layer nl, for single frequency RefractiveIndexProfile object */
  InverseLength getAbsCOLines(unsigned int nl) const
  {
    return InverseLength(imag((vv_N_COLinesPtr_[0]->at(nl))), InverseLength::UnitInverseMeter);
  }
  /** Function to retrieve CO lines Absorption Coefficient at layer nl and frequency channel nf, for RefractiveIndexProfile object with a spectral grid */
  InverseLength getAbsCOLines(unsigned int nf, unsigned int nl) const
  {
    return InverseLength(imag((vv_N_COLinesPtr_[nf]->at(nl))), InverseLength::UnitInverseMeter);
  }
  /** Function to retrieve CO lines Absorption Coefficient at layer nl, spectral window spwid and channel nf */
  InverseLength getAbsCOLines(unsigned int spwid,
                              unsigned int nf,
                              unsigned int nl) const
  {
    unsigned int j = v_transfertId_[spwid] + nf;
    return InverseLength(imag((vv_N_COLinesPtr_[j]->at(nl))), InverseLength::UnitInverseMeter);
  }



  /** Function to retrieve N2O lines Absorption Coefficient at layer nl, for single frequency RefractiveIndexProfile object */
  InverseLength getAbsN2OLines(unsigned int nl) const
  {
    return InverseLength(imag((vv_N_N2OLinesPtr_[0]->at(nl))), InverseLength::UnitInverseMeter);
  }
  /** Function to retrieve N2O lines Absorption Coefficient at layer nl and frequency channel nf, for RefractiveIndexProfile object with a spectral grid */
  InverseLength getAbsN2OLines(unsigned int nf, unsigned int nl) const
  {
    return InverseLength(imag((vv_N_N2OLinesPtr_[nf]->at(nl))), InverseLength::UnitInverseMeter);
  }
  /** Function to retrieve N2O lines Absorption Coefficient at layer nl, spectral window spwid and channel nf */
  InverseLength getAbsN2OLines(unsigned int spwid,
                               unsigned int nf,
                               unsigned int nl) const
  {
    unsigned int j = v_transfertId_[spwid] + nf;
    return InverseLength(imag((vv_N_N2OLinesPtr_[j]->at(nl))), InverseLength::UnitInverseMeter);
  }

  /** Function to retrieve NO2 lines Absorption Coefficient at layer nl, for single frequency RefractiveIndexProfile object */
  InverseLength getAbsNO2Lines(unsigned int nl) const
  {
    return InverseLength(imag((vv_N_NO2LinesPtr_[0]->at(nl))), InverseLength::UnitInverseMeter);
  }
  /** Function to retrieve NO2 lines Absorption Coefficient at layer nl and frequency channel nf, for RefractiveIndexProfile object with a spectral grid */
  InverseLength getAbsNO2Lines(unsigned int nf, unsigned int nl) const
  {
    return InverseLength(imag((vv_N_NO2LinesPtr_[nf]->at(nl))), InverseLength::UnitInverseMeter);
  }
  /** Function to retrieve NO2 lines Absorption Coefficient at layer nl, spectral window spwid and channel nf */
  InverseLength getAbsNO2Lines(unsigned int spwid,
                               unsigned int nf,
                               unsigned int nl) const
  {
    unsigned int j = v_transfertId_[spwid] + nf;
    return InverseLength(imag((vv_N_NO2LinesPtr_[j]->at(nl))), InverseLength::UnitInverseMeter);
  }


  /** Function to retrieve SO2 lines Absorption Coefficient at layer nl, for single frequency RefractiveIndexProfile object */
  InverseLength getAbsSO2Lines(unsigned int nl) const
  {
    return InverseLength(imag((vv_N_SO2LinesPtr_[0]->at(nl))), InverseLength::UnitInverseMeter);
  }
  /** Function to retrieve SO2 lines Absorption Coefficient at layer nl and frequency channel nf, for RefractiveIndexProfile object with a spectral grid */
  InverseLength getAbsSO2Lines(unsigned int nf, unsigned int nl) const
  {
    return InverseLength(imag((vv_N_SO2LinesPtr_[nf]->at(nl))), InverseLength::UnitInverseMeter);
  }
  /** Function to retrieve SO2 lines Absorption Coefficient at layer nl, spectral window spwid and channel nf */
  InverseLength getAbsSO2Lines(unsigned int spwid,
                               unsigned int nf,
                               unsigned int nl) const
  {
    unsigned int j = v_transfertId_[spwid] + nf;
    return InverseLength(imag((vv_N_SO2LinesPtr_[j]->at(nl))), InverseLength::UnitInverseMeter);
  }





  /** Function to retrieve total Dry Absorption Coefficient at layer nl, for single frequency RefractiveIndexProfile object */
  InverseLength getAbsTotalDry(unsigned int nl) const
  {
    return getAbsTotalDry(nl, 0);
  }
  /** Function to retrieve total Dry Absorption Coefficient at layer nl and frequency channel nf, for RefractiveIndexProfile object with a spectral grid */
  InverseLength getAbsTotalDry(unsigned int nf, unsigned int nl) const
  {
    return InverseLength(imag(
			      vv_N_O2LinesPtr_[nf]->at(nl)
			      + vv_N_DryContPtr_[nf]->at(nl) + vv_N_O3LinesPtr_[nf]->at(nl)
			      + vv_N_COLinesPtr_[nf]->at(nl) + vv_N_N2OLinesPtr_[nf]->at(nl)
			      + vv_N_NO2LinesPtr_[nf]->at(nl) + vv_N_SO2LinesPtr_[nf]->at(nl)), InverseLength::UnitInverseMeter);
  }
  /** Function to retrieve total Dry Absorption Coefficient at layer nl, spectral window spwid and channel nf */
  InverseLength getAbsTotalDry(unsigned int spwid,
                               unsigned int nf,
                               unsigned int nl) const
  {
    unsigned int j = v_transfertId_[spwid] + nf;
    return InverseLength(imag(
			      vv_N_O2LinesPtr_[j]->at(nl)
			      + vv_N_DryContPtr_[j]->at(nl) + vv_N_O3LinesPtr_[j]->at(nl)
			      + vv_N_COLinesPtr_[j]->at(nl) + vv_N_N2OLinesPtr_[j]->at(nl)
			      + vv_N_NO2LinesPtr_[j]->at(nl) + vv_N_SO2LinesPtr_[j]->at(nl)), InverseLength::UnitInverseMeter);
  }

  /** Function to retrieve total Wet Absorption Coefficient at layer nl, for single frequency RefractiveIndexProfile object */
  InverseLength getAbsTotalWet(unsigned int nl) const
  {
    return getAbsTotalWet(nl, 0);
  }
  /** Function to retrieve total Wet Absorption Coefficient at layer nl and frequency channel nf, for RefractiveIndexProfile object with a spectral grid */
  InverseLength getAbsTotalWet(unsigned int nf, unsigned int nl) const
  {
    return InverseLength(imag((vv_N_H2OLinesPtr_[nf]->at(nl)
        + vv_N_H2OContPtr_[nf]->at(nl))), InverseLength::UnitInverseMeter);
  }
  /** Function to retrieve total Wet Absorption Coefficient at layer nl, spectral window spwid and channel nf */
  InverseLength getAbsTotalWet(unsigned int spwid,
                               unsigned int nf,
                               unsigned int nl) const
  {
    unsigned int j = v_transfertId_[spwid] + nf;
    return InverseLength(imag((vv_N_H2OLinesPtr_[j]->at(nl)
        + vv_N_H2OContPtr_[j]->at(nl))), InverseLength::UnitInverseMeter);
  }

  Opacity getAverageO2LinesOpacity(unsigned int spwid);
  Opacity getAverageO3LinesOpacity(unsigned int spwid);
  Opacity getAverageN2OLinesOpacity(unsigned int spwid);
  Opacity getAverageNO2LinesOpacity(unsigned int spwid);
  Opacity getAverageSO2LinesOpacity(unsigned int spwid);
  Opacity getAverageCOLinesOpacity(unsigned int spwid);
  Opacity getAverageDryContOpacity(unsigned int spwid);

  /** Function to retrieve the integrated Dry Opacity along the atmospheric path for single frequency RefractiveIndexProfile object */
  Opacity getDryOpacity(){return getDryOpacity(0);}
  Opacity getDryOpacityUpTo(Length refalti){return getDryOpacityUpTo(0,refalti);}
  /** Function to retrieve the integrated Dry Opacity along the atmospheric path for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Opacity getDryOpacity(unsigned int nc);
  Opacity getDryOpacityUpTo(unsigned int nc, Length refalti);
  Opacity getDryOpacity(unsigned int spwid, unsigned int nc){ if(!spwidAndIndexAreValid(spwid, nc)) return Opacity(-999.0); return getDryOpacity(v_transfertId_[spwid] + nc);}
  Opacity getDryOpacityUpTo(unsigned int spwid, unsigned int nc, Length refalti)
  {
    if(!spwidAndIndexAreValid(spwid, nc))return Opacity(-999.0);
    return getDryOpacityUpTo(v_transfertId_[spwid] + nc, refalti);
  }
  Opacity getAverageDryOpacity(unsigned int spwid);

  /** Function to retrieve the integrated Dry Continuum Opacity along the atmospheric path for single frequency RefractiveIndexProfile object */
  Opacity getDryContOpacity();
  /** Function to retrieve the integrated Dry Continuum Opacity along the atmospheric path for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Opacity getDryContOpacity(unsigned int nc);
  Opacity getDryContOpacity(unsigned int spwid, unsigned int nc);
  /** Function to retrieve the integrated O2 Lines Opacity along the atmospheric path for single frequency RefractiveIndexProfile object */
  Opacity getO2LinesOpacity();
  /** Function to retrieve the integrated O2 Lines Opacity along the atmospheric path for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Opacity getO2LinesOpacity(unsigned int nc);
  Opacity getO2LinesOpacity(unsigned int spwid, unsigned int nc);
  /** Function to retrieve the integrated O3 Lines Opacity along the atmospheric path for single frequency RefractiveIndexProfile object */
  Opacity getO3LinesOpacity();
  /** Function to retrieve the integrated O3 Lines Opacity along the atmospheric path for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Opacity getO3LinesOpacity(unsigned int nc);
  Opacity getO3LinesOpacity(unsigned int spwid, unsigned int nc);
  /** Function to retrieve the integrated CO Lines Opacity along the atmospheric path for single frequency RefractiveIndexProfile object */
  Opacity getCOLinesOpacity();
  /** Function to retrieve the integrated CO Lines Opacity along the atmospheric path for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Opacity getCOLinesOpacity(unsigned int nc);
  Opacity getCOLinesOpacity(unsigned int spwid, unsigned int nc);
  /** Function to retrieve the integrated N2O Lines Opacity along the atmospheric path for single frequency RefractiveIndexProfile object */
  Opacity getN2OLinesOpacity();
  /** Function to retrieve the integrated N2O Lines Opacity along the atmospheric path for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Opacity getN2OLinesOpacity(unsigned int nc);
  Opacity getN2OLinesOpacity(unsigned int spwid, unsigned int nc);
  /** Function to retrieve the integrated NO2 Lines Opacity along the atmospheric path for single frequency RefractiveIndexProfile object */
  Opacity getNO2LinesOpacity();
  /** Function to retrieve the integrated NO2 Lines Opacity along the atmospheric path for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Opacity getNO2LinesOpacity(unsigned int nc);
  Opacity getNO2LinesOpacity(unsigned int spwid, unsigned int nc);
  /** Function to retrieve the integrated SO2 Lines Opacity along the atmospheric path for single frequency RefractiveIndexProfile object */
  Opacity getSO2LinesOpacity();
  /** Function to retrieve the integrated SO2 Lines Opacity along the atmospheric path for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Opacity getSO2LinesOpacity(unsigned int nc);
  Opacity getSO2LinesOpacity(unsigned int spwid, unsigned int nc);

  /** Function to retrieve the integrated Wet Opacity along the atmospheric path
   corresponding to the 1st guess water column (from AtmProfile object) or the
    for single frequency RefractiveIndexProfile object */
  Opacity getWetOpacity() {return getWetOpacity(getGroundWH2O());}
  /** Function to retrieve the integrated Wet Opacity along the atmospheric path
   for a given integratedwatercolumn */
  Opacity getWetOpacity(const Length &integratedwatercolumn);
  /** Function to retrieve the integrated Wet Opacity along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile
   object with a spectral grid */
  Opacity getWetOpacity(const Length &integratedwatercolumn, unsigned int nc);
  /** Function to retrieve the integrated Wet Opacity along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile
   object with a spectral grid */
  Opacity getWetOpacity(const Length &integratedwatercolumn, unsigned int spwid, unsigned int nc);
  /** Function to retrieve the integrated H2O Lines Opacity along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object */


  /** Function to retrieve the integrated H2O Lines Opacity along the atmospheric path
   corresponding to the 1st guess water column (from AtmProfile object)  */
  Opacity getH2OLinesOpacity() {return getH2OLinesOpacity(getGroundWH2O());}
  /** Function to retrieve the integrated H2O Lines Opacity along the atmospheric path
   corresponding to a given water vapor column */
  Opacity getH2OLinesOpacity(const Length &integratedwatercolumn);
  /** Function to retrieve the integrated H2O Lines Opacity along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Opacity getH2OLinesOpacity(const Length &integratedwatercolumn, unsigned int nc);
  /** Function to retrieve the integrated H2O Lines Opacity along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object of a spectral window */
  Opacity getH2OLinesOpacity(const Length &integratedwatercolumn, unsigned int spwid, unsigned int nc);
  /** Function to retrieve the integrated H2O Continuum Opacity along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object */

  /** Function to retrieve the integrated H2O Continuum Opacity along the atmospheric path
   corresponding to the 1st guess water column (from AtmProfile object)  */
  Opacity getH2OContOpacity() {return getH2OContOpacity(getGroundWH2O());}
  /** Function to retrieve the integrated H2O Continuum Opacity along the atmospheric path
   corresponding to a given water vapor column */
  Opacity getH2OContOpacity(const Length &integratedwatercolumn);
  /** Function to retrieve the integrated H2O Continuum Opacity along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Opacity getH2OContOpacity(const Length &integratedwatercolumn, unsigned int nc);
  /** Function to retrieve the integrated H2O Continuum Opacity along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object of a spectral window */
  Opacity getH2OContOpacity(const Length &integratedwatercolumn, unsigned int spwid, unsigned int nc);





  Opacity getAverageWetOpacity(const Length &integratedwatercolumn, unsigned int spwid);
  Opacity getAverageH2OLinesOpacity(const Length &integratedwatercolumn, unsigned int spwid);
  Opacity getAverageH2OContOpacity(const Length &integratedwatercolumn, unsigned int spwid);

  /** Function to retrieve the integrated Atmospheric Phase Delay (Dry part) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object */
  Angle getNonDispersiveDryPhaseDelay();
  /** Function to retrieve the integrated Atmospheric Path length (Dry part) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object */
  Length getNonDispersiveDryPathLength();
  /** Function to retrieve the integrated Atmospheric Phase Delay (Dry part) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Angle getNonDispersiveDryPhaseDelay(unsigned int nc);
  /** Function to retrieve the integrated Atmospheric Path length (Dry part) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Length getNonDispersiveDryPathLength(unsigned int nc);
  /** Function to retrieve the integrated Atmospheric Phase Delay (Dry part) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object with several Spectral Grids */
  Angle getNonDispersiveDryPhaseDelay(unsigned int spwid, unsigned int nc);
  /** Function to retrieve the integrated Atmospheric Path length (Dry part) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with everal Spectral Grids */
  Length getNonDispersiveDryPathLength(unsigned int spwid, unsigned int nc);

  /** Function to retrieve the integrated Atmospheric Phase Delay (Dry part) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object */
  Angle getDispersiveDryPhaseDelay();
  /** Function to retrieve the integrated Atmospheric Path length (Dry part) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object */
  Length getDispersiveDryPathLength();
  /** Function to retrieve the integrated Atmospheric Phase Delay (Dry part) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Angle getDispersiveDryPhaseDelay(unsigned int nc);
  /** Function to retrieve the integrated Atmospheric Path length (Dry part) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Length getDispersiveDryPathLength(unsigned int nc);
  /** Function to retrieve the integrated Atmospheric Phase Delay (Dry part) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object with several Spectral Grids */
  Angle getDispersiveDryPhaseDelay(unsigned int spwid, unsigned int nc);
  /** Function to retrieve the integrated Atmospheric Path length (Dry part) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with everal Spectral Grids */
  Length getDispersiveDryPathLength(unsigned int spwid, unsigned int nc);

  /** Function to retrieve the average integrated Atmospheric Path Length (Dry part) in spectral Window spwid */
  Length getAverageNonDispersiveDryPathLength(unsigned int spwid);
  /** Function to retrieve the average integrated Atmospheric Path Length (Dry part) in spectral Window spwid */
  Length getAverageDispersiveDryPathLength(unsigned int spwid);
  /** Function to retrieve the average integrated Atmospheric Phase Delay (Dry part) in spectral Window spwid */
  Angle getAverageNonDispersiveDryPhaseDelay(unsigned int spwid);
  /** Function to retrieve the average integrated Atmospheric Phase Delay (Dry part) in spectral Window spwid */
  Angle getAverageDispersiveDryPhaseDelay(unsigned int spwid);

  /** Function to retrieve the integrated Atmospheric Phase Delay (due to O2 lines) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object */
  Angle getO2LinesPhaseDelay();
  /** Function to retrieve the integrated Atmospheric Path length (due to O2 Lines) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object */
  Length getO2LinesPathLength();
  /** Function to retrieve the integrated Atmospheric Phase Delay (due to O2 Lines) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Angle getO2LinesPhaseDelay(unsigned int nc);
  /** Function to retrieve the integrated Atmospheric Path length (due to O2 Lines) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Length getO2LinesPathLength(unsigned int nc);
  /** Function to retrieve the integrated Atmospheric Phase Delay (due to O2 Lines) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object with several Spectral Grids */
  Angle getO2LinesPhaseDelay(unsigned int spwid, unsigned int nc);
  /** Function to retrieve the average integrated Atmospheric Phase Delay (due to O2 Lines) in spectral Window spwid */
  Angle getAverageO2LinesPhaseDelay(unsigned int spwid);
  /** Function to retrieve the integrated Atmospheric Path length (due to O2 Lines) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with everal Spectral Grids */
  Length getO2LinesPathLength(unsigned int spwid, unsigned int nc);
  /** Function to retrieve the average integrated Atmospheric Path Length (due to O2 Lines) in spectral Window spwid */
  Length getAverageO2LinesPathLength(unsigned int spwid);

  /** Function to retrieve the integrated Atmospheric Phase Delay (due to O3 lines) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object */
  Angle getO3LinesPhaseDelay();
  /** Function to retrieve the integrated Atmospheric Path length (due to O3 Lines) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object */
  Length getO3LinesPathLength();
  /** Function to retrieve the integrated Atmospheric Phase Delay (due to O3 Lines) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Angle getO3LinesPhaseDelay(unsigned int nc);
  /** Function to retrieve the integrated Atmospheric Path length (due to O3 Lines) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Length getO3LinesPathLength(unsigned int nc);
  /** Function to retrieve the integrated Atmospheric Phase Delay (due to O3 Lines) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object with several Spectral Grids */
  Angle getO3LinesPhaseDelay(unsigned int spwid, unsigned int nc);
  /** Function to retrieve the average integrated Atmospheric Phase Delay (due to O3 Lines) in spectral Window spwid */
  Angle getAverageO3LinesPhaseDelay(unsigned int spwid);
  /** Function to retrieve the integrated Atmospheric Path length (due to O3 Lines) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with everal Spectral Grids */
  Length getO3LinesPathLength(unsigned int spwid, unsigned int nc);
  /** Function to retrieve the average integrated Atmospheric Path Length (due to O3 Lines) in spectral Window spwid */
  Length getAverageO3LinesPathLength(unsigned int spwid);

  /** Function to retrieve the integrated Atmospheric Phase Delay (due to CO lines) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object */
  Angle getCOLinesPhaseDelay();
  /** Function to retrieve the integrated Atmospheric Path length (due to CO Lines) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object */
  Length getCOLinesPathLength();
  /** Function to retrieve the integrated Atmospheric Phase Delay (due to CO Lines) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Angle getCOLinesPhaseDelay(unsigned int nc);
  /** Function to retrieve the integrated Atmospheric Path length (due to CO Lines) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Length getCOLinesPathLength(unsigned int nc);
  /** Function to retrieve the integrated Atmospheric Phase Delay (due to CO Lines) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object with several Spectral Grids */
  Angle getCOLinesPhaseDelay(unsigned int spwid, unsigned int nc);
  /** Function to retrieve the average integrated Atmospheric Phase Delay (due to CO Lines) in spectral Window spwid */
  Angle getAverageCOLinesPhaseDelay(unsigned int spwid);
  /** Function to retrieve the integrated Atmospheric Path length (due to CO Lines) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with everal Spectral Grids */
  Length getCOLinesPathLength(unsigned int spwid, unsigned int nc);
  /** Function to retrieve the average integrated Atmospheric Path Length (due to CO Lines) in spectral Window spwid */
  Length getAverageCOLinesPathLength(unsigned int spwid);




  /** Function to retrieve the integrated Atmospheric Phase Delay (due to N2O lines) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object */
  Angle getN2OLinesPhaseDelay();
  /** Function to retrieve the integrated Atmospheric Path length (due to N2O Lines) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object */
  Length getN2OLinesPathLength();
  /** Function to retrieve the integrated Atmospheric Phase Delay (due to N2O Lines) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Angle getN2OLinesPhaseDelay(unsigned int nc);
  /** Function to retrieve the integrated Atmospheric Path length (due to N2O Lines) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Length getN2OLinesPathLength(unsigned int nc);
  /** Function to retrieve the integrated Atmospheric Phase Delay (due to N2O Lines) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object with several Spectral Grids */
  Angle getN2OLinesPhaseDelay(unsigned int spwid, unsigned int nc);
  /** Function to retrieve the average integrated Atmospheric Phase Delay (due to N2O Lines) in spectral Window spwid */
  Angle getAverageN2OLinesPhaseDelay(unsigned int spwid);
  /** Function to retrieve the integrated Atmospheric Path length (due to N2O Lines) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with everal Spectral Grids */
  Length getN2OLinesPathLength(unsigned int spwid, unsigned int nc);
  /** Function to retrieve the average integrated Atmospheric Path Length (due to N2O Lines) in spectral Window spwid */
  Length getAverageN2OLinesPathLength(unsigned int spwid);



  /** Function to retrieve the integrated Atmospheric Phase Delay (due to NO2 lines) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object */
  Angle getNO2LinesPhaseDelay();
  /** Function to retrieve the integrated Atmospheric Path length (due to NO2 Lines) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object */
  Length getNO2LinesPathLength();
  /** Function to retrieve the integrated Atmospheric Phase Delay (due to NO2 Lines) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Angle getNO2LinesPhaseDelay(unsigned int nc);
  /** Function to retrieve the integrated Atmospheric Path length (due to NO2 Lines) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Length getNO2LinesPathLength(unsigned int nc);
  /** Function to retrieve the integrated Atmospheric Phase Delay (due to NO2 Lines) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object with several Spectral Grids */
  Angle getNO2LinesPhaseDelay(unsigned int spwid, unsigned int nc);
  /** Function to retrieve the average integrated Atmospheric Phase Delay (due to NO2 Lines) in spectral Window spwid */
  Angle getAverageNO2LinesPhaseDelay(unsigned int spwid);
  /** Function to retrieve the integrated Atmospheric Path length (due to NO2 Lines) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with everal Spectral Grids */
  Length getNO2LinesPathLength(unsigned int spwid, unsigned int nc);
  /** Function to retrieve the average integrated Atmospheric Path Length (due to NO2 Lines) in spectral Window spwid */
  Length getAverageNO2LinesPathLength(unsigned int spwid);


  /** Function to retrieve the integrated Atmospheric Phase Delay (due to SO2 lines) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object */
  Angle getSO2LinesPhaseDelay();
  /** Function to retrieve the integrated Atmospheric Path length (due to SO2 Lines) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object */
  Length getSO2LinesPathLength();
  /** Function to retrieve the integrated Atmospheric Phase Delay (due to SO2 Lines) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Angle getSO2LinesPhaseDelay(unsigned int nc);
  /** Function to retrieve the integrated Atmospheric Path length (due to SO2 Lines) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Length getSO2LinesPathLength(unsigned int nc);
  /** Function to retrieve the integrated Atmospheric Phase Delay (due to SO2 Lines) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object with several Spectral Grids */
  Angle getSO2LinesPhaseDelay(unsigned int spwid, unsigned int nc);
  /** Function to retrieve the average integrated Atmospheric Phase Delay (due to SO2 Lines) in spectral Window spwid */
  Angle getAverageSO2LinesPhaseDelay(unsigned int spwid);
  /** Function to retrieve the integrated Atmospheric Path length (due to SO2 Lines) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with everal Spectral Grids */
  Length getSO2LinesPathLength(unsigned int spwid, unsigned int nc);
  /** Function to retrieve the average integrated Atmospheric Path Length (due to SO2 Lines) in spectral Window spwid */
  Length getAverageSO2LinesPathLength(unsigned int spwid);









  /** Function to retrieve the integrated Atmospheric Phase Delay (Dispersive part) along the atmospheric path
   corresponding to the 1st guess water column (from AtmProfile object)
    for single frequency RefractiveIndexProfile object */
  Angle getDispersiveH2OPhaseDelay() {return getDispersiveH2OPhaseDelay(getGroundWH2O());}
  /** Function to retrieve the integrated Atmospheric Phase Delay (Dispersive part) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object */
  Angle getDispersiveH2OPhaseDelay(const Length &integratedwatercolumn);
  /** Function to retrieve the integrated Atmospheric Path length (Dispersive part) along the atmospheric path
   corresponding to the 1st guess water column (from AtmProfile object)
    for single frequency RefractiveIndexProfile object */
  Length getDispersiveH2OPathLength() {return getDispersiveH2OPathLength(getGroundWH2O());}
  /** Function to retrieve the integrated Atmospheric Path length (Dispersive part) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object */
  Length getDispersiveH2OPathLength(const Length &integratedwatercolumn);
  /** Function to retrieve the integrated Atmospheric Phase Delay (Dispersive part) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Angle getDispersiveH2OPhaseDelay(const Length &integratedwatercolumn, unsigned int nc);
  /** Function to retrieve the integrated Atmospheric Path length (Dispersive part) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Length getDispersiveH2OPathLength(const Length &integratedwatercolumn, unsigned int nc);
  /** Function to retrieve the integrated Atmospheric Phase Delay (Dispersive part) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object with several Spectral Grids */
  Angle getDispersiveH2OPhaseDelay(const Length &integratedwatercolumn, unsigned int spwid, unsigned int nc);
  /** Function to retrieve the average integrated Atmospheric Phase Delay (Dispersive part) in spectral Window spwid */
  Angle getAverageDispersiveH2OPhaseDelay(const Length &integratedwatercolumn, unsigned int spwid);
  /** Function to retrieve the integrated Atmospheric Path length (Dispersive part) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with everal Spectral Grids */
  Length getDispersiveH2OPathLength(const Length &integratedwatercolumn, unsigned int spwid, unsigned int nc);
  /** Function to retrieve the average integrated Atmospheric Path Length (Dispersive part) in spectral Window spwid */
  Length getAverageDispersiveH2OPathLength(const Length &integratedwatercolumn, unsigned int spwid);
  /** Function to retrieve the integrated Atmospheric Phase Delay (Non-Dispersive part) along the atmospheric path
   corresponding to the 1st guess water column (from AtmProfile object)
    for single frequency RefractiveIndexProfile object */
  Angle getNonDispersiveH2OPhaseDelay() {return getNonDispersiveH2OPhaseDelay(getGroundWH2O());}
  /** Function to retrieve the integrated Atmospheric Phase Delay (Non-Dispersive part) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object */
  Angle getNonDispersiveH2OPhaseDelay(const Length &integratedwatercolumn);
  /** Function to retrieve the integrated Atmospheric Path Length (Non-Dispersive part) along the atmospheric path
   corresponding to the 1st guess water column (from AtmProfile object)
    for single frequency RefractiveIndexProfile object */
  Length getNonDispersiveH2OPathLength() {return getNonDispersiveH2OPathLength(getGroundWH2O());}
  /** Function to retrieve the integrated Atmospheric Path Length (Non-Dispersive part) along the atmospheric path
   corresponding to a given water vapor column
    for single frequency RefractiveIndexProfile object */
  Length getNonDispersiveH2OPathLength(const Length &integratedwatercolumn);
  /** Function to retrieve the integrated Atmospheric Phase Delay (Non-Dispersive part) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Angle getNonDispersiveH2OPhaseDelay(const Length &integratedwatercolumn, unsigned int nc);
  /** Function to retrieve the integrated Atmospheric Path Length (Non-Dispersive part) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with a spectral grid */
  Length getNonDispersiveH2OPathLength(const Length &integratedwatercolumn, unsigned int nc);
  /** Function to retrieve the integrated Atmospheric Phase Delay (Non-Dispersive part) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with several spectral grids */
  Angle getNonDispersiveH2OPhaseDelay(const Length &integratedwatercolumn, unsigned int spwid, unsigned int nc);
  /** Function to retrieve the average integrated Atmospheric Path Length (Non-Dispersive part) in spectral Window spwid */
  Length getAverageNonDispersiveH2OPathLength(const Length &integratedwatercolumn, unsigned int spwid);
  /** Function to retrieve the average integrated Atmospheric Phase Delay (Non-Dispersive part) in spectral Window spwid */
  Angle getAverageNonDispersiveH2OPhaseDelay(const Length &integratedwatercolumn, unsigned int spwid);
  /** Function to retrieve the integrated Atmospheric Path Length (Non-Dispersive part) along the atmospheric path
   corresponding to a given water vapor column
    for channel nc in an RefractiveIndexProfile object with several spectral grids */
  Length getNonDispersiveH2OPathLength(const Length &integratedwatercolumn, unsigned int spwid, unsigned int nc);

  void updateNewSpectralWindows();

  //@}

protected:

  vector<vector<std::complex<double> >*> vv_N_H2OLinesPtr_; //!< H2O lines absorption coefficient and delay term (rad m^-1,m^-1)
  vector<vector<std::complex<double> >*> vv_N_H2OContPtr_; //!< H2O continuum absorption coefficient and delay term  (rad m^-1,m^-1)
  vector<vector<std::complex<double> >*> vv_N_O2LinesPtr_; //!< O2 lines absorption coefficient and delay term (rad m^-1,m^-1)
  vector<vector<std::complex<double> >*> vv_N_DryContPtr_; //!< Dry continuum absorption coefficient and delay term  (rad m^-1,m^-1)
  vector<vector<std::complex<double> >*> vv_N_O3LinesPtr_; //!< O3 lines absorption coefficient and delay term  (rad m^-1,m^-1)
  vector<vector<std::complex<double> >*> vv_N_COLinesPtr_; //!< CO lines absorption coefficient and delay term  (rad m^-1,m^-1)
  vector<vector<std::complex<double> >*> vv_N_N2OLinesPtr_; //!< N2O lines absorption coefficient and delay term  (rad m^-1,m^-1)
  vector<vector<std::complex<double> >*> vv_N_NO2LinesPtr_; //!< NO2 lines absorption coefficient and delay term  (rad m^-1,m^-1)
  vector<vector<std::complex<double> >*> vv_N_SO2LinesPtr_; //!< SO2 lines absorption coefficient and delay term  (rad m^-1,m^-1)

  /* vecteur de vecteurs ???? */

  /**
   * Method to build the profile of the absorption coefficients,
   */
  void mkRefractiveIndexProfile(); //!<  builds the absorption profiles, returns error code: <0 unsuccessful
  void rmRefractiveIndexProfile(); //!<  deletes all the layer profiles for all the frequencies

  bool updateRefractiveIndexProfile(const Length &altitude,
                                    const Pressure &groundPressure,
                                    const Temperature &groundTemperature,
                                    double tropoLapseRate,
                                    const Humidity &relativeHumidity,
                                    const Length &wvScaleHeight);

  /** Method to assess if a channel frequency index is valid. A spectral grid has a certain number of
   *  frequency points. Hence the index is not valid if it is equal or greater than this number. Would
   *  the spectral grid correspond to a single spectral window this index is the channel frequency index of
   *  that window. When ther are more than a single spectral window in the spectral grid it is recommended
   *  to use the method which include a spectral window identifier.
   *  @param nc channel frequency index in the spectral grid.
   *  @ret   true if the index is smaller than the number of grid points, else false.
   *  @post  update the state of the RefractiveIndexProfile object for ALL the frequency grid points
   *         which have not yet their profiles determined, then return true, else return true
   *         directly, the object being already up-to-date.
   */
  bool chanIndexIsValid(unsigned int nc);
  /** Method to assess if a channel frequency index of a given spectral window is valid.
   *  @param spwid spectral window identifier
   *  @param nc channel frequency index relative to that spectral window (e.g. if a spectral window has
   *         1024 frequency channels then the index nc must be in the range 0 to 1023 in order to return
   *         true).
   *  @ret   true if the index is smaller than the number of grid points for that spectral window, else false.
   *  @post  update the state of the RefractiveIndexProfile object for all the frequency grid points if
   *         there was not yet profiles determined for that index, then return true, else return true
   *         directly, the object being already up-to-date.
   */
  bool spwidAndIndexAreValid(unsigned int spwid, unsigned int idx);
}; // class RefractiveIndexProfile

ATM_NAMESPACE_END

#endif /*!_ATM_REFRACTIVEINDEXPROFILE_H*/
