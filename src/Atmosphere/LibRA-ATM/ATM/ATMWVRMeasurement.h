#ifndef _ATM_WVRMEASUREMENT_H
#define _ATM_WVRMEASUREMENT_H
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
 * "@(#) $Id: ATMWVRMeasurement.h Exp $"
 *
 * who       when      what
 * --------  --------  ----------------------------------------------
 * pardo     24/03/09  created
 */

#ifndef __cplusplus
#error "This is a C++ include file and cannot be used from plain C"
#endif

#include "ATMAngle.h"
#include "ATMCommon.h"
#include "ATMLength.h"
#include "ATMTemperature.h"
#include <vector>
#include <math.h>

using std::vector;

ATM_NAMESPACE_BEGIN

/*! \brief This is an auxiliary class that allows to create objects corresponding
 *        to measurements (and their analysis) of a multichannel water vapor radiometer system.
 */
class WVRMeasurement
{
public:
  /** Class constructor with no radiometric channels */
  WVRMeasurement();
  /** The basic class constructor, with the elevation angle corresponding to the measurement, and the measurements themselves
   in the channels of the water vapor radiometer system  */
  WVRMeasurement(const Angle &elevation,
                 const vector<Temperature> &measuredSkyBrightness);
  /** A more general class constructor, with the elevation angle corresponding to the measurement, the measurements themselves
   in the channels of the water vapor radiometer system, the fitted values in those channels and the sigma of the fit. This
   constructor is not intended to be used by the user, but by water vapor retrieval procedures inside the SkyStatus class */
  WVRMeasurement(const Angle &elevation,
                 const vector<Temperature> &measuredSkyBrightness,
                 const vector<Temperature> &fittedSkyBrightness,
                 const Length &retrievedWaterVaporColumn,
                 const Temperature &sigma_fittedSkyBrightness);

  /** Destructor */
  virtual ~WVRMeasurement();

  /** Accessor to elevation */
  Angle getElevation() const { return elevation_; }
  /** Accessor to air mass */
   double getAirMass() const { return 1.0 / sin(elevation_.get()); }
  /** Accessor to measured sky brightness temperature */
   vector<Temperature> getmeasuredSkyBrightness() const { return v_measuredSkyBrightness_; }
  /** Accessor to fitted sky brightness temperatures */
   vector<Temperature> getfittedSkyBrightness() const { return v_fittedSkyBrightness_; }
  /** Setter of fitted sky brightness temperatures */
   void setfittedSkyBrightness(const vector<Temperature> &a) { v_fittedSkyBrightness_ = a; }
  /** Accessor to retrieved water vapor column */
   Length getretrievedWaterVaporColumn() const { return retrievedWaterVaporColumn_; }
  /** Setter of retrieved water vapor column */
   void setretrievedWaterVaporColumn(const Length &a) { retrievedWaterVaporColumn_ = a; }
  /** Accessor to sigma of the fit (in K) */
   Temperature getSigmaFit() const { return sigma_fittedSkyBrightness_; }
  /** Setter of sigma of the fit (in K) */
   void setSigmaFit(const Temperature &a) { sigma_fittedSkyBrightness_ = a; }

protected:
  Angle elevation_;
  vector<Temperature> v_measuredSkyBrightness_;
  vector<Temperature> v_fittedSkyBrightness_; // !< Fitted sky brightness temperatures over all WVR channels for each event
  Length retrievedWaterVaporColumn_; // !< Retrieved zenith water vapor column for each event
  Temperature sigma_fittedSkyBrightness_; // !< Sigma on the fitted sky brightness temperatures (average sigma over the WVR channels for each event).
}; // class WVRMeasurement

ATM_NAMESPACE_END

#endif /*!_ATM_WVRMEASUREMENT_H*/
