#ifndef _ATM_SKYSTATUS_H
#define _ATM_SKYSTATUS_H
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
 * "@(#) $Id: ATMSkyStatus.h Exp $"
 *
 * who       when      what
 * --------  --------  ----------------------------------------------
 * pardo     24/03/09  created
 */

#ifndef __cplusplus
#error This is a C++ include file and cannot be used from plain C
#endif

#include "ATMCommon.h"
#include "ATMRefractiveIndexProfile.h"
#include "ATMWaterVaporRadiometer.h"
#include "ATMWVRMeasurement.h"

//#include <math.h>
#include <string>
#include <vector>

using std::string;
using std::vector;

ATM_NAMESPACE_BEGIN

/*! \brief From the layerThickness and layerTemperature arrays (from AtmProfile),
 *   the RefractiveIndexProfile array, and a brightness temperature measured to
 *   the sky at the corresponding frequency, this Class retrieves the water vapor
 *   column that corresponds to the measurement.
 *
 *   The strating point will therefore be an RefractiveIndexProfile object form
 *   which the essential information will be inherited: layerThickness,
 *   layerTemperature vectors (these two inherited themselves from AtmProfile),
 *   and absTotalDry and absTotalWet arrays. The numerical methods
 *   will be an iteration wh2o to match the measured brightness temperature
 */
class SkyStatus: public RefractiveIndexProfile
{
public:

      //@{

  /** The basic constructor. It will allow to perform forward radiative transfer enquires in
   *  the spectral bands of refractiveIndexProfile, and to use some of these bands for water vapor retrievals.
   *  The basic constructor needs at least an  RefractiveIndexProfile object (it brings the spectral information via the
   *  SpectralGrid object inherited by it). Additional inputs can be the
   *  Air mass (default 1.0), sky background temperature (default 2.73 K), and user water vapor column (default 1 mm).
   *  Methods inside the class allow to retrieve the water vapor column and eventually update the user value to the retrieved one.
   *  In order to perform such retrievals, the numerical methods will need as inputs the measured T_EBBs towards the sky in
   *  some spectral bands that will be designed as water vapor radiometry channels. The sideband gains and coupling to the sky
   *  of these vapor radiometry channels should also be inputs of the retrieval methods, otherwise they will be set to default
   *  values. See documentation ef particular methods for more information.
   */
  SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile);
  /** Class constructor with additional inputs */
  SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
            double airMass);
  /** Class constructor with additional inputs */
  SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
            const Temperature &temperatureBackground);
  /** Class constructor with additional inputs */
  SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
            const Length &wh2o);
  /** Class constructor with additional inputs */
  SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
            const Temperature &temperatureBackground,
            double airMass);
  /** Class constructor with additional inputs */
  SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
            double airMass,
            const Temperature &temperatureBackground);
  /** Class constructor with additional inputs */
  SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
            const Length &wh2o,
            double airMass);
  /** Class constructor with additional inputs */
  SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
            double airMass,
            const Length &wh2o);
  /** Class constructor with additional inputs */
  SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
            const Length &wh2o,
            const Temperature &temperatureBackground);
  /** Class constructor with additional inputs */
  SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
            const Temperature &temperatureBackground,
            const Length &wh2o);
  /** Class constructor with additional inputs */
  SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
            const Temperature &temperatureBackground,
            double airMass,
            const Length &wh2o);
  /** Class constructor with additional inputs */
  SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
            const Temperature &temperatureBackground,
            const Length &wh2o,
            double airMass);
  /** Class constructor with additional inputs */
  SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
            double airMass,
            const Temperature &temperatureBackground,
            const Length &wh2o);
  /** Class constructor with additional inputs */
  SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
            double airMass,
            const Length &wh2o,
            const Temperature &temperatureBackground);
  /** Class constructor with additional inputs */
  SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
            const Length &wh2o,
            const Temperature &temperatureBackground,
            double airMass);
  /** Class constructor with additional inputs */
  SkyStatus(const RefractiveIndexProfile &refractiveIndexProfile,
            const Length &wh2o,
            double airMass,
            const Temperature &temperatureBackground);

  /** A copy constructor for deep copy   */
  SkyStatus(const SkyStatus &);

  virtual ~SkyStatus();

  //@}

  //@{


  /** Accessor to get the user water vapor column. This is the water vapor column used for forward
   radiative transfer calculations. It can be equal to wh2o_retrieved_ when this one is avalilable
   and the proper setter or update accessor (retrieveandupdateWaterVapor) has been used. */
  Length getUserWH2O() const { return wh2o_user_; }
  /** Setter for user zenith water vapor column for forward radiative transfer calculations. The user
   zenith water vapor column equals the retrieved zenith water vapor column from H2O radiometers,
   every time the last one is derived with the retrieveandupdateWaterVapor accessor.*/
  void setUserWH2O(const Length &wh2o)
  {
    if(wh2o.get() == wh2o_user_.get()) {
    } else {
      wh2o_user_ = wh2o;
    }
  }
  /** Alternative form of the setter for user zenith water vapor column for forward radiative transfer calculations. The user
   zenith water vapor column equals the retrieved zenith water vapor column from H2O radiometers,
   every time the last one is derived with the retrieveandupdateWaterVapor accessor.*/
  void setUserWH2O(double dwh2o, const string &units)
  {
    Length wh2o(dwh2o, units);
    if(wh2o.get() == wh2o_user_.get()) {
    } else {
      wh2o_user_ = wh2o;
    }
  }
  /** Accessor to get airmass */
  double getAirMass() const { return airMass_; }
  /** Setter for air mass in SkyStatus without performing water vapor retrieval */
  void setAirMass(double airMass)
  {
    if(airMass == airMass_) {
    } else {
      airMass_ = airMass;
    }
  }
  /** Accessor to get sky background temperature */
  Temperature getSkyBackgroundTemperature() const { return skyBackgroundTemperature_; }
  /** Setter for sky background temperature in SkyStatus without
   performing water vapor retrieval */
  void setSkyBackgroundTemperature(const Temperature &skyBackgroundTemperature)
  {
    if(skyBackgroundTemperature.get(Temperature::UnitKelvin) == skyBackgroundTemperature_.get(Temperature::UnitKelvin)) {
    } else {
      skyBackgroundTemperature_
          = Temperature(skyBackgroundTemperature.get(Temperature::UnitKelvin), Temperature::UnitKelvin);
    }
  }

  /** Accesor to the average Equivalent Blackbody Temperature in spectral window 0, for the current conditions
   and a perfect sky coupling */
  Temperature getAverageTebbSky()
  {
    unsigned int n = 0;
    return getAverageTebbSky(n);
  }
  /** Accesor to the average Equivalent Blackbody Temperature in spectral window spwid, for the current conditions
   and a perfect sky coupling */
  Temperature getAverageTebbSky(unsigned int spwid)
  {
    return getAverageTebbSky(spwid,
                             getUserWH2O(),
                             getAirMass(),
                             1.0,
                             getGroundTemperature());
  }
  /** Accesor to the average Equivalent Blackbody Temperature in spectral window 0, for the current
   conditions, except water column wh2o, and a perfect sky coupling */
  Temperature getAverageTebbSky(const Length &wh2o) { return getAverageTebbSky(0, wh2o); }
  /** Accesor to the average Equivalent Blackbody Temperature in spectral window spwid, for the current
   conditions, except water column wh2o, and a perfect sky coupling */
  Temperature getAverageTebbSky(unsigned int spwid, const Length &wh2o)
  {
    return getAverageTebbSky(spwid,
                             wh2o,
                             getAirMass(),
                             1.0,
                             getGroundTemperature());
  }
  /** Accesor to the average Equivalent Blackbody Temperature in spectral window 0, for the current
   conditions, except Air Mass airmass, and a perfect sky coupling */
  Temperature getAverageTebbSky(double airmass) { return getAverageTebbSky(0, airmass); }
  /** Accesor to the average Equivalent Blackbody Temperature in spectral window spwid, for the current
   conditions, except Air Mass airmass, and a perfect sky coupling */
  Temperature getAverageTebbSky(unsigned int spwid, double airmass)
  {
    return getAverageTebbSky(spwid,
                             getUserWH2O(),
                             airmass,
                             1.0,
                             getGroundTemperature());
  }
  /** Accesor to the average Equivalent Blackbody Temperature in spectral window 0, for the current
   conditions, except water column wh2o, introducing a Sky Coupling and Spill Over Temperature */
  Temperature getAverageTebbSky(const Length &wh2o,
                                double skycoupling,
                                const Temperature &Tspill)
  {
    return getAverageTebbSky(0, wh2o, skycoupling, Tspill);
  }
  /** Accesor to the average Equivalent Blackbody Temperature in spectral window spwid, for the current
   conditions, except water column wh2o, introducing a Sky Coupling and Spill Over Temperature */
  Temperature getAverageTebbSky(unsigned int spwid,
                                const Length &wh2o,
                                double skycoupling,
                                const Temperature &Tspill)
  {
    return getAverageTebbSky(spwid, wh2o, getAirMass(), skycoupling, Tspill);
  }
  /** Accesor to the average Equivalent Blackbody Temperature in spectral window 0, for the current
   conditions, except Air Mass airmass, introducing a Sky Coupling and Spill Over Temperature */
  Temperature getAverageTebbSky(double airmass,
                                double skycoupling,
                                const Temperature &Tspill)
  {
    return getAverageTebbSky(0, airmass, skycoupling, Tspill);
  }
  /** Accesor to the average Equivalent Blackbody Temperature in spectral window spwid, for the current
   conditions, except Air Mass airmass, introducing a Sky Coupling and Spill Over Temperature */
  Temperature getAverageTebbSky(unsigned int spwid,
                                double airmass,
                                double skycoupling,
                                const Temperature &Tspill)
  {
    return getAverageTebbSky(spwid, getUserWH2O(), airmass, skycoupling, Tspill);
  }
  /** Accesor to the average Equivalent Blackbody Temperature in spectral window 0 for the current
   conditions, introducing a Sky Coupling and Spill Over Temperature */
  Temperature getAverageTebbSky(double skycoupling, const Temperature &Tspill)
  {
    unsigned int n = 0;
    return getAverageTebbSky(n, skycoupling, Tspill);
  }
  /** Accesor to the average Equivalent Blackbody Temperature in spectral window spwid for the current
   conditions, introducing a Sky Coupling and Spill Over Temperature */
  Temperature getAverageTebbSky(unsigned int spwid,
                                double skycoupling,
                                const Temperature &Tspill)
  {
    return getAverageTebbSky(spwid,
                             getUserWH2O(),
                             getAirMass(),
                             skycoupling,
                             Tspill);
  }
  /** Accesor to the average Equivalent Blackbody Temperature in spectral window spwid, for Water Column
   wh2o, Air Mass airmass, Sky Coupling skycoupling, and Spill Over Temperature Tspill */
  Temperature getAverageTebbSky(unsigned int spwid,
                                const Length &wh2o,
                                double airmass,
                                double skycoupling,
                                const Temperature &Tspill);

  Temperature getAverageTebbSky(unsigned int spwid,
                                const Length &wh2o,
                                double airmass,
                                double skycoupling,
                                double signalgain,
                                const Temperature &Tspill);

  /** Accesor to the Equivalent Blackbody Temperature in spectral window 0 and channel 0, for the currnet
   (user) Water Vapor Column, the current Air Mass, and perfect Sky Coupling to the sky */
  Temperature getTebbSky()
  {
    unsigned int n = 0;
    return getTebbSky(n);
  }
  /** Accesor to the Equivalent Blackbody Temperature in spectral window 0 and channel nc, for the currnet
   (user) Water Vapor Column, the current Air Mass, and perfect Sky Coupling to the sky */
  Temperature getTebbSky(unsigned int nc) // There was reported a bug at Launchpad:
                                          // that the result did not take into account
                                          // the actual column of water. But it is not true.
                                          // The column of water is taken into account as
                                          // seen in the next accessor, that is referred by this one.
  {
    unsigned int n = 0;
    return getTebbSky(n, nc);
  }
  /** Accesor to the Equivalent Blackbody Temperature in spectral window spwid and channel nc, for the currnet
   (user) Water Vapor Column, the current Air Mass, and perfect Sky Coupling to the sky */
  Temperature getTebbSky(unsigned int spwid, unsigned int nc)
  {
    return getTebbSky(spwid,
                      nc,
                      getUserWH2O(),
                      getAirMass(),
                      1.0,
                      Temperature(100, Temperature::UnitKelvin));
  }
  /** Accesor to the Equivalent Blackbody Temperature in spectral window 0 and channel 0, for the currnet
   (user) Water Vapor Column, Air Mass airmass, and perfect Sky Coupling to the sky */
  Temperature getTebbSky(double airmass)
  {
    unsigned int n = 0;
    return getTebbSky(n, airmass);
  }
  /** Accesor to the Equivalent Blackbody Temperature in spectral window 0 and channel nc, for the currnet
   (user) Water Vapor Column, Air Mass airmass, and perfect Sky Coupling to the sky */
  Temperature getTebbSky(unsigned int nc, double airmass)
  {
    unsigned int n = 0;
    return getTebbSky(n, nc, airmass);
  }
  /** Accesor to the Equivalent Blackbody Temperature in spectral window spwid and channel nc, for the currnet
   (user) Water Vapor Column, Air Mass airmass, and perfect Sky Coupling to the sky */
  Temperature getTebbSky(unsigned int spwid, unsigned int nc, double airmass)
  {
    return getTebbSky(spwid, nc, getUserWH2O(), airmass, 1.0, Temperature(100,
                                                                          Temperature::UnitKelvin));
  }
  /** Accesor to the Equivalent Blackbody Temperature in spectral window 0 and channel 0, for Water
   Vapor Column wh2o, the current Air Mass, and perfect Sky Coupling to the sky */
  Temperature getTebbSky(const Length &wh2o)
  {
    unsigned int n = 0;
    return getTebbSky(n, wh2o);
  }
  /** Accesor to the Equivalent Blackbody Temperature in spectral window 0 and channel nc, for Water
   Vapor Column wh2o, the current Air Mass, and perfect Sky Coupling to the sky */
  Temperature getTebbSky(unsigned int nc, const Length &wh2o)
  {
    unsigned int n = 0;
    return getTebbSky(n, nc, wh2o);
  }
  /** Accesor to the Equivalent Blackbody Temperature in spectral window spwid and channel nc, for Water
   Vapor Column wh2o, the current Air Mass, and perfect Sky Coupling to the sky */
  Temperature getTebbSky(unsigned int spwid, unsigned int nc, const Length &wh2o)
  {
    return getTebbSky(spwid, nc, wh2o, getAirMass(), 1.0, Temperature(100, Temperature::UnitKelvin));
  }
  /** Accesor to the Equivalent Blackbody Temperature in spectral window 0 and channel nc, for the current
   (user) Water Vapor Column and Air Mass, Sky Coupling skycoupling, and Spill Over Temperature Tspill */
  Temperature getTebbSky(unsigned int nc,
                         double skycoupling,
                         const Temperature &Tspill)
  {
    unsigned int n = 0;
    return getTebbSky(n, nc, skycoupling, Tspill);
  }
  /** Accesor to the Equivalent Blackbody Temperature in spectral window spwid and channel nc, for the current
   (user) Water Vapor Column and Air Mass, Sky Coupling skycoupling, and Spill Over Temperature Tspill */
  Temperature getTebbSky(unsigned int spwid,
                         unsigned int nc,
                         double skycoupling,
                         const Temperature &Tspill)
  {
    return getTebbSky(spwid,
                      nc,
                      getUserWH2O(),
                      getAirMass(),
                      skycoupling,
                      Tspill);
  }
  /** Accesor to the Equivalent Blackbody Temperature in spectral window 0 and channel nc, for Water
   Vapor Column wh2o, Air Mass airmass, and perfect Sky Coupling to the sky */
  Temperature getTebbSky(unsigned int nc, const Length &wh2o, double airmass)
  {
    unsigned int n = 0;
    return getTebbSky(n, nc, wh2o, airmass);
  }
  /** Accesor to the Equivalent Blackbody Temperature in spectral window spwid and channel nc, for Water
   Vapor Column wh2o, Air Mass airmass, and perfect Sky Coupling to the sky */
  Temperature getTebbSky(unsigned int spwid,
                         unsigned int nc,
                         const Length &wh2o,
                         double airmass)
  {
    return getTebbSky(spwid, nc, wh2o, airmass, 1.0, Temperature(100, Temperature::UnitKelvin));
  }
  /** Accesor to the Equivalent Blackbody Temperature in spectral window 0 and channel nc, for Water
   Vapor Column wh2o, the current Air Mass, Sky Coupling skycoupling, and Spill Over Temperature Tspill */
  Temperature getTebbSky(unsigned int nc,
                         const Length &wh2o,
                         double skycoupling,
                         const Temperature &Tspill)
  {
    unsigned int n = 0;
    return getTebbSky(n, nc, wh2o, skycoupling, Tspill);
  }
  /** Accesor to the Equivalent Blackbody Temperature in spectral window spwid and channel nc, for Water
   Vapor Column wh2o, the current Air Mass, Sky Coupling skycoupling, and Spill Over Temperature Tspill */
  Temperature getTebbSky(unsigned int spwid,
                         unsigned int nc,
                         const Length &wh2o,
                         double skycoupling,
                         const Temperature &Tspill)
  {
    return getTebbSky(spwid, nc, wh2o, getAirMass(), skycoupling, Tspill);
  }
  /** Accesor to the Equivalent Blackbody Temperature in spectral window 0 and channel nc, for the current
   User Water Column, Air Mass airmass, Sky Coupling skycoupling, and Spill Over Temperature Tspill */
  Temperature getTebbSky(unsigned int nc,
                         double airmass,
                         double skycoupling,
                         const Temperature &Tspill)
  {
    unsigned int n = 0;
    return getTebbSky(n, nc, airmass, skycoupling, Tspill);
  }
  /** Accesor to the Equivalent Blackbody Temperature in spectral window spwid and channel nc, for the current
   User Water Column, Air Mass airmass, Sky Coupling skycoupling, and Spill Over Temperature Tspill */
  Temperature getTebbSky(unsigned int spwid,
                         unsigned int nc,
                         double airmass,
                         double skycoupling,
                         const Temperature &Tspill)
  {
    return getTebbSky(spwid, nc, getUserWH2O(), airmass, skycoupling, Tspill);
  }
  /** Accesor to the Equivalent Blackbody Temperature in spectral window 0 and channel nc, for Water Column
   wh2o, Air Mass airmass, Sky Coupling skycoupling, and Spill Over Temperature Tspill */
  Temperature getTebbSky(unsigned int nc,
                         const Length &wh2o,
                         double airmass,
                         double skycoupling,
                         const Temperature &Tspill)
  {
    unsigned int n = 0;
    return getTebbSky(n, nc, wh2o, airmass, skycoupling, Tspill);
  }
  /** Accesor to the Equivalent Blackbody Temperature in spectral window spwid and channel nc, for Water Column
   wh2o, Air Mass airmass, Sky Coupling skycoupling, and Spill Over Temperature Tspill */
  Temperature getTebbSky(unsigned int spwid,
                         unsigned int nc,
                         const Length &wh2o,
                         double airmass,
                         double skycoupling,
                         const Temperature &Tspill);


   /** Accesor to the average Rayleigh-Jeans Temperature in spectral window 0, for the current conditions
   and a perfect sky coupling */
  Temperature getAverageTrjSky()
  {
    unsigned int n = 0;
    return getAverageTrjSky(n);
  }
  /** Accesor to the average Rayleigh-Jeans Temperature in spectral window spwid, for the current conditions
   and a perfect sky coupling */
  Temperature getAverageTrjSky(unsigned int spwid)
  {
    return getAverageTrjSky(spwid,
                             getUserWH2O(),
                             getAirMass(),
                             1.0,
                             getGroundTemperature());
  }
  /** Accesor to the average Rayleigh-Jeans Temperature in spectral window 0, for the current
   conditions, except water column wh2o, and a perfect sky coupling */
  Temperature getAverageTrjSky(const Length &wh2o) { return getAverageTrjSky(0, wh2o); }
  /** Accesor to the average Rayleigh-Jeans Temperature in spectral window spwid, for the current
   conditions, except water column wh2o, and a perfect sky coupling */
  Temperature getAverageTrjSky(unsigned int spwid, const Length &wh2o)
  {
    return getAverageTrjSky(spwid,
                             wh2o,
                             getAirMass(),
                             1.0,
                             getGroundTemperature());
  }
  /** Accesor to the average Rayleigh-Jeans Temperature in spectral window 0, for the current
   conditions, except Air Mass airmass, and a perfect sky coupling */
  Temperature getAverageTrjSky(double airmass) { return getAverageTrjSky(0, airmass); }
  /** Accesor to the average Rayleigh-Jeans Temperature in spectral window spwid, for the current
   conditions, except Air Mass airmass, and a perfect sky coupling */
  Temperature getAverageTrjSky(unsigned int spwid, double airmass)
  {
    return getAverageTrjSky(spwid,
                             getUserWH2O(),
                             airmass,
                             1.0,
                             getGroundTemperature());
  }
  /** Accesor to the average Rayleigh-Jeans Temperature in spectral window 0, for the current
   conditions, except water column wh2o, introducing a Sky Coupling and Spill Over Temperature */
  Temperature getAverageTrjSky(const Length &wh2o,
                                double skycoupling,
                                const Temperature &Tspill)
  {
    return getAverageTrjSky(0, wh2o, skycoupling, Tspill);
  }
  /** Accesor to the average Rayleigh-Jeans Temperature in spectral window spwid, for the current
   conditions, except water column wh2o, introducing a Sky Coupling and Spill Over Temperature */
  Temperature getAverageTrjSky(unsigned int spwid,
                                const Length &wh2o,
                                double skycoupling,
                                const Temperature &Tspill)
  {
    return getAverageTrjSky(spwid, wh2o, getAirMass(), skycoupling, Tspill);
  }
  /** Accesor to the average Rayleigh-Jeans Temperature in spectral window 0, for the current
   conditions, except Air Mass airmass, introducing a Sky Coupling and Spill Over Temperature */
  Temperature getAverageTrjSky(double airmass,
                                double skycoupling,
                                const Temperature &Tspill)
  {
    return getAverageTrjSky(0, airmass, skycoupling, Tspill);
  }
  /** Accesor to the average Rayleigh-Jeans Temperature in spectral window spwid, for the current
   conditions, except Air Mass airmass, introducing a Sky Coupling and Spill Over Temperature */
  Temperature getAverageTrjSky(unsigned int spwid,
                                double airmass,
                                double skycoupling,
                                const Temperature &Tspill)
  {
    return getAverageTrjSky(spwid, getUserWH2O(), airmass, skycoupling, Tspill);
  }
  /** Accesor to the average Rayleigh-Jeans Temperature in spectral window 0 for the current
   conditions, introducing a Sky Coupling and Spill Over Temperature */
  Temperature getAverageTrjSky(double skycoupling, const Temperature &Tspill)
  {
    unsigned int n = 0;
    return getAverageTrjSky(n, skycoupling, Tspill);
  }
  /** Accesor to the average Rayleigh-Jeans Temperature in spectral window spwid for the current
   conditions, introducing a Sky Coupling and Spill Over Temperature */
  Temperature getAverageTrjSky(unsigned int spwid,
                                double skycoupling,
                                const Temperature &Tspill)
  {
    return getAverageTrjSky(spwid,
                             getUserWH2O(),
                             getAirMass(),
                             skycoupling,
                             Tspill);
  }
  /** Accesor to the average Rayleigh-Jeans Temperature in spectral window spwid, for Water Column
   wh2o, Air Mass airmass, Sky Coupling skycoupling, and Spill Over Temperature Tspill */
  Temperature getAverageTrjSky(unsigned int spwid,
                                const Length &wh2o,
                                double airmass,
                                double skycoupling,
                                const Temperature &Tspill);

  Temperature getAverageTrjSky(unsigned int spwid,
                                const Length &wh2o,
                                double airmass,
                                double skycoupling,
                                double signalgain,
                                const Temperature &Tspill);

  /** Accesor to the Rayleigh-Jeans Temperature in spectral window 0 and channel 0, for the currnet
   (user) Water Vapor Column, the current Air Mass, and perfect Sky Coupling to the sky */
  Temperature getTrjSky()
  {
    unsigned int n = 0;
    return getTrjSky(n);
  }
  /** Accesor to the Rayleigh-Jeans Temperature in spectral window 0 and channel nc, for the currnet
   (user) Water Vapor Column, the current Air Mass, and perfect Sky Coupling to the sky */
  Temperature getTrjSky(unsigned int nc) // There was reported a bug at Launchpad:
                                          // that the result did not take into account
                                          // the actual column of water. But it is not true.
                                          // The column of water is taken into account as
                                          // seen in the next accessor, that is referred by this one.
  {
    unsigned int n = 0;
    return getTrjSky(n, nc);
  }
  /** Accesor to the Rayleigh-Jeans Temperature in spectral window spwid and channel nc, for the currnet
   (user) Water Vapor Column, the current Air Mass, and perfect Sky Coupling to the sky */
  Temperature getTrjSky(unsigned int spwid, unsigned int nc)
  {
    return getTrjSky(spwid,
                      nc,
                      getUserWH2O(),
                      getAirMass(),
                      1.0,
                      Temperature(100, Temperature::UnitKelvin));
  }
  /** Accesor to the Rayleigh-Jeans Temperature in spectral window 0 and channel 0, for the currnet
   (user) Water Vapor Column, Air Mass airmass, and perfect Sky Coupling to the sky */
  Temperature getTrjSky(double airmass)
  {
    unsigned int n = 0;
    return getTrjSky(n, airmass);
  }
  /** Accesor to the Rayleigh-Jeans Temperature in spectral window 0 and channel nc, for the currnet
   (user) Water Vapor Column, Air Mass airmass, and perfect Sky Coupling to the sky */
  Temperature getTrjSky(unsigned int nc, double airmass)
  {
    unsigned int n = 0;
    return getTrjSky(n, nc, airmass);
  }
  /** Accesor to the Rayleigh-Jeans Temperature in spectral window spwid and channel nc, for the currnet
   (user) Water Vapor Column, Air Mass airmass, and perfect Sky Coupling to the sky */
  Temperature getTrjSky(unsigned int spwid, unsigned int nc, double airmass)
  {
    return getTrjSky(spwid, nc, getUserWH2O(), airmass, 1.0, Temperature(100,
                                                                         Temperature::UnitKelvin));
  }
  /** Accesor to the Rayleigh-Jeans Temperature in spectral window 0 and channel 0, for Water
   Vapor Column wh2o, the current Air Mass, and perfect Sky Coupling to the sky */
  Temperature getTrjSky(const Length &wh2o)
  {
    unsigned int n = 0;
    return getTrjSky(n, wh2o);
  }
  /** Accesor to the Rayleigh-Jeans Temperature in spectral window 0 and channel nc, for Water
   Vapor Column wh2o, the current Air Mass, and perfect Sky Coupling to the sky */
  Temperature getTrjSky(unsigned int nc, const Length &wh2o)
  {
    unsigned int n = 0;
    return getTrjSky(n, nc, wh2o);
  }
  /** Accesor to the Rayleigh-Jeans Temperature in spectral window spwid and channel nc, for Water
   Vapor Column wh2o, the current Air Mass, and perfect Sky Coupling to the sky */
  Temperature getTrjSky(unsigned int spwid, unsigned int nc, const Length &wh2o)
  {
    return getTrjSky(spwid, nc, wh2o, getAirMass(), 1.0, Temperature(100, Temperature::UnitKelvin));
  }
  /** Accesor to the Rayleigh-Jeans Temperature in spectral window 0 and channel nc, for the current
   (user) Water Vapor Column and Air Mass, Sky Coupling skycoupling, and Spill Over Temperature Tspill */
  Temperature getTrjSky(unsigned int nc,
                         double skycoupling,
                         const Temperature &Tspill)
  {
    unsigned int n = 0;
    return getTrjSky(n, nc, skycoupling, Tspill);
  }
  /** Accesor to the Rayleigh-Jeans Temperature in spectral window spwid and channel nc, for the current
   (user) Water Vapor Column and Air Mass, Sky Coupling skycoupling, and Spill Over Temperature Tspill */
  Temperature getTrjSky(unsigned int spwid,
                         unsigned int nc,
                         double skycoupling,
                         const Temperature &Tspill)
  {
    return getTrjSky(spwid,
                      nc,
                      getUserWH2O(),
                      getAirMass(),
                      skycoupling,
                      Tspill);
  }
  /** Accesor to the Rayleigh-Jeans Temperature in spectral window 0 and channel nc, for Water
   Vapor Column wh2o, Air Mass airmass, and perfect Sky Coupling to the sky */
  Temperature getTrjSky(unsigned int nc, const Length &wh2o, double airmass)
  {
    unsigned int n = 0;
    return getTrjSky(n, nc, wh2o, airmass);
  }
  /** Accesor to the Rayleigh-Jeans Temperature in spectral window spwid and channel nc, for Water
   Vapor Column wh2o, Air Mass airmass, and perfect Sky Coupling to the sky */
  Temperature getTrjSky(unsigned int spwid,
                         unsigned int nc,
                         const Length &wh2o,
                         double airmass)
  {
    return getTrjSky(spwid, nc, wh2o, airmass, 1.0, Temperature(100, Temperature::UnitKelvin));
  }
  /** Accesor to the Rayleigh-Jeans Temperature in spectral window 0 and channel nc, for Water
   Vapor Column wh2o, the current Air Mass, Sky Coupling skycoupling, and Spill Over Temperature Tspill */
  Temperature getTrjSky(unsigned int nc,
                         const Length &wh2o,
                         double skycoupling,
                         const Temperature &Tspill)
  {
    unsigned int n = 0;
    return getTrjSky(n, nc, wh2o, skycoupling, Tspill);
  }
  /** Accesor to the Rayleigh-Jeans Temperature in spectral window spwid and channel nc, for Water
   Vapor Column wh2o, the current Air Mass, Sky Coupling skycoupling, and Spill Over Temperature Tspill */
  Temperature getTrjSky(unsigned int spwid,
                         unsigned int nc,
                         const Length &wh2o,
                         double skycoupling,
                         const Temperature &Tspill)
  {
    return getTrjSky(spwid, nc, wh2o, getAirMass(), skycoupling, Tspill);
  }
  /** Accesor to the Rayleigh-Jeans Temperature in spectral window 0 and channel nc, for the current
   User Water Column, Air Mass airmass, Sky Coupling skycoupling, and Spill Over Temperature Tspill */
  Temperature getTrjSky(unsigned int nc,
                         double airmass,
                         double skycoupling,
                         const Temperature &Tspill)
  {
    unsigned int n = 0;
    return getTrjSky(n, nc, airmass, skycoupling, Tspill);
  }
  /** Accesor to the Rayleigh-Jeans Temperature in spectral window spwid and channel nc, for the current
   User Water Column, Air Mass airmass, Sky Coupling skycoupling, and Spill Over Temperature Tspill */
  Temperature getTrjSky(unsigned int spwid,
                         unsigned int nc,
                         double airmass,
                         double skycoupling,
                         const Temperature &Tspill)
  {
    return getTrjSky(spwid, nc, getUserWH2O(), airmass, skycoupling, Tspill);
  }
  /** Accesor to the Rayleigh-Jeans Temperature in spectral window 0 and channel nc, for Water Column
   wh2o, Air Mass airmass, Sky Coupling skycoupling, and Spill Over Temperature Tspill */
  Temperature getTrjSky(unsigned int nc,
                         const Length &wh2o,
                         double airmass,
                         double skycoupling,
                         const Temperature &Tspill)
  {
    unsigned int n = 0;
    return getTrjSky(n, nc, wh2o, airmass, skycoupling, Tspill);
  }
  /** Accesor to the Rayleigh-Jeans Temperature in spectral window spwid and channel nc, for Water Column
   wh2o, Air Mass airmass, Sky Coupling skycoupling, and Spill Over Temperature Tspill */
  Temperature getTrjSky(unsigned int spwid,
                         unsigned int nc,
                         const Length &wh2o,
                         double airmass,
                         double skycoupling,
                         const Temperature &Tspill);


  /** Accesor to get the integrated zenith Wet Opacity for the current conditions, for a single
   frequency RefractiveIndexProfile object or for the point 0 of spectral window 0 of a
   multi-window RefractiveIndexProfile object.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Opacity getWetOpacity(){unsigned int n = 0; return getWetOpacity(n);}
  Opacity getWetOpacityUpTo(Length refalti){unsigned int n = 0; return getWetOpacityUpTo(n, refalti);}
  /** Accesor to get the integrated zenith Wet Opacity for the current conditions,
   for the point nc of spectral window 0.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Opacity getWetOpacity(unsigned int nc){return getH2OLinesOpacity(nc) + getH2OContOpacity(nc);}
  Opacity getWetOpacityUpTo(unsigned int nc, Length refalti){return getH2OLinesOpacityUpTo(nc, refalti) + getH2OContOpacityUpTo(nc, refalti);}
  /** Accesor to get the integrated zenith Wet Opacity for the current conditions,
   for the point nc of spectral window spwid.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Opacity getWetOpacity(unsigned int spwid, unsigned int nc)
  {
  if(!spwidAndIndexAreValid(spwid, nc)) return (double) -999.0;
    return getWetOpacity(v_transfertId_[spwid] + nc);
  }
  Opacity getWetOpacityUpTo(unsigned int spwid, unsigned int nc, Length refalti)
  {
  if(!spwidAndIndexAreValid(spwid, nc)) return (double) -999.0;
  return getWetOpacityUpTo(v_transfertId_[spwid] + nc, refalti);
  }
  /** Accesor to get the average Wet Opacity for the current conditions,
   in spectral window spwid.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Opacity getAverageWetOpacity(unsigned int spwid)
  {
    return RefractiveIndexProfile::getAverageWetOpacity(getGroundWH2O(),spwid)
        * ((getUserWH2O().get()) / (getGroundWH2O().get()));
  }
  Opacity getAverageH2OLinesOpacity(unsigned int spwid)
  {
    return RefractiveIndexProfile::getAverageH2OLinesOpacity(getGroundWH2O(),spwid)
        * ((getUserWH2O().get()) / (getGroundWH2O().get()));
  }
  Opacity getAverageH2OContOpacity(unsigned int spwid)
  {
    return RefractiveIndexProfile::getAverageH2OContOpacity(getGroundWH2O(),spwid)
        * ((getUserWH2O().get()) / (getGroundWH2O().get()));
  }

  Opacity getDryOpacity()
  {
    return RefractiveIndexProfile::getDryOpacity();
  }
  Opacity getDryOpacity(unsigned int nc)
  {
    return RefractiveIndexProfile::getDryOpacity(nc);
  }
  Opacity getDryOpacity(unsigned int spwid, unsigned int nc)
  {
    return RefractiveIndexProfile::getDryOpacity(spwid, nc);
  }
  Opacity getAverageDryOpacity(unsigned int spwid)
  {
    return RefractiveIndexProfile::getAverageDryOpacity(spwid);
  }

  Opacity getAverageO2LinesOpacity(unsigned int spwid)
  {
    return RefractiveIndexProfile::getAverageO2LinesOpacity(spwid);
  }
  Opacity getAverageO3LinesOpacity(unsigned int spwid)
  {
    return RefractiveIndexProfile::getAverageO3LinesOpacity(spwid);
  }
  Opacity getAverageN2OLinesOpacity(unsigned int spwid)
  {
    return RefractiveIndexProfile::getAverageN2OLinesOpacity(spwid);
  }
  Opacity getAverageNO2LinesOpacity(unsigned int spwid)
  {
    return RefractiveIndexProfile::getAverageNO2LinesOpacity(spwid);
  }
  Opacity getAverageSO2LinesOpacity(unsigned int spwid)
  {
    return RefractiveIndexProfile::getAverageSO2LinesOpacity(spwid);
  }
  Opacity getAverageCOLinesOpacity(unsigned int spwid)
  {
    return RefractiveIndexProfile::getAverageCOLinesOpacity(spwid);
  }
  Opacity getAverageDryContOpacity(unsigned int spwid)
  {
    return RefractiveIndexProfile::getAverageDryContOpacity(spwid);
  }

  Opacity getTotalOpacity()
  {
    return getWetOpacity() + getDryOpacity();
  }
  Opacity getTotalOpacity(unsigned int nc)
  {
    return getWetOpacity(nc) + getDryOpacity(nc);
  }
  Opacity getTotalOpacity(unsigned int spwid, unsigned int nc)
  {
    return getWetOpacity(spwid, nc) + getDryOpacity(spwid, nc);
  }
  Opacity getAverageTotalOpacity(unsigned int spwid)
  {
    return getAverageWetOpacity(spwid) + getAverageDryOpacity(spwid);
  }
  Opacity getTotalOpacityUpTo(Length refalti)
  {
    unsigned int n = 0;
    return getTotalOpacityUpTo(n, refalti);
  }
  //  Opacity getTotalOpacityUpTo(unsigned int nc, Length refalti);
  Opacity getTotalOpacityUpTo(unsigned int nc, Length refalti)
  {
    return getWetOpacityUpTo(nc,refalti) + getDryOpacityUpTo(nc,refalti);
  }
  Opacity getTotalOpacityUpTo(unsigned int spwid, unsigned int nc, Length refalti)
  {
    Opacity wrongOp(-999.0,Opacity::UnitNeper);
    if(!spwidAndIndexAreValid(spwid, nc)) return wrongOp;
    return getTotalOpacityUpTo(v_transfertId_[spwid] + nc, refalti);
  }
  /** Accesor to get the integrated zenith H2O Lines Opacity for the current conditions,
   for a single frequency RefractiveIndexProfile object or for the point 0 of spectral
   window 0 of a multi-window RefractiveIndexProfile object.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Opacity getH2OLinesOpacity(){ unsigned int n = 0; return getH2OLinesOpacity(n);}
  Opacity getH2OLinesOpacityUpTo(Length refalti){ unsigned int n = 0; return getH2OLinesOpacityUpTo(n, refalti);}
  /** Accesor to get the integrated zenith H2O Lines Opacity for the current conditions,
   for the point nc of spectral window 0.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Opacity getH2OLinesOpacity(unsigned int nc);
  Opacity getH2OLinesOpacityUpTo(unsigned int nc, Length refalti);
  /** Accesor to get the integrated zenith H2O Lines Opacity for the current conditions,
   for the point nc of spectral window spwid.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Opacity getH2OLinesOpacity(unsigned int spwid, unsigned int nc){ if(!spwidAndIndexAreValid(spwid, nc)) return (double) -999.0; return getH2OLinesOpacity(v_transfertId_[spwid] + nc);}
  Opacity getH2OLinesOpacityUpTo(unsigned int spwid, unsigned int nc, Length refalti){ if(!spwidAndIndexAreValid(spwid, nc)) return (double) -999.0; return getH2OLinesOpacityUpTo(v_transfertId_[spwid] + nc, refalti);}
  /** Accesor to get the integrated zenith H2O Continuum Opacity for the current conditions,
   for a single frequency RefractiveIndexProfile object or for the point 0 of spectral
   window 0 of a multi-window RefractiveIndexProfile object.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Opacity getH2OContOpacity(){unsigned int n = 0; return getH2OContOpacity(n);}
  Opacity getH2OContOpacityUpTo(Length refalti){unsigned int n = 0; return getH2OContOpacityUpTo(n, refalti);}
  /** Accesor to get the integrated zenith H2O Continuum Opacity for the current conditions,
   for the point nc of spectral window 0.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Opacity getH2OContOpacity(unsigned int nc);
  Opacity getH2OContOpacityUpTo(unsigned int nc, Length refalti);
  /** Accesor to get the integrated zenith H2O Continuum Opacity for the current conditions,
   for the point nc of spectral window spwid.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Opacity getH2OContOpacity(unsigned int spwid, unsigned int nc){ if(!spwidAndIndexAreValid(spwid, nc)) return (double) -999.0; return getH2OContOpacity(v_transfertId_[spwid] + nc);}
  Opacity getH2OContOpacityUpTo(unsigned int spwid, unsigned int nc, Length refalti){ if(!spwidAndIndexAreValid(spwid, nc)) return (double) -999.0; return getH2OContOpacityUpTo(v_transfertId_[spwid] + nc, refalti);}
  /** Accesor to get the integrated zenith H2O Atmospheric Phase Delay (Dispersive part)
   for the current conditions, for a single frequency RefractiveIndexProfile object or
   for the point 0 of spectral window 0 of a multi-window RefractiveIndexProfile object.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Angle getDispersiveH2OPhaseDelay()
  {
    unsigned int n = 0;
    return getDispersiveH2OPhaseDelay(n);
  }
  /** Accesor to get the integrated zenith H2O Atmospheric Phase Delay (Dispersive part)
   for the current conditions, for a single frequency RefractiveIndexProfile object or
   for the point nc of spectral window 0.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Angle getDispersiveH2OPhaseDelay(unsigned int nc);
  /** Accesor to get the integrated zenith H2O Atmospheric Phase Delay (Dispersive part)
   for the current conditions, for a single frequency RefractiveIndexProfile object or
   for the point nc of spectral window spwid.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Angle getDispersiveH2OPhaseDelay(unsigned int spwid, unsigned int nc);

  /** Accesor to get the integrated zenith H2O Atmospheric Path length (Dispersive part)
   for the current conditions, for a single frequency RefractiveIndexProfile object or
   for the point 0 of spectral window 0 of a multi-window RefractiveIndexProfile object.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Length getDispersiveH2OPathLength()
  {
    unsigned int n = 0;
    return getDispersiveH2OPathLength(n);
  }
  /** Accesor to get the integrated zenith H2O Atmospheric Path length (Dispersive part)
   for the current conditions, for a single frequency RefractiveIndexProfile object or
   for the point nc of spectral window 0.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Length getDispersiveH2OPathLength(unsigned int nc);
  /** Accesor to get the integrated zenith H2O Atmospheric Path length (Dispersive part)
   for the current conditions, for a single frequency RefractiveIndexProfile object or
   for the point nc of spectral window spwid.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Length getDispersiveH2OPathLength(unsigned int spwid, unsigned int nc);

  /** Accesor to get the integrated zenith H2O Atmospheric Phase Delay (Non-Dispersive part)
   for the current conditions, for a single frequency RefractiveIndexProfile object or
   for the point 0 of spectral window 0 of a multi-window RefractiveIndexProfile object.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Angle getNonDispersiveH2OPhaseDelay()
  {
    unsigned int n = 0;
    return getNonDispersiveH2OPhaseDelay(n);
  }
  /** Accesor to get the integrated zenith H2O Atmospheric Phase Delay (Non-Dispersive part)
   for the current conditions, for a single frequency RefractiveIndexProfile object or
   for the point nc of spectral window 0.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Angle getNonDispersiveH2OPhaseDelay(unsigned int nc);
  /** Accesor to get the integrated zenith H2O Atmospheric Phase Delay (Non-Dispersive part)
   for the current conditions, for a single frequency RefractiveIndexProfile object or
   for the point nc of spectral window spwid.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Angle getNonDispersiveH2OPhaseDelay(unsigned int spwid, unsigned int nc);

  /** Accesor to get the integrated zenith H2O Atmospheric Path length (Non-Dispersive part)
   for the current conditions, for a single frequency RefractiveIndexProfile object or
   for the point 0 of spectral window 0 of a multi-window RefractiveIndexProfile object.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Length getNonDispersiveH2OPathLength()
  {
    unsigned int n = 0;
    return getNonDispersiveH2OPathLength(n);
  }
  /** Accesor to get the integrated zenith H2O Atmospheric Path length (Non-Dispersive part)
   for the current conditions, for a single frequency RefractiveIndexProfile object or
   for the point nc of spectral window 0.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Length getNonDispersiveH2OPathLength(unsigned int nc);
  /** Accesor to get the integrated zenith H2O Atmospheric Path length (Non-Dispersive part)
   for the current conditions, for a single frequency RefractiveIndexProfile object or
   for the point nc of spectral window spwid.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Length getNonDispersiveH2OPathLength(unsigned int spwid, unsigned int nc);

  /** Accessor to get the average integrated zenith Atmospheric Phase Delay (Dispersive part)
   for the current conditions in spectral Window spwid.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Angle getAverageDispersiveH2OPhaseDelay(unsigned int spwid);
  /** Accessor to get the average integrated zenith Atmospheric Phase Delay (Dispersive part)
   for the current conditions in spectral Window 0.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Angle getAverageDispersiveH2OPhaseDelay()
  {
    unsigned int n = 0;
    return getAverageDispersiveH2OPhaseDelay(n);
  }
  /** Accessor to get the average integrated zenith Atmospheric Phase Delay (Non-Dispersive
   part) in spectral Window spwid.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Angle getAverageNonDispersiveH2OPhaseDelay(unsigned int spwid);
  /** Accessor to get the average integrated zenith Atmospheric Phase Delay (Non-Dispersive
   part) in spectral Window 0.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Angle getAverageNonDispersiveH2OPhaseDelay()
  {
    unsigned int n = 0;
    return getAverageNonDispersiveH2OPhaseDelay(n);
  }

  /** Accessor to get the average integrated zenith Atmospheric Path Length (Dispersive part)
   in spectral Window spwid.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Length getAverageDispersiveH2OPathLength(unsigned int spwid);
  /** Accessor to get the average integrated zenith Atmospheric Path Length (Dispersive part)
   in spectral Window 0.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Length getAverageDispersiveH2OPathLength()
  {
    unsigned int n = 0;
    return getAverageDispersiveH2OPathLength(n);
  }
  /** Accessor to get the average integrated zenith Atmospheric Path Length (Non-Dispersive
   part) in spectral Window spwid.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Length getAverageNonDispersiveH2OPathLength(unsigned int spwid);
  /** Accessor to get the average integrated zenith Atmospheric Path Length (Non-Dispersive
   part) in spectral Window 0.
   There is overloading. The same accessor exists in RefractiveIndexProfile but in that
   case the returned value corresponds to the zenith water vapor column of the AtmProfile object.*/
  Length getAverageNonDispersiveH2OPathLength()
  {
    unsigned int n = 0;
    return getAverageNonDispersiveH2OPathLength(n);
  }

  Length getAverageH2OPathLength(unsigned int spwid)
  {
    return getAverageDispersiveH2OPathLength(spwid)
        + getAverageNonDispersiveH2OPathLength(spwid);
  }
  Length getAverageH2OPathLength()
  {
    unsigned int n = 0;
    return getAverageH2OPathLength(n);
  }
  double getAverageH2OPathLengthDerivative(unsigned int spwid)
  {
    return getAverageH2OPathLength(spwid).get(Length::UnitMicrons)
        / getUserWH2O().get(Length::UnitMicrons);
  } // in microns/micron_H2O
  double getAverageH2OPathLengthDerivative()
  {
    unsigned int n = 0;
    return getAverageH2OPathLengthDerivative(n);
  } // in microns/micron_H2O


  Angle getAverageH2OPhaseDelay(unsigned int spwid)
  {
    return getAverageDispersiveH2OPhaseDelay(spwid)
        + getAverageNonDispersiveH2OPhaseDelay(spwid);
  }
  Angle getAverageH2OPhaseDelay()
  {
    unsigned int n = 0;
    return getAverageH2OPhaseDelay(n);
  }

  Length getAverageDispersiveDryPathLength(unsigned int spwid)
  {
    return getAverageO2LinesPathLength(spwid)
        + getAverageO3LinesPathLength(spwid)
        + getAverageN2OLinesPathLength(spwid)
        + getAverageCOLinesPathLength(spwid)
        + getAverageNO2LinesPathLength(spwid)
        + getAverageSO2LinesPathLength(spwid);
  }
  Length getAverageDispersiveDryPathLength()
  {
    unsigned int n = 0;
    return getAverageDispersiveDryPathLength(n);
  }

  Length getAverageNonDispersiveDryPathLength(unsigned int spwid)
  {
    return RefractiveIndexProfile::getAverageNonDispersiveDryPathLength(spwid);
  }
  Length getAverageNonDispersiveDryPathLength()
  {
    unsigned int n = 0;
    return getAverageNonDispersiveDryPathLength(n);
  }
  double
      getAverageNonDispersiveDryPathLength_GroundPressureDerivative(unsigned int spwid);
  double
      getAverageNonDispersiveDryPathLength_GroundTemperatureDerivative(unsigned int spwid);
  double
      getAverageDispersiveDryPathLength_GroundPressureDerivative(unsigned int spwid);
  double
      getAverageDispersiveDryPathLength_GroundTemperatureDerivative(unsigned int spwid);

  Length getAverageO2LinesPathLength(unsigned int spwid)
  {
    return RefractiveIndexProfile::getAverageO2LinesPathLength(spwid);
  }
  Length getAverageO2LinesPathLength()
  {
    unsigned int n = 0;
    return getAverageO2LinesPathLength(n);
  }

  /** Setter for a new set of basic atmospheric parameters. Automatically updates the AtmProfile
   (if generated from this basic atmospheric parameters), RefractiveIndexProfile, and SkyStatus
   objects if the one or more value differ from the previous ones (there is overloading) */
  bool setBasicAtmosphericParameters(const Length &altitude,
                                     const Pressure &groundPressure,
                                     const Temperature &groundTemperature,
                                     double tropoLapseRate,
                                     const Humidity &humidity,
                                     const Length &wvScaleHeight);
  /** The same setter is available for only one parameter each time (except for "Length wvScaleHeight" due to a type
   conflict with "Length altitude") */
  bool setBasicAtmosphericParameters(const Length &altitude);
  /** The same setter is available for only one parameter each time (except for "Length wvScaleHeight" due to a type
   conflict with "Length altitude") */
  bool setBasicAtmosphericParameters(const Temperature &groundTemperature);
  /** The same setter is available for only one parameter each time (except for "Length wvScaleHeight" due to a type
   conflict with "Length altitude") */
  bool setBasicAtmosphericParameters(const Pressure &groundPressure);
  /** The same setter is available for only one parameter each time (except for "Length wvScaleHeight" due to a type
   conflict with "Length altitude") */
  bool setBasicAtmosphericParameters(const Humidity &humidity);
  /** The same setter is available for only one parameter each time (except for "Length wvScaleHeight" due to a type
   conflict with "Length altitude") */
  bool setBasicAtmosphericParameters(double tropoLapseRate);
  /** The same setter is available for couples of parameters. The order does not matter except for (Length altitude) and
   (Length wvScaleHeight) because they are of the same type (Length). The convention is that if a type Length is in first
   position, it represents the altitude, and if it is in second position, it represents wvScaleHeight */
  bool setBasicAtmosphericParameters(const Length &altitude,
                                     const Temperature &groundTemperature);
  /** The same setter is available for couples of parameters. The order does not matter except for (Length altitude) and
   (Length wvScaleHeight) because they are of the same type (Length). The convention is that if a type Length is in first
   position, it represents the altitude, and if it is in second position, it represents wvScaleHeight */
  bool setBasicAtmosphericParameters(const Length &altitude, const Pressure &groundPressure);
  /** The same setter is available for couples of parameters. The order does not matter except for (Length altitude) and
   (Length wvScaleHeight) because they are of the same type (Length). The convention is that if a type Length is in first
   position, it represents the altitude, and if it is in second position, it represents wvScaleHeight */
  bool setBasicAtmosphericParameters(const Length &altitude, const Humidity &humidity);
  /** The same setter is available for couples of parameters. The order does not matter except for (Length altitude) and
   (Length wvScaleHeight) because they are of the same type (Length). The convention is that if a type Length is in first
   position, it represents the altitude, and if it is in second position, it represents wvScaleHeight */
  bool setBasicAtmosphericParameters(const Length &altitude, double tropoLapseRate);
  /** The same setter is available for couples of parameters. The order does not matter except for (Length altitude) and
   (Length wvScaleHeight) because they are of the same type (Length). The convention is that if a type Length is in first
   position, it represents the altitude, and if it is in second position, it represents wvScaleHeight */
  bool setBasicAtmosphericParameters(const Length &altitude, const Length &wvScaleHeight);
  /** The same setter is available for couples of parameters. The order does not matter */
  bool setBasicAtmosphericParameters(const Temperature &groundTemperature,
                                     const Pressure &groundPressure);
  /** The same setter is available for couples of parameters. The order does not matter  */
  bool setBasicAtmosphericParameters(const Pressure &groundPressure,
                                     const Temperature &groundTemperature);
  /** The same setter is available for couples of parameters. The order does not matter  */
  bool setBasicAtmosphericParameters(const Temperature &groundTemperature,
                                     const Humidity &humidity);
  /** The same setter is available for couples of parameters. The order does not matter  */
  bool setBasicAtmosphericParameters(const Humidity &humidity,
                                     const Temperature &groundTemperature);
  /** The same setter is available for couples of parameters. The order does not matter  */
  bool setBasicAtmosphericParameters(const Temperature &groundTemperature,
                                     double tropoLapseRate);
  /** The same setter is available for couples of parameters. The order does not matter  */
  bool setBasicAtmosphericParameters(double tropoLapseRate,
                                     const Temperature &groundTemperature);
  /** The same setter is available for couples of parameters. The order does not matter except for (Length altitude) and
   (Length wvScaleHeight) because they are of the same type (Length). The convention is that if a type Length is in first
   position, it represents the altitude, and if it is in second position, it represents wvScaleHeight */
  bool setBasicAtmosphericParameters(const Temperature &groundTemperature,
                                     const Length &wvScaleHeight);
  /** The same setter is available for couples of parameters. The order does not matter  */
  bool
      setBasicAtmosphericParameters(const Pressure &groundPressure, const Humidity &humidity);
  /** The same setter is available for couples of parameters. The order does not matter  */
  bool
      setBasicAtmosphericParameters(const Humidity &humidity, const Pressure &groundPressure);
  /** The same setter is available for couples of parameters. The order does not matter  */
  bool setBasicAtmosphericParameters(const Pressure &groundPressure,
                                     double tropoLapseRate);
  /** The same setter is available for couples of parameters. The order does not matter */
  bool setBasicAtmosphericParameters(double tropoLapseRate,
                                     const Pressure &groundPressure);
  /** The same setter is available for couples of parameters. The order does not matter except for (Length altitude) and
   (Length wvScaleHeight) because they are of the same type (Length). The convention is that if a type Length is in first
   position, it represents the altitude, and if it is in second position, it represents wvScaleHeight */
  bool setBasicAtmosphericParameters(const Pressure &groundPressure,
                                     const Length &wvScaleHeight);
  /** The same setter is available for couples of parameters. The order does not matter  */
  bool setBasicAtmosphericParameters(const Humidity &humidity, double tropoLapseRate);
  /** The same setter is available for couples of parameters. The order does not matter  */
  bool setBasicAtmosphericParameters(double tropoLapseRate, const Humidity &humidity);
  /** The same setter is available for couples of parameters. The order does not matter except for (Length altitude) and
   (Length wvScaleHeight) because they are of the same type (Length). The convention is that if a type Length is in first
   position, it represents the altitude, and if it is in second position, it represents wvScaleHeight */
  bool setBasicAtmosphericParameters(const Humidity &humidity, const Length &wvScaleHeight);
  /** The same setter is available for couples of parameters. The order does not matter except for (Length altitude) and
   (Length wvScaleHeight) because they are of the same type (Length). The convention is that if a type Length is in first
   position, it represents the altitude, and if it is in second position, it represents wvScaleHeight */
  bool setBasicAtmosphericParameters(double tropoLapseRate,
                                     const Length &wvScaleHeight);

  /** Setter of new spectral windows (adds frequencies to the SpectralGrid and calculates the corresponding absorption coefficients) */
  void addNewSpectralWindow(unsigned int numChan,
                            unsigned int refChan,
                            const Frequency &refFreq,
                            const Frequency &chanSep)
  {
    RefractiveIndexProfile::add(numChan, refChan, refFreq, chanSep);
  }
  /** Setter of new spectral windows (adds frequencies to the SpectralGrid and calculates the corresponding absorption coefficients) */
  void addNewSpectralWindow(unsigned int numChan,
                            unsigned int refChan,
                            double* chanFreq,
                            Frequency::Units freqUnits)
  {
    RefractiveIndexProfile::add(numChan, refChan, chanFreq, freqUnits);
  }
  /** Setter of new spectral windows (adds frequencies to the SpectralGrid and calculates the corresponding absorption coefficients) */
  void addNewSpectralWindow(unsigned int numChan,
                            double refFreq,
                            double* chanFreq,
                            Frequency::Units freqUnits)
  {
    RefractiveIndexProfile::add(numChan, refFreq, chanFreq, freqUnits);
  }
  /** Setter of new spectral windows (adds frequencies to the SpectralGrid and calculates the corresponding absorption coefficients) */
  void addNewSpectralWindow(unsigned int numChan,
                            double refFreq,
                            const vector<double> &chanFreq,
                            Frequency::Units freqUnits)
  {
    RefractiveIndexProfile::add(numChan, refFreq, chanFreq, freqUnits);
  }
  /** Setter of new spectral windows (adds frequencies to the SpectralGrid and calculates the corresponding absorption coefficients) */
  void addNewSpectralWindow(unsigned int numChan,
                            unsigned int refChan,
                            const Frequency &refFreq,
                            const Frequency &chanSep,
                            const Frequency &intermediateFreq,
                            const SidebandSide &sbSide,
                            const SidebandType &sbType)
  {
    RefractiveIndexProfile::add(numChan,
                                refChan,
                                refFreq,
                                chanSep,
                                intermediateFreq,
                                sbSide,
                                sbType);
  }
  void addNewSpectralWindow(const vector<Frequency> &chanFreq)
  {
    RefractiveIndexProfile::add(chanFreq);
  }

  /** Zenith Water Vapor column retrieval based on fitting the vector of zenith atmospheric transmission given as input.
   Requirements: A) The FTS individual frequencies must be in spectral window 0, and B) the input zenith atmospheric
   transmission vector should match those frequencies (i.e. both vectors should have the same size),
   C) the air mass for the retrieval is the current one. If a different one is desired, it should be
   changed using the setAirMass setter. */
  Length WaterVaporRetrieval_fromFTS(const vector<double> &v_transmission)
  {
    unsigned int spwId = 0;
    Frequency f1(-999, Frequency::UnitGigaHertz);
    Frequency f2(-999, Frequency::UnitGigaHertz);
    return WaterVaporRetrieval_fromFTS(spwId, v_transmission, f1, f2);
  }
  /** Same as above but using for the retrieval only the measurements between frequencies f1 and f2>f1 */
  Length WaterVaporRetrieval_fromFTS(const vector<double> &v_transmission,
                                     const Frequency &f1,
                                     const Frequency &f2)
  {
    unsigned int spwId = 0;
    return WaterVaporRetrieval_fromFTS(spwId, v_transmission, f1, f2);
  }
  /** Zenith Water Vapor column retrieval based on fitting the vector of zenith atmospheric transmission given as input.
   Requirements: A) The FTS individual frequencies must be in spectral window spwId, and B) the input zenith atmospheric
   transmission vector should match those frequencies (i.e. both vectors should have the same size),
   C) the air mass for the retrieval is the current one. If a different one is desired, it should be
   changed using the setAirMass setter. */
  Length WaterVaporRetrieval_fromFTS(unsigned int spwId,
                                     const vector<double> &v_transmission)
  {
    Frequency f1(-999, Frequency::UnitGigaHertz);
    Frequency f2(-999, Frequency::UnitGigaHertz);
    return WaterVaporRetrieval_fromFTS(spwId, v_transmission, f1, f2);
  }
  /** Same as above but using for the retrieval only the measurements between frequencies f1 and f2>f1 */
  Length WaterVaporRetrieval_fromFTS(unsigned int spwId,
                                     const vector<double> &v_transmission,
                                     const Frequency &f1,
                                     const Frequency &f2);

  /** Zenith Water Vapor column retrieval based on fitting measured brightness temperatures of the atmosphere */
  Length WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                      const vector<Percent> &signalGain,
                                      const vector<vector<Temperature> > &vv_tebb,
                                      const vector<vector<double> > &spwId_filters,
                                      double airmass,
                                      const vector<double> &skycoupling,
                                      const vector<Temperature> &tspill);

  Length WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                      const vector<Percent> &signalGain,
                                      const vector<Temperature> &v_tebb,
                                      const vector<vector<double> > &spwId_filters,
                                      double airmass,
                                      const vector<double> &skycoupling,
                                      const vector<Temperature> &tspill);

  /* to implement in .cpp
  Length WaterVaporRetrieval_fromTEBB(unsigned int spwId,
                                      Percent signalGain,
                                      Temperature tebb,
                                      vector<double> spwId_filter,
                                      double airmass,
                                      double skycoupling,
                                      Temperature tspill); */


  Length WaterVaporRetrieval_fromTEBB(unsigned int spwId,
                                      const Percent &signalGain,
                                      const vector<Temperature> &v_tebb,
                                      const vector<double> &spwId_filter,
                                      double airmass,
                                      double skycoupling,
                                      const Temperature &tspill);

  Length WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                      const Percent &signalGain,
                                      const vector<Temperature> &v_tebb,
                                      const vector<vector<double> > &spwId_filters,
                                      double airmass,
                                      double skycoupling,
                                      const Temperature &tspill);

  Length WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                      const Percent &signalGain,
                                      const vector<Temperature> &v_tebb,
                                      double airmass,
                                      double skycoupling,
                                      const Temperature &tspill);  // equivalent eliminating the vectors already implemented, see below

  Length WaterVaporRetrieval_fromTEBB(unsigned int spwId,
                                      const vector<Temperature> &v_tebb,
                                      double skycoupling,
                                      const Temperature &tspill);

  Length WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
				      const vector<vector<Temperature> > &vv_tebb,
				      const vector<double> &skycoupling,
				      const vector<Temperature> &tspill);

  Length WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                      const vector<Temperature> &v_tebb,
                                      const vector<double> &skycoupling,
                                      const vector<Temperature> &tspill);

  /* to be implemented in the .cpp
  Length WaterVaporRetrieval_fromTEBB(unsigned int spwId,
                                      Temperature tebb,
                                      double skycoupling,
                                      Temperature tspill); */

  Length WaterVaporRetrieval_fromTEBB(unsigned int spwId,
                                      const vector<Temperature> &v_tebb,
                                      const vector<double> &spwId_filter,
                                      double skycoupling,
                                      const Temperature &tspill);

  Length WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                      const vector<vector<Temperature> > &vv_tebb,
                                      const vector<vector<double> > &spwId_filters,
                                      const vector<double> &skycoupling,
                                      const vector<Temperature> &tspill);

  Length WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                      const vector<Temperature> &v_tebb,
                                      const vector<vector<double> > &spwId_filters,
                                      const vector<double> &skycoupling,
                                      const vector<Temperature> &tspill);

  /* to be implemented in the .cpp
  Length WaterVaporRetrieval_fromTEBB(unsigned int spwId,
                                      Temperature tebb,
                                      vector<double> spwId_filters,
                                      double skycoupling,
                                      Temperature tspill); */

  Length WaterVaporRetrieval_fromTEBB(unsigned int spwId,
                                      const vector<Temperature> &v_tebb,
                                      double airmass,
                                      double skycoupling,
                                      const Temperature &tspill);

  Length WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                      const vector<vector<Temperature> > &vv_tebb,
                                      double airmass,
                                      const vector<double> &skycoupling,
                                      const vector<Temperature> &tspill);

  Length WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
				      const vector<Temperature> &v_tebb,
				      double airmass,
				      const vector<double> &skycoupling,
				      const vector<Temperature> &tspill);

  /* to be implemented in the .cpp
  Length WaterVaporRetrieval_fromTEBB(unsigned int spwId,
				      const Temperature &tebb,
				      double airmass,
				      double skycoupling,
				      const Temperature &tspill);  */

  Length WaterVaporRetrieval_fromTEBB(unsigned int spwId,
                                      const vector<Temperature> &v_tebb,
                                      const vector<double> &spwId_filter,
                                      double airmass,
                                      double skycoupling,
                                      const Temperature &tspill);

  Length WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                      const vector<vector<Temperature> > &vv_tebb,
                                      const vector<vector<double> > &spwId_filters,
                                      double airmass,
                                      const vector<double> &skycoupling,
                                      const vector<Temperature> &tspill);

  Length WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                      const vector<Temperature> &v_tebb,
                                      const vector<vector<double> > &spwId_filters,
                                      double airmass,
                                      const vector<double> &skycoupling,
                                      const vector<Temperature> &tspill);

  /* to be implemented in the .cpp
  Length WaterVaporRetrieval_fromTEBB(unsigned int spwId,
                                      const Temperature &v_tebb,
                                      const vector<double> &spwId_filter,
                                      double airmass,
                                      double skycoupling,
                                      const Temperature &tspill);   */

  Length WaterVaporRetrieval_fromTEBB(unsigned int spwId,
                                      const Percent &signalGain,
                                      const vector<Temperature> &v_tebb,
                                      double skycoupling,
                                      const Temperature &tspill);

  Length WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                      const vector<Percent> &signalGain,
                                      const vector<vector<Temperature> > &vv_tebb,
                                      const vector<double> &skycoupling,
                                      const vector<Temperature> &tspill);

  Length WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                      const vector<Percent> &signalGain,
                                      const vector<Temperature> &v_tebb,
                                      const vector<double> &skycoupling,
                                      const vector<Temperature> &tspill);

  /* to be implemented in the .cpp
  Length WaterVaporRetrieval_fromTEBB(unsigned int spwId,
                                      const Percent &signalGain,
                                      constr Temperature &v_tebb,
                                      double skycoupling,
                                      const Temperature &tspill); */

  Length WaterVaporRetrieval_fromTEBB(unsigned int spwId,
                                      const Percent &signalGain,
                                      const vector<Temperature> &v_tebb,
                                      const vector<double> &spwId_filter,
                                      double skycoupling,
                                      const Temperature &tspill);

  Length WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                      const vector<Percent> &signalGain,
                                      const vector<vector<Temperature> > &vv_tebb,
                                      const vector<vector<double> > &spwId_filters,
                                      const vector<double> &skycoupling,
                                      const vector<Temperature> &tspill);

  Length WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                      const vector<Percent> &signalGain,
                                      const vector<Temperature> &v_tebb,
                                      const vector<vector<double> > &spwId_filters,
                                      const vector<double> &skycoupling,
                                      const vector<Temperature> &tspill);

  /* to be implemented in the .cpp
  Length WaterVaporRetrieval_fromTEBB(unsigned int spwId,
                                      const Percent &signalGain,
                                      const Temperature &v_tebb,
                                      const vector<double> &spwId_filter,
                                      double skycoupling,
                                      const Temperature &tspill); */

  Length WaterVaporRetrieval_fromTEBB(unsigned int spwId,
                                      const Percent &signalGain,
                                      const vector<Temperature> &v_tebb,
                                      double airmass,
                                      double skycoupling,
                                      const Temperature &tspill);

  Length WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                      const vector<Percent> &signalGain,
                                      const vector<vector<Temperature> > &vv_tebb,
                                      double airmass,
                                      const vector<double> &skycoupling,
                                      const vector<Temperature> &tspill);

  Length WaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
				      const vector<Percent> &signalGain,
				      const vector<Temperature> &v_tebb,
				      double airmass,
				      const vector<double> &skycoupling,
				      const vector<Temperature> &tspill);

  // implemented 17/9/09
  Length WaterVaporRetrieval_fromTEBB(unsigned int spwId,
                                      const Percent &signalGain,
                                      const Temperature &tebb,
                                      double airmass,
                                      double skycoupling,
                                      const Temperature &tspill);

  double SkyCouplingRetrieval_fromTEBB(unsigned int spwId,
                                       const vector<Temperature> &v_tebb,
                                       double skycoupling,
                                       const Temperature &tspill);

  double SkyCouplingRetrieval_fromTEBB(unsigned int spwId,
                                       const vector<Temperature> &v_tebb,
                                       const vector<double> &spwId_filter,
                                       double skycoupling,
                                       const Temperature &tspill);

  double SkyCouplingRetrieval_fromTEBB(unsigned int spwId,
                                       const vector<Temperature> &v_tebb,
                                       double airmass,
                                       double skycoupling,
                                       const Temperature &tspill);

  double SkyCouplingRetrieval_fromTEBB(unsigned int spwId,
                                       const vector<Temperature> &v_tebb,
                                       const vector<double> &spwId_filter,
                                       double airmass,
                                       double skycoupling,
                                       const Temperature &tspill);

  /** (Re)setter of water vapor radiometer channels */
  void setWaterVaporRadiometer(const WaterVaporRadiometer &waterVaporRadiometer)
  {
    waterVaporRadiometer_ = waterVaporRadiometer;
  }
  /** Performs a sky coupling retrieval using WVR measurement sets between n and m
   (obtains the ratio between the current sky couplings of all WVR channels and the best fit that
   results from WVR measurement sets between n and m, keeping the other parameters unchanged)
   and uses the result to update the Skycoupling of all WVR channels */
  void updateSkyCoupling_fromWVR(vector<WVRMeasurement> &RadiometerData,
                                 unsigned int n,
                                 unsigned int m);
  /** Performs a sky coupling retrieval using WVR measurement sets between n and m
   (obtains the ratio between the current sky coupling of a single WVR channel (ichan) and the best fit that
   results from WVR measurement sets between n and m, keeping the other parameters unchanged)
   and uses the result to update the Skycoupling of than WVR channel */
  void updateSkyCouplingChannel_fromWVR(vector<WVRMeasurement> &RadiometerData,
					unsigned int ichan,
					unsigned int n,
					unsigned int m);
  void updateSkyCoupling_fromWVR(vector<WVRMeasurement> &RadiometerData,
                                 unsigned int n)
  {
    updateSkyCoupling_fromWVR(RadiometerData, n, n + 1);
  }
  void updateSkyCoupling_fromWVR(vector<WVRMeasurement> &RadiometerData,
                                 unsigned int n,
                                 unsigned int m,
                                 WaterVaporRadiometer &external_wvr)
  {
    updateSkyCoupling_fromWVR(RadiometerData, n, m);
    external_wvr = waterVaporRadiometer_;
  }
  void updateSkyCoupling_fromWVR(vector<WVRMeasurement> &RadiometerData,
                                 unsigned int n,
                                 WaterVaporRadiometer &external_wvr)
  {
    updateSkyCoupling_fromWVR(RadiometerData, n, n + 1);
    external_wvr = waterVaporRadiometer_;
  }
  /** Accessor to the sky coupling of channel n of the Water Vapor Radiometer*/
  double getWaterVaporRadiometerSkyCoupling(int n) const
  {
    return waterVaporRadiometer_.getSkyCoupling()[n];
  }
  /** Accessor to the sky coupling of channel n of the Water Vapor Radiometer*/
  Percent getWaterVaporRadiometerSignalGain(int n) const
  {
    return waterVaporRadiometer_.getsignalGain()[n];
  }
  /** Performs water vapor retrieval for WVR measurement sets between n and m */
  void WaterVaporRetrieval_fromWVR(vector<WVRMeasurement> &RadiometerData,
                                   unsigned int n,
                                   unsigned int m);
  void WaterVaporRetrieval_fromWVR(vector<WVRMeasurement> &RadiometerData,
                                   unsigned int n)
  {
    WaterVaporRetrieval_fromWVR(RadiometerData, n, n + 1);
  }
  /** Performs water vapor retrieval for one WVR measurement */
  void WaterVaporRetrieval_fromWVR(WVRMeasurement &RadiometerData);
  /** Accessor to get or check the water vapor radiometer channels */
  WaterVaporRadiometer getWaterVaporRadiometer() const
  {
    return waterVaporRadiometer_;
  }

  double
      sigmaSkyCouplingRetrieval_fromWVR(double par_fit,
                                        const WaterVaporRadiometer &wvr,
                                        vector<WVRMeasurement> &RadiometerData,
                                        unsigned int n,
                                        unsigned int m);
  double sigmaSkyCouplingRetrieval_fromWVR(double par_fit,
                                           const WaterVaporRadiometer &wvr,
                                           vector<WVRMeasurement> &RadiometerData,
                                           unsigned int n)
  {
    return sigmaSkyCouplingRetrieval_fromWVR(par_fit, wvr, RadiometerData, n, n
        + 1);
  }

  /*
    return the rms of the residuals for one channel ichan, as a function of the multiplicative
    factor par_fit to this channel's coupling efficiency.
   */
  double sigmaSkyCouplingChannelRetrieval_fromWVR(double par_fit,
						  const WaterVaporRadiometer &wvr,
						  vector<WVRMeasurement> &RadiometerData,
						  unsigned int ichan,
						  unsigned int n,
						  unsigned int m);
  Temperature getWVRAverageSigmaTskyFit(const vector<WVRMeasurement> &RadiometerData,
                                        unsigned int n,
                                        unsigned int m);
  Temperature getWVRAverageSigmaTskyFit(const vector<WVRMeasurement> &RadiometerData,
                                        unsigned int n)
  {
    return getWVRAverageSigmaTskyFit(RadiometerData, n, n + 1);
  }
  /*
     return the rms of fit residual for a single channel, in a measurement
  */
  Temperature getWVRSigmaChannelTskyFit(const vector<WVRMeasurement> &RadiometerData,
					unsigned int ichan,
					unsigned int n,
					unsigned int m);
  /*
    return the rms of water vapor retrieved values in a measurnment (n to m)
  */
  double getSigmaTransmissionFit(unsigned int spwId,
                                 const vector<double> &v_transmission,
                                 double airm,
                                 const Frequency &f1,
                                 const Frequency &f2);
  //double getSigmaTransmissionFit(){}
  Temperature getSigmaFit(unsigned int spwId,
                          const vector<Temperature> &v_tebbspec,
                          double skyCoupling,
                          const Temperature &Tspill)
  {
    return getSigmaFit(spwId, v_tebbspec, getUserWH2O(), skyCoupling, Tspill);
  }
  Temperature getSigmaFit(unsigned int spwId,
                          const vector<Temperature> &v_tebbspec,
                          const Length &wh2o,
                          double skyCoupling,
                          const Temperature &Tspill)
  {
    return getSigmaFit(spwId,
                       v_tebbspec,
                       wh2o,
                       getAirMass(),
                       skyCoupling,
                       Tspill);
  }
  Temperature getSigmaFit(unsigned int spwId,
                          const vector<Temperature> &v_tebbspec,
                          double airmass,
                          double skyCoupling,
                          const Temperature &Tspill)
  {
    return getSigmaFit(spwId,
                       v_tebbspec,
                       getUserWH2O(),
                       airmass,
                       skyCoupling,
                       Tspill);
  }
  Temperature getSigmaFit(unsigned int spwId,
                          const vector<Temperature> &v_tebbspec,
                          const Length &wh2o,
                          double airmass,
                          double skyCoupling,
                          const Temperature &Tspill);

  //@}

protected:

  double airMass_; //!< Air Mass used for the radiative transfer
  Temperature skyBackgroundTemperature_; //!< Blackbody temperature of the sky background
  Length wh2o_user_; //!< Water vapor column used for radiative transfer calculations. If not provided,
  //!< the one retrieved from the water vapor radiometer channels will be used.
  WaterVaporRadiometer waterVaporRadiometer_; // !< Identifiers, sky coupling, and sideband gain of channels corresponding to the water vapor radiometer.


  void iniSkyStatus(); //!< Basic Method initialize the class when using the constructors.
  bool updateProfilesAndRadiance(const Length &altitude,
                                 const Pressure &groundPressure,
                                 const Temperature &groundTemperature,
                                 double tropoLapseRate,
                                 const Humidity &relativeHumidity,
                                 const Length &wvScaleHeight);
  void rmSkyStatus(); //!< Resets retrieved water column to zero, or the default value.

  Length mkWaterVaporRetrieval_fromFTS(unsigned int spwId,
                                       const vector<double> &v_transmission,
                                       //double airm,     // unused parameter
                                       const Frequency &f1,
                                       const Frequency &f2);

  Length mkWaterVaporRetrieval_fromTEBB(unsigned int spwId,
                                        const Percent &signalGain,
                                        const vector<Temperature> &v_measuredSkyTEBB,
                                        double airm,
                                        const vector<double> &spwId_filter,
                                        double skycoupling,
                                        const Temperature &tspill)
  {
		return mkWaterVaporRetrieval_fromTEBB(
				vector<unsigned int>(1, spwId),
				vector<Percent>(1, signalGain),
				vector<vector<Temperature> >(1, v_measuredSkyTEBB),
				airm,
				vector<vector<double> >(1, spwId_filter),
				vector<double>(1, skycoupling),
				vector<Temperature>(1, tspill));
  }

  Length
      mkWaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                     const vector<Percent> &signalGain,
                                     const vector<vector<Temperature> > &measuredSkyTEBB,
                                     double airm,
                                     const vector<vector<double> > &spwId_filters,
                                     const vector<double> &skycoupling,
                                     const vector<Temperature> &tspill);

  Length
      mkWaterVaporRetrieval_fromTEBB(const vector<unsigned int> &spwId,
                                     const vector<Percent> &signalGain,
                                     const vector<Temperature> &measuredAverageSkyTEBB,
                                     double airm,
                                     const vector<vector<double> > &spwId_filters,
                                     const vector<double> &skycoupling,
                                     const vector<Temperature> &tspill);

  double mkSkyCouplingRetrieval_fromTEBB(unsigned int spwId,
                                         const Percent &signalGain,
                                         const vector<Temperature> &measuredSkyTEBB,
                                         double airm,
                                         const vector<double> &spwId_filter,
                                         double skycoupling,
                                         const Temperature &tspill);

  WVRMeasurement
      mkWaterVaporRetrieval_fromWVR(const vector<Temperature> &measuredSkyBrightnessVector,
                                    const vector<unsigned int> &radiometricChannels,
                                    const vector<double> &skyCoupling,
                                    const vector<Percent> &signalGain,
                                    const Temperature &spilloverTemperature,
                                    const Angle &elevation);

  double RT(double pfit_wh2o,
            double skycoupling,
            double tspill,
            unsigned int spwid,
            unsigned int nc)
  {
    return RT(pfit_wh2o, skycoupling, tspill, airMass_, spwid, nc);
  }
  double RT(double pfit_wh2o,
            double skycoupling,
            double tspill,
            unsigned int spwid)
  {
    return RT(pfit_wh2o, skycoupling, tspill, airMass_, spwid);
  }
  double RT(double pfit_wh2o,
            double skycoupling,
            double tspill,
            double airmass,
            unsigned int spwId,
            unsigned int nc);

  double RT(double pfit_wh2o,
            double skycoupling,
            double tspill,
            double airmass,
            unsigned int spwid)
  {
    double tebb_channel = 0.0;
    for(unsigned int n = 0; n < v_numChan_[spwid]; n++) {
      tebb_channel = tebb_channel + RT(pfit_wh2o,
                                       skycoupling,
                                       tspill,
                                       airmass,
                                       spwid,
                                       n) / (v_numChan_[spwid]);
    }
    return tebb_channel;
  }

  double RT(double pfit_wh2o,
            double skycoupling,
            double tspill,
            unsigned int spwid,
            const Percent &signalgain)
  {
    vector<double> spwId_filter;
    for(unsigned int n = 0; n < v_numChan_[spwid]; n++) {
      spwId_filter.push_back(1.0);
    }
    return RT(pfit_wh2o,
              skycoupling,
              tspill,
              airMass_,
              spwid,
              spwId_filter,
              signalgain);
  }

  double RT(double pfit_wh2o,
            double skycoupling,
            double tspill,
            double airmass,
            unsigned int spwid,
            const Percent &signalgain)
  {
    vector<double> spwId_filter;
    for(unsigned int n = 0; n < v_numChan_[spwid]; n++) {
      spwId_filter.push_back(1.0);
    }
    return RT(pfit_wh2o,
              skycoupling,
              tspill,
              airmass,
              spwid,
              spwId_filter,
              signalgain);
  }

  double RT(double pfit_wh2o,
            double skycoupling,
            double tspill,
            unsigned int spwid,
            const vector<double> &spwId_filter)
  {
    return RT(pfit_wh2o,
              skycoupling,
              tspill,
              airMass_,
              spwid,
              spwId_filter,
              Percent(100.0, Percent::UnitPercent));
  }

  double RT(double pfit_wh2o,
            double skycoupling,
            double tspill,
            unsigned int spwid,
            const vector<double> &spwId_filter,
            const Percent &signalgain)
  {
    return RT(pfit_wh2o,
              skycoupling,
              tspill,
              airMass_,
              spwid,
              spwId_filter,
              signalgain);
  }

  double RT(double pfit_wh2o,
            double skycoupling,
            double tspill,
            double airmass,
            unsigned int spwid,
            const vector<double> &spwId_filter,
            const Percent &signalgain);


  double RTRJ(double pfit_wh2o,
            double skycoupling,
            double tspill,
            unsigned int spwid,
            unsigned int nc)
  {
    return RTRJ(pfit_wh2o, skycoupling, tspill, airMass_, spwid, nc);
  }
  double RTRJ(double pfit_wh2o,
            double skycoupling,
            double tspill,
            unsigned int spwid)
  {
    return RTRJ(pfit_wh2o, skycoupling, tspill, airMass_, spwid);
  }
  double RTRJ(double pfit_wh2o,
            double skycoupling,
            double tspill,
            double airmass,
            unsigned int spwId,
            unsigned int nc);

  double RTRJ(double pfit_wh2o,
            double skycoupling,
            double tspill,
            double airmass,
            unsigned int spwid)
  {
    double trj_channel = 0.0;
    for(unsigned int n = 0; n < v_numChan_[spwid]; n++) {
      trj_channel = trj_channel + RTRJ(pfit_wh2o,
                                       skycoupling,
                                       tspill,
                                       airmass,
                                       spwid,
                                       n) / (v_numChan_[spwid]);
    }
    return trj_channel;
  }

  double RTRJ(double pfit_wh2o,
            double skycoupling,
            double tspill,
            unsigned int spwid,
            const Percent &signalgain)
  {
    vector<double> spwId_filter;
    for(unsigned int n = 0; n < v_numChan_[spwid]; n++) {
      spwId_filter.push_back(1.0);
    }
    return RTRJ(pfit_wh2o,
              skycoupling,
              tspill,
              airMass_,
              spwid,
              spwId_filter,
              signalgain);
  }

  double RTRJ(double pfit_wh2o,
            double skycoupling,
            double tspill,
            double airmass,
            unsigned int spwid,
            const Percent &signalgain)
  {
    vector<double> spwId_filter;
    for(unsigned int n = 0; n < v_numChan_[spwid]; n++) {
      spwId_filter.push_back(1.0);
    }
    return RTRJ(pfit_wh2o,
              skycoupling,
              tspill,
              airmass,
              spwid,
              spwId_filter,
              signalgain);
  }

  double RTRJ(double pfit_wh2o,
            double skycoupling,
            double tspill,
            unsigned int spwid,
            const vector<double> &spwId_filter)
  {
    return RTRJ(pfit_wh2o,
              skycoupling,
              tspill,
              airMass_,
              spwid,
              spwId_filter,
              Percent(100.0, Percent::UnitPercent));
  }

  double RTRJ(double pfit_wh2o,
            double skycoupling,
            double tspill,
            unsigned int spwid,
            const vector<double> &spwId_filter,
            const Percent &signalgain)
  {
    return RTRJ(pfit_wh2o,
              skycoupling,
              tspill,
              airMass_,
              spwid,
              spwId_filter,
              signalgain);
  }

  double RTRJ(double pfit_wh2o,
            double skycoupling,
            double tspill,
            double airmass,
            unsigned int spwid,
            const vector<double> &spwId_filter,
            const Percent &signalgain);


private:

  // no accessors provided for these
  double sigma_transmission_FTSfit_;
  Temperature sigma_TEBBfit_;

}; // class SkyStatus

ATM_NAMESPACE_END

#endif /*!_ATM_SKYSTATUS_H*/
