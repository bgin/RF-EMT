#ifndef _ATM_TEMPERATURE_H
#define _ATM_TEMPERATURE_H
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
 * "@(#) $Id: ATMTemperature.h Exp $"
 *
 * who       when      what
 * --------  --------  ----------------------------------------------
 * pardo     24/03/09  created
 */

#ifndef __cplusplus
#error This is a C++ include file and cannot be used from plain C
#endif

#include "ATMCommon.h"
#include <string>

using std::string;

ATM_NAMESPACE_BEGIN


/*! \brief Temperature value with units
 *
 *   This class defines physical parameters having units of Temperature.
 *   For example, the Radiance can be converted into Brightness Temperature.
 *   Default is K (International System), also used by the ATM library.
 */
class Temperature
{
public:
  enum Units {
    UnitKelvin,
    UnitMilliKelvin,
    UnitCelsius,
    UnitFahrenheit,
    NumTemperatureUnits
  };

  /** Default constructor */
  Temperature();
  /** A full constructor: Temperature value assumed by default to be in K (International System) */
  Temperature(double temperature);
  /** A full constructor: Temperature value + unit. Valid units are K [k], mK [mk], and C [c].
   *  If none of these implemented units is given, the SI value will be returned. */
  Temperature(double temperature, const string &units);
  Temperature(double temperature, Units units);

  /** Destructor */
  virtual ~Temperature() {};

  //@{
  /** Accessor to the temperature value in International System units (K) */
  double get() const { return valueIS_; }
  /** Accessor to the temperature value in specified units. Valid units are K [k], mK [mk], and C [c] */
  double get(const string &units) const;
  double get(Units units) const;
  //@}

  Temperature& operator=(const Temperature &rhs){ if(&rhs != this) valueIS_ = rhs.valueIS_; return *this; }
  Temperature& operator=(double rhs) { valueIS_ = rhs; return *this; }
  Temperature operator+(const Temperature &rhs) { return Temperature(valueIS_ + rhs.get()); }
  Temperature operator-(const Temperature &rhs) { return Temperature(valueIS_ - rhs.get()); }
  Temperature operator*(double scf) { return Temperature(valueIS_ * scf); }
  Temperature operator*(float scf) { return Temperature(valueIS_ * (double) scf); }
  Temperature operator*(int scf) { return Temperature(valueIS_ * (double) scf); }
  Temperature operator*(unsigned int scf) { return Temperature(valueIS_ * (double) scf); }
  Temperature operator/(double scf) { return Temperature(valueIS_ / scf); }
  Temperature operator/(float scf) { return Temperature(valueIS_ / (double) scf); }
  Temperature operator/(int scf) { return Temperature(valueIS_ / (double) scf); }
  Temperature operator/(unsigned int scf) { return Temperature(valueIS_ / (double) scf); }
  bool operator<(const Temperature &rhs) const { return (valueIS_ < rhs.get()); }
  bool operator>(const Temperature &rhs) const { return (valueIS_ > rhs.get()); }
  bool operator<=(const Temperature &rhs) const { return (valueIS_ <= rhs.get()); }
  bool operator>=(const Temperature &rhs) const { return (valueIS_ >= rhs.get()); }
  bool operator==(const Temperature &rhs) const { return (valueIS_ == rhs.get()); }
  bool operator!=(const Temperature &rhs) const { return (valueIS_ != rhs.get()); }

private:
  double valueIS_;
}; // class Temperature

ATM_NAMESPACE_END

#endif /*!_ATM_TEMPERATURE_H*/


