#ifndef _ATM_INVERSELENGTH_H
#define _ATM_INVERSELENGTH_H
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
* "@(#) $Id: ATMInverseLength.h Exp $"
*
* who       when      what
* --------  --------  ----------------------------------------------
* pardo     24/03/09  created
*/

#ifndef __cplusplus
#error "This is a C++ include file and cannot be used from plain C"
#endif

#include "ATMCommon.h"
#include <string>



ATM_NAMESPACE_BEGIN
/*! \brief Class for those physical parameters having dimensions of Inverse Length [L^-1]
 *
 *   This class is defined for those physical parameters that have units of Inverse Length, for example
 *   the absorption coefficients. Desfault units are m^-1 (International System).
 */
class InverseLength
{
public:
  enum Units {
    UnitInverseMeter,
    UnitInverseKiloMeter,
    UnitInverseMilliMeter,
    UnitInverseMicron,
    UnitInverseNanoMeter,
    NumInverseLengthUnits
  };

  /** Default constructor: Length value set to 0 m^-1 */
  InverseLength();
  /** A full constructor: value in default units (m^-1) */
  InverseLength(double inverseLength);
  /** A full constructor: value + units. Valid units are:  km-1 [KM-1], m-1 [M-1], mm-1 [MM-1], micron-1 [MICRON-1], nm-1 [NM-1]. */
  InverseLength(double inverseLength, const std::string &units);
  InverseLength(double inverseLength, Units units);
  /** Copy constructor */
  InverseLength (const InverseLength &inverseLength);

  /** Destructor */
  virtual ~InverseLength();

  /** Accessor to get the numerical value of inverse length (in International System units: m^-1) */
  inline double get() const { return valueIS_; }
  /** Accessor to the inverse length value in specified units. Implemented units are km-1 [KM-1], m-1 [M-1], mm-1 [MM-1], micron-1 [MICRON-1], nm-1 [NM-1].
   *  If none of these implemented units is given, the SI value will be returned. */
  inline double get(const std::string &units) const { return sget(valueIS_, units); }
  inline double get(Units units) const { return sget(valueIS_, units); }

  /** Operator "equal to a InverseLength" */
  inline InverseLength& operator=(const InverseLength &rhs) { if(&rhs != this) valueIS_ = rhs.valueIS_; return *this; }
  /** Operator "equal to a double converted to InverseLength in m-1" */
  inline InverseLength& operator=(double rhs) { valueIS_ = rhs; return *this; }
  /** Operator "addition of inverse lengths" */
  inline InverseLength operator+(const InverseLength &rhs) { return InverseLength(valueIS_ + rhs.get()); }
  /** Operator "substraction of inverse lengths" */
  inline InverseLength operator-(const InverseLength &rhs) { return InverseLength(valueIS_ - rhs.get()); }
  /** Operator "multiplication of a inverse length by a double" */
  inline InverseLength operator*(double scf) { return InverseLength(valueIS_ * scf); }
  /** Operator "multiplication of a inverse length by a float" */
  inline InverseLength operator*(float scf) { return InverseLength(valueIS_ * (double) scf); }
  /** Operator "multiplication of a inverse length by an int" */
  inline InverseLength operator*(int scf) { return InverseLength(valueIS_ * (double) scf); }
  /** Operator "multiplication of a inverse length by an unsigned int" */
  inline InverseLength operator*(unsigned int scf) { return InverseLength(valueIS_ * (double) scf); }
  /** Operator "division of a inverse length by a double" */
  inline InverseLength operator/(double scf) { return InverseLength(valueIS_ / scf); }
  /** Operator "division of a inverse length by a float" */
  inline InverseLength operator/(float scf) { return InverseLength(valueIS_ / (double) scf); }
  /** Operator "division of a inverse length by an int" */
  inline InverseLength operator/(int scf) { return InverseLength(valueIS_ / (double) scf); }
  /** Operator "division of a inverse length by an unsigned int" */
  inline InverseLength operator/(unsigned int scf) { return InverseLength(valueIS_ / (double) scf); }
  /** Operator "comparator < for two inverse lengths" */
  inline bool operator<(const InverseLength & rhs) const { return (valueIS_ < rhs.get()); }
  /** Operator "comparator > for two inverse lengths" */
  inline bool operator>(const InverseLength & rhs) const { return (valueIS_ > rhs.get()); }
  /** Operator "comparator <= for two inverse lengths" */
  inline bool operator<=(const InverseLength & rhs) const { return (valueIS_ <= rhs.get()); }
  /** Operator "comparator >= for two inverse lengths" */
  inline bool operator>=(const InverseLength & rhs) const { return (valueIS_ >= rhs.get()); }
  /** Operator "comparator == for two inverse lengths" */
  inline bool operator==(const InverseLength & rhs) const { return (valueIS_ == rhs.get()); }
  /** Operator "comparator != for two inverse lengths" */
  inline bool operator!=(const InverseLength & rhs) const { return (valueIS_ != rhs.get()); }

private:
  static double sget(double value, const std::string &units);
  static double sput(double value, const std::string &units);
  static double sget(double value, Units units);
  static double sput(double value, Units units);

private:
  double valueIS_;
}; // class InverseLength

ATM_NAMESPACE_END

#endif /*!_ATM_INVERSELENGTH_H*/
