#ifndef _ATM_LENGTH_H
#define _ATM_LENGTH_H
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
 * "@(#) $Id: ATMLength.h Exp $"
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

/*! \brief Class for those physical parameters having dimensions of Length [L]
 *
 *   This class is defined for those physical parameters that have units of Length, for example
 *   the water vapor scale height. Default units are m (International System).
 */
class Length
{
public:
  enum Units {
    UnitKiloMeter,
    UnitMeter,
    UnitMilliMeter,
    UnitMicron,
    UnitMicrons,
    UnitNanoMeter,
    NumLengthUnits
  };

  /** Default constructor: Length value set to 0 m */
  Length();
  /** A full constructor: value in default units (m) */
  Length(double length);
  /** A full constructor: value + units. Valid units are:  km [KM], m [M], mm [MM], micron [MICRON], nm [NM]. */
  Length(double length, const std::string &units);
  Length(double length, Units units);
  /** Copy constructor */
  Length(const Length &length);

  /** Destructor */
  virtual ~Length();

  /** Accessor to get the numerical value of length (in International System units: m) */
  inline double get() const { return valueIS_; }
  /** Accessor to the length value in specified units. Implemented units are km [KM], m [M], mm [MM], micron [MICRON], nm [NM].
   *  If none of these implemented units is given, the SI value will be returned. */
  inline double get(const std::string &units) const { return sget(valueIS_, units); }
  inline double get(const Units units) const { return sget(valueIS_, units); }
  /** Accessor to the length in specified units as a formatted std::string.
   *  Implemented units are km [KM], m [M], mm [MM], micron [MICRON], nm [NM].
   *  If none of these implemented units is given, the SI value will be returned. */
  std::string get(const std::string &form, const std::string &units) const;

  /** Operator "equal to a Length" */
  inline Length& operator=(const Length &rhs) { if(&rhs != this) valueIS_ = rhs.valueIS_; return *this; }
  /** Operator "equal to a double converted to Length in m" */
  inline Length& operator=(double rhs) { valueIS_ = rhs; return *this; }
  /** Operator "addition of lengths" */
  inline Length operator+(const Length &rhs) { return Length(valueIS_ + rhs.get()); }
  /** Operator "substraction of lengths" */
  inline Length operator-(const Length &rhs) { return Length(valueIS_ - rhs.get()); }
  /** Operator "multiplication of a length by a double" */
  inline Length operator*(double scf) { return Length(valueIS_ * scf); }
  /** Operator "multiplication of a length by a float" */
  inline Length operator*(float scf) { return Length(valueIS_ * (double) scf); }
  /** Operator "multiplication of a length by an int" */
  inline Length operator*(int scf) { return Length(valueIS_ * (double) scf); }
  /** Operator "multiplication of a length by an unsigned int" */
  inline Length operator*(unsigned int scf) { return Length(valueIS_ * (double) scf); }
  /** Operator "division of a length by a double" */
  inline Length operator/(double scf) { return Length(valueIS_ / scf); }
  /** Operator "division of a length by a float" */
  inline Length operator/(float scf) { return Length(valueIS_ / (double) scf); }
  /** Operator "division of a length by an int" */
  inline Length operator/(int scf) { return Length(valueIS_ / (double) scf); }
  /** Operator "division of a length by an unsigned int" */
  inline Length operator/(unsigned int scf) { return Length(valueIS_ / (double) scf); }
  /** Operator "comparator < for two lengths" */
  inline bool operator<(const Length &rhs) const { return (valueIS_ < rhs.get()); }
  /** Operator "comparator > for two lengths" */
  inline bool operator>(const Length &rhs) const { return (valueIS_ > rhs.get()); }
  /** Operator "comparator <= for two lengths" */
  inline bool operator<=(const Length &rhs) const { return (valueIS_ <= rhs.get()); }
  /** Operator "comparator >= for two lengths" */
  inline bool operator>=(const Length &rhs) const { return (valueIS_ >= rhs.get()); }
  /** Operator "comparator == for two lengths" */
  inline bool operator==(const Length &rhs) const { return (valueIS_ == rhs.get()); }
  /** Operator "comparator != for two lengths" */
  inline bool operator!=(const Length &rhs) const { return (valueIS_ != rhs.get()); }

private:
  static double sget(double value, const std::string &units);
  static double sput(double value, const std::string &units);
  static double sget(double value, Units units);
  static double sput(double value, Units units);

private:
  double valueIS_;
}; // class Length

ATM_NAMESPACE_END

#endif /*!_ATM_LENGTH_H*/
