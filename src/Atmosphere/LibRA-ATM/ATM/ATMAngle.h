#ifndef _ATM_ANGLE_H
#define _ATM_ANGLE_H
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
 * "@(#) $Id: ATMAngle.h Exp $"
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

using std::string;

ATM_NAMESPACE_BEGIN

/*! \brief Defines a Class for those parameters being angles.
 *
 *   Defines a Class for those physical parameters having units of angle. For
 *   example the Total Phase Delay can be expressed as an angle (but also as a length).
 */
class Angle
{
 public:
  enum Units {
    UnitRadian,
    UnitDegree,
    NumAngleUnits
  };

  /** Default constructor */
  Angle();
  /** A full constructor: Angle value in default units (SI: radians) */
  Angle(double angle);
  /** A full constructor: Angle value + units. Implemented units are: deg [DEG], rad [RAD] [Rad]*/
  Angle(double angle, const string &units);
  Angle(double angle, Units units);

  /** Destructor */
  virtual ~Angle();

  /** Accessor to get the angle value in SI units (rad) */
  double get() const { return valueIS_; }
  /** Accessor to the angle value in specified units. Implemented units are: deg [DEG], rad [RAD] [Rad].
   *  If none of these implemented units is given, the SI value will be returned. */
  double get(const string &units) const;
  double get(Units units) const;

  /** Operator "equal to a Angle" */
  Angle& operator=(const Angle &rhs) { if(&rhs != this) valueIS_ = rhs.valueIS_; return *this; }
  /** Operator "equal to a double converted to Angle in Hz" */
  Angle& operator=(double rhs) { valueIS_=rhs; return *this; }
  /** Operator "addition of angles" */
  Angle operator+(const Angle &rhs) { return Angle(valueIS_+rhs.get()); }
  /** Operator "substraction of angles" */
  Angle operator-(const Angle &rhs) { return Angle(valueIS_-rhs.get()); }
  /** Operator "multiplication of a angle by a double" */
  Angle operator*(double scf) { return Angle(valueIS_*scf); }
  /** Operator "multiplication of a angle by a float" */
  Angle operator*(float scf) { return Angle(valueIS_*(double)scf); }
  /** Operator "multiplication of a angle by an int" */
  Angle operator*(int scf) { return Angle(valueIS_*(double)scf); } // rhs scale factor
  /** Operator "multiplication of a angle by an unsigned int" */
  Angle operator*(unsigned int scf) { return Angle(valueIS_*(double)scf); } // rhs scale factor
  /** Operator "division of a angle by a double" */
  Angle operator/(double scf) { return Angle(valueIS_/scf); }
  /** Operator "division of a angle by a float" */
  Angle operator/(float scf) { return Angle(valueIS_/(double)scf); }
  /** Operator "division of a angle by an int" */
  Angle operator/(int scf) { return Angle(valueIS_/(double)scf); }
  /** Operator "division of a angle by an unsigned int" */
  Angle operator/(unsigned int scf) { return Angle(valueIS_/(double)scf); }
  /** Operator "comparator < for two angles" */
  bool operator<(const Angle & rhs) const  { return (valueIS_<rhs.get()); }
  /** Operator "comparator > for two angles" */
  bool operator>(const Angle & rhs) const  { return (valueIS_>rhs.get()); }
  /** Operator "comparator <= for two angles" */
  bool operator<=(const Angle & rhs) const  { return (valueIS_<=rhs.get()); }
  /** Operator "comparator >= for two angles" */
  bool operator>=(const Angle & rhs) const  { return (valueIS_>=rhs.get()); }
  /** Operator "comparator == for two angles" */
  bool operator==(const Angle & rhs) const  { return (valueIS_==rhs.get()); }
  /** Operator "comparator != for two angles" */
  bool operator!=(const Angle & rhs) const  { return (valueIS_!=rhs.get()); }

 private:
  double valueIS_;
}; // class Angle

ATM_NAMESPACE_END

#endif /*!_ATM_ANGLE_H*/
