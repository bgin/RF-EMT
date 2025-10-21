#ifndef _ATM_OPACITY_H
#define _ATM_OPACITY_H
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
 * "@(#) $Id: ATMOpacity.h Exp $"
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

/*! \brief Class for opacities [no dimensions]
 *
 *   This class is defined for opacities.
 *   Default units are np (Transmission=exp(-opacity)).
 */
class Opacity
{
public:
  enum Units {
    UnitDeciBel,
    UnitNeper,
    NumOpacityUnits
  };

  /** Default constructor: Opacity value set to 0 np */
  Opacity();
  /** A full constructor: value in default units (np) */
  Opacity(double opacity);
  /** A full constructor: value + units. Valid units are:  np [neper] [NP] [NEPER], db [DB]. */
  Opacity(double opacity, const std::string &units);
  Opacity(double opacity, Units units);
  /** Copy constructor */
  Opacity (const Opacity &opacity);

  /** Destructor */
  virtual ~Opacity();

  /** Accessor to get the numerical value of opacity (in np) */
  double get() const { return valueIS_; }
  /** Accessor to the opacity value in specified units. Implemented units are np [neper] [NP] [NEPER], db [DB].
   *  If none of these implemented units is given, the value in neper will be returned. */
  inline double get(const std::string &units) const { return sget(valueIS_, units); }
  inline double get(Units units) const { return sget(valueIS_, units); }

  /** Operator "equal to a Opacity" */
  inline Opacity& operator=(const Opacity &rhs) { if(&rhs != this) valueIS_ = rhs.valueIS_; return *this; }
  /** Operator "equal to a double converted to Opacity in m" */
  inline Opacity& operator=(const double &rhs) { valueIS_ = rhs; return *this; }
  /** Operator "addition of opacities" */
  inline Opacity operator+(const Opacity &rhs) { return Opacity(valueIS_ + rhs.get()); }
  /** Operator "substraction of opacities" */
  inline Opacity operator-(const Opacity &rhs) { return Opacity(valueIS_ - rhs.get()); }
  /** Operator "multiplication of an opacity by a double" */
  inline Opacity operator*(double scf) { return Opacity(valueIS_ * scf); }
  /** Operator "multiplication of an opacity by a float" */
  inline Opacity operator*(float scf) { return Opacity(valueIS_ * (double) scf); }
  /** Operator "multiplication of an opacity by an int" */
  inline Opacity operator*(int scf) { return Opacity(valueIS_ * (double) scf); }
  /** Operator "multiplication of an opacity by an unsigned int" */
  inline Opacity operator*(unsigned int scf) { return Opacity(valueIS_ * (double) scf); }
  /** Operator "division of a opacity by an int" */
  inline Opacity operator/(double scf) { return Opacity(valueIS_ / scf); }
  /** Operator "division of a opacity by a float" */
  inline Opacity operator/(float scf) { return Opacity(valueIS_ / (double) scf); }
  /** Operator "division of a opacity by an int" */
  inline Opacity operator/(int scf) { return Opacity(valueIS_ / (double) scf); }
  /** Operator "division of a opacity by an unsigned int" */
  inline Opacity operator/(unsigned int scf) { return Opacity(valueIS_ / (double) scf); }
  /** Operator "comparator < for two opacities" */
  inline bool operator<(const Opacity &rhs) const { return (valueIS_ < rhs.get()); }
  /** Operator "comparator > for two opacities" */
  inline bool operator>(const Opacity &rhs) const { return (valueIS_ > rhs.get()); }
  /** Operator "comparator <= for two opacities" */
  inline bool operator<=(const Opacity &rhs) const { return (valueIS_ <= rhs.get()); }
  /** Operator "comparator >= for two opacities" */
  inline bool operator>=(const Opacity &rhs) const { return (valueIS_ >= rhs.get()); }
  /** Operator "comparator == for two opacities" */
  inline bool operator==(const Opacity &rhs) const { return (valueIS_ == rhs.get()); }
  /** Operator "comparator != for two opacities" */
  inline bool operator!=(const Opacity &rhs) const { return (valueIS_ != rhs.get()); }

private:
  static double sget(double value, const std::string &units);
  static double sput(double value, const std::string &units);
  static double sget(double value, Units units);
  static double sput(double value, Units units);

private:
  double valueIS_;
}; // class Opacity

ATM_NAMESPACE_END

#endif /*!_ATM_OPACITY_H*/
