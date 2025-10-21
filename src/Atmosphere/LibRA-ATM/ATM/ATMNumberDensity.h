#ifndef _ATM_NUMBERDENSITY_H
#define _ATM_NUMBERDENSITY_H
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
 * "@(#) $Id: ATMNumberDensity.h Exp $"
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

/*! \brief Number Density value with units
 *
 *   This class defines physical parameters having units of Number Density.
 *   For example, the profiles of minor gases, such as O3, are in number density.
 *   Default is m**-3 (International System).
 */
class NumberDensity
{
 public:
  enum Units {
    UnitInverseCubicCentiMeter,
    UnitInverseCubicMeter,
    NumNumberDensityUnits
  };

  /** Default constructor */
  NumberDensity();
  /** A full constructor: Number Density value assumed by default to be in m**-3 (International System) */
  NumberDensity(double numberdensity);
  /** A full constructor: Number Density value + unit. Valid units are m**-3 [M**-3], cm**-3 [CM**-3].
   *  If none of these implented units is given, the SI value will be returned. */
  NumberDensity(double numberdensity, const string &units);
  NumberDensity(double numberdensity, Units units);

  /** Destructor */
  virtual ~NumberDensity();

  /** Accessor to the numberdensity value in International System units (K) */
  double get() const { return valueIS_; }
  /** Accessor to the numberdensity value in specified units. Valid units are K [k], mK [mk], and C [c] */
  double get(const string &units) const;
  double get(const Units units) const;

  NumberDensity& operator=(const NumberDensity &rhs) { if(&rhs != this) valueIS_ = rhs.valueIS_; return *this; }
  NumberDensity& operator=(const double &rhs) { valueIS_ = rhs; return *this; }
  NumberDensity operator+(const NumberDensity &rhs) { return NumberDensity(valueIS_ + rhs.get()); }
  NumberDensity operator-(const NumberDensity &rhs) { return NumberDensity(valueIS_ - rhs.get()); }
  NumberDensity operator*(double scf) { return NumberDensity(valueIS_ * scf); }
  NumberDensity operator*(float scf) { return NumberDensity(valueIS_ * (double) scf); }
  NumberDensity operator*(int scf) { return NumberDensity(valueIS_ * (double) scf); }
  NumberDensity operator*(unsigned int scf) { return NumberDensity(valueIS_ * (double) scf); }
  NumberDensity operator/(double scf) { return NumberDensity(valueIS_ / scf); }
  NumberDensity operator/(float scf) { return NumberDensity(valueIS_ / (double) scf); }
  NumberDensity operator/(int scf) { return NumberDensity(valueIS_ / (double) scf); }
  NumberDensity operator/(unsigned int scf) { return NumberDensity(valueIS_ / (double) scf); }
  bool operator<(const NumberDensity &rhs) const { return (valueIS_ < rhs.get()); }
  bool operator>(const NumberDensity &rhs) const { return (valueIS_ > rhs.get()); }
  bool operator<=(const NumberDensity &rhs) const { return (valueIS_ <= rhs.get()); }
  bool operator>=(const NumberDensity &rhs) const { return (valueIS_ >= rhs.get()); }
  bool operator==(const NumberDensity &rhs) const { return (valueIS_ == rhs.get()); }
  bool operator!=(const NumberDensity &rhs) const { return (valueIS_ != rhs.get()); }

 private:
  double valueIS_;
}; // class NumberDensity

ATM_NAMESPACE_END

#endif /*!_ATM_NUMBERDENSITY_H*/


