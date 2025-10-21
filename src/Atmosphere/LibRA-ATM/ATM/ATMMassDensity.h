#ifndef _ATM_MASSDENSITY_H
#define _ATM_MASSDENSITY_H
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
 * "@(#) $Id: ATMMassDensity.h Exp $"
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
/*! \brief Mass Density value with units
 *
 *   This class defines physical parameters having units of Mass Density.
 *   For example, the profiles of minor gases, such as O3, are in mass density.
 *   Default is m**-3 (International System).
 */
class MassDensity
{
public:
  enum Units {
    UnitGramPerCubicCentiMeter,
    UnitGramPerCubicMeter,
    UnitKiloGramPerCubicMeter,
    NumMassDensityUnits
  };

  /** Default constructor */
  MassDensity();
  /** A full constructor: Mass Density value assumed by default to be in kgm**-3 (International System) */
  MassDensity(double massdensity);
  /** A full constructor: Mass Density value + unit.
   *  Valid units are kgm**-3 [kg m**-3, KGM**-3, KG M**-3], gcm**-3 [g cm**-3, GCM**-3, G CM**-3].
   *  If none of these implemented units is given, the SI value will be returned. */
  MassDensity(double massdensity, const string &units);
  MassDensity(double massdensity, Units units);

  /** Destructor */
  virtual ~MassDensity();

  /** Accessor to the massdensity value in International System units (kgm**-3) */
  double get() const { return valueIS_; }
  /** Accessor to the massdensity value in specified units.
   *  Valid units are kgm**-3 [kg m**-3, KGM**-3, KG M**-3], gcm**-3 [g cm**-3, GCM**-3, G CM**-3]. */
  double get(const string &units) const;
  double get(const Units units) const;

  MassDensity& operator=(const MassDensity &rhs) { if(&rhs != this) valueIS_ = rhs.valueIS_; return *this; }
  MassDensity& operator=(double rhs) { valueIS_ = rhs; return *this; }
  MassDensity operator+(const MassDensity &rhs) { return MassDensity(valueIS_ + rhs.get()); }
  MassDensity operator-(const MassDensity &rhs) { return MassDensity(valueIS_ - rhs.get()); }
  MassDensity operator*(double scf) { return MassDensity(valueIS_ * scf); }
  MassDensity operator*(float scf) { return MassDensity(valueIS_ * (double) scf); }
  MassDensity operator*(int scf) { return MassDensity(valueIS_ * (double) scf); }
  MassDensity operator*(unsigned int scf) { return MassDensity(valueIS_ * (double) scf); }
  MassDensity operator/(double scf) { return MassDensity(valueIS_ / scf); }
  MassDensity operator/(float scf) { return MassDensity(valueIS_ / (double) scf); }
  MassDensity operator/(int scf) { return MassDensity(valueIS_ / (double) scf); }
  MassDensity operator/(unsigned int scf) { return MassDensity(valueIS_ / (double) scf); }
  bool operator<(const MassDensity & rhs) const { return (valueIS_ < rhs.get()); }
  bool operator>(const MassDensity & rhs) const { return (valueIS_ > rhs.get()); }
  bool operator<=(const MassDensity & rhs) const { return (valueIS_ <= rhs.get()); }
  bool operator>=(const MassDensity & rhs) const { return (valueIS_ >= rhs.get()); }
  bool operator==(const MassDensity & rhs) const { return (valueIS_ == rhs.get()); }
  bool operator!=(const MassDensity & rhs) const { return (valueIS_ != rhs.get()); }

 private:
  double valueIS_;
}; // class MassDensity

ATM_NAMESPACE_END

#endif /*!_ATM_MASSDENSITY_H*/
