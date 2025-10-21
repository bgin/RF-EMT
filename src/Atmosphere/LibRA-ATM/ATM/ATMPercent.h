#ifndef _ATM_PERCENT_H
#define _ATM_PERCENT_H
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
* "@(#) $Id: ATMPercent.h Exp $"
*
* who       when      what
* --------  --------  ----------------------------------------------
* pardo     2009-03-24  created
*/

#ifndef __cplusplus
#error "This is a C++ include file and cannot be used from plain C"
#endif

#include "ATMCommon.h"
#include <string>



ATM_NAMESPACE_BEGIN

/*! \brief Defines a class for quantities that may be expressed in percent.
 *
 *   This class is aimed at quantities that may be defined in percent. Default is proposed with no units
 *   (relative to one), and an alternative is considered for units in %.
 */
class Percent
{
public:

  enum Units {
    UnitPercent,
    NumPercentUnits
  };

  /** Default constructor */
  Percent();
  /** A full constructor: Amount in default units (relative to 1) */
  Percent(double percent);
  /** A full constructor: Amount + units. Valid units are: % */
  Percent(double percent, const std::string &units);

  Percent(double percent, Units units);

  ~Percent();

  //@{
  /** Accessor to the percent value in International System units */
  double get()const;
  /** Accessor to the percent value in specified units.  */
   double get(const std::string &units)const;
   double get(Units units) const;
  //@}

  /** Operator "equal to a Percent" */
  inline Percent& operator=(const Percent &rhs) { if(&rhs != this) valueIS_ = rhs.valueIS_; return *this; }
  /** Operator "equal to a double converted to Percent in %" */
  inline Percent& operator=(double rhs) { valueIS_ = rhs; return *this; }
  /** Operator "addition of percentages" */
  inline Percent operator+(const Percent &rhs) { return Percent(valueIS_ + rhs.get()); }
  /** Operator "substraction of percentages" */
  inline Percent operator-(const Percent &rhs) { return Percent(valueIS_ - rhs.get()); }
  /** Operator "multiplication of a percent by a double" */
  inline Percent operator*(double scf) { return Percent(valueIS_ * scf); }
  /** Operator "multiplication of a percent by a float" */
  inline Percent operator*(float scf) { return Percent(valueIS_ * (double) scf); }
  /** Operator "multiplication of a percent by an int" */
  inline Percent operator*(int scf) { return Percent(valueIS_ * (double) scf); }
  /** Operator "multiplication of a percent by an unsigned int" */
  inline Percent operator*(unsigned int scf) { return Percent(valueIS_ * (double) scf); }
  /** Operator "division of a percent by a double" */
  inline Percent operator/(double scf) { return Percent(valueIS_ / scf); }
  /** Operator "division of a percent by a float" */
  inline Percent operator/(float scf) { return Percent(valueIS_ / (double) scf); }
  /** Operator "division of a percent by an int" */
  inline Percent operator/(int scf) { return Percent(valueIS_ / (double) scf); }
  /** Operator "division of a percent by an unsigned int" */
  inline Percent operator/(unsigned int scf) { return Percent(valueIS_ / (double) scf); }
  /** Operator "comparator < for two percentages" */
  inline bool operator<(const Percent &rhs) const { return (valueIS_ < rhs.get()); }
  /** Operator "comparator > for two percentages" */
  inline bool operator>(const Percent &rhs) const { return (valueIS_ > rhs.get()); }
  /** Operator "comparator <= for two percentages" */
  inline bool operator<=(const Percent &rhs) const { return (valueIS_ <= rhs.get()); }
  /** Operator "comparator >= for two percentages" */
  inline bool operator>=(const Percent &rhs) const { return (valueIS_ >= rhs.get()); }
  /** Operator "comparator == for two percentages" */
  inline bool operator==(const Percent &rhs) const { return (valueIS_ == rhs.get()); }
  /** Operator "comparator != for two percentages" */
  inline bool operator!=(const Percent &rhs) const { return (valueIS_ != rhs.get()); }

private:
  double valueIS_;
}; // class Percent

ATM_NAMESPACE_END

#endif /*!_ATM_PERCENT_H*/


