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
 * "@(#) $Id: ATMNumberDensity.cpp Exp $"
 *
 * who       when      what
 * --------  --------  ----------------------------------------------
 * pardo     24/03/09  created
 */

#include "ATMNumberDensity.h"



ATM_NAMESPACE_BEGIN

NumberDensity::NumberDensity() :
  valueIS_(0.0)
{
}

NumberDensity::NumberDensity(double numberdensity) :
  valueIS_(numberdensity)
{
}

NumberDensity::NumberDensity(double numberdensity, const std::string &units)
{
  if(units == "cm**-3" || units == "CM**-3") {
    valueIS_ = 1.0E+6 * numberdensity;
  } else if(units == "m**-3" || units == "M**-3") {
    valueIS_ = numberdensity;
  } else {
    // Exception: unknown number density unit. S.I. unit (m**-3) used by default.
    valueIS_ = numberdensity;
  }
}

NumberDensity::NumberDensity(double numberdensity, NumberDensity::Units units)
{
  if(units == NumberDensity::UnitInverseCubicCentiMeter) {
    valueIS_ = 1.0E+6 * numberdensity;
  } else if(units == NumberDensity::UnitInverseCubicMeter) {
    valueIS_ = numberdensity;
  } else {
    // Exception: unknown number density unit. S.I. unit (m**-3) used by default.
    valueIS_ = numberdensity;
  }
}

NumberDensity::~NumberDensity()
{
}

double NumberDensity::get(const std::string &units) const
{
  if(units == "cm**-3" || units == "CM**-3") {
    return 1.0E-6 * valueIS_;
  } else if(units == "m**-3" || units == "M**-3") {
    return valueIS_;
  } else {
    // Exception: unknown number density unit. S.I. unit (m**-3) used by default.
    return valueIS_;
  }
}

double NumberDensity::get(NumberDensity::Units units) const
{
  if(units == NumberDensity::UnitInverseCubicCentiMeter) {
    return 1.0E-6 * valueIS_;
  } else if(units == NumberDensity::UnitInverseCubicMeter) {
    return valueIS_;
  } else {
    // Exception: unknown number density unit. S.I. unit (m**-3) used by default.
    return valueIS_;
  }
}

ATM_NAMESPACE_END
