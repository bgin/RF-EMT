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
 * "@(#) $Id: ATMAngle.cpp Exp $"
 *
 * who       when      what
 * --------  --------  ----------------------------------------------
 * pardo     24/03/09  created
 */

#include "ATMAngle.h"



ATM_NAMESPACE_BEGIN

Angle::Angle() :
  valueIS_(0.0)
{
}

Angle::Angle(double angle) :
  valueIS_(angle)
{
}

Angle::Angle(double angle, const std::string &units)
{
  if(units == "Rad" || units == "RAD" || units == "rad") {
    valueIS_ = angle;
  } else if(units == "deg" || units == "DEG") {
    valueIS_ = (angle / 360.0) * 6.2831852;
  } else {
    valueIS_ = angle;
  }
}

Angle::Angle(double angle, Angle::Units units)
{
  if(units == Angle::UnitRadian) {
    valueIS_ = angle;
  } else if(units == Angle::UnitDegree) {
    valueIS_ = (angle / 360.0) * 6.2831852;
  } else {
    valueIS_ = angle;
  }
}

Angle::~Angle()
{
}

double Angle::get(const std::string &units) const
{
  if(units == "Rad" || units == "RAD" || units == "rad") {
    return valueIS_;
  } else if(units == "deg" || units == "DEG") {
    return 360.0 * (valueIS_ / 6.2831852);
  } else {
    return valueIS_;
  }
}

double Angle::get(Angle::Units units) const
{
  if(units == Angle::UnitRadian) {
    return valueIS_;
  } else if(units == Angle::UnitDegree) {
    return 360.0 * (valueIS_ / 6.2831852);
  } else {
    return valueIS_;
  }
}

ATM_NAMESPACE_END
