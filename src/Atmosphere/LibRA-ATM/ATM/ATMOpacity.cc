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
 * "@(#) $Id: ATMOpacity.cpp Exp $"
 *
 * who       when      what
 * --------  --------  ----------------------------------------------
 * pardo     24/03/09  created
 */

#include <stdio.h>

#include "ATMOpacity.h"



ATM_NAMESPACE_BEGIN

Opacity::Opacity() :
  valueIS_(0.0)
{
}

Opacity::Opacity(double opacity) :
  valueIS_(opacity)
{
}

Opacity::Opacity(double opacity, const std::string &units)
{
  valueIS_ = sput(opacity, units);
}

Opacity::Opacity(double opacity, Opacity::Units units)
{
  valueIS_ = sput(opacity, units);
}

Opacity::Opacity(const Opacity &opacity) :
  valueIS_(opacity.valueIS_)
{
}

Opacity::~Opacity()
{
}

double Opacity::sget(double value, const std::string &units)
{
  if(units == "db" || units == "DB") {
    return value * 4.34294482;
  } else if(units == "np" || units == "NP" || units == "neper" || units == "NEPER"){
    return value;
  } else {
    // Exception: Unknown unit, neper (np) used by default)
    return value;
  }
}

double Opacity::sput(double value, const std::string &units)
{
  if(units == "db" || units == "DB") {
    return value / 4.34294482;
  } else if(units == "np" || units == "NP" || units == "neper" || units == "NEPER") {
    return value;
  } else {
    // Exception: Unknown unit, neper (np) used by default)
    return value;
  }
}

double Opacity::sget(double value, Opacity::Units units)
{
  if(units == Opacity::UnitDeciBel) {
    return value * 4.34294482;
  } else if(units == Opacity::UnitNeper){
    return value;
  } else {
    // Exception: Unknown unit, neper (np) used by default)
    return value;
  }
}

double Opacity::sput(double value, Opacity::Units units)
{
  if(units == Opacity::UnitDeciBel) {
    return value / 4.34294482;
  } else if(units == Opacity::UnitNeper) {
    return value;
  } else {
    // Exception: Unknown unit, neper (np) used by default)
    return value;
  }
}

ATM_NAMESPACE_END

