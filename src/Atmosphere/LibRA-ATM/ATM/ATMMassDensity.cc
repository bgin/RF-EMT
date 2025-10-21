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
 * "@(#) $Id: ATMMassDensity.cpp Exp $"
 *
 * who       when      what
 * --------  --------  ----------------------------------------------
 * pardo     24/03/09  created
 */

#include "ATMMassDensity.h"

ATM_NAMESPACE_BEGIN

MassDensity::MassDensity() :
  valueIS_(0.0)
{
}

MassDensity::MassDensity(double massdensity) :
  valueIS_(massdensity)
{
}

MassDensity::MassDensity(double massdensity, const std::string &units)
{
  if(units == "gcm**-3" || units == "g cm**-3" || units == "GCM**-3" ||
     units == "G CM**-3" || units == "g/cm^3") {
    valueIS_ = 1.0E+3 * massdensity;
  } else if(units == "gm**-3" || units == "g m**-3" || units == "GM**-3" ||
      units == "G M**-3" || units == "g/m^3") {
    valueIS_ = 1.0E-3 * massdensity;
  } else if(units == "kgm**-3" || units == "kg m**-3" || units == "KGM**-3" ||
      units == "KG M**-3" || units == "kg/m^3") {
    valueIS_ = massdensity;
  } else {
    // Exception: unknown number density unit. S.I. unit (kg m**-3) used by default.
    valueIS_ = massdensity;
  }
}

MassDensity::MassDensity(double massdensity, MassDensity::Units units)
{
  if(units == MassDensity::UnitGramPerCubicCentiMeter) {
    valueIS_ = 1.0E+3 * massdensity;
  } else if(units == MassDensity::UnitGramPerCubicMeter) {
    valueIS_ = 1.0E-3 * massdensity;
  } else if(units == MassDensity::UnitKiloGramPerCubicMeter) {
    valueIS_ = massdensity;
  } else {
    // Exception: unknown number density unit. S.I. unit (kg m**-3) used by default.
    valueIS_ = massdensity;
  }
}

MassDensity::~MassDensity() { };


double MassDensity::get(const std::string &units) const
{
  if(units == "gcm**-3" || units == "g cm**-3" || units == "GCM**-3" ||
      units == "G CM**-3" || units == "g/cm^3") {
    return 1.0E-3 * valueIS_;
  } else if(units == "gm**-3" || units == "g m**-3" || units == "GM**-3" ||
      units == "G M**-3" || units == "g/m^3") {
    return 1.0E+3 * valueIS_;
  } else if(units == "kgm**-3" || units == "kg m**-3" || units == "KGM**-3" ||
      units == "KG M**-3" || units == "kg/m^3") {
    return valueIS_;
  } else {
    // Exception: unknown number density unit. S.I. unit (kg m**-3) used by default.
    return valueIS_;
  }
}

double MassDensity::get(MassDensity::Units units) const
{
  if(units == MassDensity::UnitGramPerCubicCentiMeter) {
    return 1.0E-3 * valueIS_;
  } else if(units == MassDensity::UnitGramPerCubicMeter) {
    return 1.0E+3 * valueIS_;
  } else if(units == MassDensity::UnitKiloGramPerCubicMeter) {
    return valueIS_;
  } else {
    // Exception: unknown number density unit. S.I. unit (kg m**-3) used by default.
    return valueIS_;
  }
}

ATM_NAMESPACE_END
