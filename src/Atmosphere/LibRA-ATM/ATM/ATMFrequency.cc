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
 * "@(#) $Id: ATMFrequency.cpp Exp $"
 *
 * who       when      what
 * --------  --------  ----------------------------------------------
 * pardo     24/03/09  created
 */

#include <stdio.h>

#include "ATMFrequency.h"



ATM_NAMESPACE_BEGIN

Frequency::Frequency() :
  valueIS_(0.0)
{
}

Frequency::Frequency(double frequency) :
  valueIS_(frequency)
{
}

Frequency::Frequency(double frequency, const std::string &units)
{
  valueIS_ = sput(frequency, units);
}

Frequency::Frequency(double frequency, Frequency::Units units)
{
  valueIS_ = sput(frequency, units);
}

Frequency::Frequency(const Frequency &frequency) :
  valueIS_(frequency.valueIS_)
{
}

Frequency::~Frequency()
{
}

double Frequency::sget(double value, const std::string &units)
{
  if(units == "THz" || units == "THZ") {
    return 1.0E-12 * value;
  } else if(units == "GHz" || units == "GHz" || units == "ghz") {
    return 1.0E-9 * value;
  } else if(units == "MHz" || units == "MHZ" || units == "mhz") {
    return 1.0E-6 * value;
  } else if(units == "kHz" || units == "KHZ" || units == "khz") {
    return 1.0E-3 * value;
  } else if(units == "Hz" || units == "HZ" || units == "hz") {
    return value;
  } else {
    return value;
  }
}
double Frequency::sput(double freq, const std::string &units)
{
  if(units == "THz" || units == "THZ") {
    return 1.0E12 * freq;
  } else if(units == "GHz" || units == "GHZ" || units == "ghz") {
    return 1.0E9 * freq;
  } else if(units == "MHz" || units == "MHZ" || units == "mhz") {
    return 1.0E6 * freq;
  } else if(units == "kHz" || units == "KHZ" || units == "khz") {
    return 1.0E3 * freq;
  } else if(units == "Hz" || units == "HZ" || units == "hz") {
    return freq;
  } else {
    return freq;
  }
}

double Frequency::sget(double value, Frequency::Units units)
{
  if(units == Frequency::UnitTeraHertz) {
    return 1.0E-12 * value;
  } else if(units == Frequency::UnitGigaHertz) {
    return 1.0E-9 * value;
  } else if(units == Frequency::UnitMegaHertz) {
    return 1.0E-6 * value;
  } else if(units == Frequency::UnitKiloHertz) {
    return 1.0E-3 * value;
  } else if(units == Frequency::UnitHertz) {
    return value;
  } else {
    return value;
  }
}
double Frequency::sput(double freq, Frequency::Units units)
{
  if(units == Frequency::UnitTeraHertz) {
    return 1.0E12 * freq;
  } else if(units == Frequency::UnitGigaHertz) {
    return 1.0E9 * freq;
  } else if(units == Frequency::UnitMegaHertz) {
    return 1.0E6 * freq;
  } else if(units == Frequency::UnitKiloHertz) {
    return 1.0E3 * freq;
  } else if(units == Frequency::UnitHertz) {
    return freq;
  } else {
    return freq;
  }
}

ATM_NAMESPACE_END
