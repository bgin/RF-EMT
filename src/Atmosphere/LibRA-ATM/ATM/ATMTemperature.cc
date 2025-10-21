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
 * "@(#) $Id: ATMTemperature.cpp Exp $"
 *
 * who       when      what
 * --------  --------  ----------------------------------------------
 * pardo     24/03/09  created
 */

#include "ATMTemperature.h"
#include "ATMCommon.h"

ATM_NAMESPACE_BEGIN

Temperature::Temperature() :
  valueIS_(0.0)
{
}

Temperature::Temperature(double temperature) :
  valueIS_(temperature)
{
}

Temperature::Temperature(double temperature, const string &units)
{
  if(units == "mK" || units == "mk") {
    valueIS_ = 1.0E-3 * temperature;
  } else if(units == "K" || units == "k") {
    valueIS_ = temperature;
  } else if(units == "C" || units == "c") {
    valueIS_ = temperature + 273.16;
  } else if(units == "F" || units == "f") {
    valueIS_ = (temperature-32.0)*(5./9.)+273.16;
  } else {
    // Exception: Unknown temperature unit. S.I. used (Kelvin)
    valueIS_ = temperature;
  }
}

Temperature::Temperature(double temperature, const Temperature::Units units)
{
  if(units == Temperature::UnitMilliKelvin) {
    valueIS_ = 1.0E-3 * temperature;
  } else if(units == Temperature::UnitKelvin) {
    valueIS_ = temperature;
  } else if(units == Temperature::UnitCelsius) {
    valueIS_ = temperature + 273.16;
  } else if(units == Temperature::UnitFahrenheit) {
    valueIS_ = (temperature-32.0)*(5./9.)+273.16;
  } else {
    // Exception: Unknown temperature unit. S.I. used (Kelvin)
    valueIS_ = temperature;
  }
}

double Temperature::get(const string &units) const
{
  if(units == "mK") {
    return 1.0E3 * valueIS_;
  } else if(units == "K") {
    return valueIS_;
  } else if(units == "C" || units == "c") {
    return valueIS_ - 273.16;
  } else if(units == "F" || units == "f") {
    return (valueIS_ - 273.16)*(9./5.)+32.0;
  } else {
    // Exception: Unknown temperature unit. S.I. used (Kelvin)
    return valueIS_;
  }
}

double Temperature::get(Temperature::Units units) const
{
  if(units == Temperature::UnitMilliKelvin) {
    return 1.0E3 * valueIS_;
  } else if(units == Temperature::UnitKelvin) {
    return valueIS_;
  } else if(units == Temperature::UnitCelsius) {
    return valueIS_ - 273.16;
  } else if(units == Temperature::UnitFahrenheit) {
    return (valueIS_ - 273.16)*(9./5.)+32.0;
  } else {
    // Exception: Unknown temperature unit. S.I. used (Kelvin)
    return valueIS_;
  }
}

ATM_NAMESPACE_END
