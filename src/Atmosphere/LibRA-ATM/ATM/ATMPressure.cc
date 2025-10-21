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
 * "@(#) $Id: ATMPressure.cpp Exp $"
 *
 * who       when      what
 * --------  --------  ----------------------------------------------
 * pardo     24/03/09  created
 */
#include "ATMPressure.h"



ATM_NAMESPACE_BEGIN

Pressure::Pressure() :
  valueIS_(0.0)
{
}
Pressure::Pressure(double pressure) :
  valueIS_(pressure)
{
}

Pressure::Pressure(double pressure, const std::string &units)
{
  if(units == "Pa" || units == "PA") {
    valueIS_ = pressure;
  } else if(units == "hPa" || units == "HPA") {
    valueIS_ = pressure * 100.0;
  } else if(units == "bar" || units == "BAR") {
    valueIS_ = 1.0E+5 * pressure;
  } else if(units == "mb" || units == "MB") {
    valueIS_ = 1.0E+2 * pressure;
  } else if(units == "mbar" || units == "MBAR") {
    valueIS_ = 1.0E+2 * pressure;
  } else if(units == "atm" || units == "ATM") {
    valueIS_ = pressure * 101325.;
  } else {
    valueIS_ = pressure;
  }

}

Pressure::Pressure(double pressure, Pressure::Units units)
{
  if(units == Pressure::UnitPascal) {
    valueIS_ = pressure;
  } else if(units == Pressure::UnitHectoPascal) {
    valueIS_ = pressure * 100.0;
  } else if(units == Pressure::UnitBar) {
    valueIS_ = 1.0E+5 * pressure;
  } else if(units == Pressure::UnitMilliBar) {
    valueIS_ = 1.0E+2 * pressure;
  } else if(units == Pressure::UnitAtmosphere) {
    valueIS_ = pressure * 101325.;
  } else {
    valueIS_ = pressure;
  }

}

Pressure::~Pressure()
{
}

double Pressure::get(std::string const &units) const
{
  if(units == "Pa" || units == "PA") {
    return valueIS_;
  } else if(units == "hPa" || units == "HPA" || units == "hpa") {
    return 1.0E-2 * valueIS_;
  } else if(units == "bar" || units == "BAR") {
    return 1.0E-5 * valueIS_;
  } else if(units == "mb" || units == "MB") {
    return 1.0E-2 * valueIS_;
  } else if(units == "mbar" || units == "MBAR") {
    return 1.0E-2 * valueIS_;
  } else if(units == "atm" || units == "ATM") {
    return valueIS_ / 101325.;
  } else {
    return valueIS_;
  }

}

double Pressure::get(Pressure::Units units) const
{
  if(units == Pressure::UnitPascal) {
    return valueIS_;
  } else if(units == Pressure::UnitHectoPascal) {
    return 1.0E-2 * valueIS_;
  } else if(units == Pressure::UnitBar) {
    return 1.0E-5 * valueIS_;
  } else if(units == Pressure::UnitMilliBar) {
    return 1.0E-2 * valueIS_;
  } else if(units == Pressure::UnitAtmosphere) {
    return valueIS_ / 101325.;
  } else {
    return valueIS_;
  }

}

ATM_NAMESPACE_END
