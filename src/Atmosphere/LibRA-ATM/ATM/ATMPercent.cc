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
 * "@(#) $Id: ATMPercent.cpp Exp $"
 *
 * who       when      what
 * --------  --------  ----------------------------------------------
 * pardo     24/03/09  created
 */

#include "ATMPercent.h"



ATM_NAMESPACE_BEGIN

  // Constructors

  Percent::Percent() : valueIS_(0.001){}      // Default humidity 0.1 %

  //  Percent::Percent(double percent) : valueIS_(percent){}

  Percent::Percent(double percent){
    if(percent <= 0.0){
      valueIS_=0.001;
    }else{
      valueIS_=percent/100.0;
    }
  }

  Percent::Percent(double percent, const std::string &units){
    if(units == "%" || units == "percent" || units == "PERCENT"){
      valueIS_ = percent / 100.0;
    } else {
      // Exception: Unknown percent unit
      valueIS_ = percent;
    }

    if(valueIS_ <= 0.0){valueIS_=0.001;}    // default humidity 0.1 %

  }

  Percent::Percent(double percent, Percent::Units units) {
    if (units == Percent::UnitPercent) {
      valueIS_ = percent / 100.0;
    } else {
      // Exception: Unknown percent unit
      valueIS_ = percent;
    }

    if(valueIS_ <= 0.0){valueIS_=0.001;}    // default humidity 0.1 %
  }

  // Destructor
  Percent::~Percent(){}

  // Accessors
  double Percent::get()const{return valueIS_;}
  double Percent::get(const std::string &units)const{
    if(units == "%" || units == "percent" || units == "PERCENT"){
      return valueIS_ * 100.0;
    } else {
      // Exception: Unknown percent unit
      return valueIS_;
    }
  }
  double Percent::get(Percent::Units units)const{
    if(units == Percent::UnitPercent){
      return valueIS_ * 100.0;
    } else {
      // Exception: Unknown percent unit
      return valueIS_;
    }
  }

ATM_NAMESPACE_END
