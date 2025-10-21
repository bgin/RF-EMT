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
 * "@(#) $Id: ATMWaterVaporRadiometer.cpp Exp $"
 *
 * who       when      what
 * --------  --------  ----------------------------------------------
 * pardo     24/03/09  created
 */

#include "ATMWaterVaporRadiometer.h"

#include <vector>
#include <string>
#include <stdio.h>
#include <iostream>



ATM_NAMESPACE_BEGIN

WaterVaporRadiometer::WaterVaporRadiometer(const std::vector<unsigned int> &IdChannels)
{

  spilloverTemperature_ = Temperature(-999.0, Temperature::UnitKelvin);
  IdChannels_ = IdChannels;

  Percent sg(50, Percent::UnitPercent); // IF DOUBLE SIDE BAND, Default Sideband Gain is 50%

  for(unsigned int i = 0; i < IdChannels.size(); i++) {

    skyCoupling_.push_back(1.0);
    signalGain_.push_back(sg);

  }

}

WaterVaporRadiometer::WaterVaporRadiometer(const std::vector<unsigned int> &IdChannels,
                                           const std::vector<double> &skyCoupling)
{

  spilloverTemperature_ = Temperature(-999.0, Temperature::UnitKelvin);
  IdChannels_ = IdChannels;

  /*    if(IdChannels.size()!=skyCoupling.size()){throw Error();} */

  if(IdChannels.size() < skyCoupling.size()) {
    for(unsigned int i = 0; i < IdChannels.size(); i++) {
      skyCoupling_.push_back(skyCoupling[i]);
    }
  } else {
    if(IdChannels.size() == skyCoupling.size()) {
      skyCoupling_ = skyCoupling;
    } else {
      for(unsigned int i = 0; i < skyCoupling.size(); i++) {
        skyCoupling_.push_back(skyCoupling[i]);
      }
      for(unsigned int i = skyCoupling.size(); i < IdChannels.size(); i++) {
        skyCoupling_.push_back(skyCoupling[skyCoupling.size() - 1]);
      }
    }
  }

}

WaterVaporRadiometer::WaterVaporRadiometer(const std::vector<unsigned int> &IdChannels,
                                           const std::vector<Percent> &signalGain)
{

  spilloverTemperature_ = Temperature(-999.0, Temperature::UnitKelvin);
  IdChannels_ = IdChannels;

  if(IdChannels.size() < signalGain.size()) {
    for(unsigned int i = 0; i < IdChannels.size(); i++) {
      signalGain_.push_back(signalGain[i]);
    }
  } else {
    if(IdChannels.size() == signalGain.size()) {
      signalGain_ = signalGain;
    } else {
      for(unsigned int i = 0; i < signalGain.size(); i++) {
        signalGain_.push_back(signalGain[i]);
      }
      for(unsigned int i = signalGain.size(); i < IdChannels.size(); i++) {
        signalGain_.push_back(signalGain[signalGain.size() - 1]);
      }
    }
  }

}

WaterVaporRadiometer::WaterVaporRadiometer(const std::vector<unsigned int> &IdChannels,
                                           const std::vector<double> &skyCoupling,
                                           const std::vector<Percent> &signalGain)
{

  spilloverTemperature_ = Temperature(-999.0, Temperature::UnitKelvin);
  IdChannels_ = IdChannels;

  if(IdChannels.size() < skyCoupling.size()) {
    for(unsigned int i = 0; i < IdChannels.size(); i++) {
      skyCoupling_.push_back(skyCoupling[i]);
    }
  } else {
    if(IdChannels.size() == skyCoupling.size()) {
      skyCoupling_ = skyCoupling;
    } else {
      for(unsigned int i = 0; i < skyCoupling.size(); i++) {
        skyCoupling_.push_back(skyCoupling[i]);
      }
      for(unsigned int i = skyCoupling.size(); i < IdChannels.size(); i++) {
        skyCoupling_.push_back(skyCoupling[skyCoupling.size() - 1]);
      }
    }
  }

  if(IdChannels.size() < signalGain.size()) {
    for(unsigned int i = 0; i < IdChannels.size(); i++) {
      signalGain_.push_back(signalGain[i]);
    }
  } else {
    if(IdChannels.size() == signalGain.size()) {
      signalGain_ = signalGain;
    } else {
      for(unsigned int i = 0; i < signalGain.size(); i++) {
        signalGain_.push_back(signalGain[i]);
      }
      for(unsigned int i = signalGain.size(); i < IdChannels.size(); i++) {
        signalGain_.push_back(signalGain[signalGain.size() - 1]);
      }
    }
  }

}

WaterVaporRadiometer::WaterVaporRadiometer(const std::vector<unsigned int> &IdChannels,
                                           const Temperature &spilloverTemperature)
{

  spilloverTemperature_ = spilloverTemperature;
  IdChannels_ = IdChannels;

  Percent sg(50, Percent::UnitPercent); // IF DOUBLE SIDE BAND, Default Sideband Gain is 50%

  for(unsigned int i = 0; i < IdChannels.size(); i++) {

    skyCoupling_.push_back(1.0);
    signalGain_.push_back(sg);

  }

}

WaterVaporRadiometer::WaterVaporRadiometer(const std::vector<unsigned int> &IdChannels,
                                           const std::vector<double> &skyCoupling,
                                           const Temperature &spilloverTemperature)
{

  spilloverTemperature_ = spilloverTemperature;
  IdChannels_ = IdChannels;

  if(IdChannels.size() < skyCoupling.size()) {
    for(unsigned int i = 0; i < IdChannels.size(); i++) {
      skyCoupling_.push_back(skyCoupling[i]);
    }
  } else {
    if(IdChannels.size() == skyCoupling.size()) {
      skyCoupling_ = skyCoupling;
    } else {
      for(unsigned int i = 0; i < skyCoupling.size(); i++) {
        skyCoupling_.push_back(skyCoupling[i]);
      }
      for(unsigned int i = skyCoupling.size(); i < IdChannels.size(); i++) {
        skyCoupling_.push_back(skyCoupling[skyCoupling.size() - 1]);
      }
    }
  }

}

WaterVaporRadiometer::WaterVaporRadiometer(const std::vector<unsigned int> &IdChannels,
                                           const std::vector<Percent> &signalGain,
                                           const Temperature &spilloverTemperature)
{

  spilloverTemperature_ = spilloverTemperature;
  IdChannels_ = IdChannels;

  if(IdChannels.size() < signalGain.size()) {
    for(unsigned int i = 0; i < IdChannels.size(); i++) {
      signalGain_.push_back(signalGain[i]);
    }
  } else {
    if(IdChannels.size() == signalGain.size()) {
      signalGain_ = signalGain;
    } else {
      for(unsigned int i = 0; i < signalGain.size(); i++) {
        signalGain_.push_back(signalGain[i]);
      }
      for(unsigned int i = signalGain.size(); i < IdChannels.size(); i++) {
        signalGain_.push_back(signalGain[signalGain.size() - 1]);
      }
    }
  }

}

WaterVaporRadiometer::WaterVaporRadiometer(const std::vector<unsigned int> &IdChannels,
                                           const std::vector<double> &skyCoupling,
                                           const std::vector<Percent> &signalGain,
                                           const Temperature &spilloverTemperature)
{

  spilloverTemperature_ = spilloverTemperature;
  IdChannels_ = IdChannels;

  if(IdChannels.size() < skyCoupling.size()) {
    for(unsigned int i = 0; i < IdChannels.size(); i++) {
      skyCoupling_.push_back(skyCoupling[i]);
    }
  } else {
    if(IdChannels.size() == skyCoupling.size()) {
      skyCoupling_ = skyCoupling;
    } else {
      for(unsigned int i = 0; i < skyCoupling.size(); i++) {
        skyCoupling_.push_back(skyCoupling[i]);
      }
      for(unsigned int i = skyCoupling.size(); i < IdChannels.size(); i++) {
        skyCoupling_.push_back(skyCoupling[skyCoupling.size() - 1]);
      }
    }
  }

  if(IdChannels.size() < signalGain.size()) {
    for(unsigned int i = 0; i < IdChannels.size(); i++) {
      signalGain_.push_back(signalGain[i]);
    }
  } else {
    if(IdChannels.size() == signalGain.size()) {
      signalGain_ = signalGain;
    } else {
      for(unsigned int i = 0; i < signalGain.size(); i++) {
        signalGain_.push_back(signalGain[i]);
      }
      for(unsigned int i = signalGain.size(); i < IdChannels.size(); i++) {
        signalGain_.push_back(signalGain[signalGain.size() - 1]);
      }
    }
  }

}

WaterVaporRadiometer::~WaterVaporRadiometer()
{
}

ATM_NAMESPACE_END
