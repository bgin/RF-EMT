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
 * "@(#) $Id: ATMError.cpp Exp $"
 *
 * who       when      what
 * --------  --------  ----------------------------------------------
 * pardo     24/03/09  created
 */

#include "ATMError.h"

#include <stdio.h>
#include <iostream>



ATM_NAMESPACE_BEGIN

ErrorLevel Error::acceptableErrorLevel = SERIOUS; // default level to stop a process
ErrorLevel Error::errlev_ = NOERROR;
std::string Error::errorMessage_ = "";

Error::Error()
{
}

Error::Error(ErrorLevel errlev)
{
  acceptableErrorLevel = errlev;
}

Error::Error(ErrorLevel errlev, const std::string &message)
{
  errlev_ = errlev;
  errorMessage_ = message;
  printMessage(message);
  if(errlev >= acceptableErrorLevel) throw Error();
}

Error::Error(ErrorLevel errlev, char *fmt, ...)
{
  errlev_ = errlev;
  char buffer[1024];
  va_list args;
  va_start(args, fmt);
  vsprintf(buffer, fmt, args);
  errorMessage_ = buffer;
  printMessage(errorMessage_);
  if(errlev >= acceptableErrorLevel)
  {
	  va_end(args);
	  throw Error();
  }
}

Error::~Error()
{
}

void Error::notify(ErrorLevel errlev, const std::string &message)
{
  errlev_ = errlev;
  errorMessage_ = message;
  printMessage(message);
  //if(errlev>=acceptableErrorLevel)exit(-1);
  if(errlev >= acceptableErrorLevel) throw Error();
}

void Error::notify(const std::string &message)
{
  errorMessage_ = message;
  std::cout << message << std::endl;
}

void Error::setAcceptableLevel(ErrorLevel errlev)
{
  acceptableErrorLevel = errlev;
}

ErrorLevel Error::getLevel()
{
  return errlev_;
}

std::string Error::getLevelToString()
{
  if(errlev_ == MINOR) return "MINOR";
  if(errlev_ == WARNING) return "WARNING";
  if(errlev_ == SERIOUS) return "SERIOUS";
  if(errlev_ == FATAL) return "FATAL";
  return "";
}

std::string Error::getErrorMessage()
{
  std::string errorMessage = errorMessage_;
  clearMessage();
  return errorMessage;
}

void Error::clearMessage()
{
  errorMessage_ = "";
  clearErrLev();

}

void Error::clearErrLev()
{
  errlev_ = NOERROR;
}

ErrorLevel Error::getAcceptableLevel()
{
  return acceptableErrorLevel;
}

std::string Error::getAcceptableLevelToString()
{
  if(acceptableErrorLevel == MINOR) return "MINOR";
  if(acceptableErrorLevel == WARNING) return "WARNING";
  if(acceptableErrorLevel == SERIOUS) return "SERIOUS";
  if(acceptableErrorLevel == FATAL) return "FATAL";
  return "NOERROR";
}

void Error::printMessage(const std::string &message)
{
  switch(errlev_) {
  case FATAL:
    std::cout << "FATAL ERROR: " + message << std::endl;
    break;
  case SERIOUS:
    std::cout << "SERIOUS ERROR: " + message << std::endl;
    break;
  case MINOR:
    std::cout << "MINOR ERROR: " + message << std::endl;
    break;
  case WARNING:
    std::cout << "WARNING ERROR: " + message << std::endl;
    break;
  default:
    std::cout << "ERROR: " + message << std::endl;
  }
}

ATM_NAMESPACE_END
