/**************************************************************************
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
 * "@(#) $Id: ATMException.cpp Exp $"
 *
 * who       when      what
 * --------  --------  ----------------------------------------------
 * jroche    09/09/09  created
 */

#include "ATMException.h"
#include <sstream>



ATM_NAMESPACE_BEGIN

AtmException::AtmException(const char *msg) throw()
{
  std::ostringstream oss;
  oss << "";  file_m    = oss.str(); oss.str(string());
  oss << "";	routine_m = oss.str(); oss.str(string());
  oss << msg;	msg_m     = oss.str(); oss.str(string());
  line_m = 0; // The first line of a file is the line number 1.
  what_m = msg_m;
}

AtmException::AtmException(const char* file, const char* routine, int line, const char* msg) throw()
{
  std::ostringstream oss;
  oss << file;    file_m    = oss.str(); oss.str(string());
  oss << routine;	routine_m = oss.str(); oss.str(string());
  oss << msg;	    msg_m     = oss.str();
  line_m = line;
  what_m = formatMsg(file, routine, line, msg);
}

string AtmException::formatMsg(const char* file, const char* routine, int line, const char* msg)
{
  static const std::string pattern = "[%F, %R, %L]: %M";
  std::ostringstream oss;

  for(string::const_iterator it(pattern.begin()), it_max(pattern.end()); it != it_max; ++it) {
    // A pattern character is found. Search for a known pattern.
    if((*it == '%') && (it + 1 != it_max)) {
      char next = *(it + 1);
      if(next == 'F') {
        oss << file;
        ++it;
      } else if(next == 'R') {
        oss << routine;
        ++it;
      } else if(next == 'L') {
        oss << line;
        ++it;
      } else if(next == 'M') {
        oss << msg;
        ++it;
      } else {
        oss << *it;
      }
    }

    // Add the current character.
    else oss << *it;
  }
  return oss.str();
}

ATM_NAMESPACE_END
