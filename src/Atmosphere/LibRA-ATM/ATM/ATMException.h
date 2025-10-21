#ifndef _ATM_EXCEPTIONS_H
#define _ATM_EXCEPTIONS_H

/****************************************************************************
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
 * "@(#) $Id: ATMException.h Exp $"
 *
 * who       when      what
 * --------  --------  ----------------------------------------------
 * jroche    09/09/09  created
 */

#ifndef __cplusplus
#error "This is a C++ include file and cannot be used from plain C"
#endif

#include "ATMCommon.h" 
#include <exception>
#include <string>

using std::string;

ATM_NAMESPACE_BEGIN

/**
 * @class AtmException
 * @brief Manage exceptions that raises in ATM.
 *
 * @remarks Don't use a STL container (particularly a std::string) in a catch
 * block of a TelCal exception class because the container may dynamically
 * allocate the memory and so raise a std::exception if the case of
 * std::bad_alloc exception.
 *
 * To create a AtmException follow the example. The #ATM_EXCEPTION_ARGS macro
 * must be used to define the position (file, method and line) where the
 * exception is raised.
 * First include the header file for exceptions.
 *
 * @code
 *  #include "ATMException.h"
 * @endcode
 *
 * Then create an exception object in a try block. Possibly insert here a log
 * message.
 *
 * @code
 *  try
 *  {
 *    if(something wrong)
 *    {
 *      throw atm::AtmException(ATM_EXCEPTION_ARGS("A message."));
 *    }
 *  }
 * @endcode
 *
 * If you don't want to indicate the location in the log message then use the
 * message() method instead of what().
 * Finally catch the exception.
 *
 * @code
 *  catch(const atm::AtmException &e)
 *  {
 *    // do something like destroy the allocated memory in the try block...
 *    throw e;
 *  }
 * @endcode
 *
 * In this catch block you can also modify the exception before throwing it. To
 * do it remove the const keyword in the catch call.
 *
 * The message may be formated using a stream. So first include the needed
 * stream library.
 *
 * @code
 *  #include "ATMException.h"
 *  #include <sstream>
 * @endcode
 *
 * Format the message like the following example where x is a variable of
 * miscellaneous type. See following web links:
 * @li ostringstream class: http://www.cplusplus.com/reference/iostream/ostringstream
 * @li ostream::operator<< operator: http://www.cplusplus.com/reference/iostream/ostream/operator%3C%3C.html
 *
 * @code
 *  std::ostringstream oss;
 *  oss << "A message " << x << " end of the message.";
 *  throw atm::AtmException(ATM_EXCEPTION_ARGS(oss.str().c_str()));
 * @endcode
 *
 * It is recommended to catch unknown exceptions with an AtmException.
 *
 * @code
 *  catch(...)
 *  {
 *    throw atm::AtmException(ATM_EXCEPTION_ARGS("Uncaught exception!"));
 *  }
 * @endcode
 *
 * To change the format of the outputed method, change the pattern of trace
 * messages using the telcal::Trace::setTracePattern() method.
 */
class AtmException : public std::exception
{
public:
  // --------------------------------------------------
  //@{
  //! @name Constructors and destructor
  // --------------------------------------------------
  /**
    * @brief Constructor.
    * @param msg a message that describes the exception.
    */
  AtmException(const char* msg) throw();
  /**
   * @brief Constructor.
   * @param file the file name where the exception raised.
   * @param routine the routine name where the exception raised.
   * @param line the line where the exception raised.
   * @param msg a message that describes the exception.
   * @see formatMsg()
   *
   * This constructor define the AtmException::what_m attribute using the
   * output string formatted with the formatMsg() method.
   */
  AtmException(const char* file, const char* routine, int line, const char* msg) throw();
  /**
   * @brief Destructor.
   */
  virtual ~AtmException() throw()
  {
  }
  //@}

  // --------------------------------------------------
  //@{
  //! @name Get attributes
  // --------------------------------------------------
  /**
   * @brief Returns the file name where the current error raised.
   * @return A C-style character string.
   */
  virtual const char* file() const throw() { return file_m.c_str(); }
  /**
    * @brief Returns the routine name where the current error raised.
    * @return A C-style character string.
    */
  virtual const char* routine() const throw() { return routine_m.c_str(); }
  /**
    * @brief Returns the line number where the current error raised.
    * @return An integer.
    * @remarks -1 is returned if no line was given.
    */
  virtual int line() const throw() { return line_m; }
  /**
    * @brief Returns the message where the current error raised.
    * @return A C-style character string.
   */
  virtual const char* message() const throw() { return msg_m.c_str(); }
  /**
    * @brief Returns a string describing the general cause of the current error.
    * @return A C-style character string.
    */
  virtual const char* what() const throw() { return what_m.c_str(); }
  //@}

private:
  // --------------------------------------------------
  //@{
  //! @name Manage messages
  // --------------------------------------------------
  /**
   * @brief Output a message that conform to a given trace pattern.
   * @param file a file name.
   * @param routine a routine name.
   * @param line a line number.
   * @param msg a message.
   * @return The formatted message.
   */
  string formatMsg(const char* file, const char* routine, int line, const char* msg);
  //@}

private:
  string file_m; //!< The file where the exception raised.
  string routine_m; //!< The routine where the exception raised.
  int line_m; //!< The line where the exception raised.
  string msg_m; //!< The message where the exception raised.
  string what_m; //!< The error message to be displayed.
}; // class AtmException

ATM_NAMESPACE_END

/**
 * @def ATM_EXCEPTION_ARGS(msg)
 * @brief Writes arguments needed to create an AtmException object.
 * @param msg a message (that describes the exception).
 * @see AtmException
 *
 * Writes (during the C++ code parsing) the list of arguments that contains the
 * position (file, function and line) in the C++ module where this macro is
 * written and the exception description \p msg.
 *
 * Use this macro to create an AtmException object:
 *
 * @code
 *  atm::AtmException(ATM_EXCEPTION_ARGS("An exception message..."));
 * @endcode
 */
#define ATM_EXCEPTION_ARGS(msg) __FILE__, __PRETTY_FUNCTION__, __LINE__, msg

#endif /*!_ATM_EXCEPTIONS_H*/
