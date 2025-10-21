#ifndef _ATM_ERROR_H
#define _ATM_ERROR_H
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
 * "@(#) $Id: ATMError.h Exp $"
 *
 * who       when      what
 * --------  --------  ----------------------------------------------
 * pardo     24/03/09  created
 */

#ifndef __cplusplus
#error "This is a C++ include file and cannot be used from plain C"
#endif

#include "ATMCommon.h"
#include <string>
#include <stdarg.h>



ATM_NAMESPACE_BEGIN

  /** Definition of the sequence of error levels sorted by ascending severity.
 */
enum ErrorLevel
{
  NOERROR = 0,
  WARNING,
  MINOR,
  SERIOUS,
  FATAL
};

/** The Error class represents a list of exceptions that may happen when errors occur.
    To every error is associated a severity level. There is at any time an acceptable level
    of severity. If an error occurs at that level or a more severe level, an exception is thrown.
    The critical level of severity below which errors are tolerated can be modified.
*/
class Error
{
public:
  Error();
  Error(ErrorLevel errlev);
  /** Declare the occurrence of an error
   *  @pre the default level of tolerance for an acceptable error is used or a different critical
   *  level has been set.
   *  @param errlev the level of severity for the error being declared
   *  @param message Error message
   *  @post the declaration may have or may have not thrown an exception depending on the current
   *  tolerance about the severity level.
   */
  Error(ErrorLevel errlev, const std::string &message);

  /** Declare the occurrence of an error
   *  @pre the default level of tolerance for an acceptable error is used or a different critical
   *  level has been set.
   *  @param errlev the level of severity for the error being declared
   *  @param fmt Error message format as in a standard C printf
   *  @param ... Variable arguments as in printf
   *  @post the declaration may have or may have not thrown an exception depending on the current
   *  tolerance about the severity level.
   */
  Error(ErrorLevel errlev, char* fmt, ...);

  /** Destructor */
  virtual ~Error();

  void notify(ErrorLevel errlev, const std::string &message);
  void notify(const std::string &message);

  /** Setter to modify the critical current level of severity below which errors are tolerated.
   *  @param criticalErrlev critical limit: when errlev >= criticalErrlev an exception will be thrown
   *  @post  New declarations of errors will be thrown or not according to this new critical limit
   */
  static void setAcceptableLevel(ErrorLevel criticalErrlev);

  /** Accessor to the current limit below wich error are not thrown
   */
  static ErrorLevel getAcceptableLevel();

  /** Accessor to the current limit in the form of a std::string
   */
  static std::string getAcceptableLevelToString();

  /** Accessor to the severity level of the error which has just been declared
   */
  static ErrorLevel getLevel();

  /** Accessor to the severity level, in the form of a std::string, of the error which has just been declared
   */
  static std::string getLevelToString();

  /** Accessor to the error message for the error which has just been declared
   */
  static std::string getErrorMessage();

  /** Clear the stack of all the error messages which passed up to the last declaration which has
   *  finally triggered the exception
   */
  static void clearMessage();

  /** Clear the stack of all the error levels which passed up to the last declaration which has
   *  finally triggered the exception
   */
  static void clearErrLev();

private:
  void printMessage(const std::string &message); //!< Print utility

private:
  static ErrorLevel acceptableErrorLevel; //!< Current acceptable level of severity
  static std::string errorMessage_; //!< Error message
  static ErrorLevel errlev_; //!< Error level
  std::string errorLog; //!< Error log (TBD)
}; // class Error

ATM_NAMESPACE_END

/** \page testError Example showing how to use the class Error

 \dontinclude testError.cpp

 The codes using the Error class must include its header file.
 \skip include
 \until using

 To illustrate how to use this class we take here an example defining a class Array with a constraint
 about a maximum index.
 \skip const
 \until testbed

 In the test program we begin by showing the default critical limit for acceptable errors:
 \skip main()
 \until default

 An array is declared and it is populated by integers for indices from 0 to 100. The opperator []
 generate an Error throwing it if the index exceeds the maximum allowed limit, we insert the loop in a try
 section. When the index j exceeds the upper limit in the indexing, an Error is thrown unless
 this type type of errors is acceptable.
 \skip Array
 \until Error

 The out-of-range indexing in the [] operator methods having been declared to be SERIOUS, when j
 exceeds the maximum size allowed, its severity, SERIOUS (see the method []) is just at the critical
 when the error become unacceptable. Hence an Error exception is thrown.
 \skip getLevel
 \until input!

 The execution proceeds; the following test attests that this error was considered as acceptable:
 \skip if(Error
 \until proceed

 We now repeat this scenario triggering with the same error but before that we  test with a tolerence
 which is more constraining: any error will be thrown only when it severity reach the FATAL limit:
 \from setAcceptableLevel
 \until cout

 The scenario is now repeated:
 \skip Test:
 \until return

 The Error is not catched because with this new contextthe error which is still SERIOUS is considered
 as acceptable.

 */

#endif /*!ATM_ERROR_H*/
