
/* Copyright (C) Bernard Gingold, 2020-2026 
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#ifndef __GMS_ERROR_MACROS_H__
#define __GMS_ERROR_MACROS_H__ 280920191202







#include <iostream>
#include <iomanip>







#if !defined (ABORT_ON_ERROR)
#define ABORT_ON_ERROR(msg,err) \
	std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << (msg) << "\n"; \
	std::cerr << "at " << __FILE__ << ":" << __LINE__ << "\n"; \
	std::exit(err);
#endif	


#if !defined (THROW_ON_RUNTIME_ERROR)
#define THROW_ON_RUNTIME_ERROR(msg)  \
	std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << (msg) << "\n"; \
	std::cerr << "at " << __FILE__ << ":" << __LINE__ <<  "\n"; \
	throw std::runtime_error(msg);
#endif

#if !defined (THROW_ON_DOMAIN_ERROR)
#define THROW_ON_DOMAIN_ERROR(msg,value) \
	std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << (msg) << " value: " << std::fixed << std::setprecision(15) << (value) << "\n"; \
	std::cerr << "at " << __FILE__ << ":" << __LINE__ << "\n"; \
	throw std::domain_error(msg);
#endif

#if !defined (THROW_ON_OVERFLOW_ERROR)
#define THROW_ON_OVERFLOW_ERROR(msg,value) \
	std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << (msg) << " value: " << std::fixed << std::setprecision << static_cast<double>((value)) << "\n"; \
	std::cerr << "at " << __FILE__ << ":" << __LINE__ <<  "\n"; \
	throw std::overflow_error(msg);
#endif

#if !defined (THROW_ON_RANGE_ERROR)
#define THROW_ON_RANGE_ERROR(msg,value) \
	std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << (msg) << " value: " << std::fixed << std::setprecision << static_cast<double>((value)) << "\n"; \
	std::cerr << "at " << __FILE__ << ":" << __LINE__  << "\n"; \
	throw std::range_error(msg);
#endif

#if !defined (THROW_ON_INVALID_ARGUMENT_ERROR)
#define THROW_ON_INVALID_ARGUMENT_ERROR(msg) \
	std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << (msg) << "\n"; \
	std::cerr << "at " << __FILE__ << ":" << __LINE__ <<  "\n"; \
	throw std::invalid_argument(msg);
#endif

#if !defined (THROW_ON_LOGIC_ERROR)
#define THROW_ON_LOGIC_ERROR(msg) \
	std::cerr << "[" << __DATE__ << ":" << __TIME__ << "]" << (msg) << "\n"; \
	std::cerr << "at " << __FILE__ << ":" << __LINE__ <<  "\n"; \
	throw std::logic_error(msg);
#endif

#if !defined (MALLOC_FAILED)
#define MALLOC_FAILED -1
#endif

#endif /*__GMS_ERROR_MACROS_H__*/
