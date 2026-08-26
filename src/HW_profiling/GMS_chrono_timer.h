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
#ifndef __GMS_CHRONO_TIMER_H__
#define __GMS_CHRONO_TIMER_H__

// Timer based on std::chrono libarary, class -- high_performance_clock

namespace file_info {
  
      const unsigned int gGMS_CHRONO_TIMER_MAJOR = 1U;

      const unsigned int gGMS_CHRONO_TIMER_MINOR = 0U;

      const unsigned int gGMS_CHRONO_TIMER_MICRO = 0U;

      const unsigned int gGMS_CHRONO_TIMER_FULLVER = 
		1000U*gGMS_CHRONO_TIMER_MAJOR + 100U*gGMS_CHRONO_TIMER_MINOR + 10U*gGMS_CHRONO_TIMER_MICRO;

      const char * const pgGMS_CHRONO_TIMER_CREATE_DATE = "04-10-2019 20:26 +00200 (THR 04 OCT 2019 GMT+2)";
// File build date (should be set after successful build of this translation unit.
      const char * const pgGMS_CHRONO_TIMER_BUILD_DATE = "00-00-0000 00:00";

      const char * const pgGMS_CHRONO_TIMER_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com.";

      const char * const pgGMS_CHRONO_TIMER_SYNOPSIS = "std::chrono library based performance measurement class.";
}


#include <chrono>
#include <vector>
#include <string>
#include "GMS_config.h"
#if (USE_ACCURATE_IEEE754_2008_FP) == 1
#include <bfp754.h>
#endif

namespace gms {
	namespace system {

		struct HRCTimer {

			        using fptr = void(*)(void);

				using PERIOD = std::chrono::duration<double>;

				using TIME_VAL = std::chrono::time_point<std::chrono::high_resolution_clock>; // <--- STD Lib OOP bullshit

				std::size_t  m_nruns;

				std::size_t  m_datum_size;

				fptr         m_ptfunc;

				fptr         m_ptwrap;

				std::string  m_func_name;

				std::string  m_file_name;

				std::vector<double> m_timing_values;

				std::vector<double> m_overhead_values;

				std::vector<double> m_delta_values;

				static constexpr uint64_t lo_bound = 3Ui64;

				static constexpr uint64_t zero = 0Ui64;

				HRCTimer() = default;

				HRCTimer(const std::size_t,
					 const std::size_t,
					 const fptr,
					 const fptr,
					 const char *,
					 const char *);

				HRCTimer(const HRCTimer &);

				HRCTimer(HRCTimer &&) noexcept(true);

				inline ~HRCTimer() noexcept(true);

				HRCTimer & operator=(const HRCTimer &);

				HRCTimer & operator=(HRCTimer &&) noexcept(true);

				inline std::size_t get_datum_size() const noexcept(true) { return (m_nruns * m_datum_size); }

				void run_benchmark(const double) noexcept(false);

#if (USE_ACCURATE_IEEE754_2008_FP) == 1
				bool compute_stats(double &,
						   double &,
						   double &,
						   double &,
						   double &);
								
#else
				bool compute_stats(double &,
						   double &,
						   double &,
						   double &,
						   double &);
								  
#endif
				void print() const;

				static bool delta_values_eq(std::vector<bool> &,
							    const HRCTimer &,
							    const HRCTimer & ,
							    const double);

				static bool delta_values_ineq(std::vector<bool> &,
							      const HRCTimer &,
							      const HRCTimer &,
							      const double	);

				static bool delta_values_lt( std::vector<bool> &,
							     const HRCTimer &,
							     const HRCTimer &,
							     const double);

				static bool delta_values_gt( std::vector<bool> &,
							     const HRCTimer &,
							     const HRCTimer &,
							     const double );
				
		};
	}
}

#endif /*__GMS_CHRONO_TIMER_H__*/
