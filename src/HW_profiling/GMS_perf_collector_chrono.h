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
#ifndef __GMS_PERF_COLLECTOR_CHRONO_H__
#define __GMS_PERF_COLLECTOR_CHRONO_H__

// File version info
namespace file_info {
     const unsigned int gGMS_PERF_COLLECTOR_CHRONO_MAJOR = 1U;

     const unsigned int gGMS_PERF_COLLECTOR_CHRONO_MINOR = 0U;

     const unsigned int gGMS_PERF_COLLECTOR_CHRONO_MICRO = 0U;

     const unsigned int gGMS_PERF_COLLECTOR_CHRONO_FULLVER = 
	1000U*gGMS_PERF_COLLECTOR_CHRONO_MAJOR+100U*gGMS_PERF_COLLECTOR_CHRONO_MINOR+10U*gGMS_PERF_COLLECTOR_CHRONO_MICRO;

     const char * const pgGMS_PERF_COLLECTOR_CHRONO_CREATE_DATE = "06-09-2018 11:32 +00200 (THR 06 SEP 2018 GMT+2)";

     const char * const pgGMS_PERF_COLLECTOR_CHRONO_BUILD_DATE = "00-00-0000 00:00";

     const char * const pgGMS_PERF_COLLECTOR_CHRONO_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

     const char * const pgGMS_PERF_COLLECTOR_CHRONO_SYNOPSIS = "Code section performance collector class based on C++ std library chrono header file.";
}
#include <vector>
#include <string>
#include <cstdint>
#include <chrono>


#include "GMS_config.h"

#if (USE_MKL) && defined (GMS_COMPILED_BY_ICC)
#include <mkl_dfti.h>
#endif
#if (USE_ACCURATE_IEEE754_2008_FP) == 1
#include <bfp754.h>
#endif

void rdft(int, int, double *, int *, double *);

namespace gms {
	namespace system {

		struct PerfCollectorChrono{

			    using PERIOD = std::chrono::duration<double>;

			    using TIME_VAL = std::chrono::time_point<std::chrono::high_resolution_clock>; // <--- STD Lib OOP bullshit
				
				bool		  m_Iscleared;

				std::size_t   m_nsamples;

				std::string   m_funcname;

				std::string   m_filename;

				std::string   m_date;

				std::string   m_time;

				std::vector<bool> m_nvalid_values;

				std::vector<TIME_VAL> m_start_values;

				std::vector<TIME_VAL> m_stop_values;

				std::vector<double> m_delta_values;

				constexpr static size_t  lo_bound = 3Ui64;

				constexpr static size_t  lowest = 1Ui64;

				constexpr static size_t  zero = 0Ui64;

				PerfCollectorChrono();

				PerfCollectorChrono(_In_ const PerfCollectorChrono &);

				PerfCollectorChrono(_In_ PerfCollectorChrono &&) noexcept(true);

				~PerfCollectorChrono() noexcept(true);

				PerfCollectorChrono & operator=(_In_ const PerfCollectorChrono &);

				PerfCollectorChrono & operator=(_In_ PerfCollectorChrono &&) noexcept(true);

				void start() noexcept(false);

				void stop() noexcept(false);

				void clear_all();

				bool compute_delta() noexcept(false);

#if (USE_ACCURATE_IEEE754_2008_FP) == 1
				// This version relies on Intel libbfp754 and that means
				// probably slower execution speed (not confirmed yet).
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
				/*
					If using MKL for FFT pass null pointers to
					correlation_set (those pointers are used to
					call Ouru FFT rdft function)
					Usage of MKL DFTI is set as a default.
				*/
				static bool correlation_set(PerfCollectorChrono &,
							    PerfCollectorChrono &,
							    double * __restrict,
							    int32_t * __restrict,
							    double * __restrict,
							    const int32_t);

				void print() const;

				bool check_data_size() const noexcept(true);

#if (USE_ACCURATE_IEEE754_2008_FP) == 1
				// This version relies on Intel libbfp754 and that means
				// probably slower execution speed (not confirmed yet).
				bool cvrt_to_double(double * __restrict,
						    const std::size_t) noexcept(true);
#else
				bool cvrt_to_double(double * __restrict,
					            const std::size_t) noexcept(true);
#endif

                              

		};
	}
}



#endif /*__GMS_PERF_COLLECTOR_CHRONO_H__*/
