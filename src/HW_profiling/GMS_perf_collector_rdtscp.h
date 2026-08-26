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
#ifndef __GMS_PERF_COLLECTOR_RDTSCP_H__
#define __GMS_PERF_COLLECTOR_RDTSCP_H__

// File version 
namespace file_info {
      const unsigned int gGMS_PERF_COLLECTOR_RDTSCP_MAJOR = 1U;

      const unsigned int gGMS_PERF_COLLECTOR_RDTSCP_MINOR = 0U;

      const unsigned int gGMS_PERF_COLLECTOR_RDTSCP_MICRO = 0U;

      const unsigned int gGMS_PERF_COLLECTOR_RDTSCP_FULLVER = 
	1000U*gGMS_PERF_COLLECTOR_RDTSCP_MAJOR + 100U*gGMS_PERF_COLLECTOR_RDTSCP_MINOR + 10U*gGMS_PERF_COLLECTOR_RDTSCP_MICRO;

      const char * const pgGMS_PERF_COLLECTOR_RDTSCP_CREATE_DATE = "05-10-2019 11:40 +00200 (SAT 05 OCT 2019 GMT+2)";

      const char * const pgGMS_PERF_COLLECTOR_RDTSCP_BUILD_DATE = __DATE__ ":"__TIME__;

      const char * const pgGMS_PERF_COLLECTOR_RDTSCP_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

      const char * const pgGMS_PERF_COLLECTOR_RDTSCP_SYNOPSIS = "Code section performance collector class based on RDTSCP machine-code instruction.";
}


#include <vector>
#include <string>
#include <cstdint>
#include "GMS_config.h"

#if (USE_MKL) == 1 && defined (GMS_COMPILED_BY_ICC)
#include <mkl_dfti.h>
#endif
#if (USE_ACCURATE_IEEE754_2008_FP) == 1
#include <bfp754.h>
#endif

void rdft(int, int, double *, int *, double *);

namespace  gms {
	namespace system {

		struct  PerfCollectorRDTSCP {
				
				bool         m_Iscleared; 

				std::size_t  m_nsamples;

				std::string  m_func_name;

				std::string  m_file_name;

				std::string  m_date;

				std::string  m_time;
				// this member describe successful computation of delta value (success == true, otherwise false)
				std::vector<bool>  m_nvalid_values;
				// profiling samples  start  values
				std::vector<uint64_t>  m_start_values;
				// profiling samples, stop values
				std::vector<uint64_t>  m_stop_values;
				// Time samples delta
				std::vector<uint64_t>  m_delta_values;
				// Output of TSC_AUX MSR (core mapping)
				std::vector<uint32_t>  m_tscaux_start;
				// Same as above (stop measurements)
				std::vector<uint32_t>  m_tscaux_stop;

				constexpr static size_t lo_bound = 3Ui64;

				constexpr static size_t lowest = 1Ui64;

				constexpr static size_t zero = 0Ui64;

				PerfCollectorRDTSCP();

				PerfCollectorRDTSCP(const PerfCollectorRDTSCP &);

				PerfCollectorRDTSCP(PerfCollectorRDTSCP &&) noexcept(true);

				~PerfCollectorRDTSCP() noexcept(true) = default;

				PerfCollectorRDTSCP & operator=(const PerfCollectorRDTSCP &);

				PerfCollectorRDTSCP & operator=(PerfCollectorRDTSCP &&) noexcept(true);

				void start() noexcept(false);

				void stop()  noexcept(false);

				void clear_all();

				bool compute_delta() noexcept(false);

#if   (USE_ACCURATE_IEEE754_2008_FP) == 1
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
				static bool correlation_set(PerfCollectorRDTSCP &,
							    PerfCollectorRDTSCP &,
							    double * __restrict,
							     int32_t * __restrict,
							    double * __restrict,
							    const int32_t );

				void print() const;

				bool check_data_size();

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




#endif /*__LAM_PERF_COLLECTOR_RDTSCP_H__*/
