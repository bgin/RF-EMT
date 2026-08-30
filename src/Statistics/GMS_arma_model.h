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
#ifndef _GMS_ARMA_MODEL_H_
#define _GMS_ARMA_MODEL_H_



namespace file_info {


	const unsigned int gGMS_ARMA_MODEL_MAJOR = 1;

	const unsigned int gGMS_ARMA_MODEL_MINOR = 1;

	const unsigned int gGMS_ARMA_MODEL_MICRO = 0;

	const unsigned int gGMS_ARMA_MODEL_FULLVER = 
	 1000U*gGMS_ARMA_MODEL_MAJOR+100U*gGMS_ARMA_MODEL_MINOR+10U*gGMS_ARMA_MODEL_MICRO;

	const char * const pgGMS_ARMA_MODEL_CREATE_DATE = "06-10-2019 14:31 +00200 (SUN 06 OCT 2019 GMT+2)";

	const char * const pgGMS_ARMA_MODEL_BUILD_DATE = __DATE__ : __TIME__;

	const char * const pgGMS_ARMA_MODEL_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

	const char * const pgGMS_ARMA_MODEL_SYNOPSIS = "AUTOMATIC AR-MA MODEL FITTING -- SCALAR CASE"
}


#include "GMS_config.h"
#if defined (__INTEL_COMPILER)
#include <../perf_headers/c++/valarray>
#else
#include <valarray>
#endif

#include <cstdint>

namespace gms {
	namespace math {
		

		  __ATTR_ALIGN__(64) struct ArmaModel {
		  
				using VAf64 = std::valarray<double>;
				using VAi32 = std::valarray<int32_t>;
			
				/*
				Low-level computational code -- members are public.
				*/
			
				VAi32 *   __restrict    m_iq;

				VAf64 *   __restrict    m_b2;

				VAi32 *   __restrict    m_ip;

				VAf64 *   __restrict    m_a2;

				VAf64 *   __restrict    m_std;

				VAf64 *   __restrict    m_cxx2;

				VAf64 *   __restrict    m_g;

				VAf64 *   __restrict    m_saic;

				double   m_aicm;

#if (USE_STRUCT_PADDING) == 1
				PAD_TO(1,32)
#endif

				int32_t  m_newn;

				int32_t  m_nmax;

				int32_t  m_mmax;

				int32_t  m_kq;

				int32_t  m_kp;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(2,4)
#endif
				

				

				ArmaModel();

				ArmaModel(const int32_t,
					  const int32_t,
					  const int32_t,
					  const double);

				ArmaModel(const ArmaModel &);

				ArmaModel(ArmaModel &&) noexcept(true);

				~ArmaModel() = default;

				ArmaModel & operator=(const ArmaModel &);

				ArmaModel & operator=(ArmaModel &&) noexcept(true);

				void compute_arma_model(int32_t,
							int32_t,
							VAf64 &,
							int32_t,
							VAi32 &,
							VAf64 &,
							VAi32 &,
							VAf64 &,
							int32_t);
							


			};
		}
	}



#endif  /*_GMS_ARMA_MODEL_H_*/
