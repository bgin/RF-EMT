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
#ifndef _GMS_BAYESIAN_ESTIMATE_H_
#define _GMS_BAYESIAN_ESTIMATE_H__


namespace file_info {

	const unsigned int gGMS_BASYESIAN_ESTIMATE_MAJOR = 1;

	const unsigned int gGMS_BAYESIAN_ESTIMATE_MINOR =  1;

	const unsigned int gGMS_BAYESIAN_ESTIMATE_MICRO =  0;

	const unsigned int gGMS_BAYESIAN_ESTIMATE_FULLVER = 
	 1000U*gGMS_BAYESIAN_ESTIMATE_MAJOR+100U*gGMS_BAYESIAN_ESTIMATE_MINOR+10U*gGMS_BAYESIAN_ESTIMATE_MICRO;

	const char * const pgGMS_BAYESIAN_ESTIMATE_CREATE_DATE = "06-10-2019 14:31 +00200 (SUN 06 OCT 2019 GMT+2)";

	const char * const pgGMS_BAYESIAN_ESTIMATE_BUILD_DATE = __DATE__ : __TIME__;

	const char * const pgGMS_BAYESIAN_ESTIMATE_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

	const char * const pgGMS_BAYESIAN_ESTIMATE_SYNOPSIS =  " BAYESIAN ESTIMATES OF TIME SERIES MODELS ";
}


#include "GMS_config.h"
#endif
#if defined (__INTEL_COMPILER)
#include <../perf_headers/c++/valarray>
#else
#include <valarray>
#endif
#include <cstdint>

namespace gms {

	namespace math {
		
		  __ATTR_ALIGN__(64) struct BayesianEstimate {

				using VAf64 = std::valarray<double>;
				using VAi32 = std::valarray<int32_t>;

				/*
					Low-level computational code -- members are public.
				*/

				int32_t m_k;

				int32_t m_n;

				int32_t m_il;

				int32_t m_mj2;

				int32_t m_m;

				double m_zmean;

				double m_sum;

				double m_aicm;

				double m_sdm;

				double m_aicb;

				double m_sdb;

				double m_ek;

				double m_oeic;

				double m_omean;

				double m_om;

				VAi32  m_ind;

				VAf64  m_a1;

				VAf64  m_sd;

				VAf64  m_aic;

				VAf64  m_dic;

				VAf64  m_a2;

				VAf64  m_c;

				VAf64  m_c1;

				VAf64  m_c2;

				VAf64  m_b;

				VAf64  m_esum;

				VAf64  m_e;

				VAf64  m_emean;

				VAf64  m_vari;

				VAf64  m_skew;

				VAf64  m_peak;

				VAf64  m_cov;

				VAf64  m_sxx;

				BayesEstimate() = default;

				BayesEstimate(const int32_t ,
					      const int32_t ,
					      const int32_t ,
					      const int32_t ,
					      const int32_t ,
					      const double );

				BayesEstimate(const BayesEstimate &);

				BayesEstimate(BayesEstimate &&) noexcept(true);

				~BayesEstimate() = default;

				BayesEstimate & operator=(const BayesEstimate &);

				BayesEstimate & operator=(BayesEstimate &&) noexcept(true);

				void compute(VAf64 &,
					     int32_t,
					     int32_t,
					     VAi32 &,
					     VAi32 &);
					     


			};
		}
	}




#endif /*_GMS_BAYESIAN_ESTIMATE_H_*/
