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


#include "GMS_mkl_cauchyrng.h"
#include "GMS_malloc.h"
#include "GMS_simd_memops.h" 
//
//	Implementation
//

gms::math::
MKLCauchyRNG::MKLCauchyRNG() {
	
	m_rvec    = NULL;
	m_a       = 0.0;
	m_beta    = 0.0;
	m_nvalues = 0;
	m_brng    = 0;
	m_seed    = 0;
	m_error   = 1;
}


gms::math::
MKLCauchyRNG::MKLCauchyRNG(const MKL_INT nvalues,
			   const MKL_UINT brng,
			   const MKL_INT seed,
			   const double a,
			   const double beta) {
	using namespace gms::common;
        constexpr size_t align64B = 64ULL;
	m_rvec    = gms_mm_malloc(static_cast<size_t>(nvalues), align64B);

	m_a       = a;
	m_beta    = beta;
	m_nvalues = nvalues;
	m_brng    = brng;
	m_seed    = seed;
	m_error   = 1;
#if defined __AVX512F__
        avx512_init_unroll8x_pd(&m_rvec[0], static_cast<int64_t>(m_nvalues), 0.0);	
#else
	avx256_init_unroll8x_pd(&m_rvec[0], static_cast<int64_t>(m_nvalues), 0.0);
#endif
}


gms::math::
MKLCauchyRNG::MKLCauchyRNG(const MKLCauchyRNG &x) {
	using namespace gms::common;
         constexpr size_t align64B = 64ULL;
	m_rvec    = gms_mm_malloc(static_cast<size_t>(x.m_nvalues), align64B);

	m_a       = x.m_a;
	m_beta    = x.m_beta;
	m_nvalues = x.m_nvalues;
	m_brng    = x.m_brng;
	m_seed    = x.m_seed;
	m_error   = x.m_error;
#if defined __AVX512F__
    #if (USE_NT_STORES) == 1
	    avx512_uncached_memmove(&m_rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
    #else
	    avx512_cached_memmove(&m_rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
    #endif
#else
    #if (USE_NT_STORES) == 1
	     avx256_uncached_memmove(&m_rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
    #else
	     avx256_cached_memmove(&m_rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
    #endif
#endif
}


gms::math::
MKLCauchyRNG::MKLCauchyRNG(MKLCauchyRNG &&x) {
	m_rvec      = &x.m_rvec[0];
	m_a         = x.m_a;
	m_beta      = x.m_beta;
	m_nvalues   = x.m_nvalues;
	m_brng      = x.m_brng;
	m_seed      = x.m_seed;
	m_error     = x.m_error;
	x.m_rvec    = NULL;
	x.m_nvalues = 0;
}


gms::math::
MKLCauchyRNG::~MKLCauchyRNG() {
         using namespace gms::common;
	if (NULL != m_rvec) gms_mm_free(m_rvec); m_rvec = NULL;
}		
		
	


gms::math::MKLCauchyRNG &
gms::math::MKLCauchyRNG
::operator=(const MKLCauchyRNG &x) {
	using namespace gms::common;
	if (this == &x) return (*this);
	gms_mm_free(m_rvec);
	m_a       = x.m_a;
	m_beta    = x.m_beta;
	m_nvalues = x.m_nvalues;
	m_brng    = x.m_brng;
	m_seed    = x.m_seed;
	m_error   = x.m_error;
        constexpr size_t align64B = 64ULL;
	double * __restrict
	        rvec{gms_edmalloca(static_cast<size_t>(m_nvalues), align64B)};

#if defined __AVX512F__
    #if (USE_NT_STORES) == 1
	avx512_uncached_memmove(&rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
    #else
	avx512_cached_memmove(&rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
    #endif
#else
    #if (USE_NT_STORES) == 1
	    avx256_uncached_memmove(&rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
    #else
	    avx256_cached_memmove(&rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
    #endif
#endif
	m_rvec = &rvec[0];
	return (*this);
}	
	


gms::math::MKLCauchyRNG &
gms::math::MKLCauchyRNG
::operator=(MKLCauchyRNG &&x) {
        using namespace gms::common; 
	if (this == &x) return (*this);
	gms_mm_free(m_rvec);
	m_rvec = &x.m_rvec[0];
	m_a = x.m_a;
	m_beta = x.m_beta;
	m_nvalues = x.m_nvalues;
	m_brng = x.m_brng;
	m_seed = x.m_seed;
	m_error = x.m_error;
	x.m_rvec = NULL;
	x.m_nvalues = 0;
	return (*this);
}	
	


void
gms::math::
MKLCauchyRNG::compute_rand_distribution(const MKL_INT method) {
	VSLStreamStatePtr stream;
	m_error = vslNewStream(&stream,m_brng,m_seed);
	if (VSL_ERROR_OK != m_error) {
		PRINT_ERROR_VALUE("vslNewStream -- failed with an error value: ", m_error)
		DUMP_CALLSTACK_ON_ERROR
		return;
	}
	m_error = vdRngCauchy(method,stream,m_nvalues,&m_rvec[0],m_a,m_beta);
	if (VSL_ERROR_OK != m_error) {
		if (VSL_ERROR_OK != m_error) {
			PRINT_ERROR_VALUE("vdRngCauchy -- failed with an error value: ", m_error)
			DUMP_CALLSTACK_ON_ERROR
			return;
		}
	}
	m_error = vslDeleteStream(&stream);
	if (VSL_ERROR_OK != m_error) {
		PRINT_ERROR_VALUE("vslDeleteStream -- failed with an error value: ", m_error)
		DUMP_CALLSTACK_ON_ERROR
		return;
	}
}

void
gms::math::
MKLCauchyRNG::
compute_rand_distribution(VSLStreamStatePtr stream,
			  const MKL_INT method) {
	// VSLStreamStatePtr must deallocated by the caller of compute_rand_distribution.
	m_error = vslNewStream(&stream,m_brng,m_seed);
	if (VSL_ERROR_OK != m_error) {
		PRINT_ERROR_VALUE("vslNewStream -- failed with an error value: ", m_error)
		DUMP_CALLSTACK_ON_ERROR
		return;
	}
	m_error = vdRngCauchy(method, stream, m_nvalues, &m_rvec[0],m_a,m_beta);
	if (VSL_ERROR_OK != m_error) {
		if (VSL_ERROR_OK != m_error) {
			PRINT_ERROR_VALUE("vdRngCauchy -- failed with an error value: ", m_error)
			DUMP_CALLSTACK_ON_ERROR
			return;
		}
	}
	// stream is deleted by the caller.
}

std::ostream &
gms::math::operator<<(std::ostream &os,
			    const MKLCauchyRNG &x) {
	
	for (MKL_INT i = 0; i != x.m_nvalues; ++i) {
		os << std::setprecision(15) << std::showpoint <<
			"Randomly Cauchy distributed: m_rvec = " << x.m_rvec[i] << std::endl;
	}
	return (os);
}
