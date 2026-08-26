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
#include "GMS_avx512complex.h"
#include "GMS_avx512vcprocedures.h"
#include "GMS_common.h"
#include "GMS_config.h"

//
//	Implementation
//

#if (GMS_MAN_PREFETCH) == 1

#if !defined (AVX512VCPROCEDURES_PREFETCH_FROM_OBJS)
#define AVX512VCPROCEDURES_PREFETCH_FROM_OBJS(obj1,obj2,idx,off,hint) \
	_mm_prefetch(reinterpret_cast<const char*>(&(obj1).data.m_Re[(idx)+(off)]), (hint)); \
	_mm_prefetch(reinterpret_cast<const char*>(&(obj2).data.m_Re[(idx)+(off)]), (hint)); \
	_mm_prefetch(reinterpret_cast<const char*>(&(obj1).data.m_Im[(idx)+(off)]), (hint)); \
	_mm_prefetch(reinterpret_cast<const char*>(&(obj2).data.m_Im[(idx)+(off)]), (hint));
#endif

#if !defined (AVX512VCPROCEDURES_PREFETCH_FROM_OBJ_PTR)
#define AVX512VCPROCEDURES_PREFETCH_FROM_OBJ_PTR(obj,ptr,idx,off,hint)  \
	_mm_prefetch(reinterpret_cast<const char*>(&(obj).data.m_Re[(idx)+(off)]), (hint)); \
	_mm_prefetch(reinterpret_cast<const char*>(&(ptr)[(idx)+(off)]), (hint));	   \
	_mm_prefetch(reinterpret_cast<const char*>(&(obj).data.m_Im[(idx)+(off)]), (hint));
#endif	

#endif

void
gms::math::avx512vcomplex_add(AVX512VComplex1D &c,
			      const AVX512VComplex1D &a,
			      const AVX512VComplex1D &b,
			      const bool do_stream_store) {
	if(a.data.m_nsize != b.data.m_nsize || 
	   c.data.m_nsize != a.data.m_nsize) {
	   return;
	   }
	int32_t i;
	if (do_stream_store) {

	     for (i = 0; i != ROUND_TO_EIGHT(c.data.m_nsize, 8); i += 16) {
		       
		          const __m512d zmm0(_mm512_load_pd(&a.data.m_Re[i+0]));
		          const __m512d zmm1(_mm512_load_pd(&b.data.m_Re[i+0]));
				  _mm512_stream_pd(&c.data.m_Re[i + 0], _mm512_add_pd(zmm0,zmm1));
				  const __m512d zmm2(_mm512_load_pd(&a.data.m_Re[i+8]));
				  const __m512d zmm3(_mm512_load_pd(&b.data.m_Re[i+8]));
				  _mm512_stream_pd(&c.data.m_Re[i + 8], _mm512_add_pd(zmm2,zmm3));
				  const __m512d zmm4(_mm512_load_pd(&a.data.m_Im[i+0]));
				  const __m512d zmm5(_mm512_load_pd(&b.data.m_Im[i+0]));
				  _mm512_stream_pd(&c.data.m_Im[i + 0], _mm512_add_pd(zmm4,zmm5));
				  const __m512d zmm6(_mm512_load_pd(&a.data.m_Im[i+8]));
				  const __m512d zmm7(_mm512_load_pd(&b.data.m_Im[i+8]));
				  _mm512_stream_pd(&c.data.m_Im[i + 8], _mm512_add_pd(zmm6,zmm7));
		
		}
		 for (; i != c.data.m_nsize; ++i) {
			 c.data.m_Re[i] = a.data.m_Re[i] + b.data.m_Re[i];
			 c.data.m_Im[i] = a.data.m_Im[i] + b.data.m_Im[i];
		 }
	
	}
	else {
			
		for (i = 0; i != ROUND_TO_EIGHT(c.data.m_nsize, 8); i += 16) {

			const __m512d zmm0(_mm512_load_pd(&a.data.m_Re[i + 0]));
			const __m512d zmm1(_mm512_load_pd(&b.data.m_Re[i + 0]));
			_mm512_store_pd(&c.data.m_Re[i + 0], _mm512_add_pd(zmm0, zmm1));
			const __m512d zmm2(_mm512_load_pd(&a.data.m_Re[i + 8]));
			const __m512d zmm3(_mm512_load_pd(&b.data.m_Re[i + 8]));
			_mm512_store_pd(&c.data.m_Re[i + 8], _mm512_add_pd(zmm2, zmm3));
			const __m512d zmm4(_mm512_load_pd(&a.data.m_Im[i + 0]));
			const __m512d zmm5(_mm512_load_pd(&b.data.m_Im[i + 0]));
			_mm512_store_pd(&c.data.m_Im[i + 0], _mm512_add_pd(zmm4, zmm5));
			const __m512d zmm6(_mm512_load_pd(&a.data.m_Im[i + 8]));
			const __m512d zmm7(_mm512_load_pd(&b.data.m_Im[i + 8]));
			_mm512_store_pd(&c.data.m_Im[i + 8], _mm512_add_pd(zmm6, zmm7));

		}
		for (; i != c.data.m_nsize; ++i) {
			c.data.m_Re[i] = a.data.m_Re[i] + b.data.m_Re[i];
			c.data.m_Im[i] = a.data.m_Im[i] + b.data.m_Im[i];
		}
	}
}

void
gms::math::avx512vcomplex_add(AVX512VComplex1D &c,
			      const AVX512VComplex1D &a,
			      const double * __restrict b,
			      const bool do_stream_store) {
	using namespace gms::common;
	if (!Is_ptr_aligned64(b) || 
	    c.data.m_nsize != a.data.m_nsize) {
		return;
	}
	int32_t i;
	if (do_stream_store) {

			for (i = 0; i != ROUND_TO_EIGHT(c.data.m_nsize, 8); i += 16) {
				 
				const __m512d zmm0(_mm512_load_pd(&a.data.m_Re[i+0]));
				const __m512d zmm1(_mm512_load_pd(&b[i+0]));
				_mm512_stream_pd(&c.data.m_Re[i+0], _mm512_add_pd(zmm0,zmm1));
				const __m512d zmm2(_mm512_load_pd(&a.data.m_Re[i+8]));
				const __m512d zmm3(_mm512_load_pd(&b[i+0]));
				_mm512_stream_pd(&c.data.m_Re[i + 8], _mm512_add_pd(zmm2,zmm3));

		}	
			for (; i != c.data.m_nsize; ++i){
	   			c.data.m_Re[i] = a.data.m_Re[i] + b[i];
	    }
	}
	else {

		for (i = 0; i != ROUND_TO_EIGHT(c.data.m_nsize, 8); i += 16) {

			const __m512d zmm0(_mm512_load_pd(&a.data.m_Re[i + 0]));
			const __m512d zmm1(_mm512_load_pd(&b[i + 0]));
			_mm512_store_pd(&c.data.m_Re[i + 0], _mm512_add_pd(zmm0, zmm1));
			const __m512d zmm2(_mm512_load_pd(&a.data.m_Re[i + 8]));
			const __m512d zmm3(_mm512_load_pd(&b[i + 0]));
			_mm512_store_pd(&c.data.m_Re[i + 8], _mm512_add_pd(zmm2, zmm3));

		}
		for (; i != c.data.m_nsize; ++i){
			c.data.m_Re[i] = a.data.m_Re[i] + b[i];
		}
	}
}	

void
gms::math::avx512vcomplex_sub(AVX512VComplex1D &c,
			      const AVX512VComplex1D &a,
			      const AVX512VComplex1D &b,
			      const bool do_stream_store) {
	if (a.data.m_nsize != b.data.m_nsize ||
		c.data.m_nsize != a.data.m_nsize) {
		return;
	}
	int32_t i;
	if (do_stream_store) {

			for (i = 0; i != ROUND_TO_FOUR(c.data.m_nsize, 8); i += 16) {
			// Number of cache misses ~ m_nsize / 8
				const __m512d zmm0(_mm512_load_pd(&a.data.m_Re[i+0]));
				const __m512d zmm1(_mm512_load_pd(&b.data.m_Re[i+0]));
				_mm512_stream_pd(&c.data.m_Re[i + 0], _mm512_sub_pd(zmm0,zmm1));
				const __m512d zmm2(_mm512_load_pd(&a.data.m_Re[i+8]));
				const __m512d zmm3(_mm512_load_pd(&b.data.m_Re[i+8]));
				_mm512_stream_pd(&c.data.m_Re[i + 8], _mm512_sub_pd(zmm2,zmm3));
				const __m512d zmm4(_mm512_load_pd(&a.data.m_Im[i+0]));
				const __m512d zmm5(_mm512_load_pd(&b.data.m_Im[i+0]));
				_mm512_stream_pd(&c.data.m_Im[i + 0], _mm512_sub_pd(zmm4,zmm5));
				const __m512d zmm6(_mm512_load_pd(&a.data.m_Im[i+8]));
				const __m512d zmm7(_mm512_load_pd(&b.data.m_Im[i+8]));
				_mm512_stream_pd(&c.data.m_Im[i + 8], _mm512_sub_pd(zmm6,zmm7));
		
		}
			// This remainder will consume 1 cache line
			for (; i != c.data.m_nsize; ++i) {
				c.data.m_Re[i] = a.data.m_Re[i] - b.data.m_Re[i];
				c.data.m_Im[i] = a.data.m_Im[i] - b.data.m_Im[i];
			}
	}
	else {
		
		for (i = 0; i != ROUND_TO_FOUR(c.data.m_nsize, 8); i += 16) {
			// Number of cache misses ~ m_nsize / 8
			const __m512d zmm0(_mm512_load_pd(&a.data.m_Re[i + 0]));
			const __m512d zmm1(_mm512_load_pd(&b.data.m_Re[i + 0]));
			_mm512_store_pd(&c.data.m_Re[i + 0], _mm512_sub_pd(zmm0, zmm1));
			const __m512d zmm2(_mm512_load_pd(&a.data.m_Re[i + 8]));
			const __m512d zmm3(_mm512_load_pd(&b.data.m_Re[i + 8]));
			_mm512_store_pd(&c.data.m_Re[i + 8], _mm512_sub_pd(zmm2, zmm3));
			const __m512d zmm4(_mm512_load_pd(&a.data.m_Im[i + 0]));
			const __m512d zmm5(_mm512_load_pd(&b.data.m_Im[i + 0]));
			_mm512_store_pd(&c.data.m_Im[i + 0], _mm512_sub_pd(zmm4, zmm5));
			const __m512d zmm6(_mm512_load_pd(&a.data.m_Im[i + 8]));
			const __m512d zmm7(_mm512_load_pd(&b.data.m_Im[i + 8]));
			_mm512_store_pd(&c.data.m_Im[i + 8], _mm512_sub_pd(zmm6, zmm7));

		}
		
		for (; i != c.data.m_nsize; ++i) {
			c.data.m_Re[i] = a.data.m_Re[i] - b.data.m_Re[i];
			c.data.m_Im[i] = a.data.m_Im[i] - b.data.m_Im[i];
		}
	}
}	

	


void
gms::math::avx512vcomplex_sub(AVX512VComplex1D &c,
			      const AVX512VComplex1D &a,
			      const double * __restrict b,
			      const bool do_stream_store) {
	using namespace gms::common;
	if (!Is_ptr_aligned64(b) ||
		c.data.m_nsize != a.data.m_nsize) {
		return;
	}
	int32_t i;
	if (do_stream_store) {

			for (i = 0; i != ROUND_TO_EIGHT(c.data.m_nsize, 8); i += 16) {

				const __m512d zmm0(_mm512_load_pd(&a.data.m_Re[i+0]));
				const __m512d zmm1(_mm512_load_pd(&b[i+0]));
				_mm512_stream_pd(&c.data.m_Re[i + 0], _mm512_sub_pd(zmm0,zmm1));
				const __m512d zmm2(_mm512_load_pd(&a.data.m_Re[i+8]));
				const __m512d zmm3(_mm512_load_pd(&b[i+8]));
				_mm512_stream_pd(&c.data.m_Re[i + 8], _mm512_sub_pd(zmm2,zmm3));
		
		}
			for (; i != c.data.m_nsize; ++i){
				c.data.m_Re[i] = a.data.m_Re[i] - b[i];
			}
	}
	else {
		
		for (i = 0; i != ROUND_TO_EIGHT(c.data.m_nsize, 8); i += 16) {

			const __m512d zmm0(_mm512_load_pd(&a.data.m_Re[i + 0]));
			const __m512d zmm1(_mm512_load_pd(&b[i + 0]));
			_mm512_store_pd(&c.data.m_Re[i + 0], _mm512_sub_pd(zmm0, zmm1));
			const __m512d zmm2(_mm512_load_pd(&a.data.m_Re[i + 8]));
			const __m512d zmm3(_mm512_load_pd(&b[i + 8]));
			_mm512_store_pd(&c.data.m_Re[i + 8], _mm512_sub_pd(zmm2, zmm3));

		}
		for (; i != c.data.m_nsize; ++i){
			c.data.m_Re[i] = a.data.m_Re[i] - b[i];
		}
	}
}	


void
gms::math::avx512vcomplex_mul(AVX512VComplex1D &c,
			      const AVX512VComplex1D &a,
			      const AVX512VComplex1D &b,
			      const bool do_stream_store) {
	if (a.data.m_nsize != b.data.m_nsize ||
	    c.data.m_nsize != a.data.m_nsize) {
		return;
	}
	int32_t i;
	if (do_stream_store) {

		for (i = 0; i != ROUND_TO_EIGHT(c.data.m_nsize, 8); i += 16) {

			const __m512d zmm0(_mm512_load_pd(&a.data.m_Re[i+0]));
			const __m512d zmm1(_mm512_load_pd(&b.data.m_Re[i+0]));
			const __m512d zmm2(_mm512_load_pd(&a.data.m_Im[i+0]));
			const __m512d zmm3(_mm512_load_pd(&b.data.m_Im[i+0]));
			_mm512_stream_pd(&c.data.m_Re[i + 0], _mm512_sub_pd(
				_mm512_mul_pd(zmm0, zmm1), _mm512_mul_pd(zmm2, zmm3)));
			_mm512_stream_pd(&c.data.m_Im[i + 0], _mm512_add_pd(
				_mm512_mul_pd(zmm2, zmm1), _mm512_mul_pd(zmm0, zmm3)));
			const __m512d zmm4(_mm512_load_pd(&a.data.m_Re[i+8]));
			const __m512d zmm5(_mm512_load_pd(&b.data.m_Re[i+8]));
			const __m512d zmm6(_mm512_load_pd(&a.data.m_Im[i+8]));
			const __m512d zmm7(_mm512_load_pd(&b.data.m_Im[i+8]));
			_mm512_stream_pd(&c.data.m_Re[i + 8], _mm512_sub_pd(
				_mm512_mul_pd(zmm4, zmm5), _mm512_mul_pd(zmm6, zmm7)));
			_mm512_stream_pd(&c.data.m_Im[i + 8], _mm512_add_pd(
				_mm512_mul_pd(zmm6, zmm5), _mm512_mul_pd(zmm4, zmm7)));

		}
		for (; i != c.data.m_nsize; ++i) {
			c.data.m_Re[i] = (a.data.m_Re[i] * b.data.m_Re[i]) - (a.data.m_Im[i] * b.data.m_Im[i]);
			c.data.m_Im[i] = (a.data.m_Im[i] * b.data.m_Re[i]) + (a.data.m_Re[i] * b.data.m_Im[i]);
		}
	}
	else {
			
		for (i = 0; i != ROUND_TO_EIGHT(c.data.m_nsize, 8); i += 16) {

			const __m512d zmm0(_mm512_load_pd(&a.data.m_Re[i + 0]));
			const __m512d zmm1(_mm512_load_pd(&b.data.m_Re[i + 0]));
			const __m512d zmm2(_mm512_load_pd(&a.data.m_Im[i + 0]));
			const __m512d zmm3(_mm512_load_pd(&b.data.m_Im[i + 0]));
			_mm512_store_pd(&c.data.m_Re[i + 0], _mm512_sub_pd(
				_mm512_mul_pd(zmm0, zmm1), _mm512_mul_pd(zmm2, zmm3)));
			_mm512_store_pd(&c.data.m_Im[i + 0], _mm512_add_pd(
				_mm512_mul_pd(zmm2, zmm1), _mm512_mul_pd(zmm0, zmm3)));
			const __m512d zmm4(_mm512_load_pd(&a.data.m_Re[i + 8]));
			const __m512d zmm5(_mm512_load_pd(&b.data.m_Re[i + 8]));
			const __m512d zmm6(_mm512_load_pd(&a.data.m_Im[i + 8]));
			const __m512d zmm7(_mm512_load_pd(&b.data.m_Im[i + 8]));
			_mm512_store_pd(&c.data.m_Re[i + 8], _mm512_sub_pd(
				_mm512_mul_pd(zmm4, zmm5), _mm512_mul_pd(zmm6, zmm7)));
			_mm512_store_pd(&c.data.m_Im[i + 8], _mm512_add_pd(
				_mm512_mul_pd(zmm6, zmm5), _mm512_mul_pd(zmm4, zmm7)));

		}
		for (; i != c.data.m_nsize; ++i) {
			c.data.m_Re[i] = (a.data.m_Re[i] * b.data.m_Re[i]) - (a.data.m_Im[i] * b.data.m_Im[i]);
			c.data.m_Im[i] = (a.data.m_Im[i] * b.data.m_Re[i]) + (a.data.m_Re[i] * b.data.m_Im[i]);
		}
	}
	
	
	
}

void
gms::math::avx512vcomplex_mul(AVX512VComplex1D &c,
			      const AVX512VComplex1D &a,
			      const double * __restrict b,
			      const bool do_stream_store) {
	using namespace gms::common;
	if (!Is_ptr_aligned64(b) ||
		c.data.m_nsize != a.data.m_nsize) {
		return;
	}
	int32_t i;
	if (do_stream_store) {

			for (i = 0; i != ROUND_TO_EIGHT(c.data.m_nsize, 8); i += 16) {

				const __m512d zmm0(_mm512_load_pd(&a.data.m_Re[i+0]));
				const __m512d zmm1(_mm512_load_pd(&b[i+0]));
				_mm512_stream_pd(&c.data.m_Re[i + 0], _mm512_mul_pd(zmm0, zmm1));
				const __m512d zmm2(_mm512_load_pd(&a.data.m_Re[i+8]));
				const __m512d zmm3(_mm512_load_pd(&b[i+8]));
				_mm512_stream_pd(&c.data.m_Re[i + 8], _mm512_mul_pd(zmm2, zmm3));
				const __m512d zmm4(_mm512_load_pd(&a.data.m_Im[i+0]));
				_mm512_stream_pd(&c.data.m_Im[i + 0], _mm512_mul_pd(zmm4, zmm1));
				const __m512d zmm5(_mm512_load_pd(&a.data.m_Im[i+8]));
				_mm512_stream_pd(&c.data.m_Im[i + 8], _mm512_mul_pd(zmm5, zmm3));
		}
		 for (; i != c.data.m_nsize; ++i) {
				c.data.m_Re[i] = a.data.m_Re[i] * b[i];
				c.data.m_Im[i] = a.data.m_Im[i] * b[i];
			}
	}
	else {
		
		for (i = 0; i != ROUND_TO_EIGHT(c.data.m_nsize, 8); i += 16) {

			const __m512d zmm0(_mm512_load_pd(&a.data.m_Re[i + 0]));
			const __m512d zmm1(_mm512_load_pd(&b[i + 0]));
			_mm512_store_pd(&c.data.m_Re[i + 0], _mm512_mul_pd(zmm0, zmm1));
			const __m512d zmm2(_mm512_load_pd(&a.data.m_Re[i + 8]));
			const __m512d zmm3(_mm512_load_pd(&b[i + 8]));
			_mm512_store_pd(&c.data.m_Re[i + 8], _mm512_mul_pd(zmm2, zmm3));
			const __m512d zmm4(_mm512_load_pd(&a.data.m_Im[i + 0]));
			_mm512_store_pd(&c.data.m_Im[i + 0], _mm512_mul_pd(zmm4, zmm1));
			const __m512d zmm5(_mm512_load_pd(&a.data.m_Im[i + 8]));
			_mm512_store_pd(&c.data.m_Im[i + 8], _mm512_mul_pd(zmm5, zmm3));
		}
		for (; i != c.data.m_nsize; ++i) {
			c.data.m_Re[i] = a.data.m_Re[i] * b[i];
			c.data.m_Im[i] = a.data.m_Im[i] * b[i];
		}
	}
}

void
gms::math::avx512vcomplex_div(AVX512VComplex1D &c,
			      const AVX512VComplex1D &a,
			      const AVX512VComplex1D &b,
			      const bool do_stream_store) {
	if (a.data.m_nsize != b.data.m_nsize ||
		c.data.m_nsize != a.data.m_nsize) {
		return;
	}
	int32_t i;
	if (do_stream_store) {

			for (i = 0; i != ROUND_TO_EIGHT(c.data.m_nsize, 8); i += 8) {
					// Will unrolling 2x not saturate divider unit.
					// We have two parallel division so at least second
					// operation will be pipelined at divider level.
				const __m512d zmm0(_mm512_load_pd(&a.data.m_Re[i+0]));
				const __m512d zmm1(_mm512_load_pd(&b.data.m_Im[i+0]));
				const __m512d zmm2(_mm512_load_pd(&a.data.m_Im[i+0]));
				const __m512d re_term1 = _mm512_add_pd(
					_mm512_mul_pd(zmm0, zmm1), _mm512_mul_pd(zmm2, zmm1));
				const __m512d re_term2 = _mm512_add_pd(
					_mm512_mul_pd(zmm2, zmm1), _mm512_mul_pd(zmm0, zmm1));
				const __m512d zmm3(_mm512_load_pd(&b.data.m_Re[i+0]));
				const __m512d den_term = _mm512_add_pd(
					_mm512_mul_pd(zmm3, zmm3), _mm512_mul_pd(zmm1, zmm1));
				_mm512_stream_pd(&c.data.m_Re[i + 0], _mm512_div_pd(re_term1, den_term));
				_mm512_stream_pd(&c.data.m_Im[i + 0], _mm512_div_pd(re_term2, den_term));

			}
			for (; i != c.data.m_nsize; ++i) {
					const double re =  (a.data.m_Re[i] * b.data.m_Im[i]) + (a.data.m_Im[i] * b.data.m_Im[i]);
					const double im =  (a.data.m_Im[i] * b.data.m_Im[i]) - (a.data.m_Re[i] * b.data.m_Im[i]);
					const double den = (b.data.m_Re[i] * b.data.m_Re[i]) + (b.data.m_Im[i] * b.data.m_Im[i]);
					c.data.m_Re[i] = re / den;
					c.data.m_Im[i] = im / den;
				}
	}
	else {
		
		for (i = 0; i != ROUND_TO_EIGHT(c.data.m_nsize, 8); i += 8) {
			// Will unrolling 2x not saturate divider unit.
			// We have two parallel division so at least second
			// operation will be pipelined at divider level.
			const __m512d zmm0(_mm512_load_pd(&a.data.m_Re[i + 0]));
			const __m512d zmm1(_mm512_load_pd(&b.data.m_Im[i + 0]));
			const __m512d zmm2(_mm512_load_pd(&a.data.m_Im[i + 0]));
			const __m512d re_term1 = _mm512_add_pd(
				_mm512_mul_pd(zmm0, zmm1), _mm512_mul_pd(zmm2, zmm1));
			const __m512d re_term2 = _mm512_add_pd(
				_mm512_mul_pd(zmm2, zmm1), _mm512_mul_pd(zmm0, zmm1));
			const __m512d zmm3(_mm512_load_pd(&b.data.m_Re[i + 0]));
			const __m512d den_term = _mm512_add_pd(
				_mm512_mul_pd(zmm3, zmm3), _mm512_mul_pd(zmm1, zmm1));
			_mm512_store_pd(&c.data.m_Re[i + 0], _mm512_div_pd(re_term1, den_term));
			_mm512_store_pd(&c.data.m_Im[i + 0], _mm512_div_pd(re_term2, den_term));

		}
		for (; i != c.data.m_nsize; ++i) {
			const double re = (a.data.m_Re[i] * b.data.m_Im[i]) + (a.data.m_Im[i] * b.data.m_Im[i]);
			const double im = (a.data.m_Im[i] * b.data.m_Im[i]) - (a.data.m_Re[i] * b.data.m_Im[i]);
			const double den = (b.data.m_Re[i] * b.data.m_Re[i]) + (b.data.m_Im[i] * b.data.m_Im[i]);
			c.data.m_Re[i] = re / den;
			c.data.m_Im[i] = im / den;
		}
	}
}

void
gms::math::avx512vcomplex_div(AVX512VComplex1D &c,
			      const AVX512VComplex1D &a,
			      const double * __restrict b,
			      const bool do_stream_store) {
	using namespace gms::common;
	if (!Is_ptr_aligned64(b) ||
		c.data.m_nsize != a.data.m_nsize) {
		return;
	}
	int32_t i;
	if (do_stream_store) {

			for (i = 0; i != ROUND_TO_EIGHT(c.data.m_nsize, 8); i += 8) {

				const __m512d zmm0(_mm512_load_pd(&a.data.m_Re[i+0]));
				const __m512d zmm1(_mm512_load_pd(&b[i+0]));
				_mm512_stream_pd(&c.data.m_Re[i], _mm512_div_pd(zmm0, zmm1));
				const __m512d zmm2(_mm512_load_pd(&a.data.m_Im[i+0]));
				_mm512_stream_pd(&c.data.m_Im[i], _mm512_div_pd(zmm2, zmm1));
	     }
	      for (; i != c.data.m_nsize; ++i) {
		       c.data.m_Re[i] = a.data.m_Re[i] / b[i];
		       c.data.m_Im[i] = a.data.m_Im[i] / b[i];
	     }
	}
	else {
		
		for (i = 0; i != ROUND_TO_EIGHT(c.data.m_nsize, 8); i += 8) {

			const __m512d zmm0(_mm512_load_pd(&a.data.m_Re[i + 0]));
			const __m512d zmm1(_mm512_load_pd(&b[i + 0]));
			_mm512_store_pd(&c.data.m_Re[i], _mm512_div_pd(zmm0, zmm1));
			const __m512d zmm2(_mm512_load_pd(&a.data.m_Im[i + 0]));
			_mm512_store_pd(&c.data.m_Im[i], _mm512_div_pd(zmm2, zmm1));
		}
		for (; i != c.data.m_nsize; ++i) {
			c.data.m_Re[i] = a.data.m_Re[i] / b[i];
			c.data.m_Im[i] = a.data.m_Im[i] / b[i];
		}
	}
}


