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

//
//	Implementation.
//

gms::math::
AVXC2f64::AVXC2f64()
:
m_cv64(_mm256_setzero_pd()) {}

lam::math::
AVXC2f64::AVXC2f64(const double re1,
		   const double im1,
	           const double re2,
	           const double im2)
:
m_cv64(_mm256_setr_pd(re1, im1, re2, im2)) {}

gms::math::
AVXC2f64::AVXC2f64(const double re,
	           const double im)
:
m_cv64(_mm256_setr_pd(re, im, 0.0, 0.0)) {}

gms::math::
AVXC2f64::AVXC2f64(const double re)
:
m_cv64(_mm256_setr_pd(re,0.0,re,0.0)) {}

gms::math::
AVXC2f64::AVXC2f64(const double re1,
	           const double re2,
	           const double re3,
		   const double re4,
		   const bool dummy)
:
m_cv64(_mm256_setr_pd(re1, re2, re3, re4)) {}

gms::math::
AVXC2f64::AVXC2f64(const AVXC2f64 &c)
:
m_cv64(c.m_cv64) {}

gms::math::
AVXC2f64::AVXC2f64(const __m256d &v)
:
m_cv64(v) {}

const __m256d
gms::math::AVXC2f64::get_cv64() const {
	return (this->m_cv64);
}

__m256d
gms::math::AVXC2f64::get_cv64() {
	return (this->m_cv64);
}

__m128d
gms::math::AVXC2f64::complex_1() const {
	return (_mm256_extractf128_pd(m_cv64,0));
}

__m128d
gms::math::AVXC2vf64::complex_2() const {
	return (_mm256_extractf128_pd(m_cv64,1));
}

void
gms::math::AVXC2vf64::set_complex_1(const double re,
			            const double im) {
	m_cv64 = _mm256_setr_pd(re,im,0.0,0.0);
}

void
gms::math::AVXC2vf64::set_complex_2(const double re,
				    const double im) {
	m_cv64 = _mm256_setr_pd(0.0,0.0,re,im);
}

void
gms::math::AVXC2vf64::set_complex_12(const double re1,
				     const double im1,
				     const double re2,
				     const double im2) {
	m_cv64 = _mm256_setr_pd(re1,im1,re2,im2);
}

gms::math::AVXC2f64 &
gms::math::AVXC2f64::load_a(const double * __restrict mem) {
	m_cv64 = _mm256_load_pd(&mem[0]);
	return (*this);
}

gms::math::AVXC2f64 &
gms::math::AVXC2f64::load_u(const double * __restrict mem) {
	m_cv64 = _mm256_loadu_pd(&mem[0]);
	return (*this);
}

void
gms::math::AVXC2f64::store_a(double * __restrict mem) const {
	_mm256_store_pd(&mem[0],m_cv64);
}

void
gms::math::AVXC2f64::store_u(double * __restrict mem) const {
	_mm256_storeu_pd(&mem[0],m_cv64);
}

void
gms::math::AVXC2f64::stream_store(double * __restrict mem) const {
	_mm256_stream_pd(&mem[0],m_cv64);
}

double
gms::math::AVXC2f64::extract_component(const int32_t idx) const {
	__declspec(align(32))double t[4] = {};
	store_a(&t[0]);
	return (t[idx & 3]);
}

gms::math::AVXC2f64  const &
gms::math::AVXC2f64::insert(const int32_t idx,
			    const double val) {
	__m256d ymm = _mm256_broadcast_sd(&val);
	switch (idx) {
		case 0: 
			m_cv64 = _mm256_blend_pd(m_cv64,ymm,1);
	    break;
		case 1:
			m_cv64 = _mm256_blend_pd(m_cv64,ymm,2);
		break;
		case 2:
			m_cv64 = _mm256_blend_pd(m_cv64,ymm,4);
		break;
		default:
			m_cv64 = _mm256_blend_pd(m_cv64,ymm,8);
		
	}
	return (*this);
}

gms::math::AVXC2f64 &
gms::math::AVXC2f64::operator=(const AVXC2f64 &c) {
	if (this == &c) return (*this);
	m_cv64 = c.m_cv64;
	return (*this);
}

gms::math::AVXC2f64::operator __m256d () const {
	return (m_cv64);
}

double
gms::math::AVXC2f64::operator[](const uint32_t idx) const {
	return (reinterpret_cast<const double*>(&m_cv64)[idx]);
}

std::ostream &
lam::math::operator<<(std::ostream &os,
					  _In_ const AVXC2f64 &c) {
	os << std::fixed << std::showpoint << 
		std::setprecision(15) <<
		"Re1: " << "{" << c.get_cv64().m256d_f64[0] << "}" <<
		"Im1: " << "{" << c.get_cv64().m256d_f64[1] << "}" <<
		"Re2: " << "{" << c.get_cv64().m256d_f64[2] << "}" <<
		"Im2: " << "{" << c.get_cv64().m256d_f64[3] << "}" << "\n";
	return (os);
}

__m128d
gms::math::extract(AVXC2f64 &c,
		   const int32_t idx) {
	return (_mm256_extractf128_pd(c,idx));
}

gms::math::AVXC2f64 
gms::math::select(const AVXC2f64 &x,
		 const AVXC2f64 &y,
		 const __m256d &pred) {
	return (_mm256_blendv_pd(x,y,pred));
}

gms::math::AVXC2f64
gms::math::csin(const AVXC2f64 &x) {
	auto c0 = x.complex_1().m128d_f64[0];
	auto c1 = x.complex_1().m128d_f64[1];
	auto c2 = x.complex_2().m128d_f64[0];
	auto c3 = x.complex_2().m128d_f64[1];
	auto t1 = std::sin(c0)*std::cosh(c1);
	auto t2 = std::cos(c0)*std::sinh(c1);
	auto t3 = std::sin(c2)*std::cosh(c3);
	auto t4 = std::cos(c2)*std::sinh(c3);
	return (AVXC2f64{t1,t2,t3,t4});
}

gms::math::AVXC2f64
gms::math::csin(const double re,
		const double im) {
	auto t1 = std::sin(re)*std::cosh(im);
	auto t2 = std::cos(re)*std::sinh(im);
	return (AVXC2f64{t1,t2});
}

gms::math::AVXC2f64
gms::math::ccos(const AVXC2f64 &x) {
	auto c0 = x.complex_1().m128d_f64[0];
	auto c1 = x.complex_1().m128d_f64[1];
	auto c2 = x.complex_2().m128d_f64[0];
	auto c3 = x.complex_2().m128d_f64[1];
	auto t1 = std::cos(c0)*std::cosh(c1);
	auto t2 = std::sin(c0)*std::sinh(c1);
	auto t3 = std::cos(c2)*std::cosh(c3);
	auto t4 = std::sin(c2)*std::sinh(c3);
	return (AVXC2f64{t1,t2,t3,t4});
}

gms::math::AVXC2f64
gms::math::ccos(const double re,
		 const double im) {
	auto t1 = std::cos(re)*std::cosh(im);
	auto t2 = std::sin(re)*std::sinh(im);
	return (AVXC2f64{ t1, t2 });
}

gms::math::AVXC2f64
gms::math::cexp(const AVXC2f64 &x) {
	auto c0 = x.complex_1().m128d_f64[0];
	auto c1 = x.complex_1().m128d_f64[1];
	auto c2 = x.complex_2().m128d_f64[0];
	auto c3 = x.complex_2().m128d_f64[1];
	auto t1 = std::exp(c0)*std::cos(c1);
	auto t2 = std::exp(c0)*std::sin(c1);
	auto t3 = std::exp(c2)*std::cos(c3);
	auto t4 = std::exp(c2)*std::sin(c3);
	return (AVXC2f64{t1,t2,t3,t4});
}

gms::math::AVXC2f64
gms::math::cexp(const double re,
		 const double im) {
	auto t1 = std::exp(re)*std::cos(im);
	auto t2 = std::exp(re)*std::sin(im);
	return (AVXC2f64{t1,t2});
}

std::pair<double,double>
gms::math::cabs(const AVXC2f64 &x) {
	auto c0 = x.complex_1().m128d_f64[0];
	auto c1 = x.complex_1().m128d_f64[1];
	auto c2 = x.complex_2().m128d_f64[0];
	auto c3 = x.complex_2().m128d_f64[1];
	auto t1 = std::sqrt((c0*c0) + (c1*c1));
	auto t2 = std::sqrt((c2*c2) + (c3*c3));
	return (std::make_pair(t1,t2));
}

gms::math::AVXC2f64
gms::math::cpowi(const AVXC2f64 &x,
		const int32_t n) {
	// Slow version because of scalar extraction.
	double r = std::sqrt((x.extract_component(0)*x.extract_component(0)) + 
			   (x.extract_component(1)*x.extract_component(1)));
	double theta = std::atan(x.extract_component(1) / x.extract_component(0));
	return (AVXC2f64{ std::pow(r, n)*std::cos(static_cast<double>(n)*theta),
					   std::pow(r, n)*std::sin(static_cast<double>(n)*theta) });
}

gms::math::AVXC2f64
gms::math::cpowi2(const AVXC2f64  &x,
		 const int32_t n) {
	// Slow version because of scalar extraction.
	double r = std::sqrt((x.extract_component(0)*x.extract_component(0)) +
		(x.extract_component(1)*x.extract_component(1)));
	double theta = std::atan(x.extract_component(1) / x.extract_component(0));
	double r2 = std::sqrt((x.extract_component(2)*x.extract_component(2)) +
		(x.extract_component(3)*x.extract_component(3)));
	double theta2 = std::atan(x.extract_component(3) / x.extract_component(2));
	return (AVXC2vf64{ std::pow(r, n)*std::cos(static_cast<double>(n)*theta),
					   std::pow(r, n)*std::sin(static_cast<double>(n)*theta),
					   std::pow(r2, n)*std::cos(static_cast<double>(n)*theta2),
					   std::pow(r2, n)*std::sin(static_cast<double>(n)*theta2) });
}

double
gms::math::cabs(const double re,
		const double im) {
	return (std::sqrt((re*re) + (im*im)));
}

gms::math::AVXC2f64
gms::math::operator+( const AVXC2f64 x,
		      const AVXC2f64 y) {
	return (_mm256_add_pd(x,y));
}

gms::math::AVXC2f64
gms::math::operator+(const AVXC2vf64 x,
		      const double re) {
	return (x + AVXC2f64{re});
}

gms::math::AVXC2f64
gms::math::operator+(const double re,
		const AVXC2f64 x) {
	return (AVXC2f64{re} + x);
}

gms::math::AVXC2f64
gms::math::operator+=(AVXC2f64 x,
		      const AVXC2f64 y) {
	x = x + y;
	return (x);
}

gms::math::AVXC2f64
gms::math::operator+=(AVXC2f64 x,
		      const double re) {
	x = x + AVXC2f64{re};
	return (x);
}

gms::math::AVXC2f64
gms::math::operator+=(const double re,
                     AVXC2f64 x) {
	x = AVXC2f64{re} + x;
	return (x);
}

gms::math::AVXC2f64
gms::math::operator-(const AVXC2f64 x,
		     const AVXC2vf64 y) {
	return (_mm256_sub_pd(x,y));
}

gms::math::AVXC2f64
gms::math::operator-(const AVXC2f64 x,
		     const double re) {
	return (x - AVXC2f64{re});
}

gms::math::AVXC2f64
gms::math::operator-(const double re,
		const AVXC2f64 x) {
	return (AVXC2f64{re} - x);
}

gms::math::AVXC2f64
gms::math::operator-=(_In_ AVXC2f64 x,
			const AVXC2vf64 y) {
	x = x - y;
	return (x);
}

gms::math::AVXC2f64
gms::math::operator-=(AVXC2vf64 x,
			const double re) {
	x = x - AVXC2f64{re};
	return (x);
}

gms::math::AVXC2f64
gms::math::operator-=(const double re,
			AVXC2f64 x) {
	x = AVXC2f64{re} - x;
	return (x);
}

gms::math::AVXC2f64
gms::math::operator-(AVXC2f64 x) {
	x = x - AVXC2f64{0.0,0.0,0.0,0.0};
	return (x);
}

gms::math::AVXC2f64
gms::math::operator*(const AVXC2f64 x,
			const AVXC2f64 y) {
	// Assume using AVX
	__m256d ymm0 = _mm256_shuffle_pd(y,y,5);
	__m256d ymm1 = _mm256_shuffle_pd(x,x,0xF); // im part
	__m256d ymm2 = _mm256_shuffle_pd(x,x,0x0); // re part
	__m256d ymm3 = _mm256_mul_pd(ymm1,ymm0);
	return (_mm256_fmaddsub_pd(ymm2,y,ymm3));
}

gms::math::AVXC2f64
gms::math::operator*(const AVXC2f64 x,
		const double re) {
	return (x * AVXC2f64{re});
}

gms::math::AVXC2f64
gms::math::operator*(const double re,
			 const AVXC2f64 x) {
	return (AVXC2f64{re} * x);
}

gms::math::AVXC2f64
gms::math::operator*=(AVXC2f64 x,
			 const AVXC2f64 y) {
	x = x * y;
	return (x);
}

gms::math::AVXC2f64
gms::math::operator*=(AVXC2f64 x,
			 const double re) {
	x = x / AVXC2f64{re};
	return (x);
}

gms::math::AVXC2f64
gms::math::operator*=(const double re,
			 AVXC2f64 x) {
	x = AVXC2f64{re} / x;
	return (x);
}

gms::math::AVXC2f64
gms::math::operator/(const AVXC2f64 x,
	             const AVXC2f64 y) {
	// Assumes AVX is available.
	__m256d ymm0 = _mm256_shuffle_pd(x,x,0); // re
	__m256d ymm1 = _mm256_mul_pd(ymm0,y);
	__m256d ymm2 = _mm256_shuffle_pd(y,y,5);
	__m256d ymm3 = _mm256_shuffle_pd(x,x,0xF); // im
	__m256d ymm4 = _mm256_fmsubadd_pd(ymm3,ymm2,ymm0);
	__m256d ymm5 = _mm256_mul_pd(y,y);
	__m256d ymm6 = _mm256_hadd_pd(ymm5,ymm5);
	return (_mm256_div_pd(ymm4,ymm6));
}

gms::math::AVXC2f64
gms::math::operator/( const AVXC2f64 x,
		      const double re) {
	return (x / AVXC2f64{re});
}

gms::math::AVXCvf64
gms::math::operator/(const double re,
		 const AVXC2f64 x) {
	return (AVXC2f64{re} / x);
}

gms::math::AVXC2f64
gms::math::operator/=(AVXC2f64 x,
		      const AVXC2vf64 y) {
	x = x / y;
	return (x);
}

gms::math::AVXC2f64
gms::math::operator/=(AVXC2f64 x,
		      const double re) {
	x = x / AVXC2f64{re};
	return (x);
}

gms::math::AVXC2f64
gms::math::operator/=(const double re,
			 AVXC2vf64 x) {
	x = AVXC2f64{re} / x;
	return (x);
}

gms::math::AVXC2f64
gms::math::operator ~ (_In_ AVXC2f64 x) {
	__m256d t = _mm256_sub_pd(_mm256_setr_pd(0.0,0.0,0.0,0.0),x);
	return (_mm256_mul_pd(_mm256_setr_pd(-1.0,1.0,-1.0,1.0),t));
}

gms::math::AVXC2f64
gms::math::operator==(const AVXC2f64 x,
		      const AVXC2f64 y) {
	return (_mm256_cmp_pd(x,y,_CMP_EQ_OQ));
}

gms::math::AVXC2f64
gms::math::operator!=(const AVXC2vf64 x,
	              const AVXC2vf64 y) {
	return (_mm256_cmp_pd(x,y,_CMP_NEQ_OQ));
}

gms::math::AVXC2f64
gms::math::operator>(const AVXC2f64 x,
		    const AVXC2f64 y) {
	return (_mm256_cmp_pd(x,y,_CMP_GT_OQ));
}

gms::math::AVXC2f64
gms::math::operator>=(const AVXC2f64 x,
		      const AVXC2f64 y) {
	return (_mm256_cmp_pd(x,y,_CMP_GE_OQ));
}

gms::math::AVXC2f64
gms::math::operator<(const AVXC2f64 x,
		const AVXC2f64 y) {
	return (_mm256_cmp_pd(x,y,_CMP_LT_OQ));
}

gms::math::AVXC2f64
gms::math::operator<=(const AVXC2f64 x,
		     const AVXC2f64 y) {
	return (_mm256_cmp_pd(x,y,_CMP_LE_OQ));
}






