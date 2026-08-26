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
#include <iostream>
#include <iomanip>
#include "GMS_chrono_timer.h"


#include "GMS_indices.h"
#include "GMS_constants.h"
#include "GMS_common.h"


gms::system
::HRCTimer
::HRCTimer(const std::size_t nruns,
	   const std::size_t datum_size,
	   const fptr ptfunc,
	   const fptr ptwrap,
	   const char * pfuncname,
	   const char * pfilename)
:
m_nruns{ nruns },
m_datum_size{datum_size},
m_ptfunc(ptfunc),
m_ptwrap(ptwrap),
m_func_name{ pfuncname },
m_file_name{ pfilename},
m_timing_values(0.0,m_nruns * m_datum_size),
m_overhead_values(0.0, m_nruns * m_datum_size),
m_delta_values(0.0, m_nruns * m_datum_size)
{}

gms::system
::HRCTimer
::HRCTimer(const HRCTimer &x)
:
m_nruns{ x.m_nruns },
m_datum_size{ x.m_datum_size },
m_ptfunc(x.m_ptfunc),
m_ptwrap(x.m_ptwrap),
m_func_name{ x.m_func_name },
m_file_name{ x.m_file_name },
m_timing_values{ x.m_timing_values },
m_overhead_values{ x.m_overhead_values },
m_delta_values{ x.m_delta_values }
{}

gms::system
::HRCTimer
::HRCTimer(HRCTimer &&x)
:
m_nruns{ x.m_nruns },
m_datum_size{ x.m_datum_size },
m_ptfunc(x.m_ptfunc),
m_ptwrap(x.m_ptwrap),
m_func_name{ std::move(x.m_func_name) },
m_file_name{ std::move(x.m_file_name) },
m_timing_values{ std::move(x.m_timing_values) },
m_overhead_values{ std::move(x.m_overhead_values) },
m_delta_values{ std::move(x.m_delta_values) }
{}

gms::system
::HRCTimer
::~HRCTimer() {
	if (m_ptfunc != nullptr) m_ptfunc = nullptr;
	if (m_ptwrap != nullptr) m_ptwrap = nullptr;
}

gms::system::HRCTimer &
gms::system::HRCTimer
::operator=(const HRCTimer &x) {
	if (this == &x) return (*this);
			m_nruns		 = x.m_nruns;
			m_datum_size = x.m_datum_size;
			m_ptfunc     = x.m_ptfunc;
			m_ptwrap     = x.m_ptwrap;
			m_func_name.operator=(x.m_func_name);
			m_file_name.operator=(x.m_file_name);
			m_timing_values.operator=(x.m_timing_values);
			m_overhead_values.operator=(x.m_overhead_values);
			m_delta_values.operator=(x.m_delta_values);
	return (*this);
}

gms::system::HRCTimer &
gms::system::HRCTimer
::operator=(HRCTimer &&x) {
	if (this == &x) return (*this);
			m_nruns = x.m_nruns;
			m_datum_size = x.m_datum_size;
			m_ptfunc = x.m_ptfunc;
			m_ptwrap = x.m_ptwrap;
			m_func_name.operator=(std::move(x.m_func_name));
			m_file_name.operator=(std::move(x.m_file_name));
			m_timing_values.operator=(std::move(x.m_timing_values));
			m_overhead_values.operator=(std::move(x.m_overhead_values));
			m_delta_values.operator=(std::move(x.m_delta_values));
	return (*this);
}

void
gms::system
::HRCTimer
::run_benchmark(_In_ const double eps) {
	
	using namespace gms::math::constants;
	using namespace gms::common;
	
	TIME_VAL STARTW, STOPW, STARTF, STOPF;
	PERIOD tdiff1,tdiff2;
	
	// function wrapper measurement loop.
	for (std::size_t i = 0Ui64; i != m_nruns; ++i) {
		for (std::size_t j = 0Ui64; j != m_datum_size; ++j) {
				
			STARTW = std::chrono::high_resolution_clock::now();
			m_ptwrap();
			STOPW  = std::chrono::high_resolution_clock::now();
			tdiff1 = STOPW - STARTW;
			m_overhead_values.operator[](Ix2D(i, m_datum_size, j)) = tdiff1.count();
		}
	}
	// function itself measurement loop
	for (std::size_t i = 0Ui64; i != m_nruns; ++i) {
		for (std::size_t j = 0Ui64; j != m_datum_size; ++j) {

			STARTF = std::chrono::high_resolution_clock::now();
			m_ptfunc();
			STOPF = std::chrono::high_resolution_clock::now();
			tdiff2 = STOPF - STARTF;
			m_timing_values.operator[](Ix2D(i, m_datum_size, j)) = tdiff2.count();
		}
	}
	// Delta mesurement loop.
	for (std::size_t i = 0Ui64; i != m_nruns; ++i) {
		for (std::size_t j = 0Ui64; j != m_datum_size; ++j) {
#if (USE_ACCURATE_IEEE754_2008_FP) == 1
			if (__binary64_quiet_greater_binary64(m_timing_values[i],
				                                  m_overhead_values[i])) {
				m_delta_values[i] = __binary64_sub_binary64_binary64(
											m_timing_values[i], m_overhead_values[i]);
			}
			else {
				m_delta_values[i] = 0.0;
			}
#else
			if ( definitely_greaterf64(m_timing_values[i],
									 m_overhead_values[i],eps)){
			 
				m_delta_values[i] = m_timing_values[i] = m_overhead_values[i];
			}
			else {
				m_delta_values[i] = 0.0;
			}
#endif
		}
	}
}

#if (USE_ACCURATE_IEEE754_2008_FP) == 1

bool
gms::system
::HRCTimer
::compute_stats(double &mean,
		double &adev,
		double &sdev,
		double &skew,
		double &kurt){
									 
		using namespace gms::math::constants;
	
		std::size_t vlen;
		__declspec(align(64)) struct {
			double len{}, sum{ 0.0 }, var{}, t{}, t2{}, t3{},
			t4{}, t5{}, isdev{};
		}r8_loc;
		vlen = m_delta_values.size();
		r8_loc.len = __binary64_from_uint64(vlen);
		for (std::size_t i = 0Ui64; i != vlen; ++i) {
			r8_loc.sum = 
				__binary64_add_binary64_binary64(r8_loc.sum, m_delta_values[i]);
			
		}	
		
		mean = __binary64_div_binary64_binary64(r8_loc.sum,r8_loc.len);
		//  Compute average deviation and variance
		for (std::size_t i = 0Ui64; i != vlen; ++i) {
			if (__binary64_quiet_not_equal_binary64(0.0, m_delta_values[i])) {
				r8_loc.t = __binary64_abs(
					__binary64_sub_binary64_binary64(m_delta_values[i], mean));
				adev = __binary64_add_binary64_binary64(adev, r8_loc.t);
				r8_loc.t2 = __binary64_mul_binary64_binary64(
					__binary64_sub_binary64_binary64(m_delta_values[i], mean),
					__binary64_sub_binary64_binary64(m_delta_values[i], mean));
				r8_loc.var = __binary64_add_binary64_binary64(r8_loc.var, r8_loc.t2);
			}
		}
		adev = __binary64_div_binary64_binary64(adev, r8_loc.len);
		if (__binary64_quiet_less_binary64(r8_loc.var, 0.0)) {
			std::printf(" HRCTimer::compute_stats: Invalid variance: %.15f\n", r8_loc.var);
			    mean = dinf, adev = dinf,
				sdev = dinf, skew = dinf,
				kurt = dinf;
			return (false);
		}
		r8_loc.var = __binary64_div_binary64_binary64(r8_loc.var, r8_loc.len);
		sdev = __binary64_sqrt_binary64(r8_loc.var);
		r8_loc.isdev = __binary64_div_binary64_binary64(1.0, sdev);
		for (std::size_t i = 0Ui64; i != vlen; ++i) {
			if (__binary64_quiet_not_equal_binary64(0.0, m_delta_values[i])) {
				r8_loc.t3 = __binary64_mul_binary64_binary64(
					__binary64_sub_binary64_binary64(m_delta_values[i], mean), r8_loc.isdev);
				r8_loc.t4 = __binary64_mul_binary64_binary64(r8_loc.t3,
					__binary64_mul_binary64_binary64(r8_loc.t3, r8_loc.t3));
				skew = __binary64_add_binary64_binary64(skew, r8_loc.t4);
				r8_loc.t5 = __binary64_mul_binary64_binary64(r8_loc.t4, r8_loc.t3);
				kurt = __binary64_add_binary64_binary64(kurt, r8_loc.t5);
			}
		}
		skew = __binary64_div_binary64_binary64(skew, r8_loc.len);
		kurt = __binary64_div_binary64_binary64(kurt,
			__binary64_sub_binary64_binary64(r8_loc.len, 3.0));
		if (__binary64_quiet_less_equal_binary64(kurt, 1.0)) {
			std::printf(" HRCTimer::compute_stats: Invalid kurtosis: %.15f\n", kurt);
				mean = dinf, adev = dinf,
				sdev = dinf, skew = dinf,
				kurt = dinf;
			return (false);
		}
		return (true);
}
#else

bool
gms::system::HRCTimer
::compute_stats(double &mean,
		double &adev,
		double &sdev,
		double &skew,
		double &kurt) {
		using namespace gms::math::constants;
	std::size_t vlen{};
	__declspec(align(64)) struct{
		double len{}, sum{0.0}, var{}, t{}, t2{}, t3{}, tmp{},
		t4{}, t5{}, isdev{}, fracp{}, ct2{};
	}r8_loc;
	vlen = m_delta_values.size();
	r8_loc.len = static_cast<double>(vlen);
	for (std::size_t i = 0Ui64; i != vlen; ++i) {
		sum += m_delta_values[i];
	}
	mean = sum / r8_loc.len;
	//  Compute average deviation and variance
	for (size_t i = 0Ui64; i != vlen; ++i) { // <-- no auto-vectorization (no countable loop)
		if (0.0 != m_delta_values.operator[](i)) {

			r8_loc.t = std::abs(m_delta_values[i] - mean); //potential catastrophic cancellation if(tmp - mean) both are very close to each other.
			adev += r8_loc.t; // <-- here potential overflow
			r8_loc.t2 = (m_delta_values[i] - mean) * (m_delta_values[i] - mean);
			r8_loc.ct2 = r8_loc.t2; // Checks if catastrophic cancellation has occurred
			r8_loc.fracp = r8_loc.t2 - static_cast<uint64_t>(r8_loc.ct2);
			if (r8_loc.fracp <= DEPS) {
				std::printf(" HRCTimer::compute_stats: Losing a significand digits: %.16f\n", r8_loc.fracp);
			}
			r8_loc.var += r8_loc.t2; // potential overflow
		}
	}
	adev /= r8_loc.len;
	if (r8_loc.var < 0.0) {
		std::printf(" HRCTimer::compute_stats: Invalid variance: %.15f\n", r8_loc.var);
		mean = dinf, adev = dinf,
		sdev = dinf, skew = dinf,
		kurt = dinf;
		return (false);
	}
	r8_loc.var /= r8_loc.len;
	sdev = std::sqrt(r8_loc.var);
	r8_loc.isdev = 1.0 / sdev;
	r8_loc.fracp = -1.0;
	for (size_t i = 0Ui64; i != vlen; ++i) {
		if (0.0 != m_delta_values.operator[](i)) {

			r8_loc.t3 = (m_delta_values[i] - mean) * r8_loc.isdev; // <-- catastrophic cancellation here
			r8_loc.ct2 = r8_loc.t3;
			r8_loc.fracp = r8_loc.t3 - static_cast<int64_t>(r8_loc.ct2);
			if (r8_loc.fracp <= DEPS) {
				std::printf(" HRCTimer::compute_stats: Losing a significant digits: %.16f\n", r8_loc.fracp);
			}

			r8_loc.t4 = r8_loc.t3*r8_loc.t3*r8_loc.t3; // <-- // Potential overflow?
			skew += r8_loc.t4; // <-- // Potential overflow?
			r8_loc.t5 = r8_loc.t4 * r8_loc.t3;
			kurt += r8_loc.t5; // <-- // Potential overflow?
		}
	}
	skew /= r8_loc.len;
	kurt /= r8_loc.len - 3.0;
	if (kurt < 1.0) {
		std::printf(" HRCTimer::compute_stats: Invalid kurtosis: %.15f\n", kurt);
		mean = dinf, adev = dinf,
		sdev = dinf, skew = dinf,
		kurt = dinf;
		return (false);
	}
	return (true);
}

#endif

void
gms::system
::HRCTimer
::print() const {
	std::cout << "	Dumping state of  " << typeid(*this).raw_name() << "\n";
	std::cout << "============================================================\n";
	std::cout << " Number of test runs:      " << m_nruns << "\n"
		<< " Datum size per run:             " << m_datum_size << "\n"
		<< " Tested function address:        " << std::hex << m_ptfunc << "\n"
		<< " Tested wrapper  address:        " << std::hex << m_ptwrap << "\n"
		<< " Tested function name:           " << m_func_name.data() << "\n"
		<< " Located in file name:           " << m_file_name.data() << "\n";
	std::cout << "======================= Printing data set =========================\n";
	for (std::size_t i = 0Ui64; i != m_nruns; ++i) {
		std::cout << " Func sample  :    Wrapper sample  :    Delta          \n";
		for (std::size_t j = 0Ui64; j != m_datum_size; ++j) {
			std::cout << std::setprecision(15) 
				<< std::setw(4)  <<   m_timing_values[Ix2D(i,m_datum_size,j)] 
				<< std::setw(8)  <<   m_overhead_values[Ix2D(i,m_datum_size,j)]
				<< std::setw(12) << m_delta_values[Ix2D(i, m_datum_size, j)] << "\n";
		}
		std::cout << "\n";
	}
}


#if !defined (GMS_CHRONO_TIMER_CHECK_ARGS_CONFORMANCE)
#define GMA_CHRONO_TIMER_CHECK_ARGS_CONFORMANCE(vres,set1,set2)              \
	  do{								   \
          if ((vres).size() != (set1).m_delta_values.size() ||               \
		      (set1).m_delta_values.size() != (set2).m_delta_values.size()){ \
			   return (false);					  \
		}								  \
	} while (0); 
#endif

bool
gms::system::HRCTimer
::delta_values_eq(std::vector<bool> &vres,
		  const HRCTimer &set1,
		  const HRCTimer &set2,
		  const double eps) {
		using namespace gms::common;
	GMS_CHRONO_TIMER_CHECK_ARGS_CONFORMANCE(vres,set1,set2)
	
	std::size_t vlen{ set1.m_delta_values.size() };
	for (std::size_t i = 0Ui64; i != vlen; ++i) {
		if (approximately_equalf64(set1.m_delta_values[i],
				    set2.m_delta_values[i], eps)) {
			vres.operator[](i) = true;
		}
		else {
			vres.operator[](i) = false;
		}
	}
	return (true);
}

bool
gms::system::HRCTimer
::delta_values_ineq(std::vector<bool> &vres,
		    const HRCTimer &set1,
		    const HRCTimer &set2,
		    const double eps) {
	using namespace gms::common;
	GMS_CHRONO_TIMER_CHECK_ARGS_CONFORMANCE(vres, set1, set2)
	std::size_t vlen{ set1.m_delta_values.size() };
	for (std::size_t i = 0Ui64; i != vlen; ++i) {
		if (!approximately_equalf64(set1.m_delta_values[i],
									set2.m_delta_values[i], eps)) {
			vres.operator[](i) = true;
		}
		else {
			vres.operator[](i) = false;
		}
	}
	return (true);
}

bool
gms::system::HRCTimer
::delta_values_lt(std::vector<bool> &vres,
		  const HRCTimer &set1,
		  const HRCTimer &set2,
		  const double eps) {
	using namespace gms::common;
	GMS_CHRONO_TIMER_CHECK_ARGS_CONFORMANCE(vres, set1, set2)
	std::size_t vlen{ set1.m_delta_values.size() };
	for (std::size_t i = 0Ui64; i != vlen; ++i) {
		if (definitely_lessf64(set1.m_delta_values[i],
							   set2.m_delta_values[i], eps)) {
			vres.operator[](i) = true;
		}
		else {
			vres.operator[](i) = false;
		}
	}
	return (true);
}

bool
gms::system::HRCTimer
::delta_values_gt(std::vector<bool> &vres,
		  const HRCTimer &set1,
		  const HRCTimer &set2,
		  const double eps) {
	using namespace gms::common;
	GMS_CHRONO_TIMER_CHECK_ARGS_CONFORMANCE(vres, set1, set2)
	std::size_t vlen{ set1.m_delta_values.size() };
	for (std::size_t i = 0Ui64; i != vlen; ++i) {
		if (definitely_greaterf64(set1.m_delta_values[i],
								  set2.m_delta_values[i], eps)) {
			vres.operator[](i) = true;
		}
		else {
			vres.operator[](i) = false;
		}
	}
	return (true);
}
