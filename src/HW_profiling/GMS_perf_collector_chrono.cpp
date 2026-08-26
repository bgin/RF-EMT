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
#include "GMS_perf_collector_chrono.h"
#include "GMS_malloc.h"
#include "GMS_constants.h"

#if (USE_MKL) == 1 && defined (GMS_COMPILED_BY_ICC)

#if !defined (PERF_COLLECTOR_CHRONO_DFTI_FAIL_CLEANUP)
#define PERF_COLLECTOR_CHRONO_DFTI_FAIL_CLEANUP(status,data1,data2,desc)	\
	   do{																	\
			if ((status) != 0) {										    \
			    _mm_free((data1));											\
			    _mm_free((data2));										    \
				DftiFreeDescriptor(&(desc));								\
				return (false);												\
		   }																\
	   } while (0);
#endif


#if !defined (PERF_COLLECTOR_CHRONO_DFTI_FAIL_CLEANUP2)
#define PERF_COLLECTOR_CHRONO_DFTI_FAIL_CLEANUP2(status,data1,data2,desc1,desc2)		\
	do {																		\
			if ((status) != 0) {											    \
			_mm_free((data1));												    \
			_mm_free((data2));													\
			DftiFreeDescriptor(&(desc1));									    \
			DftiFreeDescriptor(&(desc2));										\
			return (false);														\
		}																		\
																				\
	} while (0);
#endif

#endif

gms::system
::PerfCollectorChrono
::PerfCollectorChrono()
:
m_Iscleared{},
m_nsamples{},
m_funcname{},
m_filename{},
m_date{},
m_time{},
m_nvalid_values(),
m_start_values(),
m_stop_values(),
m_delta_values()
{}

gms::system
::PerfCollectorChrono
::PerfCollectorChrono(const PerfCollectorChrono &x)
:
m_Iscleared{ x.m_Iscleared },
m_nsamples{ x.m_nsamples },
m_funcname{ x.m_funcname },
m_filename{ x.m_filename },
m_date{ x.m_date },
m_time{ x.m_time },
m_nvalid_values{ x.m_nvalid_values },
m_start_values{ x.m_start_values },
m_stop_values{ x.m_stop_values },
m_delta_values{ x.m_delta_values }
{}

gms::system
::PerfCollectorChrono
::PerfCollectorChrono(PerfCollectorChrono &&x)
:
m_Iscleared{ x.m_Iscleared },
m_nsamples{ x.m_nsamples },
m_funcname{ std::move(x.m_funcname) },
m_filename{ std::move(x.m_filename) },
m_date{ std::move(x.m_date) },
m_time{ std::move(x.m_time) },
m_nvalid_values{ std::move(x.m_nvalid_values) },
m_start_values{ std::move(x.m_start_values) },
m_stop_values{ std::move(x.m_stop_values) },
m_delta_values{ std::move(x.m_delta_values) }
{}

gms::system::PerfCollectorChrono &
gms::system::PerfCollectorChrono
::operator=(const PerfCollectorChrono &x){
	if (this == &x) return (*this);
			m_Iscleared = x.m_Iscleared;
			m_nsamples = x.m_nsamples;
			m_funcname.operator=(x.m_funcname);
			m_filename.operator=(x.m_filename);
			m_date.operator=(x.m_date);
			m_time.operator=(x.m_time);
			m_nvalid_values.operator=(x.m_nvalid_values);
			m_start_values.operator=(x.m_start_values);
			m_stop_values.operator=(x.m_stop_values);
			m_delta_values.operator=(x.m_delta_values);
	 return (*this);
}

gms::system::PerfCollectorChrono &
gms::system::PerfCollectorChrono
::operator=(PerfCollectorChrono &&x) {
	if (this == &x) return (*this);
			m_Iscleared = x.m_Iscleared;
			m_nsamples = x.m_nsamples;
			m_funcname.operator=(std::move(x.m_funcname));
			m_filename.operator=(std::move(x.m_filename));
			m_date.operator=(std::move(x.m_date));
			m_time.operator=(std::move(x.m_time));
			m_nvalid_values.operator=(std::move(x.m_nvalid_values));
			m_start_values.operator=(std::move(x.m_start_values));
			m_stop_values.operator=(std::move(x.m_stop_values));
			m_delta_values.operator=(std::move(x.m_delta_values));
	return (*this);
}

void
gms::system
::PerfCollectorChrono::start() {
	using namespace std::chrono;
	m_start_values.push_back(high_resolution_clock::now()); // Exception unsafe.
}

void
gms::system
::PerfCollectorChrono::stop() {
	using namespace std::chrono;
	m_stop_values.push_back(high_resolution_clock::now());
}

void
gms::system
::PerfCollectorChrono::clear_all() {
	if (!m_Iscleared) {
		if (m_start_values.size() > lowest) {
			for (std::size_t i = 0Ui64; i != m_start_values.size(); ++i) {
				m_start_values[i] = TIME_VAL{};
			}
		}
		if (m_stop_values.size() > lowest) {
			for (std::size_t i = 0Ui64; i != m_stop_values.size(); ++i) {
				m_stop_values[i] = TIME_VAL{};
			}
		}
		if (m_delta_values.size() > lowest) {
			for (std::size_t i = 0Ui64; i != m_delta_values.size(); ++i) {
				m_delta_values[i] = 0.0;
			}
		}
	}
}

bool
gms::system
::PerfCollectorChrono
::compute_delta() {
	if (!check_data_size()) {
		return (false);
	}
	double result{0.0};
	PERIOD tdiff;
	for (std::size_t i = 0Ui64; i != m_start_values.size(); ++i) {
		tdiff = m_stop_values[i] - m_start_values[i];
		if (tdiff.count() > 0.0){
		    m_delta_values.push_back(tdiff.count());
			m_nvalid_values.push_back(true);
			++m_nsamples;
		}
		else {
			m_delta_values.push_back(0.0);
			m_nvalid_values.push_back(false);
		}
	}
	return (true);
}

#if (USE_ACCURATE_IEEE754_2008_FP) == 1
	 
bool
gms::system
::PerfCollectorChrono
::compute_stats(double &mean,
	        double &adev,
		double &sdev,
		double &skew,
		double &kurt) {
	using namespace gms::math::constants;
	if (m_Iscleared || m_delta_values.size() < lo_bound) {
			mean = dinf, adev = dinf,
			sdev = dinf, skew = dinf,
			kurt = dinf;
		return (false);
	}
	std::size_t vlen{m_delta_values.size()};
	__declspec(align(64)) struct {
		double len{}, var{}, t{}, t2{}, t3{},
		t4{}, t5{}, isdev{}, sum{0.0};
	}r8_loc;
	r8_loc.len = __binary64_from_uint64(vlen);
	// Compute mean value.
	for (std::size_t i = 0Ui64; i != vlen; ++i) {
		 r8_loc.sum = __binary64_add_binary64_binary64(r8_loc.sum,
						         m_delta_values[i]);
	}
	mean = r8_loc.sum / r8_loc.len;
	//  Compute average deviation and variance
	for (std::size_t i = 0Ui64; i != vlen; ++i) {
		if (__binary64_quiet_not_equal_binary64(0.0,
							m_delta_values[i])) {
			r8_loc.t = __binary64_abs(__binary64_sub_binary64_binary64(
											m_delta_values[i],mean));
			adev = __binary64_add_binary64_binary64(adev,r8_loc.t);
			r8_loc.t2 = __binary64_mul_binary64_binary64(
						__binary64_sub_binary64_binary64(m_delta_values[i],mean),
						__binary64_sub_binary64_binary64(m_delta_values[i], mean));
			r8_loc.var = __binary64_add_binary64_binary64(r8_loc.var,r8_loc.t2);
		}
	}
	adev = __binary64_div_binary64_binary64(adev,r8_loc.len);
	if (__binary64_quiet_less_binary64(r8_loc.var, 0.0)) {
		std::printf(" PerfCollectorChrono::compute_stats: Invalid variance: %.15f\n", r8_loc.var);
		    mean = dinf, adev = dinf,
			sdev = dinf, skew = dinf,
			kurt = dinf;
		return (false);
	}
	r8_loc.var = __binary64_div_binary64_binary64(r8_loc.var,r8_loc.len);
	sdev = __binary64_sqrt_binary64(r8_loc.var);
	r8_loc.isdev = __binary64_div_binary64_binary64(1.0,sdev);
	for (std::size_t i = 0Ui64; i != vlen; ++i) {
		if (__binary64_quiet_not_equal_binary64(0.0, 
									m_delta_values[i])) {
			r8_loc.t3 = __binary64_mul_binary64_binary64(
				__binary64_sub_binary64_binary64(m_delta_values[i],mean),r8_loc.isdev);
			r8_loc.t4 = __binary64_mul_binary64_binary64(r8_loc.t3,
				__binary64_mul_binary64_binary64(r8_loc.t3, r8_loc.t3));
			skew = __binary64_add_binary64_binary64(skew,r8_loc.t4);
			r8_loc.t5 = __binary64_mul_binary64_binary64(r8_loc.t4,r8_loc.t3);
			kurt = __binary64_add_binary64_binary64(kurt,r8_loc.t5);
		}
	}
	skew = __binary64_div_binary64_binary64(skew, r8_loc.len);
	kurt = __binary64_div_binary64_binary64(kurt,
					__binary64_sub_binary64_binary64(r8_loc.len, 3.0));
	if (__binary64_quiet_less_equal_binary64(kurt, 1.0)) {
		std::printf(" PerfCollectorChrono::compute_stats: Invalid kurtosis: %.15f\n", kurt);
			mean = dinf, adev = dinf,
			sdev = dinf, skew = dinf,
			kurt = dinf;
		return (false);
	}
	return (true);
}

#else

bool
gms::system
::PerfCollectorChrono::
compute_stats(double &mean,
	      double &adev,
	      double &sdev,
	      double &skew,
	      double &kurt) {
	using namespace gms::math::constants;
	if (m_Iscleared || m_delta_values.size() < lo_bound) {
		    mean = dinf, adev = dinf,
			sdev = dinf, skew = dinf,
			kurt = dinf;
		return (false);
	}
	uint64_t vlen{};
	__declspec(align(64)) struct{
		double len{}, sum{}, var{}, t{}, t2{}, t3{}, tmp{},
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
				std::printf(" PerfCollectorChrono::compute_stats: Losing a significand digits: %.16f\n", r8_loc.fracp);
			}
			r8_loc.var += r8_loc.t2; // potential overflow
		}
	}
	adev /= r8_loc.len;
	if (r8_loc.var < 0.0) {
		std::printf(" PerfCollectorChrono::compute_stats: Invalid variance: %.15f\n", r8_loc.var);
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
				std::printf(" PerfCollectorChrono::compute_stats: Losing a significant digits: %.16f\n", r8_loc.fracp);
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
		    std::printf(" PerfCollectorChrono::compute_stats: Invalid kurtosis: %.15f\n", kurt);
			mean = dinf, adev = dinf,
			sdev = dinf, skew = dinf,
			kurt = dinf;
		return (false);
	}
	return (true);
}

#endif

bool
gms::system
::PerfCollectorChrono
::correlation_set(PerfCollectorChrono &set1,
		  PerfCollectorChrono &set2,
		  double * __restrict correlated_set,
		  int32_t * __restrict ip,
		  double * __restrict w,
		  const int32_t n) {
	if (set1.m_delta_values.size() != set2.m_delta_values.size() ||
		static_cast<size_t>(n) != set1.m_delta_values.size()) {
		return (false);
	}
#if (USE_MKL) == 1
	DFTI_DESCRIPTOR_HANDLE ddh1, ddh2;
	MKL_LONG status{};
#endif
	double t{}, invals{};
	int32_t snvals{};
	auto set_size = set1.m_delta_values.size();


	typedef double * __restrict __attribute__((aligned(64))) aligned_r8ptr;

	aligned_r8ptr dataset1 = (double*)gms::common::gms_mm_malloc(set_size,align64B);
	aligned_r8ptr dataset2 = (double*)gms::common::gms_mm_malloc(set_size,align64B);
#if (USE_ACCURATE_IEEE754_2008_FP) == 1
	auto bOk1 = set1.cvrt_to_double(&dataset1[0],set_size);
	auto bOk2 = set2.cvrt_to_double(&dataset2[0],set_size);
#else
        auto bOk1 = set1.cvrt_to_double(&dataset1[0],set_size);
	auto bOk2 = set2.cvrt_to_double(&dataset2[0],set_size);
#endif
	if (!bOk1 || !bOk2) {
	        gms_mm_free(dataset1);
		gms_mm_free(dataset2);
		return (false);
	}
	snvals = n << 1;
	invals = 1.0 / static_cast<double>(snvals);
	// MKL_DFTI
#if (USE_MKL) == 1
	status = DftiCreateDescriptor(&ddh1,DFTI_DOUBLE,DFTI_REAL,1,set_size);
	PERF_COLLECTOR_CHRONO_DFTI_FAIL_CLEANUP(status,dataset1,dataset2,ddh1)
	status = DftiCommitDescriptor(ddh1);
	PERF_COLLECTOR_CHRONO_DFTI_FAIL_CLEANUP(status,dataset1,dataset2,ddh1)
	status = DftiComputeForward(ddh1,dataset1);
	PERF_COLLECTOR_CHRONO_DFTI_FAIL_CLEANUP(status,dataset1,dataset2,ddh1)
	status = DftiCreateDescriptor(&ddh2,DFTI_DOUBLE,DFTI_REAL,1,set_size);
	PERF_COLLECTOR_CHRONO_DFTI_FAIL_CLEANUP2(status, dataset1, dataset2, ddh1, ddh2)
	status = DftiCommitDescriptor(ddh2);
	PERF_COLLECTOR_CHRONO_DFTI_FAIL_CLEANUP2(status, dataset1, dataset2, ddh1, ddh2)
	status = DftiComputeForward(ddh2, dataset2);
	PERF_COLLECTOR_CHRONO_DFTI_FAIL_CLEANUP2(status, dataset1, dataset2, ddh1, ddh2)
	for (int32_t i = 0; i != n - 1; i += 2) {
		t = dataset1[i];
		dataset1[i] = (dataset1[i] * dataset2[i] + dataset1[i + 1] * dataset2[i + 1]) * invals;
		dataset1[i + 1] = (dataset1[i + 1] * dataset2[i] - t * dataset2[i + 1]) * invals;
	}
	dataset1[0] = dataset1[0] * dataset2[0] * invals;
	dataset1[1] = dataset1[1] * dataset2[1] * invals;
	status = DftiComputeBackward(ddh1, dataset1);
	PERF_COLLECTOR_CHRONO_DFTI_FAIL_CLEANUP2(status, dataset1, dataset2, ddh1, ddh2)
	memcpy(&correlated_set[0], &dataset1[0],set_size);
	gms_mm_free(dataset1); gms_mm_free(dataset2);
	DftiFreeDescriptor(&ddh1); DftiFreeDescriptor(&ddh2);
	return (true);
#else
	ip[0] = 0;
	rdft(n, 1, &dataset1[0], &ip[0], &w[0]);
	rdft(n, 1, &dataset2[0], &ip[0], &w[0]);
	for (int32_t i = 1; i != n - 1; i += 2) {
		t = dataset1[i];
		dataset1[i] = (dataset1[i] * dataset2[i] + dataset1[i + 1] * dataset2[i + 1]) * invals;
		dataset1[i + 1] = (dataset1[i + 1] * dataset2[i] - t * dataset2[i + 1]) * invals;
	}
	dataset1[0] = dataset1[0] * dataset2[0] * invals;
	dataset1[1] = dataset1[1] * dataset2[1] * invals;
	rdft(n, -1, &dataset1[0], &ip[0], &w[0]);
	memcpy(&correlated_set[0], &dataset1[0], set_size);

	if (dataset1) gms_mm_free(dataset1);
	if (dataset2) gms_mm_free(dataset2);
	return (true);
#endif

}

void
gms::system
::PerfCollectorChrono::print() const {
	std::cout << " Dumping state of: " << typeid(*this).raw_name() << "\n";
	std::cout << "======================================================================================\n";
	std::cout << " Reset state of collector: " << std::boolalpha << m_Iscleared << "\n"
		<< "       Number of valid samples:  " << m_nsamples << "\n"
		<< "       Collected at function:    " << m_funcname.data() << "\n"
		<< "       Collected in file:        " << m_filename.data() << "\n"
		<< "       Collected at date:        " << m_date.data() << "\n"
		<< "       Collected at time:        " << m_time.data() << "\n"
		<< "       ======================================================================================\n";
	std::cout << "  Sample_state    :   delta_value    \n";
	for (std::size_t i = 0Ui64; i != m_start_values.size(); ++i) {
		std::cout << std::setw(4) << m_nvalid_values[i] 
			      << std::setw(8) << std::setprecision(16) << m_delta_values[i] << "\n";
	}
}

bool
gms::system
::PerfCollectorChrono
::check_data_size() const {
	if (!m_Iscleared) {
		if (m_start_values.size() == m_stop_values.size()){
			 
			return (true);
		}
		else { return (false); }
	}
	return (false);
}


#if (USE_ACCURATE_IEEE754_2008_FP) == 1
bool
gms::system
::PerfCollectorChrono
::cvrt_to_double(double * __restrict data,
		 const std::size_t data_len) {
	if (m_delta_values.size() != data_len) {
		return (false);
	}
	// Begin non vectorized conversion (numerically stable conversion)
	for (std::size_t i = 0Ui64; i != data_len; ++i) {
		data[i] = __binary64_from_uint64(m_delta_values[i]);
	}
	return (true);
}
#else
bool
gms::system
::PerfCollectorChrono
::cvrt_to_double(double * __restrict data,
		 const std::size_t data_len) {
       	if (m_delta_values.size() != data_len) {
		return (false);
	}
	for(std::size_t i = 0Ui64; i != data_len; ++i) {
		data[i] = static_cast<double>(m_delta_values[i]);
	}
	return (true);
}
#endif
