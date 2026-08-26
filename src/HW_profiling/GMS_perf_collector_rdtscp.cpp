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
#include "GMS_perf_collector_rdtscp.h"


#include "GMS_malloc.h"
#include "GMS_constants.h"


#if (USE_MKL) == 1 && defined (GMS_COMPILED_BY_ICC)

       #if !defined (PERF_COLLECTOR_RDTSCP_DFTI_FAIL_CLEANUP)
       #define PERF_COLLECTOR_RDTSCP_DFTI_FAIL_CLEANUP(status,data1,data2,desc)        \
	        do {														               \
		        if ((status) != 0) {									               \
		            _mm_free((data1)); _mm_free((data2));							   \
					 DftiFreeDescriptor(&(desc));										   \
		             return (false);											       \
	             }																       \
	      } while (0);
       #endif

		#if !defined (PERF_COLLECTOR_RDTSCP_DFTI_FAIL_CLEANUP2)
	    #define PERF_COLLECTOR_RDTSCP_DFTI_FAIL_CLEANUP2(status,data1,data2,desc1,desc2) \
			do {																		  \
				if ((status) != 0) {													  \
				      _mm_free((data1)); _mm_free((data2));									  \
					  DftiFreeDescriptor(&(desc1));										  \
					  DftiFreeDescriptor(&(desc2));										  \
					  return (false);													  \
				}																		  \
			} while (0);
        #endif
#endif

gms::system
::PerfCollectorRDTSCP
::PerfCollectorRDTSCP()
:
m_Iscleared{},
m_nsamples{},
m_func_name{},
m_file_name{},
m_date{},
m_time{},
m_nvalid_values{},
m_start_values{},
m_stop_values{},
m_delta_values{},
m_tscaux_start{},
m_tscaux_stop{}
{}

gms::system
::PerfCollectorRDTSCP::
PerfCollectorRDTSCP(const PerfCollectorRDTSCP &x)
:
m_Iscleared{ x.m_Iscleared },
m_nsamples{  x.m_nsamples },
m_func_name{ x.m_func_name },
m_file_name{ x.m_file_name },
m_date{      x.m_date },
m_time{      x.m_time },
m_nvalid_values{ x.m_nvalid_values },
m_start_values{ x.m_start_values },
m_stop_values{ x.m_stop_values },
m_delta_values{x.m_delta_values},
m_tscaux_start{ x.m_tscaux_start },
m_tscaux_stop{ x.m_tscaux_stop }
{}

gms::system::PerfCollectorRDTSCP
::PerfCollectorRDTSCP(PerfCollectorRDTSCP &&x)
:
m_Iscleared{ std::move(x.m_Iscleared) },
m_nsamples{ std::move(x.m_nsamples) },
m_func_name{ std::move(x.m_func_name) },
m_file_name{ std::move(x.m_file_name) },
m_date{ std::move(x.m_date) },
m_time{ std::move(x.m_time) },
m_nvalid_values{ std::move(x.m_nvalid_values) },
m_start_values{ std::move(x.m_start_values) },
m_stop_values{ std::move(x.m_stop_values) },
m_delta_values{ std::move(x.m_delta_values) },
m_tscaux_start{ std::move(x.m_tscaux_start) },
m_tscaux_stop{ std::move(x.m_tscaux_stop) }
{}

gms::system::PerfCollectorRDTSCP &
gms::system::PerfCollectorRDTSCP
::operator=(const PerfCollectorRDTSCP &x) {
	if (this == &x) return (*this);
		m_Iscleared = x.m_Iscleared;
		m_nsamples  = x.m_nsamples;
		m_func_name.operator=(x.m_func_name);
		m_file_name.operator=(x.m_file_name);
		m_date.operator=(x.m_date);
		m_time.operator=(x.m_time);
		m_nvalid_values.operator=(x.m_nvalid_values);
		m_start_values.operator=(x.m_start_values);
		m_stop_values.operator=(x.m_stop_values);
		m_delta_values.operator=(x.m_delta_values);
		m_tscaux_start.operator=(x.m_tscaux_start);
		m_tscaux_stop.operator=(x.m_tscaux_stop);
	 return (*this);
}

gms::system::PerfCollectorRDTSCP &
gms::system::PerfCollectorRDTSCP
::operator=(PerfCollectorRDTSCP &&x){
	if (this == &x) return (*this);
	    m_Iscleared = std::move(x.m_Iscleared);
	    m_nsamples  = std::move(x.m_nsamples);
	    m_func_name.operator=(std::move(x.m_func_name));
	    m_file_name.operator=(std::move(x.m_file_name));
	    m_date.operator=(std::move(x.m_date));
	    m_time.operator=(std::move(x.m_time));
	    m_nvalid_values.operator=(std::move(x.m_nvalid_values));
	    m_start_values.operator=(std::move(x.m_start_values));
	    m_stop_values.operator=(std::move(x.m_stop_values));
		m_delta_values.operator=(std::move(x.m_delta_values));
	    m_tscaux_start.operator=(std::move(x.m_tscaux_start));
	    m_tscaux_stop.operator=(std::move(x.m_tscaux_stop));
	return (*this);
}

void
gms::system
::PerfCollectorRDTSCP::start() {
	int32_t dummy1[4] = {}, dummy2{};
	uint32_t tscaux{9999};
	__cpuid(&dummy1[0],dummy2);
	m_start_values.push_back(__rdtscp(&tscaux));
	m_tscaux_start.push_back(tscaux);
}

void
gms::system
::PerfCollectorRDTSCP::stop() {
	int32_t dummy1[4] = {}, dummy2{};
	uint32_t tscaux{9999};
	m_stop_values.push_back(__rdtscp(&tscaux));
	__cpuid(&dummy1[0],dummy2);
	m_tscaux_stop.push_back(tscaux);
}

void
gms::system
::PerfCollectorRDTSCP::clear_all() {
	if (!m_Iscleared) {  // Guard against potential misaligned size (by thrown exception i.e. std::bad_alloc)
		if (m_start_values.size() > lowest) {
			for (std::size_t i = 0Ui64; i != m_start_values.size(); ++i) {
				m_start_values.operator[](i) = zero;
			}
		}
		if (m_stop_values.size() > lowest) {
			for (std::size_t i = 0Ui64; i != m_stop_values.size(); ++i) {
				m_stop_values.operator[](i) = zero;
			}
		}
		if (m_nvalid_values.size() > lowest) {
			for (std::size_t i = 0Ui64; i != m_nvalid_values.size(); ++i) {
				m_nvalid_values.operator[](i) = false;
			}
		}
		if (m_delta_values.size() > lowest) {
			for (std::size_t i = 0Ui64; i != m_delta_values.size(); ++i) {
				m_delta_values.operator[](i) = zero;
			}
		}
		if (m_tscaux_start.size() >= lowest) {
			for (size_t i = 0ULL; i != m_tscaux_start.size(); ++i) {
				m_tscaux_start.operator[](i) = 9999;
			}
		}
		if (m_tscaux_stop.size() >= lowest) {
			for (size_t i = 0ULL; i != m_tscaux_stop.size(); ++i) {
				m_tscaux_stop.operator[](i) = 9999;
			}
		}
		m_Iscleared = true;
	}
}

bool
gms::system
::PerfCollectorRDTSCP
::compute_delta() {
	if (!check_data_size()) {
		return (false);
	}
	uint64_t res{};
	for (std::size_t i = 0Ui64; i != m_start_values.size(); ++i) {
		if ((res = m_stop_values[i] - m_start_values[i]) > zero) {
			m_delta_values.push_back(res);
			m_nvalid_values.push_back(true);
			++m_nsamples;
	    }
	     else {
		   m_delta_values.push_back(zero);
		   m_nvalid_values.push_back(false);
	    }
    }
	return (true);
}

#if (USE_ACCURATE_IEEE754_2008_FP) == 1
	 
bool
gms::system
::PerfCollectorRDTSCP
::compute_stats(double &mean ,
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
	uint64_t s{}, prev{}, vlen{};
	__declspec(align(64)) struct {
		double len{}, sum{}, var{}, t{}, t2{}, t3{}, tmp{},
		t4{}, t5{}, isdev{};
	}r8_loc;
	vlen = m_delta_values.size();
	r8_loc.len = __binary64_from_uint64(vlen);
	// Compute mean guard against an overflow
	for (std::size_t i = 0Ui64; i != vlen; ++i) {
		if (s < prev) return (false);
		prev = s;
		s += m_delta_values.operator[](i);
	}
	mean = __binary64_from_uint64(s) / r8_loc.len;
	//  Compute average deviation and variance
	for (std::size_t i = 0Ui64; i != vlen; ++i) {
		if (zero != m_delta_values[i]) { 
			r8_loc.tmp = __binary64_from_uint64(m_delta_values[i]);
			r8_loc.t = __binary64_abs(__binary64_sub_binary64_binary64(r8_loc.tmp,mean));
			adev = __binary64_add_binary64_binary64(adev,r8_loc.t);
			r8_loc.t2 = __binary64_mul_binary64_binary64(
				             __binary64_sub_binary64_binary64(r8_loc.tmp,mean),
				             __binary64_sub_binary64_binary64(r8_loc.tmp,mean));
			r8_loc.var = __binary64_add_binary64_binary64(r8_loc.var,r8_loc.t2);
		}
	}
	adev = __binary64_div_binary64_binary64(adev,r8_loc.len);
	if(__binary64_quiet_less_binary64(r8_loc.var,0.0)) {
		std::printf(" PerfCollectorRDTSCP::compute_stats: Invalid variance: %.15f\n", r8_loc.var);
		    mean = dinf, adev = dinf,
			sdev = dinf, skew = dinf,
			kurt = dinf;
		return (false);
	}
	r8_loc.var = __binary64_div_binary64_binary64(r8_loc.var,r8_loc.len);
	sdev = __binary64_sqrt_binary64(r8_loc.var);
	r8_loc.isdev = __binary64_div_binary64_binary64(1.0,sdev);
	r8_loc.tmp = 0.0;
	for (std::size_t i = 0Ui64; i != vlen; ++i) {
		if (zero != m_delta_values[i]) { // Yes I know it should be computed once before.
			r8_loc.tmp = __binary64_from_uint64(m_delta_values[i]);
			r8_loc.t3  = __binary64_mul_binary64_binary64(
			 	             __binary64_sub_binary64_binary64(r8_loc.tmp,mean),r8_loc.isdev);
			r8_loc.t4  = __binary64_mul_binary64_binary64(r8_loc.t3,
						     __binary64_mul_binary64_binary64(r8_loc.t3,r8_loc.t3));
			skew = __binary64_add_binary64_binary64(skew,r8_loc.t4);
			r8_loc.t5 = __binary64_mul_binary64_binary64(r8_loc.t4,r8_loc.t3);
			kurt = __binary64_add_binary64_binary64(kurt,r8_loc.t5);
		}			
	}	
	skew = __binary64_div_binary64_binary64(skew,r8_loc.len);
	kurt = __binary64_div_binary64_binary64(kurt, 
		         __binary64_sub_binary64_binary64(r8_loc.len,3.0));
	if (__binary64_quiet_less_equal_binary64(kurt, 1.0)) {
		std::printf(" PerfCollectorRDTSCP::compute_stats: Invalid kurtosis: %.15f\n", kurt);
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
::PerfCollectorRDTSCP::
compute_stats(double &mean ,
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
	uint64_t s{}, prev{};
	__declspec(align(64)) struct{
		double len{}, sum{}, var{}, t{}, t2{}, t3{}, tmp{},
		       t4{}, t5{}, isdev{}, fracp{}, ct2{};
	}r8_loc;
	r8_loc.len = static_cast<double>(m_delta_values.size());
	// Compute mean (guard against an overflow)
	for (std::size_t i = 0Ui64; i != m_delta_values.size(); ++i) {
		if (s < prev) return (false);
		prev = s;
		s += m_delta_values.operator[](i);
	}
	mean = static_cast<double>(s) / r8_loc.len;
	//  Compute average deviation and variance
	for (std::size_t i = 0Ui64; i != m_delta_values.size(); ++i) { // <-- no auto-vectorization (no countable loop)
		if (zero != m_delta_values.operator[](i)) {
			r8_loc.tmp = static_cast<double>(m_delta_values[i]);
			r8_loc.t   = std::abs(r8_loc.tmp - mean); //potential catastrophic cancellation if(tmp - mean) both are very close to each other.
			adev += r8_loc.t; // <-- here potential overflow
			r8_loc.t2 = (r8_loc.tmp - mean) * (r8_loc.tmp - mean);
			r8_loc.ct2 = r8_loc.t2; // Checks if catastrophic cancellation has occurred
			r8_loc.fracp = r8_loc.t2 - static_cast<uint64_t>(r8_loc.ct2);
			if (r8_loc.fracp <= DEPS) {
				std::printf(" PerfCollectorRDTSCP::compute_stats: Losing a significand digits: %.16f\n", r8_loc.fracp);
			}
			r8_loc.var += r8_loc.t2; // potential overflow
		}
	}
	adev /= r8_loc.len;
	if (r8_loc.var < 0.0) {
		std::printf(" PerfCollectorRDTSCP::compute_stats: Invalid variance: %.15f\n", r8_loc.var);
		return (false);
	}
	r8_loc.var /= r8_loc.len;
	sdev = std::sqrt(r8_loc.var);
	r8_loc.isdev = 1.0 / sdev;
	r8_loc.fracp = -1.0;
	r8_loc.tmp = 0.0;
	for (std::size_t i = 0Ui64; i != m_delta_values.size(); ++i) {
		if (m_delta_values.operator[](i) != zero) {
			r8_loc.tmp = static_cast<double>(m_delta_values[i]);
			r8_loc.t3 = (r8_loc.tmp - mean) * r8_loc.isdev; // <-- catastrophic cancellation here
			r8_loc.ct2 = r8_loc.t3;
			r8_loc.fracp = r8_loc.t3 - static_cast<int64_t>(r8_loc.ct2);
			if (r8_loc.fracp <= DEPS) {
				std::printf(" PerfCollectorRDTSCP::compute_stats: Losing a significant digits: %.16f\n", r8_loc.fracp);
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
		std::printf(" PerfCollectorRDTSCP::compute_stats: Invalid kurtosis: %.15f\n", kurt);
		return (false);
	}
	return (true);
}



#endif

bool
gms::system
::PerfCollectorRDTSCP
::correlation_set(PerfCollectorRDTSCP &set1,
		  PerfCollectorRDTSCP &set2,
		  double * __restrict correlated_set,
		  int32_t * __restrict ip,
		  double * __restrict w,
		  const int32_t n) {
	if (set1.m_delta_values.size() != set2.m_delta_values.size() ||
		static_cast<std::size_t>(n) != set1.m_delta_values.size()) {
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
	aligned_r8ptr dataset1 = gms::common::gms_mm_malloc(set_size,align64B);
	aligned_r8ptr dataset2 = gms::common::gms_mm_malloc(set_size,align64B);
#if (USE_ACCURATE_IEEE754_2008_FP) == 1
	auto bOk1 = set1.cvrt_to_double(&dataset1[0],set_size);
	auto bOk2 = set2.cvrt_to_double(&dataset2[0],set_size);
	
#else
	auto bOk1 = set1.cvrt_to_double(&dataset1[0],set_size);
	auto bOk2 = set2.cvrt_to_double(&dataset2[0],set_size);
#endif
	if (!bOk1 || !bOk2) { // If you got here somehow datasets length is corrupted
	        gms_mm_free(dataset1);
		gms_mm_free(dataset2);
	 	return (false);
	}
	snvals = n << 1;
	invals = 1.0 / static_cast<double>(snvals);
	// MKL DFTI
#if (USE_MKL) == 1
	status = DftiCreateDescriptor(&ddh1,DFTI_DOUBLE,DFTI_REAL,1,set_size);
	PERF_COLLECTOR_RDTSCP_DFTI_FAIL_CLEANUP(status,dataset1,dataset2,ddh1)
	status = DftiCommitDescriptor(ddh1);
	PERF_COLLECTOR_RDTSCP_DFTI_FAIL_CLEANUP(status,dataset1,dataset2,ddh1)
	status = DftiComputeForward(ddh1,dataset1);
	PERF_COLLECTOR_RDTSCP_DFTI_FAIL_CLEANUP(status,dataset1,dataset2,ddh1)
	status = DftiCreateDescriptor(&ddh2,DFTI_DOUBLE,DFTI_REAL,1,set_size);
	PERF_COLLECTOR_RDTSCP_DFTI_FAIL_CLEANUP2(status,dataset1,dataset2,ddh1,ddh2)
	status = DftiCommitDescriptor(ddh2);
	PERF_COLLECTOR_RDTSCP_DFTI_FAIL_CLEANUP2(status,dataset1,dataset2,ddh1,ddh2)
	status = DftiComputeForward(ddh2,dataset2);
	PERF_COLLECTOR_RDTSCP_DFTI_FAIL_CLEANUP2(status,dataset1,dataset2,ddh1,ddh2)
	for (int32_t i = 0; i != n - 1; i += 2) {
		t = dataset1[i];
		dataset1[i] = (dataset1[i] * dataset2[i] + dataset1[i + 1] * dataset2[i + 1]) * invals;
		dataset1[i + 1] = (dataset1[i + 1] * dataset2[i] - t * dataset2[i + 1]) * invals;
	}
	dataset1[0] = dataset1[0] * dataset2[0] * invals;
	dataset1[1] = dataset1[1] * dataset2[1] * invals;
	status = DftiComputeBackward(ddh1,dataset1);
	PERF_COLLECTOR_RDTSCP_DFTI_FAIL_CLEANUP2(status,dataset1,dataset2,ddh1,ddh2)
	memcpy(&correlated_set[0], &dataset1[0], set_size);
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
::PerfCollectorRDTSCP
::print() const {
	std::cout << " Dumping state of: " << typeid(*this).raw_name() << "\n";
	std::cout << "======================================================================================\n";
	std::cout << " Reset state of collector: " << std::boolalpha << m_Iscleared << "\n"
		<< "       Number of valid samples:  " << m_nsamples << "\n"
		<< "       Collected at function:    " << m_func_name.data() << "\n"
		<< "       Collected in file:        " << m_file_name.data() << "\n"
		<< "       Collected at date:        " << m_date.data() << "\n"
		<< "       Collected at time:        " << m_time.data() << "\n"
		<< "       ======================================================================================\n";
	std::cout << " Invalid :   start  :    stop   :    delta    :   TSC_AUX(start)  :   TSC_AUX(stop) \n"
		      << " ========================================================================================\n";
	for (size_t i = 0Ui64; i != m_delta_values.size(); ++i) {
	    	std::cout << std::setw(4)  <<  m_nvalid_values[i] 
				      << std::setw(8)  <<  m_start_values[i]
		              << std::setw(12) <<  m_stop_values[i] 
			          << std::setw(16) <<  m_delta_values[i] 
					  << std::setw(18) <<  m_tscaux_start[i]
			          << std::setw(22) <<  m_tscaux_stop[i] << "\n";
	}
	std::cout << "=========================================================================================\n"
	          << "==============================  End of dump =============================================\n";
}


// Helpers
bool
gms::system::PerfCollectorRDTSCP
::check_data_size() {
	if (!m_Iscleared) {
		if (m_start_values.size() == m_stop_values.size()){
			return (true);
		}
		else {
			return (false);
		}
	}
	else return (false);
}

#if (USE_ACCURATE_IEEE754_2008_FP) == 1
bool 
gms::system::PerfCollectorRDTSCP
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
gms::system::PerfCollectorRDTSCP
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
