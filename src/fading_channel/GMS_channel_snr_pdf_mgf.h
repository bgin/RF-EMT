/*
 * Copyright (C) Bernard Gingold, 2020-2026 
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

#ifndef __GMS_CHANNEL_SNR_PDF_MGF_H__
#define __GMS_CHANNEL_SNR_PDF_MGF_H__ 270720261119

#include <cstdint>
#include <cmath> // bessel J0
#include <complex> // for various analytical MGF
#include "GMS_config.h"
#if (USE_OPENMP) == 1
#include <omp.h>
#endif
#include "GMS_cephes_sin_cos.h"


namespace file_info 
{

     static const unsigned int GMS_CHANNEL_SNR_PDF_MGF_MAJOR = 1;
     static const unsigned int GMS_CHANNEL_SNR_PDF_MGF_MINOR = 1;
     static const unsigned int GMS_CHANNEL_SNR_PDF_MGF_MICRO = 0;
     static const unsigned int GMS_CHANNEL_SNR_PDF_MGF_FULLVER =
       1000U*GMS_CHANNEL_SNR_PDF_MGF_MAJOR+100U*GMS_CHANNEL_SNR_PDF_MGF_MINOR+
       10U*GMS_CHANNEL_SNR_PDF_MGF_MICRO;
     static const char GMS_CHANNEL_SNR_PDF_MGF_CREATION_DATE[] = "27-07-2026 11:20AM +00200 (MON 27 JUL 2026 GMT+2)";
     static const char GMS_CHANNEL_SNR_PDF_MGF_BUILD_DATE[]    = __DATE__; 
     static const char GMS_CHANNEL_SNR_PDF_MGF_BUILD_TIME[]    = __TIME__;
     static const char GMS_CHANNEL_SNR_PDF_MGF_SYNOPSIS[]      = "Analytical PDF and MGF of the SNR for the common models of fading channels.\
	                                                          Based on the M.K Simon, M.S. Alouini: Digital Communication over Fading Channels 1st ed\
															  ISBN-13 978-0471317791";

}

namespace gms 
{

namespace fading_channel
{

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif 
__ATTR_ALWAYS_INLINE__
static inline 
float
rayleigh_fading_pdf(const float snr_instant,
                    const float snr_avg)
{
    const float inv_snr_avg{1.0f/snr_avg};
    const float snr_i_to_snr_a{snr_instant/snr_avg};
#if (USE_OPENMP) == 1
    const float exp_val{std::exp(-snr_i_to_snr_a)};
#else 
    const float exp_val{gms::math::ceph_expf(snr_i_to_snr_a)};
#endif 
    return (inv_snr_avg*exp_val);
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif 
__ATTR_ALWAYS_INLINE__
static inline 
std::complex<float>
rayleigh_fading_mgf(const std::complex<float> s,
                    const float snr_avg)
{
    std::complex<float> s_m_snr_avg{s*snr_avg};
    std::complex<float> cterm{1.0f-s_m_snr_avg};
    return (1.0f/cterm);
}

/*
Precondition: 0<=q<=1
*/
#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif 
__ATTR_ALWAYS_INLINE__
static inline 
float
nakagami_q_fading_pdf(const float q,
                      const float snr_inst,
                      const float snr_avg)
{
    const float qpow2{q*q};
    const float qpow4{qpow2*qpow2};
    const float denom{1.0f/(4.0f*qpow2*snr_avg)};
    const float num{1.0f+qpow2};
    const float denom2{1.0f/(2.0f*q*snr_avg)};
    const float exp_arg{(num*num)*snr_inst*denom};
    const float factor{num*denom2};
    const float bes_j0_arg{(1.0f-qpow4)*snr_inst*denom};
    const float exp_val{exp_arg};
    const float bes_j0_val{std::cyl_bessel_j(0,bes_j0_arg)};
    return (factor*exp_val*bes_j0_val);
}

#if (USE_OPENMP) == 1
#pragma omp declare simd simdlen(16)
#endif 
__ATTR_ALWAYS_INLINE__
static inline 
std::complex<float>
nakagami_q_fading_mgf(const std::complex<float> s,
                      const float q,
                      const float snr_avg)
{
    const std::complex<float> s_m_snr_avg{s*snr_avg};
    const std::complex<float> cfactor{1.0f-2.0f*s_m_snr_avg};
    const float qpow2{q*q};
    const float denom{(1.0f+qpow2)*(1.0f+qpow2)};
    const std::complex<float> twosq{2.0f*s_m_snr_avg};
    const float inv_den{1.0f*denom};
    const std::complex<float> cmplx_num{(twosq*twosq)*qpow2};
    const std::complex<float> cmplx_ratio{cmplx_num*inv_den};
    const std::complex<float> cmplx_res{cfactor+cmplx_ratio};
    return (1.0f/(std::sqrt(cmplx_res)));
}


} // fading_channel

} // gms

#endif /*__GMS_CHANNEL_SNR_PDF_MGF_H__*/