
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
#include "GMS_bayesian_estimate.h"
#include "GMS_timsac_iface.h"


gms::math
::BayesEstimate
::BayesEstimate(const int32_t k,
		const int32_t n,
		const int32_t il,
		const int32_t mj2,
		const int32_t ival,
		const double dval)
:
m_k{ k },
m_n{ n },
m_il{ il },
m_mj2{ mj2 },
m_m{},
m_zmean{},
m_sum{},
m_aicm{},
m_sdm{},
m_aicb{},
m_sdb{},
m_ek{},
m_oeic{},
m_omean{},
m_om{},
m_ind(ival,   static_cast<size_t>(m_k)),
m_a1(dval,    static_cast<size_t>(m_k)),
m_sd(dval,    static_cast<size_t>(m_k+1)),
m_aic(dval,   static_cast<size_t>(m_k+1)),
m_dic(dval,   static_cast<size_t>(m_k+1)),
m_a2(dval,    static_cast<size_t>(m_k)),
m_c(dval,     static_cast<size_t>(m_k)),
m_c1(dval,    static_cast<size_t>(m_k+1)),
m_c2(dval,    static_cast<size_t>(m_k)),
m_b(dval,     static_cast<size_t>(m_k)),
m_esum(dval,  static_cast<size_t>(m_k+1)),
m_e(dval,     static_cast<size_t>(m_n*m_il)),
m_emean(dval, static_cast<size_t>(m_il)),
m_vari(dval,  static_cast<size_t>(m_il)),
m_skew(dval,  static_cast<size_t>(m_il)),
m_peak(dval,  static_cast<size_t>(m_il)),
m_cov(dval,   static_cast<size_t>(m_mj2)),
m_sxx(dval,   static_cast<size_t>(121)) {}

gms::math
::BayesEstimate
::BayesEstimate(const BayesEstimate &x)
:
m_k{     x.m_k },
m_n{     x.m_n },
m_il{    x.m_il },
m_mj2{   x.m_mj2 },
m_m{     x.m_m},
m_zmean{ x.m_zmean},
m_sum{   x.m_sum},
m_aicm{  x.m_aicm},
m_sdm{   x.m_sdm},
m_aicb{  x.m_aicb},
m_sdb{   x.m_sdb},
m_ek{    x.m_ek},
m_oeic{  x.m_oeic},
m_omean{ x.m_omean },
m_om{    x.m_om},
m_ind(   x.m_ind),
m_a1(    x.m_a1),
m_sd(    x.m_sd),
m_aic(   x.m_aic),
m_dic(   x.m_dic),
m_a2(    x.m_a2),
m_c(     x.m_c),
m_c1(    x.m_c1),
m_c2(    x.m_c2),
m_b(     x.m_b),
m_esum(  x.m_esum),
m_e(     x.m_e),
m_emean( x.m_emean),
m_vari(  x.m_vari),
m_skew(  x.m_skew),
m_peak(  x.m_peak),
m_cov(   x.m_cov),
m_sxx(   x.m_sxx)  {}

gms::math
::BayesEstimate
::BayesEstimate(BayesEstimate &&x)
:
m_k{     std::move(x.m_k) },
m_n{     std::move(x.m_n) },
m_il{    std::move(x.m_il) },
m_mj2{   std::move(x.m_mj2) },
m_m{     std::move(x.m_m) },
m_zmean{ std::move(x.m_zmean) },
m_sum{   std::move(x.m_sum) },
m_aicm{  std::move(x.m_aicm) },
m_sdm{   std::move(x.m_sdm) },
m_aicb{  std::move(x.m_aicb) },
m_sdb{   std::move(x.m_sdb) },
m_ek{    std::move(x.m_ek) },
m_oeic{  std::move(x.m_oeic) },
m_omean{ std::move(x.m_omean) },
m_om{  std::move(x.m_om)},
m_ind( std::move(x.m_ind)),
m_a1(  std::move(x.m_a1)),
m_sd(  std::move(x.m_sd)),
m_aic( std::move(x.m_aic)),
m_dic( std::move(x.m_dic)),
m_a2(  std::move(x.m_a2)),
m_c(   std::move(x.m_c)),
m_c1(  std::move(x.m_c1)),
m_c2(  std::move(x.m_c2)),
m_b(   std::move(x.m_b)),
m_esum(std::move(x.m_esum)),
m_e(   std::move(x.m_e)),
m_emean(std::move(x.m_emean)),
m_vari(std::move(x.m_vari)),
m_skew(std::move(x.m_skew)),
m_peak(std::move(x.m_peak)),
m_cov(std::move(x.m_cov)),
m_sxx(std::move(x.m_sxx)) {}

gms::math::BayesEstimate &
gms::math::BayesEstimate
::operator=(const BayesEstimate &x) {

	if (this == &x) return (*this);

	BayesEstimate tmp{x};
	std::swap(*this,tmp);
	
	return (*this);
}

gms::math::BayesEstimate &
gms::math::BayesEstimate
::operator=(BayesEstimate &&x) {

	if (this == &x) return (*this);

	*this = std::move(x);
	return (*this);
}
	


void
gms::math::
::BayesEstimate
::compute(VAf64 & zs,
	  int32_t imodel,
	  int32_t lag,
	  VAi32 & lg1,
	  VAi32 & lg2) {
	 

	

	MOD_TIMSAC_mp_BSUBSTF( &zs[0],
						   &m_n,
						   &imodel,
						   &lag,
						   &m_k,
						   &m_il,
						   &lg1[0],
						   &lg2[0],
						   &m_zmean,
						   &m_sum,
						   &m_m,
						   &m_aicm,
						   &m_sdm,
						   &m_a1[0],
						   &m_sd[0],
						   &m_aic[0],
						   &m_dic[0],
						   &m_aicb,
						   &m_sdb,
						   &m_ek,
						   &m_a2[0],
						   &m_ind[0],
						   &m_c[0],
						   &m_c1[0],
						   &m_c2[0],
						   &m_b[0],
						   &m_oeic,
						   &m_esum[0],
						   &m_omean,
						   &m_om,
						   &m_e[0],
						   &m_emean[0],
						   &m_vari[0],
						   &m_skew[0],
						   &m_peak[0],
						   &m_cov[0],
						   &m_sxx[0]);

	

}



