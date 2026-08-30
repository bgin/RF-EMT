
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

#include "GMS_kmeans.h"
#include "GMS_kmeans_wrapper.h"
#include "GMS_malloc.h"
#include "GMS_indices.h"


 		  

using namespace gms::common;

gms::math
::KMEANS::KMEANS()
:
m_dim_num{},
m_samples_num{},
m_cluster_num{},
m_iter_max{},
m_iter_num{},
m_seed{},
m_cluster(NULL),
m_cluster_center(NULL),
m_cluster_population(NULL),
m_cluster_energy(NULL) {}

gms::math
::KMEANS::KMEANS(
const int32_t dim_num,
const int32_t samples_num,
const int32_t cluster_num,
const int32_t iter_max,
const int32_t seed)
:
m_dim_num{ dim_num },
m_samples_num{ samples_num },
m_cluster_num{ cluster_num },
m_iter_max{ iter_max },
m_iter_num{},
m_seed{ seed },
m_cluster(reinterpret_cast<int32_t*>(gms_mm_malloc(m_samples_num*sizeof(int32_t), 64ULL))),
m_cluster_center(reinterpret_cast<double*>(gms_mm_malloc((m_dim_num*m_cluster_num)*sizeof(double), 64ULL))),
m_cluster_population(reinterpret_cast<int32_t*>(gms_mm_malloc(m_cluster_num*sizeof(int32_t), 64ULL))),
m_cluster_energy(reinterpret_cast<double*>(gms_mm_malloc(m_cluster_num*sizeof(double), 64ULL))) {
	
#pragma simd vectorlengthfor(double) vecremainder
	for (int32_t i = 0; i != m_samples_num; ++i)
		m_cluster[i] = 0;

	for (int32_t i = 0; i != m_dim_num; ++i)
#pragma simd vectorlengthfor(double) vecremainder
	for (int32_t j = 0; j != m_cluster_num; ++j)
		m_cluster_center[Ix2D(i, m_dim_num, j)] = 0.0;

#pragma simd vectorlengthfor(double) vecremainder
	for (int32_t i = 0; i != m_cluster_num; ++i) {
		m_cluster_population[i] = 0;
		m_cluster_energy[i] = 0.0;
	}
}

gms::math
::KMEANS::KMEANS(const KMEANS &x)
:
m_dim_num{ x.m_dim_num },
m_samples_num{ x.m_samples_num },
m_cluster_num{ x.m_cluster_num },
m_iter_max{ x.m_iter_max },
m_iter_num{ x.m_iter_num },
m_seed{ x.m_seed },
m_cluster(reinterpret_cast<int32_t*>(gms_mm_malloc(m_samples_num*sizeof(int32_t), 64ULL))),
m_cluster_center(reinterpret_cast<double*>(gms_mm_malloc((m_dim_num*m_cluster_num)*sizeof(double), 64ULL))),
m_cluster_population(reinterpret_cast<int32_t*>(gms_mm_malloc(m_cluster_num*sizeof(int32_t), 64ULL))),
m_cluster_energy(reinterpret_cast<double*>(gms_mm_malloc(m_cluster_num*sizeof(double), 64ULL))) {
	
#pragma simd vectorlengthfor(double) vecremainder
	for (int32_t i = 0; i != m_samples_num; ++i)
		m_cluster[i] = x.m_cluster[i];

	for (int32_t i = 0; i != m_dim_num; ++i)
#pragma simd vectorlengthfor(double) vecremainder
	for (int32_t j = 0; j != m_cluster_num; ++j)
		m_cluster_center[Ix2D(i, m_cluster_num, j)] = x.m_cluster_center[Ix2D(i, m_cluster_num, j)];

#pragma simd vectorlengthfor(double) vecremainder
	for (int32_t i = 0; i != m_cluster_num; ++i) {
		m_cluster_population[i] = x.m_cluster_population[i];
		m_cluster_energy[i] = x.m_cluster_energy[i];
	}
}

gms::math
::KMEANS::KMEANS(KMEANS &&x)
:
m_dim_num{ x.m_dim_num },
m_samples_num{ x.m_samples_num },
m_cluster_num{ x.m_cluster_num },
m_iter_max{ x.m_iter_max },
m_iter_num{ x.m_iter_num },
m_seed{ x.m_seed },
m_cluster(NULL),
m_cluster_center(NULL),
m_cluster_population(NULL),
m_cluster_energy(NULL) {
	m_cluster = &x.m_cluster[0];
	m_cluster_center = &x.m_cluster_center[0];
	m_cluster_population = &x.m_cluster_population[0];
	m_cluster_energy = &x.m_cluster_energy[0];

	x.m_cluster_num = 0;
	x.m_dim_num = 0;
	x.m_samples_num = 0;
	x.m_iter_max = 0;
	x.m_iter_num = 0;
	x.m_seed = 0;
	x.m_cluster = NULL;
	x.m_cluster_center = NULL;
	x.m_cluster_population = NULL;
	x.m_cluster_energy = NULL;
}

gms::math
::KMEANS::~KMEANS() {
	gms_mm_free(m_cluster_energy);
        m_cluster_energy = NULL;
	gms_mm_free(m_cluster_population);
	m_cluster_population = NULL;
	gms_mm_free(m_cluster_center);
	m_cluster_center = NULL;
	gms_mm_free(m_cluster);
	m_cluster = NULL;
}

gms::math::KMEANS &
gms::math::KMEANS::operator=(const KMEANS &x) {
	if (this == &x) return (*this);
	m_dim_num = x.m_dim_num;
	m_samples_num = x.m_samples_num;
	m_cluster_num = x.m_cluster_num;
	m_iter_max = x.m_iter_max;
	m_iter_num = x.m_iter_num;
	m_seed = x.m_seed;
	_mm_free(m_cluster_energy);
	_mm_free(m_cluster_population);
	_mm_free(m_cluster_center);
	_mm_free(m_cluster);
	int32_t * __restrict cluster =
		reinterpret_cast<int32_t*>(gms_mm_malloc(m_samples_num*sizeof(int32_t), 64ULL));
	double * __restrict cluster_center =
		reinterpret_cast<double*>(gms_mm_malloc((m_dim_num*m_cluster_num)*sizeof(double), 64ULL));
	int32_t * __restrict cluster_population =
		reinterpret_cast<int32_t*>(gms_mm_malloc(m_cluster_num*sizeof(int32_t), 64ULL));
	double * __restrict cluster_energy =
		reinterpret_cast<double*>(gms_mm_malloc(m_cluster_num*sizeof(double), 64ULL));
	
#pragma simd vectorlengthfor(double) vecremainder
	for (int32_t i = 0; i != x.m_samples_num; ++i)
		cluster[i] = x.m_cluster[i];

	for (int32_t i = 0; i != x.m_dim_num; ++i)
#pragma simd vectorlengthfor(double)  vecremainder
	for (int32_t j = 0; j != x.m_cluster_num; ++j)
		cluster_center[Ix2D(i, x.m_dim_num, j)] = x.m_cluster_center[Ix2D(i, x.m_dim_num, j)];
#pragma simd vectorlengthfor(double)  vecremainder
	for (int32_t i = 0; i != x.m_cluster_num; ++i) {
		cluster_population[i] = x.m_cluster_population[i];
		cluster_energy[i] = x.m_cluster_energy[i];
	}
	m_cluster = cluster;
	m_cluster_center = cluster_center;
	m_cluster_population = cluster_population;
	m_cluster_energy = cluster_energy;
	return (*this);
}

gms::math::KMEANS &
gms::math::KMEANS::operator=(KMEANS &&x) {
	if (this == &x) return (*this);
	m_cluster_num = x.m_cluster_num;
	m_dim_num = x.m_dim_num;
	m_samples_num = x.m_samples_num;
	m_iter_max = x.m_iter_max;
	m_iter_num = x.m_iter_num;
	m_seed = x.m_seed;
	gms_mm_free(m_cluster_energy);
	gms_mm_free(m_cluster_population);
	gms_mm_free(m_cluster_center);
	gms_mm_free(m_cluster);
	m_cluster = &x.m_cluster[0];
	m_cluster_center = &x.m_cluster_center[0];
	m_cluster_population = &x.m_cluster_population[0];
	m_cluster_energy = &x.m_cluster_energy[0];

	x.m_cluster_num = 0;
	x.m_dim_num = 0;
	x.m_samples_num = 0;
	x.m_iter_max = 0;
	x.m_iter_num = 0;
	x.m_seed = 0;
	x.m_cluster = NULL;
	x.m_cluster_center = NULL;
	x.m_cluster_population = NULL;
	x.m_cluster_energy = NULL;
	return (*this);
}


void
gms::math::
compute_hmeans01(KMEANS &hm01,
		 double * __restrict sig_samples,
		 double * __restrict cluster_variance) {
		
	
	
	cluster_initialize_5(&hm01.m_dim_num, &hm01.m_samples_num, &hm01.m_cluster_num,
		&sig_samples[0], &hm01.m_seed, &hm01.m_cluster_center[0]);

	hmeans_01(&hm01.m_dim_num, &hm01.m_samples_num, &hm01.m_cluster_num, &hm01.m_iter_max,
		&hm01.m_iter_num, &sig_samples[0], &hm01.m_cluster[0], &hm01.m_cluster_center[0],
		&hm01.m_cluster_population[0], &hm01.m_cluster_energy[0]);

	cluster_variance_compute(&hm01.m_dim_num, &hm01.m_samples_num, &hm01.m_cluster_num,
		&sig_samples[0], &hm01.m_cluster[0], &hm01.m_cluster_center[0], &cluster_variance[0]);

	cluster_print_summary(&hm01.m_samples_num, &hm01.m_cluster_num, &hm01.m_cluster_population[0],
		&hm01.m_cluster_energy[0], &cluster_variance[0]);
	
}

void
gms::math::compute_hmeans02(KMEANS &hm02,
						         double * __restrict sig_samples,
							  double * __restrict cluster_variance) {
						          

	
	cluster_initialize_5(&hm02.m_dim_num, &hm02.m_samples_num, &hm02.m_cluster_num,
		&sig_samples[0], &hm02.m_seed, &hm02.m_cluster_center[0]);

	hmeans_02(&hm02.m_dim_num, &hm02.m_samples_num, &hm02.m_cluster_num, &hm02.m_iter_max,
		&hm02.m_iter_num, &sig_samples[0], &hm02.m_cluster[0], &hm02.m_cluster_center[0],
		&hm02.m_cluster_population[0], &hm02.m_cluster_energy[0], &hm02.m_seed);

	cluster_variance_compute(&hm02.m_dim_num, &hm02.m_samples_num, &hm02.m_cluster_num,
		&sig_samples[0], &hm02.m_cluster[0], &hm02.m_cluster_center[0], &cluster_variance[0]);

	cluster_print_summary(&hm02.m_samples_num, &hm02.m_cluster_num, &hm02.m_cluster_population[0],
		&hm02.m_cluster_energy[0], &cluster_variance[0]);
	
}

void
gms::math::compute_kmeans01(       KMEANS &km01,
					 double * __restrict sig_samples,
					 double * __restrict cluster_variance) {
				        

	
	cluster_initialize_5(&km01.m_dim_num,&km01.m_samples_num,&km01.m_cluster_num,
		&sig_samples[0], &km01.m_seed, &km01.m_cluster_center[0]);

	kmeans_01(&km01.m_dim_num,&km01.m_samples_num,&km01.m_cluster_num,&km01.m_iter_max,
		&km01.m_iter_num, &sig_samples[0], &km01.m_cluster[0], &km01.m_cluster_center[0],
		&km01.m_cluster_population[0], &km01.m_cluster_energy[0]);

	cluster_variance_compute(&km01.m_dim_num, &km01.m_cluster_num, &km01.m_cluster_num,
							 &sig_samples[0],&km01.m_cluster[0],&km01.m_cluster_center[0],
							 &cluster_variance[0] );

	cluster_print_summary(&km01.m_samples_num, &km01.m_cluster_num, &km01.m_cluster_population[0],
		&km01.m_cluster_energy[0], &cluster_variance[0]);
	
}

void
gms::math::compute_kmeans02(       KMEANS &km02,
				         double * __restrict sig_samples,
				         double * __restrict cluster_variance) {
				       

	
	cluster_initialize_1(&km02.m_dim_num,&km02.m_samples_num,&km02.m_cluster_num,
						 &sig_samples[0], &km02.m_cluster_center[0]);

	kmeans_02(&km02.m_dim_num,&km02.m_samples_num,&km02.m_cluster_num,&km02.m_iter_max,
		&km02.m_iter_num, &sig_samples[0], &km02.m_cluster[0], &km02.m_cluster_center[0],
		&km02.m_cluster_population[0], &km02.m_cluster_energy[0]);

	cluster_variance_compute(&km02.m_dim_num,&km02.m_samples_num,&km02.m_cluster_num,
		&sig_samples[0], &km02.m_cluster[0], &km02.m_cluster_center[0],&cluster_variance[0]);

	cluster_print_summary(&km02.m_samples_num, &km02.m_cluster_num, &km02.m_cluster_population[0],
		&km02.m_cluster_energy[0], &cluster_variance[0]);
	
}

void
gms::math::compute_kmeans03(       KMEANS &km03,
				         double * __restrict sig_samples,
				         double * __restrict cluster_variance) {
				        

	
	cluster_initialize_1(&km03.m_dim_num,&km03.m_samples_num,&km03.m_cluster_num,
							 &sig_samples[0],&km03.m_cluster_center[0]);
	
	kmeans_03(&km03.m_dim_num,&km03.m_samples_num,&km03.m_cluster_num,&km03.m_iter_max,
		&km03.m_iter_num, &sig_samples[0], &km03.m_cluster[0], &km03.m_cluster_center[0],
		&km03.m_cluster_population[0], &km03.m_cluster_energy[0]);

	cluster_variance_compute(&km03.m_dim_num,&km03.m_samples_num,&km03.m_cluster_num,
		&sig_samples[0], &km03.m_cluster[0], &km03.m_cluster_center[0], &cluster_variance[0]);

	cluster_print_summary(&km03.m_samples_num, &km03.m_cluster_num, &km03.m_cluster_population[0],
		&km03.m_cluster_energy[0], &cluster_variance[0]);

	
}

void
gms::math::
compute_hmeans01_and_kmeans01(           KMEANS &km,
				         double * __restrict sig_samples,
				         double * __restrict cluster_variance) {
				        

	

	cluster_initialize_5(&km.m_dim_num,&km.m_samples_num,&km.m_cluster_num,
		&sig_samples[0], &km.m_seed, &km.m_cluster_center[0]);

	hmeans_01(&km.m_dim_num,&km.m_samples_num,&km.m_cluster_num,&km.m_iter_max,
		&km.m_iter_num, &sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0],
		&km.m_cluster_population[0], &km.m_cluster_energy[0]);

	cluster_variance_compute(&km.m_dim_num,&km.m_samples_num,&km.m_cluster_num,
		&sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0], &cluster_variance[0]);

	cluster_print_summary(&km.m_samples_num, &km.m_cluster_num, &km.m_cluster_population[0],
		&km.m_cluster_energy[0], &cluster_variance[0]);

	kmeans_01(&km.m_dim_num,&km.m_samples_num,&km.m_cluster_num,&km.m_iter_max,&km.m_iter_num,
			  &sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0],
			  &km.m_cluster_population[0], &km.m_cluster_energy[0]);

	cluster_variance_compute(&km.m_dim_num,&km.m_samples_num,&km.m_cluster_num,
		&sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0], &cluster_variance[0]);

	cluster_print_summary(&km.m_samples_num, &km.m_cluster_num, &km.m_cluster_population[0],
		&km.m_cluster_energy[0], &cluster_variance[0]);

	
}

void
gms::math::
compute_hmeans01_and_kmeans02(  KMEANS &km,
			        double * __restrict sig_samples,
			        double * __restrict cluster_variance) {
			       

	

	cluster_initialize_5(&km.m_dim_num, &km.m_samples_num, &km.m_cluster_num,
		&sig_samples[0], &km.m_seed, &km.m_cluster_center[0]);

	hmeans_01(&km.m_dim_num, &km.m_samples_num, &km.m_cluster_num, &km.m_iter_max,
		&km.m_iter_num, &sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0],
		&km.m_cluster_population[0], &km.m_cluster_energy[0]);

	cluster_variance_compute(&km.m_dim_num, &km.m_samples_num, &km.m_cluster_num,
		&sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0], &cluster_variance[0]);

	cluster_print_summary(&km.m_samples_num, &km.m_cluster_num, &km.m_cluster_population[0],
		&km.m_cluster_energy[0], &cluster_variance[0]);

	kmeans_02(&km.m_dim_num,&km.m_samples_num,&km.m_cluster_num,&km.m_iter_max,
		&km.m_iter_num, &sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0],
		&km.m_cluster_population[0], &km.m_cluster_energy[0]);

	cluster_variance_compute(&km.m_dim_num, &km.m_samples_num, &km.m_cluster_num,
		&sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0], &cluster_variance[0]);

	cluster_print_summary(&km.m_samples_num, &km.m_cluster_num, &km.m_cluster_population[0],
		&km.m_cluster_energy[0], &cluster_variance[0]);

	
}

void
gms::math::
compute_hmeans01_and_kmeans03(  KMEANS &km,
			        double * __restrict sig_samples,
			        double * __restrict cluster_variance) {
			       

	

	cluster_initialize_5(&km.m_dim_num, &km.m_samples_num, &km.m_cluster_num,
		&sig_samples[0], &km.m_seed, &km.m_cluster_center[0]);

	hmeans_01(&km.m_dim_num, &km.m_samples_num, &km.m_cluster_num, &km.m_iter_max,
		&km.m_iter_num, &sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0],
		&km.m_cluster_population[0], &km.m_cluster_energy[0]);

	cluster_variance_compute(&km.m_dim_num, &km.m_samples_num, &km.m_cluster_num,
		&sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0], &cluster_variance[0]);

	cluster_print_summary(&km.m_samples_num, &km.m_cluster_num, &km.m_cluster_population[0],
		&km.m_cluster_energy[0], &cluster_variance[0]);

	kmeans_03(&km.m_dim_num, &km.m_samples_num, &km.m_cluster_num, &km.m_iter_max,
		&km.m_iter_num, &sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0],
		&km.m_cluster_population[0], &km.m_cluster_energy[0]);

	cluster_variance_compute(&km.m_dim_num, &km.m_samples_num, &km.m_cluster_num,
		&sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0], &cluster_variance[0]);

	cluster_print_summary(&km.m_samples_num, &km.m_cluster_num, &km.m_cluster_population[0],
		&km.m_cluster_energy[0], &cluster_variance[0]);

	

}

void
gms::math::
compute_hmeans_w_01(     KMEANS &km,
		         double * __restrict sig_samples,
		         double * __restrict cluster_variance,
		         double * __restrict weight) {
		         

	
	cluster_initialize_5(&km.m_dim_num,&km.m_samples_num,&km.m_cluster_num,
		&sig_samples[0], &km.m_seed, &km.m_cluster_center[0]);

	hmeans_w_01(&km.m_dim_num,&km.m_samples_num,&km.m_cluster_num,&km.m_iter_max,
		&km.m_iter_num, &sig_samples[0], &weight[0], &km.m_cluster[0],
		&km.m_cluster_center[0], &km.m_cluster_population[0], &km.m_cluster_energy[0]);

	cluster_variance_compute(&km.m_dim_num,&km.m_samples_num,&km.m_cluster_num,
		&sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0], &cluster_variance[0]);

	cluster_print_summary(&km.m_dim_num, &km.m_cluster_num, &km.m_cluster_population[0],
		&km.m_cluster_energy[0], &cluster_variance[0]);

	
}

void
gms::math::
compute_hmeans_w_02(     KMEANS &km,
		         double * __restrict sig_samples,
			 double * __restrict cluster_variance,
		         double * weight) {
		        


	cluster_initialize_5(&km.m_dim_num, &km.m_samples_num, &km.m_cluster_num,
		&sig_samples[0], &km.m_seed, &km.m_cluster_center[0]);

	hmeans_w_02(&km.m_dim_num, &km.m_samples_num, &km.m_cluster_num, &km.m_iter_max,
		&km.m_iter_num, &sig_samples[0], &weight[0], &km.m_cluster[0],
		&km.m_cluster_center[0], &km.m_cluster_population[0], &km.m_cluster_energy[0],&km.m_seed);

	cluster_variance_compute(&km.m_dim_num, &km.m_samples_num, &km.m_cluster_num,
		&sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0], &cluster_variance[0]);

	cluster_print_summary(&km.m_dim_num, &km.m_cluster_num, &km.m_cluster_population[0],
		&km.m_cluster_energy[0], &cluster_variance[0]);

	
}

void
gms::math::
compute_kmeans_w_01(     KMEANS &km,
		         double * __restrict sig_samples,
		         double * __restrict cluster_variance,
		         double * weight) {
		       

	

	cluster_initialize_5(&km.m_dim_num, &km.m_samples_num, &km.m_cluster_num,
		&sig_samples[0], &km.m_seed, &km.m_cluster_center[0]);

	kmeans_w_01(&km.m_dim_num,&km.m_samples_num,&km.m_cluster_num,&km.m_iter_max,
		&km.m_iter_num, &sig_samples[0], &weight[0], &km.m_cluster[0], &km.m_cluster_center[0],
		&km.m_cluster_population[0], &km.m_cluster_energy[0]);

	cluster_variance_compute(&km.m_dim_num, &km.m_samples_num, &km.m_cluster_num,
		&sig_samples[0], &km.m_cluster[0], &km.m_cluster_center[0], &cluster_variance[0]);

	cluster_print_summary(&km.m_dim_num, &km.m_cluster_num, &km.m_cluster_population[0],
		&km.m_cluster_energy[0], &cluster_variance[0]);

	
}


