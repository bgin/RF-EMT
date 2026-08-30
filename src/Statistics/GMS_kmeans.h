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
#ifndef __GMS_KMEANS_H__
#define __GMS_KMEANS_H__ 061020191431



namespace file_info {


	const unsigned int GMS_KMEANS_MAJOR = 1U;

	const unsigned int GMS_KMEANS_MINOR = 0U;

	const unsigned int GMS_KMEANS_MICRO = 0U;

	const unsigned int GMS_KMEANS_FULLVER = 
	 1000U*GMS_KMEANS_MAJOR+100U*GMS_KMEANS_MINOR+10U*GMS_KMEANS_MICRO;

	const char * const pgGMS_KMEANS_CREATE_DATE = "06-10-2019 14:31 +00200 (SUN 06 OCT 2019 GMT+2)";

	const char * const pgGMS_KMEANS_BUILD_DATE = __DATE__ ":" __TIME__;

	const char * const pgGMS_KMEANS_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

	const char * const pgGMS_KMEANS_SYNOPSIS = "High-level wrappers for Fortran 90 KMEANS library."; 
}

#include <cstdint>
#include "GMS_config.h"


namespace gms {
	namespace math {
	
		

			__ATTR_ALIGN__(64) struct KMEANS {

				
				/*
				On input, the user
				may specify an initial cluster for each point, or leave all entrie of
				CLUSTER set to 0.  On output, CLUSTER contains the index of the
				cluster to which each data point belongs.
				*/

			        int32_t * __restrict m_cluster;
				/*
				the coordinates of the cluster centers.
				*/

			        double * __restrict m_cluster_center;
				/*
				the number of
				points assigned to each cluster.
				*/

			        int32_t * __restrict m_cluster_population;
				/*
				the energy of
				the clusters.
				*/

			        double * __restrict m_cluster_energy;
			         
			         // the number of spatial dimensions.
				int32_t m_dim_num;
				// the number of signal samples (re,im) components
				int32_t m_samples_num;
				// the number of clusters.
				int32_t m_cluster_num;
				// the maximum number of iterations.
				int32_t m_iter_max;
				// the number of iterations taken.
				int32_t m_iter_num;
				// a seed for the random number generator.
				int32_t m_seed;
				PAD_TO(0,8)

				//
				//	Default Constructor implemented, although it does not have
				//  any statistical meaning.
				//
			    
				KMEANS();

				//
				//  Main class Constructor - responsible for scalar member initialization and
				//	allocation and zero-initialzation of array members.
				// 
				KMEANS(	const int32_t,const int32_t,
					const int32_t,const int32_t,
					const int32_t);

				//
				// Copy Constructor -- deep copy semantics
				//
				KMEANS(const KMEANS &);

				//
				// Move Constructor -- shallow copy semantics
				//
				KMEANS( KMEANS &&) noexcept(true);

				~KMEANS() noexcept(true);

				KMEANS & operator=(const KMEANS &);

				KMEANS & operator=(KMEANS &&) noexcept(true);



			};

			//
			//	Run computation of H-MEANS01
			//  Second argument it is radar samples signal whose number of domain is always 2
			//  and it contains real and imaginary part of signal sample.
			//  Signal samples:
			/*		Re			Im
			0.532095    0.894141
			0.189043    0.613426
			0.550977    0.415724
			0.003973    0.605760
			0.894230    0.666812
			0.073072    0.740658
			0.640180    0.926186
			0.389914    0.553149
			0.046918    0.172275
			0.820327    0.578957
			*/
			//  Third argument is storage for variance computation and must be zero-filled.
			//
			//  Warning:
			//             second argument (from left) must have size of (dim_num*samples_num)
			//			   third argument (from left) must have size of cluster_num
			//				
			//  
			void compute_hmeans01(KMEANS &,
					      double * __restrict,
					      double * __restrict);
					    

			//
			//	Run computation of H-MEANS02
			//  Second argument it is radar samples signal whose number of domain is always 2
			//  and it contains real and imaginary part of signal sample.
			//  Signal samples:
			/*		Re			Im
			0.532095    0.894141
			0.189043    0.613426
			0.550977    0.415724
			0.003973    0.605760
			0.894230    0.666812
			0.073072    0.740658
			0.640180    0.926186
			0.389914    0.553149
			0.046918    0.172275
			0.820327    0.578957
			*/
			//  Third argument is storage for variance computation and must be zero-filled.
			//
			//  Warning:
			//             second argument (from left) must have size of (dim_num*samples_num)
			//			   third argument (from left) must have size of cluster_num
			//				
			//  

			void compute_hmeans02(KMEANS &, 
					      double * __restrict,
					      double * __restrict);
					     

			//
			//	Run subroutine K-MEANS01
			//  Second argument it is radar samples signal whose number of domain is always 2
			//  and it contains real and imaginary part of signal sample.
			//  Signal samples:
			/*		Re			Im
			0.532095    0.894141
			0.189043    0.613426
			0.550977    0.415724
			0.003973    0.605760
			0.894230    0.666812
			0.073072    0.740658
			0.640180    0.926186
			0.389914    0.553149
			0.046918    0.172275
			0.820327    0.578957
			*/
			//  Third argument is storage for variance computation and must be zero-filled.
			//
			//  Warning:
			//             second argument (from left) must have size of (dim_num*samples_num)
			//			   third argument (from left) must have size of cluster_num
			//				
			//  

			void compute_kmeans01(KMEANS &, 
					      double * __restrict,
					      double * __restrict);
					      

			//
			//	Run subroutine K-MEANS02
			//  Second argument it is radar samples signal whose number of domain is always 2
			//  and it contains real and imaginary part of signal sample.
			//  Signal samples:
			/*		Re			Im
			0.532095    0.894141
			0.189043    0.613426
			0.550977    0.415724
			0.003973    0.605760
			0.894230    0.666812
			0.073072    0.740658
			0.640180    0.926186
			0.389914    0.553149
			0.046918    0.172275
			0.820327    0.578957
			*/
			//  Third argument is storage for variance computation and must be zero-filled.
			//
			//  Warning:
			//             second argument (from left) must have size of (dim_num*samples_num)
			//			   third argument (from left) must have size of cluster_num
			//				
			//  

			void compute_kmeans02(KMEANS &,
					      double * __restrict,
					      double * __restrict);
					      

			/*
			//
			//	Run subroutine K-MEANS03
			//  Second argument it is radar samples signal whose number of domain is always 2
			//  and it contains real and imaginary part of signal sample.
			//  Signal samples:
			/*		Re			Im
			0.532095    0.894141
			0.189043    0.613426
			0.550977    0.415724
			0.003973    0.605760
			0.894230    0.666812
			0.073072    0.740658
			0.640180    0.926186
			0.389914    0.553149
			0.046918    0.172275
			0.820327    0.578957
			*/
			//  Third argument is storage for variance computation and must be zero-filled.
			//
			//  Warning:
			//             second argument (from left) must have size of (dim_num*samples_num)
			//			   third argument (from left) must have size of cluster_num
			//				
			//  
	        
			void compute_kmeans03(KMEANS &,
					      double * __restrict,
					      double * __restrict);
					      

			
			/*
			//
			//	Run subroutines H-MEANS01 and K-MEANS01
			//  Second argument it is radar samples signal whose number of domain is always 2
			//  and it contains real and imaginary part of signal sample.
			//  Signal samples:
			/*		Re			Im
			0.532095    0.894141
			0.189043    0.613426
			0.550977    0.415724
			0.003973    0.605760
			0.894230    0.666812
			0.073072    0.740658
			0.640180    0.926186
			0.389914    0.553149
			0.046918    0.172275
			0.820327    0.578957
			*/
			//  Third argument is storage for variance computation and must be zero-filled.
			//
			//  Warning:
			//             second argument (from left) must have size of (dim_num*samples_num)
			//			   third argument (from left) must have size of cluster_num
			//				
			//  
			
			void compute_hmeans01_and_kmeans01(KMEANS &,
							   double * __restrict,
							   double * __restrict);
							  

			/*
			//
			//	Run subroutines H-MEANS01 and K-MEANS02
			//  Second argument it is radar samples signal whose number of domain is always 2
			//  and it contains real and imaginary part of signal sample.
			//  Signal samples:
			/*		Re			Im
			0.532095    0.894141
			0.189043    0.613426
			0.550977    0.415724
			0.003973    0.605760
			0.894230    0.666812
			0.073072    0.740658
			0.640180    0.926186
			0.389914    0.553149
			0.046918    0.172275
			0.820327    0.578957
			*/
			//  Third argument is storage for variance computation and must be zero-filled.
			//
			//  Warning:
			//             second argument (from left) must have size of (dim_num*samples_num)
			//			   third argument (from left) must have size of cluster_num
			//				
			//  

			void compute_hmeans01_and_kmeans02(KMEANS &,
							   double * __restrict,
							   double * __restrict);
							  

			/*
			//
			//	Run subroutines H-MEANS01 and K-MEANS03
			//  Second argument it is radar samples signal whose number of domain is always 2
			//  and it contains real and imaginary part of signal sample.
			//  Signal samples:
			/*		Re			Im
			0.532095    0.894141
			0.189043    0.613426
			0.550977    0.415724
			0.003973    0.605760
			0.894230    0.666812
			0.073072    0.740658
			0.640180    0.926186
			0.389914    0.553149
			0.046918    0.172275
			0.820327    0.578957
			*/
			//  Third argument is storage for variance computation and must be zero-filled.
			//
			//  Warning:
			//             second argument (from left) must have size of (dim_num*samples_num)
			//			   third argument (from left) must have size of cluster_num
			//				
			//  

			void compute_hmeans01_and_kmeans03(KMEANS &,
							   double * __restrict,
							   double * __restrict);
							 

			/*
			//
			//	Run subroutines HMEANS_W_01 
			//  Second argument it is radar samples signal whose number of domain is always 2
			//  and it contains real and imaginary part of signal sample.
			//  Signal samples:
			/*		Re			Im
			0.532095    0.894141
			0.189043    0.613426
			0.550977    0.415724
			0.003973    0.605760
			0.894230    0.666812
			0.073072    0.740658
			0.640180    0.926186
			0.389914    0.553149
			0.046918    0.172275
			0.820327    0.578957
			*/
			//  Third argument is storage for variance computation and must be zero-filled.
			//
			//  Warning:
			//             second argument (from left) must have size of (dim_num*samples_num)
			//			   third argument (from left) must have size of cluster_num
			//			   fourth argument (from left) must have size of (samples_num)
			//  

			void compute_hmeans_w_01(KMEANS &,
						 double * __restrict,
						 double * __restrict,
						 double * __restrict);
						 

			/*
			//
			//	Run subroutines HMEANS_W_02
			//  Second argument it is radar samples signal whose number of domain is always 2
			//  and it contains real and imaginary part of signal sample.
			//  Signal samples:
			/*		Re			Im
			0.532095    0.894141
			0.189043    0.613426
			0.550977    0.415724
			0.003973    0.605760
			0.894230    0.666812
			0.073072    0.740658
			0.640180    0.926186
			0.389914    0.553149
			0.046918    0.172275
			0.820327    0.578957
			*/
			//  Third argument is storage for variance computation and must be zero-filled.
			//
			//  Warning:
			//             second argument (from left) must have size of (dim_num*samples_num)
			//			   third argument (from left) must have size of cluster_num
			//			   fourth argument (from left) must have size of (samples_num)
			//  

			void compute_hmeans_w_02(KMEANS &,
						 double * __restrict,
						 double * __restrict,
						 double * __restrict);
						

			/*
			//
			//	Run subroutines KMEANS_W_01
			//  Second argument it is radar samples signal whose number of domain is always 2
			//  and it contains real and imaginary part of signal sample.
			//  Signal samples:
			/*		Re			Im
			0.532095    0.894141
			0.189043    0.613426
			0.550977    0.415724
			0.003973    0.605760
			0.894230    0.666812
			0.073072    0.740658
			0.640180    0.926186
			0.389914    0.553149
			0.046918    0.172275
			0.820327    0.578957
			*/
			//  Third argument is storage for variance computation and must be zero-filled.
			//
			//  Warning:
			//             second argument (from left) must have size of (dim_num*samples_num)
			//			   third argument (from left) must have size of cluster_num
			//			   fourth argument (from left) must have size of (samples_num)
			//  

			void compute_kmeans_w_01(KMEANS &,
						 double * __restrict,
						 double * __restrict,
						 double * __restrict);
						

		


		}
	}




#endif /*__GMS_KMEANS_H__*/
