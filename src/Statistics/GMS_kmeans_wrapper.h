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
#ifndef __GMS_KMEANS_WRAPPER_H__
#define __GMS_KMEANS_WRAPPER_H__ 151020191603

//
//	F90 kmeans library wrapper
//



namespace file_info {

  
 	const unsigned int GMS_KMEANS_WRAPPER_MAJOR = 1U;

	const unsigned int GMS_KMEANS_WRAPPER_MINOR = 0U;

	const unsigned int GMS_KMEANS_WRAPPER_MICRO = 0U;

	const unsigned int GMS_KMEANS_WRAPPER_FULLVER = 
	 1000U*GMS_KMEANS_WRAPPER_MAJOR+100U*GMS_KMEANS_WRAPPER_MINOR+10U*GMS_KMEANS_WRAPPER_MICRO;

	const char * const GMS_KMEANS_WRAPPER_CREATE_DATE = "15-10-2019 16:03 +00200 (TUE 15 OCT 2019 GMT+2)";

	const char * const GMS_KMEANS_WRAPPER_BUILD_DATE = __DATE__ " " __TIME__;

	const char * const GMS_KMEANS_WRAPPER_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

	const char * const GMS_KMEANS_WRAPPER_SYNOPSIS =   "C wrappers for Fortran 90 K-MEANS Library.";
}



#if defined (GMS_CXX_98) || defined (GMS_CXX_11) || defined (GMS_CXX_14)

	extern "C"  {

		/*
			!*****************************************************************************80
!
!! CLUSTER_ENERGY_COMPUTE computes the energy of the clusters.
!
!  Discussion:
!
!    The cluster energy is defined as the sum of the distance
!    squared from each point to its cluster center.  It is the goal
!    of the H-means and K-means algorithms to find, for a fixed number
!    of clusters, a clustering that minimizes this energy
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of data points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Input, integer ( kind = 4 ) CLUSTER(POINT_NUM), the cluster to which each
!    data point belongs.
!
!    Input, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), the 
!    centers associated with the minimal energy clustering.
!
!    Output, real ( kind = 8 ) CLUSTER_ENERGY(CLUSTER_NUM), the energy
!    associated with each cluster.
		*/

		void cluster_energy_compute(int * , int * , int * , double * __restrict, 
					    int * , double * __restrict, double * __restrict);


		/*
			 CLUSTER_INITIALIZE_1 initializes the clusters to data points.
!
!  Discussion:
!
!    The cluster centers are simply chosen to be the first data points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the coordinates 
!    of the points.
!
!    Output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the coordinates of the cluster centers.
		*/
		void cluster_initialize_1(int * ,int * , int * , double * __restrict, double * __restrict );

		/*
			 CLUSTER_INITIALIZE_2 initializes the cluster centers to random values.
!
!  Discussion:
!
!    In this case, the hyperbox containing the data is computed.
!
!    Then the cluster centers are chosen uniformly at random within
!    this hyperbox.
!
!    Of course, if the data is not smoothly distributed throughout
!    the box, many cluster centers will be isolated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the coordinates 
!    of the points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), 
!    the coordinates of the cluster centers.
		*/
		void cluster_initialize_2(int * , int *, int *, double * __restrict, int *, double * __restrict);

		/*
			CLUSTER_INITIALIZE_3 initializes the cluster centers to random values.
!
!  Discussion:
!
!    In this case, each point is randomly assigned to a cluster, and
!    the cluster centers are then computed as the centroids of the points 
!    in the cluster.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the coordinates 
!    of the points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), 
!    the coordinates of the cluster centers.
!
		*/
		void cluster_initialize_3(int *, int *, int *, double * __restrict, int *, double * __restrict);

		/*
			CLUSTER_INITIALIZE_4 initializes the cluster centers to random values.
!
!  Discussion:
!
!    In this case, each data point is divided randomly among the
!    the cluster centers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the coordinates 
!    of the points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), 
!    the coordinates of the cluster centers.
		*/
		void cluster_initialize_4(int *, int *, int *, double * __restrict, int *, double * __restrict);

		/*
			 CLUSTER_INITIALIZE_5 initializes the cluster centers to random values.
!
!  Discussion:
!
!    In this case, each cluster center is a random convex combination 
!    of the data points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the coordinates 
!    of the points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the coordinates of the cluster centers.
		*/
		void cluster_initialize_5(int *, int *, int *, double * __restrict, int *, double * __restrict);

		/*
			 CLUSTER_PRINT_SUMMARY prints a summary of data about a clustering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, integer ( kind = 4 ) CLUSTER_POPULATION(CLUSTER_NUM), the number of
!    points assigned to each cluster.
!
!    Input, real ( kind = 8 ) CLUSTER_ENERGY(CLUSTER_NUM), the energy of 
!    the clusters.
!
!    Input, real ( kind = 8 ) CLUSTER_VARIANCE(CLUSTER_NUM), the variance of 
!    the clusters.
		*/
		void cluster_print_summary(int *, int *, int * __restrict, double * __restrict, double * __restrict);

		/*
			 CLUSTER_VARIANCE_COMPUTE computes the variance of the clusters.
!
!  Discussion:
!
!    The cluster variance (from the cluster center) is the average of the 
!    sum of the squares of the distances of each point in the cluster to the 
!    cluster center.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of data points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Input, integer ( kind = 4 ) CLUSTER(POINT_NUM), the cluster to which each
!    data point belongs.
!
!    Input, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), the 
!    centers associated with the minimal energy clustering.
!
!    Output, real ( kind = 8 ) CLUSTER_VARIANCE(CLUSTER_NUM), the variance
!    associated with each cluster.
!
		*/
		void cluster_variance_compute(int *, int *, int *, double * __restrict, int * __restrict,
								      double * __restrict, double * __restrict);

		/*
			HMEANS_01 applies the H-Means algorithm.
!
!  Discussion:
!
!    The data for the H-Means problem is a set of N points X in
!    M-dimensions, and a desired number of clusters K.
!
!    The goal is to determine K points Z, called cluster centers, so that
!    if we associate each point X with its nearest Z value, we minimize
!    the standard deviation or cluster energy.  Writing CLUSTER(I) to
!    indicate the index of the nearest cluster center to point X(I), the
!    energy can be written as:
!
!      Energy = Sum ( 1 <= I <= N ) || X(I) - Z(CLUSTER(I)) ||**2
!
!    where
!
!      || X - Z ||**2 = Sum ( 1 <= J <= M ) ( X(J) - Z(J) )**2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wendy Martinez, Angel Martinez,
!    Computational Statistics Handbook with MATLAB,
!    pages 373-376,
!    Chapman and Hall / CRC, 2002.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of data points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Output, integer ( kind = 4 ) IT_NUM, the number of iterations taken.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Input/output, integer ( kind = 4 ) CLUSTER(POINT_NUM).  On input, the user 
!    may specify an initial cluster for each point, or leave all entrie of
!    CLUSTER set to 0.  On output, CLUSTER contains the index of the
!    cluster to which each data point belongs.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), the 
!    centers associated with the minimal energy clustering.
!
!    Output, integer ( kind = 4 ) CLUSTER_POPULATION(CLUSTER_NUM),
!    the populuation of each cluster.
!
!    Output, real ( kind = 8 ) CLUSTER_ENERGY(CLUSTER_NUM), the energy
!    associated with each cluster.
		*/
		void hmeans_01(int *, int *, int *, int *, int *, double * __restrict, int * __restrict,
					   double * __restrict, int * __restrict, double * __restrict);

		/*
			HMEANS_02 applies the H-Means algorithm.
!
!  Discussion:
!
!    This is a simple routine to group a set of points into K clusters,
!    each with a center point, in such a way that the total cluster 
!    energy is minimized.  The total cluster energy is the sum of the
!    squares of the distances of each point to the center of its cluster.
!
!    The algorithm begins with an initial estimate for the cluster centers:
!
!    1. The points are assigned to the nearest cluster centers.
!
!    2. The iteration stops if the total energy has not changed 
!        significantly, or we have reached the maximum number of iterations.
!
!    3. Each cluster center is replaced by the centroid of the points
!       in the cluster.
!
!    4. Return to step 1.
!
!    The algorithm may fail to find the best solution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Output, integer ( kind = 4 ) IT_NUM, the number of iterations taken.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the coordinates 
!    of the points.
!
!    Input/output, integer ( kind = 4 ) CLUSTER(POINT_NUM).  On input, the user 
!    may specify an initial cluster for each point, or leave all entrie of
!    CLUSTER set to 0.  On output, CLUSTER contains the index of the
!    cluster to which each data point belongs.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the coordinates of the cluster centers.
!
!    Output, integer ( kind = 4 ) CLUSTER_POPULATION(CLUSTER_NUM), the number of
!    points assigned to each cluster.
!
!    Output, real ( kind = 8 ) CLUSTER_ENERGY(CLUSTER_NUM), the energy of 
!    the clusters.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
		*/
		void hmeans_02(int *, int *, int *, int *, int *, double * __restrict, int * __restrict,
					   double * __restrict, int * __restrict, double * __restrict, int * );

		/*
			 HMEANS_W_01 applies the weighted H-Means algorithm. 
!
!  Discussion:
!
!    The input data for the weight H-Means problem includes:
!    * a set of N data points X in M dimensions, 
!    * a set of N nonnegative weights W,
!    * a desired number of clusters K.
!    * an initial set of cluster centers Z,
!    * an (optional) initial set of cluster assignments.
!
!    The goal is to determine K points Z, called cluster centers, and
!    to assign each point X(I) to some cluster Z(J), so that we minimize
!    the weighted standard deviation of the distance of each data point
!    to the center of its cluster.  Writing J = CLUSTER(I) to
!    indicate the index of the nearest cluster center Z(J) to the 
!    point X(I), the quantity we are trying to minimize is the sum
!    of the weighted cluster energies E(J), where:
!
!      E(J) = Sum ( 1 <= I <= N ) W(I) * || X(I) - Z(J) ||**2
!
!    Here, we assume that we are using the Euclidean norm, so that
!    
!      || X(I) - Z(J) ||**2 = Sum ( 1 <= K <= M )
!        ( X(I)(K) - Z(J)(K) )**2
!
!    In this notation, X(I)(K) is the K-th spatial component of the
!    I-th data point.
!
!    Note that this routine should give the same results as HMEANS_01
!    in any case in which all the entries of the WEIGHT vector are equal.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wendy Martinez, Angel Martinez,
!    Computational Statistics Handbook with MATLAB,
!    pages 373-376,
!    Chapman and Hall / CRC, 2002.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of data points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Output, integer ( kind = 4 ) IT_NUM, the number of iterations taken.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Input, real ( kind = 8 ) WEIGHT(POINT_NUM), the weights
!    assigned to the data points.  These must be nonnegative, and
!    at least one must be strictly positive.
!
!    Input/output, integer ( kind = 4 ) CLUSTER(POINT_NUM).  On input, the user 
!    may specify an initial cluster for each point, or leave all entrie of
!    CLUSTER set to 0.  On output, CLUSTER contains the index of the
!    cluster to which each data point belongs.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), the
!    centers associated with the minimal energy clustering.
!
!    Output, integer ( kind = 4 ) CLUSTER_POPULATION(CLUSTER_NUM),
!    the populuation of each cluster.
!
!    Output, real ( kind = 8 ) CLUSTER_ENERGY(CLUSTER_NUM), the energy
!    associated with each cluster.
!
		*/
		void hmeans_w_01(int * , int *, int *, int *, int *, double * __restrict,
						 double * __restrict, int * __restrict, double * __restrict, 
						 int * __restrict, double * __restrict);


		/*
			 HMEANS_W_02 applies the weighted H-Means algorithm.
!
!  Discussion:
!
!    The input data for the weight H-Means problem includes:
!    * a set of N data points X in M dimensions, 
!    * a set of N nonnegative weights W,
!    * a desired number of clusters K.
!    * an initial set of cluster centers Z,
!    * an (optional) initial set of cluster assignments.
!
!    The goal is to determine K points Z, called cluster centers, and
!    to assign each point X(I) to some cluster Z(J), so that we minimize
!    the weighted standard deviation of the distance of each data point
!    to the center of its cluster.  Writing J = CLUSTER(I) to
!    indicate the index of the nearest cluster center Z(J) to the 
!    point X(I), the quantity we are trying to minimize is the sum
!    of the weighted cluster energies E(J), where:
!
!      E(J) = Sum ( 1 <= I <= N ) W(I) * || X(I) - Z(J) ||**2
!
!    Here, we assume that we are using the Euclidean norm, so that
!    
!      || X(I) - Z(J) ||**2 = Sum ( 1 <= K <= M )
!        ( X(I)(K) - Z(J)(K) )**2
!
!    In this notation, X(I)(K) is the K-th spatial component of the
!    I-th data point.
!
!    Note that this routine should give the same results as HMEANS_02
!    in any case in which all the entries of the WEIGHT vector are equal.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Output, integer ( kind = 4 ) IT_NUM, the number of iterations taken.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the coordinates 
!    of the points.
!
!    Input/output, integer ( kind = 4 ) CLUSTER(POINT_NUM).  On input, the user 
!    may specify an initial cluster for each point, or leave all entrie of
!    CLUSTER set to 0.  On output, CLUSTER contains the index of the
!    cluster to which each data point belongs.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the coordinates of the cluster centers.
!
!    Output, integer ( kind = 4 ) CLUSTER_POPULATION(CLUSTER_NUM), the number of
!    points assigned to each cluster.
!
!    Output, real ( kind = 8 ) CLUSTER_ENERGY(CLUSTER_NUM), the energy of 
!    the clusters.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
		*/
	
		void hmeans_w_02(int *, int *, int *, int *, int *, double * __restrict, double * __restrict,int * __restrict,
						 double * __restrict, int * __restrict, double * __restrict, int * );

		/*
			I4_UNIFORM returns a scaled pseudorandom I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) I4_UNIFORM, a number between A and B.
!
		*/

		int i4_uniform(int *, int *, int *);

		/*
			  KMEANS_01 applies the K-Means algorithm.
!
!  Discussion:
!
!    Given a matrix of POINT_NUM observations on DIM_NUM variables, the
!    observations are to be allocated to CLUSTER_NUM clusters in such 
!    a way that the within-cluster sum of squares is minimized.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2004
!
!  Author:
!
!    Original FORTRAN77 version by David Sparks.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Sparks,
!    Algorithm AS 58: 
!    Euclidean Cluster Analysis,
!    Applied Statistics,
!    Volume 22, Number 1, 1973, pages 126-130.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Output, integer ( kind = 4 ) IT_NUM, the number of iterations taken.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the points.
!
!    Output, integer ( kind = 4 ) CLUSTER(POINT_NUM), indicates which cluster
!    each point belongs to.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the cluster centers.
!
!    Output, integer ( kind = 4 ) CLUSTER_POPULATION(CLUSTER_NUM), the number 
!    of points in each cluster.
!
!    Output, real ( kind = 8 ) CLUSTER_ENERGY(CLUSTER_NUM), the 
!    cluster energies.
		*/

		void kmeans_01(int *, int *, int *, int *, int *, double * __restrict, int * __restrict,
					   double * __restrict, int * __restrict, double * __restrict);

		/*
			KMEANS_02 applies the K-Means algorithm.
!
!  Discussion:
!
!    The routine attempts to divide POINT_NUM points in 
!    DIM_NUM-dimensional space into CLUSTER_NUM clusters so that the within 
!    cluster sum of squares is minimized.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2002
!
!  Author:
!
!    Original FORTRAN77 by John Hartigan, Manchek Wong.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Hartigan, Manchek Wong,
!    Algorithm AS 136:
!    A K-Means Clustering Algorithm,
!    Applied Statistics,
!    Volume 28, Number 1, 1979, pages 100-108.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Output, integer ( kind = 4 ) IT_NUM, the number of iterations taken.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the coordinates 
!    of the points.
!
!    Output, integer ( kind = 4 ) CLUSTER(POINT_NUM), the cluster each 
!    point belongs to.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the cluster centers.
!
!    Output, integer ( kind = 4 ) CLUSTER_POPULATION(CLUSTER_NUM), the number 
!    of points in each cluster.
!
!    Output, real ( kind = 8 ) CLUSTER_ENERGY(CLUSTER_NUM), the 
!    within-cluster sum of squares.
!
		*/

		void kmeans_02(int *, int *, int *, int *, int *, double * __restrict, int * __restrict,
					   double * __restrict, int * __restrict, double * __restrict );

		/*
			 KMEANS_02_OPTRA carries out the optimal transfer stage.
!
!  Discussion:
!
!    Each point is re-allocated, if necessary, to the cluster that
!    will induce a maximum reduction in the within-cluster sum of
!    squares.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 March 2002
!
!  Author:
!
!    Original FORTRAN77 by John Hartigan, Manchek Wong.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Hartigan, Manchek Wong,
!    Algorithm AS 136:
!    A K-Means Clustering Algorithm,
!    Applied Statistics,
!    Volume 28, Number 1, 1979, pages 100-108.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the coordinates of 
!    the points.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the cluster centers.
!
!    Input/output, integer ( kind = 4 ) CLUSTER(POINT_NUM), the cluster 
!    each point belongs to.
!
!    Input/output, integer ( kind = 4 ) CLUSTER2(POINT_NUM), the cluster 
!    to which each point is most likely to be transferred to.
!
!    Input/output, integer ( kind = 4 ) CLUSTER_POPULATION(CLUSTER_NUM), 
!    the number of points in each cluster.
!
!    Input/output, real ( kind = 8 ) AN1(CLUSTER_NUM), 
!    CLUSTER_POPULATION(L) / (CLUSTER_POPULATION(L) - 1)
!
!    Input/output, real ( kind = 8 ) AN2(CLUSTER_NUM), 
!    CLUSTER_POPULATION(L) / (CLUSTER_POPULATION(L) + 1)
!
!    Input/output, integer NCP(CLUSTER_NUM), ?
!
!    Input/output, real ( kind = 8 ) D(POINT_NUM), ?
!
!    Input/output, integer ( kind = 4 ) ITRAN(CLUSTER_NUM), 
!    1 if cluster L is updated in the quick-transfer stage,
!    0 otherwise.  Reset to zero on output.
!
!    Input/output, integer ( kind = 4 ) LIVE(CLUSTER_NUM), ?
!
!    Input/output, integer ( kind = 4 ) INDX, ?
!
		*/

		void kmeans_02_optra(int *, int *, int *, double * __restrict, double * __restrict, int * __restrict,
							 int * __restrict, int * __restrict, double * __restrict, double * __restrict,
							 int * __restrict, double * __restrict, int * __restrict, int * __restrict, int *);

		/*
			KMEANS_02_QTRAN carries out the quick transfer stage.
!
!  Discussion:
!
!    For each point I, CLUSTER(I) and CLUSTER2(I) are switched, if necessary, 
!    to reduce within-cluster sum of squares.  The cluster centers are
!    updated after each step.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2002
!
!  Author:
!
!    Original FORTRAN77 by John Hartigan, Manchek Wong.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Hartigan, Manchek Wong,
!    Algorithm AS 136:
!    A K-Means Clustering Algorithm,
!    Applied Statistics,
!    Volume 28, Number 1, 1979, pages 100-108.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the coordinates 
!    of the points.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the cluster centers.
!
!    Input/output, integer ( kind = 4 ) CLUSTER(POINT_NUM), the cluster 
!    each point belongs to.
!
!    Input/output, integer ( kind = 4 ) CLUSTER2(POINT_NUM), the cluster to 
!    which each point is most likely to be transferred to.
!
!    Input/output, integer ( kind = 4 ) CLUSTER_POPULATION(CLUSTER_NUM), 
!    the number of points in each cluster.
!
!    Input/output, real ( kind = 8 ) AN1(CLUSTER_NUM), 
!    CLUSTER_POPULATION(L) / (CLUSTER_POPULATION(L) - 1).
!
!    Input/output, real ( kind = 8 ) AN2(CLUSTER_NUM), 
!    CLUSTER_POPULATION(L) / (CLUSTER_POPULATION(L) + 1).
!
!    Input/output, integer NCP(CLUSTER_NUM), ?
!
!    Input/output, real ( kind = 8 ) D(POINT_NUM), ?
!
!    Input/output, integer ITRAN(CLUSTER_NUM), 
!    1 if cluster L is updated in the quick-transfer stage,
!    0 otherwise.
!
!    Input/output, integer ( kind = 4 ) INDX, is set to 0 if any 
!    updating occurs.
!
		*/

		void kmeans_02_qtran(int *, int *, int *, double * __restrict, double * __restrict, int * __restrict,
							 int * __restrict, int * __restrict, double * __restrict, double * __restrict,
							 int * __restrict, double * __restrict, int * __restrict, int * );

		/*
			 KMEANS_03 applies the K-Means algorithm.
!
!  Discussion:
!
!    It is possible for a straightforward K-Means algorithm to
!    halt at a non-optimal partition of the points.  This routine
!    tries to improve the input partition if possible.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 October 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wendy Martinez, Angel Martinez,
!    Computational Statistics Handbook with MATLAB,
!    pages 373-376,
!    Chapman and Hall / CRC, 2002.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of data points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Output, integer ( kind = 4 ) IT_NUM, the number of iterations taken.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Output, integer ( kind = 4 ) CLUSTER(POINT_NUM), the cluster to which
!    each point belongs.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), the 
!    centers associated with the clustering.  On output, these may 
!    have been altered.
!
!    Output, integer ( kind = 4 ) CLUSTER_POPULATION(CLUSTER_NUM), the number
!    of points in each cluster.
!
!    Output, real ( kind = 8 ) CLUSTER_ENERGY(CLUSTER_NUM), the energy of 
!    the clusters.
!
		*/

		void kmeans_03(int *, int *, int *, int *, int *, double * __restrict,
					   int * __restrict, double * __restrict, int * __restrict,
					   double * __restrict);

		/*
			 KMEANS_W_01 applies the weighted K-Means algorithm.
!
!  Discussion:
!
!    The input data for the weight K-Means problem includes:
!    * a set of N data points X in M dimensions, 
!    * a set of N nonnegative weights W,
!    * a desired number of clusters K.
!    * an initial set of cluster centers Z,
!    * an (optional) initial set of cluster assignments.
!
!    The goal is to determine K points Z, called cluster centers, and
!    to assign each point X(I) to some cluster Z(J), so that we minimize
!    the weighted standard deviation of the distance of each data point
!    to the center of its cluster.  Writing J = CLUSTER(I) to
!    indicate the index of the nearest cluster center Z(J) to the 
!    point X(I), the quantity we are trying to minimize is the sum
!    of the weighted cluster energies E(J), where:
!
!      E(J) = Sum ( 1 <= I <= N ) W(I) * || X(I) - Z(J) ||^2
!
!    Here, we assume that we are using the Euclidean norm, so that
!    
!      || X(I) - Z(J) ||^2 = Sum ( 1 <= K <= M )
!         ( X(I)(K) - Z(J)(K) )^2
!
!    In this notation, X(I)(K) is the K-th spatial component of the
!    I-th data point.
!
!    Note that this routine should give the same results as KMEANS_01
!    in any case in which all the entries of the WEIGHT vector are equal.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Sparks,
!    Algorithm AS 58: Euclidean Cluster Analysis,
!    Applied Statistics,
!    Volume 22, Number 1, 1973, pages 126-130.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Output, integer ( kind = 4 ) IT_NUM, the number of iterations taken.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the points.
!
!    Input, real ( kind = 8 ) WEIGHT(POINT_NUM), the weights.
!
!    Output, integer ( kind = 4 ) CLUSTER(POINT_NUM), indicates which cluster
!    each point belongs to.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the cluster centers.
!
!    Output, integer ( kind = 4 ) CLUSTER_POPULATION(CLUSTER_NUM), the number 
!    of points in each cluster.
!
!    Output, real ( kind = 8 ) CLUSTER_ENERGY(CLUSTER_NUM), the 
!    cluster energies.
		*/

		void kmeans_w_01(int *, int *, int *, int *, int *, double * __restrict,
					     double * __restrict, int * __restrict, double * __restrict,
						 int * __restrict, double * __restrict);

		/*
			KMEANS_W_03 applies the weighted K-Means algorithm.
!
!  Discussion:
!
!    The input data for the weight K-Means problem includes:
!    * a set of N data points X in M dimensions, 
!    * a set of N nonnegative weights W,
!    * a desired number of clusters K.
!    * an initial set of cluster centers Z,
!    * an (optional) initial set of cluster assignments.
!
!    The goal is to determine K points Z, called cluster centers, and
!    to assign each point X(I) to some cluster Z(J), so that we minimize
!    the weighted standard deviation of the distance of each data point
!    to the center of its cluster.  Writing J = CLUSTER(I) to
!    indicate the index of the nearest cluster center Z(J) to the 
!    point X(I), the quantity we are trying to minimize is the sum
!    of the weighted cluster energies E(J), where:
!
!      E(J) = Sum ( 1 <= I <= N ) W(I) * || X(I) - Z(J) ||**2
!
!    Here, we assume that we are using the Euclidean norm, so that
!    
!      || X(I) - Z(J) ||**2 = Sum ( 1 <= K <= M )
!        ( X(I)(K) - Z(J)(K) )**2
!
!    In this notation, X(I)(K) is the K-th spatial component of the
!    I-th data point.
!
!    Note that this routine should give the same results as KMEANS_01
!    in any case in which all the entries of the WEIGHT vector are equal.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wendy Martinez, Angel Martinez,
!    Computational Statistics Handbook with MATLAB,
!    pages 373-376,
!    Chapman and Hall / CRC, 2002.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of data points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Output, integer ( kind = 4 ) IT_NUM, the number of iterations taken.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Input, real ( kind = 8 ) WEIGHT(POINT_NUM), the weights.
!
!    Input/output, integer ( kind = 4 ) CLUSTER(POINT_NUM), the cluster 
!    to which each point belongs.  On output, these may have been altered.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), the
!    centers associated with the clustering.  On output, these may
!    have been altered.
!
!    Output, integer ( kind = 4 ) CLUSTER_POPULATION(CLUSTER_NUM), the number
!    of points in each cluster.
!
!    Output, real ( kind = 8 ) CLUSTER_ENERGY(CLUSTER_NUM), the energy of
!    the clusters.
		*/

		void kmeans_w_03(int *, int *, int *, int *, int *, double * __restrict,
		                 double * __restrict, int * __restrict, double * __restrict,
						 int * __restrict, double * __restrict);

		/*
			 R4_UNIFORM_01 returns a unit pseudorandom R4.
!
!  Discussion:
!
!    An R4 is a real ( kind = 4 ) value.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r4_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R4_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, 
!    which should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 4 ) R4_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
		*/

		float r4_uniform_01(int *);

		/*
			 R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r8_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
		*/

		double r8_uniform_01(int *);

		/*
			 RANDOM_INITIALIZE initializes the FORTRAN90 random number seed.
!
!  Discussion:
!
!    If you don't initialize the random number generator, its behavior
!    is not specified.  If you initialize it simply by:
!
!      call random_seed ( )
!
!    its behavior is not specified.  On the DEC ALPHA, if that's all you
!    do, the same random number sequence is returned.  In order to actually
!    try to scramble up the random number generator a bit, this routine
!    goes through the tedious process of getting the size of the random
!    number seed, making up values based on the current time, and setting
!    the random number seed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED.
!    If SEED is zero on input, then you're asking this routine to come up
!    with a seed value, which is returned as output.
!    If SEED is nonzero on input, then you're asking this routine to
!    use the input value of SEED to initialize the random number generator,
!    and SEED is not changed on output.
!
		*/

		void random_initialize(int *);

		/*
			R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of real ( kind = 8 ) values.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in 
!    the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values
		*/

		void r8vec_uniform_01(int *, int *, double * __restrict);

		/*
			R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of real ( kind = 8 ) values.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the array.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(M,N), the array of pseudorandom values.
		*/
		void r8mat_uniform_01(int *, int *, int *, double * __restrict);



} // extern "C"

#endif

#endif /*__GMS_KMEANS_WRAPPER_H__*/
