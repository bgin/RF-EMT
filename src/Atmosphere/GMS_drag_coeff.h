
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

#ifndef __GMS_DRAG_COEFF_H__
#define __GMS_DRAG_COEFF_H__



	/*
			 This code is based on Fortran 2003 version developed by Bernard Gingold
		     References:
					
							  "Fall Velocities of Hydrometeors in the Atmosphere:
                              Refinements to a Continuous Analytical Power Law "
 
                              VITALY I. KHVOROSTYANOV
                                                  Central Aerological Observatory, Dolgoprudny, Moscow, Russia
                              JUDITH A. CURRY
                                                  School of Earth and Atmospheric Sciences, Georgia Institute of Technology, Atlanta, Georgia 
	*/

#include <cstdint>

namespace gms {
	

		/*
				Computes drag coefficients per single body (hydrometeor)
		        Cd = C0(1 + delta0/Re^0.5)^2
		*/
		void compute_drag_coeff(_Out_  double * __restrict ,
							    _In_   const int32_t,
								_In_   const double * __restrict,
								_In_   const double,
								_In_   const double * __restrict );


		/*
				Calcualtion of Reynolds number (Re)
		*/
		void compute_reynolds_number(_Out_ double * __restrict, 
									 _In_  const int32_t,
									 _In_  const double * __restrict,
									 _In_  const double,
									 _In_  const double * __restrict,
									 _In_  const double * __restrict,
									 _In_  const int32_t,
									 _In_  const int32_t,
									 _In_  const int32_t );

		/*
				  Calculation of Davies (Best) number X
		*/
		void compute_davies_number(_Out_ double * __restrict,
								   _In_  const int32_t,
								   _In_  const double,
								   _In_  const double,
								   _In_  const double * __restrict,
								   _In_  const int32_t,
								   _In_  const int32_t,
								   _In_  const int32_t,
								   _In_  const double,
								   _In_  const double,
								   _In_  const double * __restrict );

		/*
				Calculation of analytical function of X known as a 'aRe'	
		*/
		void compute_aRe(_Out_ double * __restrict ,
						 _In_  double * __restrict,
						 _In_   double * __restrict,
						 _In_  const int32_t,
						 _In_  const double,
						 _In_  const double * __restrict);

		/*
				Calculation of analytical function of X known as a 'bRe'
		*/
		void compute_bRe(_Out_ double * __restrict,
						 _In_  const int32_t,
						 _In_  const double * __restrict,
						 _In_  const double * __restrict);

		/*
			      Calculation of turbulence correction for Re ~ 10.0^3.0
			      Psi - interpolation function implementated as a subroutine.
		*/
		void psi_interpolation(_Out_ double *       __restrict,
							   _In_  const double * __restrict,
							   _In_  const double * __restrict,
							   _In_  const int32_t,
							   _In_  const double,
							   _In_  const double * __restrict);

		/*
				  Calculation of bRe - turbulent correction of (delta bRe)
		*/
		void compute_delta_bRe(_Out_ double * __restrict,
							   _In_  const int32_t,
							   _In_  const double * __restrict,
							   _In_  const double * __restrict,
							   _In_  const double * __restrict,
							   _In_  const double );
	
}

#endif /*__GMS_DRAG_COEFF_H__*/
