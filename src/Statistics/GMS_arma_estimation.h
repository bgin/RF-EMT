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

#ifndef __GMS_ARMA_ESTIMATION_H__
#define __GMS_ARMA_ESTIMATION_H__




#include "GMS_config.h"


/*
      Converted to C++ from F90 implementation of the following ASA algorithms
      ALGORITHM AS 154  APPL. STATIST. (1980) VOL.29, P.311
      ALGORITHM AS 182  APPL. STATIST. (1982) VOL.31, NO.2
      Modified by Bernard Gingold on 29-11-2020 SUN 2:53PM +00200
*/

/*
   !  ALGORITHM AS 182  APPL. STATIST. (1982) VOL.31, NO.2

!  Finite sample prediction from ARIMA processes.
*/

namespace gms {

        namespace stat {

__ATTR_HOT__
__ATTR_ALIGN__(32)
void forkal(const int32_t ip,
            const int32_t iq,
            const int32_t ir,
            const int32_t np,
            const int32_t ird,
            const int32_t irz,
            const int32_t id,
            const int32_t il,
            const int32_t n,
            const int32_t nrbar,
            float * __restrict __attribute__((aligned(64))) phi,
            float * __restrict __attribute__((aligned(64))) theta,
            float * __restrict __attribute__((aligned(64))) delta,
            float * __restrict __attribute__((aligned(64))) w,
            float * __restrict __attribute__((aligned(64))) y,
            float * __restrict __attribute__((aligned(64))) amse,
            float * __restrict __attribute__((aligned(64))) a,
            float * __restrict __attribute__((aligned(64))) p,
            float * __restrict __attribute__((aligned(64))) v,
            float * __restrict __attribute__((aligned(64))) resid,
            float * __restrict __attribute__((aligned(64))) e,
            float * __restrict __attribute__((aligned(64))) xnext,
            float * __restrict __attribute__((aligned(64))) xrow,
            float * __restrict __attribute__((aligned(64))) rbar,
            float * __restrict __attribute__((aligned(64))) thetab,
            float * __restrict __attribute__((aligned(64))) store,
            int32_t &ifault ); 
/*
    !  INVOKING THIS SUBROUTINE SETS THE VALUES OF V AND PHI, AND
!  OBTAINS THE INITIAL VALUES OF A AND P.
!  THIS ROUTINE IS NOT SUITABLE FOR USE WITH AN AR(1) PROCESS.
!  IN THIS CASE THE FOLLOWING INSTRUCTIONS SHOULD BE USED FOR INITIALISATION.
!  V(1) = 1.0
!  A(1) = 0.0
!  P(1) = 1.0 / (1.0 - PHI(1) * PHI(1))
*/

__ATTR_HOT__
__ATTR_ALIGN__(32)
void starma(const int32_t ip,
            const int32_t iq,
            const int32_t ir,
            const int32_t np,
            float * __restrict __attribute__((aligned(64))) phi,
            float * __restrict __attribute__((aligned(64))) theta,
            float * __restrict __attribute__((aligned(64))) a,
            float * __restrict __attribute__((aligned(64))) p,
            float * __restrict __attribute__((aligned(64))) v,
            float * __restrict __attribute__((aligned(64))) thetab,
            float * __restrict __attribute__((aligned(64))) xnext,
            float * __restrict __attribute__((aligned(64))) xrow,
            float * __restrict __attribute__((aligned(64))) rbar,
            const int32_t nrbar,
            int32_t &ifault); 
/*
   ! N.B. Argument NRBAR has been removed.

!   ALGORITHM AS 154.3  APPL. STATIST. (1980) VOL.29, P.311

!   FORTRAN VERSION OF REVISED VERSION OF ALGORITHM AS 75.1
!   APPL. STATIST. (1974) VOL.23, P.448
!   SEE REMARK AS R17 APPL. STATIST. (1976) VOL.25, P.323
*/

__ATTR_HOT__
__ATTR_ALIGN__(32)
void inclu2(const int32_t np,
            const float weight,
            float * __restrict __attribute__((aligned(64))) xnext,
            float * __restrict __attribute__((aligned(64))) xrow,
            const float ynext,
            float * __restrict __attribute__((aligned(64))) d,
            float * __restrict __attribute__((aligned(64))) rbar,
            float * __restrict __attribute__((aligned(64))) thetab,
            float &ssqerr,
            float &recres,
            int32_t &irank,
            int32_t &ifault); 
/*
   !   ALGORITHM AS 154.4  APPL. STATIST. (1980) VOL.29, P.311

!   REVISED VERSION OF ALGORITHM AS 75.4
!   APPL. STATIST. (1974) VOL.23, P.448
!   INVOKING THIS SUBROUTINE OBTAINS BETA BY BACKSUBSTITUTION
!   IN THE TRIANGULAR SYSTEM RBAR AND THETAB.
*/

__ATTR_HOT__
__ATTR_ALIGN__(32)
void regres(const int32_t np,
            const int32_t nrbar,
            float * __restrict __attribute__((aligned(64))) rbar,
            float * __restrict __attribute__((aligned(64))) thetab,
            float * __restrict __attribute__((aligned(64))) beta); 

__ATTR_HOT__
__ATTR_ALIGN__(32)
void karma(const int32_t ip,
           const int32_t iq,
	   const int32_t ir,
	   float * __restrict __attribute__((aligned(64))) phi,
	   float * __restrict __attribute__((aligned(64))) theta,
	   float * __restrict __attribute__((aligned(64))) a,
	   float * __restrict __attribute__((aligned(64))) p,
	   float * __restrict __attribute__((aligned(64))) v,
	   const int32_t n,
	   float * __restrict __attribute__((aligned(64))) w,
	   float * __restrict __attribute__((aligned(64))) resid,
	   float &sumlog,
	   float &ssq,
	   int32_t iupd,
	   const float delta,
	   float * __restrict __attribute__((aligned(64))) e,
	   int32_t &nit); 




} // stat

} // gms






#endif /*__ARMA_ESTIMATION_H__*/
