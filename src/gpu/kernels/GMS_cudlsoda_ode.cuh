
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

#ifndef __GMS_CUDLSODA_ODE_CUH__
#define __GMS_CUDLSODA_ODE_CUH__

#include <stdint.h>

/*

     List of test ODE for CUDA lsoda integrator.
*/


struct ODE_1 {

    __device__ void operator()(int32_t *  __restrict__ neq,
                               double  *  __restrict__ t,
                               double  *  __restrict__ y,
                               double  *  __restrict__ ydot) {
        
            ydot[0] = 1.0E4*y[1]*y[2]-0.04E0*y[0];
            ydot[2] = 3.0E7*y[1]*y[1];
            ydot[1] = 1.0*(ydot[0]+ydot[2]);
    }                 
};

/*

! DEMONSTRATION PROGRAM FOR THE DVODE_F90 PACKAGE.

! The following is a simple example problem, with the coding
! needed for its solution by DVODE_F90. The problem is from
! chemical kinetics, and consists of the following three rate
! equations:
!     dy1/dt = -.04d0*y1 + 1.d4*y2*y3
!     dy2/dt = .04d0*y1 - 1.d4*y2*y3 - 3.d7*y2**2
!     dy3/dt = 3.d7*y2**2
! on the interval from t = 0.0d0 to t = 4.d10, with initial
! conditions y1 = 1.0d0, y2 = y3 = 0.0d0. The problem is stiff.
*/
struct ODE_2 {
 
   __device__ void operator()( int32_t *  __restrict__ neq,
                               double  *  __restrict__ t,
                               double  *  __restrict__ y,
                               double  *  __restrict__ ydot) {
    
          ydot[0] = -.04E0*y[0] + 1.0E4*y[1]*y[2];
          ydot[2] = 3.E7*y[1]*y[1];
          ydot[1] = -ydot[0]-ydot[2];
   } 

};


struct JAC_2 {

    __device__ void operator()(int32_t * __restrict__ neq,
                               double  * __restrict__ t,
                               double  * __restrict__ y,
                               int32_t * __restrict__ ml,
                               int32_t * __restrict__ mu,
                               double  * __restrict__ pd,
                               int32_t * __restrict__ nrpd) {
         pd[0] = -0.04;
         pd[1] = 1.0E4*y[2];
         pd[2] = 1.0E4*y[1];
         pd[3] = 0.04;
         pd[5] = -pd[2];
         pd[7] = 6.0E7*y[1];
         pd[4] = -pd[1]-pd[7];
    }
};


/*
       SUBROUTINE F2(NEQ,T,Y,YDOT)
        IMPLICIT NONE
        INTEGER NEQ
        DOUBLE PRECISION T, Y, YDOT
        DIMENSION Y(2), YDOT(2)

        YDOT(1) = Y(2)
        YDOT(2) = 100.0D0*(1.0D0-Y(1)*Y(1))*Y(2) - Y(1)
        RETURN
      END SUBROUTINE F2

      SUBROUTINE JAC2(NEQ,T,Y,ML,MU,PD,NROWPD)
        IMPLICIT NONE
        INTEGER NEQ, ML, MU, NROWPD
        DOUBLE PRECISION T, Y, PD
        DIMENSION Y(2), PD(NROWPD,2)

        PD(1,1) = 0.0D0
        PD(1,2) = 1.0D0
        PD(2,1) = -200.0D0*Y(1)*Y(2) - 1.0D0
        PD(2,2) = 100.0D0*(1.0D0-Y(1)*Y(1))
        RETURN
      END SUBROUTINE JAC2
*/


struct ODE_3 {

     __device__ void operator()(int32_t *  __restrict__ neq,
                                double  *  __restrict__ t,
                                double  *  __restrict__ y,
                                double  *  __restrict__ ydot) {
         
           ydot[0] = ydot[1];
           ydot[1] = 100.0E0*(1.0E0-y[0]*y[0])*y[1]-y[0];

     }
};


struct JAC_3 {

      __device__ void operator()(int32_t * __restrict__ neq,
                                 double  * __restrict__ t,
                                 double  * __restrict__ y,
                                 int32_t * __restrict__ ml,
                                 int32_t * __restrict__ mu,
                                 double  * __restrict__ pd,
                                 int32_t * __restrict__ nrpd) {
          
            pd[0] = 0.0E0;
            pd[1] = 1.0E0;
            pd[2] = -200.0E0*y[0]*y[1]-1.0E0;
            pd[3] = 100.0E0*(1.0E0-y[0]*y[1]);
      }
};


/*
      SUBROUTINE F1(NEQ,T,Y,YDOT)
        IMPLICIT NONE
        INTEGER NEQ
        DOUBLE PRECISION T, Y, YDOT
        DIMENSION Y(1), YDOT(1)

        YDOT(1) = ((2.0D0*LOG(Y(1))+8.0D0)/T-5.0D0)*Y(1)
        RETURN
      END SUBROUTINE F1
*/


struct ODE_3 {

      __device__ void operator()(int32_t *  __restrict__ neq,
                                 double  *  __restrict__ t,
                                 double  *  __restrict__ y,
                                 double  *  __restrict__ ydot) {
           
           ydot[0] = ((2.0E0*__logf(y[0])+8.0E0)/t-5.0E0)*y[0];
      }
};


/*
      SUBROUTINE F1(NEQ,T,Y,YDOT)
        IMPLICIT NONE
        INTEGER NEQ
        DOUBLE PRECISION T, Y, YDOT
        DIMENSION Y(2), YDOT(2)

        YDOT(1) = Y(2)
        YDOT(2) = 3.0D0*(1.0D0-Y(1)*Y(1))*Y(2) - Y(1)
        RETURN
      END SUBROUTINE F1

      SUBROUTINE JAC1(NEQ,T,Y,ML,MU,PD,NROWPD)
        IMPLICIT NONE
        INTEGER NEQ, ML, MU, NROWPD
        DOUBLE PRECISION T, Y, PD
        DIMENSION Y(2), PD(NROWPD,2)

        PD(1,1) = 0.0D0
        PD(1,2) = 1.0D0
        PD(2,1) = -6.0D0*Y(1)*Y(2) - 1.0D0
        PD(2,2) = 3.0D0*(1.0D0-Y(1)*Y(1))
        RETURN
      END SUBROUTINE JAC1
*/


struct ODE_4 {

      __device__ void operator()(int32_t *  __restrict__ neq,
                                 double  *  __restrict__ t,
                                 double  *  __restrict__ y,
                                 double  *  __restrict__ ydot) {
          
           ydot[0] = y[1];
           ydot[1] = 3.0E0*(1.0E0-y[0]*y[0])*y[1]-y[0];
     }
};


struct JAC_4 {

      __device__ void operator()(int32_t * __restrict__ neq,
                                 double  * __restrict__ t,
                                 double  * __restrict__ y,
                                 int32_t * __restrict__ ml,
                                 int32_t * __restrict__ mu,
                                 double  * __restrict__ pd,
                                 int32_t * __restrict__ nrpd) {
          
           pd[0] = 0.0E0;
           pd[1] = 1.0E0;
           pd[2] = -6.0E0*y[0]*y[1]-1.0E0;
           pd[3] = 3.0E0*(1.0E0-y[0]*y[0]);
     }
};


/*
   
     SUBROUTINE DERIVS(NEQ,T,Y,YDOT)
!     Subroutine to evaluate dy/dt for this problem
        IMPLICIT NONE
        INTEGER NEQ
        DOUBLE PRECISION T, Y, YDOT
        DIMENSION Y(NEQ), YDOT(NEQ)
        DOUBLE PRECISION PROD1, PROD2, PROD3, PROD4

        PROD1 = 7.89D-10*Y(1)
        PROD2 = 1.1D7*Y(1)*Y(3)
        PROD3 = 1.13D9*Y(2)*Y(3)
        PROD4 = 1.13D3*Y(4)
        YDOT(1) = -PROD1 - PROD2
        YDOT(2) = PROD1 - PROD3
        YDOT(4) = PROD2 - PROD4
        YDOT(3) = YDOT(2) - YDOT(4)
        RETURN
      END SUBROUTINE DERIVS

      SUBROUTINE JACD(NEQ,T,Y,ML,MU,PD,NROWPD)
!     Subroutine to define the exact Jacobian for this problem
        IMPLICIT NONE
        INTEGER NEQ, ML, MU, NROWPD
        DOUBLE PRECISION T, Y, PD
        DIMENSION Y(NEQ), PD(NROWPD,NEQ)
        DOUBLE PRECISION A, B, CM, C

        A = 7.89D-10
        B = 1.1D7
        CM = 1.13D9
        C = 1.13D3
        PD(1,1) = -A - B*Y(3)
        PD(1,2) = 0.D0
        PD(1,3) = -B*Y(1)
        PD(1,4) = 0.D0
        PD(2,1) = A
        PD(2,2) = -CM*Y(3)
        PD(2,3) = -CM*Y(2)
        PD(2,4) = 0.D0
        PD(3,1) = A - B*Y(3)
        PD(3,2) = -CM*Y(3)
        PD(3,3) = -B*Y(1) - CM*Y(2)
        PD(3,4) = C
        PD(4,1) = B*Y(3)
        PD(4,2) = 0.D0
        PD(4,3) = B*Y(1)
        PD(4,4) = -C
        RETURN
      END SUBROUTINE JACD
*/


struct ODE_5 {

     __device__ void operator()( int32_t *  __restrict__ neq,
                                 double  *  __restrict__ t,
                                 double  *  __restrict__ y,
                                 double  *  __restrict__ ydot) {
      
            double prod1,prod2,prod3,prod4;
            prod1 = 7.89E0-10.0E0*y[0];
            prod2 = 1.1E7*y[0]*y[2];
            prod3 = 1.13E9*y[1]*y[2];
            prod4 = 1.13E3*y[3];
            ydot[0] = -prod1 - prod2;
            ydot[1] = prod1 - prod3;
            ydot[3] = prod2 - prod4;
            ydot[2] = ydot[1] - ydot[3];
     }
          
};


struct JAC_5 {

      __device__ void operator()(int32_t * __restrict__ neq,
                                 double  * __restrict__ t,
                                 double  * __restrict__ y,
                                 int32_t * __restrict__ ml,
                                 int32_t * __restrict__ mu,
                                 double  * __restrict__ pd,
                                 int32_t * __restrict__ nrpd) {
             
                double A,B,CM,C;
                A = 7.89E-10;
                B = 1.1E+7;
                CM= 1.13E+9;
                C = 1.13E+3;
                pd[0] = -A - B*y[2];
                pd[1] = 0.0E+0;
                pd[2] = -B*y[0];
                pd[3] = 0.0E+0;
                pd[4] = A;
                pd[5] = -CM*y[2];
                pd[6] = -CM*y[1];
                pd[7] = 0.0E+0;
                pd[8] = A-B*y[2];
                pd[9] = -CM*y[2];
                pd[10]= -B*y[0]-CM*y[1];
                pd[11]= C;
                pd[12]= B*y[2];
                pd[13]= 0.0E+0;
                pd[14]= B*y[0];
                pd[15]= -C;
       }
};


/*

    SUBROUTINE DERIVS(NEQ,T,Y,YDOT)
!     Subroutine to evaluate dy/dt for this problem
        IMPLICIT NONE
        INTEGER NEQ
        DOUBLE PRECISION T, Y, YDOT, A, B
        DIMENSION Y(NEQ), YDOT(NEQ)

        A = 3.12121212D0
        B = 2.11111111D0
        YDOT(1) = A*Y(3)
        YDOT(2) = B*Y(4)
        YDOT(3) = -A*Y(1)
        YDOT(4) = -B*Y(2)
        RETURN
      END SUBROUTINE DERIVS
     
*/

struct ODE_6 {

      __device__ void operator()(int32_t *  __restrict__ neq,
                                 double  *  __restrict__ t,
                                 double  *  __restrict__ y,
                                 double  *  __restrict__ ydot) {
            
           const double A =  3.12121212;
           const double B =  2.11111111;
           ydot[0]        = A*y[2];
           ydot[1]        = B*y[3];
           ydot[2]        = -A*y[0];
           ydot[3]        = -B*y[1];  
      }
};


/*
     SUBROUTINE DERIVS(NEQ,T,Y,YDOT)
!     Subroutine to evaluate dy/dt for this problem
        IMPLICIT NONE
        INTEGER NEQ
        DOUBLE PRECISION T, Y, YDOT
        DIMENSION Y(NEQ), YDOT(NEQ)

        YDOT(1) = -1.71D0*Y(1) + 0.43D0*Y(2) + 8.32D0*Y(3) + 0.0007D0
        YDOT(2) = 1.71D0*Y(1) - 8.75D0*Y(2)
        YDOT(3) = -10.03D0*Y(3) + 0.43D0*Y(4) + 0.035D0*Y(5)
        YDOT(4) = 8.32D0*Y(2) + 1.71D0*Y(3) - 1.12D0*Y(4)
        YDOT(5) = -1.745D0*Y(5) + 0.43D0*(Y(6)+Y(7))
        YDOT(6) = -280D0*Y(6)*Y(8) + 0.69D0*Y(4) + 1.71D0*Y(5) - 0.43D0*Y(6) + &
          0.69D0*Y(7)
        YDOT(7) = 280D0*Y(6)*Y(8) - 1.81D0*Y(7)
        YDOT(8) = -YDOT(7)
        RETURN
      END SUBROUTINE DERIVS

      SUBROUTINE JACD(NEQ,T,Y,ML,MU,PD,NROWPD)
!     Subroutine to define the dense Jacobian for this problem
        IMPLICIT NONE
        INTEGER NEQ, ML, MU, NROWPD
        DOUBLE PRECISION T, Y, PD
        DIMENSION Y(NEQ), PD(NROWPD,NEQ)

        PD(1:NEQ,1:NEQ) = 0D0
        PD(1,1) = -1.71D0
        PD(1,2) = 0.43D0
        PD(1,3) = 8.32D0
        PD(2,1) = 1.71D0
        PD(2,2) = -8.75D0
        PD(3,3) = -10.03D0
        PD(3,4) = 0.43D0
        PD(3,5) = 0.035D0
        PD(4,2) = 8.32D0
        PD(4,3) = 1.71D0
        PD(4,4) = -1.12D0
        PD(5,5) = -1.745D0
        PD(5,6) = 0.43D0
        PD(5,7) = 0.43D0
        PD(6,4) = 0.69D0
        PD(6,5) = 1.71D0
        PD(6,6) = -280D0*Y(8) - 0.43D0
        PD(6,7) = 0.69D0
        PD(6,8) = -280D0*Y(6)
        PD(7,6) = 280D0*Y(8)
        PD(7,7) = -1.81D0
        PD(7,8) = 280D0*Y(6)
        PD(8,6) = -280D0*Y(8)
        PD(8,7) = 1.81D0
        PD(8,8) = -280D0*Y(6)
        RETURN
      END SUBROUTINE JACD

*/


struct ODE_7 {

      __device__ void operator()(int32_t *  __restrict__ neq,
                                 double  *  __restrict__ t,
                                 double  *  __restrict__ y,
                                 double  *  __restrict__ ydot) { {
         
             ydot[0] = 1.71E0*y[0]+0.43E0*y[1]+8.32E0*y[2]+0.0007E0;
             ydot[1] = 1.71E0*y[0]-8.75E0*y[1];
             ydot[2] = -10.03E0*y[2]+0.43E0*y[3]+0.035E0*y[4];
             ydot[3] = 8.32E0*y[1]+1.71E0*y[2]-1.12E0*y[3];
             ydot[4] = -1.745E0*y[4] + 0.43E0*(y[5]+y[6]);
             ydot[5] = -280E0*y[5]*y[7]+0.69E0*y[3]+1.71E0*y[4]-0.43E0*y[6]+0.69E0*y[6];
             ydot[6] = 280.0E0*y[5]*y[7]-1.81E0*y[6];
             ydot[7] = -ydot[6];
      }
};


struct JAC_7 {

      __device__ void operator()(int32_t * __restrict__ neq,
                                 double  * __restrict__ t,
                                 double  * __restrict__ y,
                                 int32_t * __restrict__ ml,
                                 int32_t * __restrict__ mu,
                                 double  * __restrict__ pd,
                                 int32_t * __restrict__ nrpd) {
          
            for(int32_t i = 0; i != neq; ++i)
                for(int32_t j = 0; j != neq; ++j)
                    pd[i+neq*j] = 0.0E+0;

            pd[0] = -1.71E+0;
            pd[1] = 0.43E+0;
            pd[2] = 8.32E+0;
            pd[8] = 1.71E+0;
            pd[9] = -8.75E+0;
            pd[18]= -10.03E+0;
            pd[19]= 0.43E+0;
            pd[20]= 0.035E+0;
            pd[24]= 8.32E+0;
            pd[25]= 7.71E+0;
            pd[26]= -1.12E+0;
            pd[36]= -1.745E+0;
            pd[37]= 0.43E+0;
            pd[38]= 0.43E+0;
            pd[43]= 0.69E+0;
            pd[44]= 1.71E+0;
            pd[45]= -280.0E+0*y[7]-0.43E+0;
            pd[46]= 0.69E+0;
            pd[47]= -280.0E+0*y[6];
            pd[53]= 280.0E+0*y[7];
            pd[54]= 1.81E+0;
            pd[55]= 280.0E+0*y[5];
            pd[61]= -280.0E+0*y[7];
            pd[62]= 1.81E+0;
            pd[63]= -280.0E+0*y[5]; 
     }
};

/*
         SUBROUTINE F1(NEQ,T,Y,YDOT)
        IMPLICIT NONE
        INTEGER NEQ
        DOUBLE PRECISION T, Y, YDOT
        DIMENSION Y(2), YDOT(2)

        YDOT(1) = Y(2)
        YDOT(2) = 3.0D0*(1.0D0-Y(1)*Y(1))*Y(2) - Y(1)
        RETURN
      END SUBROUTINE F1

      SUBROUTINE JAC1(NEQ,T,Y,ML,MU,PD,NROWPD)
        IMPLICIT NONE
        INTEGER NEQ, ML, MU, NROWPD
        DOUBLE PRECISION T, Y, PD
        DIMENSION Y(2), PD(NROWPD,2)

        PD(1,1) = 0.0D0
        PD(1,2) = 1.0D0
        PD(2,1) = -6.0D0*Y(1)*Y(2) - 1.0D0
        PD(2,2) = 3.0D0*(1.0D0-Y(1)*Y(1))
        RETURN
      END SUBROUTINE JAC1
*/


struct ODE_8 {

     __device__ void operator()(int32_t *  __restrict__ neq,
                                 double  *  __restrict__ t,
                                 double  *  __restrict__ y,
                                 double  *  __restrict__ ydot) {

             ydot[0] = y[1];
             ydot[1] = 3.0E+0*(1.0E+0-y[0]*y[0])*y[1]-y[0];
      }
};


struct JAC_8 {

     __device__ void operator()(int32_t * __restrict__ neq,
                                 double  * __restrict__ t,
                                 double  * __restrict__ y,
                                 int32_t * __restrict__ ml,
                                 int32_t * __restrict__ mu,
                                 double  * __restrict__ pd,
                                 int32_t * __restrict__ nrpd) {
          
            pd[0] = 0.0E+0;
            pd[1] = 1.0E+0;
            pd[2] = -6.0E+0*y[0]*y[1]-1.0E+0;
            pd[3] = 3.0E+0*(1.0E+0-y[0]*y[0]);
      }
};


/*
       SUBROUTINE DERIVS(NEQ,T,Y,DY)
!     Subroutine to evaluate dy/dt for this problem
        IMPLICIT NONE
        INTEGER NEQ
        DOUBLE PRECISION T, Y, DY
        DIMENSION Y(NEQ), DY(NEQ)

        DY(1) = -0.04D0*Y(1) + 1D4*Y(2)*Y(3)
        DY(2) = 0.04D0*Y(1) - 1D4*Y(2)*Y(3) - 3D7*Y(2)**2
        DY(3) = 3D7*Y(2)**2
        RETURN
      END SUBROUTINE DERIVS

      SUBROUTINE JACD(NEQ,T,Y,ML,MU,PD,NROWPD)
!     Subroutine to define the exact Jacobian for this problem
        IMPLICIT NONE
        INTEGER NEQ, ML, MU, NROWPD
        DOUBLE PRECISION T, Y, PD
        DIMENSION Y(NEQ), PD(NROWPD,NEQ)

        PD(1:NEQ,1:NEQ) = 0D0
        PD(1,1) = -0.04D0
        PD(1,2) = 1D4*Y(3)
        PD(1,3) = 1D4*Y(2)
        PD(2,1) = 0.04D0
        PD(2,2) = -1D4*Y(3) - 6D7*Y(2)
        PD(2,3) = -1D4*Y(2)
        PD(3,1) = 0D0
        PD(3,2) = 6D7*Y(2)
        PD(3,3) = 0D0
        RETURN
      END SUBROUTINE JACD
*/


struct ODE_9 {

     __device__ void operator()( int32_t *  __restrict__ neq,
                                 double  *  __restrict__ t,
                                 double  *  __restrict__ y,
                                 double  *  __restrict__ ydot) {
         
           ydot[0] = -0.04E0*y[0]+1.0E+4*y[1]*y[2];
           ydot[1] = 0.04E0*y[0]-1.0E+4*y[1]*y[2]-3.0E+7*(y[1]*y[1]);
           ydot[2] = 3.0E+7*(y[1]*y[1]);
    }
};


struct JAC_9 {

       __device__ void operator()(int32_t * __restrict__ neq,
                                 double  * __restrict__ t,
                                 double  * __restrict__ y,
                                 int32_t * __restrict__ ml,
                                 int32_t * __restrict__ mu,
                                 double  * __restrict__ pd,
                                 int32_t * __restrict__ nrpd) {
            
            pd[0] = -0.04E0;
            pd[1] = 1.0E+4*y[2];
            pd[2] = 1.0E+4*y[1];
            pd[3] = 0.04E0;
            pd[4] = -1.0E+4*y[2]-6.0E+7*y[1];
            pd[5] = -1.0E+4*y[1];
            pd[6] = 0.0E0;
            pd[7] = 6.0E+7*y[1];
            pd[8] = 0.0E+0;
      }
};


/*

      SUBROUTINE DERIVS(NEQ,T,Y,YDOT)
!     Subroutine to evaluate dy/dt for this problem
        IMPLICIT NONE
        INTEGER NEQ
        DOUBLE PRECISION T, Y, YDOT
        DIMENSION Y(NEQ), YDOT(NEQ)

        YDOT(1) = 77.27D0*(Y(2)+Y(1)*(1.D0-8.375D-6*Y(1)-Y(2)))
        YDOT(2) = (Y(3)-(1.D0+Y(1))*Y(2))/77.27D0
        YDOT(3) = 0.161D0*(Y(1)-Y(3))
        RETURN
      END SUBROUTINE DERIVS

      SUBROUTINE JACD(NEQ,T,Y,ML,MU,PD,NROWPD)
!     Subroutine to define the exact Jacobian for this problem
        IMPLICIT NONE
        INTEGER NEQ, ML, MU, NROWPD
        DOUBLE PRECISION T, Y, PD
        DIMENSION Y(NEQ), PD(NROWPD,NEQ)

        PD(1:NEQ,1:NEQ) = 0D0
        PD(1,1) = 77.27D0*(1.D0-2.D0*8.375D-6*Y(1)-Y(2))
        PD(1,2) = 77.27D0*(1.D0-Y(1))
        PD(1,3) = 0.D0
        PD(2,1) = -Y(2)/77.27D0
        PD(2,2) = -(1.D0+Y(1))/77.27D0
        PD(2,3) = 1.D0/77.27D0
        PD(3,1) = .161D0
        PD(3,2) = .0D0
        PD(3,3) = -.161D0
        RETURN
      END SUBROUTINE JACD
   
*/


struct ODE_10 {

       __device__ void operator()(int32_t *  __restrict__ neq,
                                 double  *  __restrict__ t,
                                 double  *  __restrict__ y,
                                 double  *  __restrict__ ydot) {
            
             ydot[0] = 77.27E+0*(y[1]+y[0]*(1.0E0-8.375E-6*y[0]-y[1]);
             ydot[1] = (y[2]-(1.0E+0+y[0])*y[1])/77.27E+0;
             ydot[2] = 0.161E+0*(y[0]-y[2]);
       }
};


struct JAC_10 {

        __device__ void operator()(int32_t * __restrict__ neq,
                                 double  * __restrict__ t,
                                 double  * __restrict__ y,
                                 int32_t * __restrict__ ml,
                                 int32_t * __restrict__ mu,
                                 double  * __restrict__ pd,
                                 int32_t * __restrict__ nrpd) {
           
             pd[0] = 77.27E0*(1.0-2.0*8.375E-6*y[0]-y[1]);
             pd[1] = 77.27E0*(1.0-y[0]);
             pd[2] = 0.0;
             pd[3] = -y[1]/77.27;
             pd[4] = -(1.0+y[0])/77.27;
             pd[5] = 1.0/77.27;
             pd[6] = 0.161;
             pd[7] = 0.0;
             pd[8] = -0.161;
       }
};


/*
       SUBROUTINE DERIVS(NEQ,T,YSOL,YDOT)
        IMPLICIT NONE
        INTEGER NEQ
        DOUBLE PRECISION T, Y, YDOT, YSOL
        DIMENSION YSOL(NEQ), YDOT(NEQ)

        Y = YSOL(1)
        YDOT(1) = (1.0D0/EPSILN)*((1.0D0-T)*Y-Y*Y)
        RETURN
      END SUBROUTINE DERIVS

      SUBROUTINE JACD(NEQ,T,Y,ML,MU,PD,NROWPD)
        IMPLICIT NONE
        INTEGER NEQ, ML, MU, NROWPD
        DOUBLE PRECISION T, Y, PD
        DIMENSION Y(NEQ), PD(NROWPD,NEQ)

        PD(1,1) = (1.0D0/EPSILN)*(1.0D0-T-2.0D0*Y(1))
        RETURN
      END SUBROUTINE JACD
*/


struct ODE_11 {

      __device__ void operator()(int32_t *  __restrict__ neq,
                                 double  *  __restrict__ t,
                                 double  *  __restrict__ y,
                                 double  *  __restrict__ ydot) {
          
          const double epsilon = 1.0E-6;
          double tmp             = y[0];
          ydot[0] = (1.0/epsilon)*((1.0-t)*tmp-tmp*tmp);
     }
};


struct JAC_11 {

      __device__ void operator()(int32_t * __restrict__ neq,
                                 double  * __restrict__ t,
                                 double  * __restrict__ y,
                                 int32_t * __restrict__ ml,
                                 int32_t * __restrict__ mu,
                                 double  * __restrict__ pd,
                                 int32_t * __restrict__ nrpd) {
           
           const double epsilon = 1.0E-6;
           pd[0] = (1.0/epsilon)*(1.0-t-2.0*y[0]);
     }
};
#endif /*__GMS_CUDLSODA_ODE_CUH__*/
