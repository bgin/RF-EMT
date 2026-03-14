! ************************************************************************
      subroutine rincon_m (rout,  cnb,  rin,n);
!
!     computes the inverse of the upper triangular vector stored
!     input matrix rin and store the result in rout. (rin=rout is
!     permitted) and to compute a condition number estimate.
!     cnb=frob.norm(r)*frob.norm(r**-1).

!     The frobenius norm is the square root of the sum of squares
!     of the elements. this condition number bound is used as
!     an upper bound and it acts as a lower bound on the actual
!     condition number of the problem. (see the book 'Solving least
!     squares', by lawson and hanson)
!
!     cnb/n <= condition number <= n*cnb
!
!     when rin is singular, rincon computes an inverse
!     corresponding to the the system having zero rows
!     and columns associated with the singular diagonal
!     entries.
!
!         cognizant persons:
!
!              dr. gerald j. bierman / keith h. bierman
!              factorized estimation applications
!              20 october 1981
!
!              revised 13 april 1982
!              revised 17 june  1982
!              version 0.01.00
!
!              revised 2 august 1982          version 0.01.01
!                                             version 1.00.00
!              revised 20      may 1985       version 2.00.00
!              revised  3      dec 1986       version 2.01.00
!
!     This "_m" version is a modified version of the original ESL routine
!        Modifications in the rincon_m version by B. Gibbs 8/27/09

!         1) upgraded to Fortran 90 syntax,
!         2) implemented vector operations using implicit MATLAB/Fortran 90 vector
!            operators (thus facilitating porting to MATLAB),
!         3) specifies input/output INTENT on all calling arguments,
!         4) rearranges calling arguments as “output, input/output, and input” to make
!            the subroutines appear more like functions,
!         5) uses a separate module to specify whether arithmetic operations are single
!            or double precision, and

! ******************************************************************
! *                                                                *
! *        subroutine rincon is a part of the fea. inc             *
! *              estimation subroutine library                     *
! *                                                                *
! *                         (esl)                                  *
! *                                                                *
! *             copyright 1982, 1983, 1984, 1985                   *
! *          factorized estimation applications inc.               *
! *   this library is licensed under the creative commons          *
! *   attribution 3.0 unported license.                            *
! *   to view a copy of this license,                              *
! *   visit http://creativecommons.org/licenses/by/3.0/            *
! *   or send a letter to creative commons, 171 second street,     *
! *   suite 300, san francisco, california, 94105, usa.            *
! *                                                                *
! ******************************************************************

      use filterPrec, only: real_kind1;   !precision (4 or 8)
      implicit none;

      integer(4),intent(in) :: n;                     !dimension of r matrices, n >= 1
      real(real_kind1),intent(in) :: rin(n*(n+1)/2);  !vector stored upper triangular matrix.
      real(real_kind1),intent(inout) :: cnb;          !if cnb is < 0, the condition number bound
                                                      ! calculation is skipped, and cnb is unchanged.
                                                      ! if cnb was non-negative on input then
                                                      ! cnb/n <= c <= cnb where c is the condition
                                                      ! number of rin. else it is unchanged.
      real(real_kind1),intent(out) :: rout(n*(n+1)/2);!vector stored upper triangular matrix.
                                                      ! contains inverse of rin) output.
                                                      ! if the compiler permits rout may overwrite rin.

      integer(4) i,j,k,ntot,jj,jjold,ii,ik,jm1,nbar;
      real(real_kind1) rnm,dinv,sum,one,z,temp,sumr;
!      epsilon is machine accuracy, used to test for zero divides:
!       0.6d-37 is suitable for VAX single precision,
!       1.0d-307 is limit for double precision on PC, but use 1d-100
      real(real_kind1) :: eps =1.0d-100;
      !___________________________________________

      z =0.d0;
      one =1.d0;
      ntot =(n*(n+1))/2;

!**     check to see if the rhs is included
      if (cnb >= 0.0) then;
         nbar = ntot;
         if (rin(ntot) == -1.0) nbar = ntot-n;     !???
         temp = maxval(abs(rin(1:nbar)));
         if (temp == 0.0) return;
         rnm = sum((rin(1:nbar)/temp)**2);
         cnb = temp*sqrt(rnm);
      endif;

      sumr     = rin(1);
      rout(1) = z;
      if (abs(sumr) > eps) rout(1) = one/sumr;

      if ( n > 1 ) then;
        jj=1;
        do j=2,n;
           jjold =jj;
           jj =jj+j;
           sumr = rin(jj);
           dinv = z;
           if (abs(sumr) > eps) dinv = one/sumr;
           rout(jj) =dinv;
           ii =0;
           ik =1;
           jm1 =j-1;
           do i=1,jm1;
              ii  =ii+i;
              ik  =ii;
              sumr = z;
              do k=i,jm1;
                 sumr =sumr +rout(ik)*rin(jjold+k);
                 ik =ik+k;
              enddo;
            rout(jjold+i) =-sumr*dinv;
            enddo;
          enddo;
      endif;

      if (cnb >= 0.0d0) then;
         temp = maxval(abs(rout(1:nbar)));
         rnm = sum((rout(1:nbar)/temp)**2);
         cnb = cnb*temp*sqrt(rnm);     !cnb =frob.norm(r) * frob.norm(r**-1)
      endif;

      return;
      end subroutine rincon_m;
