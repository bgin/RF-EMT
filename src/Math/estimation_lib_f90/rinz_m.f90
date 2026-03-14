! ************************************************************************
      subroutine rinz_m (x,ierr,  r,n,z);
!
!     Computes the solution vector x to the matrix equation
!     r*x=z, where r is an n by n upper-triangular vector
!     matrix. computation is minimized and accuracy enhanced
!     by not computing the matrix inverse. if r is singular,
!     with its jth diagonal element as its last diagonal zero
!     value, x(i) is computed  only for j < i <= n,
!     the other components,   x(i), are left unchanged for
!     1 <= i <= j, and the value of j is returned as ierr.
!
!     This version was converted to use Fortran 90 syntax and allows
!     variable real precision B. Gibbs, 8/24/09

! *** required subroutine:   abs    <forlib>
!
!         cognizant persons:
!              dr. gerald j. bierman / keith h. bierman
!              factorized estimation applications
!              01 may 1978
!
!              revised  3 august 1982         version 0.01.00
!                                             version 1.00.00
!              revised 20    may 1985         version 2.00.00
!
!     This "_m" version is a modified version of the original ESL routine
!        Modifications in the rinz_m version by B. Gibbs 8/27/09

!         1) upgraded to Fortran 90 syntax,
!         2) specifies input/output INTENT on all calling arguments,
!         3) rearranges calling arguments as “output, input/output, and input” to make
!            the subroutines appear more like functions,
!         4) uses a separate module to specify whether arithmetic operations are single
!            or double precision, and

! ******************************************************************
! *                                                                *
! *        subroutine rinz  is a part of the fea. inc              *
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

      use filterPrec, only: real_kind1;            !real precision (4 or 8)
      implicit none;
      integer(4),intent(in) :: n;                  !matrix vector dimension, n > 0
      real(real_kind1),intent(in) :: r(n*(n+1)/2); !vector stored upper triangular matrix
      real(real_kind1),intent(in) :: z(n);         !right hand side vector of equation r*x=z
      real(real_kind1),intent(out) :: x(n);        !solution vector that may overwrite z, compiler permitting.
      integer(4),intent(out) :: ierr;              !error return flag: 0 =normal return,
                                                   ! j = solution x valid only for components > j

      integer(4) i,j,k,ix,jj,ij;
      real(real_kind1)  rsum, s, dx;
!      epsilon is machine accuracy: 0.6d-37 is suitable for VAX single precision,
!       1.0d-307 is limit for double precision on PC, but use 1d-100
      real(real_kind1) :: eps =1.0d-100;
      real(real_kind1) :: zero = 0.0d0;
      !________________________________________________


      ierr = 0;
      jj   = n*(n+1)/2;
      ix   = n;
      rsum  = zero;
      if (n > 1) then;
         do i = 2,n;
            s  = r(jj);
            dx = z(ix) - rsum;
            if (abs(s) <= eps * abs(dx)) then;
              ierr =ix;
              return;
            endif;
            x(ix) = dx/s;
            ij    = jj-1;
            rsum   = zero;
            do j = ix, n;
               rsum  = rsum+r(ij)*x(j);
               ij   = ij+j;
            enddo;
            jj = jj-ix;
            ix = ix-1;
         enddo;
      endif;

      s = r(1);
      dx = z(1) - rsum;
      if (abs(s) <= eps * abs(dx)) then;
        ierr = ix;   ! *** error return ***
        return;
      endif;
      x(1) = dx/s;

      return;
      end subroutine rinz_m;
