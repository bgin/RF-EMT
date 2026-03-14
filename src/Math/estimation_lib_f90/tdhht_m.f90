      subroutine tdhht_m (s,v,  imaxs,irs,jcs,jstart,jstop);
!
!       tdhht transforms a rectangular double subscripted matrix s
!       to an upper triangular or partially upper triangular form
!       by the application of householder orthogonal
!       transformations. it is assumed that the first "jstart"-1
!       columns of s are already triangularized. the algorithm is
!       described in "factorization methods for discrete sequential
!       estimation" by g.j.bierman, academic press, 1977
!
!  *** required subroutines: abs,  sqrt   <forlib>
!
!         cognizant persons:
!
!              dr. gerald j. bierman / keith h. bierman
!              factorized estimation applications
!              18 october 1978
!              revised 29 may 1982        version 0.01.00
!                                         version 1.00.00
!              revised 05 may 1984        version 2.00.00
!              revised 16 may 1986        version 2.01.00
!              revised 17 oct 1986        version 2.02.00
!
!     This "_m" version is a modified version of the original ESL routine
!        Modifications in the tdhht_m version by B. Gibbs 8/27/09

!         1) upgraded to Fortran 90 syntax,
!         2) implemented many vector operations using implicit MATLAB/Fortran 90 vector
!            operators (thus facilitating porting to MATLAB),
!         3) specifies input/output INTENT on all calling arguments,
!         4) rearranges calling arguments as "output, input/output, and input" to make
!            the subroutines appear more like functions,
!         5) uses a separate module to specify whether arithmetic operations are single
!            or double precision, and
!         6) tests for matrix sparseness and only operate on non-zero elements

! ******************************************************************
! *                                                                *
! *        subroutine tdhht  is a part of the fea. inc             *
! *              estimation subroutine library                     *
! *                                                                *
! *                         (esl)                                  *
! *                                                                *
! *           copyright 1982, 1983, 1984, 1985, 1986               *
! *          factorized estimation applications inc.               *
! *   this library is licensed under the creative commons          *
! *   attribution 3.0 unported license.                            *
! *   to view a copy of this license,                              *
! *   visit http://creativecommons.org/licenses/by/3.0/            *
! *   or send a letter to creative commons, 171 second street,     *
! *   suite 300, san francisco, california, 94105, usa.            *
! *                                                                *
! ******************************************************************
!
      use filterPrec, only: real_kind1;               !precision (4 or 8)
      implicit none;

      integer(4),intent(in) :: imaxs;                 !row dimension of s
      integer(4),intent(in) :: irs;                   !number of rows in s  (irs <= imaxs and irs >= 1)
      integer(4),intent(in) :: jcs;                   !number of columns in s
      integer(4),intent(in) :: jstart;                !index of the first column to be triangularized. if
                                                      ! jstart < 1 it is assumed that jstart=1, i.e.
                                                      ! start triangularization at column 1.
      integer(4),intent(in) :: jstop;                 !index of last column to be triangularized.
                                                      ! if jstop < jstart or jstop > jcs then
                                                      ! if irs <= jcs jstop is set equal to irs-1
                                                      ! if irs > jcs jstop is set equal to jcs
                                                      ! i.e. the triangularization is completed as far as possible
      real(real_kind1),intent(inout) :: v(irs);       !work vector
      real(real_kind1),intent(inout) :: s(imaxs,jcs); !input (possibly partially) triangular matrix. the
                                                      ! output (possibly partially) triangular result
                                                      ! overwrites the input.

      integer(4) i,j,k,jp1,jstt,jstp,irs1;
      real(real_kind1) :: vsum, delta,deps,one,zero,tmpmax,tmp;

!**      epsilon is machine accuracy: 0.6d-37 is suitable for VAX single precision,
!**      1.0d-307 is limit for double precision on PC, but use 1d-100
      real(real_kind1) :: eps =1.0d-100;
      !______________________________________________

      if (irs <= 1) return;

      deps = eps;
      one  = 1.d0;
      zero = 0.0d0;
      jstt = jstart;
      jstp = jstop;
      if (jstt < 1) jstt=1;
      if (.not. (jstp >= jstt .and. (jstp <= jcs)                       &
     &                 .and. (jstp <= irs))) then;
        if (irs <= jcs) jstp=irs-1;
        if (irs > jcs) jstp=jcs;
      endif;

      do j=jstt,jstp;         !column loop
         tmpmax = zero;
         jp1    = j+1;
         if (jp1 > irs) return;

!         tmpmax = maxval(abs(s(jp1:irs,j)))
         irs1 =0;
         do i=jp1,irs;
            tmp = abs(s(i,j));
            !**  tmpmax = max(s(i,j); i=(j+1),irs)
            if (tmpmax < tmp) tmpmax = tmp;
            if (tmp /= zero) irs1 =irs;   !row of last non-zero element
         enddo;

         !**  if tmpmax = 0, column j elements below j-th are
         !**  zero and this step of the algorithm is omitted
         if (tmpmax <= eps) cycle;

         tmp = abs(s(j,j));
         if (tmp > tmpmax) tmpmax = tmp;
         vsum = zero;
         v(:) =zero;
         do i = j, irs1;
            v(i)   = s(i,j);
            s(i,j) = zero;
            vsum    = vsum+(v(i)/tmpmax)**2;
         enddo;
         vsum = tmpmax * sqrt(vsum);
         !**  vsum = 2-norm of column j

         if (v(j) > zero) vsum = -vsum;
         s(j,j)  = vsum;

         !**  normalize householder vector, v
         v(j:irs1) = v(j:irs1)/vsum;
         v(j) = v(j) - one;
         vsum  = v(j);

         if (jp1 > jcs) cycle;

         !**  the householder transformation is T =I-(1/vsum)*v*v**t,
         !**  where v is norm scaled.
         do k = jp1,jcs;
            delta = dot_product(s(j:irs1,k),v(j:irs1));
            if (abs(delta) <= deps) cycle;
            s(j,k)  = s(j,k) + delta;
            delta   = delta/vsum;
            s(jp1:irs1,k)  = s(jp1:irs1,k) +delta*v(jp1:irs1);
         enddo;
      enddo;

      return;
      end subroutine tdhht_m;
