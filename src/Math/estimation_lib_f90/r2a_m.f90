      subroutine r2a_m (a,  r,n,namr,imaxa,ja,nama);
!
!     to place the triangular vector stored matrix "r"
!     into the matrix "a" and to arrange the columns
!     to match the desired nama parameter list.
!     names in the nama list that do not correspond
!     to any name in the namr list have zero entries
!     in the corresponding "a" column
!
!     This "_m" version is a modified version of the original ESL routine
!        Modifications in the r2a_m version by B. Gibbs 8/27/09

!         1) upgrading to Fortran 90 syntax,
!         2) implementing many vector operations using implicit MATLAB/Fortran 90 vector
!            operators (thus facilitating porting to MATLAB),
!         3) specifying input/output INTENT on all calling arguments,
!         4) rearranging calling arguments as “output, input/output, and input” to make
!            the subroutines appear more like functions,
!         5) using a separate module to specify whether arithmetic operations are single
!            or double precision.

!          cognizant persons:
!
!             dr. gerald j. bierman/keith h. bierman
!             factorized estimation applications inc.
!
!             01 september  1976
!
!             revised 28 march   1982
!             revised  2 august  1982
!             revised 20 august  1982   version 0.00.02
!
!             revised 20      may 1985  version 2.00.00
!
! ******************************************************************
! *                                                                *
! *        subroutine r2a    is a part of the fea. inc             *
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
!
      use filterPrec, only: real_kind1;   !precision (4 or 8)
      implicit none;

      integer(4),intent(in) :: n;                  !dimension of r
      integer(4),intent(in) :: imaxa;              !row dimension of a
      integer(4),intent(in) :: ja;                 !number of parameter names associated
                                                   ! with the output "a" matrix
      character(*),intent(in) :: namr(n);          !parameter names associated with r
      character(*),intent(in) :: nama(ja);         !names associated with "a"
      real(real_kind1),intent(in) :: r(n*(n+1)/2); !input upper triangular vector stored array
      real(real_kind1),intent(out) :: a(imaxa,ja); !matrix containing the rearranged r matrix

      integer(4) i,j,jj;
      !____________________________________

      a(:n,:) =0.d0;

      do j=1, ja;
        do i=1,n;
          if (namr(i) == nama(j)) then;
            jj = i*(i-1)/2;
            a(1:i,j) = r(jj+1:jj+i);
            exit;
          endif;
        enddo;
      enddo;

      return;
      end subroutine r2a_m;
