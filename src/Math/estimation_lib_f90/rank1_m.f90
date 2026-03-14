      subroutine rank1_m (uout,ierr,  v,  uin,n,c);
!
!     stable  u-d factor rank 1 update
!              (uout)*dout*(uout)**t = (uin)*din*(uin)**t + c*v*v**t
!
!     This "_m" version is a modified version of the original Bierman/Nead ESP routine
!        Modifications in the rank1_m version by B. Gibbs 8/27/09
!         1) upgraded to Fortran 90 syntax,
!         2) implemented vector operations using implicit MATLAB/Fortran 90 vector
!            operators (thus facilitating porting to MATLAB),
!         3) specifies input/output INTENT on all calling arguments,
!         4) rearranges calling arguments as “output, input/output, and input” to make
!            the subroutines appear more like functions,
!         5) uses a separate module to specify whether arithmetic operations are single
!            or double precision, and
!         6) tests for matrix sparseness and only operate on non-zero elements

!**************************************
! The author grants the user a non-exclusive, worldwide, royalty-free copyright license to

! 1. reproduce and modify the software for your own purposes, and
! 2. to distribute the software provided that you give credit to the author,
!    do not make claims against the author or Wiley-Interscience,
!    and mark any modifications as your own.  If the software is incorporated in
!    a commercial product, you  agree to defend and indemnify this author and
!    Wiley-Interscience against any losses, damages and costs arising from claims,
!    lawsuits and other legal actions brought by a third party.

! The software is provided on an “as is” basis, without warranties or conditions
! of any kind, either express or implied (including any warranties or conditions
! of title, non-infringement, merchantability or fitness for a particular purpose).
! The user is solely responsible for determining the appropriateness of using and
! distributing the program and assumes all risks associated with its exercise of
! rights under this agreement.
!**************************************
      use filterPrec, only: real_kind1;      !precision (4 or 8)
      implicit none;

      integer(4),intent(in) :: n;            !matrix dimension, n >= 1
      real(real_kind1),intent(in) :: c;      !input scalar. should be non-negative
      real(real_kind1),intent(inout) :: v(n);!input vector for rank one modification.
                                             ! v is destroyed during the process
      real(real_kind1),intent(in) :: uin(n*(n+1)/2);
                                             !input vector stored positive semi-definite u-d
                                             ! array, with d elements stored on the diagonal
      real(real_kind1),intent(out) :: uout(n*(n+1)/2);
                                             !output vector stored positive (possibly)
                                             ! semi-definite u-d result.
                                             ! compiler permitting uout=uin is allowed
      integer(4),intent(out) :: ierr;        !0 =normal return:
                                             ! j = error return due to negative computed j-th diagonal
                                             ! -1 = ud result is singular

      integer(4) np2,l,i,j,jj,jm1;
      real(real_kind1) alpha, beta, s, d, zero, cj;
!      epsilon is machine accuracy: 0.6d-37 is suitable for VAX single precision,
!       1.0d-307 is limit for double precision on PC, but use 1d-100
      real(real_kind1) :: eps =1.0d-100;
      real(real_kind1) :: tst =.0625d0;      !used for rank1 algorithm switching
      !___________________________________________________

      zero    = 0.d0;
      ierr = 0;

      jj = n*(n+1)/2;
      cj = c;
      if (cj <= zero) then;
        uout(:jj) =uin(:jj);
        return;
      endif;

      if (n > 1) then;
        np2 =n+2;
        do l=2,n;
          j =np2-l;
          s =v(j);
          beta =cj*s;
          d =uin(jj)+beta*s;
          if (d <= eps) then;
            if (d < zero) then;
              ierr = j;
              return;
            endif;
            jj = jj-j;
            ierr = -1;
            uout(jj+1:jj+j) = zero;

          else;
            beta     = beta/d;
            alpha    = uin(jj)/d;
            cj       = alpha*cj;
            uout(jj) = d;
            jj       = jj-j;
            jm1      = j-1;
            if (alpha >= tst) then;
               v(1:jm1)  = v(1:jm1) -s*uin(jj+1:jj+jm1);
               uout(jj+1:jj+jm1) = beta*v(1:jm1) + uin(jj+1:jj+jm1);
            else;
              do i=1,jm1;
                d = v(i) - s*uin(jj+i);
                uout(jj+i) = alpha*uin(jj+i) + beta*v(i);
                v(i) = d;
              enddo;
            endif;
          endif;

        enddo;  !___ end l,j loop
      endif;

      uout(1) = uin(1) + (cj*v(1))*v(1); !cj*v(1)**2

      return;
      end subroutine rank1_m;
