      subroutine ud2sig_m (sig,  u,n,nsig);
!     computes standard deviations (sigmas) from u-d covariance factors
!     (B. Gibbs - 8/24/09)
!
!     This "_m" version is an independently-delevoped version of ESP/ESL routine ud2sig
!        Features in version by B. Gibbs 8/27/09

!         1) Fortran 90 syntax,
!         2) implemented vector operations using implicit MATLAB/Fortran 90 vector
!            operators (thus facilitating porting to MATLAB),
!         3) specifies input/output INTENT on all calling arguments,
!         4) rearranges calling arguments as “output, input/output, and input” to make
!            the subroutines appear more like functions,
!         5) uses a separate module to specify whether arithmetic operations are single
!            or double precision, and

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
      use filterPrec, only: real_kind1;         !precision (4 or 8)
      implicit none;

      integer(4),intent(in) :: n;               !state dimension  (n >= 1)
      integer(4),intent(in) :: nsig;            !number of sigmas to be computed (nsig <= n)
      real(real_kind1),intent(in) :: u(n*(n+1)/2);!input vector stored array containing the u-d factors.
                                                !  the d (diagonal) elements are stored on the diagonal
      real(real_kind1),intent(out) :: sig(nsig);!vector of output standard deviations.

      integer(4) i,j,k,jj,jjl,jm1;
      real(real_kind1) s;
      !_____________________________________________________

      jj    = 1;
      sig(1) =u(jj);
      if (n > 1) then;
        do j=2,n;
           jjl = jj;
           jj  = jj+j;
           s   = u(jj);
           jm1 = nsig;
           if (j <= nsig) then;
              sig(j) = s;
              jm1    = j-1;
           endif;
           sig(1:jm1) =sig(1:jm1) +s*u(jjl+1:jjl+jm1)**2;
        enddo;
      endif;

      sig(1:nsig) =sqrt(sig(1:nsig));

      return;
      end subroutine ud2sig_m;
