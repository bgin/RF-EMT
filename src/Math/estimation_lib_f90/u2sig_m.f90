      subroutine u2sig_m (sig,  u,n,sf,jstrt,nsig);

!     computes standard deviations (sigmas) from u-d covariance factors.
!     Based on Bierman/Nead's ESP routine (B. Gibbs - 8/24/09)

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
      use filterPrec, only: real_kind1;   !
      implicit none;

      integer(4),intent(in) :: n;               !state dimension  (n >= 1)
      integer(4),intent(in) :: jstrt;           !start location of sigmas that are to be computed
      integer(4),intent(in) :: nsig;            !number of sigmas to be computed (nsig+jstrt-1 <= n)
      real(real_kind1),intent(in) :: sf(n);     !scale factors. if sf(1) < 0., scale factors are not
                                                ! used (implicitly treated as 1)
      real(real_kind1),intent(in) :: u(n*(n+1)/2);!input vector stored array containing the u-d factors.
!                                                 !  the d (diagonal) elements are stored on the diagonal
      real(real_kind1),intent(out) :: sig(nsig);  !vector of output standard deviations.

      integer(4) i,j,k,jsp1,jj,jjl,km1;
      real(8) s;
      !_____________________________________________________

      jsp1  = jstrt+1;
      jj    = jstrt*jsp1/2;
      k     = 1;
      sig(1) =u(jj);

      if (jsp1 <= n) then;
        do j=jsp1,n;
           k   = k+1;
           jjl = jj+jstrt-1;
           jj  = jj+j;
           s   = u(jj);
           km1 = nsig;
           if (k <= nsig) then;
              sig(k) = s;
              km1    = k-1;
           endif;
           sig(1:km1) =sig(1:km1) +s*u(jjl+1:jjl+km1)**2;
        enddo;
      endif;

!**       convert sig variances to sigma
      sig(:nsig) = sqrt(sig(:nsig));
      if (sf(1) >= 0.0) then;
        sig(:nsig) =sig(:nsig)*sf(:nsig);
      endif;

      return;
      end subroutine u2sig_m;
