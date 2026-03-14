! *************************************************************
      subroutine cov2ud_m (u,n);
!
!        Computes the U-D factors of a positive semi-definite input matrix
!        (upper triangular by columns U). The input vector stored
!        matrix is overwritten by the output U-D factors, which are also vector stored.
!        Singular input covariances result in output matrices with zero columns

!        This is a modification of the JPL ESP routine COV2UD written by
!        G.J. Bierman & R.A.Jacobson (feb. 1977)

!        Modifications in the cov2ud_m version by B. Gibbs 8/27/09
!         1) upgrading to Fortran 90 syntax,
!         2) specifying input/output INTENT on all calling arguments,
!         3) using a separate module to specify whether arithmetic operations are single
!            or double precision

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
      use filterPrec, only: real_kind1;               !precision 4 or 8
      implicit  none;

      integer(4),intent(in) :: n;                     !matrix dimension, n >= 2
      real(real_kind1),intent(inout) :: u(n*(n+1)/2); !input vector stored covariance matrix.
                                                      ! on output it contains the vector
                                                      ! stored u-d covariance factors.

      integer(4) i,j,k,l,none,np2,jj,ij,ik,kj,kk,jm1;
      real(real_kind1) zero,one,alpha,beta;
      !_______________________________________________

      zero =0.d0;
      one =1.d0;
      none =1;

      jj =n*(n+1)/2;
      np2 =n+2;

      do l=2,n;
        j =np2-l;
        alpha =zero;
        if (u(jj) < zero) then;
          write (6,100) j,u(jj);
          u(jj) =zero;
        endif;
        if (u(jj) > zero) alpha =one/u(jj);
        jj =jj-j;
        kk =0;
        kj =jj;
        jm1 =j-1;

        do k=1,jm1;
          kj =kj+1;
          beta =u(kj);
          u(kj) =alpha*u(kj);
          ij =jj;
          ik =kk;
          do i=1,k;
            ik =ik+1;
            ij =ij+1;
            u(ik) =u(ik) -beta*u(ij);
          enddo;
          kk =kk+k;
        enddo;
      enddo;  !___ end l loop

      if (u(1) < zero) then;
        write (6,100) none, u(1);
        u(1) =zero;
      endif;

      return;

  100 format (10x,'in cov2ud at step',i4,', diagonal entry =',e12.4)
      end subroutine cov2ud_m;
