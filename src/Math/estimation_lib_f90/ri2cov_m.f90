      subroutine ri2cov_m (sig,covout,  rinv,n,krow,kcol);
!
!     Computes the standard deviations and, if desired, the
!     covariance matrix of a vector stored upper triangular square root
!     information matrix. the output covariance matrix is also vector stored.

!     This "_m" version is a modified version of the original Bierman/Nead ESP routine
!        Modifications in the ri2cov_m version by B. Gibbs 8/27/09

!         1) upgraded to Fortran 90 syntax,
!         2) specifies input/output INTENT on all calling arguments,
!         3) rearranges calling arguments as “output, input/output, and input” to make
!            the subroutines appear more like functions,
!         4) uses a separate module to specify whether arithmetic operations are single
!            or double precision.

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
      use filterPrec, only: real_kind1;                  !precision (4 or 8)
      implicit none;

      integer(4),intent(in) :: n;                        !dimension of the rinv matrix (>= 1)
      integer(4),intent(in) :: krow;                     !krow > 0: computes the covariance and sigmas
                                                         ! < 0: computes only the sigmas of the first
                                                         ! krow variables of the rinv matrix!
                                                         ! = 0:  no covariance, but all sigmas
                                                         !     (e.g., use all n rows of rinv)
      integer(4),intent(in) :: kcol;                     !number of columns of covout that are
                                                         ! computed. if kcol <= 0 then kcol=krow.
                                                         ! if krow <= 0 this input is ignored.
      real(real_kind1),intent(in) :: rinv(n*(n+1)/2);    !vector stored upper triangular
                                                         ! covariance square root.
                                                         ! (rinv = rinverse is the inverse of
                                                         ! the srif matrix)
      real(real_kind1),intent(out) :: sig(n);            ! output vector of standard deviations
      real(real_kind1),intent(out) :: covout(n*(n+1)/2); !vector stored covariance (if the compiler
                                                         ! permits covout may overwrite rinv)

      integer(4) lim,iks,i,j,k,jj,ik,nm1,kkol,ijs,jp1,imj,ijk;
      real(real_kind1) rsum;
      real(real_kind1) :: zero =0.0d0;
      !______________________________________________

      lim  = n;
      if (krow /= 0) lim  = abs(krow);

      !**  compute sigmas
      iks = 0;
      do j = 1, lim;
         iks = iks+j;
         rsum = zero;
         ik  = iks;
         do k=j, n;
            rsum = rsum +rinv(ik)**2;
            ik = ik+k;
         enddo;
         sig(j) = sqrt(rsum);
      enddo;

      if (krow <= 0) return;

      kkol = kcol;
      if (kkol <= 0) kkol = krow;
      jj = 0;

      if (n > 1) then;
         ! *** compute covariance
         nm1 = lim;
         if (krow == n) nm1 = n-1;
         do j=1, nm1;
            jj = jj+j;
            covout(jj) = sig(j)**2;
            ijs = jj+j;
            jp1 = j+1;
            do i = jp1, kkol;
               ik  = ijs;
               imj = i-j;
               rsum = zero;
               do k = i, n;
                  ijk = ik + imj;
                  rsum = rsum + rinv(ik)*rinv(ijk);
                  ik = ik+k;
               enddo;
               covout(ijs) = rsum;
               ijs = ijs+i;
            enddo;
         enddo;
      endif;

      if (krow == n) covout(jj+n) = sig(n)**2;

      return;
      end subroutine ri2cov_m;
