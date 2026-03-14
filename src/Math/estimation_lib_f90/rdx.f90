!*************************************************
      subroutine rdx (del,  r,dx,n);
!        rdx computes del =r*dx where r is the upper
!        triangular square-root information matrix

!        B. Gibbs 8/27/09

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
      use filterPrec, only: real_kind1;            !precision (4 or 8)
      implicit none;

      integer(4),intent(in) :: n;                  !state dimension
      real(real_kind1),intent(in) :: r(n*(n+1)/2); !srif r matrix
      real(real_kind1),intent(in) :: dx(n);        !vector
      real(real_kind1),intent(out) :: del(n);      !computed solution

      integer(4) i,j,k;
      real(real_kind1) rsum;
      !_______________________________________________

      k=0;
      do i=1,n;
        rsum =0.d0;
        k =(i*(i+1))/2;
        do j=i,n;
          rsum =rsum +r(k)*dx(j);
          k =k+j;
        enddo;
        del(i) =rsum;
      enddo;

      return;
      end subroutine rdx;
