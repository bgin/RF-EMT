! *****************************************************************
      subroutine symult (c,  a,b,n);
!         symult multiplies  c=a*b where b and c are vectors and
!         a is a symmetric matrix stored as upper triangular by
!         columns.  symult can be called repeatedly to form a matrix mult.

!** Author: B. Gibbs, 12/2009

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
      use filterPrec, only: real_kind1;             !precision (4 or 8)
      implicit none;
      integer(4),intent(in) :: n;                   !dimension of a
      real(real_kind1),intent(in) :: a((n*(n+1))/2);!symmetric matrix
      real(real_kind1),intent(in) :: b(n);          !column vector
      real(real_kind1),intent(out) :: c(n);         !column vector

      integer(4) j,j1,k,k1,k2;
      real(real_kind1) suma;
      !__________________________________________

      k1=0;
      do j=1,n;
        suma=0.d0;
        do k=1,j;
          k1=k1+1;
          suma =suma+a(k1)*b(k);
        enddo;

        if (j < n) then;
          j1 =j+1;
          k2 =k1+j;
          do k=j1,n;
            suma =suma+a(k2)*b(k);
            k2 =k2+k;
          enddo;
        endif;

        c(j)=suma;
      enddo;

      return;
      end subroutine symult;
