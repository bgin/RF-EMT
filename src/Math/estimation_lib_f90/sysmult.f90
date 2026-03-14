! *****************************************************************
      subroutine sysmult (c,  a,b,n);
!         symult2 multiplies  c=a*b where b and c are vectors and
!         a is a symmetric matrix stored as upper triangular by
!         columns.  b is assumed to be sparse !
!         symult can be called repeatedly to form a matrix multiply.

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

      !** local
      integer(4) i,j,k,k1,k2,nb,ib(n);
      real(real_kind1) :: sump,d(n),bb;
      !___________________________________________

      nb =0;
      do i=1,n;
        if (b(i) /= 0.d0) then;
          nb =nb+1;
          d(nb) =b(i);
          ib(nb) =i;
        endif;
      enddo;

      c(:) =0.d0;
      do i=1,nb;
        j =ib(i);  !row index of b, column index of a
        bb =d(i);
        k1 =(j*(j-1))/2;
        c(1:j) =c(1:j) +bb*a(k1+1:k1+j);
        if (j < n) then;
          k2 =k1+j+j;
          do k=j+1,n;
            c(k) =c(k) +bb*a(k2);
            k2 =k2+k;
          enddo;
        endif;
      enddo;

      return;
      end subroutine sysmult;
