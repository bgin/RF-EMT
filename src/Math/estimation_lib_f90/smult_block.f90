! ************************************************************************
      subroutine smult_block (c,  a,b, n,m,k, nra,nrb,nrc);

!***    smult_block computes matrix multiplication c = a * b
!***    using where a is stored sparse as sparse vectors a,

!** Author: B. Gibbs, 1/2010

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
      implicit none;
      integer(4),intent(in) :: n;      !number of rows in matrix a.
      integer(4),intent(in) :: m;      !input number of columns in matrix a and number of
                                       ! rows in matrix b.
      integer(4),intent(in) :: k;      !number of columns in matrix b.
      integer(4),intent(in) :: nra;    !row dimension of matrix a in the calling program
      integer(4),intent(in) :: nrb;    !row dimension of matrix b in the calling program
      integer(4),intent(in) :: nrc;    !row dimension of matrix c in the calling program
      real(8),intent(in) :: a(nra,m);  !matrix a.  size (n,m).
      real(8),intent(in) :: b(nrb,k);  !matrix b.  size (m,l).
      real(8),intent(out) :: c(nrc,k); !matrix c = a * b.  size (n,l).
                                       !matrices a and c may share the same storage locations

      !*** local
      integer(4) i,j,jr1,jr2;
      !______________________________________________________________

      do i=1,n;      !row index of matrix a
        !*** find first and last non-zero elements in row a(i,:)
        jr1 =0;
        jr2 =0;
        do j=1,m;    !column index of a
          if (abs(a(i,j)) > 1.d-60) then;
            jr2 =j;
            if (jr1 == 0) jr1 =j;
          end if;
        enddo;

        if (jr1 > 0) then;
          c(i,:) =matmul(a(i,jr1:jr2),b(jr1:jr2,:));
        else;
          c(i,:) =0.d0;
        endif;
      enddo;

      return;
      end subroutine smult_block;
