! ************************************************************************
      subroutine smult (c,  a,b, n,m,l, nra,nrb,nrc);

!***    smult computes matrix multiplication c = a * b
!***    using where a is sparse

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
      use filterPrec, only: real_kind1;   !precision (4 or 8)
      implicit none;
      integer(4),intent(in) :: n;      !number of rows in matrix a.
      integer(4),intent(in) :: m;      !input number of columns in matrix a and number of
                                       ! rows in matrix b.
      integer(4),intent(in) :: l;      !number of columns in matrix b.
      integer(4),intent(in) :: nra;    !row dimension of matrix a in the calling program
      integer(4),intent(in) :: nrb;    !row dimension of matrix b in the calling program
      integer(4),intent(in) :: nrc;    !row dimension of matrix c in the calling program
      real(real_kind1),intent(in) :: a(nra,m);  !matrix a.  size (n,m).
      real(real_kind1),intent(in) :: b(nrb,l);  !matrix b.  size (m,l).
      real(real_kind1),intent(out) :: c(nrc,l); !matrix c = a * b.  size (n,l).
                                       !matrices a and c may share the same storage locations

      !*** local
      integer(4) i,j,k,jj,ind(m);
      real(real_kind1) d(m),sumr;
      !______________________________________________________________

      do i=1,n;      !row index of matrix a
        !*** find non-zero elements in row a(i,:)
        jj =0;
        do j=1,m;    !column index of a
          if (abs(a(i,j)) > 1.d-60) then;
            jj =jj+1;
            d(jj) =a(i,j);
            ind(jj) =j;
          end if;
        enddo;

        do j=1,l;    !column index of matrix b
          sumr =0.d0;
          if (jj > 0)  then;
            do k=1,jj;
              sumr =sumr +d(k)*b(ind(k),j);
            enddo;
          endif;
          c(i,j) =sumr;
        enddo;
      enddo;

      return;
      end subroutine smult;
