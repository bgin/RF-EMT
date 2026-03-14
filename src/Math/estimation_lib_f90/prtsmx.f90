! ********************************************************************
      subroutine prtsmx (a,n,ax);

!** Subroutine prtsmx prints the upper triangular portion of a symmetric
!** square matrix, one row at a time using sparse matrix print logic.

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

      character(*),intent(in) :: ax;                  !matrix label
      integer(4),intent(in) :: n;                     !matrix dimension
      real(real_kind1),intent(in) :: a((n*(n+1))/2);  !symmetric matrix

!       local
      integer i,j,k,k1,k2,na,ia(n);
      real(real_kind1) b(n);
      !___________________________________________________

      write(6,101) ax,n;
 101  format('SYMMETRIC MATRIX ',a,' OF DIMENSION ',i6)

      k2 =0;
      na =0;

      do i=1,n;
        k=0;
        do j=1,i;
          if (a(k2+j) /= 0.d0) then;
            k =k+1;
            ia(k) =j;
            b(k) =a(k2+j);
          endif;
        enddo;

        na =na+k;
        k1 =k2+1;
        k2 =k2+i;
        if (k > i/3 .or. i <= 10) then;
          write(6,100) i,(a(j),j=k1,k2);
 100      format(1x,i4,10g12.5,9(/5x,10g12.5))
        else;
!        if (k == 0) write(6,100) i
         if (k > 0) write(6,300) i,(b(j),ia(j),j=1,k);
  300     format(/1x,i4,5(g12.5,' (',i4,')',6x),                        &
     &          6(/5x,5(g12.5,' (',i4,')',6x)))
        endif;
      enddo;

      if (na == 0) write (6,'(" *** MATRIX IS ALL ZERO")');

      return;
      end subroutine prtsmx;
