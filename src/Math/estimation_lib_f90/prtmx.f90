! **********************************************************************
!
!       **************
!       *            *
!       *   prtmx    *
!       *            *
!       **************
!
      subroutine prtmx (a,nn,n,l,ax);
!
!     subroutine prtmx prints a rectangular matrix, one row at atime
!     using sparse matrix print logic.

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

      character(*),intent(in) :: ax;   !string describing matrix a.
      integer(4),intent(in) :: nn;     !row dimension of a in the calling program.
      integer(4),intent(in) :: n;      !number of rows in matrix a.
      integer(4),intent(in) :: l;      !number of columns in matrix a.
      real(real_kind1),intent(in) :: a(nn,l);   !matrix a(n,l)

!       local
      integer(4) i,j,k,na,ia(l);
      real(4) b(l);
      !________________________________________________

      write(6,101) ax,n,l;
 101  format(/'matrix ',a,' of dimension ',i6,' by ',i6)

      if (n > nn)  then;
        write (6,'("Prtmx dimensioning error: nn=",i5,", n =",i5)') nn,n;
        call errhand ('prtmx error: dimensioning');
      endif;

      if (l <= 10)  then;
        do i=1,n;
          write(6,100) i,(a(i,j),j=1,l);
        enddo;
 100    format(1x,i4,10g12.5,5(/4x,10g12.5))
        return;
      endif;

      na =0;
      do i=1,n;
        k =0;
        do j=1,l;
          if(a(i,j) /= 0.d0) then;
            k =k+1;
            ia(k) =j;
            b(k) =a(i,j);
          endif;
        enddo;
        na =na+k;

        if (k > l/3 .or. l <= 10) then;
          write(6,200) i,(a(i,j),j=1,l);
  200     format(/1x,i4,10g12.5,9(/5x,10g12.5))
        else;
!          if (k == 0)write(6,100) i
          if (k > 0) write(6,300) i,(b(j),ia(j),j=1,k);
  300     format(/1x,i4,5(g12.5,' (',i4,')',6x),                        &
     &          6(/5x,5(g12.5,' (',i4,')',6x)))
        endif;

      enddo;

      if (na == 0) write (6,'(" *** MATRIX IS ALL ZERO")');

      return;
      end subroutine prtmx;
