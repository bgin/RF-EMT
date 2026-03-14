!**********************************************************************
      subroutine matInvGJ (ierr,  a,b,  np,n,m);

!**  matrix inversion and linear equation solver using Gauss-Jordan elimination
!**  Based on gaussj in Numerical Recipes in Fortran, Version 2

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

      integer(4),intent(in) :: np;        !row dimension of a and b (np >= n)
      integer(4),intent(in) :: n;         !number of rows and columns in a
      integer(4),intent(in) :: m;         !number of columns in b
      real(real_kind1),intent(inout) :: a(np,n);!matrix to be inverted, and inverse
      real(real_kind1),intent(inout) :: b(np,m);!right-hand-side matrix to be multiplied: b=inv(a)*b
      integer(4),intent(out) :: ierr;     !error return flag. 0=no error, 1=matrix singular

!****  local
      integer(4) i,j,k,icol,irow,indxc(n),indxr(n),ipiv(n);
      real(real_kind1) biga,atmp,pivinv,atmp1(max(n,m));
      !____________________________________________________________

      ierr =0;
      ipiv(:) =0;

      !***  compute pivot vector
      do i=1,n;
        biga =0.d0;

        do j=1,n;
          if (ipiv(j) /= 1) then;
            do k=1,n;
              if (ipiv(k) == 0) then;
                if (abs(a(j,k)) >= biga) then;
                  biga =abs(a(j,k));
                  irow =j;
                  icol =k;
                endif;
              else if (ipiv(k) > 1) then;
                ierr =1;       !stop 'matInvGJ error: singular matrix'
                return;
              endif;
            enddo;
          endif;
        enddo;

        ipiv(icol) =ipiv(icol)+1;

        !*** interchange rows if needed to put the pivot element on diagonal.
        !*** columns are relabeled, not physicaly interchanged

        if (irow /= icol) then;
          atmp1(:n) =a(irow,:);
          a(irow,:) =a(icol,:);
          a(icol,:) =atmp1(:n);
          atmp1(:m) =b(irow,:);
          b(irow,:) =b(icol,:);
          b(icol,:) =atmp1(:m);
        endif;

        indxr(i) =irow;
        indxc(i) =icol;

        if (a(icol,icol) == 0.d0) then;
          ierr =1;       !stop 'matInvGJ error: singular matrix'
          return;
        endif;
        pivinv =1.d0/a(icol,icol);
        a(icol,icol) =1.d0;
        a(icol,:) =a(icol,:)*pivinv;
        b(icol,:) =b(icol,:)*pivinv;

        !*** reduce the rows, except for pivot
        do j=1,n;
          if (j /= icol) then;
            atmp =a(j,icol);
            a(j,icol) =0.d0;
            a(j,:) =a(j,:)-a(icol,:)*atmp;
            b(j,:) =b(j,:)-b(icol,:)*atmp;
          endif;
        enddo;
      enddo;!___ end do i loop

      !** now reset column order to original
      do j=n,1,-1;
        if (indxr(j) /= indxc(j)) then;
          do k=1,n;
            atmp =a(k,indxr(j));
            a(k,indxr(j)) =a(k,indxc(j));
            a(k,indxc(j)) =atmp;
          enddo;
        endif;
      enddo;

      return;
      end subroutine matInvGJ;
