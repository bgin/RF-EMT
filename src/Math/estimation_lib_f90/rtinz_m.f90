      subroutine rtinz_m (x,ier,  r,n,y);
!** solves transpose(r)*x =y for x by forward substitution.
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
      implicit none;
      integer(4),intent(in) :: n;            !dimension of r
      real(8),intent(in) :: r(n*(n+1)/2);    !upper triangular matrix
      real(8),intent(in) :: y(n);            !RHS vector
      real(8),intent(out) :: x(n);           !solution vector
      integer(4),intent(out) :: ier;         !error status: 0 =OK, i = R is singular at row i

      integer(4) i,j,k;
      real(8) ysum;
      !_______________________________________

      ier =0;
      if (r(1) == 0.d0) then;
        ier =1;
        return;
      endif;

      x(1) =y(1)/r(1);
      k =1;
      do i=2,n;
        ysum =0.d0;
        do j=1,i-1;
          ysum =ysum +r(k+j)*x(j);
        enddo;
        k =k+i;
        if (r(k) == 0.d0) then;
          ier =i;
          return;
        endif;
        x(i) =(y(i)-ysum)/r(k);
      enddo;

      return;
      end subroutine rtinz_m;
