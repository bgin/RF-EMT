      subroutine rtinz (r,n,y,   x,ier);
!** solves transpose(r)*x =y for x by forward substitution.
!** Author: B. Gibbs, 12/2009

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
      end subroutine rtinz;
