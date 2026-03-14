!**********************************************************
      subroutine cgnr (ier,  xe,  h,y,n,m,mmax);
!** preconditioned conjugate gradient method for solving the least squares normal equations
!** to minimize the residual.  The measurement equation is y =H*xe +r.
!** Least squares normal equations are xe = (H^T*H)^(-1) * H^T*y.
!** References:
!**    1) C.T. Kelley, Iterative Methods for Linear and Nonlinear Equations,
!**       SIAM, Philadelphia (1995),
!**    2) Å Björck, Numerical Methods for Least Squares Problems, SIAM, Philadelphia (1996)

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

      integer(4),intent(in) :: n;         !column dimension of matrix H
      integer(4),intent(in) :: m;         !number of actual measurements (rows) in H and y
      integer(4),intent(in) :: mmax;      !row dimension of H
      real(8),intent(in) :: h(mmax,n);    !measurement partial matrix
      real(8),intent(in) :: y(m);         !measurement vector
      real(8),intent(inout) :: xe(n);     !state vector
      integer(4),intent(out) :: ier;      !error return: 0 =OK, 1 =not converged in 100 iterations

      integer(4) i,j;
      integer(4) :: imax =100;            !max allowed iterations
      real(8) hs(m,n),r(m),q(m),rho0,tau,taul,alpha;
      real(8) d(n),p(n),s(n);
      real(8) :: eps =1.d-12;
      !_____________________________________

      ier =0;
      !** compute pre-conditioning as diagonal
      do i=1,n;
        d(i) =sqrt(sum(h(:m,i)**2));
!        d(i) =1.d0      !### test (little difference)
        hs(:,i) =h(:m,i)/d(i);
      enddo;
      xe(:) =xe(:)*d(:);

      r(:) =y(:)-matmul(hs,xe);
      s(:) =matmul(r(:),hs);
      p(:) =s(:);
      rho0 =sqrt(sum(s(:)**2));
      taul =1.d0;

      do i=1,imax;
        tau =sum(s(:)**2);
        if (i > 1) then;
          p(:) =s(:) +(tau/taul)*p(:);
        endif;
        taul =tau;
        q(:) =matmul(hs,p);
        alpha =tau/sum(q(:)**2);
        xe(:) =xe(:) +alpha*p(:);
        s(:) =s(:) -alpha*matmul(q,hs);
!       write (6,'("cgnr: ",2i3,2g13.5)') n,i,sqrt(sum(s(:)**2)),rho0;
        if (sqrt(sum(s(:)**2)) < eps*rho0) exit;
        if (i >= imax) then;
          ier =1;
          xe(:) =xe(:)/d(:);
          return;
        endif;
      enddo;

      xe(:) =xe(:)/d(:);

      return;
      end subroutine cgnr;
