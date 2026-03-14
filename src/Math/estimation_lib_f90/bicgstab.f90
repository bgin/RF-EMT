!**********************************************************
      subroutine bicgstab (ier,  xe,  a,b,n,nd);
!** preconditioned bi-conjugate gradient stabilized method for solving A * xe =b.
!** where A is n x n.  Note version assumes that matrix A is full.  For most problems of this
!** type matrix A is very sparse, so sparse matrix storage should be used for A and the
!** internal multiplications modified accordingly.
!** Reference: C.T. Kelley, Iterative Methods for Linear and Nonlinear Equations, SIAM, Philadelphia (1995)

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

      integer(4),intent(in) :: n;         !row and column dimension of A
      integer(4),intent(in) :: nd;        !row dimension of A used in the calling routine (nd >= n)
      real(8),intent(in) :: a(nd,n);      !square n x n matrix
      real(8),intent(in) :: b(n);         !n x 1 right-hand-side vector
      real(8),intent(inout) :: xe(n);     !output n x 1 vector
      integer(4),intent(out) :: ier;      !error code: 0 =no error, 1=maximum iterations (200) reached

      integer(4) i,j;
      integer(4) :: imax=200;
      real(8) d(n),hs(n,n),r(n),r0(n),rho0,rhol,rho,alpha,beta,w;
      real(8) p(n),v(n),s(n),t(n);
      real(8) :: eps =1.d-6;              !tolerance test for convergence
      !_____________________________________

      ier =1;
      !** compute pre-conditioning as diagonal
      do i=1,n;
        d(i)=sqrt(sum(a(:,i)**2));
        hs(:,i) =a(:,i)/d(i);
      enddo;
      xe(:) =xe(:)/d(:);

      r(:) =b(:) -matmul(hs,xe);
      r0(:) =r(:);
      rho0 =sqrt(sum(b(:)**2));
      rhol =1.d0;
      alpha =1.d0;
      w =1.d0;
      v(:) =0.d0;
      p(:) =0.d0;
      rho =dot_product(r0,r);

      do i=1,imax;
        beta =(rho/rhol)*(alpha/w);
        p(:) =r(:) +beta*(p(:)-w*v(:));
        v(:) =matmul(hs,p);
        alpha =rho/dot_product(r0,v);
        s(:) =r(:) -alpha*v(:);
        t(:) =matmul(hs,s);
        w =dot_product(t,s)/sqrt(sum(t(:)**2));
        rho =w*dot_product(r0,t);
        xe(:) =xe(:) +alpha*p(:) +w*s(:);
        r(:) =s(:) -w*t(:);
        if (sqrt(sum(r(:)**2)) < eps*rho0) exit;
        if (i >= imax) then;
          ier =1;
          return;
        endif;
      enddo;

      ier =0;
      xe(:) =xe(:)*d(:);

      return;
      end subroutine bicgstab;
