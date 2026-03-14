!*******************************************************
      real(8) function condest (r,n);
!** L1 condition number estimate of upper triangular matrix using the LINPACK
!** algorithm contained in Cline, Moler, Stewart and Wilkinson (1979).
!** Reference: Algorithm A3.3.1, J.E. Dennis and R.B. Schnabel,
!** "Numerical Methods for Unconstrained Optimization and Nonlinear Equations",
!** Prentice-Hall (1983); republished by SIAM, Philadelphia (1996)

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
      integer(4),intent(in) :: n;         !matrix dimension
      real(8),intent(in) :: r(n*(n+1)/2); !upper triangular positive semi-definite matrix
                                          ! stored in a vector as upper triangular by columns

      integer(4) i,j,k,i1,ierr;
      real(8) temp,tempm,est,d(n),p(n),pm(n),x(n),rx(n),xp,xm,xnorm;
      !____________________________________

      k =1;
      est =r(1);
      d(1) =r(1);
      do i=2,n;
        temp =sum(abs(r(k+1:k+i)));
        est =max(est,temp);
        k =k+i;
        d(i) =r(k);
      enddo;

      x(1) =1.d0/d(1);
      k =1;
      do i=1,n;
        p(i) =r(k)*x(1);
        k =k+i;
      enddo;

      k =1;
      do j=2,n;
        k=k+i;
        xp =(1.d0-p(j))/d(j);
        xm =(-1.d0-p(j))/d(j);
        temp =abs(xp);
        tempm =abs(xm);
        i1 =(j*(j+1))/2;
        do i=j+1,n;
          i1 =(i*(i-1))/2;
          pm(i) =p(i)+r(j+i1)*xm;
          tempm =tempm +(abs(pm(i))/abs(d(i)));
          p(i) =p(i) +r(j+i1)*xp;
          temp =temp +(abs(p(i))/abs(d(i)));
          i1 =i1+i;
        enddo;
        if (temp >= tempm) then;
          x(j) =xp;     !e(j) =1
        else;           !e(j) =-1
          x(j) =xm;
          p(j+1:n) =pm(j+1:n);
        endif;
      enddo;

      xnorm =sum(abs(x(1:n)));
      est =est/xnorm;
      call rinz_m (rx,ierr,  r,n,x);   !rx =inv(r)*x  rinz(r,n,z,x,ierr)
      if (ierr /= 0) stop 'rinz error in condest';
      xnorm =sum(abs(rx(:n)));
      condest =est*xnorm;

      return;
      end function condest;
