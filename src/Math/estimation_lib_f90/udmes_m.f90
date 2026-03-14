! **********************************************************
      subroutine udmes1 (f,g,alpha,dz,  u,n,r,a,x,z);

!       computes estimate and u-d measurement updated
!       covariance, p=u*d*u**t
!       cognizant persons: g.j. bierman/m.w. nead (jpl, feb.1978)

!     This "_m" version is a modified version of the original ESP routine
!        Modifications in the udmes_m version by B. Gibbs 8/27/09

!         1) upgraded to Fortran 90 syntax,
!         2) implemented many vector operations using implicit MATLAB/Fortran 90 vector
!            operators (thus facilitating porting to MATLAB),
!         3) specifies input/output INTENT on all calling arguments,
!         4) rearranges calling arguments as “output, input/output, and input” to make
!            the subroutines appear more like functions,
!         5) uses a separate module to specify whether arithmetic operations are single
!            or double precision, and
!         6) tests for matrix sparseness and only operate on non-zero elements
!         7) original udmes was broken into udmes1 and udmes2 so that measurement outlier
!            test could be performed after udmes1 and udmes1 bypassed.  Also udmes1 output
!            is useful for some applications.

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
      use filterPrec, only: real_kind1;      !precision (4 or 8)
      implicit  none;
      integer(4),intent(in) :: n;            !dimension of the state estimate. n > 1
      real(real_kind1),intent(in) :: u(n*(n+1)/2); !upper triangular matrix, with d
                                             ! elements stored as the diagonal.
                                             ! u is vector stored and corresponds to
                                             ! the a priori covariance.
      real(real_kind1),intent(in) :: r;      !measurement variance
      real(real_kind1),intent(in) :: a(n);   !vector of measurement coefficients.
      real(real_kind1),intent(in) :: x(n);   !state vector
      real(real_kind1),intent(in) :: z;      !measurement

      real(real_kind1),intent(out) :: f(n);  !contains u**t*a.  f can overwrite a to save storage
      real(real_kind1),intent(out) :: g(n);  !intermediate vector.  After calling udmes2, g will be
                                             ! the unweighted kalman gain. The kalman
                                             ! gain k is equal to g/alpha
      real(real_kind1),intent(out) :: alpha; !innovations variance of the measurement residual
      real(real_kind1),intent(out) :: dz;    !measurement residual

      integer(4) j,k,l,jj,jjn,np1,np2,ntot,jm1;
      real(real_kind1) sum1,zero,one;
      !________________________________________________________

      zero=0.d0;
      one=1.d0;
      np1=n+1;
      np2=n+2;
      ntot=n*np1/2;
      dz =z-dot_product(a,x);

      jjn =ntot;
      alpha =r;
      do l=2,n;
        j =np2-l;
        jj =jjn-j;
        jm1 =j-1;
        sum1 =a(j) +dot_product(u(jj+1:jj+jm1),a(1:jm1));
        f(j) =sum1;
        g(j) =sum1*u(jjn);
        alpha =alpha +f(j)*g(j);
        jjn =jj;
      enddo;

!***     f=u**t*a and g=d*(u**t*a)
      f(1) =a(1);
      g(1) =u(1)*f(1);
      alpha =alpha +f(1)*g(1);

      return;
      end subroutine udmes1;
! ************************************************************
      subroutine udmes2 (g,u,x,  n,r,f,dz);
!***     this is the second half of esp subroutine udmeas
!***     eqn. numbers refer to bierman's 1975 cdc paper, pp. 337-346.
      use filterPrec, only: real_kind1;   !
      implicit  none;
      integer(4),intent(in) :: n;         !dimension of the state estimate. n > 1
      real(real_kind1),intent(in) :: r;            !measurement variance
      real(real_kind1),intent(in) :: f(n);         !contains u**t*a.
      real(real_kind1),intent(in) :: dz;           !measurement residual

      real(real_kind1),intent(inout) :: u(n*(n+1)/2);! updated, vector stored factors and estimate
                                          !upper triangular matrix, with d
                                          ! elements stored as the diagonal.
      real(real_kind1),intent(inout) :: x(n);      !updated state vector
      real(real_kind1),intent(inout) :: g(n);      !vector of unweighted kalman gains. the kalman
                                          ! gain k is equal to g/alpha

      integer(4) j,k,jm1,kj;
      real(real_kind1) sum1,beta,gamma,zero,one,temp,alpha,p,s,fnp1;
      !________________________________________________________

      zero =0.d0;
      one =1.d0;
      sum1 =r+g(1)*f(1);
      gamma =0.;                             !for r=0 case
      if (sum1 > zero) gamma =one/sum1;
      if (f(1) /= zero) u(1) =u(1)*r*gamma;  !d(1)

      kj=2;
      do j=2,n;
        beta =sum1;                          !beta=sum1(j-1)
        temp =g(j);
        sum1 =sum1+temp*f(j);
        p =-f(j)*gamma;                      !p=-f(j)*(1/sum1(j-1))  eqn(21)
        jm1 =j-1;
        do k=1,jm1;
          s =u(kj);
          !**  eqns(22) and (23)
          u(kj) =s+p*g(k);
          g(k) =g(k)+temp*s;
          kj =kj+1;
        enddo;

        if (temp /= zero) then;
          !**  r=0 case
          gamma =one/sum1;
          u(kj) =u(kj)*beta*gamma;           !d(j)   eqn(19)
        endif;
        kj =kj+1;
      enddo;

      alpha =sum1;
      fnp1 =dz*gamma;
      x(:) =x(:) +g(:)*fnp1;

      return;
      end subroutine udmes2;
