! **********************************************************************
      subroutine marple_spect (pwrDB,delf, nc,mmax,nplot,x,dt,tol1,tol2);
!***    marple_spect computes the power spectrum of real(4) data vector x
!***    using the Marple (1980) AR modeling method.
!**     "A New Autogregressive Spectral Analysis Algorithm", IEEE Trans on Acoustics,
!**     Speech and Signal Processing, ASSP-28, pp. 441-454 (Aug. 1980)

!**     Author: B. Gibbs, 1/2010

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

      integer(4),intent(in) :: nc;      !number of data points in c
      integer(4),intent(in) :: mmax;    !maximum AR order.  mmax <= nc/2
      integer(4),intent(in) :: nplot;   !# of pwrDB values to be computed
      real(4),intent(in) :: x(nc);      !input data vector of length n
      real(4),intent(in) :: dt;         !time step between points (sec)
      real(4),intent(in) :: tol1;       !0.01 to 0.001; use 0.01 for example
      real(4),intent(in) :: tol2;       !0.01 to 0.001; use 0.002 for example
      real(4),intent(out) :: pwrDB(nplot);!power spectrum (dB) of x
      real(4),intent(out) :: delf;      !frequency increment for pwrDB

      integer(4) i,norder,ip,status;
      real(4) dum,pmag,theta,psum,pvar,xm(nc);
      real(4) g(mmax),pm,fmax,xmean,xvar,e,e0;
      real(8) :: pi =3.141592653589793d0;
      complex z,z1,p;
      !________________________________________

      xmean =sum(x(:))/nc;
      xm(:) =x(:)-xmean;
      xvar =sum(xm(:)**2)/nc;

      !** evaluate Marple AR model
      call lstsqs (g,e,e0,norder,status,  nc,mmax,xm,tol1,tol2);
      pvar =e0/(2.d0*nc);
      pm =e/(2.d0*nc);
      write (6,'("marple_spect lstsqs status =",i3)') status;
      if (status < 4) then;
        write (6,'("marple_spect error: status =",i3)') status;
        stop 'marple_spect error';
      endif;

      write (6,'("Marple data variance =",g13.5,                        &
     &     ", prediction error variance =",g13.5)') pvar,e/(2.*nc);
      write (6,'("Marple g =",5(/10g13.5))') g(:norder);
      write (6,'("Marple optimal order =",i4)') norder;

      fmax =1.d0/(2.d0*dt);   !Nyquist frequency
      delf =fmax/nplot;
      dum =2.*pm*dt;       !factor 2 to account for negative frequencies
!     write (6,'(i3,", dt,tspan,dum =",4g13.5)')                        &
!    &                n,dt,tspan,dum;
      psum =0.;

      do ip=1,nplot;
        theta =(ip-1)*delf*2.*pi*dt;
        z =cmplx(cos(theta),-sin(theta));

        !** compute z-transform of
        z1 =cmplx(1.,0.);
        p =1.d0;  !g(1)       !note difference in AR definition from Ulrich and Bishop
        do i=1,norder;
          z1 =z1*z;
          p =p +g(i)*z1;
        enddo;
        pmag =dum/(real(p)**2 +aimag(p)**2);
        psum =psum +pmag;
        pwrDB(ip) =10.*log10(pmag);
      enddo;
!     write (6,'("pwr(:) =",10(/10g13.5))') pwr(:);
      write (6,'("data variance =",g13.5,                               &
     &        ", Marple summed power =",g13.5)') xvar,psum*delf;  !checks

      return;
      end subroutine marple_spect;
!***********************************************************************************
      subroutine lstsqs (a,e,e0,m,status,  n,mmax,x,tol1,tol2);
!**  lstsqs.f uses L. Marple's algorithms to compute parameters of an autoregressive model
!**  "A New Autogregressive Spectral Analysis Algorithm", IEEE Trans on Acoustics,
!**   Speech and Signal Processing, ASSP-28, pp. 441-454 (Aug. 1980)

      implicit none;
      integer(4),intent(in) :: n;   !number of data samples
      integer(4),intent(in) :: mmax;!max allowed model order
      real(4),intent(in) :: x(n);   !data samples
      real(4),intent(in) :: tol1;   !tolerance factor that stops recursion at order
                                    ! m when e(m)/e(0) < tol1
      real(4),intent(in) :: tol2;   !tolerance factor that stops recursion at order
                                    ! m when (e(m)-e(m-1)/e(m-1) < tol2

      integer(4),intent(out) :: m;  !model order at exit
      integer(4),intent(out) :: status;!error status: 1 =not converged at maximum order,
                                    ! 2 =, 3 = 4 =, 5 =OK
      real(4),intent(out) :: e0;     !total signal energy
      real(4),intent(out) :: e;      !prediction error energy (forward and backward) at order m
      real(4),intent(out) :: a(mmax);!array of m autoregressive coefficients

      integer(4) i,j,k,n1,nm,m1,k1,nk,m2,mk;
      real(8) c(mmax+1),d(mmax+1),r(mmax+1),f,b,h,s,u,v;
      real(8) save1,save2,save3,save4,delta,q1,q2,c1,c2,c3,c4,c5,c6;
      real(8) g,w,den,q3,q4,q5,q6,q7,alpha,eold,y1,y2,y3,y4;
      !___________________________________________________________________

!**     initialization
      e0 =2.*sum(x(:n)**2);
      q1 =1.d0/e0;
      q2 =q1*x(1);
      g =q1*x(1)**2;
      w =q1*x(n)**2;
      den =1.d0-g-w;
      q4 =1.d0/den;
      q5 =1.d0-g;
      q6 =1.d0-w;
      f =x(1);
      b =x(n);
      h =q2*x(n);
      s =q2*x(n);
      u =q1*x(n)**2;
      v =q2*x(1);
      e =e0*den;
      q1 =1.d0/e;
      c(1) =q1*x(1);
      d(1) =q1*x(n);
      save1 =0.;
      n1 =n+1;
      nm =n-1;
      do k=1,nm;
        save1 =save1 +x(k+1)*x(k);
      enddo;
      r(1) =2.d0*save1;
      a(1) =-q1*r(1);
      e =e*(1.d0-a(1)**2);

!**     loop for AR order
      m =1;
      do;
        !** prediction filter error update
        eold =e;
        m1 =m+1;
        f =x(m1);
        b =x(nm);
        do k=1,m;
          f =f+x(m1-k)*a(k);
          b =b+x(nm+k)*a(k);
        enddo;

        !** auxilliary vectors order update
        q1 =1.d0/e;
        q2 =q1*f;
        q3 =q1*b;
        do k=m,1,-1;
          k1 =k+1;
          c(k1) =c(k) +q2*a(k);
          d(k1) =d(k) +q3*a(k);
        enddo;
        c(1) =q2;
        d(1) =q3;

        !** scalar order update
        q7 =s**2;
        y1 =f**2;
        y2 =v**2;
        y3 =b**2;
        y4 =u**2;
        g =g +y1*q1 +q4*(y2*q6 +q7*q5 +2.d0*v*h*s);
        w =w +y3*q1 +q4*(y4*q5 +q7*q6 +2.d0*s*h*u);
        h =0.d0;
        s =0.d0;
        u =0.d0;
        v =0.d0;
        do k=0,m;
          k1 =k+1;
          nk =n-k;
          h =h +x(nm+k)*c(k1);
          s =s +x(nk)*c(k1);
          u =u +x(nk)*d(k1);
          v =v +x(k1)*c(k1);
        enddo;

        !** denominator update
        q5 =1.d0-g;
        q6 =1.d0-w;
        den =q5*q6 -h**2;
        if (den <= 0.d0) then;
          status =2;
          return;
        endif;

        !** time shift variables update
        q4 =1.d0/den;
        q1 =q1*q4;
        alpha =1.d0/(1.d0+(y1*q6+y3*q5+2.d0*(h*f*b))*q1);
        e =alpha*e;
        c1 =q4*(f*q6+b*h);
        c2 =q4*(b*q5+h*f);
        c3 =q4*(v*q6+h*s);
        c4 =q4*(s*q5+v*h);
        c5 =q4*(s*q6+h*u);
        c6 =q4*(u*q5+s*h);
        do k=1,m;
          k1 =k+1;
          a(k) =alpha*(a(k)+c1*c(k1)+c2*d(k1));
        enddo;
        m2 =m/2+1;
        do k=1,m2;
          mk =m+2-k;
          save1 =c(k);
          save2 =d(k);
          save3 =c(mk);
          save4 =d(mk);
          c(k) =c(k) +c3*save3+c4*save4;
          d(k) =d(k) +c5*save3+c6*save4;
          if (mk /= k) then;
            c(mk) =c(mk) +c3*save1 +c4*save2;
            d(mk) =d(mk) +c5*save1 +c6*save2;
          endif;
        enddo;

        if (m >= mmax) then;
          status =1;
          return;
        endif;

        !** order update
        m =m+1;
        nm =n-m;
        m1 =m-1;
        delta =0.d0;
        c1 =x(n1-m);
        c2 =x(m);
        do k=m1,1,-1;
          r(k+1) =r(k) -x(n1-k)*c1 -x(k)*c2;
          delta =delta +r(k+1)*a(k);
        enddo;
        save1 =0.d0;
        do k=1,nm;
          save1 =save1 +x(k+m)*x(k);
        enddo;
        r(1) =2.d0*save1;
        delta =delta+r(1);
        q2 =-delta/e;
        a(m) =q2;
        m2 =m/2;
        do k=1,m2;
          mk =m-k;
          save1 =a(k);
          a(k) =a(k) +q2*a(mk);
          if (k /= mk) a(mk) =a(mk) +q2*save1;
        enddo;

        y1 =q2**2;
        e =e*(1.d0-y1);

        if (y1 >= 1.d0) then;
          status =3;
          return;
        endif;
        if (e < e0*tol1) then;
          status =4;    !normal return
          return;
        endif;
        if ((eold-e) < eold*tol2) then;
          status =5;    !normal return
          return;
        endif;
      enddo;   !___ end order loop

      return;
      end subroutine lstsqs;
