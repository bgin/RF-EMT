! **********************************************************************
      subroutine YW_spect (pwrDB,delf,g,  x,nc,dt,nlag,nplot);
!***    YW_spect computes the power spectrum of real(4) data vector x
!***    using the Yule-Walker maximum entropy method (MEM).  The Yule-Walker normal
!***    equations, containing sample autocorrelation coefficients, are solved
!***    to estimate coefficients of an autoregressive (AR) prediction error filter.
!***    YW_spect first fits models and computes the Akaike final prediction error
!***    (FPE) for all model orders up to the maximum npmax.  It then re-computes the
!***    model for the order with minimum FPE and uses the AR coefficients to compute
!***    the power spectral density (PSD).

!***    Subroutine ywpr is adapted from Ulyrich and Bishop, Maximum Entropy Spectral Analysis and
!**     Autoregressive Decomposition", Rev. Geophysics and Space Physics, V13, pp. 183-200, Feb 1975

!**     Author: B. Gibbs, 12/2009

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
      integer(4),intent(in) :: nlag;    !# of lags for autocorrelation function.  nlag <= nc/2
      integer(4),intent(in) :: nplot;   !# of pwrDB values to be computed
      real(4),intent(in) :: x(nc);      !input data vector of length n
      real(4),intent(in) :: dt;         !time step between points (sec)
      real(4),intent(out) :: pwrDB(nplot);!power spectrum (dB) of x
      real(4),intent(out) :: delf;      !frequency increment for pwrDB

      integer(4) n,k,i,j,isn,lext,lg,norder,ip;
      real(4) freq,dum,tspan,pmag,theta,psum;
      real(4) g(nc),phi(nc),fpe(nc),pm,fmax,xmean,xvar;
      real(8) :: pi =3.141592653589793d0;
      complex z,z1,p;
      !________________________________________

      xmean =sum(x(:))/nc;
      xvar =sum((x(:)-xmean)**2)/nc;

      !** evaluate Yule-Walker AR model and FPE versus order
      lg =nlag;     !max possible AR order
      lext =max(lg+1,nc/5);   !# of requested output autocorrelation coefficients in phi
      call ywpr (g,phi,fpe,pm,  lg,lext,nc,x);

      write (6,'("YW FPE =",5(/10g13.5))') fpe(:lg);
      write (6,'("YW g =",5(/10g13.5))') g(:lg);
      write (6,'("YW phi =",5(/10g13.5))') phi(:nlag+1);
      norder =minloc(fpe(2:lg),1) +1;
      write (6,'("YW-MEM optimal order =",i4)') norder;

      !** recompute Yule-Walker AR model using the minimum FPE model order
      call ywpr (g,phi,fpe,pm,  norder,lext,nc,x);

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
        p =g(1);
        do i=2,norder;
          z1 =z1*z;
          p =p +g(i)*z1;
        enddo;
        pmag =dum/(real(p)**2 +aimag(p)**2);
        psum =psum +pmag;
        pwrDB(ip) =10.*log10(pmag);
      enddo;
!     write (6,'("pwr(:) =",10(/10g13.5))') pwr(:);
      write (6,'("data variance =",g13.5,", YW summed power =",g13.5)') &
     &            xvar,psum*delf;  !checks

      return;
      end subroutine YW_spect;
! **********************************************************************
      subroutine ywpr (g,phi,fpe,pm,  lg,lext,n,x);
!***   compute Yule-Walker power spectrum of x using AIC FPE to
!***   select model order (see Ulrych and Bishop, 1975)
      implicit none;

      integer(4),intent(in) :: n;      !number of data points in x
      integer(4),intent(in) :: lg;     !# of lags used for computing direct autocorrelation
                                       ! function and fpe, i.e., max order of AR model
      integer(4),intent(in) :: lext;   !number of output correlation coeff (lext > lg)
      real(4),intent(in) :: x(n);      !input data vector of length n

      real(4),intent(out) :: g(lg);    !prediction error coefficients
      real(4),intent(out) :: fpe(lg);  !final prediction error versus AR order
      real(4),intent(out) :: phi(lext);!computed autocorrelation coefficients up to lag lext (lext > lg)
      real(4),intent(out) :: pm;       !updated variance

      integer(4) i,j,k,nn;
      real(8) h(n),dphi(n),suma,vp,dp,ftemp,phi8(lext),g8(lg),xmean;
      !__________________________________________________________

      if (lext < lg) stop 'ywpr error: lext < lg';

      xmean =sum(x(:n))/n;

      !** compute autocorrelation function phi
      do i=1,lg;
        j =min(n+i-1,n) -i+1;
        phi8(i) =dot_product(x(i:i+j-1)-xmean,x(1:j)-xmean) /n;
      enddo;
      dphi(:lg) =phi8(:lg);

      g8(1) =1.d0;
      g8(2) =-dphi(2)/dphi(1);
      ftemp =(real(n+1)/real(n-1)) *phi8(1);     !scaled variance of x
      fpe(1) =0.;

      !** recursively compute AR coefficients and fpe vs model order
      do nn=2,lg;
        vp =dot_product(g8(1:nn),dphi(1:nn));
        dp =0.;
        do i=1,nn;
          dp =dp +g8(nn+1-i)*dphi(i+1);   !dot_product(g(nn:1:-1),dphi(2:n+1))
        enddo;
        pm =vp;
        if (n /= nn) then;
          fpe(nn) =(real(n+nn)/real(n-nn)) *vp/ftemp;
          fpe(nn) =log10(fpe(nn));
        endif;

        if (nn == lg) exit;
        g8(nn+1) =-dp/vp;
        forall (i=2:nn) h(i) =g8(i) +g8(nn+1)*g8(nn+2-i);
        g8(2:nn) =h(2:nn);
      enddo;

      !** compute extended autocorrelation function
      do j=lg+1,lext;
        suma =0.;
        do i=2,lg;
          suma =suma -dphi(j+1-i)*g8(i);
        enddo;
        dphi(j) =suma;
        phi8(j) =suma;
      enddo;

      phi(:) =phi8(:);
      g(:) =g8(:);

      return;
      end subroutine ywpr;
