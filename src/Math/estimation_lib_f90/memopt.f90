! *****************************************************************
      subroutine memopt (delf,pwrDB,                                    &
     &                   f,m,npmin,npmax,nlag,dt,iunit,junit,nplt);
!** memopt computes the power spectrum of real(4) data vector f
!** using the Burg maximum entropy method (MEM).  The Burg MEM computes autoregressive
!** coefficients that minimize the squared prediction error of forward and backward
!** prediction filters.  Memopt first fits models and computes the Akaike final
!** prediction error (FPE) for all model orders up to the maximum npmax.  It then
!** re-computes the model for the order with minimum FPE and uses the AR coefficients
!** to compute the power spectral density (PSD).
!** Reference: T.J. Ulrych and T.N. Bishop, "Maximum Entropy Spectral Analysis and
!** Autoregressive Decomposition", Rev. Geophysics and Space Phys, 13, 183-200 (Feb 1975).
!** Reprinted in Modern Spectrum Analysis, IEEE Press (1978)

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
      integer(4),intent(in) :: m;       !# data points
      integer(4),intent(in) :: npmin;   !min fit order
      integer(4),intent(in) :: npmax;   !max fit order
      integer(4),intent(in) :: nlag;    !# lags for autocorrelation fcn (nlag >= npmax,
                                        ! but nlag should be < m/2)
      integer(4),intent(in) :: iunit;   !plot unit
      integer(4),intent(in) :: junit;   !debug unit
      integer(4),intent(in) :: nplt;    !number of points to compute in pwrDB
      real(4),intent(in) :: f(m);       !time series
      real(4),intent(in) :: dt;         !time step between pts
      real(4),intent(out) :: pwrDB(nplt);!output power
      real(4),intent(out) :: delf;      !increment of frequency for pwrDB

      integer(4) i,n,nord;
      real(8) g(npmax+1);
      real(4) acov(nlag),fpe(npmax+1),pm,fmin,pwr,fm(m);
      !_____________________________________________

      if (nlag < npmax) stop 'memopt error: nlag < npmax';
      fm(:m) =f(:m);       !make local copy of f because mempr removes mean

!      write(0,420);
      write(6,420);
  420 format(/' beginning execution of Burg MEM')
      call mempr (g,fpe,acov,pwr,  fm,  m,npmax,nlag,junit);

      fmin=1.e28;
      do i=npmin,npmax;
        if (fpe(i) < fmin) then;
          nord =i;
          fmin =fpe(i);
        endif;
      enddo;

      write(6,1100) nord,fmin,10.d0**(fmin/2.d0);
1100  format(/"Burg MEM optimal order =",i4,                            &
     &        ", log10 (final prediction error) =",                     &
     &        g13.5,", normalized RMS =",g13.5)

      if ( nord < npmax ) then;
        call mempr (g,fpe,acov,pwr,  fm,  m,nord,nlag,junit);
      endif;
      call memplt (delf,pwrDB, nord,g,pwr,dt,m,iunit,junit,nplt);

      return;
      end subroutine memopt;
! *****************************************************************
      subroutine mempr (g,fpe,phi,pwr,  f,  m,lg,lext,junit);     !(m,f,lg,g,phi,lext,fpe,pwr,junit)
!***    mempr computes the Burg estimates of the AR prediction error filter
!**     coefficients and the Akaike final prediction error (FPE).
!**     Adapted from Ulyrich and Bishop, Maximum Entropy Spectral Analysis and
!**     Autoregressive Decomposition", Rev. Geophysics and Space Physics, V13, pp. 183-200, Feb 1975

!        note: mempr removes the mean from f

      implicit none;
      integer(4),intent(in) :: m;       !# data points
      integer(4),intent(in) :: lg;      !max fit order
      integer(4),intent(in) :: lext;    !# lags for autocorrelation fcn (lext >= lg,
                                        ! but lext should be < m/2)
      integer(4),intent(in) :: junit;   !debug unit
      real(4),intent(inout) :: f(m);    !time series: mempr removes the mean from f
      real(8),intent(out) :: g(lg);     !AR coefficients
      real(4),intent(out) :: fpe(lg);   !FPE vs order
      real(4),intent(out) :: phi(lext); !autocovariance vs lag
      real(4),intent(out) :: pwr;       !power of input white noise

      real(8) per(m),pef(m),h(lg),pm,dm,sumf,sn,sd,amean,ftemp;
      integer(4) i,j,k,n,nn,jj;
      !______________________________________

!     write(0,'(/" beginning execution of subroutine mempr"/)');
      if (lext < lg) then;
         write (0,'("mempr error lext < m")');
         stop;
      endif;

      amean =sum(f(:m))/m;
      sumf =sum(f(:m)**2);

      f(:) =f(:)-amean;
      pm =sum(f(:m)**2)/m;  !more accurate than sumf/m-amean**2
      dm =pm;
      phi(1) =pm;
      ftemp =(m+1)*pm/(m-1);
      fpe(1) =0.;
!     write(0,1010) lg, pm, sqrt(pm);
      write(6,1010) lg, pm, sqrt(pm);
 1010 format(/' mempr beginning lag loop; nn=2 to ',i4,                 &
     &        ', initial power =',g13.5,", RMS =",g13.5)

      do nn=2,lg;
!       write(0,'(/" nn=",i5/)') nn;
        n =nn-2;
        if (n == 0) then;
          pef(1:m) =0.d0;
          per(1:m) =0.d0;
        endif;

        sn=0.;
        sd=0.;
        jj =m-n-1;
        do j=1,jj;
          sn =sn-2.*(f(j+n+1)+pef(j))*(f(j)+per(j));
          sd =sd+(f(j+n+1)+pef(j))**2+(f(j)+per(j))**2;
        enddo;

        g(nn) =sn/sd;
        if (n > 0) then;
          do j=1,n;
            k =n-j+2;
            h(j+1) =g(j+1) +g(nn)*g(k);
          enddo;
          do j=1,n;
            g(j+1) =h(j+1);
          enddo;
          jj =jj-1;
        endif;

        do j=1,jj;
          per(j) =per(j) +g(nn)*(pef(j)+f(j+nn-1));
          pef(j) =pef(j+1) +g(nn)*(per(j+1)+f(j+1));
        enddo;
        sumf =0.d0;
        do j=2,nn;
          sumf =sumf -phi(nn+1-j)*g(j);
        enddo;
        phi(nn) =sumf;
        dm =(1.0d0-g(nn)**2)*dm;
        pm =dm;
        if (nn /= m) then;
          fpe(nn) =(m+nn)*pm/((m-nn)*ftemp);
          fpe(nn) =log10(fpe(nn));
        endif;
      enddo;

      g(1)=1.;
      do j=lg+1,lext;
        sumf=0.;
        do i=2,lg;
          sumf =sumf -phi(j+1-i)*g(i);
        enddo;
        phi(j) =sumf;
      enddo;

      write(6,1100) lg,pm,sqrt(pm);
 1100 format(/"max tested order (lg)=",i4,", fit residual power (pm) =",&
     &      g13.5,", RMS =",g13.5)
      write(6,1000) 'g',g;
!     write(6,1000) 'f',f
      write(6,1000) 'phi (covariance)',phi;
      write(6,1000) 'fpe (prediction error)',fpe;
 1000 format(/a," array",100(/(1x,10g13.5)))
 1001 format(/a," array",100(/(1x,10i10)))

!     write(junit+5,1100) lg, pm, sqrt(pm);
!     write(junit+5,'(/" phi =")');
!     write(junit+5,1000) 'phi (covariance)',phi;
!     write(junit+5,'(/" fpe =")');
!     write(junit+5,1000) 'fpe (prediction error)',fpe;
      pwr =pm;

      return;
      end subroutine mempr;
! *****************************************************************
      subroutine memplt (delf,pwrDB,                                    &
     &                   n,g,pm,dt,m,iunit,junit,nplt);
      implicit none;
      integer(4),intent(in) :: n;      !# data points
      integer(4),intent(in) :: m;      !not used
      integer(4),intent(in) :: iunit;  !plot unit
      integer(4),intent(in) :: junit;  !debug unit
      integer(4),intent(in) :: nplt;   !number of points to plot
      real(8),intent(in) :: g(n);  !time series
      real(4),intent(in) :: pm;    !white noise PSD ?
      real(4),intent(in) :: dt;    !time step between pts
      real(4),intent(out) :: pwrDB(nplt);
      real(4),intent(out) :: delf;  !increment of frequency for pwrDB

      integer(4) i,j,ix,m2,i1,imin;
!      complex(8) z,sumg,d1   !cmplx
      real(8) z(2),sumg(2),d1(2);  !cmplx
      real(4) dmin,dmax,pwr,range,scale,f;
      real(8) domg,dum,thet,psum;
      !______________________________________

      delf =1./(2*(nplt-1)*dt);        !Nyquist frequency/(nplt-1)
      domg =3.141592653589793d0/(nplt-1);
      dum =2.*pm*dt;
      dmax =-1e28;
      psum =0.d0;
      do i=1,nplt;
        thet =(i-1)*domg;
        z(:) =(/ cos(thet), -sin(thet) /);!complex
        sumg(:) = (/ g(1), 0.d0 /);
        d1(:) =(/ 1.d0, 0.d0 /);
!        z =cmplx(cos(thet),-sin(thet))
!        sumg =g(1)
!        d1 =1.
        do j=2,n;
          d1(:) =(/ d1(1)*z(1)-d1(2)*z(2), d1(1)*z(2)+d1(2)*z(1) /);
          sumg(:) =sumg(:) +g(j)*d1(:);
!          d1 =d1*z
!          sumg =sumg+cmplx(g(j))*d1
        enddo;
        pwr =dum/(sumg(1)**2+sumg(2)**2);
        psum =psum +pwr;
!        pwr =dum/(real(sumg)**2+aimag(sumg)**2)
        pwrDB(i) =10.*log10(pwr);
        dmax =max(dmax,pwrDB(i));
      enddo;

      write (6,'("Burg MEM summed power =",g13.5)') psum*delf;

      range=50.;
      i1 =dmax*0.1+0.999+1.e5;
      dmax =10*(i1-100000);
      dmin =dmax-range;
      scale =100./range;
!     write(iunit,200)dmax
  200 format('0 decibel plot. max db =',f12.0)
      imin=dmin;
!      if (iplot >= 0 ) call LPplot(0,iunit,0,imin,10,0,0.,0.)

      do i=1,nplt;
        f =(i-1)*delf;
!        ix =scale*(pwrDB(i+1)-dmin)+0.5
!        if (iplot >= 0 ) call LPplot(1,iunit,i,ix,0,0,f,pwrDB(i+1))
!        if (iplot <= 0 ) write(junit,321) f,pwrDB(i+1)
        write(junit,'(2f15.2)') f,pwrDB(i);
      enddo;

      return;
      end subroutine memplt;
! *****************************************************************
      subroutine LPplot (np,iout,ii,ix,io,ip,f,db);
!**     line printer plot
      real(4) aplt(101);
      integer(4) iaxs(7);
      data aplt /101*' ' /;
      data ab,a1,ax,ao,ap /' ','1','x','o','+' /;
      !__________________________________________

      if (np <= 0) then;
        do i=1,6;
          iaxs(i) =ix+(i-1)*io;
        enddo;
        write(iout,100) iaxs(1:6);
  100   format('1',20x,5(i5,15x),i5/                                    &
     &     ' ',23x,'1',5('---------1---------1'))
        return;
      endif;

      do i=1,6;
        aplt(20*(i-1)+1) =a1;
      enddo;
      ix =max(min(ix+1,101),1);
      aplt(ix) =ax;

      if (np > 1) then;
        io =max(min(io+1,101),1);
        aplt(io) =ao;
        if (np > 2) then;
          ip =max(min(ip+1,101),1);
          aplt(ip) =ap;
        endif;
      endif;

      write(iout,200)ii,f,db,aplt;
  200 format(1h ,i4,f10.4,f8.2,1x,101a1)
      aplt(ix) =ab;
      if (np >= 2) aplt(io) =ab;
      if (np >= 3) aplt(ip) =ab;

      return;
      end subroutine LPplot;
