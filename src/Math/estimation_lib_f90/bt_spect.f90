! **********************************************************************
      subroutine BT_spect (pwr,delf,n3,  c,nc,dt,nlag,iWindow);
!***   compute Blackman-Tukey correlogram power spectrum of data vector c using FFT with optional tapered windows
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

      integer(4),intent(in) :: nc;      !number of data points in c
      integer(4),intent(in) :: nlag;    !# of lags for autocorrelation function.  nlag <= nc/10 is recommended
      integer(4),intent(in) :: iWindow; !window function: 1=Bartlett, 2=Hamming (not reliable), >2 =Bartlett
      real(4),intent(in) :: c(nc);      !input data vector of length n
      real(4),intent(in) :: dt;         !time step between points (for sec
      real(4),intent(out) :: pwr(nc/2+1);!power spectrum of c (length nlag/2+1)
      real(4),intent(out) :: delf;      !frequency increment for d
      integer(4),intent(out) :: n3;     !# of samples in d

      integer(4) n,k,i,j,isn;
      real(4) freq,fm,dum,tspan,cm(nc),cmean,cvar,psum;
      real(4) a(2*nlag+1),b(2*nlag+1),c1(2*nlag+1),d1(2*nlag+1);
      real(4) sumx,scalew,sumw,r(-nlag:nlag),w(-nlag:nlag);
      !________________________________________

      fm =3.141593/nlag;
      cmean =sum(c(:nc))/nc;
      cm(:) =c(:)-cmean;         !remove mean
      cvar =sum(cm(:)**2)/nc;
      sumw =0.;
      r(:) =0.d0;

      !** compute unbiased autocorrelation function scaled by window function
      do k=0,nlag;
        sumx =0.d0;
        do j=1,nc-k;
          sumx =sumx +cm(j)*cm(j+k);
        enddo;

        if (iWindow == 2) then;
          !** Hamming window function (not reliable - scale is incorrect)
          w(k) =0.54+0.46*cos(fm*k);
        else;
          w(k) =1.0-abs(k/nlag);      !defaults to Bartlett window function
        endif;
        sumw =sumw +w(k);

        r(k) =sumx*w(k)/(nc-k);
        if (k > 0) then;
          r(-k) =r(k);
          sumw =sumw +w(k);
        endif;
      enddo;
!     write (6,'("r(0:nlag) =",10(/10g13.5))') r(0:nlag);

      scalew =(2*nlag+1)/sumw;
      r(:) =r(:)*scalew;
      n =2*nlag+1;
      b(:) =0.;
      a(1) =r(0);
      do k=1,nlag;
        a(k+1) =r(k);
        a(n+1-k) =r(-k);
      enddo;
!     write (6,'("a(:) =",10(/10g13.5))') a(:);

      !** now compute Fourier transform of autocorrelation function (could replace with fft)
      !** Note: dft does not scale c1,d1 by 1/n for isn > 0
      isn =+1;
      call dft (c1,d1,   isn,n,a,b);
!     write (6,'("c1(:) =",10(/10g13.5))') c1(:);
!     write (6,'("d1(:) =",10(/10g13.5))') d1(:);

      tspan =dt*n;
      delf =1./tspan;
      n3 =n/2+1;
      dum =2.*tspan;    !factor 2 accounts for negative frequencies
!     write (6,'(i3,", dt,tspan,dum =",4g13.5)')                        &
!    &                n,dt,tspan,dum;
      pwr(1) =0.5*dum*sqrt(c1(1)**2+d1(1)**2);
      psum =pwr(1);
      pwr(1) =10.*log10(0.5*dum*sqrt(c1(1)**2+d1(1)**2));
      do i=2,n3;
        freq =(i-1)*delf;
        pwr(i) =dum*sqrt(c1(i)**2+d1(i)**2);
        psum =psum +pwr(i);
        pwr(i) =10.*log10(pwr(i));
      enddo;
!     write (6,'("d(:) =",10(/10g13.5))') d(:);
      write (6,'("BT_spect data variance =",g13.5,                      &
     &       ", summed power =",g13.5)') cvar,psum*delf;  !checks

      return;
      end subroutine BT_spect;
