! **********************************************************************
      subroutine spectd (c,d,nc,nplt,iunit,junit,dt,iplot,iWindow);

!*** subroutine spect computes the power spectrum of real(4) data vector c via the
!*** periodogram using a DFT with optional tapered windows.  Subroutine spectd is
!*** essentially the same as subroutine spect except that it uses the DFT (with
!*** any number of data points) rather than FFT  Variable iWindow specifies
!*** use of either rectangular (0), 10% cosine taper (1), Hamming (2), or Welch 50%
!*** overlapping segments using the Bartlett window (3).

!***  programmed by B. Gibbs, 12/2009

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
      integer(4),intent(in) :: nplt;    !number of points of d which are to be plotted
      integer(4),intent(in) :: iunit;   !plot unit #
      integer(4),intent(in) :: junit;   !unit # for printout
      integer(4),intent(in) :: iplot;   !0 = no plot, >0 =plot
      integer(4),intent(in) :: iWindow; !window function: 0=none, 1=10% cosine taper
      real(4),intent(in) :: c(nc);      !input data vector of length n
      real(4),intent(in) :: dt;         !time step between points (for sec
      real(4),intent(out) :: d(nc/2+1); !power spectrum of c (length n/2)

      integer(4) n,n2,n3,n5,nd,k,i,i1,i2,isn,imin,ix,np;
      real(4) fn,sumc,freq,pwr,f5,dum,tspan,dmin,dmax,range;
      real(4) a(nc),b(nc),cc(nc),dd(nc),f,fr,scalep,scalew;
      !________________________________________

      n =nc;
      isn =0;

      do while (isn == 0);
        write(6,'(/"start execution of subroutine spect: n=",i5)') n;
        n2 =n/2;
        n3 =n2+1;
        fn =float(n);
        n5 =n/10+1;
        f5 =3.141593/(0.1*fn);
        sumc =sum(c(1:n))/fn;
!       write(0,'(/"mean from spect :",g13.4)') sumc;
        k=1;
        do i=1,n2;
          a(i) =(c(k)-sumc)/fn;
          b(i) =(c(k+1)-sumc)/fn;
          k =k+2;
        enddo;
        nd =n+1-n5;

!.......multiply by 10% cosine taper window function
        if (iwindow == 1) then;
          do i=1,n5;
            dum =0.5*(1.-cos(f5*float(i-1)));
            i2 =i/2;
            if ( i /= 2*i2 ) then;
              a(i2+1)=dum*a(i2+1);
            else;
              b(i2)=dum*b(i2);
            endif;
            i1=n-i+1;
            i2=i1/2;
            if ( i1 /= 2*i2 ) then;
              a(i2+1)=dum*a(i2+1);
            else;
              b(i2)=dum*b(i2);
            endif;
          enddo;
          scalew =1./0.8625;
        else;
          scalew =1.;
        endif;

!       write(0,'(/"calling fft from spect")');
        isn =1;
        call fft (a,b,n2,n2,n2,isn);
        if (isn == 0 ) then;
          n =n-2;
          write (6,'("fft cannot handle n: try n=",i5)') n;
        endif;
      enddo;

      write(0,'(/" calling dft from spect")');
!      call realtr(a,b,n2,1)
      a(:n) =(c(:n)-sumc);
      b(:n) =0.;
      isn =1;        !forward transform
      call dft (cc,dd,  isn,n,a,b);
      write(0,'(/"dft completed execution")');

      tspan =fn*dt;
      dum =scalew*tspan;

      dmax =-1.e28;
      pwr =cc(1)**2*dum;
      d(1) =10.*log10(pwr);
      do i=2,n3;
        freq=(i-1)/tspan;
        pwr =(cc(i)**2 +dd(i)**2)*dum*4.;  !### 2 or 4 ?
        d(i) =10.*log10(pwr);
        if (d(i) > dmax ) dmax=d(i);
        if ( i > nplt ) cycle;
        write(6,300) i,freq,pwr,d(i);
  300   format(1x,i5,2g13.6,f14.6)
      enddo;

      if ( iplot > 0 ) then;
        range =50.;
        i1 =0.1*dmax+0.999+1.e5;
        dmax =10*(i1-100000);
        dmin =dmax-range;
        scalep =100./range;
        write(iunit,200) dmax;
  200   format('0  decibel plot : max db=',f12.6)
        imin =dmin;
!        call plot (0,iunit,0,imin,10,0,0.,0.)
        np =min0(nplt,n2);
        do i=1,np;
          f =i/tspan;
          ix =floor(scalep*(d(i+1)-dmin)+0.5);
!          call plot (1,iunit,i,ix,0,0,f,d(i+1))
        enddo;

      else;
        write(junit,310) 0.,d(1);
  310   format(2g13.5)
        do i=2,n2;
          fr =(i-1)/(fn*dt);
          if (i <= nplt) write(junit,310) fr,d(i);
        enddo;
        if (nplt >= n3) write(junit,310) n2/(fn*dt),d(n3);
      endif;

   88 continue;
!     write(0,'(/" exiting spect")');
      return;
      end subroutine spectd;
