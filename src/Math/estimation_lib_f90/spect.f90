! **********************************************************************
      subroutine spect (d,delf,n3,  c,nc,nplt,junit,dt,nseg,iWindow);

!*** subroutine spect computes the power spectrum of real(4) data vector c via the
!*** periodogram using the FFT with optional tapered windows. The FFT is Singleton's
!*** mixed-radix algorithm that works with most even numbers of data points nc.
!*** If the initial number nc is not acceptable for use in the FFT, spect will reduce nc
!*** by increments of 2 until an acceptable number is found.  Variable iWindow specifies
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
      integer(4),intent(in) :: junit;   !unit # for printout
      integer(4),intent(in) :: nseg;    !# of segments n which to average nc data points, i.e., ~nc/nseg per FFT
      integer(4),intent(in) :: iWindow; !window function: 0=none, 1=10% cosine taper, 2=Hamming,
                                        ! 3=Welch 50% overlapping segments using Bartlett window
      real(4),intent(in) :: c(nc);      !input data vector of length n
      real(4),intent(in) :: dt;         !time step between points (for sec
      real(4),intent(out) :: d(nplt);   !power spectrum of c (length nc/2)
      real(4),intent(out) :: delf;      !frequency increment for d
      integer(4),intent(out) :: n3;      !# of samples in d

      integer(4) n,n2,n5,k,i,j,isn,ierr,iseg,nseg1;
      logical(4) :: lfft =.true.;
      real(4) fn,sumc,freq,pwr(nc/2+1),f5,dum,tspan,dmax,cmean,cvar;
      real(4) a(nc/2+1),b(nc/2+1),c1(nc),d1(nc),fr,scalew,scalep,sumw;
      !________________________________________

      if (iWindow < 3) then;
        n =2*(nc/(2*nseg));   !make n even
      else; !50% overlapping windows (nseg must be odd)
        n =2*(nc/(2*(nseg+1)/2));
      endif;

      if (n <= 10) then;
        write (6,'("spect error: less than 10 pts per segment")');
        stop 'spect error: < 10 points/seg';
      endif;

      isn =1;      !forward transform

      !** check if n is OK, and if not reduce
      do j=1,10;
!       write(6,'(/"start execution of subroutine spect: n=",i5)') n;
        fn =real(n);
        n2 =n/2;
        sumc =sum(c(:n))/fn;
        k =1;
!       write(0,'(/"mean from spect :",g13.4)') sumc;
        do i=1,n2;
          a(i) =(c(k)-sumc)/fn;
          b(i) =(c(k+1)-sumc)/fn;
          k =k+2;
        enddo;
        call fft (ierr,  a,b,n2,n2,n2,isn);
        if (ierr == 0 ) exit;
        n =n-2;
        write (6,'("fft cannot handle n: try n=",i5)') n;
        if (j > 9) stop 'spect error: cannot factor n';
      enddo;

      if (iWindow < 3) then;
        nseg1 =nc/n;
      else;
        nseg1 =nc/(0.5*n)-1;
      endif;
      delf =1.d0/(n*dt);
!      n2 =n/2
      n3 =n2+1;
      pwr(:) =0.;
      k =1;

      do iseg =1,nseg1;
!       write (6,'("** iseg =",i4,", n =",i4,", k =",i4)') iseg,n,k;
        c1(:n) =c(k:k+n-1);
        sumc =sum(c1(:n))/fn;   !remove mean from each segment, or overall ??
        c1(:n) =c1(:n) -sumc;
!       write(0,'(/"mean from spect :",g13.4)') sumc;

        if (iwindow == 1) then;
          !** multiply by 10% cosine taper window function
          n5 =n/10+1;
!          nd =n+1-n5
          f5 =3.141593/(0.1*fn);
          sumw =n-2*n5;
          do i=1,n5;
            dum =0.5*(1.-cos(f5*float(i-1)));
            sumw =sumw +2.*dum**2;
            c1(i) =c1(i)*dum;
            c1(n+1-i) =c1(n+1-i)*dum;
          enddo;
          scalew =1./(sumw/n);   !1./0.8625
!          write (6,'("10% taper scale =",g13.5)') sumw/n  !0.89796 for n=500, 0.87296 if dum**2
          k =k+n;
        else if (iwindow == 2) then;
          !** multiply by Hamming window function #
          f5 =3.141593/n2;
          sumw =0.;
          do i=1,n2;
            dum =0.54+0.46*cos(f5*(n2-i));
            sumw =sumw +2.*dum**2;
            c1(i) =c1(i)*dum;
            c1(n+1-i) =c1(n+1-i)*dum;
          enddo;
          scalew =1./(sumw/n);   !###
          k =k+n;
        else if (iwindow == 3) then;
          !** Welch 50% overlapping windows with Bartlett window function
          !** multiply by 50% Bartlett window function
          sumw =0.;
          do i=1,n2;
            dum =abs(i/(0.5*fn));
            sumw =sumw +2.*dum**2;
            c1(i) =c1(i)*dum;
            c1(n+1-i) =c1(n+1-i)*dum;
          enddo;
          scalew =1./(sumw/n);   !###
          !** reset k for 50% overlap
          k =k+n2;
        else;
          !** rectangular window
          scalew =1.;
          k =k+n;
        endif;

        if (lfft) then;
          !** copy windowed c1 data to real/imaginary data for efficient FFT evaluation
          j=1;
          do i=1,n2;
            a(i) =c1(j)/(2.*fn);
            b(i) =c1(j+1)/(2.*fn);
            j =j+2;
          enddo;

!         write(0,'(/"calling fft from spect")');
          isn =1;      !forward transform
          call fft (ierr,  a,b,n2,n2,n2,isn);
          if (ierr /= 0 ) then;
            write (6,'("fft cannot handle n=",i5)') n;
            stop;
          endif;

!         write(0,'(/" calling realtr from spect")');
          isn =1;
          call realtr(a,b,n2,isn);
!         write(0,'(/" realtr completed execution")');
        else;
          d1(:) =0.;
          isn =1;
          call dft (a,b,  isn,n,c1,d1);
        endif;

        tspan =fn*dt;
        scalep =scalew*tspan;

        dmax =-1.e28;
!        pwr(1) =max(a(1)**2*scalep,1.e-36)
!        write (0,'("pwr0 =",g13.5)') pwr
!        d(1) =10.*log10(pwr(1))

        pwr(1) =max(pwr(1),1.e-36);
        do i=2,n3;
          freq =(i-1)/tspan;
          pwr(i) =pwr(i) +(a(i)**2+b(i)**2)*scalep;
        enddo;
      enddo;      !___ end segment loop

      pwr(1) =pwr(1)/nseg1;
      pwr(2:n3) =2.*pwr(2:n3)/nseg1;    !factor 2 account for negative frequencies
      cmean =sum(c(:nc))/nc;
      cvar =sum((c(:nc)-cmean)**2)/nc;
      write (6,'("FFT spect data variance =",g13.5,                     &
     &       ", summed power =",g13.5)') cvar,sum(pwr(:n3))*delf; !checks
      d(1) =-50.;

!      write (0,'("pwr0 =",g13.5)') pwr
      do i=2,n3;
        freq =(i-1)/tspan;
        d(i) =10.*log10(pwr(i));
        if (d(i) > dmax ) dmax =d(i);
        if ( i > nplt ) cycle;
        write(6,300) i,freq,pwr(i),d(i);
  300   format(1x,i5,2g13.6,f14.6)
      enddo;

      if (junit > 0) then;
        write (junit,1000);
 1000   format ('variables ="frequency", "power (db)"')
        write(junit,310) 0.,d(1);
  310   format(2g13.5)
        do i=1,n2;
          fr =(i-1)/(fn*dt);
          if (i > nplt) exit;
          write(junit,310) fr,d(i);
        enddo;
        if (nplt >= n3) write(junit,310) n2/(fn*dt),d(n3);
      endif;

   88 continue;
!     write(0,'(/" exiting spect")');
      return;
      end subroutine spect;
