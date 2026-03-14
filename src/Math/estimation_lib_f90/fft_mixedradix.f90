!** fft_mixedRadix.f is a Fortran 90 updated version of Singleton's mixed-radix
!** Fast Fourier transform (FFT) as published in

!** 1) Singleton, R.C. (1969), “An algorithm for computing the mixed radix Fast Fourier
!**    Transform,” IEEE Trans. on Audio Electroacoustics, AU-17, 93-103,
!**    June 1969; reprinted in Digital Signal Processing, IEEE Press, New York, 1972
!** 2) Singleton, R.C. (1979), “Mixed radix Fast Fourier Transform,”
!**    Programs for Digital Signal Processing, IEEE Press, New York

!** This version was modified by B. Gibbs 10/24/2009 to 1) partially
!** break up code segments as subroutines, 2) minimize use of "go to"
!** statements, 3) use Fortran 90 features (intent, do while, etc.), and
!** locate common variables in module fftcom.

!**************************************
! The author grants the user a non-exclusive, worldwide, royalty-free copyright license to

! 1. reproduce and modify the software for your own purposes, and
! 2. to distribute the software provided that you give credit to the author,
!    do not make claims against the author or Wiley-Interscience,
!    and mark any modifications as your own.  If the software is incorporated in
!    a commercial product, you agree to defend and indemnify this author and
!    Wiley-Interscience against any losses, damages and costs arising from claims,
!    lawsuits and other legal actions brought by a third party.

! The software is provided on an “as is” basis, without warranties or conditions
! of any kind, either express or implied (including any warranties or conditions
! of title, non-infringement, merchantability or fitness for a particular purpose).
! The user is solely responsible for determining the appropriateness of using and
! distributing the program and assumes all risks associated with its exercise of
! rights under this agreement.
!**************************************

!**     Using Intel compiler
!**     Include in ifort.cfg:
!**   /check:bounds /traceback /include:\kf_book\distribSWm\library /link /nodefaultlib:msvcrt.lib
!**     \kf_book\distribSWm\library\libkf_bnd_ifort.lib

!**     To compile
!**   set ifortcfg=\kf_book\distribSWm\ifort.cfg
!**   ifort fft_mixedRadix.f

!************************************************************
      module fftcom;
!***   common variables for subroutine fft and realtr
      implicit none;
      integer(4),parameter :: maxf =23;   !maxf must be >= the maximum prime factor of n
      integer(4),parameter :: maxp =209;  !maxp must be > the number of prime factors of n
                                          !In addition, if the square-free portion k of n has two or
                                          !  more prime factors, then maxp must be >= k-1.
      integer(4),parameter :: mmax =11;   !maximum number of factors allowed
      integer(4),parameter :: real_kind =4;  !4 for single precision, 8 for double precision

      integer(4) nn;                      !#
      integer(4) jc;                      !#
      integer(4) nt;                      !initially inc*ntot, but decremented by kspnn at each stage
      integer(4) inc;                     !increment: +1 for forward transform, -1 for inverse
      real(real_kind) at(maxf),bt(maxf);  !array storage for maximum prime factor of 23
      real(real_kind) ck(maxf),sk(maxf);  !
      real(real_kind) sd,cd;              !sin & cos of ?

      end module fftcom;
!*****************************************************
      subroutine fft (ierr,  a,b,ntot,n,nspan,isn);
!  multivariate complex fourier transform, computed in place
!    using mixed-radix fast fourier transform algorithm.
!  by r. c. singleton, stanford research institute, oct. 1968

!** this version was modified by B. Gibbs 10/24/2009 to 1) partially
!** break up code segments as subroutines, 2) minimize use of "go to"
!** statements, 3) use Fortran 90 features (intent, do while, etc.), and
!** locate common variables in module fftcom.

!  arrays a and b originally hold the real and imaginary
!    components of the data, and return the real and
!    imaginary components of the resulting fourier coefficients.
!  multivariate data is indexed according to the fortran
!    array element successor function, without limit
!    on the number of implied multiple subscripts.
!    the subroutine is called once for each variate.
!    the calls for a multivariate transform may be in any order.
!  ntot is the total number of complex data values.
!  n is the dimension of the current variable.
!  nspan/n is the spacing of consecutive data values
!    while indexing the current variable.
!  the sign of isn determines the sign of the complex
!    exponential, and the magnitude of isn is normally one.
!  a tri-variate transform with a(n1,n2,n3), b(n1,n2,n3)
!    is computed by
!      call fft(a,b,n1*n2*n3,n1,n1,1)
!      call fft(a,b,n1*n2*n3,n2,n1*n2,1)
!      call fft(a,b,n1*n2*n3,n3,n1*n2*n3,1)
!  for a single-variate transform,
!    ntot = n = nspan = (number of complex data values), f.g.
!      call fft(a,b,n,n,n,1)
!  the data may alternatively be stored in a single complex
!    array a, then the magnitude of isn changed to two to
!    give the correct indexing increment and a(2) used to
!    pass the initial address for the sequence of imaginary
!    values, e.g.
!      call fft(a,a(2),ntot,n,nspan,2)

!  arrays at(maxf), ck(maxf), bt(maxf), sk(maxf), and np(maxp)
!    are used for temporary storage.  if the available storage
!    is insufficient, the program returns isn =0.

!   *** to convert program to double precision, change real_kind to 8 in fftmod

!    maxf must be >= the maximum prime factor of n.
!    maxp must be > the number of prime factors of n.
!    In addition, if the square-free portion k of n has two or
!    more prime factors, then maxp must be >= k-1.

      use fftcom;
      implicit none;
      integer(4),intent(in) :: ntot;      !total number of complex data points
      integer(4),intent(in) :: n;         !dimension of the current variable.
      integer(4),intent(in) :: nspan;     !nspan/n is the spacing of consecutive data values
                                          !  while indexing the current variable
      integer(4),intent(in) :: isn;       !The sign of isn determines the sign of the complex
                                          !   exponential, and the magnitude of isn is normally one.
                                          !If the given ntot cannot be factored, isn will be set to 0 on output.
      real(real_kind),intent(inout) :: a(ntot);   !real component of input data, or output FFT
      real(real_kind),intent(inout) :: b(ntot);   !imaginary component of input data, or output FFT
      integer(4),intent(out) :: ierr;     !error return: 0 =no error, 1=could not factor n

!**  array storage in nfac for a maximum of 11 factors of n.
!**  if n has more than one square-free factor, the product of the
!**    square-free factors must be <= 210

      integer(4) nfac(mmax),np(maxp),maxf1,ier;
      integer(4) i,j,k,jj,ks,kspan,jf,m,kt,kk;
      integer(4) kspnn,k1,k2,k3;
      real(real_kind) ak,bk,c1,c2,s1,s2,s72,c72,s120,rad,radf;
      !_____________________________________________________________

      ierr =0;
!     write(0,435);
  435 format(/' beginning execution of subroutine fft')
      if (n < 2) return;

      maxf1 =maxf;

      inc =isn;
      rad =6.28318530717958647692528676655900576d0;
      s72 =rad/5.0;
      c72 =cos(s72);
      s72 =sin(s72);
      s120 =0.86602540378443864676372317075293618d0;     !=sqrt(3)/2

      if (isn < 0) then;
        s72 =-s72;
        s120 =-s120;
        rad =-rad;
        inc =-inc;
      endif;

      nt =inc*ntot;
      ks =inc*nspan;
      kspan =ks;
      nn =nt-inc;
      jc =ks/n;
      radf =rad*float(jc)*0.5;
      i =0;
      jf =0;

!  determine the factors of n
      call get_factor (nfac,m,kt,  n);
      write (6,'("n =",i5,", factors =",11i4)') n,nfac(:m);

!  compute fourier transform (loop for processing factors)

  100 sd =radf/float(kspan);
      cd =2.0*sin(sd)**2;
      sd =sin(sd+sd);
      kk =1;
      i =i+1;

      if (nfac(i) == 2) then;
        !**  transform for factor of 2 (including rotation factor)
        call transform2 (a,b,kk,kspan,  ntot);   !  nn,jc,nt,inc,ntot,sd,cd)
        if (kk > kspan) go to 800;
        go to 100;

      else if (nfac(i) == 4) then;
        !**  transform for factor of 4
        call transform4 (a,b,kk,kspan,  isn,ntot);
        if (kspan == jc) go to 800;
        go to 100;
      endif;

  600 k =nfac(i);
      kspnn =kspan;
      kspan =kspan/k;

      if (k == 3) then;
        !**  transform for factor of 3 (optional code)
        call transform3 (a,b,kk,kspan,  ntot,s120);

      else if (k == 5) then;
        !**  transform for factor of 5 (optional code)
        call transform5 (a,b,kk,kspan,  ntot,s72,c72);
      else;
        !**  transform for odd factors
        call transformOdd (ier,  a,b,kk,kspan,jf,k,                     &
     &                      ntot,kspnn,rad);
        if (ier /= 0) go to 998;
      endif;

!**  multiply by rotation factor (except for factors of 2 and 4)
  700 if (i == m) go to 800;
      kk =jc+1;

      do while (kk <= jc+jc);
        c2=1.0-cd;
        s1=sd;

        do while (kk <= kspan);
          c1 =c2;
          s2 =s1;
          kk =kk+kspan;

          do;
            ak =a(kk);
            a(kk) =c2*ak-s2*b(kk);
            b(kk) =s2*ak+c2*b(kk);
            kk =kk+kspnn;
            if (kk <= nt) cycle;    !go to 730
            ak =s1*s2;
            s2 =s1*c2+c1*s2;
            c2 =c1*c2-ak;
            kk =kk-nt+kspan;
            if (kk > kspnn) exit;   !not: go to 730
          enddo;

          c2 =c1-(cd*c1+sd*s1);
          s1 =s1+(sd*c1-cd*s1);
          !**  the following three statements compensate for truncation
          !**  error. if rounded arithmetic used, they may be deleted.
          c1 =0.5/(c2**2+s1**2)+0.5;
          s1 =c1*s1;
          c2 =c1*c2;
          kk =kk-kspnn+jc;
        enddo;
!        if (kk <= kspan) go to 720

        kk=kk-kspan+jc+inc;
      enddo;
!      if (kk <= jc+jc) go to 710
      go to 100;

!**  permute the results to normal order---done in two stages
!**  permutation for square factors of n
  800 np(1)=ks;
      if (kt == 0) go to 890;

      k =kt+kt+1;
      if (m < k) k=k-1;
      j=1;
      np(k+1) =jc;
  810 np(j+1) =np(j)/nfac(j);
      np(k) =np(k+1)*nfac(j);
      j =j+1;
      k =k-1;
      if (j < k) go to 810;

      k3 =np(k+1);
      kspan =np(2);
      kk =jc+1;
      k2 =kspan+1;
      j =1;
      if (n /= ntot) go to 850;

!**  permutation for single-variate transform (optional code)
!      do
  820   ak =a(kk);
        a(kk) =a(k2);
        a(k2) =ak;
        bk =b(kk);
        b(kk) =b(k2);
        b(k2) =bk;
        kk =kk+inc;
        k2 =kspan+k2;
        if (k2 < ks) go to 820;

  830   k2 =k2-np(j);
        j =j+1;
        k2 =np(j+1)+k2;
        if (k2 > np(j)) go to 830;

        j=1;
  840   if (kk < k2) go to 820;

      kk =kk+inc;
      k2 =kspan+k2;
      if (k2 < ks) go to 840;
      if (kk < ks) go to 830;
      jc =k3;
      go to 890;

!**  permutation for multivariate transform
  850 k =kk+jc;

      do while (kk < k);
        ak =a(kk);
        a(kk) =a(k2);
        a(k2) =ak;
        bk =b(kk);
        b(kk) =b(k2);
        b(k2) =bk;
        kk =kk+inc;
        k2 =k2+inc;
      enddo;
!      if (kk < k) go to 860

      kk =kk+ks-jc;
      k2 =k2+ks-jc;
      if (kk < nt) go to 850;
      k2 =k2-nt+kspan;
      kk =kk-nt+jc;
      if (k2 < ks) go to 850;

  870 k2 =k2-np(j);
      j =j+1;
      k2 =np(j+1)+k2;
      if (k2 > np(j)) go to 870;
      j =1;
  880 if (kk < k2) go to 850;

      kk =kk+jc;
      k2 =kspan+k2;
      if (k2 < ks) go to 880;
      if (kk < ks) go to 870;
      jc=k3;

  890 if (2*kt+1 >= m) return;
      kspnn=np(kt+1);

!**  permutation for square-free factors of n
      j =m-kt;
      nfac(j+1) =1;
  900 nfac(j) =nfac(j)*nfac(j+1);
      j =j-1;
      if (j /= kt) go to 900;
      kt =kt+1;
      nn =nfac(kt)-1;
      if (nn > maxp) go to 998;
      jj =0;
      j =0;
      go to 906;

  902 jj =jj-k2;
      k2 =kk;
      k =k+1;
      kk =nfac(k);

  904 jj =kk+jj;
      if (jj >= k2) go to 902;

      np(j) =jj;
  906 k2 =nfac(kt);
      k =kt+1;
      kk =nfac(k);
      j =j+1;
      if (j <= nn) go to 904;

!  determine the permuation cycles of length greater than 1
      j=0;
      go to 914;

  910 k =kk;
      kk =np(k);
      np(k) =-kk;
      if (kk /= j) go to 910;
      k3 =kk;
  914 j =j+1;
      kk =np(j);
      if (kk < 0) go to 914;
      if (kk /= j) go to 910;
      np(j) =-j;
      if (j /= nn) go to 914;
      maxf1 =inc*maxf;
      !  reorder a and b, following the permutation cycles
      go to 950;

  924 j =j-1;
      if (np(j) < 0) go to 924;

      jj =jc;
  926 kspan =jj;
      if (jj > maxf1) kspan =maxf1;
      jj =jj-kspan;
      k =np(j);
      kk =jc*k+i+jj;
      k1 =kk+kspan;
      k2 =0;

      do while (k1 /= kk);
        k2 =k2+1;
        at(k2) =a(k1);
        bt(k2) =b(k1);
        k1 =k1-inc;
      enddo;
!      if (k1 /= kk) go to 928

  932 k1 =kk+kspan;
      k2 =k1-jc*(k+np(k));
      k =-np(k);

      do while (k1 /= kk);
        a(k1) =a(k2);
        b(k1) =b(k2);
        k1 =k1-inc;
        k2 =k2-inc;
      enddo;
!      if (k1 /= kk) go to 936

      kk =k2;
      if (k /= j) go to 932;

      k1 =kk+kspan;
      k2 =0;
      do while (k1 /= kk);
        k2 =k2+1;
        a(k1) =at(k2);
        b(k1) =bt(k2);
        k1 =k1-inc;
      enddo;
!      if (k1 /= kk) go to 940

      if (jj /= 0) go to 926;
      if (j /= 1) go to 924;

  950 j =k3+1;
      nt =nt-kspnn;
      i =nt-inc+1;
      if (nt >= 0) go to 924;
      return;

!  error finish, insufficient array storage
  998 ierr =1;
      write (6,999);
  999 format(/'array bounds exceeded within subroutine fft')
      return;
      end subroutine fft;
! *********************************************************************
      subroutine realtr (a,b,n,isn);
!  if isn=1, this subroutine completes the fourier transform
!    of 2*n real data values, where the original data values are
!    stored alternately in arrays a and b, and are first
!    transformed by a complex fourier transform of dimension n.
!    the cosine coefficients are in a(1),a(2),...a(n+1) and
!    the sine coefficients are in b(1),b(2),...b(n+1).
!    a typical calling sequence is
!     call fft(a,b,n,n,n,1)
!     call realtr(a,b,n,1)
!    the results should be multiplied by 0.5/n to give the
!    usual scaling of coefficients.
!  if isn=-1, the inverse transformation is done, the first step
!    in evaluating a real fourier series.
!    a typical calling sequence is
!      call realtr(a,b,n,-1)
!      call fft(a,b,n,n,n,-1)
!    the results should be multiplied by 0.5 to give the usual
!    scaling, and the time domain results alternate in arrays a
!    and b, i.e. a(1),b(1),a(2),b(2),...a(n),b(n).
!  the data may alternatively be stored in a single complex
!    array a, then the magnitude of isn changed to two to
!    give the correct indexing increment and a(2) used to
!    pass the initial address for the sequence of imaginary
!    values, e.g.
!      call fft(a,a(2),n,n,n,2)
!      call realtr(a,a(2),n,2)
!    in this case, the cosine and sine coefficients alternate in a.
!  by r.c. singleton, stanford research institute, oct. 1968

!   *** to convert program to double precision, change real_kind to 8

      use fftcom, only:real_kind;
      implicit none;
      integer(4),intent(in) :: n;
      integer(4),intent(in) :: isn;
      real(real_kind),intent(inout) :: a(2*n+2),b(2*n+2);

      integer(4) inc,nk,nh,j,k;
      real(real_kind) im,sd,cd,sn,cn,aa,ab,ba,bb,re;
      !____________________________________________

      inc =abs(isn);
      nk =n*inc+2;
      nh =nk/2;
      sd =3.14159265358979323846264338327950288d0/float(2*n);
      cd =2.0*sin(sd)**2;
      sd =sin(sd+sd);
      sn =0.0;
      if (isn < 0) then;
        cn =-1.0;
        sd =-sd;
      else;
        cn =1.0;
        a(nk-1) =a(1);
        b(nk-1) =b(1);
      endif;

      do j=1,nh,inc;
        k =nk-j;
        aa =a(j)+a(k);
        ab =a(j)-a(k);
        ba =b(j)+b(k);
        bb =b(j)-b(k);
        re =cn*ba+sn*ab;
        im =sn*ba-cn*ab;
        b(k) =im-bb;
        b(j) =im+bb;
        a(k) =aa-re;
        a(j) =aa+re;
        aa =cn-(cd*cn+sd*sn);
        sn =(sd*cn-cd*sn)+sn;
        !  the following three statements compensate for truncation
        !  error.  if rounded arithmetic is used, substitute
        !  cn=aa
        cn =0.5/(aa**2+sn**2)+0.5;
        sn =cn*sn;
        cn =cn*aa;
      enddo;

      return;
      end subroutine realtr;
!********************************************************
      subroutine get_factor (nfac,m,kt,  n);
!**      compute factors of n
      use fftcom, only:mmax;
      implicit none;
      integer(4),intent(in) :: n;            !number fo be factored
!      integer(4),intent(in) :: mmax          !maximum number of factors allowed
      integer(4),intent(out) :: nfac(mmax);  !factors of n
      integer(4),intent(out) :: m;           !number of factors found
      integer(4),intent(out) :: kt;          !remaining factors

      integer(4) k,j,jj;
      !___________________________________

      m =0;
      k =n;

      do while (k == (k/16)*16);
        m =m+1;
        if (m > mmax) go to 100;
        nfac(m) =4;
        k =k/16;
      enddo;

      j =3;
      jj =9;

      do while (jj <= k);
        do while (mod(k,jj) == 0);
          m =m+1;
          if (m > mmax) go to 100;
          nfac(m) =j;
          k =k/jj;
        enddo;
        j =j+2;
        jj =j**2;
      enddo;

      if (k <= 4) then;
        kt =m;
        nfac(m+1) =k;
        if (k /= 1) m =m+1;

      else;
        if (k == (k/4)*4) then;
          m =m+1;
          if (m > mmax) go to 100;
          nfac(m) =2;
          k =k/4;
        endif;

        kt =m;
        j =2;
        do while (j <= k);
          if (mod(k,j) == 0) then;
            m =m+1;
            if (m > mmax) go to 100;
            nfac(m) =j;
            k =k/j;
          endif;
          j =((j+1)/2)*2+1;
        enddo;
      endif;

      if (kt == 0) return;
      j =kt;

      do while (j /= 0);
        m =m+1;
        if (m > mmax) go to 100;
        nfac(m) =nfac(j);
        j =j-1;
      enddo;
      return;

  100 write (6,'("get_factor error: more than 11 factors found",        &
     &   " for n =",i5)') n;
      m =0;
      return;
      end subroutine get_factor;
!********************************************************************
      subroutine transform2 (a,b,kk,kspan,  ntot);
!  transform for factor of 2 (including rotation factor)
      use fftcom, only:real_kind,nn,jc,nt,inc,sd,cd;        !all input
      implicit none;
      integer(4),intent(in) :: ntot;
!      real(real_kind),intent(in) :: sd,cd
      integer(4),intent(inout) :: kk,kspan;
      real(real_kind),intent(inout) :: a(ntot),b(ntot);

      integer(4) k1,k2;
      real(real_kind) ak,bk,c1,s1;
      !___________________________________________________________

      kspan =kspan/2;
      k1 =kspan+2;

  210 k2 =kk+kspan;
      ak =a(k2);
      bk =b(k2);
      a(k2) =a(kk)-ak;
      b(k2) =b(kk)-bk;
      a(kk) =a(kk)+ak;
      b(kk) =b(kk)+bk;
      kk =k2+kspan;
      if (kk <= nn) go to 210;
      kk =kk-nn;
      if (kk <= jc) go to 210;
      if (kk > kspan) return; !go to 800
  220 c1 =1.0-cd;
      s1 =sd;

  230 k2 =kk+kspan;
      ak =a(kk)-a(k2);
      bk =b(kk)-b(k2);
      a(kk) =a(kk)+a(k2);
      b(kk) =b(kk)+b(k2);
      a(k2) =c1*ak-s1*bk;
      b(k2) =s1*ak+c1*bk;
      kk =k2+kspan;
      if (kk < nt) go to 230;
      k2 =kk-nt;
      c1 =-c1;
      kk =k1-k2;
      if (kk > k2) go to 230;
      ak =c1-(cd*c1+sd*s1);
      s1 =(sd*c1-cd*s1)+s1;
!  the following three statements compensate for truncation
!    error. if rounded arithmetic is used, substitue
!     c1=ak
      c1 =0.5/(ak**2+s1**2)+0.5;
      s1 =c1*s1;
      c1 =c1*ak;
      kk =kk+jc;
      if (kk < k2) go to 230;
      k1 =k1+inc+inc;
      kk =(k1-kspan)/2+jc;
      if (kk <= jc+jc) go to 220;

      return;
      end subroutine transform2;
!********************************************************************
      subroutine transform3 (a,b,kk,kspan,  ntot,s120);
!  transform for factor of 3 (optional code)
      use fftcom, only:real_kind,nn;
      implicit none;
      integer(4),intent(in) :: ntot;
      real(real_kind),intent(in) :: s120;
      integer(4),intent(inout) :: kk,kspan;
      real(real_kind),intent(inout) :: a(ntot),b(ntot);

      integer(4) k1,k2;
      real(real_kind) ak,bk,aj,bj;
      !___________________________________________________________

  320 k1 =kk+kspan;
      k2 =k1+kspan;
      ak =a(kk);
      bk =b(kk);
      aj =a(k1)+a(k2);
      bj =b(k1)+b(k2);
      a(kk) =ak+aj;
      b(kk) =bk+bj;
      ak =-0.5*aj+ak;
      bk =-0.5*bj+bk;
      aj =(a(k1)-a(k2))*s120;
      bj =(b(k1)-b(k2))*s120;
      a(k1) =ak-bj;
      b(k1) =bk+aj;
      a(k2) =ak+bj;
      b(k2) =bk-aj;
      kk =k2+kspan;
      if (kk < nn) go to 320;
      kk =kk-nn;
      if (kk <= kspan) go to 320;

      return;
      end subroutine transform3;
!******************************************************************
      subroutine transform4 (a,b,kk,kspan,  isn,ntot);
!  transform for factor of 4
      use fftcom, only:real_kind,jc,nt,inc,sd,cd;     !all input
      implicit none;
      integer(4),intent(in) :: ntot,isn;     !jc,nt,inc
!      real(real_kind),intent(in) :: sd,cd
      integer(4),intent(inout) :: kk,kspan;
      real(real_kind),intent(inout) :: a(ntot),b(ntot);

      integer(4) k1,k2;
      real(real_kind) c1,s1,kspnn,akp,akm,k3,ajp,ajm;
      real(real_kind) bkp,bkm,bjp,bjm,c2,s2,c3,s3;
      !___________________________________________________________

      kspnn=kspan;
      kspan=kspan/4;
  410 c1=1.0;
      s1=0.;

  420 k1=kk+kspan;
      k2=k1+kspan;
      k3=k2+kspan;
      akp=a(kk)+a(k2);
      akm=a(kk)-a(k2);
      ajp=a(k1)+a(k3);
      ajm=a(k1)-a(k3);
      a(kk)=akp+ajp;
      ajp=akp-ajp;
      bkp=b(kk)+b(k2);
      bkm=b(kk)-b(k2);
      bjp=b(k1)+b(k3);
      bjm=b(k1)-b(k3);
      b(kk)=bkp+bjp;
      bjp=bkp-bjp;
      if (isn < 0) go to 450;

      akp=akm-bjm;
      akm=akm+bjm;
      bkp=bkm+ajm;
      bkm=bkm-ajm;
      if (s1 == 0.0) go to 460;

  430 a(k1)=akp*c1-bkp*s1;
      b(k1)=akp*s1+bkp*c1;
      a(k2)=ajp*c2-bjp*s2;
      b(k2)=ajp*s2+bjp*c2;
      a(k3)=akm*c3-bkm*s3;
      b(k3)=akm*s3+bkm*c3;
      kk=k3+kspan;
      if (kk <= nt) go to 420;

  440 c2=c1-(cd*c1+sd*s1);
      s1=(sd*c1-cd*s1)+s1;
!  the following three statements compensate for truncation
!    error. if rounded arithmetic is used, substitute
!     c1=c2
      c1=0.5/(c2**2+s1**2)+0.5;
      s1=c1*s1;
      c1=c1*c2;
      c2=c1**2-s1**2;
      s2=2.0*c1*s1;
      c3=c2*c1-s2*s1;
      s3=c2*s1+s2*c1;
      kk=kk-nt+jc;
      if (kk <= kspan) go to 420;
      kk=kk-kspan+inc;
      if (kk <= jc) go to 410;
      if (kspan == jc) return;  !go to 800
      return;  !go to 100

  450 akp=akm+bjm;
      akm=akm-bjm;
      bkp=bkm-ajm;
      bkm=bkm+ajm;
      if (s1 /= 0.0) go to 430;

  460 a(k1)=akp;
      b(k1)=bkp;
      a(k2)=ajp;
      b(k2)=bjp;
      a(k3)=akm;
      b(k3)=bkm;
      kk=k3+kspan;
      if (kk <= nt) go to 420;
      go to 440;

      return;
      end subroutine transform4;
!****************************************************************
      subroutine transform5 (a,b,kk,kspan,  ntot,s72,c72);
!  transform for factor of 5 (optional code)
      use fftcom, only:real_kind,nn;
      implicit none;
      integer(4),intent(in) :: ntot;
      real(real_kind),intent(in) :: c72,s72;
      integer(4),intent(inout) :: kk,kspan;
      real(real_kind),intent(inout) :: a(ntot),b(ntot);

      integer(4) k1,k2,k3,k4;
      real(real_kind) ak,bk,aj,bj,c2,s2,aa;
      real(real_kind) akp,akm,bkp,bkm,ajp,ajm,bjp,bjm,bb;
      !___________________________________________________________

  510 c2 =c72**2-s72**2;
      s2 =2.0*c72*s72;

  520 k1 =kk+kspan;
      k2 =k1+kspan;
      k3 =k2+kspan;
      k4 =k3+kspan;
      akp =a(k1)+a(k4);
      akm =a(k1)-a(k4);
      bkp =b(k1)+b(k4);
      bkm =b(k1)-b(k4);
      ajp =a(k2)+a(k3);
      ajm =a(k2)-a(k3);
      bjp =b(k2)+b(k3);
      bjm =b(k2)-b(k3);
      aa =a(kk);
      bb =b(kk);
      a(kk) =aa+akp+ajp;
      b(kk) =bb+bkp+bjp;
      ak =akp*c72+ajp*c2+aa;
      bk =bkp*c72+bjp*c2+bb;
      aj =akm*s72+ajm*s2;
      bj =bkm*s72+bjm*s2;
      a(k1) =ak-bj;
      a(k4) =ak+bj;
      b(k1) =bk+aj;
      b(k4) =bk-aj;
      ak =akp*c2+ajp*c72+aa;
      bk =bkp*c2+bjp*c72+bb;
      aj =akm*s2-ajm*s72;
      bj =bkm*s2-bjm*s72;
      a(k2) =ak-bj;
      a(k3) =ak+bj;
      b(k2) =bk+aj;
      b(k3) =bk-aj;
      kk =k4+kspan;
      if (kk < nn) go to 520;
      kk =kk-nn;
      if (kk <= kspan) go to 520;

      return;
      end subroutine transform5;
!***********************************************************************
      subroutine transformOdd (ier,  a,b,kk,kspan,jf,k,                 &
     &                         ntot,kspnn,rad);
!  transform for odd factors
      use fftcom, only:real_kind,nn,maxf,at,bt,ck,sk;   !note: at,bt,ck,sk are inout
      implicit none;
      integer(4),intent(in) :: ntot,kspnn;
      real(real_kind),intent(in) :: rad;
      integer(4),intent(inout) :: kk,kspan,jf,k;
      real(real_kind),intent(inout) :: a(ntot),b(ntot);
!      real(real_kind),intent(inout) :: at(maxf),bt(maxf),ck(maxf),sk(maxf)
      integer(4),intent(out) :: ier;

      integer(4) k1,k2,jj,j;
      real(real_kind) ak,bk,aj,bj,s1,c1,aa,bb;
      !___________________________________________________________

      ier =0;

      if (k /= jf) then;
        jf=k;
        s1=rad/float(k);
        c1=cos(s1);
        s1=sin(s1);
        if (jf > maxf) then;
          ier =1;   !go to 998
          return;
        endif;

        ck(jf)=1.0;
        sk(jf)=0.0;
        j=1;
        do while (j < k);
          ck(j)=ck(k)*c1+sk(k)*s1;
          sk(j)=ck(k)*s1-sk(k)*c1;
          k=k-1;
          ck(k)=ck(j);
          sk(k)=-sk(j);
          j=j+1;
        enddo;
!        if (j < k) go to 630
      endif;

  640 k1=kk;
      k2=kk+kspnn;
      aa=a(kk);
      bb=b(kk);
      ak=aa;
      bk=bb;
      j=1;
      k1=k1+kspan;

      do while (k1 < k2);
        k2=k2-kspan;
        j=j+1;
        at(j)=a(k1)+a(k2);
        ak=at(j)+ak;
        bt(j)=b(k1)+b(k2);
        bk=bt(j)+bk;
        j=j+1;
        at(j)=a(k1)-a(k2);
        bt(j)=b(k1)-b(k2);
        k1=k1+kspan;
      enddo;
!      if (k1 < k2) go to 650

      a(kk)=ak;
      b(kk)=bk;
      k1=kk;
      k2=kk+kspnn;
      j=1;

      do while (j < k);
        k1=k1+kspan;
        k2=k2-kspan;
        jj=j;
        ak=aa;
        bk=bb;
        aj=0.0;
        bj=0.0;
        k=1;

        do while (k < jf);
          k=k+1;
          ak=at(k)*ck(jj)+ak;
          bk=bt(k)*ck(jj)+bk;
          k=k+1;
          aj=at(k)*sk(jj)+aj;
          bj=bt(k)*sk(jj)+bj;
          jj=jj+j;
          if (jj > jf) jj=jj-jf;
        enddo;
!       if (k < jf) go to 670

        k=jf-j;
        a(k1)=ak-bj;
        b(k1)=bk+aj;
        a(k2)=ak+bj;
        b(k2)=bk-aj;
        j=j+1;
      enddo;
!      if (j < k) go to 660

      kk =kk+kspnn;
      if (kk <= nn) go to 640;
      kk =kk-nn;
      if (kk <= kspan) go to 640;

      return;
      end subroutine transformOdd;

