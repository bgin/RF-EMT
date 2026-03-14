!**      subroutine cva performs canonical variate analysis using
!**      the method outlined by w. larimore (1983, 1987).
!**      Programmed by B. Gibbs, 9/30/2009

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
!**************************************************************

      module cvacom;
!**     common varibles used in the CVA computations
      implicit none;

      integer(4) :: ny;           !number of meas  at each sample
      integer(4) :: nu;           !number of controls at each sample
      integer(4) :: nyu;          !=ny+nu
      integer(4) :: m;            !total number of measurement times
      integer(4) :: nlag;         !number of lags to be used in computing sample covariances
      integer(4) :: nlago;        !optimal number of lags to be used for cva
      integer(4) :: nynlag;       !=ny*nlag
      integer(4) :: nyunlg;       !=nyu*nlag
      integer(4) :: nyuyunlg;     !=nyu*(nlag+1)
      integer(4) :: nrana;        !rank of
      integer(4) :: nrowpp;       !row dimension of sigpp, sigpf, jp
      integer(4) :: nrowff;       !row dimension of sigff, jf

      real(8),allocatable :: sigff(:,:);  !sample correlation of yu(:,it+nlag-i)*yu(:,it+nlag-j)^T for i,j=1:nlag
      real(8),allocatable :: sigpf(:,:);  !sample correlation of yu(:,it-i)*yu(:,it+nlag-j+nlag)^T
      real(8),allocatable :: sigpp(:,:);  !sample correlation of yu(:,it-i)*yu(:,it-j)^T
      real(8),allocatable :: sigp1p1(:,:);!sample correlation of yu(:,it-i)*yu(:,it-j)^T for i,j=0:nlag
      real(8),allocatable :: s(:);
!      real(8),allocatable :: jp(:,:),jf(:,:),covyu_yup(:,:)

      real(8) :: pi =3.141592653589793d0;

      end module cvacom;
!*********************************************************************
      subroutine cva (phi,h,q,rn,bw,gu,ac,n,                            &
     &             nyi,nyui,mi,minlag,maxlag,ndim,yu,eps,dt,directUfeed);

!**      subroutine cva performs canonical variate analysis using
!**      the method outlined by w. larimore (1983, 1987).
!**      programmed by b. gibbs, 9/30/2009
!**
!**      State space Model is
!**         y(i) =Ac*u(i) +H*x(i) +Bw*w(i) +v(i)
!**         x(i+1) =Phi*x(i) +Gu*u(i) +w(i)
!**      where
!**         u(i) is a known control
!**         w(i) is state noise
!**         v(i) is measurement noise
!**
!**      If directUfeed = .false., A is constrained to be 0
!**
!**     TTD:
!**       2) aic for ar model
!**       3) recursive state_space1 (?)
!**       4) biases

      use cvacom;
      implicit none;
      integer(4),intent(in) :: nyi;          !number of measurements at each time
      integer(4),intent(in) :: nyui;         !number of measurements+controls at each time (nyu =ny+nu)
      integer(4),intent(in) :: mi;           !total number of measurement times
      integer(4),intent(in) :: minlag;       !minimum number of lags used to compute cva model
      integer(4),intent(in) :: maxlag;       !maximum number of lags used to compute cva model. maxlag must
                                             ! be greater than the actual number of states in system.
      integer(4),intent(in) :: ndim;         !row dimension of phi and q in calling program
      real(8),intent(in) :: yu(nyui,mi);     !input data: measurements (ny) and control (nu) versus time (m)
      real(8),intent(in) :: eps;             !tolerance used to determine system rank (e.g. eps=1d-12)
      real(8),intent(in) :: dt;              !time step (sec)
      logical(4),intent(in) :: directUfeed;  !true allows direct feed of control u to measurement y (ac)
      real(8),intent(out) :: phi(ndim,ndim); !computed system transition matrix
      real(8),intent(out) :: h(nyi,ndim);    !computed system measurement matrix
      real(8),intent(out) :: q(ndim,ndim);   !computed state noise covariance matrix
      real(8),intent(out) :: rn(nyi,nyi);    !computed measurement noise covariance matrix
      real(8),intent(out) :: bw(nyi,ndim);   !computed direct feed of process noise q to meas y
      real(8),intent(out) :: gu(ndim,nyui-nyi);!computed feed of control to state x(i+1)
      real(8),intent(out) :: ac(ny,nyui-nyi);!computed direct feed of control to measurements (only if directUfeed =.true.)
      integer(4),intent(out) :: n;           !computed state order

      integer(4) :: smodel =1;   !1 is more efficient state space model computation,
                                 !2 is for testing and can later be used for recursive computation of Rn

      integer(4) ierr,nlag2,nnlag2;    !,nynlag,nyunlg,nrana,nyuyunlg
      integer(4) kmin,i,j,k,kk,it,i2,lag;    !,nrowpp,nrowff
      integer(4) nyunlg1,nynlag1,lagm;
      integer(4) :: nplot =100;
      real(8),allocatable :: jp(:,:),jf(:,:),covyu_yup(:,:),uc(:,:);
!      For nlag=50, nyu=ny=2: sigff(10000),sigpf(10000),sigpp(10000)
      character(40) plotFile;
      real(8) :: aic(ndim),aicm,logDet;
      !______________________________________________________

      ny =nyi;
      nyu =nyui;
      nu =nyu-ny;
      m =mi;           !total number of measurement times
      nlag =maxlag;
      nlag2 =nlag*2;
      nnlag2 =(nlag2*(nlag2+1))/2;
      nynlag =ny*nlag;
      nyunlg =nyu*nlag;
      nyuyunlg =nyu*(nlag+1);
      nrowpp =nyunlg;
      nrowff =nynlag;

      allocate (sigff(nrowff,nynlag),sigpf(nrowpp,nynlag),              &
     &          sigpp(nrowpp,nyunlg),jp(nrowpp,nyunlg),                 &
     &          jf(nrowff,nynlag),covyu_yup(nyu,nyu+nyunlg),            &
     &          s(nyunlg),sigp1p1(nyuyunlg,nyuyunlg),uc(nu,m),          &
     &          stat=ierr);
      if (ierr /= 0) stop 'cva error allocating sigff, sigpf, ...';
      uc(:nu,:) =yu(ny+1:ny+nu,:);

      !** call cvacorr0 computes the autocovariances needed for cva using the
      !** brute force approach.  cvacorr0 should only be used to check
      !** the efficient subroutine (cvacorr).  Output is sigff,sigpf,sigpp.
!      call cvacorr0 (yu)     !sigff,sigpf,sigpp,  yu,ny,nyu,m,nlag,nynlag,nyunlg)
      call cvacorr1 (yu);     !sigff,sigpf,sigpp,  yu,ny,nyu,m,nlag,nynlag,nyunlg)
      write (0,'("completed cvacorr")');

      !** calculate cov [(yi,ui)^T * (ui^T, p^T)]  and cov (yu * yu^T) (redundant with sigp1p1)
      covyu_yup(:,:) =0.d0;
      do it =nlag+1, m-nlag+1;
        forall (i=1:nyu)                                                &
     &          covyu_yup(:,i) =covyu_yup(:,i) +yu(:nyu,it)*yu(i,it);
        do lag=1,nlag;
          i2 =lag*nyu;
          do i=1,nyu;
            covyu_yup(:,i2+i) =covyu_yup(:,i2+i)                        &
     &                        +yu(:nyu,it)*yu(i,it-lag);       !cov (yi,ui)^T * p^T
          enddo;   !___ end i loop
        enddo;   !___ end lag loop
      enddo;   !___ end it loop
      covyu_yup(:,:) =covyu_yup(:,:)/(m-2*nlag+1);

!      if (smodel ==2)
      !** compute cov [(yu, p)^T, (yu, p)] by brute force (not efficient)
!      call get_sigp1p1 (sigp1p1,  yu)     !for testing alternate state space model and in arOrder

      !** the following is a more efficient method for computing sigp1p1
      sigp1p1(:nyu,:nyuyunlg) =covyu_yup(:nyu,:nyuyunlg);
      sigp1p1(:nyuyunlg,:nyu) =transpose(covyu_yup(:nyu,:nyuyunlg));
      sigp1p1(nyu+1:nyuyunlg,nyu+1:nyuyunlg) =sigpp(:,:);

      write (6,'(/"last columns of sigpf")');
      do i=1,ny;
        write (6,'(/i2,20(/10g13.5))') i,sigpf(:,nynlag-ny+i);   !checks with covyu_p
      enddo;

      !** compute optimal ar order using sigpp,sigff,sigpf,sigp1p1
      !** This selects optimal lags too low, but rest seems to work except for plotPsdY2
      lagm =nlag;    !default lags
      nlago =nlag;
      call arOrder (lagm,jp,jf,s,nrana,  minlag,eps);
!      lagm =20    !test - works ###
!      lagm =max(lagm,15)  test
      nlago =lagm;

      !** reset covariance sizes for optimal lags and re-compute cva
      nyunlg =nyu*lagm;
      nynlag =ny*lagm;

!***      canonical variate analysis:
!***      Use sigpf,sigff,sigpp to compute jp, jf
!***         output J = u**T*sigpp**(-1/2)  (stored in jp)
!***         output L = v**t*sigff**(-1/2)  (stored in jf)
      call canon (jp,jf,s,nrana,    nynlag,nyunlg,eps);

!     call prtmx (jf,nrowff,nynlag,nrana,"L");
!     call prtmx (jp,nrowpp,nrana,nyunlg,"J");
      write (0,'("completed canon")');

!***   evaluate AIC and find optimal model order using sigpp, jf, jp, yu, covyu_yup
      call maic (kmin,  jf,jp,yu,covyu_yup,directUfeed);
      write (0,'("completed maic")');

      !** generate plot file for cva model using sigpp
      plotFile ='ypsd_corr1.plt';
      call plotPsdY2 (nplot,dt,plotFile);    !ny,nyu,nyunlg,nlag,sigpp,)

!***   form state space arrays
      n =min(kmin,ndim);
      write (6,'("state space model ",i2)') smodel;

      if (smodel == 1) then;
        !** compute state space model from output J (jp), sigpp, and covyu_yup
        call state_space (phi,q,h,rn,bw,gu,ac,                          &
     &                    jp,covyu_yup,n,ndim,directUfeed);
      else;
        !** compute state space model from output J (jp), sigp1p1
        call state_space1 (phi,q,h,rn,bw,gu,ac,                         &
     &                     jp,n,ndim,directUfeed);
      endif;
      write (0,'("completed state_space")');

!** alternate computation of AIC using state space model
      aicm =1.d10;
      aic(:) =0.d0;
      kmin =1;

      do n =ny,min(nrana,ndim);
        write (6,'(/"*** testing model order ",i3)') n;

        if (smodel == 1) then;
          !** compute state space model from output J (jp), sigpp, and covyu_yup
          call state_space (phi,q,h,rn,bw,gu,ac,                        &
     &                      jp,covyu_yup,n,ndim,directUfeed);
        else;
          !** compute state space model from output J (jp), sigp1p1
          call state_space1 (phi,q,h,rn,bw,gu,ac,                       &
     &                       jp,n,ndim,directUfeed);
        endif;

        if (n == 4) then;       !### for true order
          call plotPsd2 (nplot,n,ndim,ny,nu,phi,q,h,rn,bw,gu,ac,        &
     &                uc,m,dt,'ypsd_cva4.plt');    !include controls
        endif;

        !** compute log-determinant of sig_yres_yres from a Kalman filter
        !** using the state space model
        call KF_logLike (logDet,  phi,q,h,rn,bw,gu,ac,yu,n,ndim);
!     &                            n,ndim,ny,nu,m,nlag)

        aic(n) =(m-nlag2+1)*(ny*(1+log(2.d0*pi)) +logdet)               &
     &            +2*(n*(2*ny+nu) +ny*(ny+1)/2);
        if (directUfeed) aic(n) =aic(n) +2*ny*nu;
        write (6,'("order ",i3,", KF logDet(Pz) =",g13.5,               &
     &           ", aic =",g13.5)') n,logDet,aic(n);
        if (aic(n) < aicm) then;
          aicm =aic(n);
          kmin =n;
        endif;
      enddo;

      write (6,'("aic =",2(/10g13.5))') aic(:min(nrana,ndim));
      write (6,'("CVA AIC optimal order =",i3)') kmin;

      n =kmin;
      if (smodel == 1) then;
        !** compute state space model from output J (jp), sigpp, and covyu_yup
        call state_space (phi,q,h,rn,bw,gu,ac,                          &
     &                    jp,covyu_yup,n,ndim,directUfeed);
      else;
        !** compute state space model from output J (jp), sigp1p1
        call state_space1 (phi,q,h,rn,bw,gu,ac,                         &
     &                     jp,n,ndim,directUfeed);
      endif;

      !** generate plot file of power spectrum using state-space model and
      !** history of controls
      call plotPsd2 (nplot,n,ndim,ny,nu,phi,q,h,rn,bw,gu,ac,            &
     &                uc,m,dt,'ypsd_cva.plt');

      deallocate (sigff,sigpf,sigpp,jp,jf,covyu_yup,s,sigp1p1,uc);
      return;
      end subroutine cva;
! ******************************************************************
      subroutine cvacorr0 (yu);   !sigff,sigpf,sigpp,  ny,nyu,m,nlag,nynlag,nyunlg)
!        cvacorr0 computes the autocovariances needed for cva using the
!        brute force approach. This subroutine should only be used to check the efficient
!        subroutine (cvacorr).  Output is sigff,sigpp, sigpf
      use cvacom, only: nyunlg,nynlag,ny,nu,nyu,m,nlag,sigpp,sigff,sigpf;
      implicit none;
!      integer(4),intent(in) :: ny,nyu        !number of meas and controls at each sample
!      integer(4),intent(in) :: m             !total number of measurement times
!      integer(4),intent(in) :: nlag          !number of lags to be used
!      integer(4),intent(in) :: nynlag        !=ny*nlag
!      integer(4),intent(in) :: nyunlg        !=nyu*nlag
      real(8),intent(in) :: yu(nyu,m);       !measurement+control data at each sample
!      real(8),intent(out) :: sigff(nynlag,nynlag)   !sample correlation of yu(:,it+nlag-i)*yu(:,it+nlag-j)^T
!      real(8),intent(out) :: sigpf(nyunlg,nynlag)   !sample correlation of yu(:,it-i)*yu(:,it+nlag-j+nlag)^T
!      real(8),intent(out) :: sigpp(nyunlg,nyunlg)   !sample correlation of yu(:,it-i)*yu(:,it-j)^T

      integer(4) i,j,it,jp,jf,jc,ip,if,ia,ib,nlag2;
      real(8) cons;
      logical(4) :: lprint =.false.;
      !__________________________________________________

!      if (nyu /= ny) stop 'cvacorr0 error: nyu /= ny'
      sigff(:,:) =0.d0;
      sigpf(:,:) =0.d0;
      sigpp(:,:) =0.d0;
      nlag2 =2*nlag;
      cons =1.d0/(m-nlag2+1);

      do j=1,nlag;
        jp =(j-1)*nyu;
        jf =(j-1)*ny;
        do i=1,nlag;
          ip =(i-1)*nyu;
          if =(i-1)*ny;
          ia =ip+1;
          ib =ip+nyu;
          do it =nlag+1, m-nlag+1;
            if (i <= j) then;
              do jc=1,ny;
                sigff(if+1:if+ny,jf+jc) =sigff(if+1:if+ny,jf+jc)+       &
     &                               yu(1:ny,it+nlag-i)*yu(jc,it+nlag-j);
              enddo;
              do jc=1,nyu;
                sigpp(ia:ib,jp+jc) =sigpp(ia:ib,jp+jc)+                 &
     &                               yu(1:nyu,it-i)*yu(jc,it-j);
              enddo;   !___ end jc loop
            endif;
            do jc=1,ny;
              sigpf(ia:ib,jf+jc) =sigpf(ia:ib,jf+jc)+                   &
     &                             yu(1:nyu,it-i)*yu(jc,it+nlag-j);
            enddo;
          enddo;   !___ end it loop
        enddo;    !___ end i loop
      enddo;   !___ end j loop

      do j=2,nynlag;
        do i=1,j-1;
          sigff(j,i) =sigff(i,j);
        enddo;
      enddo;
      do j=2,nyunlg;
        do i=1,j-1;
          sigpp(j,i) =sigpp(i,j);
        enddo;
      enddo;

      sigff(:,:) =cons*sigff(:,:);
      sigpf(:,:) =cons*sigpf(:,:);
      sigpp(:,:) =cons*sigpp(:,:);

      if (lprint) then;
        call prtmx(sigff,nynlag,nynlag,nynlag,'sigff0');
        call prtmx(sigpf,nyunlg,nyunlg,nynlag,'sigpf0');
        call prtmx(sigpp,nyunlg,nyunlg,nyunlg,'sigpp0');
      endif;

      return;
      end subroutine cvacorr0;
! *******************************************************************
      subroutine cvacorr1 (yu);     !sigff,sigpf,sigpp,  ny,nyu,m,nlag,nynlag,nyunlg)
!**      cvacorr computes the sample autocorrelation for cva using efficient methods.
!**      Output is sigff,sigpp, sigpf
      use cvacom, only: nyunlg,nynlag,ny,nu,nyu,m,nlag,sigpp,sigff,sigpf;
      implicit none;
      real(8),intent(in) :: yu(nyu,m);       !measurement+control data at each sample

      integer(4) i,j,k,lag1,lag,ic,k1,k2,jf,jp,if,ke;
      integer(4) nlag2,nnlag2,ierr,it,i1,i2,ip,ijd,j1,j2;
      logical(4) :: lprint =.false.;
      real(8) cons;
      real(8),allocatable :: b(:,:,:,:),r(:,:,:);
      !__________________________________________________

      nlag2 =2*nlag;
      nnlag2 =(nlag2*(nlag2+1))/2;
      allocate (b(nyu,nyu,nnlag2,nlag2),r(nyu,nyu,-nlag2+1:nlag2-1),    &
     &          stat=ierr);
      if (ierr /= 0) stop 'cvacorr error allocating b,r';

      r(:,:,:) =0d0;
      b(:,:,:,:) =0.d0;

!**      computation of autocorrelation r
      do lag=0,nlag2-1;
        do it=1,m-lag;    !time loop
          do i=1,nyu;     !measurement/control loop
            r(:nyu,i,lag) =r(:nyu,i,lag) +yu(:nyu,it+lag)*yu(i,it);
          enddo;
        enddo;
      enddo;
      do lag=1,nlag2-1;
        r(:,:,-lag) =transpose(r(:,:,lag));
      enddo;

!**      compute beginning effect - b
      do ic =0,nlag2-2;
        do j =nlag2-1,1+ic,-1;
          i =j-ic;
          if (i < 1) exit;
          do i1=1,nyu;
            b(1:nyu,i1,i,j) =b(1:nyu,i1,i+1,j+1)                        &
     &                       +yu(1:nyu,nlag2-i)*yu(i1,nlag2-j);
          enddo;
        enddo;
      enddo;
      do j=2,nlag2;
        do i=1,j-1;
          b(:,:,j,i) =transpose(b(:,:,i,j));
        enddo;
      enddo;

!**      now form the sample covariances including beginning effects
      cons =1.d0/float(m-nlag2+1);
      nu =nyu-ny;

      do j=1,nlag;
        jf =(j-1)*ny;
        jp =(j-1)*nyu;
        do i=1,nlag;
          if =(i-1)*ny;
          ip =(i-1)*nyu;
          ijd =(j-i);
          do j2=1,ny;
            do i2=1,ny;
              sigff(if+i2,jf+j2) =cons*(r(i2,j2,ijd) -b(i2,j2,i,j));
            enddo;
          enddo;

          ijd =(j-i-nlag);
          do j2=1,ny;
            do i2=1,nyu;
              sigpf(ip+i2,jf+j2) =cons*(r(i2,j2,ijd) -b(i2,j2,i+nlag,j));
            enddo;
          enddo;

          ijd =((j+nlag)-(i+nlag));
          k =k+nlag;
          do j2=1,nyu;
            do i2=1,nyu;
              sigpp(ip+i2,jp+j2) =cons*(r(i2,j2,ijd)                    &
     &                              -b(i2,j2,i+nlag,j+nlag));
            enddo;
          enddo;
        enddo;
      enddo;

!***      now compute end effects (e) re-using b array.
!***      first zero first row of b
      ke =1;
      do lag1=1,nlag2;
        b(:,:,1,lag1) =0.d0;
        ke =ke+lag1;
      enddo;
!      b(:,:,:,:) =0.d0

      do ic=0,nlag2-2;
        do i=2,nlag2-ic;
          j =i+ic;
          k1 =((j-1)*j)/2+i;
          k2 =k1-j;
          do i1=1,nyu;
            do j1=1,nyu;
              b(j1,i1,i,j) =b(j1,i1,i-1,j-1)                            &
     &                      +yu(j1,m+2-i)*yu(i1,m+2-j);
            enddo;
          enddo;
        enddo;
      enddo;
      do j=2,nlag2;
        do i=1,j-1;
          b(:,:,j,i) =transpose(b(:,:,i,j));
        enddo;
      enddo;

!***       add end effects to sig arrays
      do j=1,nlag;
        jf =(j-1)*ny;
        jp =(j-1)*nyu;
        do i=1,nlag;
          if =(i-1)*ny;
          ip =(i-1)*nyu;
          do j2=1,ny;
            do i2=1,ny;
              sigff(if+i2,jf+j2) =sigff(if+i2,jf+j2) -cons*b(i2,j2,i,j);
            enddo;
          enddo;

          do j2=1,ny;
            do i2=1,nyu;
              sigpf(ip+i2,jf+j2) =sigpf(ip+i2,jf+j2)                    &
     &                            -cons*b(i2,j2,i+nlag,j);
            enddo;
          enddo;

          do j2=1,nyu;
            do i2=1,nyu;
              sigpp(ip+i2,jp+j2) =sigpp(ip+i2,jp+j2)                    &
     &                           -cons*b(i2,j2,i+nlag,j+nlag);
            enddo;
          enddo;
        enddo;
      enddo;

      deallocate (b,r);
      if (lprint) then;
        call prtmx(sigff,nynlag,nynlag,nynlag,'sigff');
        call prtmx(sigpf,nyunlg,nyunlg,nynlag,'sigpf');
        call prtmx(sigpp,nyunlg,nyunlg,nyunlg,'sigpp');
      endif;

      return;
      end subroutine cvacorr1;
! ******************************************************************
      subroutine get_sigp1p1 (sigp1p1,  yu);    !,nyu,m,nlag,nyunlg,nyuyunlg)
!        get_sigp1p1 computes sample cov ((y,u,p) (y,u,p)^T) using the
!        brute force approach. This subroutine should only be used for testing more efficient methods
      use cvacom, only: nyunlg,nyuyunlg,nyu,m,nlag;
      implicit none;
      real(8),intent(in) :: yu(nyu,m);       !measurement+control data at each sample
      real(8),intent(out) :: sigp1p1(nyuyunlg,nyuyunlg);  !sample correlation of yu(:,it-i)*yu(:,it-j)^T for i,j over 0 to nlag

      integer(4) i,j,it,jp,jc,ip,ia,ib,nlag2,ntot;
      real(8) cons;
      !__________________________________________________

      sigp1p1(:,:) =0.d0;
      nlag2 =2*nlag;
      cons =1.d0/(m-nlag2+1);
      ntot =nyuyunlg;

      do j=0,nlag;
        jp =j*nyu;
        do i=0,nlag;
          ip =i*nyu;
          ia =ip+1;
          ib =ip+nyu;
          do it =nlag+1, m-nlag+1;
            if (i <= j) then;
              do jc=1,nyu;
                sigp1p1(ia:ib,jp+jc) =sigp1p1(ia:ib,jp+jc)+             &
     &                               yu(1:nyu,it-i)*yu(jc,it-j);
              enddo;   !___ end jc loop
            endif;
          enddo;   !___ end it loop
        enddo;    !___ end i loop
      enddo;   !___ end j loop

      do j=2,ntot;
        do i=1,j-1;
          sigp1p1(j,i) =sigp1p1(i,j);
        enddo;
      enddo;
      sigp1p1(:,:) =cons*sigp1p1(:,:);

!     call prtmx(sigp1p1,nyunlg+nyu,nyunlg+nyu,nyunlg+nyu,'sigp1p1');

      return;
      end subroutine get_sigp1p1;
!*************************************************************************
      subroutine arOrder (lagm,jp,jf,s,nrana,  minlag,eps);
!**      subroutine arOrder computes the optimal AR order using the AIC
      use cvacom, only: nyunlg,nynlag,nyuyunlg,nrowpp,nrowff,           &
     &                  ny,nu,nyu,m,nlag,sigp1p1,pi;
      implicit none;
      integer(4),intent(in) :: minlag;             !minimum # lags to be tested
      real(8),intent(in) :: eps;                   !tolerance for testing singular values for rank (e.g., 1d-12)
      real(8),intent(out) :: jp(nrowpp,nyunlg);    !mapping from canonical states to past
      real(8),intent(out) :: jf(nrowff,nynlag);    !mapping from canonical states to future
      real(8),intent(out) :: s(nyunlg);            !singular values of sigpp^(-1/2)*sigpf*sigff**(-1/2)
      integer(4),intent(out) :: nrana;             !rank of jf
      integer(4),intent(out) :: lagm;              !optimal lags

      integer(4) lag,i,j,k,kk,ierr;
      integer(4) nyunlg1,nynlag1;
      real(8) logdet,aicLg(nlag),aicm,d1(ny,nyuyunlg),ypre(ny,nyunlg);
      real(8) sigres(ny*(ny+1)/2),eps1;
!      real(8) :: pi =3.141592653589793d0
      !________________________________________________

      write (6,'(//"**********************************"/                &
     &           "*** start AIC testing of AR order")');

      !** test for optimal # lags
      aicm =1.d20;
      lagm =0;
      aicLg(:) =0.;

      do lag=minlag,nlag;
        nyunlg1 =nyu*lag;
        nynlag1 =ny*lag;
        !**  canonical variate analysis: use sigpf,sigff,sigpp to compute jp, jf
        call canon (jp,jf,s,nrana,  nynlag1,nyunlg1,eps);
        kk =min(nynlag1,nrana);
        if (nynlag1 > nrana)                                            &
     &      write (6,'("for lag ",i2,", nynlag1, nrana =",2i4)')        &
     &            lag,nynlag1,nrana;
        ypre(:ny,:nyunlg1) =                                            &
     &            matmul(jf(nynlag1-ny+1:nynlag1,:kk),jp(:kk,:nyunlg1));
        d1(:ny,:nyu+nyunlg1) =matmul(ypre(:ny,:nyunlg1),                &
     &                          sigp1p1(nyu+1:nyu+nyunlg1,:nyu+nyunlg1));
        k =0;
        do i=1,ny;
          do j=1,i;
            k =k+1;
            sigres(k) =sigp1p1(i,j) -(d1(i,j)+d1(j,i))                  &
     &         +dot_product(d1(i,nyu+1:nyu+nyunlg1),ypre(j,:nyunlg1));
          enddo;
        enddo;

        !** get cholesky factor of sigres so that log determinant can be computed
        eps1 =1.d-4;
        call cfactor (ierr,  sigres,ny,eps1);
        if (ierr < 0) then;
          write (6,'("cva error for lag",i3,": sigres is singular")')lag;
          cycle;
        endif;
        logdet =0.d0;
        do j=1,ny;
          kk =(j*(j+1))/2;
          logdet =logdet +2.d0*log(sigres(kk));
        enddo;
        write (6,'("For lag ",i3,", logdet(sigres) =",g13.5)')lag,logdet;

        aicLg(lag) =(m-2*nlag+1)*(ny*(1+log(2.d0*pi)) +logdet)          &
     &            +2*(nrana*(ny+nu) +ny*(ny+1)/2);        !this probably isn't correct ###
!     &            +2*(nrana*(2*ny+nu) +ny*(ny+1)/2)         !this probably isn't correct ###
!        if (directUfeed) aicLg(lag) =aicLg(lag) +2*ny*nu
        write (6,'("for lag =",i3,", aic =",g13.5)') lag,aicLg(lag);
        if (aicLg(lag) < aicm) then;
          aicm =aicLg(lag);
          lagm =lag;
        else;
          if (lagm > 0 .and. lag > lagm+20) exit;
        endif;
      enddo;

      write (6,'("ar-Order aicLG =",10(/10g13.5))') aicLg(:);
      write (6,'("optimal AR lag =",i3)') lagm;

      return;
      end subroutine arOrder;
! ******************************************************************
      subroutine canon (jp,jf,s,nrana,  nynlag,nyunlg,eps);
!**      subroutine canon performs the canonical variate analysis
!**      on the sigma arrays
      use cvacom, only: nrowpp,nrowff,sigpp,sigff,sigpf;
      implicit none;
      integer(4),intent(in) :: nynlag;             !=ny*nlag
      integer(4),intent(in) :: nyunlg;             !=(ny+nu)*nlag
      real(8),intent(in) :: eps;                   !tolerance for testing singular values for rank (e.g., 1d-12)
      real(8),intent(out) :: jp(nrowpp,nyunlg);    !mapping from canonical states to past
      real(8),intent(out) :: jf(nrowff,nynlag);    !mapping from canonical states to future
      real(8),intent(out) :: s(nyunlg);            !singular values of sigpp^(-1/2)*sigpf*sigff**(-1/2)
      integer(4),intent(out) :: nrana;             !rank of jf

      integer(4) i,nranf,nranp,lwork,info;
      real(8) s1,si,suma;    !e(nyunlg),work(nyunlg) for Linpack
      real(8),allocatable :: sigpp2i(:,:),sigff2i(:,:),sigff2(:,:);
      real(8),allocatable :: sigffs(:,:), sigpfs(:,:);   !shifted versions of sigff and sigpf
      real(8),allocatable :: work1(:),d(:,:),u(:,:),vt(:,:);
      integer(4) iwork(10*nyunlg);        !LAPACK SVD (only need 8*nyunlg)
      logical(4) :: check =.true., linv =.true.;
      !__________________________________________________

      lwork =nyunlg*(10*nyunlg+4);   !nyunlg*(7*nyunlg+4)
          ! 3*min(M,N)*min(M,N) + max(max(M,N),4*min(M,N)*min(M,N)+4*min(M,N))
      allocate (work1(lwork),sigpp2i(nyunlg,nyunlg),                    &
     &          sigff2i(nynlag,nynlag),sigff2(nynlag,nynlag),           &
     &          d(nyunlg,nyunlg),u(nyunlg,nyunlg),                      &
     &          sigffs(nynlag,nynlag), sigpfs(nyunlg,nynlag),           & !
     &          vt(nyunlg,nyunlg),stat=info);
      if (info /= 0) stop 'allocate work1,... failed';


      !** shift future if sigff and sigpf were computed with more lags than now used
      i =nrowff -nynlag;
      sigffs(:nynlag,:nynlag) =sigff(i+1:nrowff,i+1:nrowff);
      sigpfs(:nyunlg,:nynlag) =sigpf(:nyunlg,i+1:nrowff);

!**     compute svd on sigpp using LAPACK : output is s, u and v
       jp(:nyunlg,:) =sigpp(:nyunlg,:);         !jp will be destroyed in dgesdd, output is V^T
       call DGESDD("A", nyunlg, nyunlg, jp, nrowpp, s, U, nyunlg,       & !divide-and-conquer algorithm
     &             vt, nyunlg, WORK1, lwork, IWORK, info );
!       call DGESVD( "A", "A", nyunlg, nyunlg, jp, nrowpp, s,           !iteration
!     &         U, nyunlg, VT, nyunlg, WORK1, lwork, info )
!      call dsvdc (sigpp,nyunlg,nyunlg,nyunlg,s,e,u,nyunlg,v,nyunlg,    !Linpack (doesn't work correctly)
!     &           work,21,info)
      if (info /= 0) then;
        write(6,1000) info;
 1000   format('error in computing svd of sigpp. info =',i4)
        stop;
      endif;
!     call prtmx (u,nyunlg,nyunlg,nyunlg,"sigpp-u");
!     write (6,'("vt row 1",20(/10g12.5))') vt(1,:nyunlg);
!     call prtmx (vt,nyunlg,nyunlg,nyunlg,"sigpp-vt");
      call prtmx (s,1,1,nyunlg,"sigpp-singVal");

!**     test for positive definite sigpp and rank
      nranp =1;
      s1 =s(1);
      do i=1,nyunlg;
        if (s(i) > eps*s1) then;
          nranp =i;
          suma =dot_product(u(:nyunlg,i),vt(i,:nyunlg)); !v(:nyunlg,i))
          if (suma < 0.d0) then;
            write(6,1020) i,suma;
 1020       format('error - matrix sigpp is not pos. semidefinite at',  &
     &          i2,g13.5)
            call prtmx (u(:nyunlg,i),1,1,nyunlg,"u-col1");
            call prtmx (vt(i,:nyunlg),1,1,nyunlg,"vt-row1");
            write (6,'("dot product =",g13.5)')                         &
     &           dot_product(u(:nyunlg,i),vt(i,:nyunlg));
            stop;
          endif;
        endif;
      enddo;

      if (check) then;
        !** CHECK SVD (sigpp =u*s*v^T)     LINPACK dsvdc doesn't work !!!
        write (6,'("sigpp: nranp, nyunlg =",2i4)') nranp,nyunlg;
        forall (i=1:nyunlg) d(:nyunlg,i) =u(:nyunlg,i)*s(i);
        jp(:nyunlg,:) =matmul(d,vt) -sigpp(:nyunlg,:);
!        call prtmx (jp,nrowpp,nyunlg,nyunlg,"sigpp-u*s*vT")
        call checkZ (jp,nrowpp,nyunlg,"U*S*VT-sigPP");
      endif;

!**     scale u by 1/sqrt(s)
      do i=1,nranp;
        si =1.d0/sqrt(s(i));
        u(:nyunlg,i) =u(:nyunlg,i)*si;
      enddo;

!**     sigpp2i =sigpp^(-1/2) = U * S^(-1/2) * V^T
      sigpp2i(:nyunlg,:) =matmul(u(:nyunlg,:nranp),vt(:nranp,:));

!**     compute svd on sigff
      jf(:nynlag,:nynlag) =sigffs(:nynlag,:nynlag);         !jf will be destroyed in dgesdd
      call DGESDD("A", nynlag, nynlag, jf, nrowff, s, U, nyunlg,        &
     &             vt, nyunlg, WORK1, lwork, IWORK, info );      !output is V^T
!      call DGESVD( "A", "A", nynlag, nynlag, jf, nrowff, s,
!     &         U, nyunlg, VT, nyunlg, WORK1, lwork, info )
!      call dsvdc (sigffs,nynlag,nynlag,nynlag,s,e,u,nyunlg,v,nyunlg,
!     1           work,21,info)
      if (info /= 0) then;
        write(6,1030) info;
 1030   format('error in computing svd of sigff. info =',i4)
        stop;
      endif;
      call prtmx (s,1,1,nynlag,"sigff-singVal");

!**      test for positive definite sigff and rank
      nranf =1;
      s1 =s(1);
      do i=1,nynlag;
        if (s(i) > eps*s1) then;
          nranf =i;
          suma =dot_product(u(:nynlag,i),vt(i,:nynlag));
          if (suma < 0.d0) then;
            write(6,1040) i;
 1040       format('error - matrix sigff is not pos. semidefinite',i3)
            stop;
          endif;
        endif;
      enddo;

      if (check) then;
        !** CHECK SVD (sigff =u*s*v^T)
        write (6,'("sigff: nranf, nynlag =",2i4)') nranf,nynlag;
        forall (i=1:nynlag) d(:nynlag,i) =u(:nynlag,i)*s(i);
        jf(:nynlag,:nynlag) =                                           &
     &               matmul(d(:nynlag,:nynlag),vt(:nynlag,:nynlag))     &
     &                       -sigffs(:nynlag,:nynlag);
!        call prtmx (jf,nrowff,nynlag,nynlag,"sigff-u*s*vT")
        call checkZ (jf,nrowff,nynlag,"U*S*VT-sigFF");
      endif;

!**      scale u by 1./sqrt(s) and form d =U*s^(1/2)
      do i=1,nranf;
        si =1.d0/sqrt(s(i));
        d(:nynlag,i) =u(:nynlag,i)*sqrt(s(i));
        u(:nynlag,i) =u(:nynlag,i)*si;
      enddo;

!**     sigff2i =sigff^(-1/2) =U * S^(-1/2) * V^T
!**     sigff2 =sigff^(1/2) = U * S^(1/2) * V^T
      sigff2i(:nynlag,:nranf) =                                         &
     &            matmul(u(:nynlag,:nranf),vt(:nranf,:nranf));
      sigff2(:nynlag,:nranf) =                                          &
     &            matmul(d(:nynlag,:nranf),vt(:nranf,:nranf));

!**     vt =sigpp2i*sigpf = sigpp^(-1/2)*sigpf
      vt(:nyunlg,:nynlag) =matmul(sigpp2i(:,:),sigpfs(:nyunlg,:));

!**     jp =vt*sigff^(-T/2) =sigpp^(-1/2)*sigpf*sigff**(-T/2)
!**     (jp is temporary storage - will be destroyed in dgesdd)
      jp(:nyunlg,:nynlag) =matmul(vt(:nyunlg,:nynlag),sigff2i(:,:));

!**     compute svd of matrix A (stored in jp)         !jp will be destroyed in dgesdd
      call DGESDD("A", nyunlg, nynlag, jp, nrowpp, s, U, nyunlg,        &
     &             vt, nyunlg, WORK1, lwork, IWORK, info );      !output is V^T
!      call DGESVD( "A", "A", nyunlg, nynlag, jp, nyunlg, s,
!     &         U, nyunlg, VT, nyunlg, WORK1, lwork, info )
!      call dsvdc (sigpfs,nyunlg,nyunlg,nynlag,s,e,u,nyunlg,v,nyunlg,
!     1            work,21,info)
      if (info /= 0) then;
        write(6,'("error in computing svd of a. info =",i4)') info;
        stop 'error in computing svd of a';
      endif;
      call prtmx (s,1,1,nynlag,                                         &
     &              "sigpp^(-1/2)*sigpf*sigff**(-T/2) -singVal");

!**      find rank
      nrana =1;
      s1 =s(1);
      do i=1,nynlag;
        if (s(i) > eps*s1) nrana =i;
      enddo;
      write(6,'("rank of sigpp^(-1/2)*sigpf*sigff**(-T/2) is",i3)')nrana;

!***    past multiplier: Jp = u^T*sigpp**(-1/2)   (store in jp)
      do i=1,nrana;
        jp(i,:nyunlg) =matmul(u(:nyunlg,i),sigpp2i(:nyunlg,:nyunlg));
      enddo;

!***    future multiplier: L = v^t*sigff**(-1/2)    (defined in Larimore papers, not used in remaining code)
!***    future multiplier: jf = sigff**(1/2)*v*S    (store in jf)
      if (linv) then;
        !** output L^(-1)*S =sigff**(1/2)*v*S
        do i=1,nrana;
          jf(:nynlag,i) =matmul(sigff2(:,:),vt(i,:nynlag))*s(i);
        enddo;
      else;       !output L
        jf(:nynlag,:nynlag) =                                           &
     &              matmul(vt(:nynlag,:nynlag),sigff2i(:nynlag,:nynlag));
      endif;

      if (check) then;
        !*** check JP * sigPP * JP^T =I
        d(:nrana,:nrana) =matmul(jp(:nrana,:nyunlg),                    &
     &     matmul(sigpp(:nyunlg,:nyunlg),transpose(jp(:nrana,:nyunlg))));
        forall (i=1:nrana) d(i,i) =d(i,i) -1.d0;
!        call prtmx (u,nyunlg,nyunlg,nyunlg,"JP * sigPP * JP^T -I")
        call checkZ (d,nyunlg,nrana,"Jp * sigPP * Jp^T -I");

        d(:,:) =matmul(u,transpose(u));
        forall (i=1:nyunlg) d(i,i) =d(i,i) -1.d0;
        call checkZ (d,nyunlg,nyunlg,"U * U^T -I");

        if (.not.linv) then;
          !*** check jf * sigFF * jf^T =I
          write (6,'("A: nrana,nynlag =",2i4)') nrana,nynlag;
          d(:nynlag,:nynlag) =matmul(sigffs(:nynlag,:nynlag),           &
     &                               transpose(jf(:nynlag,:nynlag)));
          u(:nynlag,:nynlag) =                                          &
     &                  matmul(jf(:nynlag,:nynlag),d(:nynlag,:nynlag));
          forall (i=1:nynlag) u(i,i) =u(i,i) -1.d0;
!          call prtmx (u,nyunlg,nynlag,nynlag,"jf * sigff * Jf^T -I")
          call checkZ (u,nyunlg,nynlag,"Jf * sigFF * Jf^T -I");
        endif;
      endif;

      deallocate (work1,sigpp2i,sigff2i,sigff2,sigffs,sigpfs,d);
      return;
      end subroutine canon;
!************************************************************
      subroutine maic (kmin,  jf,jp,yu,covyu_yup,directUfeed);  !sigpp,
      !** computes optimal state order using Akaike information criteria.
      !** other variables from cvacom:
      !**   sigpp(nrowpp,nyunlg) =sample correlation matrix for past-past
      !**   s(nyunlg) = singular values of sigpp^(-1/2)*sigpf*sigff**(-1/2) (not used)

      use cvacom, only: nyunlg,nynlag,ny,nu,nyu,m,nlag,nlago,nrana,     &
     &                  nrowpp,nrowff,sigpp,s,pi;
      implicit none;
      logical(4),intent(in) :: directUfeed;
      real(8),intent(in) :: jf(nrowff,nrana);      !mapping from canonical states to future
      real(8),intent(in) :: jp(nrowpp,nyunlg);     !mapping from canonical states to past
      real(8),intent(in) :: yu(nyu,m);             !measurements and controls vs time
      real(8),intent(in) :: covyu_yup(nyu,nyu+nyunlg);
      integer(4),intent(out) :: kmin;              !computed model order for minimum AIC

      integer(4) i,j,k,it,lag,ip,kk,ierr,nlag2;
      logical(4) :: alternateAIC =.false.;
      real(8) sigres(ny*(ny+1)/2),sigres1(ny*(ny+1)/2),x(nyunlg);
      real(8) covyu_p(nyu,nyunlg),yubar(nyu,nyunlg),ubar(nu,nyunlg);
      real(8) siguu(nu,nu),siguui(nu,nu),sigyu(ny,nu),siguy(nu,ny);
      real(8) syyu(ny,ny),yres(ny),eps1,logdet,logdet1,aic(nyunlg),aicm;
!      real(8) d1(nyu+nyunlg,nyu+nyunlg),d2(ny,nyu+nyunlg),d3(ny,ny)
      !__________________________________________________

      write (6,'(//"**********************************"/                &
     &         "**** starting CVA AIC order testing")');

      nlag2 =2*nlag;

      covyu_p(:nyu,:nyunlg) =covyu_yup(:nyu,nyu+1:nyu+nyunlg);
      call prtmx (covyu_p,nyu,nyu,nyunlg,"covyu_p");
      call prtmx (covyu_yup,nyu,ny,ny,"covyy");

      !** yubar =E(yu,p^T)*jp^T =covyu_p * jp^T
      do i=1,nrana;
        yubar(:nyu,i) =matmul(covyu_p(:nyu,:nyunlg),jp(i,:nyunlg));
      enddo;
      call prtmx (yubar,nyu,nyu,nrana,"yubar");

      if (nu > 0 .and. directUfeed) then;
        ubar(:nu,:nrana) =yubar(ny+1:nyu,:nrana);
        siguu(:,:) =covyu_yup(ny+1:nyu,ny+1:nyu);
      endif;

      aicm =1.d20;
      kmin =0;

      !**  loop for order testing
      do k=1,nrana;

        if (nu > 0 .and. directUfeed) then;
          siguui(:,:) =covyu_yup(ny+1:nyu,ny+1:nyu)                     &
     &                 -matmul(ubar(:nu,:k),transpose(ubar(:nu,:k)));
          sigyu(:ny,:nu) =covyu_yup(:ny,ny+1:nyu)                       &
     &                 -matmul(yubar(:ny,:k),transpose(ubar(:nu,:k)));
          siguy(:nu,:ny) =transpose(sigyu(:,:));
          !** form siguy =(siguu)^(-1) * siguy
          call matInvGJ (ierr,  siguui,siguy,  nu,nu,ny);
          if (ierr /= 0) stop 'maic error: inverting siguu';
          syyu(:ny,:ny) =matmul(sigyu(:ny,:nu),siguy(:nu,:ny));
!         call prtmx(syyu,ny,ny,ny,"syyu");
        else;
          syyu(:,:) =0.d0;
        endif;

        kk =0;
        do i=1,ny;
          do j=1,i;
            kk =kk+1;
            sigres(kk) =covyu_yup(i,j)                                  &
     &                 -dot_product(yubar(i,:k),yubar(j,:k)) -syyu(i,j);
          enddo;
!          resvar(i) =sigres(kk)
        enddo;
        write (6,'(/"order ",i3,", sigres =",10g13.5)') k,sigres(:kk);

        logdet1 =0.d0;
        if (alternateAIC) then;
          !** alternate calculation of sigres from Larimore memo 20 feb, 1988
          !** (needs more work)
!          d1(:nyu,:nyu+nyunlg) =covyu_yup(:nyu,:nyu+nyunlg)
!          d1(nyu+1:nyu+nyunlg,:nyu)
!     &        =transpose(covyu_yup(:nyu,nyu+1:nyu+nyunlg))
!          d1(nyu+1:nyu+nyunlg,nyu+1:nyu+nyunlg) =sigpp(:,:)
!          d2(:ny,:k) =matmul(jf(nynlag-ny+1:nynlag,:k),d1)
!          d3(:ny,:ny) =
!     &        matmul(d2(:ny,:k),transpose(jf(nynlag-ny+1:nynlag,:k))
!          kk =0
!          do i=1,ny
!            do j=1,i
!              kk =kk+1
!              sigres(kk) =d3(i,j)
!            enddo
!          enddo

          !** alternate calculation of sigres using actual prediction residual.
          !** This is slow and matches sigres when directUfeed =.false.
          !** It doesn't change the order when directUfeed=T, so avoid!
          sigres1(:) =0.d0;
          do it =nlag+1, m-nlag+1;       !time to predict to
            x(:k) =0.d0;
            do lag=1,nlago;
              ip =(lag-1)*nyu;
              x(:k) =x(:k) +matmul(jp(:k,ip+1:ip+nyu),yu(:nyu,it-lag));
            enddo;
            kk =min(k,nrana);
            yres(:ny) =yu(:ny,it)                                       &
     &                 -matmul(jf(nynlag-ny+1:nynlag,:kk),x(:kk)); !need extra term for direct feed ###
            kk =0;
            do i=1,ny;
              do j=1,i;
                kk =kk+1;
                sigres1(kk) =sigres1(kk) +yres(i)*yres(j);
              enddo;
            enddo;
          enddo;
          sigres1(:kk) =sigres1(:kk)/(m-2*nlag+1);

          !** get cholesky factor of sigres so that log determinant can be computed
          eps1 =1.d-4;
          call cfactor (ierr,  sigres1,ny,eps1);
          if (ierr < 0) stop 'maic error: sigres1 is singular';
          logdet1 =0.d0;
          do j=1,ny;
            kk =(j*(j+1))/2;
            logdet1 =logdet1 +2.d0*log(sigres1(kk));
          enddo;
          write (6,'("order",i3,", alternate sigres =",10g13.5)')       &
     &              k,sigres1(:kk);    !same as sigres
        endif;    !___ end alternateAIC

!        jfp(1:ny,:nyunlg) =matmul(jf(nynlag-ny+1:nynlag,:k),jp(:k,:nyungp))     !### for psd

        !** get cholesky factor of sigres so that log determinant can be computed
        eps1 =1.d-4;
        call cfactor (ierr,  sigres,ny,eps1);
        if (ierr < 0) stop 'maic error: sigres is singular';
        logdet =0.d0;
        do j=1,ny;
          kk =(j*(j+1))/2;
          logdet =logdet +2.d0*log(sigres(kk));
        enddo;
        write (6,'("maic logdet(sigres) =",2g13.5)') logdet,logdet1;

        aic(k) =(m-nlag2+1)*(ny*(1+log(2.d0*pi)) +logdet)               &
     &            +2*(k*(2*ny+nu) +ny*(ny+1)/2);
        if (directUfeed) aic(k) =aic(k) +2*ny*nu;
        if (aic(k) < aicm) then;
          aicm =aic(k);
          kmin =k;
        endif;
      enddo;

      write (6,'("aic =",2(/10g13.5))') aic(:nrana);
      write (6,'("AIC optimal CVA order =",i3)') kmin;

      return;
      end subroutine maic;
!************************************************************
      subroutine maic1 (kmin,  jf,jp,yu,covyu_yup,directUfeed);  !sigpp,
      !** computes optimal state order using Akaike information criteria.
      !** Alternate method (BPG).  This gives same results as subroutine maic, but is slower
      !** other variables from cvacom:
      !**   sigpp(nrowpp,nyunlg) =sample correlation matrix for past-past
      !**   s(nyunlg) = singular values of sigpp^(-1/2)*sigpf*sigff**(-1/2) (not used)

      use cvacom, only: nyunlg,nynlag,ny,nu,nyu,m,nlag,nlago,nrana,     &
     &                  nrowpp,nrowff,sigpp,s,pi;
      implicit none;
      logical(4),intent(in) :: directUfeed;
      real(8),intent(in) :: jf(nrowff,nrana);      !mapping from canonical states to future
      real(8),intent(in) :: jp(nrowpp,nyunlg);     !mapping from canonical states to past
      real(8),intent(in) :: yu(nyu,m);             !measurements and controls vs time
      real(8),intent(in) :: covyu_yup(nyu,nyu+nyunlg);
      integer(4),intent(out) :: kmin;              !computed model order for minimum AIC

      integer(4) i,j,k,it,lag,ip,kk,ierr,nlag2;
      real(8) sigres(ny*(ny+1)/2),sigres1(ny*(ny+1)/2),x(nyunlg);
      real(8) covyu_p(nyu,nyunlg),yubar(nyu,nyunlg),ubar(nu,nyunlg);
      real(8) siguu(nu,nu),siguui(nu,nu),sigyu(ny,nu),siguy(nu,ny);
      real(8) syyu(ny,ny),yres(ny),eps1,logdet,logdet1,aic(nyunlg),aicm;
      real(8) d1(nu+nrowpp,nu+nrowpp),d2(nu+nrowpp,ny),d3(ny,nu+nrowpp);
      !__________________________________________________

      write (6,'(//"**********************************"/                &
     &         "**** starting CVA AIC1 order testing")');

      nlag2 =2*nlag;
      covyu_p(:nyu,:nyunlg) =covyu_yup(:nyu,nyu+1:nyu+nyunlg);
      call prtmx (covyu_p,nyu,nyu,nyunlg,"covyu_p");
      call prtmx (covyu_yup,nyu,ny,ny,"covyy");

      !** yubar =E(yu,p^T)*jp^T =covyu_p * jp^T
      do i=1,nrana;
        yubar(:nyu,i) =matmul(covyu_p(:nyu,:nyunlg),jp(i,:nyunlg));
      enddo;
      call prtmx (yubar,nyu,nyu,nrana,"yubar");

      if (nu > 0 .and. directUfeed) then;
        ubar(:nu,:nrana) =yubar(ny+1:nyu,:nrana);
        siguu(:,:) =covyu_yup(ny+1:nyu,ny+1:nyu);
      endif;

      aicm =1.d20;
      kmin =0;

      !**  loop for order testing
      do k=1,nrana;

        if (nu > 0 .and. directUfeed) then;
          d1(:nu,:nu) =siguu(:,:);
          d1(:nu,nu+1:nu+k) =ubar(:nu,:k);
          d1(nu+1:nu+k,:nu) =transpose(ubar(:nu,:k));
          d1(nu+1:nu+k,nu+1:nu+k) =0.d0;
          forall (i=nu+1:nu+k) d1(i,i) =1.d0;
          d3(:ny,:nu) =covyu_yup(:ny,ny+1:nyu);
          d3(:ny,nu+1:nu+k) =yubar(:ny,:k);
          d2(:nu+k,:ny) =transpose(d3(:ny,:nu+k));
          !** form d2 =d1^(-1) * d2
          call matInvGJ (ierr,  d1,d2,  nu+nrowpp,nu+k,ny);
          if (ierr /= 0) stop 'maic error: inverting siguu';
          syyu(:ny,:ny) =matmul(d3(:ny,:nu+k),d2(:nu+k,:ny));
!         call prtmx(syyu,ny,ny,ny,"syyu");
        else;
          syyu(:,:) =matmul(yubar(:ny,:k),transpose(yubar(:ny,:k)));
        endif;

        kk =0;
        do i=1,ny;
          do j=1,i;
            kk =kk+1;
            sigres(kk) =covyu_yup(i,j) -syyu(i,j);
          enddo;
        enddo;
        write (6,'(/"order ",i3,", sigres =",10g13.5)') k,sigres(:kk);

        !** get cholesky factor of sigres so that log determinant can be computed
        eps1 =1.d-4;
        call cfactor (ierr,  sigres,ny,eps1);
        if (ierr < 0) stop 'maic error: sigres is singular';
        logdet =0.d0;
        do j=1,ny;
          kk =(j*(j+1))/2;
          logdet =logdet +2.d0*log(sigres(kk));
        enddo;
        write (6,'("maic logdet(sigres) =",2g13.5)') logdet;

        aic(k) =(m-nlag2+1)*(ny*(1+log(2.d0*pi)) +logdet)               &
     &            +2*(k*(2*ny+nu) +ny*(ny+1)/2);
        if (directUfeed) aic(k) =aic(k) +2*ny*nu;
        if (aic(k) < aicm) then;
          aicm =aic(k);
          kmin =k;
        endif;
      enddo;

      write (6,'("aic =",2(/10g13.5))') aic(:nrana);
      write (6,'("AIC optimal CVA order =",i3)') kmin;

      return;
      end subroutine maic1;
!*****************************************************
      subroutine checkZ (a,ndim,n,aname);
      !** checks that the a matrix is close to zero and prints elements that are not
      implicit none;
      integer(4),intent(in) :: ndim,n;
      real(8),intent(in) :: a(ndim,n);
      character(*),intent(in) :: aname;

      integer(4) i,j;
      !___________________________________________

      do i=1,n;
        do j=1,n;
          if (abs(a(i,j)) > 2.d-12) then;
            write (6,'("error: ",a," (",2i4," =",g13.5)')               &
     &           aname,i,j,a(i,j);
          endif;
        enddo;
      enddo;

      return;
      end subroutine checkZ;
! ******************************************************************
      subroutine state_space (phi,q,h,rn,bw,gu,ac,                      &
     &                   jp,covyu_yup,n,ndim,directUfeed); !sigpp,nyunlg,ny,nu,
!**      subroutine state_space computes the state space model from the
!**      jp,sigpp,covyu_yup arrays.

!**      Model is
!**         y(i) =Ac*u(i) +H*x(i) +Bw*w(i) +v(i)
!**         x(i+1) =Phi*x(i) +Gu*u(i) +w(i)
!**      where
!**         u(i) is a known control
!**         w(i) is state noise
!**         v(i) is measurement noise
!**
!**      If directUfeed = .false., Ac is constrained to be 0
      use cvacom, only: nyunlg,ny,nu,nyu,nlag,nlago,nrowpp,sigpp;
      implicit none;
      integer(4),intent(in) :: n;               !model order (# states in x)
      integer(4),intent(in) :: ndim;            !row dimensioning of phi,q
      real(8),intent(in) :: jp(nrowpp,nyunlg);  !CVA-computed past-to-state matrix
      real(8),intent(in) :: covyu_yup(ny+nu,ny+nu+nyunlg); !sample covariance of (y,u) (y,u,p)^T
      logical(4),intent(in) :: directUfeed;      !if false, Ac is constrained to be 0
      real(8),intent(out) :: phi(ndim,n);       !computed state transition matrix
      real(8),intent(out) :: h(ny,n);           !computed measurement matrix
      real(8),intent(out) :: q(ndim,n);         !computed state transition matrix
      real(8),intent(out) :: rn(ny,ny);         !computed meas. noise covariance
      real(8),intent(out) :: bw(ny,n);          !computed feed of state noise to measurements
      real(8),intent(out) :: gu(ndim,nu);       !computed feed of control to state
      real(8),intent(out) :: ac(ny,nu);         !computed feed of control to measurements

      integer(4) i,ierr,lwork,info,nplot,i1,i2;
      logical(4) :: check =.true., dprint;
      real(8),allocatable :: d1(:,:),d2(:,:),d3(:,:),d4(:,:);
      real(8),allocatable :: covpi1_pi(:,:),covyu_up(:,:);
      real(8),allocatable :: covpi1_uipi(:,:),covym1_um(:,:);
      real(8),allocatable :: covym1_ym1(:,:),covum_um(:,:),s(:,:);
      real(8),allocatable :: covpi1_pi1(:,:),h1(:,:),covium_um(:,:);
!      real(8),allocatable :: covu_uup(:,:)
!      real(8),allocatable :: work(:),sv(:),u(:,:),vt(:,:)
      complex,allocatable :: phic(:,:),workc(:),wc(:),vl(:,:),vr(:,:);
      real(4),allocatable :: rwork(:);
      !__________________________________________________

      write (6,'(//"**********************************"/                &
     &         "**** starting state_space array calculations")');

      allocate (covpi1_pi(nyunlg,nyunlg),covpi1_uipi(nyunlg,nu+nyunlg), &
     &          d1(ny+n,nyunlg),d2(ny+n,nyu+nyunlg),d3(ny+n,nu+n),      &
     &          d4(ny+n,ny+n),covium_um(nu+n,nu+n),covym1_um(ny+n,nu+n),&
     &          covyu_up(ny+nu,nu+nyunlg),covpi1_pi1(nyunlg,nyunlg),    &
     &          covym1_ym1(ny+n,ny+n),covum_um(nu+n,nu+n),h1(ny,n),     &
     &          s(ny+n,ny+n),stat=ierr);     !covu_uup(nu,(nu+1)*nlag),
      if (ierr /= 0) stop 'allocate covpi1_pi,... failed';

      !** copy cov [(y,u) * (u,p)^T]
      covyu_up(:nyu,:nu+nyunlg) =covyu_yup(:nyu,ny+1:nyu+nyunlg);

      !** covpi1_pi =cov [pi+1 * pi^T]
      covpi1_pi(:nyu,:nyunlg) =covyu_up(:nyu,nu+1:nu+nyunlg);
      covpi1_pi(nyu+1:nyunlg,:nyunlg) =sigpp(:nyunlg-nyu,:nyunlg);

      !** covpi1_pi1 =cov [pi+1 * pi+1^T]
      covpi1_pi1(:nyu,:nyunlg) =covyu_yup(:nyu,:nyunlg);
      covpi1_pi1(nyu+1:nyunlg,:nyu) =                                   &
     &             transpose(covyu_yup(:nyu,nyu+1:nyunlg));
      covpi1_pi1(nyu+1:nyunlg,nyu+1:nyunlg) =                           &
     &              sigpp(:nyunlg-nyu,:nyunlg-nyu);

      dprint =check .and. (n == 4);    !###
      if (dprint) then;
        call prtmx(sigpp,nyunlg,nyunlg,nyunlg,"sigpp");
        call prtmx(covpi1_pi1,nyunlg,nyunlg,nyunlg,"covpi1_pi1");
      endif;

      !** covpi1_uipi =cov [pi+1, (Ui^T , pi^T)] = (cov [pi+1, Ui^T] ,covpi1_pi)
      covpi1_uipi(:nyu,:nu) =covyu_up(:nyu,:nu);
      covpi1_uipi(nyu+1:nyunlg,:nu) =                                   &
     &                  transpose(covyu_up(ny+1:nyu,nu+1:nu+nyunlg-nyu));
      covpi1_uipi(:nyunlg,nu+1:nu+nyunlg) =covpi1_pi(:nyunlg,:nyunlg);

      d1(:ny,:) =0.d0;
      forall (i=1:ny) d1(i,i) =1.d0;
      d1(ny+1:,:) =jp(:n,:);

      d2(:,:nu+nyunlg) =matmul(d1,covpi1_uipi);
      covym1_um(:ny+n,:nu) =d2(:,:nu);
      covym1_um(:ny+n,nu+1:nu+n) =                                      &
     &                matmul(d2(:,nu+1:nu+nyunlg),transpose(jp(:n,:)));

      d2(:,:nyunlg) =matmul(d1,covpi1_pi1);
      covym1_ym1(:ny+n,:ny+n) =matmul(d2(:,:nyunlg),transpose(d1));

      !** cov [(u,m) * (U,m)^T)]  where m =Jp * p and cov (m * m^T] =I
      covum_um(:,:) =0.d0;
      forall (i=nu+1:nu+n) covum_um(i,i) =1.d0;
      covum_um(:nu,:nu) =covyu_up(ny+1:nyu,:nu);
      d2(:nu,:n) =                                                      &
     &     matmul(covyu_up(ny+1:nyu,nu+1:nu+nyunlg),transpose(jp(:n,:)));
      covum_um(:nu,nu+1:) =d2(:nu,:n);
      covum_um(nu+1:,:nu) =transpose(d2(:nu,:n));
      if (n == 4) then;
        call prtmx(covum_um,nu+n,nu+n,nu+n,"covum_um");
      endif;

      !** invert covum_um
      covium_um(:,:) =covum_um(:,:);    !save before inverting
      call matInvGJ (ierr,  covium_um,d2,  nu+n,nu+n,0);
      if (ierr /= 0) stop 'state_space error: inverting covum_um';

      d3(:ny+n,:nu+n) =matmul(covym1_um,covium_um);
      if (nu > 0) then;
        gu(:n,:nu) =d3(ny+1:ny+n,:nu);
        ac(:ny,:nu) =d3(:ny,:nu);
      endif;
      phi(:n,:n) =d3(ny+1:ny+n,nu+1:nu+n);
      h(:ny,:n) =d3(:ny,nu+1:nu+n); !also equal cov (y,p^T)*Jp^T

      call prtmx (phi,ndim,n,n,"phi");
      call prtmx (gu,ndim,n,nu,"gu");
      call prtmx (h,ny,ny,n,"h");
      call prtmx (ac,ny,ny,nu,"ac");

      !** now form S to get Q and Rn
      if (directUfeed) then;
        d4(:ny+n,:ny+n) =matmul(d3(:ny+n,:nu+n),transpose(covym1_um));
        s(:,:) =covym1_ym1(:,:) -d4(:,:);
      else;    !no directUfeed from control to measurements
        ac(:ny,:nu) =0.d0;
        h1(:ny,:n) =covym1_um(:ny,nu+1:nu+n);
!        h1(:ny,:n) =
!     &        matmul(covyu_up(:ny,nu+1:nu+nyunlg),transpose(jp(:n,:)))     !same as above, but not necessarily same as h
        call prtmx (h1,ny,ny,n,"h1");
        d3(:ny,:nu) =0.d0;                !explicitly zero A (direct feed) term
        d3(:ny,nu+1:nu+n) =h(:ny,:n);
        d4(:ny+n,:ny+n) =matmul(d3(:ny+n,:nu+n),transpose(covym1_um));
        s(:,:) =covym1_ym1(:,:) -(d4(:,:) +transpose(d4));
        d1(:ny+n,:nu+n) =matmul(d3(:ny+n,:nu+n),covium_um);
        s(:,:) =s(:,:)                                                  &
     &         +matmul(d1(:ny+n,:nu+n),transpose(d3(:ny+n,:nu+n)));
      endif;

      q(:n,:n) =s(ny+1:,ny+1:);
      !** invert Q  (should use pseudo-inverse for Q as it may not be positive definite)
      d1(:n,:n) =q(:n,:n);
      call matInvGJ (ierr,  d1,bw,  ny+n,n,0);
      if (ierr /= 0) stop 'state_space error: inverting q';
      bw(:ny,:n) =matmul(s(:ny,ny+1:ny+n),d1(:n,:n));
      rn(:ny,:ny) =s(:ny,:ny)                                           &
     &         -matmul(bw(:ny,:n),transpose(s(:ny,ny+1:ny+n)));

      call prtmx (q,ndim,n,n,"q");
      call prtmx (bw,ny,ny,n,"bw");
      call prtmx (rn,ny,ny,ny,"rn");

      if (dprint) then;
        !** check calculations
        d3(:ny,:nu) =ac(:ny,:nu);
        d3(:ny,nu+1:nu+n) =h(:ny,:n);
        d3(ny+1:ny+n,:nu) =gu(:n,:nu);
        d3(ny+1:ny+n,nu+1:nu+n) =phi(:n,:n);

        d1(:ny+n,:nu+n) =matmul(d3(:ny+n,:nu+n),covum_um(:nu+n,:nu+n));
        d2(:ny+n,:ny+n) =                                               &
     &             matmul(d1(:ny+n,:nu+n),transpose(d3(:ny+n,:nu+n)));
        d1(:ny,:ny) =matmul(matmul(bw(:ny,:n),q(:n,:n)),                &
     &                      transpose(bw(:ny,:n)));
        call prtmx(d1,ny+n,ny,ny,"bw*q*bw^T");
        d1(:ny,:ny) =d1(:ny,:ny) +rn(:ny,:ny);
        call prtmx(d1,ny+n,ny,ny,"rn +bw*q*bw^T");
        d2(:ny,:ny) =d2(:ny,:ny) +d1(:ny,:ny)                           &
     &                 -covyu_yup(:ny,:ny);
        call prtmx(d2,ny+n,ny,ny,"covyy-err");
        d2(:ny,ny+1:ny+n) =d2(:ny,ny+1:ny+n)+matmul(bw(:ny,:n),q(:n,:n))&
     &            -covym1_ym1(:ny,ny+1:ny+n);
        call prtmx(d2(:ny,ny+1:ny+n),ny,ny,n,"covym-err");
        d1(:n,:nyunlg) =                                                &
     &       matmul(jp(:n,:nyunlg),covpi1_pi1(:nyunlg,:nyunlg));
        d1(:n,:n) =matmul(d1(:n,:nyunlg),transpose(jp(:n,:nyunlg)));
        d2(ny+1:ny+n,ny+1:ny+n) =d2(ny+1:ny+n,ny+1:ny+n) +q(:n,:n)      &
     &         -d1(:n,:n);
!        forall (i=ny+1:ny+n) d2(i,i) =d2(i,i) -1.d0   !d1 should be I
        call prtmx(d2(ny+1:ny+n,ny+1:ny+n),n,n,n,"covmm-err");
        call prtmx(d1,ny+n,n,n,"covJpi1_pi1Jt");   !~= I to 3 digits

        !** compute complex right eigenvalues
        lwork =10*n;       !should only need 4*n but more helps
        allocate (phic(n,n),workc(lwork),wc(n),                         &
     &            vl(n,n),vr(n,n),rwork(10*n),stat=info);
        if (info /= 0) stop 'allocate phic,... failed';
        phic(:,:) =phi(:,:);
        call CGEEV("N", "V", N, phic, n, wc, VL, n, VR, n,              &
     &                  WORKc, LWORK, RWORK, INFO );
        write (6,'("returned from cgeev: ",i3)') info;
        write (6,'("phi complex eigenvalues",10(/2g13.5))') wc(:n);
        deallocate (phic,workc,wc,vl,vr,rwork);
      endif;

      deallocate (covpi1_pi,covpi1_uipi,d1,d2,d3,d4,covym1_um,          &
     &            covyu_up,covym1_ym1,covum_um,h1,s); !,covu_uup)
      return;
      end subroutine state_space;
! ******************************************************************
      subroutine state_space1 (phi,q,h,rn,bw,gu,ac,                     &
     &                         jp,n,ndim,directUfeed);  !nyunlg,ny,nu,sigp1p1,nyuyunlg,
!**      subroutine state_space computes the state space model from the
!**      jp,sigpp,covyu_yup arrays.  This version uses the innovations
!**      form of the state space model, which allows recursive computation
!**      of the measurement noise covariance.  This is more accurate when
!**      measurement noise is small.

!**      As of 10/20/09, this version without the recursive computation of Rn
!**      checks with subroutine state_space (at least for directUfeed)

!**      Model is
!**         y(i) =Ac*u(i) +H*x(i) +v(i)               !v includes Bw*q =0
!**         x(i+1) =Phi*x(i) +Gu*u(i) +K*v(i) +w(i)   !this is not a complete model as K*v(i) is rank ny
!**      where
!**         u(i) is a known control
!**         w(i) is state noise
!**         v(i) is measurement noise
!**
!**      If directUfeed = .false., Ac is constrained to be 0
!**      other variables from cvacom: sigp1p1(nyuyunlg,nyuyunlg) =sample covariance of past+1,past+1, i.e., (y,u,p) (y,u,p)

      use cvacom, only: nyunlg,ny,nu,nyu,nlag,nyuyunlg,nrowpp,sigp1p1;
      implicit none;
      integer(4),intent(in) :: n;               !model order (# states in x)
      integer(4),intent(in) :: ndim;            !row dimensioning of phi,q
      real(8),intent(in) :: jp(nrowpp,nyunlg);  !CVA-computed past-to-state matrix
      logical(4),intent(in) :: directUfeed;      !if false, Ac is constrained to be 0
      real(8),intent(out) :: phi(ndim,n);       !computed state transition matrix
      real(8),intent(out) :: h(ny,n);           !computed measurement matrix
      real(8),intent(out) :: q(ndim,n);         !computed state transition matrix
      real(8),intent(out) :: rn(ny,ny);         !computed meas. noise covariance
      real(8),intent(out) :: bw(ny,n);          !computed feed of state noise to measurements
      real(8),intent(out) :: gu(ndim,nu);       !computed feed of control to state
      real(8),intent(out) :: ac(ny,nu);         !computed feed of control to measurements

      integer(4) i,i1,i2,ierr,lwork,info,nplot;
      integer(4) nup,num,nyup,nyum,nuyunlg;
      logical(4) :: check =.true., dprint;
      real(8),allocatable :: d1(:,:),d2(:,:),d3(:,:),d4(:,:);
      real(8),allocatable :: covup_up(:,:),covy_up(:,:);
      real(8),allocatable :: covum_um(:,:),covium_um(:,:);
      real(8),allocatable :: covvum_vum(:,:),covivum_vum(:,:);
      real(8),allocatable :: covy_um(:,:),covp1_vum(:,:),k(:,:);
!      real(8),allocatable :: covu_uup(:,:)
!      real(8),allocatable :: work(:),sv(:),u(:,:),vt(:,:)
      complex,allocatable :: phic(:,:),workc(:),wc(:),vl(:,:),vr(:,:);
      real(4),allocatable :: rwork(:);
      !__________________________________________________

      write (6,'(//"**********************************"/                &
     &         "**** starting state_space1 array calculations")');

      nup =nu +nyunlg;
      nyup =nyu+nyunlg;
      nuyunlg =nu+nyunlg;
      num =nu+n;
      nyum =nyu+n;
!      write (6,'(10i4)') nyu,nlag,nup,nyup,nuyunlg,num,nyum
      dprint =check .and. (n == 4);    !###

      allocate (covup_up(nup,nup),covy_up(ny,nup),                      &
     &          d1(nuyunlg,nuyunlg),d2(nyup,nyuyunlg),d3(nyunlg,nyum),  &
     &          d4(ny+n,ny+n),covum_um(num,num),covium_um(num,num),     &
     &          covy_um(ny,num),covvum_vum(nyum,nyum),                  &
     &          covivum_vum(nyum,nyum),covp1_vum(nyunlg,nyum),          &
     &          k(n,ny),stat=ierr);    !covu_uup(nu,(nu+1)*nlag),
      if (ierr /= 0) then;
        write (6,'("allocate ier =",i6)') ierr;
        stop 'allocate covup_up,... failed';
      endif;

      covup_up(:nup,:nup) =sigp1p1(ny+1:ny+nup,ny+1:ny+nup);
      covy_up(:ny,:nup) =sigp1p1(:ny,ny+1:ny+nup);

      !** form (I O / 0 J) and multiply to get cov [(u m) (u m)^T]
      d1(:,:) =0.d0;
      forall (i=1:nu) d1(i,i) =1.d0;
      d1(nu+1:num,nu+1:nup) =jp(:n,:nyunlg);
      covum_um(:num,:num) =                                             &
     &         matmul(matmul(d1(:num,:nup),covup_up(:nup,:nup)),        &
     &                                   transpose(d1(:num,:nup)));
      covy_um(:ny,:num) =                                               &
     &         matmul(covy_up(:ny,:nup),transpose(d1(:num,:nup)));

!      if (dprint) then
!        call prtmx(covum_um,nu+n,nu+n,nu+n,"covum_um")
!      endif

      !** invert covum_um(:num,:num)
      covium_um(:num,:num) =covum_um(:num,:num);
      d2(:,:) =0.d0;
      call matInvGJ (ierr,  covium_um,d2,  num,num,0);
      if (ierr /= 0) stop 'state_space error: inverting covum_um';

      if (nu > 0) then;
        if (directUfeed) then;
          ac(:ny,:nu) =matmul(covy_um(:ny,:num),covium_um(:num,:nu));
        else;
          ac(:,:) =0.d0;
        endif;
      endif;
      h(:ny,:n) =matmul(covy_um(:ny,:num),covium_um(:num,nu+1:num));

      d2(:,:) =0.d0;
      forall (i=1:nyu) d2(i,i) =1.d0;
      if (nu > 0) d2(:ny,ny+1:nyu) =-ac(:ny,:nu);
      d2(:ny,nyu+1:nyuyunlg) =-matmul(h(:ny,:n),jp(:n,:nyunlg));
      d2(nyu+1:nyum,nyu+1:nyuyunlg) =jp(:n,:nyunlg);

      !** form cov [(v,u,m) (v,u,m)^T) and cov [p_i+1 (v,u,m)^T)
      covvum_vum(:nyum,:nyum) =                                         &
     &  matmul(matmul(d2(:nyum,:nyuyunlg),sigp1p1(:nyuyunlg,:nyuyunlg)),&
     &                 transpose(d2(:nyum,:nyuyunlg)));
      covp1_vum(:nyunlg,:nyum) =matmul(sigp1p1(:nyunlg,:nyuyunlg),      &
     &                              transpose(d2(:nyum,:nyuyunlg)));

      !** invert covvum_vum (save in covivum_vum)
      covivum_vum(:nyum,:nyum) =covvum_vum(:nyum,:nyum);
      call matInvGJ (ierr,  covivum_vum,d2,  nyum,nyum,0);
      if (ierr /= 0) stop 'state_space error: inverting covvum_vum';

      d3(:nyunlg,:nyum) =matmul(covp1_vum(:nyunlg,:nyum),               &
     &                     covivum_vum(:nyum,:nyum));
      k(:n,:ny) =matmul(jp(:n,:nyunlg),d3(:nyunlg,:ny));
      if (nu > 0) gu(:n,:nu)=matmul(jp(:n,:nyunlg),d3(:nyunlg,ny+1:nyu));
      phi(:n,:n) =matmul(jp(:n,:nyunlg),d3(:nyunlg,nyu+1:nyum));

      call prtmx (phi,ndim,n,n,"phi");
      call prtmx (gu,ndim,n,nu,"gu");
      call prtmx (h,ny,ny,n,"h");
      call prtmx (ac,ny,ny,nu,"ac");

      !** approximate calculation (should replace with recursive)
      !** This includes both measurement noise and the effects of q,
      !** and matches Rn +Bw*Q*Bw^T of subroutine state_space
      rn(:ny,:ny) =sigp1p1(:ny,:ny) -                                   &
     &    matmul(matmul(covp1_vum(:ny,ny+1:nyum),covium_um(:num,:num)), &
     &                 transpose(covp1_vum(:ny,ny+1:nyum)));
      !** this gives same answer as eqn above
!      rn(:ny,:ny) =covvum_vum(:ny,:ny) -
!     &    matmul(matmul(covvum_vum(:ny,ny+1:nyum),covium_um(:num,:num)),
!     &                 transpose(covvum_vum(:ny,ny+1:nyum)))

!      q(:n,:n) =matmul(matmul(k(:n,:ny),rn(:ny,:ny)),       !this is wrong ###
!     &                transpose(k(:n,:ny)))
      q(:n,:n) =matmul(matmul(jp(:n,:nyunlg),sigp1p1(:nyunlg,:nyunlg)), &
     &                transpose(jp(:n,:nyunlg)));
      d1(:n,:nu) =gu(:n,:nu);
      d1(:n,nu+1:nu+n) =phi(:n,:n);
      q(:n,:n) =q(:n,:n)-matmul(matmul(d1(:n,:num),covum_um(:num,:num)),&
     &                          transpose(d1(:n,:num)));

      !** still need bw ###
      d1(:ny,:nu) =ac(:ny,:nu);
      d1(:ny,nu+1:num) =h(:ny,:n);
      d2(:ny,:num) =matmul(d1(:ny,:nu+n),covum_um(:num,:num));
      d1(:nu,:n) =transpose(gu(:n,:nu));
      d1(nu+1:num,:n) =transpose(phi(:n,:n));
      bw(:ny,:n)=matmul(sigp1p1(:ny,:nyunlg),transpose(jp(:n,:nyunlg))) &
     &            -matmul(d2(:ny,:num),d1(:num,:n));
      !** invert q
      d1(:n,:n) =q(:n,:n);
      call matInvGJ (ierr,  d1,d3,  nuyunlg,n,0);
      if (ierr /= 0) stop 'state_space error: inverting q';
      bw(:ny,:n) =matmul(bw(:ny,:n),d1(:n,:n));

      !** compute rn without Bw*q*Bw^T
      rn(:ny,:ny) =rn(:ny,:ny) -                                        &
     &      matmul(matmul(bw(:ny,:n),q(:n,:n)),transpose(bw(:ny,:n)));

      call prtmx (q,ndim,n,n,"q");
      call prtmx (bw,ny,ny,n,"bw");
      call prtmx (rn,ny,ny,ny,"rn");

      if (dprint) then;
        !** check calculations
        d3(:ny,:nu) =ac(:ny,:nu);
        d3(:ny,nu+1:nu+n) =h(:ny,:n);
        d3(ny+1:ny+n,:nu) =gu(:n,:nu);
        d3(ny+1:ny+n,nu+1:nu+n) =phi(:n,:n);

        d1(:ny+n,:nu+n) =matmul(d3(:ny+n,:nu+n),covum_um(:nu+n,:nu+n));
        d2(:ny+n,:ny+n) =                                               &
     &             matmul(d1(:ny+n,:nu+n),transpose(d3(:ny+n,:nu+n)));
!        d1(:ny,:ny) =matmul(matmul(bw(:ny,:n),q(:n,:n)),
!     &                      transpose(bw(:ny,:n)))
!        call prtmx(d1,ny+n,ny,ny,"bw*q*bw^T")

        !** compute complex right eigenvalues
        lwork =10*n;       !should only need 4*n but more helps
        allocate (phic(n,n),workc(lwork),wc(n),                         &
     &            vl(n,n),vr(n,n),rwork(10*n),stat=info);
        if (info /= 0) stop 'allocate phic,... failed';
        phic(:,:) =phi(:,:);
        call CGEEV("N", "V", N, phic, n, wc, VL, n, VR, n,              &
     &                  WORKc, LWORK, RWORK, INFO );
        write (6,'("returned from cgeev: ",i3)') info;
        write (6,'("phi complex eigenvalues",10(/2g13.5))') wc(:n);
        deallocate (phic,workc,wc,vl,vr,rwork);
      endif;

      deallocate (covup_up,covy_up,d1,d2,d3,d4,covum_um,covium_um,      &
     &          covy_um,covvum_vum,covivum_vum,covp1_vum,k);
!     &          covu_uup)

      return;
      end subroutine state_space1;
!******************************************************************
      subroutine KF_logLike (logDet,  phi,q,h,rn,bw,gu,ac,yu,n,ndim);
!     &                            n,ndim,ny,nu,m,nlag)
!**  compute log-determinant of sample sig_yres-yres from a Kalman filter
!**  using the CVA state space model

!**      Model is
!**         y(i) =Ac*u(i) +H*x(i) +Bw*w(i) +v(i)
!**         x(i+1) =Phi*x(i) +Gu*u(i) +w(i)
!**      where
!**         u(i) is a known control
!**         w(i) is state noise
!**         v(i) is measurement noise
!**
      use cvacom, only: ny,nu,nyu,m,nlag,nlago;
      implicit none;
      integer(4),intent(in) :: n;               !model order (# states in x)
      integer(4),intent(in) :: ndim;            !row dimensioning of phi,q
!      logical(4),intent(in) :: directUfeed       !if false, Ac is constrained to be 0
      real(8),intent(in) :: phi(ndim,n);        !computed state transition matrix
      real(8),intent(in) :: h(ny,n);            !computed measurement matrix
      real(8),intent(in) :: q(ndim,n);          !computed state transition matrix
      real(8),intent(in) :: rn(ny,ny);          !computed meas. noise covariance
      real(8),intent(in) :: bw(ny,n);           !computed feed of state noise to measurements
      real(8),intent(in) :: gu(ndim,nu);        !computed feed of control to state
      real(8),intent(in) :: ac(ny,nu);          !computed feed of control to measurements
      real(8),intent(in) :: yu(ny+nu,m);        !measurements
      real(8),intent(out) :: logDet;            !log determinant of sample residual covariance

      integer(4) i,ierr,it,j,kk;
      logical(4) :: dprint =.false.;
      real(8) :: x(n),p(n,n);
      real(8) sigres(ny*(ny+1)/2),yres(ny),eps1,uc(nu);
      real(8) hp(ny,n),c0(ny,ny),c(ny,ny),ci(ny,ny),kg(n,ny);
      real(8) bq(ny,n),yp(ny),dqr(n,ny),dqrh(n,n),dqrdqt(n,n);
      real(8) d1(n,n),bwqbwt(ny,ny);!,hb(ny,ny)
      !__________________________________________________

!      nyu =ny+nu

      !** initialize arrays
      sigres(:) =0.d0;
      x(:) =0.d0;
      p(:n,:n) =0.d0;       !state error covariance
      forall (i=1:n) p(i,i) =1.d0;   !from cov [m * m*T] =I

      bq(:ny,:n) =matmul(bw(:ny,:n),q(:n,:n));
      bwqbwt(:ny,:ny) =matmul(bq(:ny,:n),transpose(bw(:ny,:n)));
      ci(:ny,:ny) =rn(:ny,:ny) +bwqbwt(:ny,:ny);
      call matInvGJ (ierr,  ci,kg,  ny,ny,0);
      if (ierr /= 0) stop 'KF_logLike error: inverting rn+bwqbwt';

      !** dqr =q*bw^T*rn^(-1), dqrh =dqr*h, dqrdqt =dqr*rn*dqr^T
      dqr(:n,:ny) =matmul(transpose(bq(:ny,:n)),ci(:ny,:ny));
      dqrh(:n,:n) =matmul(dqr(:n,:ny),h(:ny,:n));
      dqrdqt(:n,:n) =matmul(dqr(:n,:ny),bq(:ny,:n));

      do it =1, m-nlag+1;       !time to predict to
        uc(:nu) =yu(ny+1:nyu,it);
        yp(:ny) =matmul(h(:ny,:n),x(:n));
        if (nu > 0) yp(:ny) =yp(:ny) +matmul(ac(:ny,:nu),uc(:nu));
        yres(:ny) =yu(:ny,it) -yp(:ny);
!        write (6,'(i3,", y =",10g13.5)') it,yu(:ny,it)
!        write (6,'(2x,"yres =",10g13.5)') yres(:ny)
        !** form hp =h*p and residual covariance c
        hp(:ny,:n) =matmul(h(:ny,:n),p(:n,:n));
        c0(:ny,:ny) =                                                   &
     &          matmul(hp(:ny,:n),transpose(h(:ny,:n))) +rn(:ny,:ny);
        c(:ny,:ny) =c0(:ny,:ny) +bwqbwt(:ny,:ny);   !add direct state noise contribution
        ci(:ny,:ny) =c(:ny,:ny);
        call matInvGJ (ierr,  ci,kg,  ny,ny,0);
        if (ierr /= 0) stop 'KF_logLike error: inverting c';

        !** measurement update
        kg(:n,:ny) =matmul(transpose(hp(:ny,:n)),ci(:ny,:ny));
        x(:n) =x(:n) +matmul(kg(:n,:ny),yres(:ny));
        p(:n,:n) =p(:n,:n) -matmul(kg(:n,:ny),hp(:ny,:n));

        if (dprint .and. n == 3) then;
          call prtmx (hp,ny,ny,n,"hp");
          call prtmx (c0,ny,ny,ny,"c0");
          call prtmx (c,ny,ny,ny,"c0+bw*q*bwT");
!          ci(:,:) =matmul(c,ci)
!          call prtmx (ci,ny,ny,ny,"ci*c")
          call prtmx (kg,n,n,ny,"kg");
          call prtmx (p,n,n,n,"P+");
          write (12,'(10g13.5)') real(it),yu(:ny,it),                   &
     &              yp(:ny),yres(:ny),(sqrt(c(i,i)),i=1,ny);
        endif;

        if (it > nlag) then;
          kk =0;
          do i=1,ny;
            do j=1,i;
              kk =kk+1;
              sigres(kk) =sigres(kk) +yres(i)*yres(j);
            enddo;
          enddo;
        endif;

        !** time update for next meas
        yres(:ny) =yu(:ny,it) -matmul(h(:ny,:n),x(:n));     !compute a posteriori measurement residual  ### this makes test less sensitive ??
        x(:n) =matmul(phi(:n,:n),x(:n)) +matmul(dqr(:n,:ny),yres(:ny));
        if (nu > 0) x(:n) =x(:n) +matmul(gu(:n,:nu),uc(:nu));
        d1(:n,:n) =phi(:n,:n) -dqrh(:n,:n);
        p(:n,:n) =matmul(matmul(d1(:n,:n),p(:n,:n)),                    &
     &              transpose(d1(:n,:n))) +q(:n,:n) -dqrdqt(:n,:n);
      enddo;   !end time loop

      sigres(:kk) =sigres(:kk)/(m-2*nlag+1);
      call prtmx (sigres,1,1,kk,'sigres');

      !** get cholesky factor of sigres so that log determinant can be computed
      eps1 =1.d-4;
      call cfactor (ierr,  sigres,ny,eps1);
      if (ierr < 0) then;
        write (6,'("KF_logLike error: sigres is singular")');
        logdet =1.d10;
        return;
      endif;

      logdet =0.d0;
      do j=1,ny;
        kk =(j*(j+1))/2;
        logdet =logdet +2.d0*log(sigres(kk));
      enddo;
      write (6,'("logdet =",2g13.5)') logdet;

      return;
      end subroutine KF_logLike;
!****************************************************************
      subroutine plotPsdY (nplot,n,jp,jf,rn,dt,plotFile);
!**      generate a plot file of measurement power spectral density for
!**      white process noise input: given the past and future mapping matrices
!**      and measurement noise covariance
      use cvacom, only:ny,nyu,nyunlg,nynlag,nlago,nrowpp,nrowff;
      implicit none;
      character(*),intent(in) :: plotFile;   !plot file name
      real(8),intent(in) :: dt;              !time step (sec)
      integer(4),intent(in) :: n;            !state order
!      integer(4),intent(in) :: ny            !# of measurements at each time
!      integer(4),intent(in) :: nyu           !# of measurements + controls at each time
      integer(4),intent(in) :: nplot;        !# of points to be plotted
!      integer(4),intent(in) :: nyunlg        !=nyu*nlag
!      integer(4),intent(in) :: nynlag        !=ny*nlag
!      integer(4),intent(in) :: nlag          !# of lags used in jp, jf
      real(8),intent(in) :: jp(nrowpp,nyunlg); !past to model state mapping
      real(8),intent(in) :: jf(nrowff,n);    !model state to future mapping
      real(8),intent(in) :: rn(ny,ny);       !measurement noise covariance (not used?)

      integer(4) i,j,lag,ip;
      real(8) ypsd(ny),ypmag,delf;
      real(8) :: pi=3.141592653589793d0;
      real(4) theta;
      complex x(n,ny),y(ny,ny),z,z1;
      !____________________________________

      open (10,file=plotfile,status='unknown');
      write (10,1000);
 1000 format('variables ="frequency", "y1-psd", "y2-psd"')

      delF =1.d0/(2*nplot*dt);   !factor 2 for negative frequencies
      do i=0,nplot;
        theta =i*2.*pi*delf;
        z =cmplx(cos(theta),sin(theta));
        z1 =z;
        x(:n,:) =0.0;
        do lag=1,nlago;
          ip =(lag-1)*nyu;
          z1 =z1*z;
          x(:n,:) =x(:n,:) +jp(:n,ip+1:ip+ny)/z1;
        enddo;
        forall (j=1:ny)                                                 &
     &     y(j,:) =matmul(jf(nynlag-ny+j,:n),x(:n,:));
!        forall (j=1:n) vr(:ny,j) =vr(:ny,j)*sqrt(q(j,j))     !assume q independent
        do j=1,ny;
          ypmag =dot_product(y(j,:),conjg(y(j,:)));
          if (i > 0) ypmag =2.*ypmag;     !factor 2 for negative frequencies
          ypsd(j) =10.*log10(ypmag);
        enddo;
        write (10,'(4g13.5)') i*delf,ypsd(:ny);
      enddo;

      close(10);
      return;
      end subroutine plotPSDY;
!****************************************************************
      subroutine plotPsdY2 (nplot,dt,plotFile);
!**      generate a plot file of measurement power spectral density for
!**      white process noise input given the past-past covariance matrix.
!**      This doesn't give accurate PSD when # lags is low, even though CVA model is good.

      use cvacom, only:ny,nyu,nyunlg,nlago,sigp1p1;
      implicit none;
      character(*),intent(in) :: plotFile;   !plot file name
      real(8),intent(in) :: dt;              !time step (sec)
      integer(4),intent(in) :: nplot;        !# of points to be plotted
!      integer(4),intent(in) :: ny            !# of measurements at each time
!      integer(4),intent(in) :: nyu           !# of measurements + controls at each time
!      integer(4),intent(in) :: nyunlg        !=nyu*nlago
!      integer(4),intent(in) :: nlag          !# of lags used in jp, jf
!      real(8),intent(in) :: sigpp(nyunlg,nyunlg)   !sample covariance of past-past

      integer(4) i,j,lag,ip;
      real(8) ypsd(ny),delf;
      real(8) :: pi=3.141592653589793d0;
      real(4) theta;
      complex yp(ny,nyunlg),yy(ny,ny),z,z1,za(nlago);
      !____________________________________

      open (10,file=plotfile,status='unknown');
      write (10,1000);
 1000 format('variables ="frequency (Hz)", "y1-psd", "y2-psd"')

      delF =1.d0/(2*nplot*dt);
      do i=0,nplot;
        theta =i*2.*pi*delf;
        z =cmplx(cos(theta),-sin(theta));    !unit delay
!        za(1) =1.
!        do j=2,nlag
!          za(j) =za(j-1)*z
!        enddo

        z1 =1.;
        yp(:,:) =sigp1p1(:ny,:ny);     !0.0
        do lag=1,nlago;
          ip =(lag-1)*nyu;
          z1 =z1*z;
          yp(:,:) =yp(:,:) +sigp1p1(ip+1:ip+ny,:nyunlg)*z1;
        enddo;

        z1 =1.;
        yy(:,:) =0.0;
        do lag=1,nlago;
          ip =(lag-1)*nyu;
          z1 =z1*z;
          yy(:,:) =yy(:,:) +yp(:,ip+1:ip+ny)*conjg(z1);
        enddo;
!        write (6,'(i3,", yy-diag =",4g13.5)') i,(yy(j,j),j=1,ny)

        do j=1,ny;
          ypsd(j) =2.*dt*sqrt(real(yy(j,j))**2+aimag(yy(j,j))**2)/nlago;   !factor of 2 for negative power,
          ypsd(j) =10.*log10(ypsd(j));
        enddo;
        write (10,'(6g13.5)') i*delf, ypsd(:ny);
      enddo;

      close(10);
      return;
      end subroutine plotPSDY2;
!****************************************************************
      subroutine plotPsd2 (nplot,n,ndim,ny,nu,phi,q,h,rn,bw,gu,ac,uc,mu,&
     &                    dt,plotFile);
!**      generate a plot file of measurement power spectral density for
!**      white process noise input: given the state transition,
!**      state noise, measurementsensitivity matrix.
!**      Also plot contribution from control inputs ?
!**      This version does not use Cholesky factorization

!**      State space Model is
!**         y(i) =Ac*u(i) +H*x(i) +Bw*w(i) +v(i)
!**         x(i+1) =Phi*x(i) +Gu*u(i) +w(i)
!**      where
!**         u(i) is a known control
!**         w(i) is state noise
!**         v(i) is measurement noise

      implicit none;
      character(*),intent(in) :: plotFile;   !plot file name
      integer(4),intent(in) :: ny;           !# of measurements at each time
      integer(4),intent(in) :: nu;           !# of controls at each time
      integer(4),intent(in) :: nplot;        !# of points to be plotted (e.g., 100)
      integer(4),intent(in) :: n;            !state order
      integer(4),intent(in) :: ndim;         !row dimension of phi,q (ndim >= n)
      integer(4),intent(in) :: mu;           !# of time samples in uc
      real(8),intent(in) :: dt;              !time step (sec)
      real(8),intent(in) :: phi(ndim,n);     !state transition matrix
      real(8),intent(in) :: q(ndim,n);       !process noise covariance
      real(8),intent(in) :: h(ny,n);         !meas sensitivity matrix
      real(8),intent(in) :: rn(ny,ny);       !meas noise covariance
      real(8),intent(in) :: bw(ny,ndim);     !state noise mapping to measurements
      real(8),intent(in) :: gu(ndim,nu);     !control mapping to states
      real(8),intent(in) :: ac(ny,nu);       !control mapping to measurements
      real(8),intent(in) :: uc(nu,mu);       !time samples of the controls

      integer(4) ip,i,j,info,nrhs,iwork(n);
      real(8) ypsd(ny),delF,pmag(nu),pusum(nu),psum(ny),scalu;
      real(8) :: pi=3.141592653589793d0;
      real(4) theta;
      complex vl(n,n),vr(n,n),zimphi(n,n),z,z1,uz(nu);
      complex hx(ny,n),hxc(n,ny),pv(ny,nu),pw(ny,ny),d1(ny,n);
      !____________________________________

      open (10,file=plotfile,status='unknown');
      write (10,1000);
 1000 format('variables ="frequency (Hz)", "y1-psd", "y2-psd"')
      delF =1.d0/(2*nplot*dt);
      scalu =sqrt(1.d0/mu);
      pusum(:) =0.d0;
      psum(:) =0.d0;

      do ip=0,nplot;
        theta =ip*2.*pi*delF*dt;
        z =cmplx(cos(theta),sin(theta));

        if (nu > 0) then;
          !** compute z-transform of control data u
          z1 =cmplx(1.,0.);
          uz(:) =0.;
          do i=mu,1,-1;
            uz(:) =uz(:) +uc(:,i)*conjg(z1);
            z1 =z1*z;
          enddo;
          uz(:) =uz(:)*scalu;
          pmag(:) =2.*dt*(real(uz(:))**2 +aimag(uz(:))**2);
          pusum(:) =pusum(:) +pmag(:);
!          write (6,'("f =",g13.5,", uz =",2g13.5)')
!          write (6,'(5g13.5)') ip*delF, pmag   !###
        endif;

        zimphi(:n,:n) =-phi(:n,:n);
        forall (j=1:n) zimphi(j,j) =zimphi(j,j) +z;
        !** invert z*I-phi
        nrhs =n;
        vl(:,:) =0.;
        do j=1,n;
          vl(j,j) =cmplx(1.,0.);
        enddo;
        vr(:,:) =zimphi(:,:);
        !** vl = complex matrix inversion of (zI-phi)
        call CGESV( n,nrhs, zimphi, n, iwork, vl, n, info);
        if (info /= 0) stop 'error inverting z*I-phi';

        !** check that zimphi(:,:) * vl =I
!       vr(:,:) =matmul(vr,vl);
!       do i=1,n;
!         vr(i,i) =vr(i,i) -1.0;
!         do j=1,n;
!           if (real(vr(i,j)) > 1.d-6 .or. aimag(vr(i,j)) > 1.d-6) then;
!             write (6,'("error: ",a," (",2i4," =",2g13.5)')            &
!    &           "zimphi*inv",i,j,vr(i,j);
!           endif;
!         enddo;
!       enddo;

        !** transfer function is hx =h*(zI-phi)^(-1)
        !** y power spectrum is (hx+bw)*Q*(hx+bw)^*T +Rn +(hx*g+ac)*|F(u)|^2)*(hx*g+ac)^*T

        hx(:ny,:) =matmul(h(:ny,:n),vl(:,:));
        hxc(:n,:ny) =transpose(conjg(hx(:ny,:n)) +bw(:ny,:n));
        d1(:ny,:n) =matmul(hx(:ny,:n)+bw(:ny,:n),q(:n,:n));
        pw(:ny,:ny) =matmul(d1(:,:),hxc(:n,:ny)) +rn(:,:);
        forall (j=1:nu)                                                 &
     &     pv(:ny,j) =(matmul(hx(:ny,:n),gu(:n,j)) +ac(:ny,j)) *uz(j);
        do j=1,ny;
          ypsd(j) =sqrt(real(pw(j,j))**2+aimag(pw(j,j))**2);   !imaginary component should be 0
          if (nu > 0 .and. mu > 0)                                      &
     &             ypsd(j) =ypsd(j) +sum(pv(j,:)*conjg(pv(j,:)));
          ypsd(j) =2.*dt*ypsd(j);   !factor 2 accounts for negative frequencies (should use 1.0 for ip=0, but this gives discontinuity)
          psum(j) =psum(j) +ypsd(j);
          ypsd(j) =10.*log10(ypsd(j));
        enddo;
        write (10,'(6g13.5)') ip*delf, ypsd(:ny);
      enddo;

      write (6,'("plotPSD2 summed power =",10g13.5)') psum(:)*delf;

      if (nu > 0) write (6,'("average U power =",4g13.5)')              &
     &                       pusum(:)*delf;
      close(10);
      return;
      end subroutine plotPSD2;
