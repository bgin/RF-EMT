!*****************************************************
      subroutine KF (xsig,  xf,pf,                                      &
     &               phit,qt,ht,y,sigr,t,m,mdim,n,ndim,ndyn,sparse);

!**   Subroutine kf is a general-purpose Kalman filter that processes a single time
!**   and measurement update in a single call.  The state noise covariance Q and
!**   state covariance P matrices are stored as square for ease of use, even though
!**   only the upper triangular partition is used.  [The computed P has the lower
!**   triangle set equal to the upper triangle).  The subroutine uses variables from two
!**   modules: filterPrec and filterMod.  The precision of the filter can be set using
!**   variable real_kind1 in filterPrec. Other variables in filterMod control options
!**   or supply array storage.  For example, variable KFmethod determines whether the
!**   filter measurement update is the short Kaplan form, Bierman's version of the
!**   Joseph update, or the full Joseph update.

!** Author: B. Gibbs, 12/2009

!***************************************
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
      use filterPrec, only: real_kind1;   !precision (4 or 8)
      use filterMod, only:first,KFmethod,nstate,mmaxi,xstate,rsum,      &
     &                    debugp,pfa,xfa;
      implicit none;
      integer(4),intent(in) :: m;                  !number of measurements in vector y
      integer(4),intent(in) :: mdim;               !row dimension of matrix h (mdim >= m)
      integer(4),intent(in) :: n;                  !number of states in xfs
      integer(4),intent(in) :: ndim;               !row dimension of phit and qt (ndim >= n)
      integer(4),intent(in) :: ndyn;               !number of dynamic (non-bias) states (ndyn <= n)
      logical(4),intent(in) :: sparse;             !if true, the hc matrix input to thhc_m is assumed to
                                                   ! be sparse and thus sparse logic is used in processing
      real(real_kind1),intent(in) :: t;            !time of current measurement y
      real(real_kind1),intent(in) :: phit(ndim,n); !state transition matrix from last measurement to t
      real(real_kind1),intent(in) :: qt(ndim,ndyn);!state process noise matrix from last measurement to t
      real(real_kind1),intent(in) :: sigr(m);      !1-sigma of measurement noise in y
      real(real_kind1),intent(in) :: y(m);         !measurement vector at time t
      real(real_kind1),intent(in) :: ht(mdim,n);   !measurement partial derivatives: h =del-y/del-x
      real(real_kind1),intent(inout) :: xf(n);     !filter state vector at last time step
      real(real_kind1),intent(inout) :: pf(n,n);   !covariance matrix (symmetric)
      real(real_kind1),intent(out) :: xsig(n);     !output filter 1-sigma uncertainty from covariance matrix

      integer(4),save :: i,j,k,ier,nf;
      real(real_kind1) eps,epsm,h(mdim,n),xa(n),vmin;
      real(real_kind1),allocatable,save :: phip(:,:);    !pa(:,:),xa(:),
      real(real_kind1),allocatable,save :: g(:,:),hp(:,:),dy(:);
      real(real_kind1),allocatable,save :: py(:),py2(:,:),a(:,:),b(:,:);
      !_______________________________

      if (ndim /=nstate) call errhand ("KF error: ndim /= nstate");
      if (n /= nstate) call errhand ("KF error: n /= nstate");
      if (mdim > mmaxi) call errhand ("KF error: mdim > mmaxi");

      if (first(1)) then;
        first(1) =.false.;
        if (.not.allocated(g)) then;
          allocate (g(n,mmaxi),phip(n,n),                               & !pa(n,n),xa(n),
     &       hp(mmaxi,n),dy(mmaxi),py(mmaxi*(mmaxi+1)/2),               &
     &       py2(mmaxi,mmaxi),a(n,n),b(n,mmaxi),                        &
     &       stat=ier);
          if (ier /= 0) call errhand ("tst_KF error: allocate");
        endif;
        pf(:,:) =0.d0;
        forall (i=1:n) pf(i,i) =xstate(i)%sig0**2;
        xf(:) =0.d0;
        nf =n;
!        do i=1,n
!          if (xstate(i)%status >= 2) nf =nf+1
!        enddo
      endif;

      !** filter time update
      phip(:ndyn,:nf) =matmul(phit(:ndyn,:nf),pf(:nf,:nf));
      pfa(:ndyn,:ndyn) =                                                &
     &   matmul(phip(:ndyn,:nf),transpose(phit(:ndyn,:nf))) +qt(:nf,:nf);
      pfa(:ndyn,ndyn+1:nf) =phip(:ndyn,ndyn+1:nf);
      pfa(ndyn+1:nf,:ndyn) =transpose(phip(:ndyn,ndyn+1:nf));
      pfa(ndyn+1:nf,ndyn+1:nf) =pf(ndyn+1:nf,ndyn+1:nf);
      xfa(:ndyn) =matmul(phit(:ndyn,:nf),xf(:nf));
      xfa(ndyn+1:nf) =xf(ndyn+1:nf);

    !** filter time update
!      nf1 =ubound(pae,1)
!      write (6,'("pa diff =",10(/10g13.5))') (pa(i,i)-pae(i,i),i=1,nf1)
!      write (6,'("xa diff =",10(/10g13.5))') (xa(i)-xfa(i),i=1,nf1)
      h(:,:) =ht(:,:);
      hp(:m,:nf) =matmul(h(:m,:nf),pfa(:nf,:nf));
!     call prtmx (hp,mmaxi,m,nf,"hp_tst");

                           !### need to write arrays for smoothing

      !** form residual covariance
      k =0;
      do i=1,m;
        do j=1,i;
          k=k+1;
          py(k) =dot_product(hp(i,:nf),h(j,:nf));
        enddo;
        py(k) =py(k) +sigr(i)**2;
      enddo;
!     call prtsmx (py,m,"py");

      !** invert residual covariance
      eps =1.d-4;
      call sinv (ier,epsm,  py,  m,eps,.true.);
      if (ier /= 0) call errhand ("KF error: sinv failed");
      k =0;
      do i=1,m;
        do j=1,i;
          k=k+1;
          py2(i,j) =py(k);
          py2(j,i) =py(k);
        enddo;
      enddo;

      g(:nf,:m) =matmul(transpose(hp(:m,:nf)),py2(:m,:m));
!      do i=1,nf
!        if (xstate(i)%status == 2) g(i,:m) =0.d0      !zero gain rows for consider states
!      enddo
      xa(:) =xfa(:);          !double precision
      dy(:m) =y(:m) -matmul(ht(:m,:nf),xa(:nf));   !form residual in double precision
      rsum(1) =rsum(1) +dot_product(matmul(dy(:m),py2(:m,:m)),dy(:m));
      xf(:nf) =xa(:nf) +matmul(g(:nf,:m),dy(:m));

      if (KFmethod == 1) then;
        !** update state covariance partitions
        !** Pf =(I-G*H)*Pa
        pf(:nf,:nf) =pfa(:nf,:nf) -matmul(g(:nf,:m),hp(:m,:nf));   !doesn't match filt_meas_update (diff states 7-9 are 0.8)

      else if (KFmethod == 2) then;
        !** bierman's version of joseph
        a(:nf,:nf) =pfa(:nf,:nf) -matmul(g(:nf,:m),hp(:m,:nf));
        b(:nf,:m) =-matmul(a(:nf,:nf),transpose(ht(:m,:nf)));
!        forall (j=1:m) b(:nf,j) =b(:nf,j) +g(:nf,j)*sigr(j)**2
        pf(:nf,:nf) =a(:nf,:nf) +matmul(b(:nf,:m),transpose(g(:nf,:m))); !matches filt_meas_update within 1e-7
        forall (j=1:m) b(:nf,j) =g(:nf,j)*sigr(j)**2;
        pf(:nf,:nf) =pf(:nf,:nf) +matmul(b(:nf,:m),transpose(g(:nf,:m)));

      else if (KFmethod == 3) then;
        !** full joseph
        !** Pf =(I-G*H)*Pa*(I-G*H)^T +G*Rmeas*G^T
        !**    =A +[-(A*H^T) +G*Rmeas]*G^T    where A=(I-G*H)*Pa
        a(:nf,:nf) = -matmul(g(:nf,:m),ht(:m,:nf));
        forall (i=1:nf) a(i,i) =a(i,i) +1.d0;
        pf(:nf,:nf) =matmul(matmul(a(:nf,:nf),pfa(:nf,:nf)),            &
     &              transpose(a(:nf,:nf)));
        forall (j=1:m) b(:nf,j) =g(:nf,j)*sigr(j)**2;
        pf(:nf,:nf) =pf(:nf,:nf) +matmul(b(:nf,:m),transpose(g(:nf,:m)));   !matches filt_meas_update within 1e-7
      endif;

      !** make symmetric: averaging matches filt_meas_update within 1e-6
      !** making symmetric (without avg) gives about the same accuracy: 1e-6
      do i=2,nf;
        do j=1,i-1;
          pf(i,j) =0.5d0*(pf(i,j)+pf(j,i));
          pf(j,i) =pf(i,j);
        enddo;
      enddo;

      !** if using single precision, optionally constrain variances to be greater than vmin  ###
      !** (mainly used for testing different algorithms using IMU data of Chapter 10)
      if (real_kind1 == 4) then;
        vmin =1.e-4;
!        do i=1,n
!          if (pf(i,i) < vmin) then
!            write (*,'("pf(",i2,",",i2,") =",g13.5)') i,i,pf(i,i)
!            if (i <= 6) pf(i,i) =vmin
!          endif
!!          pf(i,i) =max(pf(i,i),1e-4)      !1e-5 for method 3  ###
!        enddo
      endif;
      forall (i=1:n) xsig(i) =sqrt(pf(i,i));


!      write (6,'("pf diff =",10(/10g13.5))') (pf(i,i)-pfe(i,i),i=1,nf1)
!      write (6,'("xf diff =",10(/10g13.5))') (xf(i)-xfe(i),i=1,nf1)
!      write (20,'(10g13.5)') t, 304.8*(xf(:6)-xtr(:6)),xf(7:9)-xtr(7:9)

      return;
      end subroutine kf;
