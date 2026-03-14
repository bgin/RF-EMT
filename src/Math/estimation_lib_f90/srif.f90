!*************************************************************************
      subroutine srif (xsig,  xfs,rsrif,                                &
     &                phit,qt,ht,y,sigr,t,m,mdim,n,ndim,ndyn,sparse);

!** computes one time update and measurements update of the Square-Root Information Filter (SRIF).
!**      Author B. Gibbs, 12/2009

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
      use filterPrec, only: real_kind1;   !precision (4 or 8)
      use filterMod, only:first,nstate,mmaxi,xstate,rsum,debugp;
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
      real(real_kind1),intent(inout) :: xfs(n);    !filter state vector at last time step
      real(real_kind1),intent(inout) :: rsrif(n*(n+3)/2);!square root information matrix (upper triangular) and vector (last column)
      real(real_kind1),intent(out) :: xsig(n);     !output filter 1-sigma uncertainty from covariance matrix

      !** local
      integer(4) i,j,k,ier,ic,nbias;
      integer(4),save :: ns2,nnd;
      real(real_kind1) c,rsos,cnb;
      real(real_kind1),allocatable,save :: xa(:),phiInv(:,:),d(:,:);
      real(real_kind1),allocatable,save :: dy(:),qsym(:),qf(:,:);
      real(real_kind1),allocatable,save :: w(:,:),hc(:,:),v(:),g(:);
      real(real_kind1),allocatable,save :: rphi(:,:),rinv(:),dx(:);
      real(real_kind1),allocatable,save :: covout(:);
      !_______________________________

      if (ndim /=nstate) call errhand ("tst_KF error: ndim /= nstate");
      if (n /= nstate) call errhand ("tst_KF error: n /= nstate");
      if (mdim > mmaxi) call errhand ("tst_KF error: mdim > mmaxi");

      if (first(3)) then;
        first(3) =.false.;
        ns2 =n*(n+1)/2;
        nnd =n+ndyn;
        if (.not.allocated(xa)) then;
          allocate (xa(n),phiInv(n,n),d(n,n),dx(n),covout(n*(n+1)/2),   &
     &       rphi(n,n),dy(mmaxi),g(n),qsym(ns2),qf(n,ndyn),             &
     &       v(nnd),w(nnd,nnd+1),hc(mdim,n+1),rinv(n*(n+1)/2),stat=ier);
          if (ier /= 0) call errhand ("tst_udKF error: allocate");
        endif;
        rsrif(:) =0.d0;
        do i=1,n;
          rsrif(i*(i+1)/2) =1.d0/xstate(i)%sig0;
        enddo;
        xfs(:) =0.d0;
        rsrif(ns2+1:ns2+n) =xfs(:n)/xstate(:n)%sig0;
        if (debugp) call prtmx (phit,ndim,n,n,"phit1");
      endif;

!** srif time update
      w(:,:) =0.d0;
      forall (i=1:ndyn) w(i,i) =1.d0;

      qsym(:) =0.d0;
      k =0;
      do i=1,ndyn;
        do j=1,i;
          k =k+1;
          qsym(k) =qt(j,i);
        enddo;
      enddo;
!        call prtsmx (qsym,n,"sqrt-qd")
      !** factor state noise matrix Q as qsyn*qsyn^T
      call cov2ud_m (qsym,ndyn);

      !** scale columns of qsym
      k =1;
      qf(:,:) =0.d0;
      qf(1,1) =sqrt(qsym(1));
      do i=2,ndyn;
        c =sqrt(qsym(k+i));
        do j=1,i-1;
          qf(j,i) =qsym(k+j)*c;
        enddo;
        k =k+i;
        qf(i,i) =c;
      enddo;

      !** compute inv(phi)
      phiInv(:n,:n) =phit(:n,:n);   !copy phi (also sets lower "bias" identity)
      nbias =n-ndyn;
      call matInvGJ (ier,  phiInv,phiInv(:,ndyn+1:n),                   &
     &                       n,ndyn,nbias);
      phiInv(:ndyn,ndyn+1:n) =-phiInv(:ndyn,ndyn+1:n);

      !** check matrix inversion
!      d(:,:) =matmul(phiinv(:n,:n),phit(:n,:n))
!      do i=1,n
!        d(i,i) =d(i,i) -1.d0
!        do j=1,n
!          !** largest errors are 3d-14 in rows 7-9, col 26, 27
!          if (abs(d(i,j)) > 1.d-14) write (*,'("phi*inv(phi)-I: (",i3,
!     &         ",",i3,") =",g13.5)') i,j,d(i,j)
!        enddo
!      enddo

      k=0;
      rphi(:,:) =0.d0;
      do i=1,n;   !r-column loop
        do j=1,i; !r-row loop
          k =k+1;
          do ic=1,n;
            rphi(j,ic) =rphi(j,ic) +rsrif(k)*phiInv(i,ic);
          enddo;
        enddo;
      enddo;
!      call prtmx (phiinv,n,n,n,"phiinv")
!      call prtmx (rphi,n,n,n,"rphi")

      do i=1,n;         !row loop
        do j=1,ndyn;    !
          w(ndyn+i,j) =-dot_product(rphi(i,:j),qf(:j,j));
        enddo;
        w(ndyn+i,ndyn+1:ndyn+n) =rphi(i,:n);
      enddo;
      w(ndyn+1:ndyn+n,nnd+1) =rsrif(ns2+1:ns2+n);

      !** triangularize w: output rsrif is in w(ndyn+1:ndyn+n,ndyn+1:ndyn+n)
      call tdhht_m (w,v,  nnd,nnd,nnd+1, 1,nnd);
      !** copy output rsrif         !### need to write arrays for smoothing
      k=0;
      do i=1,n;
        do j=1,i;
          k=k+1;
          rsrif(k) =w(ndyn+j,ndyn+i);
        enddo;
      enddo;
      rsrif(ns2+1:ns2+n) =w(ndyn+1:ndyn+n,nnd+1);
      xa(:n) =matmul(phit(:n,:n),xfs(:n));

!** srif measurement update
      dy(:m) =y(:m);!-matmul(ht(:m,:n),xa(:n))
      do i=1,m;
        hc(i,:n) =ht(i,:n)/sigr(i);
        hc(i,n+1) =dy(i)/sigr(i);
      enddo;
      rsos =0.d0;
      call thhc_m (rsrif,hc,rsos,  n,mdim,m,1,sparse);
      rsum(3) =rsum(3) +rsos**2;

      !** compute state estimate xfs =inv(rsrif)*v
      v(:n) =rsrif(ns2+1:ns2+n);
      call rinz_m (dx,ier,  rsrif,n,v);
      if (ier /= 0) then;
        write (*,'("srif error: rinz_m singular at row ",i3)') ier;
        call errhand ("srif error: rinz_m singular");
      endif;
      xfs(:n) =dx(:n);
!      xfs(:n) =xa(:n) +dx(:n)

      !** ### check rsrif*xfs -v  !get errors of 1e8 at end, but v is 1d23
!      call rdx (g,  rsrif,dx,n)
!      if (debugp) write (*,'(g13.5,", r*x-z =",3(/10g13.5))')
!     &     t,(g(:n)-v(:n))    !/abs(v(:n))

      !** compute xsig     ###
      call rincon_m (rinv,  cnb,  rsrif,n);     !compute rinv =inv(rsrif)
!      call ri2cov_m (xsig,covout,  rinv,n,n,0)    !compute full covariance
      call ri2cov_m (xsig,covout,  rinv,n,0,n);   !compute only sigmas

      return;
      end subroutine srif;
