!*****************************************************
      subroutine udKF (xsig,  xfu,ud,                                   &
     &                phit,qt,ht,y,sigr,t,m,mdim,n,ndim,ndyn,sparse);
!** computes one time update and measurements update of the U-D filter.
!**   The state noise covariance Q and state covariance P matrices are stored as square
!**   for ease of use, even though only the upper triangular partition is used.
!**   The subroutine uses variables from two modules: filterPrec and filterMod.
!**   The precision of the filter can be set using variable real_kind1 in filterPrec.

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
      use filterPrec, only: real_kind1;   !precision (4 or 8)
      use filterMod, only:first,KFmethod,nstate,mmaxi,xstate,rsum,      &
     &                    debugp,method_udtime;
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
      real(real_kind1),intent(inout) :: xfu(n);    !filter state vector at last time step
      real(real_kind1),intent(inout) :: ud(n*(n+1)/2);!square root information matrix (upper triangular) and vector (last column)
      real(real_kind1),intent(out) :: xsig(n);     !output filter 1-sigma uncertainty from covariance matrix

      integer(4) i,j,k,ier,ns2,nnd;
      real(real_kind1) z,Pz,dz,rvar,c;
      real(real_kind1),allocatable,save :: uda(:),xa(:),phiua(:,:);
      real(real_kind1),allocatable,save :: g(:),dy(:),pfu(:),psym(:);
      real(real_kind1),allocatable,save :: d(:),w(:,:),dw(:),v(:),vd(:);
      real(real_kind1),allocatable,save :: f(:),h(:);
      !_______________________________

      if (ndim /=nstate) call errhand ("tst_KF error: ndim /= nstate");
      if (n /= nstate) call errhand ("tst_KF error: n /= nstate");
      if (mdim > mmaxi) call errhand ("tst_KF error: mdim > mmaxi");

      if (first(2)) then;
        first(2) =.false.;
        ns2 =n*(n+1)/2;
        if (.not.allocated(uda)) then;
          allocate (uda(ns2+n),pfu(ns2),xa(n),g(n),phiua(ndim,n),       &
     &       dy(mmaxi),d(n),psym(ns2),f(n),h(n),                        &
     &       v(2*n),vd(n),w(n,2*n),dw(2*n),stat=ier);
          if (ier /= 0) call errhand ("tst_udKF error: allocate");
        endif;
        ud(:) =0.d0;
        do i=1,n;
          ud(i*(i+1)/2) =xstate(i)%sig0**2;
        enddo;
        xfu(:) =0.d0;
      endif;

!     write (6,'("starting U-D time update")');
      !** UD filter time update              !### recode for biases
      call phiu_m (phiua,  phit,ndim,n,n,ud,n,n);   !*** only do for ndyn rows
      k=0;
      do i=1,n;
        k =k+i;
        d(i) =ud(k);
      enddo;

      w(:n,:n) =phiua(:n,:n);
      dw(:n) =d(:n);
      psym(:) =0.d0;
      k =0;
      do i=1,ndyn;
        do j=1,i;
          k =k+1;
          psym(k) =qt(j,i);
        enddo;
      enddo;
!        call prtsmx (psym,n,"sqrt-qd")
      call cov2ud_m (psym,ndyn);
      nnd =n+ndyn;

      if (method_udtime == 1) then;
        !** modified weighted gram-schmidt time update   !### recode for biases
        w(:n,n+1:2*n) =0.d0;
        k =0;
        do i=1,ndyn;
          do j=1,i;
            k =k+1;
            w(j,i+n) =psym(k);
          enddo;
          w(i,i+n) =1.d0; !replaces diagonal
          dw(i+n) =psym(k);
        enddo;
!        call prtmx (w,n,n,2*n,"w")
!        call prtmx (dw,1,1,2*n,"dw")

        !*** weighted gram-schmidt orthogonalization
        call wgs_m (ud,v,  w,  n,n,nnd,dw,sparse);

      else;
        !*** weighted gram-schmidt orthogonalization for phiu  !### recode for biases
        !*** followed by rank-1 updates for q
        call wgs_m (ud,v,  w,  n,n,n,dw,sparse);
        !** rank1 updates for process noise
        k =0;
        do i=1,ndyn;
          v(:i) =psym(k+1:k+i);
          v(i) =1.d0;
          v(i+1:ndyn) =0.d0;
          c =psym(k+i);
          call rank1_m (ud,ier,  v,  ud,ndyn,c);
          k =k+i;
        enddo;
      endif;

!      call prtmx (w,n,n,nnd,"wgs-w")
!      call prtsmx (ud,n,"wgs-ud")
!      call ud2cov_m (ud,pfu,n)

      xa(:) =matmul(phit(:n,:n),xfu(:n));
      xfu(:) =xa(:);
!      write (6,'("xf diff =",10(/10g13.5))') xfu(:nest) -xfa(:nest)
!      write (6,'("pa diff =",10(/10g13.5))')
!     &      (pfu(i*(i+1)/2)-pfa(i,i),i=1,nest)
                       !### need to write arrays for smoothing

      !** process meas one at a time
      do i=1,m;
        !** UD filter measurement update
        z =y(i);
        rvar =sigr(i)**2;
        h(:n) =ht(i,:n);
        call udmes1 (f,g,Pz,dz,  ud,n,rvar,h,xfu,z);
        dz =y(i) -dot_product(h(:m),xfu(:n));   !recompute dz using double precision y
        rsum(2) =rsum(2) +dz**2/pz;

!        chiSq =dz**2/Pz
!        if (chiSq > editThresh**2) then
!          write (*,1100) t,z,chiSq
!          return
!        endif
        call udmes2 (g,ud,xfu,  n,Rvar,f,dz);
      enddo;

      !** compute covariance
!      call ud2cov_m (pfu,  ud,n)         !this is slow
!      do i=1,n
!        k =(i*(i+1))/2
!        xsig(i) =sqrt(pfu(k))
!      enddo
      call ud2sig_m (xsig,  ud,n,n); !this is faster

      !** write double precision output file for future comparison
      if (debugp) write (40) 2*n+1, t, xfu(:n),xsig(:n);

      return;
      end subroutine udkf;
