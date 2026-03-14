! *********************************************************************
!**   description
!         Subroutine optl_LM computes the state estimate xest minimizing a
!         residual sum-of-squares cost function.  The optimal state estimate
!         is obtained using Gauss-Newton iterations with Levinburg-Marquard
!         optimization logic to handle nonlinearities.

!         A user-supplied function
!              funct (cost,pf,ainf,b,dxn,nmeas,  partial,xest,nest,debugp)

!         must be provided where (partial,xest,nest,debugp) are inputs and
!              (cost,pf,ainf,b,dxn,nmeas) are outputs.
!            real(8) ainf(nest*(nest+1)/2) is the upper triangular by columns portion
!                                          of the least squares information matrix
!            real(8) b(nest) is the least squares information vector
!            logical(4) partial denotes whether the cost partial derivatives are to be computed,
!            real(8) dxn(nest) is the computed Gauss-Newton step,
!            logical(4) debugp denoted whether debug print is generated on unit 6

!         Author B. Gibbs  12/2009

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
!----------------------------------------------------------------------
      subroutine optl_LM (cost,nmeas,istat,  xest,pf,                   &
     &                    nest,niter,convtst,lambda0,maxLambda,lprnt_o);

      use optlDat, only:fstat,srif;
      implicit none;

      integer(4),intent(in) :: nest;         !number of states
      integer(4),intent(in) :: lprnt_o;      !debug print flag for optimization
      integer(4),intent(in) :: niter;        !maximum number of newton iterations
      real(8),intent(in) :: convtst;         !required test metric for convergence (e.g, 0.001)
      real(8),intent(in) :: lambda0;         !starting and minimum LM scaling parameter.
                                             ! Normally lambda=1d-4 works well but occasionally need 1d-6
      real(8),intent(in) :: maxLambda;       !maximum allowed value of lambda. normally =100.
      real(8),intent(inout) :: xest(nest);   !state vector
      real(8),intent(inout) :: pf(nest*(nest+1)/2); !upper triangular by columns portion of the covariance matrix
      real(8),intent(out) :: cost;           !least-squares cost function (0.5*SOS)
      integer(4),intent(out) :: nmeas;        !# of scalar measurements
      integer(4),intent(out) :: istat;       !completion status:
                                             ! 0 =converged and OK
                                             ! -1 =not converged
                                             ! -2 =funct failed
                                             ! -3 =inversion of information matrix failed
      !** local
      logical(4) :: partial,lconvg;
      integer(4) i,j,k,ier,iter,iter2,maxIter2,nn2,nstrt;
      real(8) eps,dx(nest),dxn(nest),b(nest),ainf(nest*(nest+3)/2);
      real(8) bL(nest),aL(nest*(nest+3)/2),costL,xestL(nest);
      real(8) h(nest+1),d(nest),rrmsh(niter);
      real(8) lambda,converg,rsos,errm;
      real(8) :: sigx(nest);     !state 1-sigma uncertainty
      !___________________________________________________

      eps =1.d-6;
      costL =1.d20;
      xestL(:) =0.d0;
      iter2 =0;
      lambda =lambda0;
      iter =0;
      istat =0;
      maxIter2 =8;
      nn2 =(nest*(nest+1))/2;
      lconvg =.false.;

      do while (iter < niter);
        iter =iter +1;
        fstat%iter =iter;
        write (6,'(/5("**********")/"Iteration ",i2)') iter;
        partial =.true.;

        !** evaluate least squares cost function, covariance matrix (pf),
        !** gradient (b) and Gauss-Newton step (dxn = -pf*b = -inv(ainf)*b)
        call funct (cost,pf,ainf,b,dxn,nmeas,  partial,xest,nest,.true.);    !funct pf=inv(ainf)
        if (nmeas < 0) then;
          istat =-2;
          write (6,'("optl_LM: least-squares funct failed")');
          return;
        endif;

        xestL(:) =xest(:);
        xest(:) =xest(:) +dxn(:);  !candidate Gauss-Newton step

        converg =sqrt(dot_product(-b,dxn)/nest);
        rrmsh(iter) =fstat%rrms;
        k =0;
        do i=1,nest;
          k =k+i;
          sigx(i) =sqrt(pf(k));
        enddo;

        write (6,'(/"Iteration ",i2,", cost =",g16.8,                   &
     &       ", converg =",g13.5,", xestL =",2(/6g16.8))')              &
     &       iter, cost, converg, xestL(:);
        write (*,'("Iteration ",i2,", cost =",g16.8,"(",i5," meas)",    &
     &       ", converg =",g13.5)') iter, cost, nmeas, converg;

        if (lprnt_o > 0) then;
          write (6,'(/" dxn =",2(/6g16.8))') dxn(:);
          write (6,'(/" sigx =",2(/6g16.8))') sigx(:);
          write (6,'(/"correlation matrix =")');
          k =0;
          do i=1,nest;
            write (6,'(i2,20f7.3)') i,(pf(k+j)/(sigx(i)*sigx(j)),j=1,i);
            k=k+i;
          enddo;
!          if (.not.srif) write (6,'(/"a-diag =",2(/10g13.5))')
!     &                       (ainf(i*(i+1)/2),i=1,nest)
        endif;

        !** step increase (Levinburg-Marquard)
!        lambda = lambda0
        lambda =max(lambda*0.01, lambda0);
        costL =cost;
        aL(:nn2) =ainf(:nn2);
        bL(:) =-b(:);
        iter2 =0;

        do;     !inner loop for Levinburg-Marquard iterations
          partial =.false.;
          !** evaluate least-squares cost function only at candidate step
          call funct (cost,pf,ainf,b,dx,nmeas,                          &
     &                partial,xest,nest,.false.);
          if (nmeas < 0) then;
            istat =-2;
            write (6,'("optl_LM: least-squares funct failed")');
            cost =1.d20;
!            return
          endif;

          if (iter2 == 0) then;
            write (6,1000) iter2, 0., cost;
          else;
            write (6,1000) iter2, lambda,cost;
          endif;
 1000     format(/"Inner iteration ",i2," with lambda =",g13.5,         &
     &         ": at new step cost =",g16.8)

          if (cost < costL) then;
            maxIter2 =5;
            exit;
          endif;

          if (iter2 >= maxIter2 .or. lambda > maxLambda) then;
            if (converg < convtst) then;
              write (6,                                                 &
     &        '(/"Inner iterations have not improved cost function but",&
     &         " last GN step met convergence test: accept last step")');
              xest(:) =xestL(:);   !reset estimate to last accepted value
            endif;
            maxIter2 =maxIter2 -1;    !reduce max inner step for next iteration
            exit;           !accept step anyway
          endif;

          !** step reduction (Levinburg-Marquard)
          lambda =lambda*10.;

          !** compute Newton step: ai =inv(a), step =ai*b
          ainf(:nn2) =aL(:nn2);

          if (srif) then;
            !** compute diagonals of R^T*R where R is the
            !** upper triangular square-root information matrix
            k =0;
            d(:) =0.d0;
            do i=1,nest;
              do j=1,i;
                k=k+1;
                d(i) =d(i) +ainf(k)**2;
              enddo;
            enddo;
            write (6,'("ainf-diag =",10g13.5)') d(:);

!            b(:) =ainf(nn2+1:nn2+nest)      !save RHS
            !** now add one-component-at-a-time lambda*d to diagonals of R^T*R,
            !** and recalculate R factor (all done using orthonormal Householder
            !** transformations for accuracy)
            do i=1,nest;
              h(:) =0.d0;
              h(i) =sqrt(lambda*d(i));
              nstrt =1;
              rsos =0.d0;
              call thhc_m (ainf,h,rsos,  nest,1,1,nstrt,.false.); !was thhcs (ainf,nest,h,1,1,rsos,nstrt)
            enddo;

            !** calculate modified Gauss-Newton step dx =R^(-1)*z where z
            !** is information vector stored in last column of R
            call rinz_m (dx,ier,  ainf,nest,ainf(nn2+1:nn2+nest));
            if (ier < 0) then;
              cost =1.d10;
              istat =-3;
              write (6,'(" optl_LM matrix inversion failed")');
              return;
            endif;

          else;
            !** Levinburg-Marquard using standard least squares normal equations:
            !** add lambda*diagonals to ainf information matrix
            k =0;
            do i=1,nest;
              k=k+i;
              d(i) =ainf(k);
              ainf(k) =ainf(k)*(1.d0+lambda);
            enddo;
            write (6,'("ainf-diag =",10g13.5)') d(:);
            !** compute modified Gauss-Newton step: first invert ainf
            call sinv (ier,errm,  ainf,  nest,eps,.true.);
            if (ier < 0) then;
              cost =1.d10;
              istat =-3;
              write (6,'(" optl_LM matrix inversion failed")');
              return;
            endif;
            call symult (dx,  ainf,bL,nest);  !dx =[ainf^(-1)]*bL using symmetric storage
          endif;
          write (6,'("cost =",g13.5,", lambda =",f9.5,                  &
     &       ", xest =",2(/10g13.5))') cost,lambda,xest(:);
          xest(:) =xestL(:) +dx(:);
          iter2 =iter2 +1;
          write (6,'("step reduction iter ",i2,", cost =",g13.5,        &
     &       ", lambda =",f9.5,", dx =",2(/10g13.5))')                  &
     &       iter2,cost,lambda,dx(:);
!          write (6,'(" xest =",2(/10g13.5))') xest(:)
        enddo;

!        write (6,'("optl_LM iter ",i3,", xest =",6g13.5,", a-diag =",
!     &            6g13.5)') iter,xest(:),(a(i*(i+1)/2),i=1,nest)
        if (converg < convtst) then;
          lconvg =.true.;
          exit;
        else if (iter > 3 .and. converg < 2.*convtst) then;
          if (rrmsh(iter-2)/rrmsh(iter) < 1.01) then;
            lconvg =.true.;
            exit;
          endif;
        endif;
      enddo;   !____ end do while (iter < 20) loop

      if (.not.lconvg) then;
         write (6,'("*** Not converged"/)');
         write (6,'(/"Last iteration cost =",g13.5,                     &
     &       ", converg =",g13.5,", xestL =",2(/6g16.8))')              &
     &        costL, converg, xestL(:);
         write (6,'(/"current cost =",g13.5,", xestL =",2(/6g16.8))')   &
     &        cost, xest(:);
         istat =-1;
      else;
!         write (6,'("converged")')
      endif;

      if (nmeas > 0) then;
        if (partial) write(6,'("nmeas =",i4,", cost =",g13.5,           &
     &            ", rms =",g13.5)') nmeas, cost, fstat%rrms;
      endif;

      return;
      end subroutine optl_LM;
