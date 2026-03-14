! *********************************************************************
      subroutine optl_back (cost,r,g,nmeas,irflag,  xtot,               &
     &                      nxtot,niter,xmin,xmax,convtst,lprnt_o);

!         Subroutine optl_back computes the state estimate xtot minimizing a
!         residual sum-of-squares cost function.  The optimal state estimate
!         is obtained using Gauss-Newton iterations with back-tracking to handle
!         nonlinearities.  It also handles minimum/maximum state constraints

!         A user-supplied function
!              funct (cost,pf,r,g,delx,nmeas,   partial,xtot,nxtot,debugp)

!         must be provided where (partial,xtot,nxtot,debugp) are inputs and
!              (cost,pf,r,g,delx,nmeas) are outputs.
!            real(8) Pf(nxtot*(nxtot+1)/2) is the upper triangular by columns
!                    portion of the covariance matrix,
!            logical(4) partial denotes whether the cost partial derivatives are to be computed,
!            real(8) delx(nxtot) is the computed Gauss-Newton step,
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
      use optlDat, only:fstat,srif;
      implicit none;

      integer(4),intent(in) :: nxtot;              !total number of states
      integer(4),intent(in) :: lprnt_o;            !debug print flag for optimization: 0=none, 1=level 1, 2=max
      integer(4),intent(in) :: niter;              !maximum number of newton iterations
      real(8),intent(in) :: xmin(nxtot);           !minimum allowable values for xtot
      real(8),intent(in) :: xmax(nxtot);           !maximum allowable values for xtot
      real(8),intent(in) :: convtst;               !convergence test (e.g. 0.05)
      real(8),intent(inout) :: xtot(nxtot);        !initial and final estimate of all states
      integer(4),intent(out) :: irflag;            !flag indicating convergence status
                                                   !  0 - successful convergence
                                                   !  1 - the maximum iteration count was reached
                                                   !  2 - more than 4 back-trackings required
      integer(4),intent(out) :: nmeas;             !number of measurements processed
      real(8),intent(out) :: g(nxtot);             !gradient of cost wrt xtot
      real(8),intent(out) :: cost;                 !least squares cost function (also fstat%cost)
      real(8),intent(out) :: r(nxtot*(nxtot+3)/2); !srif information array

      !***  local variables
      integer(4) i,npnp1,ncons,jback,nback,npnp2,nhalf,iter;
      logical(4) partial,lcons;
      integer icons(nxtot);
      real(8) z(nxtot),sigx(nxtot),xl(nxtot),xnewt(nxtot);
      real(8) delx(nxtot),scalenewton;
      real(8) pf(nxtot*(nxtot+1)/2),rtmp(nxtot*(nxtot+3)/2);
      real(8) cost0,costl,cost1,cost2,step,step1,step2,testc,djds;
      !________________________________________________________

      costl =1.d20;
      irflag =-1;
      npnp2 =(nxtot*(nxtot+3))/2;
      npnp1 =(nxtot*(nxtot+1))/2;
      jback =0;

!****   iteration loop
      do iter=1,niter;
        fstat%iter =iter;
        nmeas =0;
        nback =0;
        nhalf =0;
        xl(:) =xtot(:);      !save last state used to compute grad, fis

        write (6,1010) iter;
 1010   format(//1x,6('**********')/                                    &
     &   ' ITERATION ',i3/' CALLING FUNCTION WITH PARTIALS')

        !**  evaluate cost function and (optionally) srif R-matrix)
        partial=.true.;
        call funct (cost,pf,r,g,delx,nmeas,   partial,xtot,nxtot,.true.);

        cost0 =fstat%cost;
        if (iter == 1) costl =fstat%cost;

        !**  compute fraction of newton step based on residual rms
        scalenewton =20.;    !3.0
!        step=min(scalenewton/sqrt(fstat%rrms),1.d0)  !use scalenewton=3
        step=1.d0;
        write (6,1400) iter,fstat%testc,fstat%cost;
        write (0,1400) iter,fstat%testc,fstat%cost;
!        write (1,1400) iter,fstat%testc,fstat%cost
 1400   format(/' AT G-N ITERATION ',i3,', CONVERG. TEST =',g15.5,      &
     &        ', COST =',g12.5)
        write (6,'(" NEWTON STEP FACTOR =",F8.3/)') step;
        testc =fstat%testc;

        !**  compute newton point and limited newton step and
        xnewt(:) =xl(:)+delx(:); !1:nxtot
        delx(:) =step*delx(:);
        xtot(:) =xtot(:)+delx(:);
        icons(:) =0;            !default is no constraints

        !_____________________
   10   continue;              !start of back-tracking logic loop

        !***  back-tracking logic
        !***  check for constrained states
        call constrx (lcons,  icons,xtot,  xmin,xmax,nxtot);

        do while (lcons);
          if (lprnt_o > 1) then;
             write (6,'(/" ICONS=",30I4,200(/7X,30I4))') icons(:nxtot);
             call prtmx (xtot,1,1,nxtot,'XTOTCI  ');
          endif;

          !** compute constrained minimum
          !** find the constrained minimum on the surface defined by
          !** the constraint vector.
          z(:) =r(npnp1+1:npnp1+nxtot);
          if (.not.srif)                                                &
     &           stop 'optl_back: cannot use constraints if not srif';
          call consop1 (z, r, icons, nxtot, npnp2, xmax, xmin, xl,      &
     &         xnewt, rtmp,     xtot);

          !** reset delx to 'STEP' times constrained step for
          !**  unconstrained states
          do i=1,nxtot;
            if (icons(i) == 0) then;
              xtot(i) =xl(i)+step*(xtot(i)-xl(i));
!           else     !no change
            endif;
          enddo;

          !**  re-compute convergence test without constrained states
          z(:) =r(npnp1+1:npnp1+nxtot);
          testc =sqrt(sum(z(:)**2)/nxtot);
          write (6,'(/" Constrained convergence test =",g12.5)')        &
     &            testc;

          !**  re-constrain states
          call constrx (lcons,  icons,xtot,  xmin,xmax,nxtot);
        enddo;    !____ end "DO WHILE (LCONS)" loop

        ncons=0;
        do i=1,nxtot;
          if (icons(i) /=0) ncons =ncons+1;
        enddo;
        write (6,'(/" ",I6," STATES WERE CONSTRAINED")') ncons;

        if (lprnt_o > 0) then;
          if (ncons > 0) then;
            do i=1,nxtot;
              if (icons(i) /= 0) then;
                write (6,'(" STATE ",I6," CONSTRAINT FLAG ",I3,         &
     &                ", VALUE =",G13.5)') i,icons(i),xtot(i);
              endif;
            enddo;
          else;
            write (6,'(/" NO STATES WERE CONSTRAINED")');
          endif;
        endif;

        !**  compute dot product of gradient and step
        !**  (xtot-xl may not equal delx because of constraints)
        if (nback == 0 .and. nhalf == 0)                                &
     &      djds =dot_product(xtot(:)-xl(:),g(:));

        if (lprnt_o > 1) call prtmx (xtot,1,1,nxtot,'XTOT    ');

        rtmp(:) =r(:);  !save r

        !** compute cost function at new point
        !** (funct must not set delx or g !!!!!!!!!)
        partial =.false.;
        call funct (cost,pf,r,g,delx,nmeas,   partial,xtot,nxtot,.true.);
        r(:) =rtmp(:);  !restore r

        write (6,1420) iter,nback,step,fstat%cost;
        write (0,1410) iter,nback,step,fstat%cost;
!        write (1,1410) iter,nback,step,fstat%cost
 1410   format(' AT G-N ITERATION ',i3,', BACK-TRACKING #',i2,          &
     &        ' STEP=',f8.5,', COST =',g12.5)
 1420   format(/' AT G-N ITERATION ',i3,', BACK-TRACKING #',i2,         &
     &        ' STEP=',f8.5,', COST =',g12.5)

        !**  test for reduction in cost function
        if (fstat%rrms > 1.d15 .or. fstat%cost > cost0) then;

          if (fstat%rrms > 1.d15)then;
            !** if get invalid model, cut interval in half
            step =0.5d0*step;
            nhalf =nhalf+1;

            if (nhalf > 9) then;
              irflag =2;
              write (0,1440);
!              write (1,1440)
              write (6,1440);
 1440         format(' EXIT OPTIMIZATION: MORE THAN 9 STEP-HALVINGS ',  &
     &                  'DUE TO INVALID MODEL')
              go to 995;
            endif;

          else if (fstat%cost > cost0) then;
            !** cost function is worse
            !** accept step if converged and cost is only slightly high
            if (testc < convtst .and. fstat%cost < 1.01*cost0) then;
              irflag =0;
              write (6,                                                 &
     &         '("converged: accept step even though slightly large")');
              go to 995;
            endif;

            !**  if cost function is worse, perform line minimization
            !**  assuming quadratic or cubic with gradient at xl known.
            !**  (cost= cost0+ djdx*dx+ dj2dx2*dx**2 where dx=1 for full step
            cost1 =cost;
            step1 =step;
!            call stephalf (step1,step2,cost2,nback,
            call backtrack (step1,step2,cost2,nback,                    &
     &                      cost0,cost1,djds,lprnt_o);
            step=step1;
          endif;

          write (6,'(/" BACK-TRACKING NUMBER ",I3)') nback;

          if (nback == 2) then;
            jback =jback+1;
          else if (step < 0.01 .and. jback >= 3) then;
            write (6,'(//" ######### OPTIMIZATION IS STUCK: "/          &
     &        " ACCEPT STEP AND ALLOW INCREASE IN COST FUNCTION")');
            jback =0;
            go to 20;
          endif;

          if (nback > 5) then;
            irflag =2;
            write (0,1430);
!            write (1,1430)
            write (6,1430);
 1430       format(' EXIT OPTIMIZATION: MORE THAN 5 BACK-TRACKS')
            go to 995;
          endif;

          xtot(:) =xl(:)+step*delx(:);   !1:nxtot
          icons(:) =0;

          go to 10;
        endif;

   20   costl =fstat%cost;
        if (nback == 0) jback =0;

        if (testc < convtst) then;
          irflag =0;
          write (6,'("converged")');
          go to 995;
        endif;
      enddo;!_______end of iteration loop

      irflag=1;

  995 continue;   !end of iteration loop
      fstat%testc =testc;

      return;
      end subroutine optl_back;
! *********************************************************************
      subroutine backtrack (step1,step2,cost2,nback,                    &
     &                      cost0,cost1,djds,lprnt_o);

!       backtrack finds the minimum along a line (line minimization)
!       assuming either a quadratic or cubic model given the cost and
!       gradient at the initial point, and the cost at 1 or 2 other points.
!       reference: numerical recipes in fortran, 1992, section 9.7

!       if nback =0 (only know the cost at the initial point and one other)
!          model: cost= cost0+ djds*step+ a*step**2
!                   given cost1 at step1.

!       if nback >0 (know the cost at the initial point and two others)
!          model: cost= cost0+ djds*step+ b*step**2 +a*step**3
!                   given cost1 at step1 and cost2 at step2.

!       inputs:
!         cost0 - cost function at starting point (step=0)
!         cost1 - cost function at step1 (most current value)
!         cost2 - cost function at step2 (only if nback > 0)
!         djds - partial deriv. of cost wrt step evaluated at step=0
!         step1 - current step
!         step2 - previous step if nback > 0. (step2 > step1)
!         nback - number of backtracking steps
!       outputs
!         step1 - updated step for minimum of function
!         step2 - previous step1
!         cost2 - previous cost1
!         nback - updated number of backtracking steps (nback=nback+1)

      implicit none;

      real(8),intent(in) :: cost0,cost1,djds;
      integer(4),intent(in) :: lprnt_o;     !print flag
      integer(4),intent(inout) :: nback;
      real(8),intent(inout) :: step1,step2,cost2;

      !** local
      real(8) a,b,d,dc1,dc2,det,djdx,step;
      !___________________________________________________


      if (nback == 0) then;
        !** only have one cost value available - use quadratic model
        djdx =step1*djds;
        a =cost1-cost0-djdx;
        if (a <= 0.d0) then;                 !should not happen
          if (lprnt_o > 0) write (6,'(" BACKTRACK ERROR: A <= 0.",      &
     &                    " HALVING STEP")');
          step =0.5*step1;
        else;   !compute optimum step for quadratic
          step =-djdx/(2.*a);
        endif;

      else if (nback > 0) then;
        !**  have two cost values available - use cubic model
        dc1 =cost1-djds*step1-cost0;
        dc2 =cost2-djds*step2-cost0;
        det =1.d0/(step1-step2);
        a =det*(dc1/step1**2-dc2/step2**2);
        b =det*(-dc1*step2/step1**2+dc2*step1/step2**2);
        d =b**2-3.*a*djds;
        if (d < 0.d0) then;
          if (lprnt_o > 0)                                              &
     &       write (6,'(" BACKTRACK ERROR: B**2-3*A*DJDS < 0: ",        &
     &           "  USING QUADRATIC MODEL")');
          a =cost1-cost0-djds*step1;
          step =-djds*step1/(2.*a);
        else;
          step =(-b+sqrt(d))/(3.*a);
        endif;

      else;
        stop 'BACKTRACK ERROR: NBACK < 0';
      endif;

      !***  limit step to range 0.20 to 0.70*step1
      step =min(step,0.7d0*step1);
      step =max(step,0.20d0*step1);

      a =cost1-cost0-djds*step1;
      if (lprnt_o > 0) then;
        write (6,1500) cost0,djds,step1,a,step;
 1500   format(/' Back-tracking along (constrained ?) Newton ',         &
     &   'direction.'/' Initial cost (c0) at step=o: ',g13.5/           &
     &   ' gradient wrt step (dcds) at step=o: ',g13.5/                 &
     &   ' previous step (step1)=',f8.5/                                &
     &   ' difference in cost from linear model: c1-(c0+dcds*step1) =', &
     &   g13.5/' computed new line step =',f8.5/)
      endif;

      !**  shift cost, step and back-tracking count
      cost2 =cost1;
      step2 =step1;
      step1 =step;
      nback =nback+1;

      return;
      end subroutine backtrack;
! ****************************************************************
      subroutine stephalf (step1,step2,cost2,nback,                     &
     &                      cost0,cost1,djds,lprnt_o);

!       stephalf cuts the current step in half

!       inputs:
!         cost0 - cost function at starting point (step=0)
!         cost1 - cost function at step1 (most current value)
!         cost2 - cost function at step2 (only if nback > 0)
!         djds - partial deriv. of cost wrt step evaluated at step=0
!         step1 - current step
!         step2 - previous step if nback > 0. (step2 > step1)
!         nback - number of backtracking steps
!       outputs
!         step1 - updated step for minimum of function
!         step2 - previous step1
!         cost2 - previous cost1
!         nback - updated number of backtracking steps (nback=nback+1)

      implicit none;

      real(8),intent(in) :: cost0,cost1,djds;
      integer(4),intent(in) :: lprnt_o;     !print flag
      integer(4),intent(inout) :: nback;
      real(8),intent(inout) :: step1,step2,cost2;

      !** local
      !___________________________________________________

      step2 =step1;
      step1 =0.5d0*step1;
      cost2 =cost1;
      nback =nback+1;

      end subroutine stephalf;
! ****************************************************************
      subroutine constrx (lcons,  icons,x,  xmin,xmax,n);
!          constrx constrains x to be within xmin to xmax and sets
!          icons for constrained states. lcons is set to true if any
!          states are constrained
      implicit none;

      integer(4),intent(in) :: n;
      real(8),intent(in) :: xmin(n),xmax(n);

      integer(4),intent(inout) :: icons(n);
      real(8),intent(inout) :: x(n);
      logical(4),intent(out) :: lcons;

      !** local
      integer(4) i;
      !________________________________________________


      lcons=.false.;
      do i=1,n;
!        icons(i)=0        !not reset !!!!!!!!!!!!!!!!

        if (x(i).gt.xmax(i)) then;
          icons(i)=1;
          x(i)=xmax(i);
          lcons=.true.;
        endif;

        if (x(i).lt.xmin(i)) then;
          icons(i)=-1;
          x(i)=xmin(i);
          lcons=.true.;
        endif;

      enddo;

      return;
      end subroutine constrx;
!***********************************************************************
!                                                        *
!                            consop1                     *
!                                                        *
!*********************************************************
         subroutine consop1                                             &
     &                   (z     , r     , icons  , n , npnp2,           &
     &                    tmax  , tmin  , x, xnewt,                     &
     &                    rtmp, y     );

!   purpose
!        finds the minimum of the quadratic function on the surface
!        defined by the constraints.
!   description
!        consop1 checks the constraint vector icons to determine if any
!        constraints exist. if not, the minimum of the quadratic func-
!        tion is set to the computed newton point. if constraints exist,
!        the constrained minimum of the quadratic function is computed.

!   inputs:
!                  z     (r) = square-root information vector
!                  r     (r) = square-root information matrix
!                  icons (i) = constraint vector
!                  n     (i) = number of unknown parmeters
!                  npnp2 (i) = (n*(n+3))/2
!                  tmax  (r) = maximum allowed value of state
!                  tmin  (r) = minimum allowed value of state
!                  x     (r) = parameter estimates at which r,z were computed
!                  xnewt (r) = newton point
!   input/output
!                  y     (r) = newton point or constraint minimum of
!                              the quadratic function.  in input,
!                              y must be the current constrained estimate.

!   associated subprograms
!        called by:  optl_back
!        calls to:   thhcs, rinzc, errhand
!   language:   fortran 90

      implicit real(8) (a-h,o-z);

      integer,intent(in) ::  n, npnp2, icons(n);
      real(8),intent(in) :: z(n), r(npnp2),tmax(n),tmin(n),x(n),xnewt(n);
      real(8),intent(out) :: y(n), rtmp(npnp2);

!            local variables
      integer(4)  i , j,  l, iflag, map(1);
      real(8)  dx(n), c(n+1);
      !__________________________________________________


      nnp1 =(n*(n+1))/2;
      if (npnp2 < nnp1+n) stop 'CONSOP1 ERROR: NPNP2 < N*(N+3)/2';

      !** determine if any constraints exist.
      iflag = 0;
      do i = 1, n;
        if (icons(i) /= 0) iflag = 1;
      end do;

      if (iflag == 0) then;
          !** constraints do not exist. set minimum of the quadratic
          !**function to the computed newton point.
          y(:) = xnewt(:);
      else;          !constraints do exist.
          !**  copy r|z to rtmp
          rtmp(1:nnp1) =r(1:nnp1);
          rtmp(nnp1+1:nnp1+n) = z(:);
!          call prtmx (z,1,1,n,'Z       ');

          !**  first set constrained states to boundary, adjust all
          !**  elements of the z vector for contrained states, and zero column.
          do i=1,n;
            dx(i) =0.d0;
            if (icons(i) /= 0) then;
              if (icons(i) > 0) then;
                  dx(i) = tmax(i)-x(i);
              else if (icons(i) < 0) then;
                  dx(i) = tmin(i)-x(i);
              endif;

              ! print *,i,n,npnp2,dx(i)
              l=(i*(i+1))/2-i;
              do j=1,i;
                rtmp(nnp1+j) =rtmp(nnp1+j)-rtmp(l+j)*dx(i);
                rtmp(l+j) =0.d0;
              enddo;
            endif;
          enddo;

          !**    next delete all constrained rows and re-triangularize
          rsos=0.;
          do i=1,n;
            if (icons(i) /= 0) then;
              ! print *,i,n,npnp2,dx(i)
              c(1:i) =0.d0;

              l=(i*(i+1))/2+i;
              do j=i+1,n+1;
                c(j) =rtmp(l);
                rtmp(l) =0.d0;
                l =l+j;
              enddo;
              nstrt =i+1;
              call thhc_m (rtmp,c,rsos,  n,1,1,nstrt,.false.); !was thhcs (rtmp,n,c,1,1,rsos,nstrt)
              ! call prtsmx (rtmp,n+1,'RTMP   ')
            endif;
          enddo;

          !**  backsolve r*dx=z for dx,
          !**  bypassing rows with constrained states
          call rinzc (rtmp,rtmp(nnp1+1),dx,n,icons);

          do i = 1, n;
            if (icons(i) > 0) then;
              !**  set y exactly to boundary for constrained states
              y(i) = tmax(i);
            else if (icons(i) < 0) then;
              y(i) = tmin(i);
            else;
              y(i) = x(i) + dx(i);
            endif;
          end do;
      end if;

      return;
      end subroutine consop1;
!*********************************************************************
!                                                        *
!                  rinzc                                 *
!                                                        *
!*********************************************************

      subroutine rinzc (r,z,dx,n,icons);
!        rinz backsolves r*dx=z for dx where r is the upper
!        triangular square-root information matrix.  constrained
!        states are bypassed and must be set on input !!!!!!!!!!!!!


!   communications
!        parameter list
!             given arguments
!                  z(n)     (n) = srif information vector
!                  n        (i) = state dimension
!                  r(*)     (r) = srif r matrix
!                  icons(n) (i) = constraint vector
!                  dx(n)    (r) = on input, dx must be set to y(i)-x(i) for
!                                 constrained states
!             both arguments
!                  dx(n)    (r) = computed solution

      implicit real*8 (a-h,o-z);
      dimension r(*),dx(n),z(n),icons(n);


!        backsolve r*dx=y for dx, bypassing constrained parameters
      l1=(n*(n+1))/2;

      do i=n,1,-1;
        if (icons(i).eq.0) then;
           sum=z(i);
           l=l1+i;
           if (i.lt.n) then;
              do j=n,i+1,-1;
                l=l-j;
                sum=sum-dx(j)*r(l);
              enddo;
           endif;
           l=l-i;
           dx(i)=sum/r(l);
        else;
!              set dx as step to boundary  (must be set in input)
!          dx(i)=y(i)-x(i)
        endif;
      enddo;

      return;
      end subroutine rinzc;
