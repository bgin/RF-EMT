! ************************************************************************
!
!       **************
!       *            *
!       *   dricss   *
!       *            *
!       **************
!
      subroutine dricss (p,e,t1,s,  ed,  nd2,nd,n,lprnt,delta);
!
!       solve the steady-state discrete algebraic ricatti equation (DARE)
!           0 =-P +Phi*(P -P*H^T*(H*P*H^T +R)^(-1)*H*P)*Phi^T +Q
!
!          ED = | Phi^(-1)                -Phi^(-1)*G*Q*G^T                    |
!               | -H^T*inv(R)*H*Phi^(-1)  Phi^T +H^T*inv(R)*H*Phi^(-1)*G*Q*G^T |

!       Reference: A.N. Beavers and E.D. Denman white paper (date ?)
!       "Steady-State Solution to the Discrete Linear Regulator Filter and Liapunov Equation"

!       Probably similar to E. Denman and A. Beavers,
!       "The matrix sign function and computations in systems",
!       Appl. Math. Comp., 2, pp. 63-94 (1976)

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

!**     Note: The Beavers and Denman test case works correctly with this code using delta =1.d-5
!**      B&D test case
!      write (6,'(/"Start of Beavers and Denman test case")')
!      p(:,:) =0.d0
!      b(:,:) =0.d0
!      q(:,:) =0.d0
!      q(3,3) =0.8d0
!      hrht(1,:) =(/ 1.d0/0.095d0, 1.d0/0.095d0, 0.d0 /)
!      hrht(2,:) =(/ 1.d0/0.095d0, 2.d0/0.095d0, 0.d0 /)
!      hrht(3,:) =0.d0
!      phi(1,:) =(/ 0.036d0, 3.0d0, 0.d0 /)
!      phi(2,:) =(/ 10.d0, 6.964d0, 0.1d0 /)
!      phi(3,:) =(/ 50.d0, -7.278d0, -0.001d0 /)
!      d1(:,:) =phi(:,:)
!      d2(:,:) =q(:,:)

!      !**  compute d1 =inv(phi), d2 =inv(phi)*q
!      call matInvGJ (ier,  d1,d2,  nd,nd,nd)
!      if (ier /= 0) stop 'gaussj error on phi'

!      b(:nd,:nd) =d1(:,:)
!      b(:nd,nd1:nd2) =-d2(:,:)
!      b(nd1:nd2,:nd) =-matmul(hrht(:,:),d1)
!      b(nd1:nd2,nd1:nd2) =transpose(phi(:,:)) +matmul(hrht,d2)

!      delta =1.d-5
!      call dricss (p,e,t1,s,  b,  nd2,nd,nd,1,delta)
!      call prtmx (p,nd,nd,nd,"p-dricss")

      implicit none;
      integer(4),intent(in) :: nd2,nd,n,lprnt;
      real(8),intent(in) :: delta;              !eigenvalue shift (e.g. 1d-5)
      real(8),intent(inout) :: ed(nd2,2*n);        !input array =[
      real(8),intent(out) :: e(nd2,2*n),t1(nd2,2*n),s(nd2,2*n); !temporary arrays
      real(8),intent(out) :: p(nd,n);          !steady-state covariance

      integer(4) nr(n),n2,i,j,k,iter,ier;
      real(8) s1,s2,sume,sold,d2(nd2);
      !________________________________________


      n2 =2*n;
      if (nd < n) stop 'dricss error: nd < n';
      if (nd2 < n2) stop 'dricss error: nd2 < n2';

!**      apply eigenvalue shift
      do i=1,n;
        ed(i,i) =ed(i,i) +delta;
        ed(i+n,i+n) =ed(i+n,i+n) -delta;
      enddo;
!      if (lprnt > 0) call prtmx(ed,nd2,n2,n2,'ed      ')

!**      t1 =ed+i
      t1(:,:) =ed(:,:);
      forall (i=1:nd2) t1(i,i) =t1(i,i) +1.d0;
      if (lprnt > 0) call prtmx(t1,nd2,n2,n2,'t1      ');

!**      invert t1
      d2(:) =0.d0;
!      call gaussj (ier,t1,n2,nd2,d2,1,1)
      call matInvGJ (ier,  t1,d2,  nd2,n2,1);
      if (ier /= 0) stop 'matInvGJ error: t1 inversion failed';

!**      e =(ed-i)*inv(ed+i)
      s1 =0.d0;
      do i=1,n2;
        do j=1,n2;
          sume =-t1(j,i) +dot_product(ed(j,:n2),t1(:n2,i));
          s1 =s1+sume**2;
          e(j,i) =sume;
          s(j,i) =sume;
        enddo;
      enddo;

      if (lprnt > 0) call prtmx(e,nd2,n2,n2,'e       ');
      s(:,:) =e(:,:);
      d2(:) =0.d0;

      !** iteration loop
      do iter =1,25;
        !**  invert e
        call matInvGJ (ier,  e,d2,  nd2,n2,1);
        if (ier /= 0) stop 'matInvGJ error: e inversion failed';
        s1 =0.d0;
        s2 =0.d0;
        do i=1,n2;
          do j=1,n2;
            sold =s(j,i);
            s(j,i) =0.5d0*(s(j,i)+e(j,i));
            if (abs(s(j,i)) < 1.d-15) s(j,i)=0.d0;
            e(j,i) =s(j,i);
            s2 =s2 +(s(j,i)-sold)**2;
            s1 =s1 +s(j,i)**2;
          enddo;
        enddo;

        !**  test for convergence
!       if(lprnt > 0) call prtmx(s,nd2,n2,n2,'s       ')
        if (s2/s1 < 1.d-20) exit;
        if (iter >= 25) stop 'dricss failed to converge';
      enddo;

      if (lprnt > 0) then;
        write(6,100) iter;
  100   format('# iter =',i5)
        call prtmx(s,nd2,n2,n2,'dricss s');
      endif;

!**      add jr
      do i=1,n;
        s(i,i) =s(i,i) +1.d0; !sign(1.d0,s(i,i))
        s(i+n,i+n) =s(i+n,i+n) -1.d0;  !+sign(1.d0,s(i+n,i+n))
      enddo;

!**      invert s+jr
      call matInvGJ (ier,  s,d2,  nd2,n2,1);
      if (ier /= 0) stop 'matInvGJ error: s inversion failed';
      if (lprnt > 0) call prtmx(s,nd2,n2,n2,'(s+jr)-1');

      t1(:n,:n) =2.d0*s(n+1:n2,:n);
      p(:n,:n) =-2.d0*s(:n,n+1:n2);

      if (lprnt == 0) return;
      call prtmx(t1,nd2,n,n,'t1      ');
      call prtmx(p,nd,n,n,'p       ');

      return;
      end subroutine dricss;
