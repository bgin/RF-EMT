! ************************************************************************
!
!       **************
!       *            *
!       *   cricss   *
!       *            *
!       **************
!
      subroutine cricss (p,e,t1,s,  ec,nd2,nd,n,lprnt);
!
!       solve the steady-state continuous algebraic ricatti equation (CARE)
!          Pdot =G*Q*G^T +F*P +P*F^T -P*H^T*inv(R)*H*P
!
!          EC = | -F              -G*Q*G^T |   (Hamiltonian matrix
!               | -H^T*inv(R)*H    F^T     |

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
      implicit none;
      integer(4),intent(in) :: nd2;
      integer(4),intent(in) :: nd;
      integer(4),intent(in) :: n;
      integer(4),intent(in) :: lprnt;
      real(8),intent(in) :: ec(nd2,2*n);     !input array
      real(8),intent(out) :: e(nd2,2*n);
      real(8),intent(out) :: t1(nd2,2*n);
      real(8),intent(out) :: s(nd2,2*n); !temporary arrays
      real(8),intent(out) :: p(nd,n);          !steady-state covariance

      integer(4) nr(n),n2,i,j,k,iter,ier;
      real(8) s1,s2,sume,sold,d2(nd2);
      !________________________________________

      n2 =2*n;
      if (nd < n) stop 'cricss error: nd < n';
      if (nd2 < n2) stop 'cricss error: nd2 < n2';

      e(:,:) =ec(:,:);
      if (lprnt > 0) call prtmx(e,nd2,n2,n2,'ec');

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
!       if (lprnt > 0) call prtmx(s,nd2,n2,n2,'s')
        if (s2/s1 < 1.d-20) exit;
        if (iter >= 25) stop 'cricss failed to converge';
      enddo;

      if (lprnt > 0) then;
        write(6,100) iter;
  100   format('# iter =',i5)
        call prtmx(s,nd2,n2,n2,'s');
      endif;

!**      add jr
      do i=1,n;
        s(i,i) =s(i,i) +1.d0; !sign(1.d0,s(i,i))
        s(i+n,i+n) =s(i+n,i+n) -1.d0;  !+sign(1.d0,s(i+n,i+n))
      enddo;

!**      invert s+jr
      call matInvGJ (ier,  s,d2,  nd2,n2,1);
      if (ier /= 0) stop 'matInvGJ error: s inversion failed';
      if (lprnt > 0) call prtmx (s,nd2,n2,n2,'(s+jr)-1');

      t1(:n,:n) =2.d0*s(n+1:n2,:n);
      p(:n,:n) =-2.d0*s(:n,n+1:n2);

      if (lprnt == 0) return;
      call prtmx (t1,nd2,n,n,'t1      ');
      call prtmx (p,nd,n,n,'p       ');

      return;
      end subroutine cricss;
