! **********************************************************************
      subroutine dft (c,d,   isn,n,a,b);
!**   discrete Fourier transform (brute-force method, not fast)
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
      integer(4),intent(in) :: n;         !number of points in a,b
      integer(4),intent(in) :: isn;       !+1 =forward transform, -1=inverse
      real(4),intent(in) :: a(n),b(n);    !input real and imaginary data
      real(4),intent(out) :: c(n),d(n);   !transform output real and imaginary parts

      integer(4) i,j;
      real(4) fn,dw,er,ei,cd,sd,dcd,del,csum,ssum;
      real(4) dtx;
      !___________________________________________

      fn =1.0/float(n);
      dw =6.2831853*fn;
      if (isn > 0) dw =-dw;
      dtx =0.;

      do i=1,n;
        er =-2.*sin(dtx/2.)**2;
        ei =sin(dtx);
        cd =1.;
        sd =0.;
        csum =a(1);
        ssum =b(1);
        do j=2,n;
          dcd =cd;
          !** compute sin and cosine by recursion
          cd =cd +er*cd -ei*sd;
          sd =sd +ei*dcd+er*sd;
          !** the following three statements compensate for truncation
          !** arithmetic. if rounded arithmetic is used, they may be removed
          del =0.5+0.5/(cd**2+sd**2);
          cd =del*cd;
          sd =del*sd;
          !** test (both are close)
!          cd =cos((j-1)*dtx)
!          sd =sin((j-1)*dtx)
          csum =csum +cd*a(j)-sd*b(j);
          ssum =ssum +cd*b(j)+sd*a(j);
        enddo;
        c(i) =csum;
        d(i) =ssum;
        dtx =dtx+dw;
      enddo;

      if (isn < 0) return;
      c(:n) =c(:n)*fn;
      d(:n) =d(:n)*fn;

      return;
      end subroutine dft;
