! ***********************************************************
      subroutine phiu_m (pu,  phi,nph,nr,nc,u,nu,npu);
!         subroutine phiu computes pu =phi*u where phi is a
!         sparse, 2-dimensional matrix and u is unit upper
!         triangular

!     Although this routine has a similar function as Bierman's ESP routine phiu,
!     it was developed independently.

!        Features in phiu_m by B. Gibbs 8/27/09
!         1) Fortran 90 syntax,
!         2) implements vector operation using implicit MATLAB/Fortran 90 vector
!            operators (thus facilitating porting to MATLAB),
!         3) specifies input/output INTENT on all calling arguments,
!         4) arranges calling arguments as “output, input/output, and input” to make
!            the subroutines appear more like functions,
!         5) uses a separate module to specify whether arithmetic operations are single
!            or double precision, and
!         6) tests for non-zero start/stop elements in each phi row and only operate on
!            the non-zero elements

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
      implicit none;
      integer(4),intent(in) :: nph;  !row dimension of phi
      integer(4),intent(in) :: nr;   !# rows of phi
      integer(4),intent(in) :: nc;   !# columns of phi
      integer(4),intent(in) :: nu;   !dimension of u - not used
      integer(4),intent(in) :: npu;  !row dimension of pu
      real(real_kind1),intent(in) :: phi(nph,nc);
      real(real_kind1),intent(in) :: u(nc*(nc+1)/2);
      real(real_kind1),intent(out) :: pu(npu,nc);

      integer(4) i,j,k,l,in1,in2,jm1;
      real(real_kind1) t(nc),sum1;
      !__________________________________________________________-

!***      row loop of phi
      do i=1,nr;
        !** store row of a in t. Also find first and last non-zero element in row
        in1 =0;
        in2 =0;
        do j=1,nc;
          t(j) =phi(i,j);
          if (t(j) /= 0.d0) then;
            if (in1 == 0) in1 =j;
            in2 =j;
          endif;
        enddo;

        !**  now multiply t*u
        l =1;
        pu(i,1) =t(1);
        do j=2,nc;
          sum1 =t(j);
          jm1 =min(j-1,in2);
          if (in1 > 0 .and. jm1 >= in1) then;
            sum1 =sum1 +dot_product(t(in1:jm1),u(l+in1:l+jm1));
          endif;
          l =l+j;
          pu(i,j) =sum1;
        enddo;
      enddo;  !___ end i loop

      return;
      end subroutine phiu_m;
