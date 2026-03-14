! ************************************************************
      subroutine wgs_m (u,v,  w,  imaxw,iw,jw,dw,sparse);
!
!     modified weighted gram-schmidt algorithm for reducing w*dw*w^T to u*d*u^T
!     form where u is a vector stored triangular matrix with the
!     resulting d elements stored on the diagonal.
!
!     (see book: ' factorization methods for discrete sequential estimation ',
!                      by g.j.bierman)
!
!     cognizant persons: g.j.bierman/m.w.nead   (jpl, feb.1978)
!     This "_m" version is a modified version of the original ESP routine
!        Modifications in the wgs_m version by B. Gibbs 8/27/09

!         1) upgraded to Fortran 90 syntax,
!         2) implemented many vector operations using implicit MATLAB/Fortran 90 vector
!            operators (thus facilitating porting to MATLAB),
!         3) specifies input/output INTENT on all calling arguments,
!         4) rearranges calling arguments as “output, input/output, and input” to make
!            the subroutines appear more like functions,
!         5) uses a separate module to specify whether arithmetic operations are single
!            or double precision, and
!         6) tests for matrix sparseness and only operate on non-zero elements

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
      use filterPrec, only: real_kind1;   !precision 4 or 8
      implicit  none;
      integer(4),intent(in) :: imaxw;     !row dimension of w matrix
      integer(4),intent(in) :: iw;        !no. rows of w matrix, dimension of u
      integer(4),intent(in) :: jw;        !no. cols of w matrix
      logical(4),intent(in) :: sparse;    !sparse operations are to be used
      real(real_kind1),intent(in) :: dw(jw);  !vector of non-negative weights for the
                                          ! orthogonalization process. the d's are unchanged
                                          ! by the calculation.
      real(real_kind1),intent(inout) :: w(imaxw,jw);!input matrix to be reduced to triangular form.
                                          ! this matrix is destroyed by the calculation
                                          ! iw <= imaxw .and. iw > 1
      real(real_kind1),intent(out) :: u(iw*(iw+1)/2);!upper triangular vector stored output
      real(real_kind1),intent(out) :: v(jw);       !work vector

      integer(4) i,j,k,l,kl,ki,ij,jm1,in(jw),numc;
      real(real_kind1) sum1,zero,dinv;
      !________________________________________________________

      zero =0.d0;
      if (iw > imaxw) call errhand ("invalid row index in wgs_m");

      do j=iw,2,-1;    !inverse row loop
        numc =0;
        v(:jw) =w(j,:jw);
        u(:jw) =dw(:jw)*v(:jw);
        sum1 =dot_product(v(:jw),u(:jw));
        do k=1,jw;   !column loop
          if (v(k) /= 0.d0) then;
            numc =numc+1;
            in(numc) =k;
          endif;
        enddo;

        w(j,j) =sum1;
        dinv =sum1;
        jm1 =j-1;

        if (sum1 <= zero) then;
          !**   w(j,.)=0. when dinv=0 (dinv=norm(w(j,.)**2))
          w(j,1:jm1) =zero;
!          write (6,'("zero row ",i4)') j
        else;
          do k=1,jm1;
            sum1 =zero;
            if (numc > 0) then;
              if (sparse) then;
                !** sum1 =dot_product(w(k,:),u(:))/dinv
                do ki=1,numc;
                  i =in(ki);
                  sum1 =sum1 +w(k,i)*u(i);
                enddo;
                sum1 =sum1/dinv;
                !**  divide here (in place of reciprocal multiply) to avoid
                !**  possible overflow
                !** w(k,:) =w(k,:) -sum1*v(:)
                do ki=1,numc;
                  i =in(ki);
                  w(k,i) =w(k,i)-sum1*v(i);
                enddo;
              else;
                sum1 =dot_product(w(k,:),u(:jw))/dinv;
                w(k,:) =w(k,:) -sum1*v(:);
              endif;
            endif;
            w(j,k) =sum1;
          enddo;
        endif;
      enddo;   !___ end l loop

!***      the left lower part of w is u transpose
      u(1) =dot_product(dw(1:jw),w(1,1:jw)**2);
      ij =1;
      do j=2,iw;
        do i=1,j;
          ij =ij+1;
          u(ij) =w(j,i);
        enddo;
      enddo;

      return;
      end subroutine wgs_m;
