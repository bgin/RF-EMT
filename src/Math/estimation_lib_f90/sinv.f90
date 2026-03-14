! **********************************************************************
!
!       **************
!       *            *
!       *    sinv    *
!       *            *
!       **************
!
      subroutine sinv (ier,errm,  a,  n,eps,debugP);
!
!        sinf.f inverts a symmetric positive-definite matrix. It calculates inverse by
!        factoring matrix a using Cholesky decomposition. Thus, the input
!        matrix a must be positive definite. Note that error flag ier
!        should always be checked.
!
!**      B. Gibbs, 9/2010

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

      integer(4),intent(in) :: n;         !input order of matrix a.
      real(real_kind1),intent(in) :: eps; !input loss-of-precision tolerance for error return.
                                          ! eps = 1.d-4 is typical.
      logical(4),intent(in) :: debugP;    !if true, debug print may be output on unit 6
      real(real_kind1),intent(inout) :: a(n*(n+1)/2);!input/output matrix a(n,n), stored as a
                                          ! vector containing the upper triangular
                                          ! partition of the n x n symmetric matrix A, stored by
                                          ! columns.
      real(real_kind1),intent(out) :: errm;!minimum remaining precision (approximate)
      integer(4),intent(out) :: ier;      !output error flag.
                                          ! = 0              no error.
                                          ! = k with k > 0.  loss of precision ocurred.  the
                                          !                  radicand formed during factorization
                                          !                  of the k+1 row was still positive
                                          !                  but no longer greater than
                                          !                  abs(eps*a(k+1,k+1).
                                          ! = k with k < 0.  The radicand formed during factorization
                                          !                  of the k+1 row was negative,
                                          !                  indicating that matrix A is not
                                          !                  positive definite. The output A is meaningless.

      integer(4) j,k,ipiv,icoln,ihor,iver,irow,icol;
      real(real_kind1) dpiv,asum;
      !_______________________________________________


      !** compute Cholesky upper-triangular factor (U) of input A =U*U^T =L*L^T.
      !** The factor U is returned in A.
      !** Also compute minimum remaining precision (approximate)

      call cfactor2 (ier,errm,  a,  n,eps,debugP);
      if (ier < 0) return;

!***     invert upper triangular matrix U
      ipiv =n*(n+1)/2;
      a(ipiv) =1.d0/a(ipiv);
      icoln =ipiv-1;
      ipiv =ipiv -n;

!***     inversion loop
      do irow =n-1,1,-1;         !row loop
        dpiv =1.d0/a(ipiv);
        a(ipiv) =dpiv;

        j =icoln;
        do icol=1,n-irow;        !column-loop
          ihor =ipiv;
          iver =j;
          asum =0.d0;
          do k =irow,n-icol;     !inner loop
            iver =iver+1;
            ihor =ihor+k;
            asum =asum+a(iver)*a(ihor);
          enddo;

          a(j) =-asum*dpiv;
          j =j-(n-icol);
        enddo;        !___ end of icol loop

        ipiv =ipiv-irow;
        icoln =icoln-1;
      enddo;        !___ end of irow loop


!***     calculate inverse. inv(a)=inv(U)*transpose(inv(U))

      do icol=1,n;             !column loop
        ipiv =ipiv+icol;
        j =ipiv;
        do irow =icol,n;       !row loop
          ihor =j;
          asum =0.d0;

          do k =irow,n;        !start inner loop
            iver =ihor+irow-icol;
            asum =asum+a(ihor)*a(iver);
            ihor =ihor+k;
          enddo;

          a(j) =asum;
          j =j+irow;
        enddo;    !____ end row loop
      enddo;

      return;
      contains;
!*********************************************************************
      subroutine cfactor2 (ier,errm,  a,  n,eps,debugP);

!** subroutine cfactor computes the Cholesky factor of a symmetric
!** positive definite matrix A, i.e., A =U^T * U where U is upper triangular.
!** Matrix A is stored as upper triangular by columns and the output U is stored in A.
!** If matrix A is singular at row i, the same row of U will be set to zero and
!** ier will be set to the first row that is found to be singular.  cfactor also tests
!** for a loss of precision or singularity and returns an error code indicating
!** at which row the loss occurred.

!** Author: B. Gibbs, 9/2010

      use filterPrec, only: real_kind1;   !precision (4 or 8)
      implicit none;

      integer(4),intent(in) :: n;               !Dimension of A
      real(real_kind1),intent(in) :: eps;       !Tolerance for loss of precision, e.g. eps =1d-4
      logical(4),intent(in) :: debugP;          !if true, debug print may be output on unit 6
      real(real_kind1),intent(inout) :: a((n*(n+1))/2);  !Symmetric matrix to be factored, or output factor
      real(real_kind1),intent(out) :: errm;     !minimum remaining precision
      integer(4),intent(out) :: ier;            !Error flag: ier = -1 = n < 0, ier = -i means matrix
                                                ! is singular at row i, ier = +i, loss of precision
                                                ! exceeds eps at row i

      integer(4) i,k,kpiv,ki,iflg,km1;
      real(real_kind1) tol,dsum,dpiv,wmax;
      !____________________________________________________

      iflg =0;
!      if (ier == -2) iflg=1

      !**  test for invalid n
      if (n < 1) then;
        ier =-1;
        return;
      endif;

      ier =0;
      kpiv =0;          !initialize diagonal index
      wmax =0.d0;

      do k=1,n;         !row index
        kpiv =kpiv+k;   !index of diagonal element of row
        ki =kpiv;
        km1 =k-1;
        !**  calculate tolerance for testing loss of significance
        tol =abs(eps*a(kpiv));

        !**   start factorization-loop over k-th row
        do i=k,n;    !i is column index, ki is index of (k,i) element
          if (k > 1) then;
            dsum =dot_product(a(kpiv-km1:kpiv-1),a(ki-km1:ki-1));
          else;
            dsum =0.d0;
          endif;

          !**  compute difference for a(ki)
          dsum =a(ki)-dsum;

          if (i == k) then;
            !**  diagonal element: test for negative pivot element and for loss of significance
            wmax =max(wmax,abs(a(kpiv)/dsum));   !(a(kpiv)/dsum)**1.3 matches actual errors - WHY ??
            if (dsum <= tol) then;
              if (dsum <= 0.d0) then;
                if (debugP)                                             &
     &                write (6,'(/"MATRIX IS SINGULAR AT ROW ",I3)') i;
                if (ier >= 0) ier =-i;    !set for first row that is singular
                dpiv =0.d0;               !set dpiv to zero elements in row k
                a(kpiv) =0.d0;
                !dsum =1.d40  !when called by sinv, set diagonal to big number to get "pseudo-inverse" ?
                !return
              else;
                if (debugP) write (6,100) k,log10(abs(a(kpiv)/dsum));
  100           format(/'AT ROW',i5,',',f7.1,                           &
     &            ' DIGITS WERE LOST IN MATRIX FACTORIZATION')
                if (ier == 0) ier =k-1;
              endif;
            endif;

            if (dsum > 0.d0) then;
              !** compute pivot element
              a(kpiv) =sqrt(dsum);
              dpiv =1.d0/a(kpiv);
            endif;

          else;      !**  calculate terms in row
            a(ki) =dsum*dpiv;
          endif;

          ki =ki+i;
        enddo;       !____ end i loop
      enddo;      !____ end k loop

      if (real_kind1 == 8) then;
!        errm =1.1d-16*sqrt(wmax)      !little difference between using max vs RSS
        errm =1.1d-16*wmax;     !1.1d-16 is mantissa LSB (53) of IEEE T_floating on PC
      else if (real_kind1 == 4) then;
        errm =6.0e-8*wmax;     !6.0e-8 is mantissa LSB (24) of IEEE S_floating on PC
      else;
        errm =0.;
      endif;

      return;
      end subroutine cfactor2;

      end subroutine sinv;
