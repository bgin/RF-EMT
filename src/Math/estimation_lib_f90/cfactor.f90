      subroutine cfactor (ier,  a,  n,eps);

!** subroutine cfactor computes the Cholesky factor of a symmetric
!** positive definite matrix A, i.e., A =U^T * U where U is upper triangular.
!** Matrix A is stored as upper triangular by columns and the output U is stored in A.
!** If matrix A is singular at row i, the same row of U will be set to zero and
!** ier will be set to the first row that is found to be singular.  cfactor also tests
!** for a loss of precision or singularity and returns an error code indicating
!** at which row the loss occurred.

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

!      use filterPrec, only: real_kind1    !precision (4 or 8)
      implicit none;

      integer(4),intent(in) :: n;               !Dimension of A
      real(8),intent(in) :: eps;                !Tolerance for loss of precision, e.g. eps =1d-4
      real(8),intent(inout) :: a((n*(n+1))/2);  !Symmetric matrix to be factored, or output factor
!      real(8),intent(out) :: errm               !minimum remaining precision
      integer(4),intent(out) :: ier;            !Error flag: ier = -1 = n < 0, ier = -i means matrix
                                                ! is singular at row i, ier = +i, loss of precision
                                                ! exceeds eps at row i

      integer(4) i,k,kpiv,ki,iflg,km1;
      real(8) tol,dsum,dpiv,work,wmax;
      real(8) errm;   !minimum remaining precision
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
                write (6,'(/"MATRIX IS SINGULAR AT ROW ",I3)') i;
                if (ier >= 0) ier =-i;    !set for first row that is singular
                dpiv =0.d0;               !set dpiv to zero elements in row k
                a(kpiv) =0.d0;
                !dsum =1.d40  !when called by sinv, set diagonal to big number to get "pseudo-inverse" ?
                !return
              else;
                work =log10(abs(a(kpiv)/dsum));
                write (6,100) k,work;
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

!      errm =1.1d-16*sqrt(wmax)      !little difference between using max vs RSS
      errm =1.1d-16*wmax;     !1.1d-16 is mantissa LSB (53) of IEEE S_floating on PC

      return;
      end subroutine cfactor;
