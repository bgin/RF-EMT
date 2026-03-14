! ****************************************************************
      subroutine triphu (u,v,  a,vd,   nda,n,d);

!         Form u such that         u * d1 * u**t = a * d * a**t.
!         Also a output such that       a * a**t = a * d * a**t  (input)

!         subroutine triphu triangularizes the input a matrix via
!         postmultiplication by a Householder transformation
!         i.e., a*d*a**t = (ap*t)*(ap*t)**t where ap=a*d**(1/2) and
!         ap*t is upper triangular.  The final step is conversion to u*d1*u**t where u is unit
!         upper triangular.  Note. a will usually be the result of phi*ui.
!         The transformation is defined as T =I-c*v*v**t and a**t*a = b**t*b
!         where vector a is a row of matrix A which is to be transformed to vector b

!         Note that triphu uses a separate module to specify whether arithmetic operations are single
!         or double precision.

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
      use filterPrec, only: real_kind1;   !precision 4 or 8
      implicit none;
      integer(4),intent(in) :: n;         !input size of input matrix a and output matrix u.
      integer(4),intent(in) :: nda;       !input row dimension of matrix a in the calling routine.
      real(8),intent(in) :: d(n);         !input vector of elements to represent diagonal matrix d.

      real(8),intent(inout) :: a(nda,n);  !input matrix a.  size (n,n).  row dimension nda.
                                          ! output matrix a such that on output a is defined as
                                          ! a(out) * a(out)**t = a(in) * d * a(in)**t.
      real(8),intent(inout) :: vd(n);     !unused

      real(8),intent(out) :: u(n*(n+1)/2);!output u-d factorization of the matrix product
                                          ! a * d * a**t.
      real(8),intent(out) :: v(n);        !required work vector

      integer(4) i,j,k,l,li,nm1,j0,ib,im1;
      real(8) c,vsum,diag,bi;

!      epsilon is machine accuracy: 0.6d-37 is suitable for VAX single precision,
!       1.0d-307 is limit for double precision on PC, but use 1d-100
      real(real_kind1) :: eps =1.0d-100;
      real(8) :: sqrEps =1.d-50;
      !________________________________________________

!***      multiply columns of a by sqrt of d
      do i=1,n;
        c =sqrt(d(i));
        a(:n,i) =c*a(:n,i);
      enddo;

!***      row loop of matrix a starting at bottom
      nm1=n-1;
      do ib=1,nm1;
        i =n+1-ib;
        im1 =i-1;
        j0 =0;
        !***  compute elementary householder vectors v,vd, and find
        !***  first non-zero element in row
        vsum =0.;
        do j=1,im1;
          v(j) =a(i,j);
          if (v(j) /= 0.d0 .and. j0 == 0) j0=j;
          a(i,j) =0.;
          vsum =vsum+v(j)**2;
        enddo;

        if (j0 == 0) cycle;   !row of a is all 0

        bi =sqrt(vsum+a(i,i)**2);
        if (a(i,i) > 0.d0) bi=-bi;
        v(i) =a(i,i)-bi;
        a(i,i) =bi;
        vsum =vsum+v(i)**2;
        if (vsum >= eps) then;
          c =2.d0/vsum;
          !***   now transform rows above row i
          do k=1,im1;
            vsum =c*dot_product(a(k,j0:i),v(j0:i));
            a(k,j0:i) =a(k,j0:i) -vsum*v(j0:i);
          enddo;
        endif;

      enddo;  !___ end ib (row) loop

!      call prtmx (a,nda,n,n,'a before')

!         now convert a and d to unit upper triangular u,d1
!         -------------------------------------------------
      l=1;
      u(1) =a(1,1)**2;
      do i=2,n;
        li =l+i;
        diag =a(i,i);
        u(li) =diag**2;
        im1 =i-1;
        if (u(li) < eps)  then;
          u(li) =eps;
          diag =sign(sqrEps,diag);
        endif;
        do j=1,im1;
          l =l+1;
          u(l) =a(j,i)/diag;
        enddo;
        l =li;
      enddo;

      return;
      end subroutine triphu;
