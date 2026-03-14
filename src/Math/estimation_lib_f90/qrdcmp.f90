! *******************************************************************
      subroutine QRdcmp (hc,rsos,   mmax,n,m);
!       QRdcmp is a general purpose routine for transforming
!       a rectangular 2-d array into an upper triangular array
!       using orthogonal householder transformations.

!**     QRdcmp is about an order of magnitude less accurate than THHCS_m,
!**     probably because it does not modify elements of hc above the diagonal
!**     at each step    !!!!!!

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
      integer(4),intent(in) :: n;            !number of rows in r
      integer(4),intent(in) :: mmax;         !row dimension of Hc
      integer(4),intent(in) :: m;            !number of rows in hc  (m <= mmax)
      real(8),intent(inout) :: hc(mmax,n+1); !measurement matrix with measurements in
                                             ! last column if rsos >= 0.  On output
                                             ! hc is the triangularized matrix + RHS
      real(8),intent(inout) :: rsos;         !output residual sum-of-squares if rsos >= 0 on input.

!       local
      integer(4) np1,ic,j;
      real(8) a,sumh,sumr,vic,c,c1,scaleh;
      real(8) :: eps =1.d-60;
      !_______________________________________________

      if (m > mmax) stop 'QRdcmp ERROR: M > MMAX';
      np1 =n;
      if (rsos >= 0.d0) np1 =n+1;

!**      loop for columns
      do ic =1,n-1;
        !** form householder vector and coefficient
        scaleh =maxval(abs(hc(ic:m,ic)));

        if (scaleh > eps) then;
          hc(ic:m,ic) =hc(ic:m,ic)/scaleh;   !doesn't seem to improve result much
          a =hc(ic,ic);
          sumh =sum(hc(ic:m,ic)**2);

          if (sumh > eps) then;
            !print *,'QRdcmp IC,SUM,NMAP',ic,sum,nmap
            sumr =sign(sqrt(sumh),a);
            vic =a +sumr;
            c =1.d0/(sumh +a*sumr);
!            c =1.d0/(sumr*(a+sumr))
            hc(ic,ic) =vic;   !vic | hc(ic+1:m,ic) is the householder vector

            !** loop thru columns
            do j=ic+1,np1;
              !**  form v*col
              sumh =dot_product(hc(ic:m,ic),hc(ic:m,j));
              c1 =c*sumh;
              if (c1 /= 0.d0) then;
                hc(ic:m,j) =hc(ic:m,j) -c1*hc(ic:m,ic);
              endif;
            enddo;
            hc(ic,ic) =-sumr*scaleh;    !reset diagonal of column ic
          endif;
        else;
          hc(ic,ic) =0.d0;
        endif;
      enddo;

!**   calculate a posteriori residual root sum-of-squares
      if (rsos >= 0.d0) then;
        rsos =sqrt(sum(hc(:m,np1)**2));
      endif;

      return;
      end subroutine QRdcmp;
