! *******************************************************************
      subroutine thhc_m (r,hc,rsos,  n,mmax,m,nstrt,sparse);
!       thhcs is a general purpose routine for transforming an upper
!       triangular plus a rectangular 2-d array into an upper
!       triangular using orthogonal householder transformations.
!       thhcs is a modification of Bierman's (ESP) thhc which takes
!       advantage of sparseness to avoid calculations.

!     This "_m" version is a modified version of the original ESP routine writtten by Bierman and Nead.
!        Modifications in the thhc version by B. Gibbs 8/27/09

!         1) upgraded to Fortran 90 syntax,
!         2) implemented many vector operations using implicit MATLAB/Fortran 90 vector
!            operators (thus facilitating porting to MATLAB),
!         3) specifies input/output INTENT on all calling arguments,
!         4) rearranges calling arguments as "output, input/output, and input" to make
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
      use filterPrec, only: real_kind1;               !precision (4 or 8)
      implicit none;
      integer(4),intent(in) :: n;                     !number of rows in r
      integer(4),intent(in) :: mmax;                  !row dimension of hc in calling program
      integer(4),intent(in) :: m;                     !number of rows in hc  (m < mmax)
      integer(4),intent(in) :: nstrt;                 !index of first nonzero column in hc
      logical(4),intent(in) :: sparse;                !if true sparse logic is to be used
      real(real_kind1),intent(inout) :: r((n*(n+3))/2);!triangular matrix to be triangularized
                                                      !(if rsos >= 0, n+1 column contains information vector)
      real(real_kind1),intent(inout) :: hc(mmax,n+1); !measurement matrix with measurements in
                                                      ! last column (if rsos >= 0)
      real(real_kind1),intent(inout) :: rsos;         !residual sum-of-squares

!**     local
      integer(4) np1,ic,j,k,l,lr,nstart,k1;
      integer(4) jsparse,itot,jtot,isparse,nmap;
      real(real_kind1) a,sumh,sumr,vic,c,c1,fsparsei,fsparsej;

!**      epsilon is machine accuracy: 0.6d-37 is suitable for VAX single precision,
!**      1.0d-307 is limit for double precision on PC, but use 1d-100
      real(real_kind1) :: eps =1.0d-100;
      integer(4) :: mapt(mmax);
      !_______________________________________________

      if (m > mmax) stop 'THHCS ERROR: M > MMAX';
      isparse =0;
      jsparse =0;
      itot =0;
      jtot =0;
      np1 =n;
      if (rsos >= 0.d0) np1=n+1;

      nstart =min(max(nstrt,1),n);
      l =(nstart*(nstart-1))/2;

      !**   loop for columns
      do ic =nstart,n;
        !**  form householder vector and coefficient
        l=l+ic;
        a =r(l);
        sumh =sum(hc(1:m,ic)**2);
        nmap =0;
        if (sparse) then;
          !** row loop
          do j=1,m;
            if (hc(j,ic) /= 0.d0) then;
              nmap =nmap+1;
              mapt(nmap) =j;
            endif;
          enddo;
          itot =itot+m+1;
          isparse =isparse+nmap+1;
        endif;

        if (sumh > eps) then;
          sumh =sumh +a**2;
!            print *,'THHCS IC,SUM,NMAP',ic,sum,nmap
          sumr =sign(sqrt(sumh),a);
          vic =a+sumr;
          c =1.d0/(sumh+a*sumr);
          !**  vic | hc(*,ic) is the householder vector
          r(l) =-sumr;
          lr =l+ic;
          jtot =jtot+(np1-ic)*(nmap+1);

          !** loop thru columns
          if (sparse) then;
            !**  form v*col
            do j=ic+1,np1;
              sumh =vic*r(lr);
              do k1=1,nmap;
                k =mapt(k1);
                sumh =sumh +hc(k,ic)*hc(k,j);
              enddo;
              c1 =-c*sumh;
              if (c1 /= 0.d0) then;
                r(lr) =r(lr) +c1*vic;
                jsparse =jsparse+(nmap+1);
                do k1=1,nmap;
                  k =mapt(k1);
                  hc(k,j) =hc(k,j) +c1*hc(k,ic);
                enddo;
              endif;
              lr =lr+j;
            enddo;
          else;
            do j=ic+1,np1;
              sumh =vic*r(lr) +dot_product(hc(:m,ic),hc(:m,j));
              c1 =-c*sumh;
              if (c1 /= 0.d0) then;
                r(lr) =r(lr) +c1*vic;
                hc(:m,j) =hc(:m,j) +c1*hc(:m,ic);
              endif;
              lr =lr+j;
            enddo;
          endif;
        endif;
      enddo;

!**       calculate a posteriori residual root sum-of-squares
      if (rsos >= 0.d0) then;
        rsos =rsos**2 +sum(hc(:m,np1)**2);
        rsos =sqrt(rsos);
      endif;

!     if (sparse) then;
!       if (itot == 0 .or. jtot == 0)                                   &
!    &      print '(" thhcs itot,jtot=",2i5)',itot,jtot;
!       fsparsei =1.d0;
!       if (itot > 0) fsparsei=float(isparse)/float(itot);
!       fsparsej=1.d0;
!       if (jtot > 0) fsparsej=float(jsparse)/float(jtot);
!       if (fsparsei < 0.8 .or. fsparsej < 0.8)                         &
!    &  write (6,'(" fractional row/column sparseness =",2f10.3)')      &
!    &    fsparsei, fsparsej;
!     endif;

      return;
      end subroutine thhc_m;
