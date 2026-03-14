! *************************************************************************
      subroutine ud2cov_m (pout,  uin,n);
!
!     computes a covariance from its u-d factorization. both matrices
!     are vector stored and the output covariance can overwrite the
!     input u-d array. uin=u-d is related to pout via pout=udu(**t)
!
!     cognizant persons: g.j.bierman/m.w.nead  (jpl, feb. 1977)
!
!     This "_m" version is a modified version of the original ESP routine
!        Modifications in the rank1_m version by B. Gibbs 8/27/09

!         1) upgraded to Fortran 90 syntax,
!         2) specifies input/output INTENT on all calling arguments,
!         3) rearranges calling arguments as “output, input/output, and input” to make
!            the subroutines appear more like functions,
!         4) uses a separate module to specify whether arithmetic operations are single
!            or double precision, and

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
      use filterPrec, only: real_kind1;              !precision (4 or 8)
      implicit  none;
      integer(4),intent(in) :: n;                    !dimension of the matrices involved, n > 1
      real(real_kind1),intent(in) :: uin(n*(n+1)/2); !u-d factors, vector stored with the d
                                                     ! entries stored on the diagonal of uin
      real(real_kind1),intent(out) :: pout(n*(n+1)/2);!output covariance, vector stored.
                                                      ! (pout=uin is permitted)

      integer(4) i,j,k,ii,ik,jj,jjl,jm1;
      real(real_kind1) s,alpha;
      !________________________________________________________

      pout(1) =uin(1);
      jj =1;
      do j=2,n;
        jjl =jj;
        jj =jj+j;
        pout(jj) =uin(jj);
        s =pout(jj);
        ii =0;
        jm1 =j-1;
        do i=1,jm1;
          ii =ii+i;
          alpha =s*uin(jjl+i);
          ik =ii;
          do k=i,jm1;
            pout(ik) =pout(ik) +alpha*uin(jjl+k); ! jjl+k =(k,j)
            ik =ik+k;
          enddo;
          pout(jjl+i) =alpha;
        enddo;
      enddo;

      return;
      end subroutine ud2cov_m;
