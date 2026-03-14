! ******************************************************************
      subroutine smop (c,  a,b,ip,iq,maxp,maxq);
!         computes the matrix function c =a*b*a**t where b anc c
!         symmetric stored upper triangular by columns

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
      use filterPrec, only: real_kind1;      !precision (4 or 8)
      implicit none;
      integer(4),intent(in) :: ip;           !row dimension of matrix a (and c)
      integer(4),intent(in) :: iq;           !column dimension of matrix a (and b)
      integer(4),intent(in) :: maxp;         !actual row dimension of a (maxp >= ip)
      integer(4),intent(in) :: maxq;         !not used
      real(real_kind1),intent(in) :: a(maxp,iq);      !input rectangular matrix
      real(real_kind1),intent(in) :: b(iq*(iq+1)/2);  !input symmetric matrix
      real(real_kind1),intent(out) :: c(ip*(ip+1)/2); !output symmetric matrix

      integer(4) i,j,k,nm,j1,k1,k2,l,li;
      real(real_kind1) asum;
      !______________________________________________
!
      nm =(ip*(ip+1))/2;
      c(:) =0.d0;

      do i=1,ip;
        li=(i*(i+1))/2;
        k1=0;
        do j=1,iq;
          asum =0.d0;
          do k=1,j;
            k1 =k1+1;
            asum =asum+a(i,k)*b(k1);
          enddo;

          if (j /= iq) then;
            j1=j+1;
            k2=k1+j;
            do k=j1,iq;
              asum =asum+a(i,k)*b(k2);
              k2=k2+k;
            enddo;
          endif;

          l=li;
          do k=i,ip;
            c(l) =c(l)+asum*a(k,j);
            l=l+k;
          enddo;
        enddo;
      enddo;

      return;
      end;
