!*******************************************************
      real(8) function l2normM (c,nd,n);
      !** compute the Frobenius L2 norm of matrix C

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
      integer(4),intent(in) :: nd;    !row dimension of matrix c
      integer(4),intent(in) :: n;     !actual dimension of c (n <= nd)
      real(8),intent(in) :: c(nd,n);  !matrix


      integer(4) i;
      !____________________________________

      l2normM =0.d0;
      do i=1,n;
        l2normM =l2normM +sum(c(i,:n)**2);
      enddo;
      l2normM =sqrt(l2normM); !/(n**2))

      return;
      end function l2normM;
