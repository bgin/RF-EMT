      module optlDat;

!** module optllDat includes model info for the funct subroutine called by
!** the optimization routines.  Other variables may be added to optlDat
!** as needed for specific models.

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
      type fstats;
        integer(4) iter,nmeas;
        real(8) rrms,cost,testc;
      end type;
      type (fstats) fstat;

      logical srif;
      end module optlDat;
