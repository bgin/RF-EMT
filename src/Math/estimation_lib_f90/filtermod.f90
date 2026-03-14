      module filterMod;
!** filterMod contains variables used in various Kalman filter and smoother software.
!** It is only necessary to set variables that are actually used in the filter code.
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
      use filterPrec, only: real_kind1;   !precision (4 or 8)

      !** contains global variables used in defining Kalman filter, UD and SRIF
      implicit none;
      integer(4) nstate;        !total number of states in model (set by model)
      integer(4) ncore;         !number of core states in model
      integer(4) nMarkov;       !number of markov states in model
      integer(4) nbias;         !number of bias states in model

      integer(4) nest;          !total number of estimated states
      integer(4) nConsider;     !total number of consider parameters
      integer(4) nUnadj;        !total number of unadjusted-analyze parameters
      integer(4) ncoreE;        !number of estimated core states
      integer(4) nMarkE;        !number of estimated Markov states
      integer(4) nbiasE;        !number of estimated bias states
      integer(4) nDynE;         !number of estimated core+markov states

      integer(4) nmeas;         !total number of vector meas in input file
      integer(4) :: method_udtime =1;!method used for ud time update: 1=wgs, 2=rank-1, 3=square phi*u*d*u^T*phi^T  ###
      integer(4) :: KFmethod =2;   !Kalman filter (covariance) measurement update method: 1=short, 2=Bierman;s Joseph, 3=full Joseph ###

      integer(4),parameter :: mmaxi =10;  !maximum number of allowed scalar measurements at one time point (10)

      integer(4),allocatable :: mapxe(:);!mapping from estimated states (nest) to total model states
      integer(4),allocatable :: mapxu(:);!mapping from unadjusted analyze parameter (nUnadj) to total model states
      integer(4),allocatable :: mapxei(:);!mapping from total model states (nstate) to estimated states
      integer(4),allocatable :: mapxui(:);!mapping from total model states (nstate) to estimated states

!      character(80) setupFile   !name of file containing data for structure xstate
!      character(80) measFile    !name of file containing measurement data
      character(20),allocatable :: xnames(:); !names of total states as set by model routines

      logical(4) :: sparse =.true.;  !true = use sparse multiplies  !###
      logical(4) :: debugp =.true.;  !generates debug print
      logical(4) :: first(3) =.true.;!initializes filters
      logical(4) :: runFilt(3) =.true.; !denotes whether Kalman covariance, U-D and SRIF filters are to be run

      real(8) tepoch;           !epoch time of x0 and p0
      real(8) tend;             !end time for filter processing
      real(8) rsum(3);          !residual sum-of-squares for 3 filters

      real(real_kind1),allocatable :: xfu(:),ud(:),xfk(:),pf(:,:),xfs(:);
      real(real_kind1),allocatable :: xfa(:),pfa(:,:),rsrif(:);
      real(8),allocatable :: dxsum(:,:),sigsum(:,:);

!      real(8) rnoise(mmaxi,mmaxi)           !meas noise covariance

      type states;
        integer(4) status;      !0=ignore, 1=unadjusted analyze, 2=consider, 3=estimate
        integer(4) xtype;       !1=core state, 2=1 or 2 order Markov, 3=bias
        character(20) xname;    !name of state that must match name in hard-coded model
        real(8) x0;             !a priori estimate at epoch time
        real(8) sig0;           !1-sigma uncertainty of x0 at epoch time
        real(8) qs;             !PSD of driving white noise
        real(8) const(2);       !known constants associated with this state in dynamic or measurement model
      end type;
      type (states),allocatable :: xstate(:);

      type meass;
        integer(4) mtype;       !meas type (must match types in model)
        integer(4) m;           !number of scalar meas included in this meas
        character(20) mname;    !name of meas for output (may include sensor name)
        real(8) t;              !meas time
        real(8) y(mmaxi);        !up to mmax scalar meas values
        real(8) sigr(mmaxi);     !1-sigma meas noise for mmax scalar meas
      end type;
!      type (meass),allocatable :: meas(:)

      end module;
