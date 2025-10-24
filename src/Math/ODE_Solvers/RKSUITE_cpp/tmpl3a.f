C
C                             Template 3a
C
C  This template comes in two variants, 3a and 3b, in files tmpl3a.for and
C  tmpl3b.for, respectively.  They illustrate different aspects of the RK
C  suite on the same problem.  The task could be handled more simply with
C  UT; the purpose of the template is to illustrate how to use CT when more
C  control over the integration is required.  Specifically, variant 3a
C  illustrates how to use CT with METHOD = 2 to advance the integration a
C  step at a time and how to use INTRP with this METHOD to get answers at
C  specific points by interpolation.
C
C  Problem:  Integrate a forced Duffing equation that models an electric
C            circuit with a nonlinear inductor.  The behavior of solutions
C            depends on the value of a parameter B that governs the size of
C            the forcing term.  This problem is discussed at length in F.C.
C            Moon, Chaotic Vibrations, John Wiley & Sons, Inc., New York,
C            1987;  numerical results are reported on p. 272.  Written as
C            a first order system, the second order equation is
C
C               y1' = y2,
C               y2' = - 0.1 y2 - y1**3 + B cos(t).
C
C            A Poincare map is used to study the qualitative behavior of the
C            solutions.  This means that solution values are to be computed
C            at t = 2*pi*k for k = 1,2,...,n and y2(t) plotted against y1(t).
C
C            Exploration of the qualitative behavior of solutions could lead
C            to integration of the differential equation for a great many
C            initial values.  Here the single set of initial values
C                         y1(0) = 3.3,        y2(0) = -3.3
C            is considered. It is planned that the integration extend for
C            many multiples of the period 2*pi of the forcing function.
C            The template takes n = 100.  To have some accuracy at the end
C            of a "long" integration, it is necessary to specify a moderately
C            stringent tolerance.  In this variant TOL is taken to be 1.0D-5.
C            The general scaling of the initial values and the anticipated
C            behavior of the solution components suggest that threshold
C            values of 1.0D-10 on each component might be reasonable.
C
C            At this tolerance METHOD = 2 is a reasonable choice.  Output
C            is desired at specific points.  By using INTRP, output this
C            infrequently has no impact on the cost of the integration.
C
C            The character of the problem depends on the value of B. According
C            to Moon, the solution is periodic when B = 9.8 and becomes
C            chaotic when B is about 10.0.  The behavior of the global error
C            is very different in these circumstances.  By definition the
C            solution is very sensitive to changes in the initial conditions
C            when there is chaos.  This means that it is very difficult to
C            compute accurately a specific solution curve -- the global error
C            must get large rather quickly. Variant 3a takes B = 9.8 so that
C            the solution is periodic.  Also included in tmpl3.out is the
C            output when B = 11.0, a chaotic solution.
C
C            With both variants of the template the integrations are
C            sufficiently expensive that if MESAGE is set .TRUE., there are
C            a considerable number of returns reporting that 5000 function
C            evaluations have been made since the last message.  To suppress
C            these messages, MESAGE is set .FALSE.  Another purpose of this
C            template is to illustrate the processing of error messages and
C            the decisions that are made at each step of an integration.
C
C  NOTES:    Working storage must be provided in the array WORK(*) of
C            length LENWRK.  Although storage is no problem with only
C            NEQ = 2 equations, LENWRK is taken here to be 14*NEQ,
C            the least that can be used in CT with METHOD = 2 when global
C            error assessment is not done (ERRASS = .FALSE.).  Storage must
C            also be supplied for interpolation in the array WRKINT(*) of
C            length NINT.  The  amount depends on the number NWANT of
C            solution components desired.  In this template, all the
C            components are approximated so NWANT = NEQ.  NINT is taken
C            here to be the minimum needed with this METHOD, namely
C            NEQ+MAX(NEQ,5*NWANT) = 6*NEQ.
C
C     .. Parameters ..
      INTEGER           NEQ, METHOD, NWANT, NINT, 
     &                  LENWRK, OUTFIL, MAXPER
      PARAMETER         (NEQ=2,METHOD=2,NWANT=NEQ,NINT=6*NEQ,
     &                  LENWRK=14*NEQ,OUTFIL=9,MAXPER=100)
      DOUBLE PRECISION  ZERO, ONE, FOUR
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0,FOUR=4.0D0)
C     .. Local Scalars ..
      DOUBLE PRECISION  HNEXT, HSTART, PI, TEND, TNOW,
     &                  TOL, TOUT, TSTART, TWOPI, WASTE
      INTEGER           CFLAG, KOUNTR, L, STPCST, STPSOK, TOTF
      LOGICAL           ERRASS, MESAGE
C     .. Local Arrays ..
      DOUBLE PRECISION  THRES(NEQ), WORK(LENWRK), WRKINT(NINT),
     &                  YNOW(NEQ), YOUT(NEQ), YPNOW(NEQ), YPOUT(1),
     &                  YSTART(NEQ)
C     .. External Subroutines ..
      EXTERNAL          CT, F, INTRP, SETUP, STAT
C     .. Scalars in Common ..
      DOUBLE PRECISION  B
C     .. Intrinsic Functions ..
      INTRINSIC         ATAN
C     .. Common blocks ..
      COMMON            /BPARAM/B
C     .. Executable Statements ..
C
      PI = FOUR*ATAN(ONE)
      TWOPI = PI + PI
C
C  Set the parameter defining the character of the solutions.
C
      B = 9.8D0
C
C  Initialize output.
C
      WRITE (*,'(A,I3,A,F5.1/)')
     &' Template 3a with METHOD = ',METHOD,' and B = ', B
C
C  Set the initial conditions.  TEND is well past the last point
C  of interest, MAXPER*2*PI, so that the code can step past this
C  point and obtain the answer there by interpolation.
C
      TSTART = ZERO
      YSTART(1) =  3.3D0
      YSTART(2) = -3.3D0
      TEND = TSTART + MAXPER*TWOPI + TWOPI
C
C  Because the solution components may be computed at many points, they are
C  written to a file on unit OUTFIL.  Subsequently, the results can be studied,
C  plotted, or manipulated as necessary.  Naming of files is system-dependent,
C  so the name here, PLOT.DAT, may need to be changed.  Begin by writing the
C  initial values to PLOT.DAT in a form suitable for plotting y2 against y1.
C
      OPEN (UNIT=OUTFIL,FILE='PLOT.DAT')
C
      WRITE (OUTFIL,'(2E13.5)') YSTART(1), YSTART(2)
C
C  Set error control parameters.
C
      TOL = 1.0D-5
      DO 20 L = 1, NEQ
         THRES(L) = 1.0D-10
   20 CONTINUE
C
C  Call the setup routine. Set MESAGE = .FALSE. because this will be a
C  a relatively expensive integration and it is annoying to be told
C  repeatedly that 5000 function evaluations have been required.  The
C  global error is not assessed. Automatic selection of an initial step
C  size is specified by setting HSTART to zero.
C
      MESAGE = .FALSE.
      ERRASS = .FALSE.
      HSTART = ZERO
      CALL SETUP(NEQ,TSTART,YSTART,TEND,TOL,THRES,METHOD,'Complex Task',
     &           ERRASS,HSTART,WORK,LENWRK,MESAGE)
C
C  KOUNTR counts the number of output points.
C
      KOUNTR = 1
      TOUT = TSTART + TWOPI
   40 CONTINUE
      CALL CT(F,TNOW,YNOW,YPNOW,WORK,CFLAG)
C
C  Loop until an output point is reached, or the code gets into
C  trouble.  Interrogate CFLAG to decide on an appropriate course
C  of action.  Because MESAGE = .FALSE., it is necessary to write 
C  out an explanation of any trouble except for a catastrophe when an
C  explanation is provided automatically.
C
C  Start of "CASE" statement implemented with COMPUTED GO TO.
      GO TO (60,60,100,120,140,160) CFLAG
C
   60 CONTINUE
C
C  CFLAG = 1 or 2.  Successful step.  Have we reached an output point?
   80 CONTINUE
      IF (TNOW.LT.TOUT) THEN
C
C  Take another step.
         GO TO 40
      ELSE
C
C  Reached an output point.  Call INTRP to obtain the answer at TOUT.
C  Specify 'Solution only' and all the components of the solution
C  (NWANT = NEQ).  The values are output in YOUT(*).  It is necessary
C  to provide storage only for the number of components requested and
C  here YPOUT is a dummy.  In general more than one output point might
C  be obtained by interpolation after a single step by CT.  That is the
C  reason for the jump to the statement labelled 80.  For many problems 
C  it is very important to take advantage of this possibility.
C
         CALL INTRP(TOUT,'Solution only',NWANT,YOUT,YPOUT,F,WORK,
     &              WRKINT,NINT)
         WRITE (OUTFIL,'(2E13.5)') YOUT(1), YOUT(2)
         IF (KOUNTR.LT.MAXPER) THEN
            KOUNTR = KOUNTR + 1
            TOUT = TOUT + TWOPI
            GO TO 80
         END IF
         GO TO 180
      END IF
C
  100 CONTINUE
C
C  CFLAG = 3.  Although a considerable amount of work has been done, continue.
      GO TO 40
C
  120 CONTINUE
C
C  CFLAG = 4.  The problem appears to be stiff.
      WRITE (*,'(A/)') ' The problem appears to be stiff. '
      STOP
C
  140 CONTINUE
C
C  CFLAG = 5.  The accuracy request is too stringent.
      WRITE (*,'(A/)') ' The accuracy request is too stringent. '
      STOP
C
  160 CONTINUE
C
C  CFLAG = 6.  The global error estimation scheme broke down.  This 
C  cannot occur because we did not ask for global error estimation.
C
      WRITE (*,'(A/)') ' Global error estimation broke down. '
      STOP
C
  180 CONTINUE
C  End of "CASE" statement implemented with COMPUTED GO TO.
C
      CLOSE (OUTFIL)
C
C  The subroutine STAT is used to obtain some information about the progress
C  of the integration. TOTF is the total number of calls to F made so far
C  in the integration; it is a machine-independent measure of work.  At present
C  the integration is finished, so the value printed out refers to the overall
C  cost of the integration.  STPCST is the cost in calls to F of a typical
C  step with this method.  WASTE is the fraction of steps after the first that
C  were rejected.  A "large" value of WASTE is usually the result of a mistake
C  in the coding of F.  It can also be the result of a solution component that
C  is not smooth -- singular or with many places where a low order derivative
C  is not continuous.  STPSOK is the number of accepted steps.
C
      CALL STAT(TOTF,STPCST,WASTE,STPSOK,HNEXT)
      WRITE (*,'(/A,1PE10.2,0P/A,I10/A,I10/A,F10.2/A,I10)')
     &' The integration reached                      ', TNOW,
     &' The cost of the integration in calls to F was', TOTF,
     &' The number of calls to F per step is         ', STPCST,
     &' The fraction of failed steps was             ', WASTE,
     &' The number of accepted steps was             ', STPSOK
C
      STOP
      END
      SUBROUTINE F(T,Y,YP)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
C     .. Array Arguments ..
      DOUBLE PRECISION  Y(*), YP(*)
C     .. Parameters ..
      DOUBLE PRECISION  K
      PARAMETER         (K=0.1)
C     .. Scalars in Common ..
      DOUBLE PRECISION  B
C     .. Intrinsic Functions ..
      INTRINSIC         COS
C     .. Common blocks ..
      COMMON            /BPARAM/B
C     .. Executable Statements ..
      YP(1) = Y(2)
      YP(2) = -K*Y(2) - Y(1)**3 + B*COS(T)
      RETURN
      END
