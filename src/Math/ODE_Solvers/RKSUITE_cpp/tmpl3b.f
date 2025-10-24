C
C                             Template 3b
C
C  This template comes in two variants, 3a and 3b, in files tmpl3a.for and
C  tmpl3b.for, respectively.  They illustrate different aspects of the RK
C  suite on the same problem.  The task could be handled more simply with
C  UT; the purpose of the template is to illustrate how to use CT when more
C  control over the integration is required.  Specifically, variant 3b
C  illustrates how to use CT with METHOD = 3 to advance the integration a
C  step at a time and how to get answers at specific points.  It also
C  illustrates how to assess the global error with GLBERR.
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
C            stringent tolerance.  In this variant TOL is taken to be 1.0D-7.
C            The general scaling of the initial values and the anticipated
C            behavior of the solution components suggest that threshold
C            values of 1.0D-10 on each component might be reasonable.
C
C            The tolerance has been reduced substantially from the value
C            used in Variant 3a.  At this tolerance METHOD = 3 is preferred.
C            Although output is desired at specific points, they do not
C            occur frequently enough to impact the integration -- the
C            tolerance is stringent enough that there are a number of steps
C            between output points.  One purpose of this template is to show
C            how to obtain output at specific points with METHOD = 3.
C
C            Another purpose of this template is to illustrate global error
C            assessment. GLBERR can be called after any call to CT to assess
C            the global error.  If the global error gets too large to be
C            meaningful or if the code has evidence that the scheme for
C            assessing the global error might be unreliable, the integration
C            is terminated by CT.  In the template the global error is reported
C            at the end of the run or whenever the assessment breaks down.
C
C            The character of the problem depends on the value of B. According
C            to Moon, the solution is periodic when B = 9.8 and becomes
C            chaotic when B is about 10.0.  The behavior of the global error
C            is very different in these circumstances.  By definition the
C            solution is very sensitive to changes in the initial conditions
C            when there is chaos.  This means that it is very difficult to
C            compute accurately a specific solution curve -- the global error
C            must get large rather quickly. Variant 3b takes B = 11.0 so that 
C            chaos is present, and the global error assessment will eventually
C            break down.  Also included in tmpl3.out is the output when 
C            B = 9.8, a periodic solution.
C
C            With the codes of RKSUITE, precisely the same results are 
C            obtained whether or not global error assessment is requested.
C            When METHOD = 2 or 3, a run costs approximately three times as 
C            many calls to the routine defining the differential equations 
C            when ERRASS is .TRUE. as when it is .FALSE., and when METHOD = 1,
C            it costs approximately four times as many.  For this reason the
C            facility is normally used only for spot checks.  A run with
C            B = 9.8 is relatively expensive; you might want to reduce the
C            the run time by reducing MAXPER.
C
C            With both variants of the template the integrations are
C            sufficiently expensive that if MESAGE is set .TRUE., there are
C            a number of returns reporting that 5000 function evaluations
C            have been made since the last message.  To suppress these
C            messages, MESAGE is set .FALSE.  Another purpose of this
C            template is to illustrate the processing of error messages and
C            the decisions that are made at each step of an integration.
C
C  NOTES:    As may be seen in tmpl3.out, with a specific computer, precision,
C            compiler, ... , the global error assessment breaks down at about
C            t = 121 after 34348 function evaluations have been made.  The
C            maximum weighted error seen at any step prior to this point is
C            about 0.07!  A more realistic assessment of the accuracy is
C            provided by the individual RMS weighted errors.  At this point
C            the larger is about 0.003.  Obviously it is difficult to follow
C            a solution in the presence of chaos even with a stringent
C            tolerance.  If the global error assessment is turned off, the
C            integration will proceed to completion at t = 200*pi.  It costs
C            65414 function evaluations.  It is important to understand what
C            all this means.  The integration is successful because the code
C            controls the local error at each step -- what it is supposed to
C            do.  However, when the problem itself is not stable, control of
C            local error does not imply control of the global error that
C            really interests us.  The global error assessment brings to our
C            attention that we are not getting what we want.  The situation
C            is different when B = 9.8 and the solution is periodic.  As the
C            results for this case in tmpl3.out show, this solution has an
C            acceptable accuracy for the whole interval of integration.
C
C            For the typical problem that is moderately stable the global
C            error is comparable to the local tolerance, and users come to
C            expect this.  It is unfortunate that the global error assessment
C            facility is too expensive for routine use because there is no 
C            doubt that many solutions of initial value problem are not nearly
C            as accurate as the users of the software think they are.
C
C            Working storage must be provided in the array WORK(*) of
C            length LENWRK.  Although storage is no problem with only
C            NEQ = 2 equations, LENWRK is taken here to be 21*NEQ,
C            the least that can be used in CT with METHOD = 3 when global
C            error assessment is done (ERRASS = .TRUE.).
C
C     .. Parameters ..
      INTEGER           NEQ, METHOD, LENWRK, OUTFIL, MAXPER
      PARAMETER         (NEQ=2,METHOD=3,LENWRK=21*NEQ,OUTFIL=9,
     &                  MAXPER=100)
      DOUBLE PRECISION  ZERO, ONE, FOUR
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0,FOUR=4.0D0)
C     .. Local Scalars ..
      DOUBLE PRECISION  ERRMAX, HNEXT, HSTART, PI, TEND, TERRMX,
     &                  TNOW, TOL, TSTART, TWOPI, WASTE
      INTEGER           CFLAG, KOUNTR, L, STPCST, STPSOK, TOTF
      LOGICAL           ERRASS, MESAGE
C     .. Local Arrays ..
      DOUBLE PRECISION  RMSERR(NEQ), THRES(NEQ), WORK(LENWRK),
     &                  YNOW(NEQ), YPNOW(NEQ), YSTART(NEQ)
C     .. External Subroutines ..
      EXTERNAL          CT, ENVIRN, F, GLBERR, RESET, SETUP, STAT
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
      B = 11.0D0
C
C  Initialize output.
C
      WRITE (*,'(A,I3,A,F5.1/)')
     &' Template 3b with METHOD = ',METHOD,' and B = ', B
C
C  Set the initial conditions.  TEND is set to the first output point
C  because with METHOD = 3, output is obtained by integrating to TEND
C  and then resetting TEND to obtain the next output point.
C
      TSTART = ZERO
      YSTART(1) =  3.3D0
      YSTART(2) = -3.3D0
      TEND = TSTART + TWOPI
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
      TOL = 1.0D-7
      DO 20 L = 1, NEQ
         THRES(L) = 1.0D-10
   20 CONTINUE
C
C  Call the setup routine. Set MESAGE = .FALSE. because this will be a
C  a relatively expensive integration and it is annoying to be told
C  repeatedly that 5000 function evaluations have been required.  The
C  global error is assessed. Automatic selection of an initial step
C  size is specified by setting HSTART to zero.
C
      MESAGE = .FALSE.
      ERRASS = .TRUE.
      HSTART = ZERO
      CALL SETUP(NEQ,TSTART,YSTART,TEND,TOL,THRES,METHOD,'Complex Task',
     &           ERRASS,HSTART,WORK,LENWRK,MESAGE)
C
C  KOUNTR counts the number of output points.
C
      KOUNTR = 1
   40 CONTINUE
      CALL CT(F,TNOW,YNOW,YPNOW,WORK,CFLAG)
C
C  The code is not permitted to go past TEND, the end of the interval of
C  integration, so it will shorten the step when necessary so as to produce
C  an answer at TEND exactly.  To obtain answers at specific points with
C  METHOD = 3, which has no interpolation capability, this fact is exploited
C  by stepping to TEND, and then resetting TEND to the next outpoint point
C  with the subroutine RESET.  KOUNTR counts the number of output points.
C
C  Loop until an output point is reached, or the code gets into trouble.
C  Interrogate CFLAG to decide on an appropriate course of action.
C  Because MESAGE = .FALSE., it is necessary to write out an explanation
C  of any trouble except for a catastrophe when an explanation is provided
C  automatically.
C
C  Start of "CASE" statement implemented with COMPUTED GO TO.
      GO TO (60,60,80,100,120,140) CFLAG
C
   60 CONTINUE
C
C  CFLAG = 1 or 2.  Successful step.  Have we reached an output point?
      IF (TNOW.LT.TEND) THEN
C
C  Take another step.
         GO TO 40
      ELSE
C
C  Reached an output point.  Save the solution and reset the end point to
C  continue on. After MAXPER output points, jump to where the global error
C  assessment is obtained.
C
         WRITE (OUTFIL,'(2E13.5)') YNOW(1), YNOW(2)
         IF (KOUNTR.LT.MAXPER) THEN
            KOUNTR = KOUNTR + 1
            TEND = TEND + TWOPI
            CALL RESET(TEND)
            GO TO 40
         ELSE
            GO TO 160
         END IF
      END IF
C
   80 CONTINUE
C
C  CFLAG = 3.  Although a considerable amount of work has been done, continue.
      GO TO 40
C
  100 CONTINUE
C
C  CFLAG = 4.  The problem appears to be stiff.
      WRITE (*,'(A/)') 
     &' The problem appears to be stiff. '
      STOP
C
  120 CONTINUE
C
C  CFLAG = 5.  The accuracy request is too stringent.
      WRITE (*,'(A/)')
     &' The accuracy request is too stringent. '
      STOP
C
  140 CONTINUE
C
C  CFLAG = 6.  The global error estimation scheme broke down.  Report this
C              fact and then jump to where the assessment is obtained.
C
      WRITE (*,'(A/)') 
     &' Global error estimation broke down. '
      GO TO 160
C
  160 CONTINUE
C  End of "CASE" statement implemented with computed GO TO.
C
      CLOSE (OUTFIL)
C
C  Obtain the global error assessment by a call to GLBERR and report it. To
C  facilitate experimentation with the cost of global error assessment, test
C  whether ERRASS was set .TRUE. so that a call to GLBERR makes sense.
C
      IF (ERRASS) THEN
         CALL GLBERR(RMSERR,ERRMAX,TERRMX,WORK)
         WRITE (*,'(/A,1PE9.2)')
     &' The tolerance was ',TOL
         WRITE (*,'(A,1PE9.2,A,1PE9.2,A)')
     &' The worst global error observed was ',ERRMAX,
     &'.  (It occurred at ',TERRMX,'.)'
         WRITE (*,'(A)')
     &' The RMS errors in the individual components were:'
         DO 180 L = 1,NEQ
            WRITE (*,'(50X,1PE9.2)') RMSERR(L)
  180    CONTINUE
      END IF
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
