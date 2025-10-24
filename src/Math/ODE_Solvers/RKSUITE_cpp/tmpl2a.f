C
C                            Template 2a
C
C  This template comes in two variants, 2a and 2b, in files tmpl2a.for and
C  tmpl2b.for, respectively.  The variant 2b differs only in that the true,
C  or global, error of the integration in 2a is assessed.  The output to the
C  standard output channel from both variants run with both METHODs 2 and 3 
C  is found in the file tmpl2.out.
C
C  Problem:  Integrate a two body problem.  The equations for the coordinates
C            (x(t),y(t)) of one body as functions of time t in a suitable
C            frame of reference are
C
C               x'' = - x/r**3,
C               y'' = - y/r**3,   where   r = SQRT(x**2 + y**2).
C
C            The initial conditions lead to elliptic motion with eccentricity
C            ECC.  This parameter will be taken to be 0.9.
C
C               x(0) = 1-ECC,     x'(0) = 0,
C               y(0) = 0,         y'(0) = SQRT((1+ECC)/(1-ECC)).
C
C            An accurate solution that shows the general behavior of the
C            orbit is desired.  The coordinates will be returned at every
C            time step in [0,20].  This is a standard test problem for which
C            there is a semi-analytical solution.  It will be compared to the
C            computed solution at the end of the interval.
C
C  Solution: Substitute y1 = x, y2 = y, y3 = x', and y4 = y' to obtain
C
C               y1' = y3,            y1(0) = 1-ECC,
C               y2' = y4,            y2(0) = 0,
C               y3' = - y1/r**3,     y3(0) = 0,
C               y4' = - y2/r**3,     y4(0) = SQRT((1+ECC)/(1-ECC)),
C            where
C               r = SQRT(y1**2 + y2**2).
C
C            Since it is the general behavior of the solution that is desired,
C            it is best to let the integrator choose where to provide answers.
C            It will produce answers more frequently where the solution
C            changes more rapidly.  Because the solution is inspected at
C            every step, the task is a "Complex Task" that is solved with CT.
C
C            On physical grounds the solution is expected to be somewhat
C            unstable when one body approaches the other. To obtain an
C            accurate solution, a  stringent relative error tolerance should
C            be imposed -- TOL = 1.0D-10 will be used.  At a tolerance this
C            stringent the highest order pair, METHOD = 3, is likely to be
C            the most efficient choice. This method is inefficient when it
C            is to produce answers at a great many specific points. It is
C            most effective when used as in this template.  The solution
C            components are expected to be of order 1, so threshold values
C            THRES(*) = 1.0D-13 are reasonable.  When a solution component is
C            smaller in magnitude than this threshold, the code will control
C            the local error to be no more than TOL*THRES(L) = 1.0D-23.  The
C            reasonableness of this choice will be monitored by printing out
C            the maximum value seen for each solution component in the course
C            of the integration.  Error and warning messages, if any, will be
C            printed out.
C
C            This is the standard test problem D5 of T.E. Hull, W.H. Enright,
C            B.M. Fellen, and A.E. Sedgwick, "Comparing Numerical Methods
C            for Ordinary Differential Equations," SIAM Journal on Numerical
C            Analysis, Vol. 9, pp. 603-637, 1972.  The analytical solution
C            in terms of the numerical solution of Kepler's equation can be
C            found there as well as in most discussions of the two body
C            problem.  The results for the particular choice of eccentricity,
C            initial conditions, and interval of this template are provided
C            here in a DATA statement.
C
C  NOTES:    Typically problems are solved for a number of tolerances, initial
C            conditions, intervals, parameter values, ... .  A prudent person
C            would make at least one run with global error assessment as a
C            spot check on the reliability of the results.  Variant 2b shows
C            how to do this.
C
C            For TOL in this range, METHOD = 3 is generally the most efficient
C            choice.  Indeed, for this specific problem and tolerance (and a
C            specific computer, precision, compiler, ... ), the results found
C            in tmpl2.out show that the cost with METHOD = 3 is 4631 calls to
C            the subroutine F.  With METHOD = 2 the cost is 8091 calls, so
C            the higher order pair is quite advantageous at a tolerance as
C            stringent as TOL = 1.0D-10.  Of course, one should not even
C            consider using METHOD = 1 at tolerances this stringent.
C
C            With METHOD = 3, the template writes to an output file results
C            at 344 steps, which is more than adequate to see how the solution
C            components behave. The global error at the end of the run is about
C            1.2D-9, rather bigger than the local error tolerance TOL.  This
C            illustrates the point that at best one can anticipate global
C            errors comparable to the tolerance.  In point of fact, this
C            problem is unstable at some points of the integration and the
C            global error assessment of variant 2b reveals that the worst
C            global error is considerably worse than the error at the end 
C            -- an example of the value of the global error assessment
C            capability.
C
C            Working storage must be provided in the array WORK(*) of length
C            LENWRK.  Because storage is no problem with only NEQ = 4
C            equations, LENWRK is taken here to be 32*NEQ, enough to handle
C            all three methods with and without global error assessment.
C
C     .. Parameters ..
      INTEGER           NEQ, METHOD, LENWRK, OUTFIL
      PARAMETER         (NEQ=4,METHOD=3,LENWRK=32*NEQ,OUTFIL=9)
      DOUBLE PRECISION  ECC, ZERO, ONE
      PARAMETER         (ECC=0.9D0,ZERO=0.0D0,ONE=1.0D0)
C     .. Local Scalars ..
      DOUBLE PRECISION  ERROR, HNEXT, HSTART, TEND, TNOW,
     &                  TOL, TSTART, WASTE, WEIGHT
      INTEGER           CFLAG, CFLAG3, L, STPCST, STPSOK, TOTF
      LOGICAL           ERRASS, MESAGE
C     .. Local Arrays ..
      DOUBLE PRECISION  THRES(NEQ), TRUY(NEQ), WORK(LENWRK), YMAX(NEQ),
     &                  YNOW(NEQ), YPNOW(NEQ), YSTART(NEQ)
C     .. External Subroutines ..
      EXTERNAL          CT, F, SETUP, STAT
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SQRT
C     .. Data statements ..
      DATA              TRUY/-1.29526625098758D0,  0.400393896379232D0,
     &                       -0.67753909247075D0, -0.127083815427869D0/
C     .. Executable Statements ..
C
C  Initialize output.
C
      WRITE (*,'(A,I10/)') ' Template 2a with METHOD = ', METHOD
C
C  Set initial conditions and the end of the interval of interest.
C
      TSTART = ZERO
      YSTART(1) = ONE - ECC
      YSTART(2) = ZERO
      YSTART(3) = ZERO
      YSTART(4) = SQRT((ONE+ECC)/(ONE-ECC))
      TEND = 20.0D0
C
C  Because the solution components may be computed at many points, they are
C  written to a file on unit OUTFIL.  Subsequently, the results can be studied,
C  plotted, or manipulated as necessary.  Naming of files is system-dependent,
C  so the name here, RESULT, may need to be changed.  For simplicity, list-
C  directed output is used.  The data could be manipulated later as required
C  by, e.g., a plot package, or the WRITE statements in this template could be
C  altered so as to output the results in the desired form.  Begin by writing
C  the initial values to RESULT.
C
      OPEN (UNIT=OUTFIL,FILE='RESULT')
      WRITE (OUTFIL,*) TSTART, YSTART(1), YSTART(2), YSTART(3),
     &  YSTART(4)
C
C  To monitor the reasonableness of the choice of THRES(L), the maximum
C  magnitude YMAX(L) of solution component L is computed and returned.
C  It is initialized here.
C
      DO 20 L = 1, NEQ
         YMAX(L) = ABS(YSTART(L))
   20 CONTINUE
C
C  Set error control parameters.
C
      TOL = 1.0D-10
      DO 40 L = 1, NEQ
         THRES(L) = 1.0D-13
   40 CONTINUE
C
C  Call the setup routine. Because messages are requested, MESAGE = .TRUE.,
C  there is no need later to test values of flags and print out explanations.
C  In this variant no error assessment is done, so ERRASS is set .FALSE..
C  By setting HSTART = ZERO, the code is told to find a starting (initial)
C  step size automatically .
C
      MESAGE = .TRUE.
      ERRASS = .FALSE.
      HSTART = ZERO
      CALL SETUP(NEQ,TSTART,YSTART,TEND,TOL,THRES,METHOD,'Complex Task',
     &           ERRASS,HSTART,WORK,LENWRK,MESAGE)
C
C  Advance the integration one step.  Note that messages are written to 
C  the standard output channel rather than to the file where the results
C  are being accumulated.  The code will not go past the specified TEND.
C  Because it will step to exactly TEND, this fact can be used to terminate
C  the run.
C
      CFLAG3 = 0
   60 CONTINUE
      CALL CT(F,TNOW,YNOW,YPNOW,WORK,CFLAG)
C
C  Update the maximum magnitudes of the solution components.
      DO 80 L = 1, NEQ
         YMAX(L) = MAX(YMAX(L),ABS(YNOW(L)))
   80 CONTINUE
C
      IF (CFLAG.LE.3) THEN
C
C  Success.  Record the results and go on towards TEND.  The combination of
C  the tolerance and the length of the interval is sufficiently demanding that
C  the code might return with CFLAG = 3 to report that 5000 calls have been
C  made to the subroutine F.  The evaluations in F are not expensive, so
C  three such returns are permitted before terminating the run.
C
         WRITE (OUTFIL,*) TNOW, YNOW(1), YNOW(2), YPNOW(1), YPNOW(2)
         IF (CFLAG.EQ.3) CFLAG3 = CFLAG3 + 1
         IF (TNOW.LT.TEND .AND. CFLAG3.LT.3) GO TO 60
      END IF
      CLOSE (OUTFIL)
C
C  At the end of the run, the maximum magnitudes of the solution components
C  are reported so that the reasonableness of the threshold values THRES(*)
C  can be checked.
C
      WRITE (*,'(A/)') '             YMAX(L) '
      DO 100 L = 1, NEQ
         WRITE (*,'(13X,1PE8.2)') YMAX(L)
  100 CONTINUE
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
C  For this particular problem a semi-analytical solution is available to
C  study the relationship between local error tolerance and true (global)
C  error.  The true solution at t = 20 is provided in TRUY(*).  If the
C  integration reached this point, the error of the computed solution is
C  computed and returned.  To interpret the results properly, the error
C  must be measured in the same way that it is measured in the code.
C
      IF (TNOW.EQ.TEND) THEN
         ERROR = ZERO
         DO 120 L = 1, NEQ
            WEIGHT = MAX(ABS(YNOW(L)),THRES(L))
            ERROR = MAX(ERROR,ABS(TRUY(L)-YNOW(L))/WEIGHT)
  120    CONTINUE
         WRITE (*,'(/A,1PE9.2)')
     &' At t = 20, the error is',ERROR
         WRITE (*,'(A,1PE9.2)') ' The tolerance is       ', TOL
      END IF
C
      STOP
      END
      SUBROUTINE F(T,Y,YP)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
C     .. Array Arguments ..
      DOUBLE PRECISION  Y(*), YP(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  R
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      R = SQRT(Y(1)**2+Y(2)**2)
      YP(1) = Y(3)
      YP(2) = Y(4)
      YP(3) = -Y(1)/R**3
      YP(4) = -Y(2)/R**3
      RETURN
      END
