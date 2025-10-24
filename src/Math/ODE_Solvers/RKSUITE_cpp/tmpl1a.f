C
C                             Template 1a
C
C  This template comes in two variants, 1a and 1b, in files tmpl1a.for and
C  tmpl1b.for, respectively.  The variant 1b differs only in that the true,
C  or global, error of the integration in 1a is assessed.  Output from both
C  variants run with all three methods is found in the file tmpl1.out.
C
C  Problem:  Compute about four correct figures in the solution of
C
C               y'' = -y
C
C            on the range [0,2*pi] at intervals of length pi/4, given
C            that y(0)=0,  y'(0)=1.
C
C  Solution: Let y1 = y and y2 = y' to obtain the first order system
C
C               y1' =   y2     with initial values   y1 = 0
C               y2' = - y1                           y2 = 1
C
C            This is a "Usual Task" that is appropriately solved with UT.
C            Although the code controls the local error rather than the true
C            error, it is "tuned" so that the true error will be comparable to
C            the local error tolerance for typical problems. Thus a relative
C            error tolerance TOL of 5.0D-5 is appropriate.  In this range of
C            tolerances, METHOD = 2 is likely to be the most efficient choice.
C            The solution components are expected to get as large as 1.0D0.
C            With this in mind, solution components smaller than, say, 1.0D-10
C            are not very interesting, and requiring five correct figures then
C            is not worth the cost. For this reason, the threshold values are
C            specified as THRES(L) = 1.0D-10 for L = 1,2.  When solution
C            component L is smaller than this threshold, the code will control
C            the local error to be no more than TOL*THRES(L) = 5.0D-15.  Error
C            and warning messages, if any, will be printed out.  Answers will
C            be computed at equally spaced points, and the true values of y
C            and y' will be printed out for comparison.
C
C  NOTES:    Typically problems are solved for a number of tolerances, initial
C            conditions, intervals, parameter values, ... .  A prudent person
C            would make at least one run with global error assessment as a
C            spot check on the reliability of the results.  Variant 1b shows
C            how to do this.
C
C            For TOL in this range, METHOD = 2 is generally the most efficient
C            choice.  Indeed, for this specific problem and tolerance (and a
C            specific computer, precision, compiler, ... ), the results found
C            in tmpl1.out show that the cost with METHOD = 2 is 109 calls to
C            the subroutine F, with METHOD = 1 it is 292 calls, and with
C            METHOD = 3 it is 105 calls. At relaxed tolerances, METHOD = 1 is
C            generally the most efficient choice, and at stringent tolerances,
C            METHOD = 3 is generally the most efficient.
C
C            In typical use of UT, the cost is scarcely affected by the
C            number of answers, but when a "large" number of answers is
C            required, this is not true of METHOD = 3.  In such a situation
C            its cost is proportional to the number of answers.
C
C            Working storage must be provided in the array WORK(*) of length
C            LENWRK.  Because storage is no problem with only NEQ = 2
C            equations, LENWRK is taken here to be 32*NEQ, enough to handle
C            all three methods with and without global error assessment.
C
C     .. Parameters ..
      INTEGER           NEQ, LENWRK, METHOD
      PARAMETER         (NEQ=2,LENWRK=32*NEQ,METHOD=2)
      DOUBLE PRECISION  ZERO, ONE, TWO, FOUR
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,FOUR=4.0D0)
C     .. Local Scalars ..
      DOUBLE PRECISION  HNEXT, HSTART, PI, T, TEND, TINC,
     &                  TLAST, TOL, TSTART, TWANT, WASTE
      INTEGER           L, NOUT, STPCST, STPSOK, TOTF, UFLAG
      LOGICAL           ERRASS, MESAGE
C     .. Local Arrays ..
      DOUBLE PRECISION  THRES(NEQ), WORK(LENWRK), Y(NEQ), YMAX(NEQ),
     &                  YP(NEQ), YSTART(NEQ)
C     .. External Subroutines ..
      EXTERNAL          F, SETUP, STAT, UT
C     .. Intrinsic Functions ..
      INTRINSIC         ATAN, COS, SIN
C     .. Executable Statements ..
C
C  Set the initial conditions.  Note that TEND is taken well past
C  the last output point, TLAST.  When this is possible, and it
C  usually is, it is good practice.
C
      TSTART = ZERO
      YSTART(1) = ZERO
      YSTART(2) = ONE
      PI = FOUR*ATAN(ONE)
      TLAST = TWO*PI
      TEND = TLAST + PI
C
C  Initialize output.
C
      WRITE (*,'(A,I10)')  ' Template 1a with METHOD = ', METHOD
      WRITE (*,'(/A/)')    '    t           y        true y   '
      WRITE (*,'(1X,F6.3,3X,F9.4,3X,F9.4)')
     &                              TSTART, YSTART(1), SIN(TSTART)
      WRITE (*,'(1X,9X,F9.4,3X,F9.4/)') YSTART(2), COS(TSTART)
C
C  Set error control parameters.
C
      TOL = 5.0D-5
      DO 20 L = 1, NEQ
         THRES(L) = 1.0D-10
   20 CONTINUE
C
C  Call the setup routine. Because messages are requested, MESAGE = .TRUE.,
C  there is no need later to test values of flags and print out explanations.
C  In this variant no error assessment is done, so ERRASS is set .FALSE..
C  By setting HSTART to zero, the code is told to find a starting (initial)
C  step size automatically .
C
      MESAGE = .TRUE.
      ERRASS = .FALSE.
      HSTART = ZERO
      CALL SETUP(NEQ,TSTART,YSTART,TEND,TOL,THRES,METHOD,'Usual Task',
     &           ERRASS,HSTART,WORK,LENWRK,MESAGE)
C
C  Compute answers at NOUT equally spaced output points. It is good
C  practice to code the calculation of TWANT so that the last value
C  is exactly TLAST.
C
      NOUT = 8
      TINC = (TLAST-TSTART)/NOUT
C
      DO 40 L = 1, NOUT
         TWANT = TLAST + (L-NOUT)*TINC
         CALL UT(F,TWANT,T,Y,YP,YMAX,WORK,UFLAG)
C
         IF (UFLAG.GT.2) GO TO 60
C
C  Success. T = TWANT. Output computed and true solution components.
         WRITE (*,'(1X,F6.3,3X,F9.4,3X,F9.4)') T, Y(1), SIN(T)
         WRITE (*,'(1X,9X,F9.4,3X,F9.4/)')        Y(2), COS(T)
   40 CONTINUE
C
C  The integration is complete or has failed in a way reported in a
C  message to the standard output channel.
   60 CONTINUE
C
C  YMAX(L) is the largest magnitude computed for the solution component
C  Y(L) in the course of the integration from TSTART to the last T.  It
C  is used to decide whether THRES(L) is reasonable and to select a new
C  THRES(L) if it is not.
C
      WRITE (*,'(A/)') '             YMAX(L) '
      DO 80 L = 1, NEQ
         WRITE (*,'(13X,1PE8.2)')    YMAX(L)
   80 CONTINUE
C
C  The subroutine STAT is used to obtain some information about the progress
C  of the integration. TOTF is the total number of calls to F made so far
C  in the integration; it is a machine-independent measure of work.  At present
C  the integration is finished, so the value printed out refers to the overall
C  cost of the integration.
C
      CALL STAT(TOTF,STPCST,WASTE,STPSOK,HNEXT)
      WRITE (*,'(/A,I10)')
     &  ' The cost of the integration in evaluations of F is', TOTF
C
      STOP
      END
      SUBROUTINE F(T,Y,YP)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
C     .. Array Arguments ..
      DOUBLE PRECISION  Y(*), YP(*)
C     .. Executable Statements ..
      YP(1) =  Y(2)
      YP(2) = -Y(1)
      RETURN
      END
