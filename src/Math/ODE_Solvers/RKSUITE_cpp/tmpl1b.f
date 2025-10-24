C
C                             Template 1b
C
C  This template comes in two variants, 1a and 1b, in files tmpl1a.for and
C  tmpl1b.for, respectively.  The variant 1b differs only in that the true,
C  or global, error of the integration in 1a is assessed.  Output from both
C  variants run with all three methods is found in the file tmpl1.out.
C
C  In this variant there is no need to repeat the statement of the problem
C  and the formulation of its solution since it differs only in that global
C  error assessment is done.  To emphasize what is different in this variant,
C  all comments not specifically related to the matter have been removed (a
C  practice not to be imitated).
C
C  Although you might well prefer an estimate of the global error in the
C  approximate solution at TWANT, and the code does have such an estimate at
C  its disposal, it is not made available to you because it is very difficult
C  to obtain a reliable estimate of the global error at a single point.  This
C  is easy to understand near a change of sign of the error, but at crude
C  tolerances and at limiting precision and when the problem is not smooth,
C  methods for estimation of global error may break down. Because an
C  assessment of the general size of the global error in the course of an
C  integration can be obtained with much greater reliability at acceptable
C  cost than an estimate at any one point, this is what is provided.  A
C  reliable global error assessment is relatively expensive, particularly for
C  METHOD = 1 because it is normally used only at crude tolerances and it
C  is difficult to get a reliable assessment then.  The assessment costs 
C  roughly double the cost of the integration itself for METHOD = 2 or 3, 
C  and triple for METHOD = 1.
C
C  After any return from UT with the integration advanced to T, the subroutine
C  GLBERR can be called to obtain an assessment of the global error.  Two kinds
C  of assessment are provided.  An array RMS(*) of length NEQ contains the
C  estimated root-mean-square error in each solution component.  This average
C  is over all solution values from the given initial value through the last
C  computed value. (UT may have advanced the integration past your last TWANT
C  and obtained a result there by interpolation, so the average can be over
C  values past TWANT.)  A scalar ERRMAX provides the largest error seen in any
C  solution component at any step of the integration so far.  The scalar
C  TERRMAX reports where the value ERRMAX was attained.  All errors are
C  measured with the same weights used in the integration so that they may
C  be compared with the local error tolerance TOL.
C
C  For the simple example problem of this template, the true solution is
C  readily available.  The output of variant 1a compares the computed to the
C  true solutions at selected points.  Examination of the file tmpl1.out
C  provided shows that the global error assessments of variant 1b are in
C  agreement with these true errors at individual output points.  With
C  METHOD = 2, the global error is quite closely related to the local error
C  tolerance TOL.  METHOD = 1 has a global error that is somewhat larger than
C  TOL.  METHOD = 3 has a global error that is considerably smaller than TOL.
C  Unlike the other two methods, this one does not have an interpolation
C  capability.  This particular problem is sufficiently easy for the method
C  and the output points are sufficiently close to one another that the step
C  size must be shortened to obtain solutions at the specified output points.
C  Reducing the step size for this reason results in more accuracy (and a
C  greater cost).
C
C  Global error assessment requires some storage.  Here the amount of
C  storage, LENWRK, allocated in WORK(*) is the same as in variant 1a
C  because it was taken there to be enough for all three methods with
C  and without global error assessment.
C
      INTEGER           NEQ, LENWRK, METHOD
      PARAMETER         (NEQ=2,LENWRK=32*NEQ,METHOD=2)
      DOUBLE PRECISION  ZERO, ONE, TWO, FOUR
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,FOUR=4.0D0)
      DOUBLE PRECISION  HNEXT, HSTART, PI, T, TEND, TINC,
     &                  TLAST, TOL, TSTART, TWANT, WASTE
      INTEGER           L, NOUT, STPCST, STPSOK, TOTF, UFLAG
      LOGICAL           ERRASS, MESAGE
      DOUBLE PRECISION  THRES(NEQ), WORK(LENWRK), Y(NEQ), YMAX(NEQ),
     &                  YP(NEQ), YSTART(NEQ)
C
C  Some additional quantities must be declared for global error assessment.
C
      DOUBLE PRECISION  ERRMAX, TERRMX
      DOUBLE PRECISION  RMSERR(NEQ)
C     
      EXTERNAL          F, SETUP, STAT, UT
      INTRINSIC         ATAN, COS, SIN
C
      TSTART = ZERO
      YSTART(1) = ZERO
      YSTART(2) = ONE
      PI = FOUR*ATAN(ONE)
      TLAST = TWO*PI
      TEND = TLAST + PI
C
      WRITE (*,'(A,I10)')  ' Template 1b with METHOD = ', METHOD
      WRITE (*,'(/A/)')    '    t           y        true y   '
      WRITE (*,'(1X,F6.3,3X,F9.4,3X,F9.4)')
     &                              TSTART, YSTART(1), SIN(TSTART)
      WRITE (*,'(1X,9X,F9.4,3X,F9.4/)') YSTART(2), COS(TSTART)
C
      TOL = 5.0D-5
      DO 20 L = 1, NEQ
         THRES(L) = 1.0D-10
   20 CONTINUE
      MESAGE = .TRUE.
C
C  In this variant error assessment is desired, so ERRASS is set .TRUE.
C  in the call to the setup routine.
C
      ERRASS = .TRUE.
      HSTART = ZERO
      CALL SETUP(NEQ,TSTART,YSTART,TEND,TOL,THRES,METHOD,'Usual Task',
     &           ERRASS,HSTART,WORK,LENWRK,MESAGE)
C
      NOUT = 8
      TINC = (TLAST-TSTART)/NOUT
C
      DO 40 L = 1, NOUT
         TWANT = TLAST + (L-NOUT)*TINC
         CALL UT(F,TWANT,T,Y,YP,YMAX,WORK,UFLAG)
C
C  A global error assessment is useless unless it is reliable, so the code
C  monitors computation of the assessment for credibility.  If the code has
C  any doubts about the reliability of the assessment, it will return  with
C  a message to this effect and UFLAG set to 6.  An attempt to continue
C  integrating after such a return is a fatal error.  As coded, the
C  integration is terminated when UFLAG = 6 and the global error is assessed
C  up to last point where the code believes the assessment to be reliable.
C
         IF (UFLAG.GT.2) GO TO 60
C
         WRITE (*,'(1X,F6.3,3X,F9.4,3X,F9.4)') T, Y(1), SIN(T)
         WRITE (*,'(1X,9X,F9.4,3X,F9.4/)')        Y(2), COS(T)
   40 CONTINUE
C
   60 CONTINUE
C
      WRITE (*,'(A/)') '             YMAX(L) '
      DO 80 L = 1, NEQ
         WRITE (*,'(13X,1PE8.2)')    YMAX(L)
   80 CONTINUE
C
C  The assessment is obtained by a call to GLBERR.
C
      CALL GLBERR(RMSERR,ERRMAX,TERRMX,WORK)
      WRITE (*,'(/A,1PE9.2)')
     &' The tolerance was ',TOL
      WRITE (*,'(A,1PE9.2,A,1PE9.2,A)')
     &' The worst global error observed was ',ERRMAX,
     &'.  (It occurred at ',TERRMX,'.)'
      WRITE (*,'(A)')
     &' The RMS errors in the individual components were:'
      DO 100 L = 1,NEQ
         WRITE (*,'(50X,1PE9.2)') RMSERR(L)
  100 CONTINUE
C
      CALL STAT(TOTF,STPCST,WASTE,STPSOK,HNEXT)
      WRITE (*,'(/A,I10)')
     &' The cost of the integration in evaluations of F was', TOTF
C
      STOP
      END
      SUBROUTINE F(T,Y,YP)
      DOUBLE PRECISION  T
      DOUBLE PRECISION  Y(*), YP(*)
      YP(1) =  Y(2)
      YP(2) = -Y(1)
      RETURN
      END
