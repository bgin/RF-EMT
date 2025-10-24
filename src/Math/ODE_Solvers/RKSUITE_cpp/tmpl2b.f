C
C                            Template 2b
C
C  This template comes in two variants, 2a and 2b, in files tmpl2a.for and
C  tmpl2b.for, respectively.  The variant 2b differs only in that the true,
C  or global, error of the integration in 2a is assessed.  The output to the
C  standard output channel from both variants run with both METHOD 2 and 3 is
C  found in the file tmpl2.out.
C
C  In this variant there is no need to repeat the statement of the problem
C  and the formulation of its solution since it differs only in that global
C  error assessment is done.  To emphasize what is different in this variant,
C  all comments not specifically related to the matter have been removed (a
C  practice not to be imitated).
C
C  Although you might well prefer an estimate of the global error in the
C  approximate solution at TNOW, and the code does have such an estimate at
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
C  After any return from CT with the integration advanced to TNOW, the
C  subroutine GLBERR can be called to obtain an assessment of the global error.
C  Two kinds of assessment are provided.  An array RMS(*) of length NEQ
C  contains the estimated root-mean-square error in each solution component.
C  This average is over all solution values from the given initial value
C  through the last computed value. A scalar ERRMAX provides the largest error
C  seen in any solution component at any step of the integration so far.  The
C  scalar TERRMX reports where the value ERRMAX was attained.  All errors are
C  measured with the same weights used in the integration so that they may
C  be compared with the local error tolerance TOL.
C
C  The output of variant 2a compares the computed to the true solution
C  components at the end of the interval of integration.  The error at this
C  point is found to be somewhat bigger than the tolerance.  However,
C  examination of the file tmpl2.out provided shows that the RMS errors are
C  rather bigger than the error at this one point, and moreover, the worst
C  global error is considerably bigger.  This illustrates the fact that
C  global errors can decrease as well as increase in the course of an
C  integration.
C
C  Global error assessment requires some storage.  Here the amount of
C  storage, LENWRK, allocated in WORK(*) is the same as in variant 2a
C  because it was taken there to be enough for all three methods with
C  and without global error assessment.
C
      INTEGER           NEQ, METHOD, LENWRK, OUTFIL
      PARAMETER         (NEQ=4,METHOD=3,LENWRK=32*NEQ,OUTFIL=9)
      DOUBLE PRECISION  ECC, ZERO, ONE
      PARAMETER         (ECC=0.9D0,ZERO=0.0D0,ONE=1.0D0)
      DOUBLE PRECISION  ERROR, HNEXT, HSTART, TEND, TNOW,
     &                  TOL, TSTART, WASTE, WEIGHT
      INTEGER           CFLAG, CFLAG3, L, STPCST, STPSOK, TOTF
      LOGICAL           ERRASS, MESAGE
      DOUBLE PRECISION  THRES(NEQ), TRUY(NEQ), WORK(LENWRK), YMAX(NEQ),
     &                  YNOW(NEQ), YPNOW(NEQ), YSTART(NEQ)
C
C  Some additional quantities must be declared for global error assessment.
C
      DOUBLE PRECISION  ERRMAX,TERRMX
      DOUBLE PRECISION  RMSERR(NEQ)
C
      EXTERNAL          CT, F, SETUP, STAT
      INTRINSIC         ABS, MAX, SQRT
      DATA              TRUY/-1.29526625098758D0,  0.400393896379232D0,
     &                       -0.67753909247075D0, -0.127083815427869D0/
C
      WRITE (*,'(A,I10/)') ' Template 2b with METHOD = ', METHOD
C
      TSTART = ZERO
      YSTART(1) = ONE - ECC
      YSTART(2) = ZERO
      YSTART(3) = ZERO
      YSTART(4) = SQRT((ONE+ECC)/(ONE-ECC))
      TEND = 20.0D0
C
      OPEN (UNIT=OUTFIL,FILE='RESULT')
      WRITE (OUTFIL,*) TSTART, YSTART(1), YSTART(2), YSTART(3),
     &  YSTART(4)
C
      DO 20 L = 1, NEQ
         YMAX(L) = ABS(YSTART(L))
   20 CONTINUE
C
      TOL = 1.0D-10
      DO 40 L = 1, NEQ
         THRES(L) = 1.0D-13
   40 CONTINUE
C
      MESAGE = .TRUE.
C
C  In this variant error assessment is desired, so ERRASS is set .TRUE. in
C  the call to the setup routine.
C
      ERRASS = .TRUE.
      HSTART = ZERO
      CALL SETUP(NEQ,TSTART,YSTART,TEND,TOL,THRES,METHOD,'Complex Task',
     &           ERRASS,HSTART,WORK,LENWRK,MESAGE)
C
      CFLAG3 = 0
   60 CONTINUE
      CALL CT(F,TNOW,YNOW,YPNOW,WORK,CFLAG)
C
      DO 80 L = 1, NEQ
         YMAX(L) = MAX(YMAX(L),ABS(YNOW(L)))
   80 CONTINUE
C
C  A global error assessment is useless unless it is reliable, so the code
C  monitors the computation of the assessment for credibility.  If the code
C  has any doubts about the reliability of the assessment, it will return
C  with a message to this effect and CFLAG set to 6.  An attempt to continue
C  integrating after such a return is a fatal error.  When CFLAG = 6 the
C  integration is terminated and the global error is assessed through the
C  last point where the code believes the assessment to be reliable.
C
      IF (CFLAG.LE.3) THEN
C
         WRITE (OUTFIL,*) TNOW, YNOW(1), YNOW(2), YPNOW(1), YPNOW(2)
         IF (CFLAG.EQ.3) CFLAG3 = CFLAG3 + 1
         IF (TNOW.LT.TEND .AND. CFLAG3.LT.3) GO TO 60
      END IF
      CLOSE (OUTFIL)
C
      WRITE (*,'(A/)') '             YMAX(L) '
      DO 100 L = 1, NEQ
         WRITE (*,'(13X,1PE8.2)') YMAX(L)
  100 CONTINUE
C
      CALL STAT(TOTF,STPCST,WASTE,STPSOK,HNEXT)
      WRITE (*,'(/A,1PE10.2,0P/A,I10/A,I10/A,F10.2/A,I10)')
     &' The integration reached                      ', TNOW,
     &' The cost of the integration in calls to F was', TOTF,
     &' The number of calls to F per step is         ', STPCST,
     &' The fraction of failed steps was             ', WASTE,
     &' The number of accepted steps was             ', STPSOK
C
      IF (TNOW.EQ.TEND) THEN
         ERROR = ZERO
         DO 120 L = 1, NEQ
            WEIGHT = MAX(ABS(YNOW(L)),THRES(L))
            ERROR = MAX(ERROR,ABS(TRUY(L)-YNOW(L))/WEIGHT)
  120    CONTINUE
         WRITE (*,'(/A,1PE9.2)') 
     &' At t = 20, the error is', ERROR
         WRITE (*,'(A,1PE9.2)')  ' The tolerance is       ', TOL
      END IF
C
C  The assessment is obtained by a call to GLBERR.
C
      CALL GLBERR(RMSERR,ERRMAX,TERRMX,WORK)
      WRITE (*,'(A,1PE9.2,A,1PE9.2,A)')
     &' The worst global error observed was ',ERRMAX,
     &'.  (It occurred at ',TERRMX,'.)'
      WRITE (*,'(A)')
     &' The RMS errors in the individual components were:'
      DO 140 L = 1,NEQ
         WRITE (*,'(50X,1PE9.2)') RMSERR(L)
  140 CONTINUE
C
      STOP
      END
      SUBROUTINE F(T,Y,YP)
      DOUBLE PRECISION  T
      DOUBLE PRECISION  Y(*), YP(*)
      DOUBLE PRECISION  R
      INTRINSIC         SQRT
      R = SQRT(Y(1)**2+Y(2)**2)
      YP(1) = Y(3)
      YP(2) = Y(4)
      YP(3) = -Y(1)/R**3
      YP(4) = -Y(2)/R**3
      RETURN
      END
